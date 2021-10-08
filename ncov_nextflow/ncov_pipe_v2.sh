# Last edit : 20210512
# FROM WITHIN THE "batch_folder"

NEXTFLOW_PIPE='/home/workflows/nCoV19_pipeline/viral_refbased.nf'
NEXTFLOW_PIPE_CONFIG='/home/workflows/nCoV19_pipeline/viral_refbased.singularity.config'
NEXTFLOW_CLUSTER='/home/workflows/nCoV19_pipeline/cluster_analysis.nf'
NEXTFLOW_CLUSTER_CONFIG='/home/workflows/nCoV19_pipeline/cluster_analysis.singularity.config'
# -----------------------------------------------------------------------------------
# MERGING DUPLICATE-SEQUENCED READS AND MOVE READS TO NECESSARY LOCATION FOR PIPELINE
# -----------------------------------------------------------------------------------

# We often sequence in duplicate to bolster genome coverage & quality. This takes all
# samples sequenced in duplicate (-1 and -2) and merges the L and R read files.
# echo -e "\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \nConcatenating L and R reads of samples sequenced in duplicate\n"

READ_LOC='Sequencing_reads/Raw'
# remove undetermined reads
undet_reads=$(ls Undetermined* 2> /dev/null | wc -l)
if [ $undet_reads != 0 ]; then rm Undetermined*; fi

# Test if reads are present
num_reads_cwd=$(ls *R1*fastq.gz 2> /dev/null | wc -l)
num_reads_readdir=$(ls $READ_LOC/*R1*fastq.gz 2> /dev/null | wc -l)
if [ $num_reads_cwd != 0 ]; then
        mkdir -p $READ_LOC
        # if duplicate reads are found, concatenate them
        num_dup_reads=$(ls *-[1-4]*fastq.gz 2> /dev/null | wc -l)
        if [ $num_dup_reads != 0 ]; then
                mkdir -p pre-merged_reads
                for read in $(echo R1 && echo R2); do for samp in $(ls *-[1-4]_*fastq.gz | sed 's/-.*//' | sort -u); do cat $samp-*$read*fastq.gz > ${samp}_S00_L001_${read}_001.fastq.gz; done; done
                mv *-[1-4]_*fastq.gz pre-merged_reads
        fi
        # move all remaining reads to "reads" directory
        num_reads_cwd=$(ls *R1*fastq.gz 2> /dev/null | wc -l)
        if [ $num_reads_cwd != 0 ]; then mkdir -p $READ_LOC && mv *fastq.gz $READ_LOC; fi
        echo -e "\nConcatenation complete\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n"
elif [ $num_reads_readdir != 0 ]; then
        echo -e "\nReads found in '$READ_LOC' directory. Continuing with analysis."
else
        echo -e "\n\n********************* ERROR: No reads found in current directory *********************\n\n"
fi

num_reads_readdir=$(ls $READ_LOC/*R1*fastq.gz 2> /dev/null | wc -l)
if [ $num_reads_readdir != 0 ]; then
        # -----------------------------------------------------------------------------------------
        # MOVE READS TO NECESSARY LOCATION FOR PIPELINE AND GENERATE INPUT "covid_samples.txt" FILE
        # -----------------------------------------------------------------------------------------

        accession=$(head -n1 Metadata_*csv | tr ',' '\n' | grep -n "Accession" | cut -f 1 -d ':' | head -n1)
        masphl_ID=$(head -n1 Metadata_*csv | tr ',' '\n' | grep -n "MASPHL ID" | cut -f 1 -d ':' | head -n1)
        collection_date=$(head -n1 Metadata_*csv | tr ',' '\n' | grep -n "Collection" | cut -f 1 -d ':' | head -n1)
        cut -f $accession,$masphl_ID,$collection_date -d ',' Metadata*.csv | grep -v "Accession*" | sed "s/,/\t/g" > covid_samples.txt

        # ---------------------------
        # GENERATE CONSENSUS SEQUENCE
        # ---------------------------
        # Here, we use an adjusted version of the Cecret pipeline from UPHL. Following, all
        # passing assemblies are moved to the "passing_consensus" folder for a cluster analysis.
        echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
        echo
        echo "Starting pipeline for consensus calling"
        echo
        nextflow run $NEXTFLOW_PIPE -c $NEXTFLOW_PIPE_CONFIG
        mkdir -p passing_consensus
        for i in $(grep PASS run_results.csv | awk -F ',' '{print $1}'); do cp covid/consensus/$i.consensus.fa passing_consensus/$i.consensus.fasta; done
        echo
        echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
fi
# ---------------------------
# RUN A CLUSTER ANALYSIS
# ---------------------------
# echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
# echo
# echo "Starting pipeline for cluster analysis"
# echo

# if [ -f run_results.csv ]; then
#       echo -e "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \nStarting pipeline for cluster analysis\n"
#     STATUS_COLUMN=$(head -n 1 run_results.csv | tr ',' '\n' | grep -n assembly_status | cut -f 1 -d ":" | head -n 1 )
#       cut -f 1,$STATUS_COLUMN -d ',' run_results.csv | sed 's/,/\t/g' > analysis_run_cluster.tsv
#       nextflow run $NEXTFLOW_CLUSTER -c $NEXTFLOW_CLUSTER_CONFIG
#       echo -e "\nCluster analysis complete\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
# else
#       echo -e "\nNo run_results.csv file found from viral pipeline. Cannot perform cluster analysis.\n"
# fi

# ------------------------------
# GENERATE AGGREGATED DATA FILES
# ------------------------------

if [ -f run_results.csv ] && [ -f Metadata*.csv ]; then

        echo "Found run_results and Metadata files. Aggregating data for BtB & MAVEN"
        # ------------
        # Capture accession by provider for cluster analyses
        # -----------

        # generating list of facilities that requested sequencing, and extracting specimens for each facility (outputs to add_to_cluster_FACILITY.txt)
        if grep -q Sequencing Metadata*csv; then
                mkdir -p sample_by_provider
                IFS=$'\n'
                PROVIDER_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n Provider | cut -f 1 -d ":" | head -n 1 )
                for i in $(grep Sequencing Metadata*csv | cut -f $PROVIDER_COL -d ',' | sort -u); do grep "$i" Metadata*csv | cut -f1 -d ',' > sample_by_provider/add_to_cluster_`echo $i | sed "s/ /_/g;s/'//g;s/\//-/g"`.txt; done
        fi

        # -------------------------------------------
        # Generate Headers for GISAID & GENBANK files
        # -------------------------------------------

        mkdir -p pub_repo_files
        cp covid/submission_files/*fasta pub_repo_files

        GISAID_OUTFILE=$(echo pub_repo_files/`pwd | rev | cut -f1 -d '/' | rev`.gisaid_metadata.csv)
    echo submitter,fn,covv_virus_name,covv_type,covv_passage,covv_collection_date,covv_location,covv_add_location,covv_host,covv_add_host_info,covv_sampling_strategy,covv_gender,covv_patient_age,covv_patient_status,covv_specimen,covv_outbreak,covv_last_vaccinated,covv_treatment,covv_seq_technology,covv_assembly_method,covv_coverage,covv_orig_lab,covv_orig_lab_addr,covv_provider_sample_id,covv_subm_lab,covv_subm_lab_addr,covv_subm_sample_id,covv_authors,covv_comment,comment_type > $GISAID_OUTFILE
    echo Submitter,FASTA filename,Virus name,Type,Passage details/history,Collection date,Location,Additional location information,Host,Additional host information,Sampling Strategy,Gender,Patient age,Patient status,Specimen source,Outbreak,Last vaccinated,Treatment,Sequencing technology,Assembly method,Coverage,Originating lab,Address,Sample ID given by the sample provider,Submitting lab,Address,Sample ID given by the submitting laboratory,Authors,Comment,Comment Icon >> $GISAID_OUTFILE

        GENBANK_OUTFILE=$(echo pub_repo_files/`pwd | rev | cut -f1 -d '/' | rev`.genbank_metadata.csv)
        echo Sequence_ID,Organism,collection-date,country,host,isolate,isolation-source,BioProject > $GENBANK_OUTFILE

        # Attributes for submission to public repositories
        SUBMITTER='Alang'
        BIOPROJECT='PRJNA686883'
        PASSAGE='Original'
        CONTINENT='North America'
        COUNTRY='USA'
        HOST='Human'
        GENDER='unknown'
        AGE='unknown'
        PT_STATUS='unknown'
        SEQ_PLATFORM='Illumina Miseq'
        ASSEMBLY_METHOD='BWA 0.7.17; iVar 1.2.2'
        ORIG_LAB='Massachusetts State Public Health Laboratory'
        LAB_ADDY='305 South St, Boston, MA 02130'
        SUBMITTING_LAB='Massachusetts State Public Health Laboratory'
        SUBMITTING_ADDY='305 South St, Boston, MA 02130'
        AUTHORS='Andrew Lang, Timelia Fink, Glen Gallagher, Sandra Smole'
        GENBANK_ORGANISM='Severe acute respiratory syndrome coronavirus 2'
        FASTA_Name=`basename covid/submission_files/*gisaid_submission.fasta`
        # ----------------------------
        # Generate BtB and MAVEN files
        # ----------------------------

        IFS=$'\n'

        # Determining column numbers of Metadata file
        COLLECTION_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n Collection | cut -f 1 -d ":" | head -n 1 )
        MASPHL_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n MASPHL | cut -f 1 -d ":" | head -n 1 )
        REASON_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n Reason | cut -f 1 -d ":" | head -n 1 )
        PROVIDER_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n Provider | cut -f 1 -d ":" | head -n 1 )
        LASTNM_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n LastName | cut -f 1 -d ":" | head -n 1 )
        FIRSTNM_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n FirstName | cut -f 1 -d ":" | head -n 1 )
        DOB_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n "Birth" | cut -f 1 -d ":" | head -n 1 )
        STREET_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n "PatientAddress" | cut -f 1 -d ":" | head -n 1 )
        CITY_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n "City" | cut -f 1 -d ":" | head -n 1 )
        STATE_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n "State" | cut -f 1 -d ":" | head -n 1 )
        ZIP_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n "PostalCode" | cut -f 1 -d ":" | head -n 1 )
        COUNTY_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n "County" | cut -f 1 -d ":" | head -n 1 )
        SOURCE_COL=$(head -n 1 Metadata*csv | tr ',' '\n' | grep -n "Source" | cut -f 1 -d ":" | head -n 1 )

        # determining column numbers of run_results file
        ACCESSION_COLUMN=$(head -n 1 run_results.csv | tr ',' '\n' | grep -n sample_id | cut -f 1 -d ":" | head -n 1 )
        STATUS_COLUMN=$(head -n 1 run_results.csv | tr ',' '\n' | grep -n assembly_status | cut -f 1 -d ":" | head -n 1 )
        PANGO_COLUMN=$(head -n1 run_results.csv | tr ',' '\n' | grep -n pangolin_lineage | cut -f 1 -d ':' | head -n 1)
        PANGO_VERS_COLUMN=$(head -n1 run_results.csv | tr ',' '\n' | grep -n pango_version | cut -f 1 -d ':' | head -n 1)
        SeqDATE_COL=$(head -n1 run_results.csv | tr ',' '\n' | grep -n seq_date | cut -f 1 -d ':' | head -n 1)
        AA_SUBS_COL=$(head -n1 run_results.csv | tr ',' '\n' | grep -n AA_substitutions | cut -f 1 -d ':' | head -n 1)
        AA_DELS_COL=$(head -n1 run_results.csv | tr ',' '\n' | grep -n AA_deletions | cut -f 1 -d ':' | head -n 1)

        # FLAT FILE FOR MAVEN (before BtB and MAVEN can talk)
        MAVEN_OUTFILE=$(echo NGS_SC2_MAVEN_`date +%Y-%m-%d`.csv)
        echo 'Accession,MASPHL_ID,CollectionDate,Reason For Sequencing,Provider,Last_Name,First_Name,Date_of_Birth,Street,City,State,ZipCode,County,AssemblyStatus,Lineage,MAVEN_Lineage' > $MAVEN_OUTFILE

        # FILE FOR UPLOAD TO BtB
        BtB_OUTFILE=$(echo BtB_upload_`pwd | rev | cut -f1 -d '/' | rev`_`date +%Y-%m-%d`.csv)
#       ORGANISM="SARS-Coronavirus 2"
        ORGANISM="SARS-CoV 2"
        BtB_TEST="SARS-CoV-2 Sequencing"
        echo 'sample_id,assembly_status,tool_lineage,lineage_to_maven,pango_version,organism,test' > $BtB_OUTFILE

        # FILE FOR TABLEAU DASHBOARD
        TABLEAU_OUTFILE=$(echo NGS_tableau_vis_`date +%Y-%m-%d`.csv)
        echo 'sample_id,assembly_status,collection_date,seq_date,lineage,County,ZipCode,Reason For Sequencing,AA_substitutions,AA_deletions' > $TABLEAU_OUTFILE

        # Merge metadata & results aspects for each specimen in results file
        for id in `cut -f $ACCESSION_COLUMN -d ',' run_results.csv | grep VX`; do

                # grab metadata aspects for specimen - I grab with a trailing comma, to ensure it's not grabbing 21VX111 when I want 21VX11
                # added the sort -u aspect to handle duplicate entries in metadta file
        COLLECTION=$(grep $id, Metadata*csv | cut -f $COLLECTION_COL -d ',' | sed "s/\t/,/g" | sort -u)
        MASPHL=$(grep $id, Metadata*csv | cut -f $MASPHL_COL -d ',' | sed "s/\t/,/g"| sort -u)
        REASON=$(grep $id Metadata*csv | cut -f $REASON_COL -d ',' | sed "s/\t/,/g")
        PROVIDER=$(grep $id, Metadata*csv | cut -f $PROVIDER_COL -d ',' | sed "s/\t/,/g"| sort -u)
        LASTNM=$(grep $id, Metadata*csv | cut -f $LASTNM_COL -d ',' | sed "s/\t/,/g"| sort -u)
        FIRSTNM=$(grep $id, Metadata*csv | cut -f $FIRSTNM_COL -d ',' | sed "s/\t/,/g"| sort -u)
        DOB=$(grep $id, Metadata*csv | cut -f $DOB_COL -d ',' | sed "s/\t/,/g" | sort -u)
        STREET=$(grep $id, Metadata*csv | cut -f $STREET_COL -d ',' | sed "s/\t/,/g"| sort -u)
        CITY=$(grep $id, Metadata*csv | cut -f $CITY_COL -d ',' | sed "s/\t/,/g"| sort -u)
        STATE=$(grep $id, Metadata*csv | cut -f $STATE_COL -d ',' | sed "s/\t/,/g"| sort -u)

        # Access cannot output leading 0 for zip code. Auto adding here
        # If zip contains only digits
        ZIP=$(grep $id, Metadata*csv | cut -f $ZIP_COL -d ',' | sed "s/\t/,/g" | sort -u)
        if [[ $ZIP =~ ^[0-9]+$ ]]; then
                #if zip is not already 5 digits, add leading 0
                zip_chars=$(awk -F '[0-9]' '{print NF-1}' <<< "$ZIP")
                if [[ $zip_chars == 4 ]]; then
                        ZIP=`echo 0$ZIP`; fi
                else ZIP="";fi

        COUNTY=$(grep $id, Metadata*csv | cut -f $COUNTY_COL -d ',' | sed "s/\t/,/g" | sort -u | tr [a-z] [A-Z])
                # grab results aspects for specimen
                STATUS=$(grep $id, run_results.csv | cut -f $STATUS_COLUMN -d ',' | tail -n1 | sort -u)
                PANGO_LINEAGE=$(grep $id, run_results.csv | cut -f $PANGO_COLUMN -d ',' | tail -n1 | sort -u)
                PANGO_VERSION=$(grep $id, run_results.csv | cut -f $PANGO_VERS_COLUMN -d ',' | tail -n1 | sort -u)

                SEQ_DATE=$(grep $id, run_results.csv | cut -f $SeqDATE_COL -d ',' | tail -n1 | sort -u)
                AA_SUBS=$(grep $id, run_results.csv | cut -f $AA_SUBS_COL -d ',' | tail -n1 | sort -u)
                AA_DELS=$(grep $id, run_results.csv | cut -f $AA_DELS_COL -d ',' | tail -n1 | sort -u)

                # Capture all VOCs/VOIs as defined by the CDC, all else captured as "other lineage"
                # ---VOCs---            ---------VOIs---------
                # B.1.1.7                       B.1.525          B.1.617.1
                # B.1.351                       B.1.526          B.1.617.2
                # P.1                           B.1.526.1        B.1.617.3
                # B.1.427                       P.2
                # B..1.429                      B.1.617

                if [[ $STATUS == PASS ]]; then
                        STATUS_OUT="PASS"
                        if [[ $PANGO_LINEAGE == "B.1.1.7" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.351" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "P.1" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "P.2" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.427" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.429" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.526" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.526.1" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.525" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.617" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.617.1" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.617.2" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        elif [[ $PANGO_LINEAGE == "B.1.617.3" ]]; then
                                LINEAGE_TO_MAVEN=$PANGO_LINEAGE

                        # elif [[ $PANGO_LINEAGE == "B.1.351.2" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "B.1.351.3" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.1" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.2" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.3" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.4" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.5" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.6" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.7" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.8" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.9" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.10" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.11" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "AY.12" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "P.1.1" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE
                        # elif [[ $PANGO_LINEAGE == "P.1.2" ]]; then
                        #       LINEAGE_TO_MAVEN=$PANGO_LINEAGE

                        else
                                LINEAGE_TO_MAVEN="Other Lineage"
                        fi




                # Only adjust lineage informaition for data passed to maven when a specimen fails sequencing
                elif [[ $STATUS == FAIL ]]; then
                STATUS_OUT="FAIL"
                LINEAGE_TO_MAVEN="INVALID"
                PANGO_LINEAGE="INVALID"
                else
                        STATUS_OUT="UNKNOWN"
                        LINEAGE_TO_MAVEN="UNKNOWN"
                        PANGO_LINEAGE="UNKNOWN";
                fi

                # add specimen data to MAVEN file
                echo $id,$MASPHL,$COLLECTION,$REASON,$PROVIDER,$LASTNM,$FIRSTNM,$DOB,$STREET,$CITY,$STATE,$ZIP,$COUNTY,$STATUS_OUT,$PANGO_LINEAGE,$LINEAGE_TO_MAVEN >> $MAVEN_OUTFILE;

                # add specimen data to BtB file
                echo $id,$STATUS_OUT,$PANGO_LINEAGE,$LINEAGE_TO_MAVEN,$PANGO_VERSION,$ORGANISM,$BtB_TEST >> $BtB_OUTFILE;

                # add specimen data to Tableau file
                echo $id,$STATUS_OUT,$COLLECTION,$SEQ_DATE,$PANGO_LINEAGE,$COUNTY,$ZIP,,$AA_SUBS,$AA_DELS >> $TABLEAU_OUTFILE

        # Specimen-specific attributes for submission to public repository
        SOURCE=$(grep $id Metadata*csv | cut -f $SOURCE_COL -d ',' | sed "s/\t/,/g" | sort -u)
                COLLECTION_YEAR=$(echo $COLLECTION | cut -f1 -d '-')
                ISOLATE=$(echo SARS-CoV-2/$HOST/$COUNTRY/$MASPHL/$COLLECTION_YEAR)
                SPEC_STATE_ID=$(echo $MASPHL | cut -f1 -d '-')
                if [[ $SPEC_STATE_ID == 'MA' ]];then
                        SPEC_STATE='Massachusetts'
                elif [[ $SPEC_STATE_ID == 'NH' ]];then
                        SPEC_STATE='New Hampshire'
                elif [[ $SPEC_STATE_ID == 'VT' ]];then
                        SPEC_STATE='Vermont'
                elif [[ $SPEC_STATE_ID == 'RI' ]];then
                        SPEC_STATE='Rhode Island'
                elif [[ $SPEC_STATE_ID == 'CT' ]];then
                        SPEC_STATE='Connecticut'
                else SPEC_STATE=''
                fi

                if [ -f covid/submission_files/$MASPHL.gisaid.fa ]; then
                        #echo "$SUBMITTER","$MASPHL.gisaid.fa","hCoV-19/$COUNTRY/$MASPHL/$COLLECTION_YEAR","betacoronavirus","$PASSAGE","$COLLECTION","$CONTINENT / $COUNTRY / $STATE",,"$HOST",,,"$GENDER","$AGE","$STATUS","$SOURCE",,,,"SEQ_PLATFORM","ASSEMBLY_METHOD",,"$ORIG_LAB","$LAB_ADDY",,"$SUBMITTING_LAB","$SUBMITTING_ADDY",,"$AUTHORS" >> $GISAID_OUTFILE
                        echo '"'$SUBMITTER'"','"'$FASTA_Name'"','"'hCoV-19/$COUNTRY/$MASPHL/$COLLECTION_YEAR'"','"'betacoronavirus'"','"'$PASSAGE'"','"'$COLLECTION'"','"'$CONTINENT / $COUNTRY / $SPEC_STATE'"',,'"'$HOST'"',,,'"'$GENDER'"','"'$AGE'"','"'$PT_STATUS'"','"'$SOURCE'"',,,,'"'$SEQ_PLATFORM'"','"'$ASSEMBLY_METHOD'"',,'"'$ORIG_LAB'"','"'$LAB_ADDY'"',,'"'$SUBMITTING_LAB'"','"'$SUBMITTING_ADDY'"',,'"'$AUTHORS'"',, >> $GISAID_OUTFILE
                fi

                if [ -f covid/submission_files/$MASPHL.genbank.fa ]; then
                        echo '"'$MASPHL'"','"'$GENBANK_ORGANISM'"','"'$COLLECTION'"','"'$COUNTRY':'$SPEC_STATE'"','"'$HOST'"','"'$ISOLATE'"','"'$SOURCE'"','"'$BIOPROJECT'"' >> $GENBANK_OUTFILE
                        #echo "$MASPHL",$ORGANISM","$COLLECTION","$COUNTRY","$HOST","$ISOLATE","$SOURCE","$BIOPROJECT" >> $GENBANK_OUTFILE
                fi
        done
else
        echo "Cannot find run_results.csv file. It appears the pipeline did not complete successfully."
fi
