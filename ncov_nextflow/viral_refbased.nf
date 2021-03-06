#!/usr/bin/env nextflow

println("Currently using the Viral Reference Based Assembly workflow for use with artic-Illumina hybrid library prep on MiSeq")
println("Version: v.20210510")

//# To be used with the ivar container staphb/ivar:1.2.2_artic20200528, this includes all artic and reference files, plus the index files already exist

params.artic_version = 'V3'
params.year = '2021'
params.reads = workflow.launchDir + '/Sequencing_reads/Raw'

maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${maxcpus}")
if ( maxcpus < 5 ) {
  medcpus = maxcpus
} else {
  medcpus = 5
}

// param that coincides with the staphb/seqyclean:1.10.09 container run with singularity
params.seqyclean_contaminant_file="/Adapters_plus_PhiX_174.fasta"

// params that coincide with the staphb/ivar:1.2.2_artic20200528 container run with singularity
// when not using the container, the reference genome will need to be indexed for use with bwa
params.primer_bed = file("/home/workflows/nCoV19_pipeline/nCoV-2019.primer.V3.bed")
params.reference_genome = file("/artic-ncov2019/primer_schemes/nCoV-2019/${params.artic_version}/nCoV-2019.reference.fasta")
params.gff_file = file("/reference/GCF_009858895.2_ASM985889v3_genomic.gff")
params.amplicon_bed = file("/artic-ncov2019/primer_schemes/nCoV-2019/${params.artic_version}/nCoV-2019_amplicon.bed")

// param that coincides with the staphb/kraken2:2.0.8-beta_hv container run with singularity
params.kraken2_db="/kraken2-db"

// This is where the results will be
params.outdir = workflow.launchDir
println("The files and directory for results is " + params.outdir)
params.log_directory = params.outdir + '/logs'
println("A table summarizing results will be created: ${workflow.launchDir}/run_results.txt")

// this sample file contains metadata for renaming files and adding collection dates to submission files
// The columns should be file_name\tsubmission_id\tcollection_date
params.sample_file = file(params.outdir + '/covid_samples.txt' )

samples = []
if (params.sample_file.exists()) {
  println("List of COVID19 samples: " + params.sample_file)
  params.sample_file
    .readLines()
    .each { samples << it.split('\t')[0] }
  }
  else {
    println("${params.sample_file} could not be found!")
    println("\tFor file submission renaming, please include a file named\n\t${params.outdir}/covid_samples.txt\n\twith columns of sample_id\tsubmission_id\tcollection_date")
    }

Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq.gz",
                  "${params.reads}/*_{1,2}.fastq*",
                  "${params.reads}/*.R{1,2}.fastq*"], size: 2 )
  .map{ reads -> [reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]] }
  .ifEmpty{
    println("No paired fastq or fastq.gz files were found at ${params.reads}")
    exit 1
  }
  .into { fastq_reads; fastq_reads2; fastq_reads3; fastq_reads4 }

process ids {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/bedtools logs/bedtools'

  input:
  set val(sample), file(reads) from fastq_reads3

  when:
  if (params.sample_file.exists()) { return true }

  output:
  tuple sample, env(sample_id), env(submission_id), env(collection_date) into submission_ids, submission_ids2

  shell:
  '''
    log_file=logs/bedtools/multicov.!{workflow.sessionId}.log
    err_file=logs/bedtools/multicov.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    sample_id='NA'
    submission_id='NA'
    collection_date='NA'

    while read line
    do
      lab_accession=$(echo $line | awk '{print $1}' )
      if [[ "!{sample}" == *"$lab_accession"* ]]
      then
        sample_id=$lab_accession
        submission_id=$(echo $line | awk '{print $2}' )
        collection_date=$(echo $line | awk '{print $3}' )
      fi
    done < !{params.sample_file}

    if [ -z "$sample_id" ] ; then sample_id=!{sample} ; fi
    if [ -z "$submission_id" ] ; then submission_id=!{sample} ; fi
    if [ -z "$collection_date" ] ; then collection_date='NA' ; fi
  '''
}

process seqyclean {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p Sequencing_reads/QCed logs/seqyclean'

  input:
  set val(sample), file(reads) from fastq_reads

  output:
  tuple sample, file("Sequencing_reads/QCed/${sample}_clean_PE{1,2}.fastq") into clean_reads, clean_reads2, clean_reads3
  file("Sequencing_reads/QCed/${sample}_clean_SE.fastq")
  file("Sequencing_reads/QCed/${sample}_clean_SummaryStatistics.{txt,tsv}")
  file("logs/seqyclean/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(pairskept), env(perc_kept) into seqyclean_results

  shell:
  '''
    log_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.log
    err_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "seqyclean version: $(seqyclean -h | grep Version)" >> $log_file

    seqyclean -minlen 25 -qual -c !{params.seqyclean_contaminant_file} -1 !{reads[0]} -2 !{reads[1]} -o Sequencing_reads/QCed/!{sample}_clean 2>> $err_file >> $log_file
    pairskept=$(cut -f 58 Sequencing_reads/QCed/!{sample}_clean_SummaryStatistics.tsv | grep -v "PairsKept" | head -n 1)
    perc_kept=$(cut -f 59 Sequencing_reads/QCed/!{sample}_clean_SummaryStatistics.tsv | grep -v "Perc_Kept" | head -n 1)
    if [ -z "$pairskept" ] ; then pairskept="0" ; fi
    if [ -z "$perc_kept" ] ; then perc_kept="0" ; fi
  '''
}

process bwa {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus maxcpus

  beforeScript 'mkdir -p covid/bwa logs/bwa_covid'

  input:
  set val(sample), file(reads) from clean_reads

  output:
  tuple sample, file("covid/bwa/${sample}.sorted.bam") into bams, bams2, bams3
  file("covid/bwa/${sample}.sorted.bam.bai") into bais
  file("logs/bwa_covid/${sample}.${workflow.sessionId}.log")
  file("logs/bwa_covid/${sample}.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/bwa_covid/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bwa_covid/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
    samtools --version >> $log_file

    # bwa mem command
    bwa mem -t !{maxcpus} !{params.reference_genome} !{reads[0]} !{reads[1]} 2>> $err_file | \
      samtools sort 2>> $err_file | \
      samtools view -F 4 -o covid/bwa/!{sample}.sorted.bam 2>> $err_file >> $log_file

    # indexing the bams
    samtools index covid/bwa/!{sample}.sorted.bam 2>> $err_file >> $log_file
  '''
}

process ivar_trim {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/ivar_trim logs/ivar_trim'

  input:
  set val(sample), file(bam) from bams

  output:
  tuple sample, file("covid/ivar_trim/${sample}.primertrim.bam") into trimmed_bams
  file("logs/ivar_trim/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    ivar version >> $log_file

    wget https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.primer.bed

    # trimming the reads
    ivar trim -e -i !{bam} -b nCoV-2019.primer.bed -p covid/ivar_trim/!{sample}.primertrim 2>> $err_file >> $log_file
  '''
}

process samtools_sort {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/trimmed logs/samtools_sort_trimmed'

  input:
  set val(sample), file(bam) from trimmed_bams

  output:
  tuple sample, file("covid/trimmed/${sample}.primertrim.sorted.bam") into sorted_bams, sorted_bams2, sorted_bams3, sorted_bams4
  file("covid/trimmed/${sample}.primertrim.sorted.bam.bai") into sorted_bais
  file("logs/samtools_sort_trimmed/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/samtools_sort_trimmed/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_sort_trimmed/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    # sorting and indexing the trimmed bams
    samtools sort !{bam} -o covid/trimmed/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
    samtools index covid/trimmed/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
  '''
}

process ivar_variants {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/variants logs/ivar_variants'

  input:
  set val(sample), file(bam) from sorted_bams

  output:
  file("covid/variants/${sample}.variants.tsv")
  file("logs/ivar_variants/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(variants_num) into variants_results

  shell:
  '''
    log_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file

    samtools mpileup -A -d 600000 -B -Q 0 --reference !{params.reference_genome} !{bam} 2>> $err_file | \
      ivar variants -p covid/variants/!{sample}.variants -q 20 -t 0.6 -m 10 -r !{params.reference_genome} -g !{params.gff_file} 2>> $err_file >> $log_file

    variants_num=$(grep "TRUE" covid/variants/!{sample}.variants.tsv | wc -l)
    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
  '''
}

process ivar_consensus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/consensus logs/ivar_consensus'

  input:
  set val(sample), file(bam) from sorted_bams2

  output:
  tuple sample, file("covid/consensus/${sample}.consensus.fa") into consensus, consensus2, consensus3
  file("logs/ivar_consensus/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(num_N), env(num_ACTG), env(num_degenerate), env(num_total) into consensus_results, consensus_results2

  shell:
  '''
    log_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file

    samtools mpileup -A -d 6000000 -B -Q 0 --reference !{params.reference_genome} !{bam} 2>> $err_file | \
      ivar consensus -t 0.6 -p covid/consensus/!{sample}.consensus -n N 2>> $err_file >> $log_file

    num_N=$(grep -v ">" covid/consensus/!{sample}.consensus.fa | grep -o 'N' | wc -l)
    if [ -z "$num_N" ] ; then num_N="0" ; fi

    num_ACTG=$(grep -v ">" covid/consensus/!{sample}.consensus.fa | grep -o -E "C|A|T|G" | wc -l)
    if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi

    num_degenerate=$(grep -v ">" covid/consensus/!{sample}.consensus.fa | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z"  | wc -l)
    if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi

    num_total=$(grep -v ">" covid/consensus/!{sample}.consensus.fa | grep -o -E '[A-Z]' | wc -l)

    sed -i 's/Consensus_//;s/.consensus_threshold.*//' covid/consensus/!{sample}.consensus.fa
  '''
}

fastq_reads2
  .combine(clean_reads2, by: 0)
  .set { raw_clean_reads }

process fastqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  beforeScript 'mkdir -p fastqc logs/fastqc'

  input:
  set val(sample), file(raw), file(clean) from raw_clean_reads

  output:
  file("fastqc/*.{html,zip}")
  tuple sample, env(raw_1), env(raw_2), env(clean_1), env(clean_2) into fastqc_results
  file("logs/fastqc/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/fastqc/!{sample}.!{workflow.sessionId}.log
    err_file=logs/fastqc/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastqc --version >> $log_file

    fastqc --outdir fastqc --threads !{task.cpus} !{sample}*.fastq* 2>> $err_file >> $log_file

    raw_1=$(unzip -p fastqc/!{raw[0].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    raw_2=$(unzip -p fastqc/!{raw[1].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    clean_1=$(unzip -p fastqc/!{clean[0].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    clean_2=$(unzip -p fastqc/!{clean[1].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
  '''
}

bams3
  .combine(sorted_bams3, by: 0)
  .into { combined_bams; combined_bams2; combined_bams3 }

process samtools_stats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/samtools_stats/bwa covid/samtools_stats/trimmed logs/samtools_stats'

  input:
  set val(sample), file(bam), file(sorted_bam) from combined_bams

  output:
  file("covid/samtools_stats/bwa/${sample}.stats.txt")
  file("covid/samtools_stats/trimmed/${sample}.stats.trim.txt")
  file("logs/samtools_stats/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/samtools_stats/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_stats/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools stats !{bam} > covid/samtools_stats/bwa/!{sample}.stats.txt 2>> $err_file
    samtools stats !{sorted_bam} > covid/samtools_stats/trimmed/!{sample}.stats.trim.txt 2>> $err_file
  '''
}

process samtools_coverage {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/samtools_coverage/bwa covid/samtools_coverage/trimmed logs/samtools_coverage'

  input:
  set val(sample), file(bwa), file(sorted) from combined_bams2

  output:
  file("covid/samtools_coverage/bwa/${sample}.cov.{txt,hist}")
  file("covid/samtools_coverage/trimmed/${sample}.cov.trim.{txt,hist}")
  file("logs/samtools_coverage/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(coverage), env(depth), env(coverage_trim), env(depth_trim), env(meanbaseq_trim), env(meanmapq_trim) into samtools_coverage_results

  shell:
  '''
    log_file=logs/samtools_coverage/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_coverage/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools coverage !{bwa} -m -o covid/samtools_coverage/bwa/!{sample}.cov.hist 2>> $err_file >> $log_file
    samtools coverage !{bwa} -o covid/samtools_coverage/bwa/!{sample}.cov.txt 2>> $err_file >> $log_file
    samtools coverage !{sorted} -m -o covid/samtools_coverage/trimmed/!{sample}.cov.trim.hist 2>> $err_file >> $log_file
    samtools coverage !{sorted} -o covid/samtools_coverage/trimmed/!{sample}.cov.trim.txt 2>> $err_file >> $log_file

    coverage=$(cut -f 6 covid/samtools_coverage/bwa/!{sample}.cov.txt | tail -n 1)
    depth=$(cut -f 7 covid/samtools_coverage/bwa/!{sample}.cov.txt | tail -n 1)
    if [ -z "$coverage" ] ; then coverage="0" ; fi
    if [ -z "$depth" ] ; then depth="0" ; fi

    coverage_trim=$(cut -f 6 covid/samtools_coverage/trimmed/!{sample}.cov.trim.txt | tail -n 1)
    depth_trim=$(cut -f 7 covid/samtools_coverage/trimmed/!{sample}.cov.trim.txt | tail -n 1)
    meanbaseq_trim=$(cut -f 8 covid/samtools_coverage/trimmed/!{sample}.cov.trim.txt | tail -n 1)
    meanmapq_trim=$(cut -f 9 covid/samtools_coverage/trimmed/!{sample}.cov.trim.txt | tail -n 1)
    if [ -z "$coverage_trim" ] ; then coverage_trim="0" ; fi
    if [ -z "$depth_trim" ] ; then depth_trim="0" ; fi
    if [ -z "$meanbaseq_trim" ] ; then meanbaseq_trim="0" ; fi
    if [ -z "$meanmapq_trim" ] ; then meanmapq_trim="0" ; fi
  '''
}

process samtools_flagstat {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/samtools_flagstat/bwa covid/samtools_flagstat/trimmed logs/samtools_flagstat'

  input:
  set val(sample), file(bwa), file(trimmed) from combined_bams3

  output:
  file("covid/samtools_flagstat/{bwa,trimmed}/${sample}.flagstat.txt")
  file("logs/samtools_flagstat/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/samtools_flagstat/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_flagstat/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools flagstat !{bwa} > covid/samtools_flagstat/bwa/!{sample}.flagstat.txt 2>> $err_file
    samtools flagstat !{trimmed} > covid/samtools_flagstat/trimmed/!{sample}.flagstat.txt 2>> $err_file
  '''
}

process kraken2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus maxcpus

  beforeScript 'mkdir -p covid/kraken2 logs/kraken2'

  input:
  set val(sample), file(clean) from clean_reads3

  output:
  file("covid/kraken2/${sample}_kraken2_report.txt")
  file("logs/kraken2/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(percentage_human), env(percentage_cov) into kraken2_results

  shell:
  '''
    log_file=logs/kraken2/!{sample}.!{workflow.sessionId}.log
    err_file=logs/kraken2/!{sample}.!{workflow.sessionId}.err
    date | tee -a $log_file $err_file > /dev/null
    kraken2 --version >> $log_file
    kraken2 --paired \
      --classified-out cseqs#.fq \
      --threads !{task.cpus} \
      --db !{params.kraken2_db} \
      !{clean[0]} !{clean[1]} \
      --report covid/kraken2/!{sample}_kraken2_report.txt \
      2>> $err_file >> $log_file
    percentage_human=$(grep "Homo sapiens" covid/kraken2/!{sample}_kraken2_report.txt | awk '{print $1}')
    percentage_cov=$(grep "Severe acute respiratory syndrome coronavirus 2" covid/kraken2/!{sample}_kraken2_report.txt | awk '{print $1}')
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    if [ -z "$percentage_cov" ] ; then percentage_cov="0" ; fi
  '''
}

//kraken2_results2.view()


process bedtools {
  publishDir "${params.outdir}", mode: 'copy'
  tag "bedtools"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/bedtools logs/bedtools'

  input:
  file(bams) from bams2.collect()
  file(bais) from bais.collect()
  file(trimmed_bams) from sorted_bams4.collect()
  file(trimmed_bais) from sorted_bais.collect()

  output:
  file("covid/bedtools/multicov.txt") into bedtools_results
  file("logs/bedtools/multicov.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/bedtools/multicov.!{workflow.sessionId}.log
    err_file=logs/bedtools/multicov.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    bedtools --version >> $log_file

    echo "primer" $(ls *bam) | tr ' ' '\t' > covid/bedtools/multicov.txt
    bedtools multicov -bams $(ls *bam) -bed !{params.amplicon_bed} | cut -f 4,6- 2>> $err_file >> covid/bedtools/multicov.txt
  '''
}

process pangolin {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus medcpus

  beforeScript 'mkdir -p covid/pangolin logs/pangolin'

  input:
  set val(sample), file(fasta) from consensus

  output:
  file("covid/pangolin/${sample}/lineage_report.csv")
  tuple sample, env(pangolin_lineage), env(pangolin_stats), env(pangolin_vers) into pangolin_results
  file("logs/pangolin/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/pangolin/!{sample}.!{workflow.sessionId}.log
    err_file=logs/pangolin/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    pangolin --version >> $log_file
    #pangolin -lv >> $log_file
    pangolin -pv >> $log_file

    pangolin_vers=$(pangolin --version | sed 's/ /-/')_$(pangolin -pv | sed 's/ /-/')

    pangolin --threads !{task.cpus} --outdir covid/pangolin/!{sample} !{fasta} 2>> $err_file >> $log_file
    pangolin_lineage=$(tail -n 1 covid/pangolin/!{sample}/lineage_report.csv | cut -f 2 -d "," | grep -v "lineage" )

    pangolin_stats=$(tail -n 1 covid/pangolin/!{sample}/lineage_report.csv | cut -f 3 -d "," )
    if [ -z "$pangolin_stats" ]; then pangolin_stats="None"; fi
  '''
}

process nextclade {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus medcpus

  beforeScript 'mkdir -p covid/nextclade logs/nextclade'

  input:
  set val(sample), file(fasta) from consensus3

  output:
  file("covid/nextclade/${sample}.nextclade.tsv")
  file("covid/nextclade/${sample}.nextclade.auspice.json")
  file("covid/nextclade/${sample}.nextclade.json")
  tuple sample, env(nextclade_lineage), env(aa_subs), env(aa_dels) into nextclade_results
  file("logs/nextclade/${sample}.${workflow.sessionId}.{log,err}")
  file("covid/nextclade/${sample}.nextcladeout.txt")

  shell:
  '''
    log_file=logs/nextclade/!{sample}.!{workflow.sessionId}.log
    err_file=logs/nextclade/!{sample}.!{workflow.sessionId}.err

    nextclade.js --version >> $log_file
    nextclade.js \
      --input-fasta !{fasta} \
      --output-json !{sample}.nextclade.json \
      --output-tsv  !{sample}.nextclade.tsv \
      --output-tree !{sample}.nextclade.auspice.json

    cp !{sample}.nextclade.json covid/nextclade/!{sample}.nextclade.json
    cp !{sample}.nextclade.tsv covid/nextclade/!{sample}.nextclade.tsv
    cp !{sample}.nextclade.auspice.json covid/nextclade/!{sample}.nextclade.auspice.json

    lineage_col=$(head -n 1 !{sample}.nextclade.tsv | tr '\t' '\n' | grep -n clade | cut -f 1 -d ":" | head -n 1 )
    nextclade_lineage=$(cut -f $lineage_col !{sample}.nextclade.tsv | grep -v clade)
    if [ -z "$nextclade_lineage" ]; then nextclade_lineage="none"; fi

    aa_subs_col=$(head -n 1 !{sample}.nextclade.tsv | tr '\t' '\n' | grep -n aaSubstitutions | cut -f 1 -d ":" | head -n 1 )
    aa_subs=$(cut -f $aa_subs_col !{sample}.nextclade.tsv | grep -v aaSubstitutions | sed 's/,/|/g')
    if [ -z "$aa_subs" ]; then aa_subs="none";fi

    aa_dels_col=$(head -n 1 !{sample}.nextclade.tsv | tr '\t' '\n' | grep -n aaDeletions | cut -f 1 -d ":" | head -n 1 )
    aa_dels=$(cut -f $aa_dels_col !{sample}.nextclade.tsv | grep -v aaDeletions | sed 's/,/|/g')
    if [ -z "$aa_dels" ]; then aa_dels="none";fi

    echo $nextclade_lineage >> covid/nextclade/!{sample}.nextcladeout.txt
    echo $aa_dels >> covid/nextclade/!{sample}.nextcladeout.txt
    echo $aa_subs >> covid/nextclade/!{sample}.nextcladeout.txt
  '''
}
seqyclean_results
  // tuple sample, env(pairskept), env(perc_kept) into seqyclean_results
  .combine(variants_results, by: 0)
  // tuple sample, env(variants_num) into variants_results
  .combine(consensus_results, by: 0)
  // tuple sample, env(num_N), env(num_ACTG), env(num_degenerate), env(num_total) into consensus_results
  .combine(fastqc_results, by: 0)
  // tuple sample, env(raw_1), env(raw_2), env(clean_1), env(clean_2) into fastqc_results
  .combine(samtools_coverage_results, by: 0)
  // tuple sample, env(coverage), env(depth), env(coverage_trim), env(depth_trim), env(meanbaseq_trim), env(meanmapq_trim) into samtools_coverage_results
  .combine(kraken2_results, by: 0)
  // tuple sample, env(percentage_human), env(percentage_cov) into kraken2_results
  .combine(pangolin_results, by: 0)
  // tuple sample, env(pangolin_lineage), env(pangolin_stats), env(pangolin_vers) into pangolin_results
  .combine(nextclade_results, by: 0)
  // tuple sample, env(nextclade_lineage), env(aa_subs), env(aa_dels) into nextclade_results
  .combine(submission_ids2, by: 0)
  // tuple sample, env(sample_id), env(submission_id), env(collection_date)
  .set { results }

process summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/summary logs/summary'

  input:
  set val(sample), val(pairskept), val(perc_kept),
    val(variants_num),
    val(num_N), val(num_ACTG), val(num_degenerate), val(num_total),
    val(raw_1), val(raw_2), val(clean_1), val(clean_2),
    val(coverage), val(depth), val(coverage_trim), val(depth_trim),
    val(meanbaseq_trim), val(meanmapq_trim),
    val(percentage_human), val(percentage_cov),
    val(pangolin_lineage), val(pangolin_stats), val(pangolin_vers),
    val(nextclade_lineage), val(aa_subs), val(aa_dels),
    val(sample_id), val(submission_id), val(collection_date) from results
  file(multicov) from bedtools_results.collect()

  output:
  file("covid/summary/${sample}.summary.txt") into summary
  file("logs/summary/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/summary/!{sample}.!{workflow.sessionId}.log
    err_file=logs/summary/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    bedtools_column=$(head -n 1 !{multicov} | tr '\t' '\n' | grep -n !{sample} | grep -v primertrim | cut -f 1 -d ":" | head -n 1 )
    amp_fail=$(cut -f $bedtools_column !{multicov} | awk '{ if ( $1 < 20 ) print $0 }' | wc -l )
    if [ -z "$amp_fail" ] ; then amp_fail=0 ; fi

    sample_id=$(echo !{sample} | cut -f 1 -d "_" )

    cov=$(echo "!{coverage} >= 95.00" | bc)
    if (($cov == 0)); then cov_status="Avg coverage < 95%;"; else cov_status="";fi

    baseq=$(echo "!{meanbaseq_trim} >= 30.00" | bc)
    if (($baseq == 0)); then baseq_status="Mean base quality < 30;"; else baseq_status=""; fi

    mapq=$(echo "!{meanmapq_trim} >= 30.00" | bc)
    if (($mapq == 0)); then mapq_status="Mean map quality < 30"; else mapq_status=""; fi

    #if (($cov == 1)) && (($baseq == 1)) && (($mapq == 1)); then assembly_status="PASS"; else assembly_status=$(echo "WARNING: ${cov_status}${baseq_status}${mapq_status}"); fi
    if (($cov == 1)) && (($baseq == 1)) && (($mapq == 1)); then assembly_status="PASS"; else assembly_status="FAIL"; fi


    seq_date=$(date +%Y-%m-%d)

    echo -e "sample_id,deidentified_id,collection_date,seq_date,pangolin_lineage,pangolin_conflict,nextclade_lineage,AA_substitutions,AA_deletions,fastqc_raw_reads_1,fastqc_raw_reads_2,fastqc_clean_reads_PE1,fastqc_clean_reads_PE2,pairs_kept_after_cleaning,percent_kept_after_cleaning,depth_before_trimming,depth_after_trimming,coverage_before_trimming,coverage_after_trimming,%_human_reads,%_SARS-COV-2_reads,num_failed_amplicons,num_variants,num_N,num_degenerate,num_ACTG,num_total,meanbaseq_trim,meanmapq_trim,assembly_status,pango_version" > covid/summary/!{sample}.summary.txt
    echo -e "${sample_id},!{submission_id},!{collection_date},${seq_date},!{pangolin_lineage},!{pangolin_stats},!{nextclade_lineage},!{aa_subs},!{aa_dels},!{raw_1},!{raw_2},!{clean_1},!{clean_2},!{pairskept},!{perc_kept},!{depth},!{depth_trim},!{coverage},!{coverage_trim},!{percentage_human},!{percentage_cov},$amp_fail,!{variants_num},!{num_N},!{num_degenerate},!{num_ACTG},!{num_total},!{meanbaseq_trim},!{meanmapq_trim},$assembly_status,!{pangolin_vers}" >> covid/summary/!{sample}.summary.txt

  '''
}

process combined_summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "summary"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/submission_files logs/summary'

  input:
  file(summary) from summary.collect()

  output:
  file("covid/summary.txt")
  file("run_results.csv")
  file("logs/summary/summary.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/summary/summary.!{workflow.sessionId}.log
    err_file=logs/summary/summary.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    cat *summary.txt | grep "assembly_status" | uniq > covid/summary.txt 2>> $err_file
    cat *summary.txt | grep -v "pangolin_lineage" | sort | uniq >> covid/summary.txt 2>> $err_file

    cp covid/summary.txt run_results.csv
  '''
}

fastq_reads4
  .combine(consensus2, by: 0)
  // tuple sample, file("covid/consensus/${sample}.consensus.fa") into consensus, consensus2
  .combine(consensus_results2, by:0)
  // tuple sample, env(num_N), env(num_ACTG), env(num_degenerate), env(num_total) into consensus_results
  .combine(submission_ids, by:0)
  // tuple val(sample), env(sample_id), env(submission_id), env(collection_date) into submission_ids
  .set{ reads_consensus }

process file_submission {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/submission_files logs/submission'

  input:
  set val(sample), file(reads),
    file(consensus),
    val(num_n), val(num_ACTG), val(num_degen), val(num_total),
    val(sample_id), val(submission_id), val(collection_date) from reads_consensus

  when:
  for(int i =0; i < samples.size(); i++) {
    if(sample.contains(samples[i])) { return true }
  }

  output:
  file("covid/submission_files/${submission_id}.{R1,R2}.fastq.gz")
  file("covid/submission_files/${submission_id}.consensus.fa")
  file("covid/submission_files/${submission_id}.{genbank,gisaid}.fa") optional true into submission_files
  file("logs/submission/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
  log_file=logs/submission/!{sample}.!{workflow.sessionId}.log
  err_file=logs/submission/!{sample}.!{workflow.sessionId}.err

  date | tee -a $log_file $err_file > /dev/null

  # getting the consensus fasta file
  # changing the fasta header
  echo ">!{submission_id}" > covid/submission_files/!{submission_id}.consensus.fa 2>> $err_file
  grep -v ">" !{consensus} >> covid/submission_files/!{submission_id}.consensus.fa 2>> $err_file

  year=$(echo !{collection_date} | cut -f 1 -d '-')

  SPEC_STATE_ID=$(echo !{submission_id} | cut -f1 -d '-')
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

  if [ "!{num_n}" -lt 15000 ] && [ "!{num_total}" -gt 25000 ]
  then
    # removing leading Ns, folding sequencing to 75 bp wide, and adding metadata for genbank submissions
    echo ">!{submission_id} [organism=Severe acute respiratory syndrome coronavirus 2][isolate=SARS-CoV-2/Human/USA/!{submission_id}/$year][host=Human][country=USA:$SPEC_STATE][collection_date=!{collection_date}]" > covid/submission_files/!{submission_id}.genbank.fa  2>> $err_file
    grep -v ">" !{consensus} | sed 's/^N*N//g' | fold -w 75 >> covid/submission_files/!{submission_id}.genbank.fa  2>> $err_file
    if [ "!{num_ACTG}" -gt 25000 ]
    then
      echo ">hCoV-19/USA/!{submission_id}/$year" > covid/submission_files/!{submission_id}.gisaid.fa
      grep -v ">" !{consensus} >> covid/submission_files/!{submission_id}.gisaid.fa
      echo "!{sample} had !{num_n} Ns and is part of the genbank and gisaid submission fasta" >> $log_file
    else
      echo "!{sample} had !{num_n} Ns and is part of the genbank submission fasta, but not gisaid" >> $log_file
    fi
  else
    echo "!{sample} had !{num_n} Ns and is not part of the genbank or the gisaid submission fasta" >> $log_file
  fi

  # copying fastq files and changing the file name
  cp !{reads[0]} covid/submission_files/!{submission_id}.R1.fastq.gz  2>> $err_file
  cp !{reads[1]} covid/submission_files/!{submission_id}.R2.fastq.gz  2>> $err_file
  '''
}

process multifasta_submission {
  publishDir "${params.outdir}", mode: 'copy'
  tag "multifasta"
  echo false
  cpus 1

  beforeScript 'mkdir -p covid/submission_files logs/multifasta_submission'

  input:
  file(fastas) from submission_files.collect()

  when:
  if (params.sample_file.exists()) { return true }

  output:
  file("covid/submission_files/*.{gisaid_submission,genbank_submission}.fasta")
  file("logs/multifasta_submission/multifasta_submission.${workflow.sessionId}.{log,err}")

  shell:
  '''
  log_file=logs/multifasta_submission/multifasta_submission.!{workflow.sessionId}.log
  err_file=logs/multifasta_submission/multifasta_submission.!{workflow.sessionId}.err

  date | tee -a $log_file $err_file > /dev/null

  run_id=$(echo "!{params.outdir}" | rev | cut -f 1 -d '/' | rev )
  if [ -z "$run_id" ] ; then run_id="Submission" ; fi
  cat *gisaid.fa > covid/submission_files/$run_id.gisaid_submission.fasta 2>> $err_file
  cat *genbank.fa > covid/submission_files/$run_id.genbank_submission.fasta 2>> $err_file
  '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("A summary of results can be found in a tab-delimited file: ${workflow.launchDir}/run_results.txt")
    if (params.sample_file.exists()) { println("SRA, GenBank, and GISAID submission-ready files are located at ${workflow.launchDir}/covid/submission_files") }
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
