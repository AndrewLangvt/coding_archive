#!/usr/bin/env python3

# brr_maineffor.py
# Andrew S. Lang
# Created: 28OCT2019
# Last Modified: 10JAN2019

import os
import re
import sys
import shutil
import glob
import argparse
import traceback
import subprocess
import brr_foundation
from shutil import which
from linecache import getline
from time import gmtime, strftime
from brr_foundation import make_directory
from brr_foundation import remove
from brr_foundation import run_command_logger
from brr_foundation import logger
from brr_foundation import confirm_complete
from brr_foundation import confirm_present
import time
import datetime
from subprocess import Popen, PIPE, STDOUT
sys.path.append(os.path.abspath(os.path.dirname(__file__)))                     # appending location of github repo to sys path to pull in other python modules

parser = argparse.ArgumentParser(description='Paramters for running pipeline')
parser = argparse.ArgumentParser(
    epilog='''QC AND ORGANISMAL TYPING PIPELINE. This pipeline is maintained by Andrew Lang
    at the Massachusetts State Public Health Lab. Additional information pertaining to pipeline usage can
    be found at https://gitlab.com/ma_ngs/brr_pipeline. For any issues with program operation, please
    contact Andrew at Andrew.Lang@massmail.state.ma.us for assistance.''')
parser._action_groups.pop()
optional = parser.add_argument_group('Optional Arugments')
optional.add_argument('-cfg', '--config', help='location of config file (Default= ./brr_config.ini)', default='./brr_config.ini')

def pull_SRRs(d, p):
    """Pulls SRA files from NCBI based upon IDs provided in a file called SRR"""

    logger('SEQUENCE FILE LOCATION\n')
    logger('Looking for SRR file')
    time.sleep(0.5)
    d.raw = f'{d.base}/raw_seqs'
    make_directory(d.raw)
    if not os.path.isfile(f'{d.base}/SRR'):
        logger('WARNING: No SRR file found.')
        p.files = "NoSRR"                                                       # This will allow for later portion of script to determine if no SRR AND no fastqs present
    else:
        logger('SRR file found.')
        make_directory(f'{d.base}/tmp-dir')
        SRRFile = open(f'{d.base}/SRR', 'r')
        logger('Checking for SRR sequence files and downloading those absent.')
        p.files = "SRRpresent"
        time.sleep(0.5)

        for line in SRRFile:
            if len(line.strip()) > 0:                                           # handles blanklines in SRR file
                ID = line.strip()
                if glob.glob(f'{d.raw}/{ID}_1.fastq.gz') and glob.glob(f'{d.raw}/{ID}_2.fastq.gz'):
                    logger(f'{ID} Zipped L and R files exist in raw_seqs.')
                elif glob.glob(f'{d.base}/{ID}_1.fastq.gz') and glob.glob(f'{d.base}/{ID}_2.fastq.gz'):
                    logger(f'{ID} Zipped L and R files exist in raw_seqs.')
                elif glob.glob(f'{d.base}/{ID}_1.fastq') and glob.glob(f'{d.base}/{ID}_2.fastq'):
                    logger(f'{ID} Unzipped L and R files exist in base directory.')
                else:
                    # FIRST, PREFETCHING
                    logger(f'Prefetching {ID}...')
                    cmd = f'singularity exec {p.singularity_loc}/{p.sratoolkit_sif} prefetch -O {d.base} {ID}'
                    run_command_logger(cmd)
                    confirm_complete(f'{d.base}/{ID}.sra')

                    # NEXT, PULLING DATASET AND SPLITTING INTO L AND R READS
                    logger(f'Pulling {ID} file from NCBI. \nThis may take a few minutes...')
                    cmd = f'singularity exec {p.singularity_loc}/{p.sratoolkit_sif} fasterq-dump \
                        {ID} --skip-technical --split-files -t {d.base}/tmp-dir -e {p.cores} -p'
                    run_command_logger(cmd)
                    try:
                        confirm_complete(f'{d.base}/{ID}_1.fastq')
                        confirm_complete(f'{d.base}/{ID}_2.fastq')
                        logger(f'{ID} L and R files successfully downloaded.')
                        remove(f'{d.base}/{ID}.sra')
                    except:
                        logger(f'Issue with {ID} SRA files. Confirm IDs are correct.')
            else:
                pass
        remove(f'{d.base}/tmp-dir')

def move_fastqs(d, p):
    """Moves sequence files (.fastq, .fq, .fastq.gz, .fq.gz) to raw_seqs."""

    logger('\nChecking for local sequence files.')
    time.sleep(0.5)
    file_types = ('*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz')                   # Looks for at least one of these file types in CWD
    files_present = []
    for type in file_types:
        files_present.extend(glob.glob(os.path.join(d.base, type)))
    if len(files_present) > 0:
        logger('Found sequence files.')
        make_directory(d.raw)

        for file in os.listdir(d.base):                                         # If files are fq or fq.gz, rename to fastq or fastq.gz respectively
            if file.endswith('fq'):
                newfile = f"{file[:2]}fastq"                                    # converting fq to fastq
                shutil.move(file, newfile)
                confirm_complete(newfile)
                logger(f'renamed {file} to {newfile} for consistency')
            elif file.endswith('fq.gz'):
                newfile = f"{file[:-5]}fastq.gz"                                # converting fq.gz to fastq.gz
                shutil.move(file, newfile)
                confirm_complete(newfile)
                logger(f'renamed {file} to {newfile} for consistency')
            else:
                pass

        for file in os.listdir(d.base):
            if file.endswith('.fastq'):                                         # If files are not zipped, zip and move to raw_seqs
                try:
                    logger(f'Zipping {file}')
                    if which('pigz') is not None:
                        subprocess.run(f'pigz {file}', shell=True)              # If 'pigz' DNE on current system, use gzip- (pigz is faster)
                        confirm_complete(f'{d.base}/{file}.gz')
                    else:
                        subprocess.run(f'gzip {file}', shell=True)
                        confirm_complete(f'{d.base}/{file}.gz')
                    logger(f'Moving {file}.gz to {d.base}/raw_seqs')
                    shutil.move(f'{d.base}/{file}.gz', f'{d.raw}/{file}.gz')
                    confirm_complete(f'{d.raw}/{file}.gz')
                    time.sleep(0.5)
                except BaseException:
                    logger(f'\nERROR: Issue with {file}. Quitting this analysis.')
                    exit()

            elif file.endswith('.fastq.gz'):                                    # If files are already zipped, move to raw_seqs
                try:
                    logger(f'{file} zipped previously, moving to {d.base}/raw_seqs')
                    shutil.move(f'{d.base}/{file}', f'{d.raw}/{file}')
                    confirm_complete(f'{d.raw}/{file}')
                    time.sleep(0.5)
                except BaseException:
                    logger(f'\nERROR: Issue moving {file} to {d.raw}. Quitting this analysis.')
                    exit()
            else:
                pass
    else:
        logger(f'No additional fastq files found.')
        time.sleep(0.5)
        if p.files == "NoSRR":
            logger('\n***ERROR: You do not have any local fastq files, nor do you have a SRR file with IDs for download. \nNo data for this analysis. Quitting.')
            quit()

def trim_reads(d, p):
    """Trims L and R read files and trims using SeqyClean."""

    # TRIM RAW DATA- BOTH ADAPTER AND QUALITY TRIMMING
    logger('TRIMMING READ DATA\n')
    make_directory(d.trim)
    for file in os.listdir(d.raw):
        if file.endswith('_R1_001.fastq.gz'):
            base = file[:-16]                                                   # removing "_R1_001.fastq.gz"
            if base.endswith('_L001') or base.endswith('_L'):                   # Illumina platform includes this LOO1 check in *most* of it's runs, but not sure if it does for
                base = base[:-5]                                                # all, so I am adding this switch for naming "base"
            else:
                pass
            isolate_r1 = file
            isolate_r2 = file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
        elif file.endswith('_1.fastq.gz'):
            base = file[:-11]                                                   # removing '_1.fastq.gz' suffix
            isolate_r1 = file
            isolate_r2 = f'{base}_2.fastq.gz'
        else:
            continue

        if confirm_present(f'{d.trim}/{base}_PE1.fastq.gz') and confirm_present(f'{d.trim}/{base}_PE2.fastq.gz'):
            logger(f'{base} has been trimmed previously. Moving to next sample.')
        else:
            logger(f'\nTrimming {base}')
            flags = ''
            if p.seqyclean_qualTrim.upper().strip() == 'Y':                     # switch to allow USER to turn off quality trimming
                flags = f'-minlen {p.seqyclean_minlen} -qual {p.seqyclean_qual} -c {p.seqyclean_adapters} {p.seqyclean_add}'
            elif p.seqyclean_qualTrim.upper().strip() == 'N':
                flags = f'-minlen {p.seqyclean_minlen} -c {p.seqyclean_adapters} {p.seqyclean_add}'
            else:
                logger('Unsure if you want to execute quality trimming or not with SeqyClean. Please input "Y" or "N" for seqyclean_qualTrim in config file.\nQuitting.')
                exit()

            cmd = f'singularity exec {p.singularity_loc}/{p.seqyclean_sif} seqyclean {flags} \
                -1 {d.raw}/{isolate_r1} -2 {d.raw}/{isolate_r2} -o {d.trim}/{base}'
            run_command_logger(cmd)
            confirm_complete(f'{d.trim}/{base}_PE1.fastq')
            confirm_complete(f'{d.trim}/{base}_PE2.fastq')
            logger(f'Zipping trimmed files for {base}')
            if which('pigz') is not None:
                subprocess.run(f'pigz {d.trim}/{base}*', shell=True)            # If 'pigz' DNE on current system, use gzip- (pigz is faster)
            else:
                subprocess.run(f'gzip {d.trim}/{base}*', shell=True)
            confirm_complete(f'{d.trim}/{base}_PE1.fastq.gz')
            confirm_complete(f'{d.trim}/{base}_PE2.fastq.gz')
    logger('\nAll reads have been trimmed and zipped. \nTrimming stage complete.')

def shuffle_reads(d, p):
    """Shuffles L and R reads into single, interleaved file"""

    # INTERLEAVE L AND R READFILES
    logger('INTERLEAVING TRIMMED READS\n')
    make_directory(d.shuffle)
    for sample in d.sample_list:
        if confirm_present(f'{d.shuffle}/{sample}.shuffled.fastq.gz'):
            logger(f'{sample}.shuffled.fastq.gz already exists. Moving to next sample.')
        else:
            logger(f'\nRunning shufflereads.pl from Lyve-SET container on {sample}')
            cmd = f'singularity exec {p.singularity_loc}/{p.lyveset_sif} run_assembly_shuffleReads.pl \
            {d.trim}/{sample}_PE1.fastq.gz {d.trim}/{sample}_PE2.fastq.gz > {d.shuffle}/{sample}.shuffled.fastq'
            run_command_logger(cmd)

            logger(f'Zipping interleaved file for {sample}')
            if which('pigz') is not None:
                subprocess.run(f'pigz {d.shuffle}/{sample}.shuffled.fastq', shell=True)              # If 'pigz' DNE on current system, use gzip- (pigz is faster)
            else:
                subprocess.run(f'gzip {d.shuffle}/{sample}.shuffled.fastq', shell=True)
            confirm_complete(f'{d.shuffle}/{sample}.shuffled.fastq.gz')
    logger('\nRead shuffling complete.')

def run_fastqc(d, p):
    """Assess quality of raw and trimmed reads."""

    # RUN FASTQC ON ALL RAW FILES
    logger('QC CHECK ON READ DATA\n')
    make_directory(d.fastqc)
    raw_files = ('*_R1_001.fastq.gz', '*_R2_001.fastq.gz', '*_1.fastq.gz', '*_2.fastq.gz')
    for f_type in raw_files:
        readfiles = glob.glob(os.path.join(d.raw, f_type))
        for file in readfiles:
            file_qc = os.path.basename(file).replace('.fastq.gz', '_fastqc.html')
            if confirm_present(f'{d.fastqc}/{file_qc}'):
                logger(f'{file_qc} exists. Moving to next readset.')
            else:
                logger(f'\nRunning FASTQC on {file}')
                cmd = f'singularity exec {p.singularity_loc}/{p.fastqc_sif} fastqc {file} -t {p.threads} {p.fastqc_add} -o {d.fastqc}'
                run_command_logger(cmd)
                confirm_complete(f'{d.fastqc}/{file_qc}')
    # RUN FASTQC ON ALL TRIMMED FILES
    trim_files = ('*_PE1.fastq.gz', '*_PE2.fastq.gz')
    for f_type in trim_files:
        readfiles = glob.glob(os.path.join(d.trim, f_type))
        for file in readfiles:
            file_qc = os.path.basename(file).replace('.fastq.gz', '_fastqc.html')
            if confirm_present(f'{d.fastqc}/{file_qc}'):
                logger(f'{file_qc} exists. Moving to next readset.')
            else:
                logger(f'\nRunning FASTQC on {file}')
                cmd = f'singularity exec {p.singularity_loc}/{p.fastqc_sif} fastqc {file} -t {p.threads} {p.fastqc_add} -o {d.fastqc}'
                run_command_logger(cmd)
                confirm_complete(f'{d.fastqc}/{file_qc}')
    logger('\nRead QC complete.')

def get_sample_list(dir):
    """Generates a list of samples in provided directory. Will remove 'R1_001.fastq.gz' or '_1.fastq.gz' extensions."""

    # GETS LIST OF SAMPLES FOR FUTURE ITERATIVE FUNCTIONS
    sample_list = []
    for file in os.listdir(dir):
        if file.endswith('_PE1.fastq.gz'):
            base = file[:-13]                                                   # removing "_PE1.fastq.gz" suffix
            sample_list.append(base)
    sample_list.sort()
    return sample_list

def run_metaphlan2(d, p):
    """Runs MetaPhlAn2 to assess metagenomic construction of sample. Here, we use it as a contamination check."""

    logger('RUNNING METAPHLAN2 TO CHECK SAMPLES FOR CONTAMINATION\n')
    make_directory(d.metaphlan2)

    # RUN METAPHLAN@ ON ALL SAMPLES
    for sample in d.sample_list:
        if confirm_present(f'{d.metaphlan2}/{sample}_profile.txt'):
            logger(f"{sample} has already been Phlan'd. Moving to next sample.")
        else:
            logger(f'\nRunning MetaPhlAn2 to assess contamination of {sample}.')
            cmd = f'singularity exec -B {p.metaphlan2_loc}:/opt/conda/bin/metaphlan_databases \
            {p.singularity_loc}/{p.metaphlan2_sif} metaphlan2.py {d.shuffle}/{sample}.shuffled.fastq.gz\
            --bowtie2out {d.metaphlan2}/{sample}.bowtie2out.txt --input_type fastq --nproc {p.cores} {p.metaphlan_add}> {d.metaphlan2}/{sample}_profile.txt'
            run_command_logger(cmd)
            confirm_complete(f'{d.metaphlan2}/{sample}_profile.txt')

    # MERGING METAPHLAN2 OUTPUT TO A COMPILED TABLE
    logger('Removing any existing MetaPhlAn summary files and generating with all samples in this run.')
    remove(f'{d.metaphlan2}/merged_metaphlan_table.txt')
    remove(f'{d.metrics}/metagenomic_composition.tsv')

    file_list = glob.glob(os.path.join(d.metaphlan2, '*_profile.txt'))          # simply passing '*_profile.txt' to metaphlan via singularity results in
    file_list.sort()                                                            # samples being out of order.
    file_str = ' '.join(file_list)

    logger('\nMerging all samples into MetaPhlAn abundance table for parsing.')
    cmd = f'singularity exec -B {p.metaphlan2_loc}:/opt/conda/bin/metaphlan_databases \
    {p.singularity_loc}/{p.metaphlan2_sif} merge_metaphlan_tables.py {file_str} >\
    {d.metaphlan2}/merged_metaphlan_table.txt'
    run_command_logger(cmd)
    confirm_complete(f'{d.metaphlan2}/merged_metaphlan_table.txt')

    # PARSING COMPILED TABLE TO GENERATE TSV WITH Genus/Species AND ABUNDANCE
    meta_f = open(f'{d.metaphlan2}/merged_metaphlan_table.txt', 'r')
    logger('Generating new file containing metagenomic composition of all samples in this run.')
    meta_summary_format = open(f'{d.metrics}/metagenomic_composition.tsv', 'w')
    abundances = {}
    for line in meta_f:
        if line.startswith('#'):
            database = line.replace('#', '')
        elif line.startswith('clade_name'):
            samp_headers = line.strip().split('\t')[2:]                         # removing the NCBI_tax_id column as it is not needed here for the summary
            new_ID_header = ['Organism']
            for ID in samp_headers:
                new_ID = ID.replace('_profile', '')
                new_ID_header.append(new_ID)
            meta_summary_format.write('\t'.join(new_ID_header) + '\n')
        else:
            outline = []
            columns = line.strip().split('\t')                                  # each identifier is prefaced by a clade type identifier. ex. 'k__Bacteria' would be Kingdom:Bacteria
            clades = columns[0].split('|')
            if len(clades) == 7:                                                # Only grabbing those entries that contain species
                id = clades[6]
                species = id.split('__')[1]
                kingdom = clades[0].split('__')[1]                              # separating by kingdom - Bacteria, virus, animalia, etc.
                kingdom = kingdom.upper()
                if kingdom not in abundances.keys():
                    abundances[kingdom] = [[species, *columns[2:]]]             # Need to classify this as a list of lists, otherwise appending to it will append to the first list of abundances entered
                else:
                    abundances[kingdom].append([species, *columns[2:]])

    t_line = []
    for kingdom in abundances.keys():
        meta_summary_format.write(f'#{kingdom} IDENTIFIED\n')
        for species_line in abundances[kingdom]:
            meta_summary_format.write('\t'.join(species_line) + '\n')
            t_line.append(species_line)
    meta_summary_format.write(f'\n#DATABASE_USED:{database}')
    meta_f.close()
    meta_summary_format.close()
    confirm_complete(f'{d.metrics}/metagenomic_composition.tsv')

    # METAPHLAN OUTPUT HAS GENUS/SPECES AS ROW IDS AND SAMPLES AS COL HEADERS. TRANSPOSING FOR MULTIQC
    for index, header in enumerate(new_ID_header):
        if header == 'Organism':
            new_ID_header[index] = 'Sample'
    transposed = zip(*t_line)
    transposed_l = []
    for line in transposed:
        transposed_l.append(line)
    meta_summary_format_t = open(f'{d.metaphlan2}/metagenomic_composition_multiqc.tsv', 'w')
    outlines_sorted_dic = {}                                                     # after transposing, the TSV is no longer sorted like all the other outputs. will sort keys of dict to achieve this
    for index, sample in enumerate(new_ID_header):
        if sample == 'Sample':
            header = '\t'.join(transposed_l[index])
            firstline = f'Sample\t{header}\n'
            meta_summary_format_t.write(firstline)
        else:
            transposed_line = '\t'.join(transposed_l[index])
            outlines_sorted_dic[sample] = transposed_line

    for sample in sorted(outlines_sorted_dic.keys()):                           # Sorting dictionary keys, which are all sample names
        line = outlines_sorted_dic[sample]
        outline = f'{sample}\t{line}\n'
        meta_summary_format_t.write(outline)
    meta_summary_format_t.close()
    confirm_complete(f'{d.metaphlan2}/metagenomic_composition_multiqc.tsv')
    logger('\nMetagenomic contamination assessment with MetaPhlAn2 is complete.')

def run_shovill(d, p):
    """Assembles reads using Shovill"""

    # RUN SHOVILL ON ALL SAMPLES
    logger('GENOME ASSEMBLY WITH SHOVILL\n')
    make_directory(d.shovill)
    for sample in d.sample_list:
        if confirm_present(f'{d.shovill}/{sample}/{sample}.fa'):
            logger(f'{sample} assembly exists. Moving to next sample.')
        else:
            isolate_r1 = f'{sample}_PE1.fastq.gz'
            isolate_r2 = f'{sample}_PE2.fastq.gz'
            if confirm_present(f'{d.trim}/{isolate_r1}') and confirm_present(f'{d.trim}/{isolate_r2}'):
                logger(f'\nNow assembling genome for {sample} using Shovill.')
                flags = f'--cpus {p.cores} --ram {p.ram} --minlen {p.shovill_minlen} --mincov {p.shovill_mincov} {p.shovill_add}'
                cmd = f'singularity exec {p.singularity_loc}/{p.shovill_sif} shovill {flags} \
                    --outdir {d.shovill}/{sample} --R1 {d.trim}/{isolate_r1} --R2 {d.trim}/{isolate_r2}'
                run_command_logger(cmd)
                confirm_complete(f'{d.shovill}/{sample}/contigs.fa')
                shutil.move(f'{d.shovill}/{sample}/contigs.fa', f'{d.shovill}/{sample}/{sample}.fa')
                confirm_complete(f'{d.shovill}/{sample}/{sample}.fa')
    logger('\nGenome assembly complete.')

def run_quast(d, p):
    """Runs Quast on all shovill assemblies to check quality"""

    # RUN QUAST ON ALL SAMPLES
    logger('QUALITY ASSESSMENT OF GENOMES WITH QUAST\n')
    make_directory(d.quast)
    for sample in d.sample_list:
        samp_str = []
        if confirm_present(f'{d.quast}/{sample}/report.tsv'):
            logger(f'QUAST has already been run on {sample} assembly. Skipping to next sample.')
        else:
            logger(f'\nRunning QUAST on {sample} assembly.')
            flags = []
            if p.quast_nocheck.upper().strip() == 'Y':
                flags.append('--no-check')
            if p.quast_noplots.upper().strip() == 'Y':
                flags.append('--no-plots')
            if p.quast_nohtml.upper().strip() == 'Y':
                flags.append('--no-html')
            if p.quast_noicarus.upper().strip() == 'Y':
                flags.append('--no-icarus')
            if p.quast_nosnps.upper().strip() == 'Y':
                flags.append('--no-snps')
            if p.quast_nogc.upper().strip() == 'Y':
                flags.append('--no-gc')
            if p.quast_nosv.upper().strip() == 'Y':
                flags.append('--no-sv')
            if p.quast_nogzip.upper().strip() == 'Y':
                flags.append('--no-gzip')
            if p.quast_noreadstats.upper().strip() == 'Y':
                flags.append('--no-read-stats')

            cmd = f'singularity exec {p.singularity_loc}/{p.quast_sif} quast.py \
                {" ".join(flags)} -t {p.threads} {p.quast_add} \
                {d.shovill}/{sample}/{sample}.fa -o {d.quast}/{sample}'
            run_command_logger(cmd)                                             # if you dont want GC content, just run this with --fast flag instead of all the `--no` flags
            confirm_complete(f'{d.quast}/{sample}/report.tsv')
    logger(f'\nQuality assessment of assemblies with QUAST is complete on all samples.')

    logger('Removing any existing QUAST summary file and generating new with all samples in this run.')
    remove(f'{d.metrics}/quast_metrics.tsv')

    logger('Generating compiled Quast metric file.')
    quast_metrics = open(f'{d.metrics}/quast_metrics.tsv', 'w')
    quast_list = ['Sample\tContigs\tLargest_Contig\tTotal_Lentgh\tGC\tN50']

    # PARSE EACH OUTPUT FILE FOR DESIRED METRICS
    for sample in d.sample_list:
        sample_metrics = open(f'{d.quast}/{sample}/report.tsv')
        for line in sample_metrics:
            cols = line.split('\t')[0]
            value = line.rstrip().split('\t')[-1]                               # acutal metric is always the last element of the tab-separated line
            if cols == '# contigs':
                contigs = value
            elif cols == 'Largest contig':
                lg_ctg = value
            elif cols == 'Total length (>= 0 bp)':
                tot_length = value
            elif cols == 'GC (%)':
                gc_content = value
            elif cols == 'N50':
                n50 = value
        samp_str = [sample, contigs, lg_ctg, tot_length, gc_content, n50]
        quast_list.append('\t'.join(samp_str))
    for line in quast_list:
        quast_metrics.write(f'{line}\n')
    quast_metrics.close()
    confirm_complete(f'{d.metrics}/quast_metrics.tsv')
    logger('\nAssembly QC complete.')

def run_bwa(d, p):
    """Maps trimmed reads to shovil assembly"""

    # RUN BWA ON ALL SAMPLES
    logger('READ ALIGNMENT WITH BWA\n')
    make_directory(d.bwa)
    make_directory(f'{d.bwa}/temp_files')
    for sample in d.sample_list:
        if confirm_present(f'{d.bwa}/{sample}.alignment.sorted.bam'):
            logger(f'{sample} reads have already been aligned to assembly. Skipping to next sample.')
        else:
            logger(f'\nIndexing {sample} assembly.')                             # First create index of assembly
            cmd = f'singularity exec {p.singularity_loc}/{p.bwa_sif} \
                bwa index {p.bwa_index_add} {d.shovill}/{sample}/{sample}.fa'
            run_command_logger(cmd)

            isolate_r1 = f'{sample}_PE1.fastq.gz'                              # Next aligning reads
            isolate_r2 = f'{sample}_PE2.fastq.gz'
            cmd = f'singularity exec {p.singularity_loc}/{p.bwa_sif} \
                bwa mem -t {p.threads} {p.bwa_mem_add} {d.shovill}/{sample}/{sample}.fa \
                {d.trim}/{isolate_r1} {d.trim}/{isolate_r2} > {d.bwa}/{sample}.alignment.sam'
            logger(f'Aligning {sample} reads to assembly.')
            run_command_logger(cmd)
            confirm_complete(f'{d.bwa}/{sample}.alignment.sam')
                                                                                # Finally, sorting SAM file and converting to BAM
            cmd = f'singularity exec {p.singularity_loc}/{p.samtools_sif} samtools sort \
                {d.bwa}/{sample}.alignment.sam --threads {p.threads} -T {d.bwa}/temp_files {p.samtools_add} \
                -O BAM -o {d.bwa}/{sample}.alignment.sorted.bam'
            logger(f'Sorting SAM file for {sample}.')
            run_command_logger(cmd)
            confirm_complete(f'{d.bwa}/{sample}.alignment.sorted.bam')

    for sample in d.sample_list:
        if confirm_present(f'{d.bwa}/{sample}.stats'):
            logger(f'Summary file already generated for {sample} BAM.')
        else:
            logger(f'Generating summary file for {sample} BAM')

            # GENERATING STATS FILE THAT I WILL PULL SUMMARY DATA FROM LATER
            cmd = f'singularity exec {p.singularity_loc}/{p.samtools_sif} samtools \
            stats {d.bwa}/{sample}.alignment.sorted.bam > {d.bwa}/{sample}.stats'
            run_command_logger(cmd)
            confirm_complete(f'{d.bwa}/{sample}.stats')

    remove(f'{d.bwa}/temp_files')
    logger('\nBWA alignment complete.')

def run_mash(d, p):
    """Metagenomic distance estimation using Min-Hash method"""

    logger('METAGENOMIC DISTANCE ESTIMATION WITH MASH\n')
    make_directory(d.mash)

    # RUN MASH ON ALL SAMPLES
    for sample in d.sample_list:
        if confirm_present(f'{d.mash}/{sample}.distance.tab'):
            logger(f'{sample} has already been MASHed. Skipping to next sample.')
        else:
            if not os.path.isfile(f'{d.mash}/{sample}.distance.tab'):  # This switch isnt SUPER necessary, but allows for re-start if only sketches completed before a failure for some reason
                logger(f"\nSketching profile of the lovely {sample}. Rembrandt ain't got nothin...")
                cmd = f'singularity exec {p.singularity_loc}/{p.mash_sif} mash sketch \
                    -m {p.mash_sketch_m} {p.mash_sketch_add} {d.shuffle}/{sample}.shuffled.fastq.gz'
                run_command_logger(cmd)
                    # QUESTION: -m 2 excludes single-copy kmers. If remove this, replace with -r to indicate input is reads.
                shutil.move(f'{d.shuffle}/{sample}.shuffled.fastq.gz.msh', f'{d.mash}/{sample}.shuffled.fastq.gz.msh')
                confirm_complete(f'{d.mash}/{sample}.shuffled.fastq.gz.msh')

                logger(f'Calculating distance for {sample}.')                   # decided not to retain the STDOUT in memory, as the files are quite large (i.e. could have piped to stdout and written at end rather than writing `tab` files)
                cmd = f'singularity exec {p.singularity_loc}/{p.mash_sif} mash dist \
                    {p.mash_dist_db} {p.mash_dist_add} {d.mash}/{sample}.shuffled.fastq.gz.msh > {d.mash}/{sample}.distance.tab'
                run_command_logger(cmd)
                confirm_complete(f'{d.mash}/{sample}.distance.tab')

    # GENERATE FILE CONTAINING TOP TEN HITS FOR EACH SAMPLE
    logger('Removing any existing MASH summary files and generating new with all samples in this run')
    remove(f'{d.metrics}/top_hits_MASH.tsv')
    remove(f'{d.mash}/sample_organismIDs.tsv')

    logger('\nGenerating top hits of MASH. Hopefully "Monster Mash" will be included... what a classic.')
    mash_samp_ID = open(f'{d.mash}/sample_organismIDs.tsv', 'w')
    mash_samp_ID.write(f'SAMPLE\tORGANISM\tDISTANCE\tPROBABILITY\n')
    mash_top_hits = open(f'{d.metrics}/top_hits_MASH.tsv', 'w')
    mash_hits_l = []
    for sample in d.sample_list:
        mash_hits = []
        mash_sample = open(f'{d.mash}/{sample}.distance.tab', 'r')              # Sort file by 3rd column and hold the top 10 in memory for writing to "top hits" file
        for line in mash_sample:
            # pval = float(line.split('\t')[3])
            # if pval < 0.01:                                                   # only pulls those with a pvalue of 0. I've commented this out, as we'll want to see what else is in there, even if only for edification purposes
            mash_hits.append(line)
        mash_hits.sort(key = lambda x: (x.split('\t')[3], x.split('\t')[2]))    # [reference-ID, query-ID, distance, p-value, shared-hashes]. I will sort by p-value and then by distance

        mash_hits_l.append(f'{sample}\n')                                       # mash_hits_l will write out the "top 10" hits for each sample into single file
        for hit in mash_hits[:10]:
            mash_hits_l.append(hit)
        mash_hits_l.append('\n')
        mash_sample.close()

        top_hit_cols = mash_hits[0].split('\t')                                 # generating 'sample_organismIDs.tsv' file containing top hit for each organism.
        sample = top_hit_cols[1]
        sample_name = os.path.basename(sample)[:-18]                            # removing '.shuffled.fastq.gz' suffix
        organism = top_hit_cols[0]
        organism_name = re.sub('.*-', '', organism)[:-4]                        # isolate the genus-species name only
        distance = top_hit_cols[2]
        probability = top_hit_cols[3]
        mash_samp_ID.write(f'{sample_name}\t{organism_name}\t{distance}\t{probability}\n')
    mash_samp_ID.close()

    for entry in mash_hits_l:                                               # writing top 10 hits for each sample to file
        mash_top_hits.write(entry)
    mash_top_hits.close()
    logger('\nMish MASH complete.')
    confirm_complete(f'{d.mash}/sample_organismIDs.tsv')
    confirm_complete(f'{d.metrics}/top_hits_MASH.tsv')

def organism_switch(d, p):
    """Determines what MASH identified as the "top" ID, and then runs Serotypefinder if E. coli, or seqsero2 and SISTR if Salmonella"""

    logger('SEROTYPE PREDICTION BASED ON MASH OUTPUT\n')
    if not os.path.isfile(f'{d.mash}/sample_organismIDs.tsv'):
        logger('Cannot find MASH ID file to direct running of Serotypefinder vs seqsero2 & SISTR. Check MASH completed successfully.')
    else:
        id_file = open(f'{d.mash}/sample_organismIDs.tsv', 'r')
        e_coli_l = []
        sal_l = []
        streppy_l = []
        streppn_l = []
        leg_l = []
        for line in id_file:
            if line.startswith('SAMPLE'):
                pass
            else:
                samp_name = line.split('\t')[0]
                organism = line.split('\t')[1]
                if "Escherichia_coli" in organism:
                    logger(f'{samp_name} has organsim ID {organism} which contains "E coli"')
                    e_coli_l.append(samp_name)
                elif "Salmonella" in organism:
                    logger(f'{samp_name} has organsim ID {organism} which contains "Salmonella"')
                    sal_l.append(samp_name)
                elif "Streptococcus_pyogenes" in organism:
                    logger(f'{samp_name} has organsim ID {organism} which contains "Streptococcus pyogenes"')
                    streppy_l.append(samp_name)
                elif "Streptococcus_pneumoniae" in organism:
                    logger(f'{samp_name} has organsim ID {organism} which contains "Streptococcus pneumoniae"')
                    streppn_l.append(samp_name)
                elif "Legionella_pneumophila" in organism:
                    logger(f'{samp_name} has organsim ID {organism} which contains "Legionella pneumophila"')
                    leg_l.append(samp_name)
                else:
                    logger(f'{samp_name} has organism ID {organism} which does not contain E coli, Salmonella, Streptococcus pyogenes, Streptococcus pneumoniae, nor Legionella pneumophila')
                    pass
        id_file.close()

    if len(e_coli_l) > 0:
        run_serotypefinder(e_coli_l, d, p)
    if len(sal_l) > 0:
        run_seqsero2(sal_l, d, p)
    if len(streppy_l) > 0:
        run_emmTyping(streppy_l, d, p)
    if len(streppn_l) > 0:
        run_seroba(streppn_l, d, p)
    if len(leg_l) > 0:
        run_legsta(leg_l, d, p)
    #     run_lpserogroup(leg_l, d, p)
    logger('\nSerotyping of all samples complete')

def run_serotypefinder(samples, d, p):
    """Runs serotypefinder on all isolates provided (should only be E. coli samples)."""

    logger('')
    logger(f'Samples {", ".join(samples)} have been identified as E. coli. Will now run Serotypefinder on these samples.')
    make_directory(d.serotypefinder)

    # RUN SEROTYPEFINDER ON ALL SAMPLES
    for sample in samples:
        if confirm_present(f'{d.serotypefinder}/{sample}/results_table.txt'):
            logger(f'Serotypefinder has already been run on {sample}. Moving to next sample.')
        else:
            logger(f'\nRunning SerotypeFinder on {sample}.')
            cmd = f'singularity exec {p.singularity_loc}/{p.serotypefinder_sif} serotypefinder.pl \
                -d {p.serotypefinder_db} -i {d.shovill}/{sample}/{sample}.fa  \
                -b {p.serotypefinder_blast} -s {p.serotypefinder_species} -k {p.serotypefinder_threshold} \
                -l {p.serotypefinder_minleng} -o {d.serotypefinder}/{sample} {p.serotypefinder_add}'
            run_command_logger(cmd)
            confirm_complete(f'{d.serotypefinder}/{sample}/results_table.txt')

    # GENERATE FILE OF COMPILED SEROTYPEFINDER OUTPUTS
    logger('Removing any existing serotype summary files and generating new with all samples in this run')
    remove(f'{d.metrics}/e_coli_serotypes.tsv')
    ecoli_sero = open(f'{d.metrics}/e_coli_serotypes.tsv', 'w')
    ecoli_sero.write(f'SAMPLE\tO_type:H_type\n')

    for sample in samples:
        results_f = open(f'{d.serotypefinder}/{sample}/results_table.txt', 'r')
        results = results_f.readlines()
        if results[1] == 'No serotype predicted.\n':
            H_type = "No_H"                                                 # If H serotype is not succesfully predicted, line 3 is blank. This will catch that.
        else:
            H_type = results[2].split('\t')[5]                              # H serotype is on 3rd line of file, 6th column

        if results[5] == 'No serotype predicted.\n':
            O_type = "No_O"                                                 # If H serotype is not succesfully predicted, line 3 is blank. This will catch that.
        else:
            O_type = results[6].split('\t')[5]                              # H serotype is on 7th line of file, 6th column
        ecoli_sero.write(f'{sample}\t{O_type}:{H_type}\n')
        results_f.close()

    ecoli_sero.close()
    confirm_complete(f'{d.metrics}/e_coli_serotypes.tsv')
    logger('Serotyping of E. coli with serotypefinder complete.\n')

def run_seqsero2(samples, d, p):
    """Runs seqsero2 on all isolates provided (should only be Salmonella samples)."""

    logger(f'Samples {", ".join(samples)} have been identified as Salmonella. Will now run seqsero2 and SISTR on these samples.')

    make_directory(d.seqsero2)
    for sample in samples:
        if confirm_present(f'{d.seqsero2}/{sample}/Seqsero_result.txt'):
            logger(f'SeqSero2 has already been run on {sample}. Skipping to next sample.')
        else:
            logger(f'\nSerotyping {sample} with SeqSero2.')
            # Seqsero v1
#            cmd = f'singularity exec {p.singularity_loc}/{p.seqsero2_sif} seqsero.py \
#            -m4 -i {d.shovill}/{sample}/{sample}.fa -d {d.seqsero2}/{sample}'

    # THIS IS FOR seqsero2 as opposed to V1
            cmd = f'singularity exec {p.singularity_loc}/{p.seqsero2_sif} SeqSero2_package.py \
            -t4 -m {p.seqsero2_mode} -i {d.shovill}/{sample}/{sample}.fa -p {p.threads} \
            -b {p.seqsero2_mapping} -d {d.seqsero2}/{sample}'
            run_command_logger(cmd)
            confirm_complete(f'{d.seqsero2}/{sample}/Seqsero_result.txt')
    # GENERATE FILE OF COMPILED seqsero2 OUTPUTS
    logger('Removing any existing seqsero2 summary files and generating new with all samples in this run')
    remove(f'{d.metrics}/seqsero2_serotypes.tsv')

    seqsero2_types = open(f'{d.metrics}/seqsero2_serotypes.tsv', 'w')
    for sample in samples:
        seqsero2_types.write(f'{sample}\n')
        sample_type = open(f'{d.seqsero2}/{sample}/Seqsero_result.txt', 'r')
        for line in sample_type:
            seqsero2_types.write(line)
        sample_type.close()
    seqsero2_types.close()
    confirm_complete(f'{d.metrics}/seqsero2_serotypes.tsv')
    logger('Serotyping of Salmonella with seqsero2 complete.')

def run_emmTyping(samples, d, p):
    """Runs the emm_typing.py script developed by the CDC to serotype Streptococcus pyogenes samples"""

    logger(f'Samples {", ".join(samples)} have been identified as Streptococcus. Will now run emm_typing on these samples.')

    make_directory(d.emm_sero)
    for sample in samples:
        if os.path.isfile(f'{d.emm_sero}/emm_{sample}/ComponentComplete.txt'):
            logger(f'EMM typing has already been run on {sample}. Skipping to next sample.')
        else:
            logger(f'\nSerotyping {sample} with EMM_typing.')

            cmd = f'singularity exec -e -B {d.base}:/EMBOSS-6.6.0/emboss/.libs/ {p.singularity_loc}/{p.emmtyping_sif} emm_typing.py \
            --profile_file_directory {p.emm_db} \
            --fastq_1 {d.trim}/{sample}_PE1.fastq.gz \
            --fastq_2 {d.trim}/{sample}_PE2.fastq.gz \
            --output_directory {d.emm_sero}/emm_{sample} {p.emm_add}'
            run_command_logger(cmd)
            remove(f'{d.base}/lt-seqret')

            if os.path.isfile(f'{d.emm_sero}/emm_{sample}/ComponentComplete.txt'):  # Cant use confirm_compete here as this is an empty "completion" file. I.e., run is completed, but file is blank
                pass
            else:
                raise BaseException

    logger('Removing any existing emm summary files and regenerating with all samples in this run')
    remove(f'{d.metrics}/emm_serotypes.tsv')
    emm_types = open(f'{d.metrics}/emm_serotypes.tsv', 'w')
    emm_type_l = [('Sample', 'EMM')]
    for sample in samples:
        samp_emm_f = open(f'{d.emm_sero}/emm_{sample}/{sample}_PE1.results.xml', 'r')
        emm = 'NA'
        emm_nonval = 'NA'
        for line in samp_emm_f:
            if 'Final_EMM_type' in line:
                value = line.strip().split('value="')[1]
                emm = value.replace('value="', '').replace('.sds"/>', '')
                emm_type_l.append((sample, emm))
     # BELOW CODE WILL ALLOW FOR INCORPORATION OF NONVALIDATED EMM TYPES, AS WELL.
     # TO CHANGE THIS, COMMENT OUT ABOVE LINE AND UNCOMMENT FOLLOWING 5 LINES
        #     if 'EMM_nonValidated' in line:
        #         value = line.strip().split('value="')[1]
        #         emm_nonval = value.replace('value="', '').replace('.sds"/>', '')
        # emms = f'{emm}:{emm_nonval}'
        # emm_type_l.append((sample, emms))

    for entry in emm_type_l:
        emm_line = '\t'.join(entry)
        emm_types.write(f'{emm_line}\n')
    emm_types.close()
    confirm_complete(f'{d.metrics}/emm_serotypes.tsv')
    logger('EMM typing of Streptococcus with emm_typing.py complete.\n')

def run_seroba(samples, d, p):
    """Runs SEROBA for Strep pneumoniae. More details found here: https://sanger-pathogens.github.io/seroba/"""

    make_directory(d.seroba)
    os.chdir(d.seroba)
    for sample in samples:
        if confirm_present(f'{d.seroba}/{sample}/pred.tsv'):
            logger(f'Seroba has already been run on {sample}. Skipping to next S. pneumoniae sample.')
        else:
            cmd = f'singularity exec {p.singularity_loc}/{p.seroba_sif} seroba runSerotyping \
            {p.seroba_db} {d.trim}/{sample}_PE1.fastq.gz {d.trim}/{sample}_PE2.fastq.gz \
            {sample} --coverage {p.seroba_cvg} {p.seroba_add}'
            run_command_logger(cmd)
            confirm_complete(f'{d.seroba}/{sample}/pred.tsv')
    os.chdir(d.base)

    remove(f'{d.metrics}/seroba_serotypes.tsv')
    seroba_seros = open(f'{d.metrics}/seroba_serotypes.tsv', 'w')
    header = 'Sample\tSerotype\n'
    sero_list = []
    for sample in samples:
        cur_sero_f = open(f'{d.seroba}/{sample}/pred.tsv', 'r')
        for line in cur_sero_f:
            if line.startswith(sample):
                sero_list.append(line)
        cur_sero_f.close()
    seroba_seros.write(header)
    for sero_entry in sero_list:
        seroba_seros.write(sero_entry)
    seroba_seros.close()
    confirm_complete(f'{d.metrics}/seroba_serotypes.tsv')

def run_legsta(samples, d, p):
    """Runs legsta to serotype Legionella pneumophila samples"""

    logger(f'Samples {", ".join(samples)} have been identified as Legionella pneumophila. Will now run lpserogroup_prediction on these samples.')

    make_directory(d.leg_sero)
    for sample in samples:
        if confirm_present(f'{d.leg_sero}/legsta_{sample}.tsv'):
            logger(f'Legsta already completed for {sample}. Skipping to next sample.')
        else:
            logger(f'\nSerotyping {sample} with legsta.')
            cmd = f'singularity exec {p.singularity_loc}/{p.legsta_sif} legsta {p.legsta_add} {d.shovill}/{sample}/{sample}.fa > {d.leg_sero}/legsta_{sample}.tsv'
            run_command_logger(cmd)
            confirm_complete(f'{d.leg_sero}/legsta_{sample}.tsv')

    remove(f'{d.metrics}/legsta_serotypes.tsv')
    leg_seros = open(f'{d.metrics}/legsta_serotypes.tsv', 'w')
    sero_list = []
    for sample in samples:
        cur_sero_f = open(f'{d.leg_sero}/legsta_{sample}.tsv', 'r')
        for line in cur_sero_f:
            if line.startswith('FILE'):
                header = line
            else:
                sero_list.append(line)
        cur_sero_f.close()
    leg_seros.write(header)
    for sero_entry in sero_list:
        leg_seros.write(sero_entry)
    leg_seros.close()
    confirm_complete(f'{d.metrics}/legsta_serotypes.tsv')

# CURRENTLY, lpserogroup_prediction cannot be run with Singularity
# def run_lpserogroup(samples, d, p):
#     """Runs lpserogroup_prediction pipeline developed by the CDC to serotype Legionella pneumophila samples"""
#
#     logger(f'Samples {", ".join(samples)} have been identified as Legionella pneumophila. Will now run lpserogroup_prediction on these samples.')
#
#     make_directory(d.leg_sero)
#     for sample in samples:
#         if os.path.isfile(f'{d.leg_sero}/'):
#             logger(f'Lp serogroup prediction already completed for {sample}. Skipping to next sample.')
#         else:
#             logger(f'\nSerotyping {sample} with lpserogroup_prediction.')
#
#             # FIRST GETTING THE REFERENCE AND GFF FROM WITHIN THE container
#             # using glob so they can change it within the container and not break this pipeline
#             ref_loc = '/opt/pipeline/reference'
#             cmd = f'singularity exec {p.singularity_loc}/{p.legtyping_sif} ls {ref_loc}'
#             run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
#             for line in run.stdout:
#                 if 'fna' in line:
#                     ref_f = line.strip()
#                     lp_reference = f'{ref_loc}/{ref_f}'
#                 elif 'gff' in line:
#                     gff_f = line.strip()
#                     lp_gff = f'{ref_loc}/{gff_f}'
#
#             make_directory(f'{d.leg_sero}/{sample}')
#             make_directory(f'{d.leg_sero}/{sample}/output')
#             cmd = f'singularity exec -B {d.leg_sero}/{sample}/output:/output -B {d.trim}:/data {p.singularity_loc}/{p.legtyping_sif} /opt/pipeline/pipeline.sh \
#             --reference={lp_reference} --gff={lp_gff} \
#             --r1=/data/{sample}_PE1.fastq.gz --r2=/data/{sample}_PE2.fastq.gz \
#             --isolate={sample} --output=/output'
#             run_command_logger(cmd)
#
#     logger('Serotyping of Legionella pneumophila with lpserogroup_prediction complete.\n')

def run_abricate(d, p):
    """Runs Abricate to identify Antibiotic Resistance in each sample"""

    logger('ANTIBIOITIC RESISTANCE AND VIRULENCE IDENTIFICATION WITH ABRICATE')
    make_directory(d.abricate)
    db_l = []
    cmd = f'singularity exec {p.singularity_loc}/{p.abricate_sif} abricate --list'
    run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    for line in run.stdout:                                                     # Generating a list of the DBs present in abricate container
        db_name = line.split('\t')[0]                                           # DB name is the first element of each line (though very first line is "DATABASE")
        if db_name != 'DATABASE':                                               # ISSUE: the 'abricate' db appears to not have been established...
            db_l.append(db_name)
    for DB in db_l:
        if confirm_present(f'{d.abricate}/{DB}_summary'):
            logger('Removing any existing Abricate summary files and regenerating will all samples in this run')
            remove(f'{d.metrics}/abricate_{DB}_summary')

        else:
            logger(f'\nNow querying Database {DB}')
            # RUN ALL SAMPLES THROUGH CURRENT DATABASE IN ABRICATE
            for sample in d.sample_list:
                if confirm_present(f'{d.abricate}/{sample}_{DB}.tab'):
                    logger(f'Abra-kad-abricate! {sample} has already been run. Skipping to next sample.')
                else:
                    logger(f'Database {DB} is being queried for sample {sample} to fabricate results... I mean abricate results!')
                    cmd = f'singularity exec {p.singularity_loc}/{p.abricate_sif} abricate \
                        --minid {p.abricate_minid} --mincov {p.abricate_mincov} --threads {p.threads} \
                        -db {DB} {p.abricate_add} {d.shovill}/{sample}/{sample}.fa > {d.abricate}/{sample}_{DB}.tab'
                    run_command_logger(cmd)
                    confirm_complete(f'{d.abricate}/{sample}_{DB}.tab')

            # GENERATE SUMMARY FILE OF CURRENT DATABASE RUN
            if confirm_present(f'{d.metrics}/abricate_{DB}_summary'):
                logger(f'{DB} database summary has already been created. Skipping to next.')
            else:
                logger(f'Generating summary file for query of {DB} database.')
                cmd = f'singularity exec {p.singularity_loc}/{p.abricate_sif} abricate \
                    --summary {d.abricate}/*_{DB}.tab > {d.metrics}/abricate_{DB}_summary'
                run_command_logger(cmd)
                confirm_complete(f'{d.metrics}/abricate_{DB}_summary')

    logger(f'Antibiotic resistance and virulence detection complete. Summary files in {d.metrics}')

def QualCov_metrics(d, p):
    """Generates files containing quality metrics for readsets"""

    # ASSESS METRICS WITH TRIMMED READS
    logger('READ DATA METRICS\n')
    ftype = 'trimmed'

    logger(f'Obtaining metrics using trimmed read files')
    for sample in d.sample_list:
        if confirm_present(f'{d.trim}/{sample}_readMetrics.tsv'):
            logger(f'Read Metrics already gathered for {sample}. Skipping to next sample.')
        else:
            logger(f'Parsing Read Metrics for {sample}.')
            cmd = f'singularity exec {p.singularity_loc}/{p.lyveset_sif} \
                run_assembly_readMetrics.pl {d.trim}/{sample}*PE*fastq.gz --fast -numcpus {p.cores} > {d.trim}/{sample}_readMetrics.tsv'
            run_command_logger(cmd)
            confirm_complete(f'{d.trim}/{sample}_readMetrics.tsv')

    remove(f'{d.metrics}/QC-metrics_trimmed.tsv')
    trim_metrics = open(f'{d.metrics}/QC-metrics_trimmed.tsv', 'w')             # First, get metrics using trimmed reads
    header = ''
    metrics_l = []
    for sample in d.sample_list:
        readMet_f = open(f'{d.trim}/{sample}_readMetrics.tsv', 'r')
        for line in readMet_f:                                                  # reading all results into single list for sorting
            if line.startswith('File'):
                header = line
            else:
                sample = line.split('\t')[0].split('/')[-1]
                metrics_l.append(line)
    metrics_l.sort(key = lambda x: int(x.split('\t')[2]))                       # Sorting by number of bases. Outputs as string, need to sort it as an INT to avoid 10000 being "smaller" than 900
    logger(f'Writing cleaned reads metrics to QC-metrics_trimmed.tsv')
    trim_metrics.write(header)
    for metric in metrics_l:
        trim_metrics.write(metric)
    trim_metrics.close()
    confirm_complete(f'{d.metrics}/QC-metrics_trimmed.tsv')
    logger('\nCalculation of read metrics complete.')

def compile_metrics(d, p):
    """Generates final output file containing metrics from trimming, assembly, Species OPTIMIZATION (METAPHLAN/MIDAS), and Serotype"""

    logger('COMPILING METRICS INTO SINGLE FILE')
    logger('This will take a few moments...')

    # OBTAIN METRICS FROM BAM FILES
    read_metrics = {}
    for sample in d.sample_list:
        sample_l = []
        if confirm_present(f'{d.bwa}/{sample}.stats'):
            bamf = open(f'{d.bwa}/{sample}.stats', 'r')
        reads = 0
        mapped = 0
        insert_avg = 0
        insert_stdev = 0
        for line in bamf:
            if line.startswith('SN'):                                           # Only care about the "Summary Numbers" lines
                metric = line.rstrip().split('\t')[1]
                value = line.rstrip().split('\t')[2]
                if metric == 'average length:':                                 # Generates a list of tuples containing my metrics for each sample
                    sample_l.append(('avgReadLength', value))
                elif metric == 'bases mapped (cigar):':
                    sample_l.append(('basesMapped', value))
                elif metric == 'maximum length:':
                    sample_l.append(('maxReadLength', value))
                elif metric == 'average quality:':
                    sample_l.append(('avgReadQual', value))
                elif metric == 'sequences:':
                    sample_l.append(('numReads', value))
                    reads = float(value)
                elif metric == 'reads mapped:':
                    sample_l.append(('mappedReads', value))
                    mapped = float(value)
                elif metric == 'insert size average:':
                    insert_avg = float(value)
                elif metric == 'insert size standard deviation:':
                    insert_stdev = float(value)
                elif metric == 'reads properly paired:':                        # If any pairs are detected, it is PE
                    if int(value) > 0:
                        sample_l.append(('PE', 'YES'))
                    else:
                        sample_l.append(('PE', 'NO'))
        read_metrics[sample] = sample_l
        if reads != 0 and mapped != 0:
            perc_mapped = mapped/reads * 100
            value = str(round(perc_mapped, 2))
            sample_l.append(('%ReadsMap', value))
        else:
            sample_l.append(('%ReadsMap', '-'))
        if insert_avg != 0 and insert_stdev != 0:
            insert_l = round(insert_avg - insert_stdev)
            insert_r = round(insert_avg + insert_stdev)
            insert_avg_stdev = f'{round(insert_avg)}[{insert_l},{insert_r}]'
            sample_l.append(('avgInsertSize', insert_avg_stdev))
        else:
            sample_l.append(('avgInsertSize', '-'))

    # ADD AVG GENOME COVERAGE ABOVE X (provided in config file)
    coverage_d = {}
    dots = ''
    for sample in d.sample_list:
        dots += '.'
        logger(dots)
        bases_covered = 0                                                       # Samtools mpileup will filter out some reads including anamalous pairs.
                                                                                # -A flag includes ALL reads
                                                                                # -q0 and -Q0 accept all reads and bases with this minimum PHRED
        cmd = f'singularity exec {p.singularity_loc}/{p.samtools_sif} \
        samtools mpileup -A -q0 -Q0 {d.bwa}/{sample}.alignment.sorted.bam > {d.bwa}/{sample}.mpileup'
        if confirm_present(f'{d.bwa}/{sample}.mpileup'):
            pass
        else:
            subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pileupfile = open(f'{d.bwa}/{sample}.mpileup')
        sample_totaldepth = 0
        for line in pileupfile:
            if line.startswith('[mpileup]'):                                    # Skipping first line of stdout
                pass
            else:
                coverage = float(line.split()[3])                               # fourth column of mpileup contains reads covering given site
                sample_totaldepth = sample_totaldepth + coverage                # summing coverage of each base to calculate average base coverage by genome size
                if coverage >= float(p.coverage):
                    bases_covered += 1
        coverage_d[sample] = [f'bases_cov_above_X:{bases_covered}', f'totalDepth:{sample_totaldepth}']
        pileupfile.close()

    # ADD QUAST METRICS
    quast_metrics = {}
    assembly_metrics_f = open(f'{d.metrics}/quast_metrics.tsv', 'r')
    for line in assembly_metrics_f:
        if line.startswith('Sample'):
            pass
        else:
            sample = line.strip().split('\t')[0]                                # First element of each line in the quast_metrics.tsv file is the sample name
            contigs = line.strip().split('\t')[1]
            lgst_contig = line.strip().split('\t')[2]
            genome_l = int(line.strip().split('\t')[3])
            gc = line.strip().split('\t')[4]
            N50 = line.strip().split('\t')[5]
            for entry in coverage_d[sample]:
                if entry.split(':')[0] == 'totalDepth':
                    totalDepth = float(entry.split(':')[1])
                elif entry.split(':')[0] == 'bases_cov_above_X':
                    bases_cov_above_X = float(entry.split(':')[1])
                else:
                    logger('ERROR with calculating genome coverage')
            average_coverage = totalDepth / genome_l                            # Calculating avg coverage as total depth over all bases covered, by total number of bases in genome
            average_coverage = str(round(average_coverage, 2))
            perc_cov_above_X = bases_cov_above_X / genome_l * 100               # User inputs desired coverage value. This allows to see how spread the coverage is
            perc_cov_above_X = str(round(perc_cov_above_X, 2))
            quast_metrics[sample] = [('Contigs', contigs), ('N50', N50), ('GC_content', gc), ('genomeLength', str(genome_l)), ('totalDepth', str(totalDepth)), ('avgCoverage', average_coverage), ('coverageAboveX', perc_cov_above_X)]
            for element in line.rstrip().split('\t')[1:]:
                quast_metrics[sample].append(element)
    assembly_metrics_f.close()

    # ADD ORGANISM ID FROM MASH
    sero_metrics = {}

    if confirm_present(f'{d.mash}/sample_organismIDs.tsv'):
        mash_ids = open(f'{d.mash}/sample_organismIDs.tsv', 'r')
        for line in mash_ids:
            organism = 'none'
            probability = '-'
            o_type = '-'
            h1_type = '-'
            h2_type = '-'
            antigenic = '-'
            sero = '-'
            if line.startswith('SAMPLE'):
                pass
            else:
                sample = line.split()[0]
                organism = line.split()[1].split('_subsp')[0]                   # Salmonella_enterica_subsp._enterica_serovar_Oranienburg_str._S_76 will return Salmonella_enterica
                distance = line.split()[2]
                probability = line.split()[3]
                if "Escherichia_coli" in organism:                              # IF E Coli, will parse serotypefinder output for serotype info
                    # ADD SEROTYPE FROM SEROTYPEFINDER
                    if confirm_present(f'{d.serotypefinder}/{sample}/results_table.txt'):
                        serotypefinder_out = open(f'{d.serotypefinder}/{sample}/results_table.txt').read().splitlines()
                        for index, line in enumerate(serotypefinder_out):
                            if line.startswith('H_type'):
                                if serotypefinder_out[index + 1] == 'No serotype predicted.':
                                    h1_type = 'nonePredicted'
                                    h2_type = 'nonePredicted'
                                else:
                                    sero_line = serotypefinder_out[index + 2]   # depending on which gene is listed will dictate which H type is being described
                                    gene = sero_line.split()[0]
                                    if gene == 'fliC':
                                        h1_type = sero_line.split('\t')[5]
                                    elif gene == 'fljB':
                                        h2_type = sero_line.split('\t')[5]
                            elif line.startswith('O_type'):
                                if serotypefinder_out[index + 1] == 'No serotype predicted.':
                                    o_type = 'none'
                                else:
                                    sero_line = serotypefinder_out[index + 2]
                                    o_type = sero_line.split('\t')[5]
                            else:
                                pass
                    else:
                        pass

                elif "Salmonella" in organism:                                  # If salmonella, will parse seqsero2 output for serotype info
                    # ADD SEROTYPE FROM seqsero2
                    if confirm_present(f'{d.seqsero2}/{sample}/Seqsero_result.txt'):
                        seqsero2_out = open(f'{d.seqsero2}/{sample}/Seqsero_result.txt', 'r')
                        for line in seqsero2_out:
                            try:
                                description = line.split('\t')[0].strip()
                                value = line.split('\t')[1].strip()
                                if description == 'O antigen prediction:':
#                                    o_type = value.split('-')[1]                # o types are listed as O-number and only want the number
                                    o_type = value                              # seqsero22 slightly modified output. No longer o-1,3. It is now just 1,3
                                elif description == 'H1 antigen prediction(fliC):':
                                    h1_type = value
                                elif description == 'H2 antigen prediction(fljB):':
                                    h2_type = value
                                elif description == 'Predicted antigenic profile:':
                                    antigenic = value
                                elif description == 'Predicted serotype:':
                                    sero = value
                                else:
                                    pass
                            except:
                                pass
                        seqsero2_out.close()
                    else:
                        pass

                elif "Streptococcus_pyogenes" in organism:                      # If Strep pyogenes, will parse emm_type output for serotype info
                    # ADD SEROTYPE FROM EMM_typing
                    if confirm_present(f'{d.emm_sero}/emm_{sample}/{sample}_PE1.results.xml'):
                        samp_emm_f = open(f'{d.emm_sero}/emm_{sample}/{sample}_PE1.results.xml', 'r')
                        for line in samp_emm_f:
                            if 'Final_EMM_type' in line:
                                value = line.strip().split('value="')[1]
                                sero = value.replace('value="', '').replace('.sds"/>', '')
                            else:
                                pass
                        samp_emm_f.close()
                    else:
                        pass
                elif "Streptococcus_pneumoniae" in organism:                      # If Strep pneumoniae, will parse emm_type output for serotype info
                    # ADD SEROTYPE FROM SEROBA
                    if confirm_present(f'{d.seroba}/{sample}/pred.tsv'):
                        samp_seroba_f = open(f'{d.seroba}/{sample}/pred.tsv', 'r')
                        for line in samp_seroba_f:
                            cols = line.strip().split('\t')
                            samp = cols[0]
                            type = cols[1]
                            if str(samp) == str(sample):
                                sero = type
                            else:
                                pass
                        samp_seroba_f.close()
                    else:
                        pass

                elif "Legionella_pneumophila" in organism:                      # If legionella, will parse legsta output for serotype info
                    # ADD SEROTYPE FROM LEGSTA
                    if confirm_present(f'{d.leg_sero}/legsta_{sample}.tsv'):
                        samp_legsta_f = open(f'{d.leg_sero}/legsta_{sample}.tsv', 'r')
                        for line in samp_legsta_f:
                            if line.startswith('FILE'):
                                pass
                            else:
                                cols = line.strip().split('\t')
                                samp = cols[0]
                                samp_name = samp.split('/')[-1].replace('.fa','')
                                stype = ';'.join(cols[1:])
                                if str(samp_name) == str(sample):
                                    sero = stype
                                else:
                                    pass
                        samp_legsta_f.close()
                    else:
                        pass

                sero_metrics[sample] = [('Organism', organism), ('Distance', distance), ('Probability', probability), ('O', o_type), ('H1', h1_type), ('H2', h2_type), ('AntigenicProfile', antigenic), ('Serotype', sero)]
    else:
        logger(f'{d.mash}/sample_organismIDs.tsv is missing...')

    compiled_file = open(f'{d.base}/compiled_metrics.tsv', 'w')
#    sero_l = ['Organism', 'O', 'H1', 'H2', 'AntigenicProfile', 'Serotype']         # TOTAL OPTIONS FOR sero_l
    sero_l = ['Organism', 'Distance', 'Probability', 'O', 'H1', 'H2', 'Serotype']
#    quast_l = ['Contigs', 'N50', 'GC_content', 'genomeLength', 'totalDepth', 'avgCoverage', 'coverageAboveX']              # TOTAL OPTIONS FOR quast_l
    quast_l = ['Contigs', 'N50', 'GC_content', 'genomeLength', 'avgCoverage', 'coverageAboveX']
#    readmetric_l = ['Sample', 'avgReadLength', 'avgReadQual', 'maxReadLength', 'basesMapped', 'numReads', 'mappedReads', '%ReadsMap', 'avgInsertSize', 'PE']    # TOTAL OPTIONS FOR readmetric_l
    readmetric_l = ['avgReadLength', 'avgReadQual', 'numReads', 'mappedReads', '%ReadsMap', 'avgInsertSize', 'PE']

    all_headers = ['Sample', *sero_l, *quast_l, *readmetric_l]
    for n, header in enumerate(all_headers):                                    # Changing coverage threshold to reflect input from user
        if header == 'coverageAboveX':
            all_headers[n] = f'coverageAbove{p.coverage}X'
    final_header = '\t'.join(all_headers)
    compiled_file.write(f'{final_header}\n')

    for sample in d.sample_list:
        values = [sample]
        for heading in sero_l:
            for entry in sero_metrics[sample]:
                if entry[0] == heading:
                    values.append(entry[1])
        for heading in quast_l:
            for entry in quast_metrics[sample]:
                if entry[0] == heading:
                    values.append(entry[1])
        for heading in readmetric_l:
            for entry in read_metrics[sample]:
                if entry[0] == heading:                                         # first element of each tuple is the metric name, second is the value
                    values.append(entry[1])
        line_format = "\t".join(values)
        compiled_file.write(f'{line_format}\n')
    compiled_file.close()

def warning_flags(d, p):
    """Generates a TSV containing PASS/FAIL for samples based upon user-defined thresholds in config file"""

    if confirm_present(f'{d.base}/compiled_metrics.tsv'):
        compiled = open(f'{d.base}/compiled_metrics.tsv', 'r')
        flags_dict = {}
        for line in compiled:
            if line.startswith('Sample'):
                pass
            else:
                columns = line.strip().split('\t')
                sample = columns[0]
                organism = columns[1]
                mash_prob = float(columns[3])
                n50 = float(columns[9])
                genome_len = float(columns[11])
                avg_coverage = float(columns[12])
                read_quality = float(columns[15])
                reads_mapped = float(columns[18])

                if mash_prob < p.thresh_mash_prob:
                    mash_flag = 'PASS'
                else:
                    mash_flag = 'FAIL'

                if n50 >= p.thresh_n50:
                    n50_flag = 'PASS'
                else:
                    n50_flag = 'FAIL'

                if p.genome_specific_compare.upper() == 'YES':                  # this switch incorporates expected genome size by organism type from the config file
                    genomes_dict = {}
                    sizes = p.specific_genome_sizes.strip().split(';')
                    for entry in sizes:
                        entry_org = entry.split(',')[0].strip()
                        genomesizes = entry.split(',')[1].strip()
                        genomes_dict[entry_org] = genomesizes

                    genome_min = int(p.thresh_genome_sz.strip().split(':')[0])       # Default is going to be that there isnt an entry for this specific organism as some organisms have _OtypeHtype extensions in MASH output
                    genome_max = int(p.thresh_genome_sz.strip().split(':')[1])
                    for org in genomes_dict.keys():                             # if an entry in the dictionary is found to be within the current organism, it wil adjust the genome size to reflect this.
                        if org in organism:
                            genome_exp_sizes = genomes_dict[org]
                            genome_min = float(genome_exp_sizes.split(':')[0])
                            genome_max = float(genome_exp_sizes.split(':')[1])
                            matchfound = True
                else:                                                           # Will compare all genomes to the same size parameters
                    genome_min = int(p.thresh_genome_sz.strip().split(':')[0])
                    genome_max = int(p.thresh_genome_sz.strip().split(':')[1])

                if genome_min <= genome_len and genome_len <= genome_max:
                    genome_flag = 'PASS'
                else:
                    genome_flag = 'FAIL'

                if avg_coverage >= p.thresh_avg_coverage:
                    coverage_flag = 'PASS'
                else:
                    coverage_flag = 'FAIL'

                if read_quality >= p.thresh_read_quality:
                    readqual_flag = 'PASS'
                else:
                    readqual_flag = 'FAIL'

                if reads_mapped >= p.thresh_mapping:
                    readmap_flag = 'PASS'
                else:
                    readmap_flag = 'FAIL'

                flags_dict[sample] = [mash_flag, genome_flag, n50_flag, coverage_flag, readmap_flag, readqual_flag]
        compiled.close()

        flags = open(f'{d.metrics}/warning_flags.tsv', 'w')
        flags.write(f'Sample\tMashPvalue<{p.thresh_mash_prob}\tGenomeSize\tn50>{round(p.thresh_n50/1000)}KB\tAvgCoverage>{p.thresh_avg_coverage}X\tReadMapping>{p.thresh_mapping}%\tReadQuality>{p.thresh_read_quality}\n')
        for sample in d.sample_list:
            values = flags_dict[sample]
            cols = [sample, *values]
            line_form = '\t'.join(cols)
            flags.write(f'{line_form}\n')
        flags.close()

        flags = open(f'{d.metrics}/warning_flags.tsv', 'r')
        flags_mqc = open(f'{d.for_multiqc}/warning_flags_multiqc.tsv', 'w')
        for line in flags:
            flags_mqc.write(line.replace('PASS', '0.1').replace('FAIL', '1'))
        flags.close()
        flags_mqc.close()

    else:
        logger(f'Cannot find compiled_metrics.tsv file in {d.base}. Quitting.')
        quit()

def extract_4_multiqc(d, p):
    """Extracts data from numerous locations and compiles in 'multiqc_input'"""

    make_directory(d.for_multiqc)

    # multiqc directory did not exist during metaphlan run. moving composition tsv to multiqc dir now
    if confirm_present(f'{d.metaphlan2}/metagenomic_composition_multiqc.tsv'):
        shutil.move(f'{d.metaphlan2}/metagenomic_composition_multiqc.tsv', f'{d.for_multiqc}/metagenomic_composition_multiqc.tsv')
    elif not confirm_present(f'{d.metaphlan2}/metagenomic_composition_multiqc.tsv') and confirm_present(f'{d.for_multiqc}/metagenomic_composition_multiqc.tsv'):
        pass
    elif not confirm_present(f'{d.metaphlan2}/metagenomic_composition_multiqc.tsv') and not confirm_present(f'{d.for_multiqc}/metagenomic_composition_multiqc.tsv'):
        logger('Cannot find metagenomic_composition_multiqc.tsv file. Make sure MetaPhlAn2 completed successfully.')
        quit()

    # GENERATING A ORGANISM DESCRIPTION FOR MULTIQC
    org_description = open(f'{d.for_multiqc}/organismal_specs.tsv', 'w')
    assembly_description = open(f'{d.for_multiqc}/assembly_specs.tsv', 'w')
    compiled = open(f'{d.base}/compiled_metrics.tsv', 'r')
    for line in compiled:                                                       # Extracting sample name, H and O types, and serotype into one file
        columns = line.split('\t')
        sampName_species = columns[:2]
        if columns[0] == 'Sample':                                              # not adjusting header line
            oh12Type_Sero = columns[4:8]
        else:
            oh12Type_Sero = []
            oh12Type_Sero_1st = [value.replace(',', ';') for value in columns[4:8]] # Replacing commas with semicolons so MultiQC recognizes them, otherwise 1 2 is less visible. Now it will be 1;2 if there are two H, O, or serotypes predicted
            for element in oh12Type_Sero_1st:
                if element != '-':
                    oh12Type_Sero.append(f'{element};')                         # appending semicolon to the end of each element. this keeps multiqc from formatting numbers with a colored abundance bar behind it. (i.e. SBT for Lp may be 36.0)
                else:
                    oh12Type_Sero.append(element)
        decriptions = [*sampName_species, *oh12Type_Sero]
        descriptions_out = '\t'.join(decriptions)
        org_description.write(f'{descriptions_out}\n')
        assembly_mets = [columns[0], *columns[9:]]                              # Extracting sample name and assembly calculations into another file.
        assembly_met_l = '\t'.join(assembly_mets)
        assembly_met_l += '\n'
        assembly_description.write(assembly_met_l)
    org_description.close()
    assembly_description.close()

    # GENERATING MASH OUTPUT FOR MULTIQC
    all_orgs = []
    all_counts = {}                                                             # Will be a nested dictionary. First key is sample ID, second key is genus_species, value is count
    for sample in d.sample_list:
        sample_orgs = {}
        mash_f = open(f'{d.mash}/{sample}.distance.tab', 'r')
        for line in mash_f:
            pval = float(line.split('\t')[3])
            if pval == 0:                                                       # only want those that we are cetrtain of (pval = o)
                name_str = line.split('\t')[0]                                  # isolating first column ( refseq-NC-1001582-PRJNA224116-SAMN02603189-GCF_000204275.1-pMC1-Bacillus_amyloliquefaciens_LL3.fna )
                fna = name_str.split('-')[-1][:-4]                              # Elements are separated by '-', and want to remove .fna from the end of the last element
                genus = fna.split('_')[0]
                if fna.split('_')[1] =='sp.':                                   # this captures the specificity of KTE92 if fna is Klebsiella_sp._KTE92.fna
                    species = 'sp'                                              # whereas if the fna is Klebsiella_pneumoniae_KCTC_2242.fna, it will capture Klebsiella_pneumoniae
                else:
                    species = fna.split('_')[1]

                gen_spec = '_'.join([genus, species])                           # First two elements are Genus species, Isolating them and joining back with underscore
                if gen_spec not in sample_orgs.keys():                          # If gen_species is not already in dictionary, add it. Otherwise, add to existing count.
                    sample_orgs[gen_spec] = 1
                else:
                    sample_orgs[gen_spec] += 1
                if gen_spec not in all_orgs:
                    all_orgs.append(gen_spec)
        all_counts[sample] = sample_orgs

    mash_multiqc = open(f'{d.for_multiqc}/mash_abundance.tsv', 'w')               # Writing to file for multiqc (Tab separated abundance file)
    header = ['Sample', *all_orgs]
    header = '\t'.join(header)
    header += '\n'
    mash_multiqc.write(header)
    for sample in d.sample_list:
        line = [sample]
        org_cts = all_counts[sample]
        for query_org in all_orgs:
            if query_org in org_cts.keys():
                count = org_cts[query_org]
            else:
                count = 0
            line.append(str(count))
        line = '\t'.join(line)
        sample_line = f'{line}\n'
        mash_multiqc.write(sample_line)
    mash_multiqc.close()

    # GENERATING ABRICATE OUTPUT FOR MULTIQC
    abricate_dbs = []
    cmd = f'singularity exec {p.singularity_loc}/{p.abricate_sif} abricate --list'
    run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    for line in run.stdout:                                                     # Generating a list of the DBs present in abricate container
        db_name = line.split('\t')[0]                                           # DB name is the first element of each line (though very first line is "DATABASE")
        if db_name != 'DATABASE':                                               # ISSUE: the 'abricate' db appears to not have been established...
            abricate_dbs.append(db_name)

    abricate_db_line = '\t'.join(abricate_dbs)
    abricate_counts = open(f'{d.for_multiqc}/abricate_total_counts.tsv', 'w')
    counts_d = {}
    for DB in abricate_dbs:
        ab_summary = open(f'{d.metrics}/abricate_{DB}_summary', 'r')
        abricate_multiqc = open(f'{d.for_multiqc}/abricate_{DB}_multiqc.tsv', 'w')
        for line in ab_summary:
            if line.startswith('#'):
                cols = line.split('\t')
                first = cols[0].replace('#', '')
                rest = cols[2:]                                                 # Second column contains total number of AR genes found. Going to use this for another output
                outline_l = [first, *rest]
                outline = '\t'.join(outline_l)
            else:
                cols = line.strip().split('\t')                                 # Very last element has newline, so without stripping it, I wasnt getting the == '.' function below to work
                sample = cols[0]
                remaining_cols = cols[2:]
                sample_new = sample.split('/')[-1].replace(f'_{DB}.tab', '')
                outline_l = [sample_new, *remaining_cols]
                for index, value in enumerate(outline_l):                       # If abricate doesnt find a gene, the entry is '.', which I need to replace with 0 for the heatmap
                    if value == '.':
                        outline_l[index] = '0'
                    elif ';' in value:                                          # numerous places have several hits for different contigs ('12.44;100.00;19.49'), so this takes the
                        numbers = value.split(';')                              # highest value, as all the user will want from the multiqc output is which genes are found, not location
                        numbers = list(map(float, numbers))
                        highest = max(numbers)
                        outline_l[index] = str(highest)
                outline = '\t'.join(outline_l)
                outline += '\n'

                total_num = cols[1]
                if sample_new not in counts_d.keys():
                    counts_d[sample_new] = [(f'{DB}', total_num)]
                else:
                    counts_d[sample_new].append((f'{DB}', total_num))
            abricate_multiqc.write(outline)
        abricate_multiqc.close()
    counts_header_l = ['Sample', *abricate_dbs]
    counts_header = '\t'.join(counts_header_l)
    abricate_counts.write(f'{counts_header}\n')
    for sample in counts_d.keys():
        line = [sample]
        for DB in abricate_dbs:
            for db_ct in counts_d[sample]:                                      # Each sample now has a list of tuples that are (DB, count)
                count_db = db_ct[0]
                if count_db == DB:
                    count = db_ct[1]
                    line.append(count)
                else:
                    pass
        line_out = '\t'.join(line)
        abricate_counts.write(f'{line_out}\n')
    abricate_counts.close()

    # GENERATING HEATMAPS FOR ABRICATE FOR PHENOTYPE RESISTANCE
    for DB in abricate_dbs:
        DB_pheno_counts = {}
        resistance_l = []
        for sample in d.sample_list:
            sample_pheno_counts = {}
            samp_abricate = open(f'{d.abricate}/{sample}_{DB}.tab', 'r')
            for line in samp_abricate:
                if line.startswith('#FILE'):
                    pass
                else:
                    resistances = []
                    if len(line.strip().split('\t')) == 15:                     # several of the DBs dont regularly ID resistance phenotype. This captures that
                        resistance = line.strip().split('\t')[14]               # resistance phenotype is last element of line
                        if ';' in resistance:
                            resistances = resistance.split(';')                 # if multiple resistances are given, they are separated by semi-colon
                        elif '/' in resistance:
                            resistances = resistance.split('/')
                        else:
                            resistances = [resistance]
                    else:
                        resistances = ['NonePredicted']
                    for resistance in resistances:
                        if resistance in sample_pheno_counts.keys():
                            sample_pheno_counts[resistance] += 1
                        else:
                            sample_pheno_counts[resistance] = 1
                        if resistance not in resistance_l:
                            resistance_l.append(resistance)
            DB_pheno_counts[sample] = sample_pheno_counts
        db_pheno_out = open(f'{d.for_multiqc}/{DB}_drug_resistance.tsv', 'w')
        headers = [*['SampleName'], *resistance_l]
        first_line = '\t'.join(headers) + '\n'
        db_pheno_out.write(first_line)
        for sample in d.sample_list:
            line = [sample]
            samp_counts_dir = DB_pheno_counts[sample]
            for resistance in resistance_l:
                if resistance in samp_counts_dir.keys():
                    count = samp_counts_dir[resistance]
                    line.append(str(count))
                else:
                    line.append(str(0))
            line_format = '\t'.join(line) + '\n'
            db_pheno_out.write(line_format)
        db_pheno_out.close()

def run_multiqc(d, p):
    """Runs multiqc and several custom scripts to obtain visual representation of pipeline findings"""

    logger('RUNNING MULTIQC TO COMPILE METRICS AND GENERATE VISUALS')
    if p.multiqc_title == 'Default':                                            # Default will name MultiQC report at the DTG it was generated
        p.multiqc_title = str(datetime.datetime.now()).replace(' ','_')
    if p.multiqc_fname == 'Default':                                            # Default will name MultiQC report at the DTG it was generated
        mqc_prefix = str(datetime.datetime.now()).replace(' ','_').split('.')[0]
        p.multiqc_fname = f'{mqc_prefix}_multiqcReport.html'
    cmd = f'singularity exec {p.singularity_loc}/{p.multiqc_sif} multiqc \
    -c {d.base}/BRR_multiqc_config.yaml --exclude general_stats --title {p.multiqc_title} \
    --comment {p.multiqc_comment} --filename {p.multiqc_fname} {p.multiqc_add} \
    {d.fastqc} {d.trim} {d.quast} {d.for_multiqc}'
    run_command_logger(cmd)

def run_pipeline(dirs, cfg, cfg_f, VERSION):
    """This function runs the main pipeline"""

    if os.path.isfile(f'{dirs.base}/logfile.txt'):
        print('Logfiles exist from previous run. Preappending "prev_" to those and generating new logfile\n')
        os.rename(f'{dirs.base}/logfile.txt', f'{dirs.base}/prev_logfile.txt')
        os.rename(f'{dirs.base}/logfile_full.txt', f'{dirs.base}/prev_logfile_full.txt')
    else:
        pass

    file_full = open('logfile_full.txt', 'a')
    logger(f'''\n****************************************************************************************************
You are running QC and Typing analysis version {VERSION}. This pipeline is maintained by Andrew Lang
at the Massachusetts State Public Health Lab. Additional information pertaining to pipeline usage can
be found at https://gitlab.com/ma_ngs/brr_pipeline. For any issues with program operation, please
contact Andrew at Andrew.Lang@massmail.state.ma.us for assistance.
****************************************************************************************************\n''')
    time.sleep(1)
    logger('\n~ - ~ - ~ - ~ - ~ - CONFIG FILE PARAMETERS FOR THIS ANALYSIS ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ -\n')
    config = open(cfg_f, 'r')
    for line in config:
        logger(line.rstrip())
    file_full.close()
    logger('\n\n~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - CONFIG END ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~\n\n')
    time.sleep(1)

    start = datetime.datetime.now()
    logger(f'Beginning analysis at {str(start).split(".")[0]}')
    logger('\n...INFRASTRUCTURE CHECK...\n')
    if os.path.isfile(f'{dirs.base}/BRR_multiqc_config.yaml'):
        logger('Found MultiQC yaml config file already present in current directory. Will use for MultiQC metric output.\n')
    else:
        brr_foundation.generate_multiqcConfig(dirs)
        logger(f'No .yaml config file for MultiQC found in {dirs.base}. This file has now been created.\n')

    file = open('logfile.txt', 'w')
    brr_foundation.get_memory(dirs, cfg)
    brr_foundation.singularity_check(cfg)
    brr_foundation.database_check(cfg)
    make_directory(dirs.metrics)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    pull_SRRs(dirs, cfg)
    move_fastqs(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    trim_reads(dirs, cfg)
    dirs.sample_list = get_sample_list(dirs.trim)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    shuffle_reads(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    run_fastqc(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    run_metaphlan2(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    run_shovill(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    run_quast(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    run_bwa(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    shuffle_reads(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    run_mash(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    organism_switch(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    run_abricate(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    QualCov_metrics(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    compile_metrics(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    extract_4_multiqc(dirs, cfg)
    warning_flags(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    run_multiqc(dirs, cfg)
    logger('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
    stop = datetime.datetime.now()
    logger(f'Analysis completed at {str(stop).split(".")[0]}')
    duration = stop - start
    logger(f'Duration: {str(duration).split(".")[0]}')

if __name__ == '__main__':

#search for ISSUE, OPTIMIZATION, QUESTION at the end of this

    logger('')
    args = parser.parse_args()
    config_f = args.config
    dirs = brr_foundation.gen_directories()

    if not os.path.isfile(config_f):
        print(f'Cannot find your specified config file at {config_f}')
        brr_foundation.gen_config(dirs)
        print('New config file generated in current directory. Confirm parameters are correct before running analysis.')
        exit()
    else:
       try:
           cfg = brr_foundation.parse_config(config_f)
           VERSION = '10JAN2019'
           run_pipeline(dirs, cfg, config_f, VERSION)
       except Exception:
           traceback.print_exc()
           exit()
