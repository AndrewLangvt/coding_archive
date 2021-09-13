#!/usr/bin/env python3

# - brr_foundation.py
# - Andrew S. Lang
# - Created: 20OCT2019
# - Last Modified: 24FEB2020

import os
import sys
import glob
import argparse
import configparser
from shutil import which
import subprocess

class Parameters:
    """Defines all aspects of config file as a class, enabling passing to any function."""

    def __init__(self):
        self.singularity_loc = ''
        self.cores = ''
        self.sample_list = []

def gen_config(d):
    """Writes new config file if not present in specified directory. Writes to CWD."""

    config = open(f'{d.base}/brr_config.ini', 'w')
    config.write("""
[Locations]

# - Location for singularity image files
singularity_loc = ~/singularity_files
# - Location for MetaPhlAn2 databses
metaphlan2_loc = ~/metaphlan2_db

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

[Thresholds]
# - Coverage: when calculating percent genome coverage, determine what level is desired. I.e. "20" will calculate
# - the proportion of genome covered by at least 20 reads. NOTE: this is different from the below "WarningFlags" as
# - this is part of the compiled_metrics, while below will output a heatmap of samples that do/do not meet those thresholds
coverage = 20

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

[WarningFlags]

# - GENOME SIZE - base pair size range of expected genome (for all organisms not in below "SpecificGenomeSizes" section)
thresh_genome_sz = 1800000:2200000
# - MAPPING RATE - percent reads mapped to assembled genome
thresh_mapping = 99.50
# - MASH PROBABILITY - p-value that the predicted organism accurately reflects your sample
thresh_mash_prob = 0.05
# - AVERAGE COVERAGE - average read coverage of assembled genome
thresh_avg_coverage = 20.00
# - AVERAGE READ QUALITY - average quality of trimmed reads
thresh_read_quality = 30.00
# - N50 Value for assembly - N50 is the size of the contig which, when added with all of the larger contigs, encompass 50% of the total genome size.
thresh_n50 = 250000

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

[SpecificGenomeSizes]
genome_specific_compare = YES
# - These are expected genome sizes for specific organisms. The pipeline will identify the organism with MASH, and compare it to these
# - The user can add to this list, but must do so in a specific format. To do this, add a semicolon and your organism (as identified by MASH)
# - along with the genome size- separated by a comma. Usually, MASH outputs organisms in the Genus_species format. i.e. to add in an
# - expected genome size for My_genome at 100kb, I would add "; My_genome, 100000" to the end of the list
specific_genome_sizes = Salmonella_enterica, 4400000:5200000; Salmonella_bongori, 4200000:4600000; Streptococcus_pyogenes, 1600000:2000000; Streptococcus_pneumoniae, 1800000:2200000; Escherichia_coli, 4400000:4800000; Legionella_pneumophila, 3200000:3600000; Campylobacter_jejuni, 1500000:1700000; Klebsiella_pneumoniae, 5000000:5400000
# - For organisms not included in this list, the pipeline will default to the thresh_genome_sz parameter in "WarningFlags" above
# - If you would prefer to simply set the same threshold for "warnings" for all assemblies, the genome_specific_compare parameter can be adjusted to "NO" and
# - the expected genome size range will be 1800000:2200000

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

[Containers]

sratoolkit = staphb/sratoolkit:2.9.2
seqyclean = staphb/seqyclean:1.10.09
fastqc = staphb/fastqc:0.11.8
metaphlan2 = andrewlangvt/metaphlan2:2
shovill = staphb/shovill:1.0.4-cv2
quast = staphb/quast:5.0.2
bwa = staphb/bwa:0.7.17
samtools = staphb/samtools:1.9
lyveset = staphb/lyveset:1.1.4f
mash = staphb/mash:2.1
serotypefinder = staphb/serotypefinder:1.1
seqsero2 = staphb/seqsero2:1.0.0
sistr = staphb/sistr:1.0.2
emmtyping = staphb/emm-typing-tool:0.0.1
# legtyping = smorrison42/lpserogroup_prediction:0.2
seroba = staphb/seroba:1.0.0
legsta = staphb/legsta:0.3.7-cv1
shigatyper = andrewlangvt/shigatyper:1
abricate = staphb/abricate:0.9.8-cv1
multiqc = ewels/multiqc:1.7
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

[RunSpecs]

# - Cores and RAM will be automatically calculated to 75%
# - of machines capability if you set them equal to "0"
cores = 0
RAM = 0
threads = 1
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

[Commands]

# - This section contains options to modify commands. Elements contained with "<>" are modifiable by the lines following.
# - Elements contained within "{}" are not modifiable to ensure correct I/O flow through pipeline.

# - FOR ADVANCED USERS WISHING TO FURTHER MODIFY COMMANDS: At the end of each block for a given analysis tool, there is a "tool_add"
# - line wherein the user can add desired parameters to the analysis. MODIFY AT YOUR OWN RISK. It is highly likely that
# - alterations to output file type (i.e. changing CSV output to TSV) will negatively impact downstream metric compiling
# - and data visualization. But hey, live life on the edge. You do you.

# - ANALYSIS TOOL THRESHOLDS

# -------------------------------------------------------------------------------------------------------------------
# - SEQYCLEAN
# - seqyclean -minlen <seqyclean_minlen> -qual <seqyclean_qual> -c <seqyclean_adapters> -1 {d.raw}/{isolate_r1} -2 {d.raw}/{isolate_r2} -o {d.trim}/{base}
seqyclean_minlen = 25
# - turns off/on quality trimming. Replace with "N" to turn off.
seqyclean_qualTrim = Y
# - Set specific error levels. Error boundaries: max_average_error (default=20 Phred), max_error_at_ends (default=20 Phred).
seqyclean_qual = 20 20
seqyclean_adapters = /Adapters_plus_PhiX_174.fasta
seqyclean_add =
# -------------------------------------------------------------------------------------------------------------------
# - FASTQC
# - fastqc <file> -t <threads> -o <fastqc_dir>
fastqc_add =
# -------------------------------------------------------------------------------------------------------------------
# - METAPHLAN2
# - metaphlan2.py {d.shuffle}/{sample}.shuffled.fastq.gz --bowtie2out {d.metaphlan2}/{sample}.bowtie2out.txt --input_type fastq --nproc {p.cores} <metaphlan_add> > {d.metaphlan2}/{sample}_profile.txt
metaphlan_add =
# -------------------------------------------------------------------------------------------------------------------
# - SHOVILL
# - shovill --cpus {p.cores} --ram {p.ram} --minlen <shovill_minlen> --mincov <shovill_mincov> <shovill_add> --outdir {d.shovill}/{sample} --R1 {d.trim}/{isolate_r1} --R2 {d.trim}/{isolate_r2}
# - Minimum contig length <0=AUTO>
shovill_minlen = 0
# - Minimum contig coverage <0=AUTO> (default: 2)
shovill_mincov = 2
shovill_add =
# -------------------------------------------------------------------------------------------------------------------
# - QUAST
# - quast.py <quast_flags> -t {p.threads} <quast_add> {d.shovill}/{sample}/{sample}.fa -o {d.quast}/{sample}
quast_nocheck = Y
quast_noplots = Y
quast_nohtml = Y
quast_noicarus = Y
quast_nosnps = Y
quast_nogc = N
quast_nosv = Y
quast_nogzip = Y
quast_noreadstats = N
quast_add =
# -------------------------------------------------------------------------------------------------------------------
# - BWA
# - bwa index <bwa_index_add> {d.shovill}/{sample}/{sample}.fa
bwa_index_add =
# - bwa mem -t {p.threads} <bwa_mem_add> {d.shovill}/{sample}/{sample}.fa {d.trim}/{isolate_r1} {d.trim}/{isolate_r2} > {d.bwa}/{sample}.alignment.sam
bwa_mem_add =
# -------------------------------------------------------------------------------------------------------------------
# - SAMTOOLS
# - samtools sort {d.bwa}/{sample}.alignment.sam --threads {p.threads} -T {d.bwa}/temp_files <samtools_add> -O BAM -o {d.bwa}/{sample}.alignment.sorted.bam
samtools_add =
# -------------------------------------------------------------------------------------------------------------------
# - MASH
# - mash sketch -m <mash_sketch_m> <mash_sketch_add> {d.shuffle}/{sample}.shuffled.fastq.gz
# - minimum copies of each k-mer required. 2 will filter out unique
mash_sketch_m = 2
mash_sketch_add =
# - mash dist <mash_dist_db> <mash_dist_add> {d.mash}/{sample}.shuffled.fastq.gz.msh > {d.mash}/{sample}.distance.tab
mash_dist_db = /db/RefSeqSketchesDefaults.msh
mash_dist_add =
# -------------------------------------------------------------------------------------------------------------------
# - SEROTYPEFINDER
# - serotypefinder.pl -d <serotypefinder_db> -i {d.shovill}/{sample}/{sample}.fa -b <sertypefinder_blast> -o {d.sero}/{sample} -s <serotypefinder_species> -k <serotypefinder_threshold> -l <serotypefinder_minleng> <serotypefinder_add>
serotypefinder_db = /serotypefinder/database
serotypefinder_blast = /blast-2.2.26
# - currently, only ecoli is an option for species
serotypefinder_species = ecoli
serotypefinder_threshold = 95.00
serotypefinder_minleng = 0.60
serotypefinder_add =
# -------------------------------------------------------------------------------------------------------------------
# - SEQSERO2
# - SeqSero2_package.py -t4 -m <seqsero2_mode> -i {d.shovill}/{sample}/{sample}.fa -p {p.threads} -b <seqsero2_mapping> -d {d.seqsero2}/{sample}'
# - mode (kmer or allele)
seqsero2_mode = k
# - mapping algorithm (sam or mem)
seqsero2_mapping = mem
# -------------------------------------------------------------------------------------------------------------------
# - EMM_TYPING
# emm_typing.py --profile_file_directory <emm_db> --fastq_1 {d.trim}/{sample}_PE1.fastq.gz --fastq_2 {d.trim}/{sample}_PE2.fastq.gz --output_directory {d.emm_sero}/emm_{sample} <emm_add>'
emm_db = /db/
emm_add =
# -------------------------------------------------------------------------------------------------------------------
# - SEROBA
# - seroba runSerotyping <seroba_db> {d.trim}/{sample}_PE1.fastq.gz {d.trim}/{sample}_PE2.fastq.gz {sample} --coverage <seroba_cvg> <seroba_add>
seroba_db = /seroba-1.0.0/database
seroba_cvg = 20
seroba_add =
# -------------------------------------------------------------------------------------------------------------------
# - LEGSTA
# - legsta <legsta_add> {d.shovill}/{sample}/{sample}.fa > {d.leg_sero}/legsta_{sample}.tsv
legsta_add =
# -------------------------------------------------------------------------------------------------------------------
# - SHIGATYPER
# - there currently are no additional flags that can be added by the user
# -------------------------------------------------------------------------------------------------------------------
# - ABRICATE
# - abricate --minid <abricate_minid> --mincov <abricate_mincov> --threads {p.threads} -db {DB} {d.shovill}/{sample}/{sample}.fa > {d.abricate}/{sample}_{DB}.tab
# - NOTE: there is no option ot alter the abricate DB, as the pipeline automatically queries all DBs included with abricate
abricate_minid = 75
abricate_mincov = 0
abricate_add =
# -------------------------------------------------------------------------------------------------------------------
# - MULTIQC
# - multiqc -c {d.base}/BRR_multiqc_config.yaml --exclude general_stats ----title <multiqc_title> --comment <multiqc_comment> <multiqc_add> {d.fastqc} {d.trim} {d.quast} {d.for_multiqc}'
# - multiqc filename Default: Date_time_multiqcReport.html
multiqc_fname = Default
# - Report title. Default: uses Date-Time group for titling
multiqc_title = Default
# - Include a custom comment at top of report. NOTE: MUST be encompassed by single quotes. Default: MultiQC Report
multiqc_comment = 'MultiQC Report'
multiqc_add =
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------
""")
    config.close()

def parse_config(config_f):
    """Parse a config file, assigning values to the class attributes"""

    logger(f'Parsing config file {config_f}')
    config_parsed = Parameters()
    parser = configparser.ConfigParser()
    parser.read(config_f)
    config_parsed.singularity_loc = os.path.expanduser(parser.get('Locations', 'singularity_loc'))
    config_parsed.metaphlan2_loc = os.path.expanduser(parser.get('Locations', 'metaphlan2_loc'))
    config_parsed.sratoolkit = parser.get('Containers', 'sratoolkit')
    config_parsed.seqyclean = parser.get('Containers', 'seqyclean')
    config_parsed.fastqc = parser.get('Containers', 'fastqc')
    config_parsed.metaphlan2 = parser.get('Containers', 'metaphlan2')
    config_parsed.shovill = parser.get('Containers', 'shovill')
    config_parsed.quast = parser.get('Containers', 'quast')
    config_parsed.bwa = parser.get('Containers', 'bwa')
    config_parsed.samtools = parser.get('Containers', 'samtools')
    config_parsed.lyveset = parser.get('Containers', 'lyveset')
    config_parsed.mash = parser.get('Containers', 'mash')
    config_parsed.seqsero2 = parser.get('Containers', 'seqsero2')
    config_parsed.emmtyping = parser.get('Containers', 'emmtyping')
    config_parsed.legsta = parser.get('Containers', 'legsta')
    config_parsed.serotypefinder = parser.get('Containers', 'serotypefinder')
    config_parsed.seroba = parser.get('Containers', 'seroba')
    config_parsed.shigatyper = parser.get('Containers', 'shigatyper')
    config_parsed.abricate = parser.get('Containers', 'abricate')
    config_parsed.multiqc = parser.get('Containers', 'multiqc')
    config_parsed.cores = parser.get('RunSpecs', 'cores')
    config_parsed.ram = parser.get('RunSpecs', 'RAM')
    config_parsed.threads = parser.get('RunSpecs', 'threads')

    # Generating names for each sif. They will contain the version as well. I.e. docker://sratoolkit:2.9.2 will be named sratoolkit_2.9.2.sif
    config_parsed.sratoolkit_sif = f"{config_parsed.sratoolkit.split('/')[1].replace(':', '_')}.sif"
    config_parsed.seqyclean_sif = f"{config_parsed.seqyclean.split('/')[1].replace(':', '_')}.sif"
    config_parsed.fastqc_sif = f"{config_parsed.fastqc.split('/')[1].replace(':', '_')}.sif"
    config_parsed.metaphlan2_sif = f"{config_parsed.metaphlan2.split('/')[1].replace(':', '_')}.sif"
    config_parsed.shovill_sif = f"{config_parsed.shovill.split('/')[1].replace(':', '_')}.sif"
    config_parsed.quast_sif = f"{config_parsed.quast.split('/')[1].replace(':', '_')}.sif"
    config_parsed.bwa_sif = f"{config_parsed.bwa.split('/')[1].replace(':', '_')}.sif"
    config_parsed.samtools_sif = f"{config_parsed.samtools.split('/')[1].replace(':', '_')}.sif"
    config_parsed.lyveset_sif = f"{config_parsed.lyveset.split('/')[1].replace(':', '_')}.sif"
    config_parsed.mash_sif = f"{config_parsed.mash.split('/')[1].replace(':', '_')}.sif"
    config_parsed.seqsero2_sif = f"{config_parsed.seqsero2.split('/')[1].replace(':', '_')}.sif"
    config_parsed.emmtyping_sif = f"{config_parsed.emmtyping.split('/')[1].replace(':', '_')}.sif"
    config_parsed.legsta_sif = f"{config_parsed.legsta.split('/')[1].replace(':', '_')}.sif"
    config_parsed.serotypefinder_sif = f"{config_parsed.serotypefinder.split('/')[1].replace(':', '_')}.sif"
    config_parsed.seroba_sif = f"{config_parsed.seroba.split('/')[1].replace(':', '_')}.sif"
    config_parsed.shigatyper_sif = f"{config_parsed.shigatyper.split('/')[1].replace(':', '_')}.sif"
    config_parsed.abricate_sif = f"{config_parsed.abricate.split('/')[1].replace(':', '_')}.sif"
    config_parsed.multiqc_sif = f"{config_parsed.multiqc.split('/')[1].replace(':', '_')}.sif"

    config_parsed.coverage = parser.get('Thresholds', 'coverage')

    config_parsed.thresh_genome_sz = parser.get('WarningFlags', 'thresh_genome_sz')
    config_parsed.thresh_mapping = float(parser.get('WarningFlags', 'thresh_mapping'))
    config_parsed.thresh_mash_prob = float(parser.get('WarningFlags', 'thresh_mash_prob'))
    config_parsed.thresh_avg_coverage = float(parser.get('WarningFlags', 'thresh_avg_coverage'))
    config_parsed.thresh_read_quality = float(parser.get('WarningFlags', 'thresh_read_quality'))
    config_parsed.thresh_n50 = float(parser.get('WarningFlags', 'thresh_n50'))

    config_parsed.specific_genome_sizes = parser.get('SpecificGenomeSizes', 'specific_genome_sizes')
    config_parsed.genome_specific_compare = parser.get('SpecificGenomeSizes', 'genome_specific_compare')

    config_parsed.seqyclean_minlen = parser.get('Commands', 'seqyclean_minlen')
    config_parsed.seqyclean_qualTrim = parser.get('Commands', 'seqyclean_qualTrim')
    config_parsed.seqyclean_qual = parser.get('Commands', 'seqyclean_qual')
    config_parsed.seqyclean_adapters = parser.get('Commands', 'seqyclean_adapters')
    config_parsed.seqyclean_add = parser.get('Commands', 'seqyclean_add')

    config_parsed.fastqc_add = parser.get('Commands', 'fastqc_add')

    config_parsed.metaphlan_add = parser.get('Commands', 'metaphlan_add')

    config_parsed.shovill_minlen = parser.get('Commands', 'shovill_minlen')
    config_parsed.shovill_mincov = parser.get('Commands', 'shovill_mincov')
    config_parsed.shovill_add = parser.get('Commands', 'shovill_add')

    config_parsed.bwa_index_add = parser.get('Commands', 'bwa_index_add')
    config_parsed.bwa_mem_add = parser.get('Commands', 'bwa_mem_add')
    config_parsed.samtools_add = parser.get('Commands', 'samtools_add')

    config_parsed.quast_nocheck = parser.get('Commands', 'quast_nocheck')
    config_parsed.quast_noplots = parser.get('Commands', 'quast_noplots')
    config_parsed.quast_nohtml = parser.get('Commands', 'quast_nohtml')
    config_parsed.quast_noicarus = parser.get('Commands', 'quast_noicarus')
    config_parsed.quast_nosnps = parser.get('Commands', 'quast_nosnps')
    config_parsed.quast_nogc = parser.get('Commands', 'quast_nogc')
    config_parsed.quast_nosv = parser.get('Commands', 'quast_nosv')
    config_parsed.quast_nogzip = parser.get('Commands', 'quast_nogzip')
    config_parsed.quast_noreadstats = parser.get('Commands', 'quast_noreadstats')
    config_parsed.quast_add = parser.get('Commands', 'quast_add')

    config_parsed.mash_sketch_m = parser.get('Commands', 'mash_sketch_m')
    config_parsed.mash_sketch_add = parser.get('Commands', 'mash_sketch_add')
    config_parsed.mash_dist_db = parser.get('Commands', 'mash_dist_db')
    config_parsed.mash_dist_add = parser.get('Commands', 'mash_dist_add')

    config_parsed.serotypefinder_db = parser.get('Commands', 'serotypefinder_db')
    config_parsed.serotypefinder_blast = parser.get('Commands', 'serotypefinder_blast')
    config_parsed.serotypefinder_species = parser.get('Commands', 'serotypefinder_species')
    config_parsed.serotypefinder_threshold= parser.get('Commands', 'serotypefinder_threshold')
    config_parsed.serotypefinder_minleng = parser.get('Commands', 'serotypefinder_minleng')
    config_parsed.serotypefinder_add = parser.get('Commands', 'serotypefinder_add')

    config_parsed.seqsero2_mode = parser.get('Commands', 'seqsero2_mode')
    config_parsed.seqsero2_mapping = parser.get('Commands', 'seqsero2_mapping')

    config_parsed.emm_db = parser.get('Commands', 'emm_db')
    config_parsed.emm_add = parser.get('Commands', 'emm_add')

    config_parsed.seroba_db = parser.get('Commands', 'seroba_db')
    config_parsed.seroba_cvg = parser.get('Commands', 'seroba_cvg')
    config_parsed.seroba_add = parser.get('Commands', 'seroba_add')

    config_parsed.legsta_add = parser.get('Commands', 'legsta_add')

    config_parsed.abricate_minid = parser.get('Commands', 'abricate_minid')
    config_parsed.abricate_mincov = parser.get('Commands', 'abricate_mincov')
    config_parsed.abricate_add = parser.get('Commands', 'abricate_add')

    config_parsed.multiqc_fname = parser.get('Commands', 'multiqc_fname')
    config_parsed.multiqc_title = parser.get('Commands', 'multiqc_title')
    config_parsed.multiqc_comment = parser.get('Commands', 'multiqc_comment')
    config_parsed.multiqc_add = parser.get('Commands', 'multiqc_add')

    logger('Parsing complete.\n')
    return config_parsed

def generate_multiqcConfig(d):
    """Generates yml file for MultiQC run to incorporate custom metrics"""

    if os.path.isfile(f'{d.base}/BRR_multiqc_config.yaml'):
        pass
    else:
        mqc_config = open(f'{d.base}/BRR_multiqc_config.yaml', 'w')
        mqc_config.write("""custom_content:
  order:
    - warnings
    - organismal_specs
    - assembly_specs
    - mash_out
    - metaphlan_out
    - abricate_out
    - abricate_argannot
    - abricate_ecoh
    - abricate_plasmidfinder
    - abricate_vfdb
    - abricate_card
    - abricate_ecoli_vf
    - abricate_ncbi
    - abricate_resfinder
    - card_drugResistance
    - ncbi_drugResistance
    - resfinder_drugResistance

custom_data:
  organismal_specs:
    id: 'organismal_specs'
    section_name: 'ORGANISMAL INFO'
    description: 'contains information about each sample pertaining to organism and serotype'
    plot_type: 'table'
    format: 'tsv'
    pconfig:
      id: 'organismal_specs'
      title: "Organism Information by Sample"
  assembly_specs:
    id: 'assembly_specs'
    section_name: 'ASSEMBLY SPECS'
    description: 'for assembly qualtiy, composition, and read alignment'
    plot_type: 'table'
    format: 'tsv'
    scale: False
    pconfig:
      id: 'assembly_specs'
      title: "Assembly Information by Sample"
  warnings:
    id: 'warnings'
    section_name: 'WARNINGS'
    description: 'Exhibits where your data did or did not meet the thresholds input via config file'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'warngings'
      title: 'Pass(grey) / Fail(black) of your data'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [1.0, '#000000'], ]
      use_legend: FALSE
  mash_out:
    id: 'mash_out'
    section_name: 'MASH'
    description: 'uses the Min-hash method to identify composition of your sample. I may be helpful to visualize this as "percents" rather than hard counts.'
    plot_type: 'bargraph'
    format: 'tsv'
    pconfig:
      id: 'mash_abundance'
      title: "Mash Organisms"
      ylab: 'Abundance'
  metaphlan_out:
    id: 'metaphlan_out'
    section_name: 'METAPHLAN'
    description: 'characterizes the metagenomic composition of your sample. Values represent what percentage of your sample is from a given genus sp.'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'metaphlan_out'
      title: 'Sample Metagenomic Composition: Organism Percent Abundance'
      xTitle: 'Organism'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  abricate_out:
    id: 'abricate_out'
    section_name: 'ABRICATE: TOTALS'
    description: 'counts of AR genes identified by Abricate to asses which DB gives you the best resolution for AR profiling'
    plot_type: 'bargraph'
    format: 'tsv'
    pconfig:
      id: 'abricate_total_counts'
      title: 'Total AR Gene Counts from Abricate by Database'
      xTitle: 'Database'
  abricate_argannot:
    id: 'abricate_argannot'
    section_name: 'ABRICATE: ARGANNOT'
    description: 'Values represent % coverage of each gene. I.e. 100% means 100% of that gene is found in that genome'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'abricate_argannot'
      title: 'AR Gene Counts from Argannot DB'
      xTitle: 'Gene'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  abricate_ecoh:
    id: 'abricate_ecoh'
    section_name: 'ABRICATE: ECOH'
    description: 'Values represent % coverage of each gene. I.e. 100% means 100% of that gene is found in that genome'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'abricate_ecoh'
      title: 'AR Gene Counts from ECOH DB'
      xTitle: 'Gene'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  abricate_plasmidfinder:
    id: 'abricate_plasmidfinder'
    section_name: 'ABRICATE: PLASMIDFINDER'
    description: 'Values represent % coverage of each gene. I.e. 100% means 100% of that gene is found in that genome'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'abricate_plasmidfinder'
      title: 'AR Gene Counts from PlasmidFinder DB'
      xTitle: 'Gene'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  abricate_vfdb:
    id: 'abricate_vfdb'
    section_name: 'ABRICATE: VFDB'
    description: 'Values represent % coverage of each gene. I.e. 100% means 100% of that gene is found in that genome'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'abricate_vfdb'
      title: 'AR Gene Counts from VF DB'
      xTitle: 'Gene'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  abricate_card:
    id: 'abricate_card'
    section_name: 'ABRICATE: CARD'
    description: 'Values represent % coverage of each gene. I.e. 100% means 100% of that gene is found in that genome'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'abricate_card'
      title: 'AR Gene Counts from Card DB'
      xTitle: 'Gene'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  abricate_ecoli_vf:
    id: 'abricate_ecoli_vf'
    section_name: 'ABRICATE: EcoliVF'
    description: 'Values represent % coverage of each gene. I.e. 100% means 100% of that gene is found in that genome'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'abricate_ecoli_vf'
      title: 'AR Gene Counts from Ecoli VF DB'
      xTitle: 'Gene'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  abricate_ncbi:
    id: 'abricate_ncbi'
    section_name: 'ABRICATE: NCBI'
    description: 'Values represent % coverage of each gene. I.e. 100% means 100% of that gene is found in that genome'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'abricate_ncbi'
      title: 'AR Gene Counts from NCBI DB'
      xTitle: 'Gene'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  abricate_resfinder:
    id: 'abricate_resfinder'
    section_name: 'ABRICATE: RESFINDER'
    description: 'Values represent % coverage of each gene. I.e. 100% means 100% of that gene is found in that genome'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'abricate_resfinder'
      title: 'AR Gene Counts from Resfinder DB'
      xTitle: 'Gene'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  card_drugResistance:
    id: 'card_drugResistance'
    section_name: 'Drug Resistance: CARD'
    description: 'shows the number of times resistance to a particular drug was predicted by Abricate using the CARD database. i.e. a count of 5 means there were 5 genes that Abricate predics will confer AR resistance to that drug.'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'card_drugResistance'
      title: 'Putative Drug Resistance: CARD'
      xTitle: 'Drug'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  ncbi_drugResistance:
    id: 'ncbi_drugResistance'
    section_name: 'Drug Resistance: NCBI'
    description: 'shows the number of times resistance to a particular drug was predicted by Abricate using the NCBI database. i.e. a count of 5 means there were 5 genes that Abricate predics will confer AR resistance to that drug.'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'ncbi_drugResistance'
      title: 'Putative Drug Resistance: NCBI'
      xTitle: 'Drug'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]
  resfinder_drugResistance:
    id: 'resfinder_drugResistance'
    section_name: 'Drug Resistance: Resfinder'
    description: 'shows the number of times resistance to a particular drug was predicted by Abricate using the Resfinder database. i.e. a count of 5 means there were 5 genes that Abricate predics will confer AR resistance to that drug.'
    plot_type: 'heatmap'
    format: 'tsv'
    pconfig:
      id: 'resfinder_drugResistance'
      title: 'Putative Drug Resistance: Resfinder'
      xTitle: 'Drug'
      square: False
      sortRows: True
      colstops: [ [0, '#FFFFFF'], [0.33, '#d0d1e6'], [0.75, '#045a8d'], [1.0, '#000000'], ]

sp:
  organismal_specs:
    fn: 'organismal_specs.tsv'
  assembly_specs:
    fn: 'assembly_specs.tsv'
  warnings:
    fn: 'warning_flags_multiqc.tsv'
  mash_out:
    fn: 'mash_abundance.tsv'
  metaphlan_out:
    fn: 'metagenomic_composition_multiqc.tsv'
  abricate_out:
    fn: 'abricate_total_counts.tsv'
  abricate_argannot:
    fn: 'abricate_argannot_multiqc.tsv'
  abricate_ecoh:
    fn: 'abricate_ecoh_multiqc.tsv'
  abricate_plasmidfinder:
    fn: 'abricate_plasmidfinder_multiqc.tsv'
  abricate_vfdb:
    fn: 'abricate_vfdb_multiqc.tsv'
  abricate_card:
    fn: 'abricate_card_multiqc.tsv'
  abricate_ecoli_vf:
    fn: 'abricate_ecoli_vf_multiqc.tsv'
  abricate_ncbi:
    fn: 'abricate_ncbi_multiqc.tsv'
  abricate_resfinder:
    fn: 'abricate_resfinder_multiqc.tsv'
  card_drugResistance:
    fn: 'card_drug_resistance.tsv'
  ncbi_drugResistance:
    fn: 'ncbi_drug_resistance.tsv'
  resfinder_drugResistance:
    fn: 'resfinder_drug_resistance.tsv'""")
        mqc_config.close()

def run_command_logger(command):
    """ Runs a command and captures STDOUT to logfile as well as printing to screen."""

    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    file_full = open('logfile_full.txt', 'a')
    file_full.write(f'RUNNING COMMAND: {command}\n')
    file = open('logfile.txt', 'a')
    file.write(f'RUNNING COMMAND: {command}\n')
    sys.stdout.write(f'RUNNING COMMAND: {command}\n')
    for line in p.stdout:
        sys.stdout.write(line)
        file_full.write(line)
    p.wait()
    file_full.close()

def logger(string):
    """Prints string to screen, as well as appending it to logfile."""

    file = open('logfile.txt', 'a')
    file_full = open('logfile_full.txt', 'a')
    print(string)
    file.write(f'{string}\n')
    file_full.write(f'{string}\n')
    file.close()
    file_full.close()

def make_directory(dirname):
    """Makes a directory if it does not exist"""

    try:
        os.mkdir(dirname)
        logger(f'{dirname} folder created.')
    except FileExistsError:
        logger(f'{dirname} already exists. Will not overwrite')

def remove(item):
    """Removes file or empty directory"""

    if os.path.isfile(item):
        os.remove(item)
        logger(f'{item} has been removed')
    elif os.path.isdir(item):
        os.rmdir(item)
        logger(f'{item} has been removed')
    else:
        pass

def get_memory(d, p):
    """Determines maximum amount of RAM on system, and calculates 75% of that for program execution"""

    run = subprocess.Popen('free --giga', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    for line in run.stdout:
        if line.startswith('Mem'):
            RAM = line.split()[6]
        else:
            pass
    RAM = int(RAM)
    ram_use = round(RAM * 0.75)
    set_ram = int(p.ram)

    if set_ram == 0:                                                            # if RAM was unset, calculate to use 75% of available RAM
        p.ram = ram_use
    elif set_ram >= RAM:                                                        # check to confirm the RAM specified is not more than what is on the machine
        p.ram = ram_use
    logger(f'Total RAM available on machine is {RAM}, RAM was set to {set_ram}, this analysis will use {p.ram}.')

    run = subprocess.Popen('nproc', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    CPU = int(run.stdout.readlines()[0])
    cpu_use = round(CPU * 0.75)
    set_cpu = int(p.cores)

    if set_cpu == 0:                                                            # if CPU was unset, calculate to use 75% of available RAM
        p.cores = cpu_use
    elif set_cpu >= CPU:                                                        # check to confirm the CPU specified is not more than what is on the machine
        p.cores = cpu_use
    logger(f'Total CPU available on machine is {CPU}, CPU was set to {set_cpu}, this analysis will use {p.cores}.')

    run = subprocess.Popen('df -h', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    for line in run.stdout:
        if '/dev/' in line and 'da1' in line:
            space_avail = line.split()[3]                                       # Isolating the number from available memory. Using Sda1/XVda1 drive
            if space_avail[-1:] == 'M':
                space = float(space_avail[:-1]) * 1000
            if space_avail[-1:] == 'G':
                space = float(space_avail[:-1]) * 1000000
            if space_avail[-1:] == 'T':
                space = float(space_avail[:-1]) * 1000000000
        else:
            pass

    file_types = ('*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz')                   # Looking for any sequence files in CWD and will count them, divide by 2 to get number of datasets
    files_present = []
    for type in file_types:
        files_present.extend(glob.glob(os.path.join(d.base, type)))
    samp_num = len(files_present) / 2
    samp_num = round(samp_num)
    if os.path.isfile(f'{d.base}/SRR'):
        srr_f = open(f'{d.base}/SRR', 'r')
        for line in srr_f:
            if len(line.strip()) > 0:
                samp_num += 1
    needed_gb = samp_num * 4.5 + 3
    needed_bytes = needed_gb * 1000000
    if space < needed_bytes :
        logger(f'WARNING: Not enough space for this analysis. Base memory requirement is 3Gb, with an additional 4.5 for each isolate. \nYour analysis of {samp_num} samples requires {needed_gb}Gb and will likely error out. Increase storage allocation, or decrease number of samples')
    else:
        gb_avail = float(space/1000000)
        gb_avail = round(gb_avail, 2)
        logger(f'The base infrastructure for this analysis requires 3Gb of space. Each sample requires 4.5Gb. \nYou are attempting to run {samp_num} samples, requiring {needed_gb}Gb. \nYou have {gb_avail}Gb space currently available. Proceeding.')

def confirm_present(file):
    """Confirms file exists before program has run. If file DNE, will not throw error (so next step of pipeline can generate it). If file exists and is empty, will throw error."""

    if os.path.isfile(file) and os.stat(file).st_size != 0:
        return True
    elif os.path.isfile(file) and os.stat(file).st_size == 0:
        logger(f'File {file} is empty. Please delete and rerun pipeline.')
        quit()
    else:
        pass

def confirm_complete(file):
    """Confirms file exists after program has run. If not, throws error and stops analysis"""

    if os.path.isfile(file) and os.stat(file).st_size != 0:
        pass
    else:
        logger(f'Issue creating {file}. Please check logfile to diagnose issue. Quitting analysis')
        exit()

def image_check(location, image, image_n):
    """Confirm Singularity Image File is available"""

    if not os.path.isfile(f'{location}/{image_n}'):
        logger(f'{image_n} not found locally. Building.')
        logger('This may take a few minutes.')
        cmd = f'singularity build {location}/{image_n} docker://{image}'
        run_command_logger(cmd)
        if not os.path.isfile(f'{location}/{image_n}'):
            logger(f'Issue building singularity image {image_n}. Please confirm you have internet access and singularity is fully functional. Quitting analysis.')
            quit()
    else:
        logger(f'{image_n} found.')

def singularity_check(cfg):
    """Confirms presence of singularity, as well as a directory containing all singularity files"""

    logger('Confirming singularity is installed and locating or downloading required image files.')
    if which('singularity') is not None:
        program_nameVers_std = subprocess.run('singularity --version', shell=True, stdout=subprocess.PIPE, universal_newlines=True).stdout
        program_nameVers = program_nameVers_std.rstrip('\n')
        logger(f'{program_nameVers} exists in your PATH')
    else:
        logger('ERROR: Cannot locate singularity in your PATH. Confirm it exists.\nSee https://sylabs.io/docs/ for installation instructions.\n')
        exit()

    if not os.path.isdir(cfg.singularity_loc):
        logger('Specified directory containing singularity image files does not exist. Creating directory {cfg.singularity_loc}')
        make_directory(cfg.singularity_loc)
    logger('Found Singularity image file location. Confirming desired images are present.')
    image_check(cfg.singularity_loc, cfg.sratoolkit, cfg.sratoolkit_sif)
    image_check(cfg.singularity_loc, cfg.seqyclean, cfg.seqyclean_sif)
    image_check(cfg.singularity_loc, cfg.fastqc, cfg.fastqc_sif)
    image_check(cfg.singularity_loc, cfg.metaphlan2, cfg.metaphlan2_sif)
    image_check(cfg.singularity_loc, cfg.shovill, cfg.shovill_sif)
    image_check(cfg.singularity_loc, cfg.quast, cfg.quast_sif)
    image_check(cfg.singularity_loc, cfg.bwa, cfg.bwa_sif)
    image_check(cfg.singularity_loc, cfg.samtools, cfg.samtools_sif)
    image_check(cfg.singularity_loc, cfg.lyveset, cfg.lyveset_sif)
    image_check(cfg.singularity_loc, cfg.mash, cfg.mash_sif)
    image_check(cfg.singularity_loc, cfg.seqsero2, cfg.seqsero2_sif)
    image_check(cfg.singularity_loc, cfg.serotypefinder, cfg.serotypefinder_sif)
    image_check(cfg.singularity_loc, cfg.emmtyping, cfg.emmtyping_sif)
    image_check(cfg.singularity_loc, cfg.legsta, cfg.legsta_sif)
    image_check(cfg.singularity_loc, cfg.seroba, cfg.seroba_sif)
    image_check(cfg.singularity_loc, cfg.shigatyper, cfg.shigatyper_sif)
    image_check(cfg.singularity_loc, cfg.abricate, cfg.abricate_sif)
    image_check(cfg.singularity_loc, cfg.multiqc, cfg.multiqc_sif)
    logger('')

def database_check(cfg):
    """Confirms local install of databases, and will install if not present."""
    logger('Confirming install of MetaPhlAn database')
    if os.path.isdir(cfg.metaphlan2_loc):
        logger(f'Intended location for MetaPhlAn2 database: {cfg.metaphlan2_loc}')
    else:
        logger('No local install of MetaPhlAn2 database found. Building locally now')
        make_directory(cfg.metaphlan2_loc)
        cmd = f'singularity exec -B {cfg.metaphlan2_loc}:/opt/conda/bin/metaphlan_databases \
        {cfg.singularity_loc}/{cfg.metaphlan2_sif} metaphlan2.py --install'
        run_command_logger(cmd)

    # - CONFRIMING SUCCESSFUL INSTALL
    cmd = f'singularity exec -B {cfg.metaphlan2_loc}:/opt/conda/bin/metaphlan_databases \
    {cfg.singularity_loc}/{cfg.metaphlan2_sif} metaphlan2.py --install'
    run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    status = str(run.stdout.readlines()[0]).rstrip()
    if status == "The database is installed":
        logger(f'{status}')
    else:
        logger(f'Issue with metaphlan database install in {cfg.metaphlan2_loc}. \nPlease remove the metaphlan directory ({cfg.metaphlan2_loc}) and rerun pipeline. \nQuitting.')
        quit()

class Directories:
    """Generates strings for all directory locations required in this analysis. NOTE- does NOT make the directories, this simply allows for easy future manipulation of I/O"""

    def __init__(self):
        self.base = ''

def gen_directories():

    dirs = Directories()
    dirs.base = os.getcwd()
    dirs.metrics = f'{dirs.base}/output_metrics'
    dirs.raw = f'{dirs.base}/raw_seqs'
    dirs.trim = f'{dirs.base}/trimmed'
    dirs.shuffle = f'{dirs.base}/shuffle'
    dirs.metaphlan2 = f'{dirs.base}/metaphlan2'
    dirs.fastqc = f'{dirs.base}/fastqc'
    dirs.shovill = f'{dirs.base}/shovill'
    dirs.bwa = f'{dirs.base}/bwa'
    dirs.quast = f'{dirs.base}/quast'
    dirs.mash = f'{dirs.base}/mash'
    dirs.serotypefinder = f'{dirs.base}/serotypefinder'
    dirs.seqsero2 = f'{dirs.base}/sero_seqsero2'
    dirs.emm_sero = f'{dirs.base}/sero_emm'
    dirs.leg_sero = f'{dirs.base}/sero_leg'
    dirs.seroba = f'{dirs.base}/seroba'
    dirs.shigatyper = f'{dirs.base}/sero_shigatyper'
    dirs.abricate = f'{dirs.base}/abricate'
    dirs.for_multiqc = f'{dirs.base}/for_multiqc'
    dirs.repo = repo_dir = os.path.dirname(__file__)
    return dirs
