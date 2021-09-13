#!/usr/bin/env python3

# brr_tree_building.py
# Andrew S. Lang
# Created: 20FEB2020
# Last Modified: 28FEB2020

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
    epilog='''Phylogenetic Tree Building. This pipeline is maintained by Andrew Lang
    at the Massachusetts State Public Health Lab. Additional information pertaining to pipeline usage can
    be found at https://gitlab.com/ma_ngs/brr_pipeline. For any issues with program operation, please
    contact Andrew at Andrew.Lang@massmail.state.ma.us for assistance.''')
#parser._action_groups.pop()             # allows me to list required arguments before optional
#required = parser.add_argument_group('Required Arguments')
#required.add_argument('-n', '--noAssemble', action='store_true', help='will not assemble genomes.')
optional = parser.add_argument_group('Optional Arugments')
optional.add_argument('-cfg', '--config', help='location of config file (Default= ./brr_tree_config.ini)', default='./brr_tree_config.ini')

# NOTE: YOUR READ FILES SHOULD BE ZIPPED.

class Directories:
    """Generates strings for all directory locations required in this analysis. NOTE- does NOT make the directories, this simply allows for easy future manipulation of I/O"""

    def __init__(self):
        self.base = ''

def gen_directories():

    dirs = Directories()
    dirs.base = os.getcwd()
    dirs.raw = f'{dirs.base}/raw_seqs'
    dirs.roary = f'{dirs.base}/roary'
    dirs.assemblies = f'{dirs.base}/assemblies'
    dirs.gffs = f'{dirs.base}/gffs'
    dirs.lyveset = f'{dirs.base}/lyveset'
    dirs.repo = repo_dir = os.path.dirname(__file__)
    return dirs

def get_sample_list(dir):
    """Generates a list of samples in provided direcotry. Will remove '.fa' or '.fasta' extensions."""
    # GETS LIST OF SAMPLES FOR FUTURE ITERATIVE FUNCTIONS
    sample_list = []
    for file in os.listdir(dir):
        if file.endswith('.fasta'):
            base = file[:-6]                                                   # removing ".fasta" suffix
            sample_list.append(base)
        elif file.endswith('.fa'):
            base = file[:-3]
            sample_list.append(base)
    sample_list.sort()
    return sample_list

def hamming(seq1, seq2):
    """Find the hamming distance (SNP differences) between core genomes in aln file. Will only process files that
    actually contain data, and requires that the sequences be of the same length. Roary ouputs satisfy these
    specifications if completed successfully."""

    while len(seq1) > 0:
        if len(seq1) != len(seq2):
            print('Core genome sequence lengths are not equal. \nCheck that your alignment completed successfully.')
            exit()
        elif len(seq1) == len(seq2):
#            return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))
            return sum(ch1 != ch2 and ch1 != '-' and ch2 != '-' for ch1, ch2 in zip(seq1, seq2)) #this will not count the "-" in alignments as SNPs

def run_prokka(d, p):
    """Runs Prokka to annotate genome assemblies"""
    prokka --cpu 1 --compliant --centre --URF --mincontiglen 500 --outdir Prokka/sample --locustag locus_tag --prefix sample --genus ${mash_result[0]} --species ${mash_result[1]} --force shovill_result/sample/contigs.fa

def run_roary(d, p):
    """Runs Roary to generate core genome file"""

    roary_cmd = f'roary -p {cpu} -f {d.roary} -e -n {d.treebld}/*gff'

def generate_SNP_MATRIX(dir):
    """UNTESTED: Generates SNP Matrix from Core Genome Alignment""")

    if os.path.isfile(f'{d.treebld}/pairwise_SNP_matrix.tsv'):
        print(f'SNP matrix already exists in {d.treebld}. \nPlease delete or rename previous matrix if you would like to generate a new one.\n')
        exit()
    else:
        pass

    # PARSE FASTA TO DETERMINE HAMMING DISTANCE
    print('Generating SNP Matrix')
    mfasta = open(roary_out_f, 'r')
    fasta_d_ord = OrderedDict()
    fasta_l = []
    hamming_d = {}

    print('Parsing sequence data...')
    for record in SeqIO.parse(mfasta, 'fasta'):                                 # Generates dictionary (isolate_ID: Sequence)
        seq = str(record.seq.upper())
        fasta_d_ord[record.id] = seq

    for id1 in fasta_d_ord.keys():                                              # Iterate over each entry and compare to all others
        for id2 in fasta_d_ord.keys():                                          # Basically doing double work of counting here, but
            seq1 = fasta_d_ord[id1]                                             # it makes the matrix generation much easier
            seq2 = fasta_d_ord[id2]
            ham_dist = hamming(seq1, seq2)
            hamming_d[id1 + ':' + id2] = str(ham_dist)

    print('..........')
    # GENERATE PAIRWISWE SNP MATRIX AND AVERAGE SNP DIFFERENCE OF EACH ISOLATE
    mtrx_fname = f'{d.treebld}/pairwise_SNP_matrix.tsv'
    mtrx_f = open(mtrx_fname, 'w')
    mtrx_out_str = '.'                                                          # Placeholder for formatting of the first line of SNP Matrix
    for id1 in fasta_d_ord.keys():
        mtrx_out_str += f'\t{id1}'
    mtrx_out_str += '\n'

    avg_snp_d = {}
    num_isolates = len(fasta_d_ord.keys())
    for id1 in fasta_d_ord.keys():
        line = id1
        snp_tot = 0
        for id2 in fasta_d_ord.keys():
            snp_count = hamming_d[id1 + ':' + id2]
            line += f'\t{snp_count}'
            snp_tot = snp_tot + int(snp_count)
        mtrx_out_str += f'{line}\n'
        avg_snp = str(round(snp_tot/num_isolates, 1))                           # Rounds to one decimal place
        avg_snp_d[id1] = avg_snp                                                # dictionary of average snp counts (isolateID: avgSNP)
    mtrx_out_str.rstrip('\n')
    mtrx_f.write(mtrx_out_str)

    print('..........')
    avg_ct_fname = f'{d.treebld}/average_SNP_counts.tsv'
    avg_ct_f = open(avg_ct_fname, 'w')
    sorted_l = sorted(avg_snp_d.items(), key=lambda x: x[1])
    avg_out_str = ''
    for isolate in sorted_l:
        iso_nm = isolate[0]
        avg_ct = isolate[1]
        avg_out_str += f'{iso_nm}\t{avg_ct}\n'
    avg_ct_f.write(avg_out_str)
    avg_ct_f.close()
    print('Complete.')

isolate best reference


def run_lyveset(d, cpu):
    """Run Lyve-SET from container"""

    if os.path.isdir(d.lyveset):
        print('You have already run Lyve-SET. Please rename "lyveset" directory, and re-enter command')
        exit()
    else:
        if os.path.isfile(f'{d.base}/samples_lyveset.txt'):
            print('Found samples list to include in this run. Creating lyveset directory structure.')
            subprocess.run(f'singularity exec {singularity_loc}/lyveset.sif set_manage.pl --create {d.lyveset}', shell=True)
            shutil.move(f'{d.base}/samples_lyveset.txt', f'{d.lyveset}/samples_lyveset.txt')
            samples = open(f'{d.lyveset}/samples_lyveset.txt', 'r')
            print('')
            for line in samples:
                if line.startswith('reference'):
                    ref = line.strip().split('\t')[1]
                    print(f'Using {ref} as reference for this Lyve-SET run.')
                elif line.startswith('Y'):
                    isolate = line.rstrip().split('\t')[1]
                    print(f'{isolate} will be included in this analysis')
                else:
                    pass
            samples.close()

            try:
                os.mkdir(f'{d.shuffle}')
            except FileExistsError:
                pass

            print('')
            samples = open(f'{d.lyveset}/samples_lyveset.txt', 'r')
            for line in samples:
                if line.startswith('reference'):
                    ref_seq_id = line.rstrip().split('\t')[1]
                    ref_seq = os.path.basename(glob.glob(f'{d.base}/ALL_assembled/{ref_seq_id}*')[0])        # This extracts just the filename from the full path (glob returns list, so I'm grabbing the first, and only, element of that list)
                    os.symlink(f'{d.base}/ALL_assembled/{ref_seq}', f'{d.lyveset}/reference/{ref_seq}')
                    if os.path.islink(f'{d.lyveset}/reference/{ref_seq}'):
                        print(f'Found reference in Lyve-SET directory structure.')
                elif line.startswith('Y'):
                    isolate = line.rstrip().split('\t')[1]
                    isolate_r1 = glob.glob(f'{d.raw}/{isolate}*_R1_*')[0]                    # Pulling R1 and R2 filename
                    isolate_r2 = glob.glob(f'{d.raw}/{isolate}*_R2_*')[0]

                    if os.path.isfile(f'{d.lyveset}/reads/{isolate}.cleaned.fastq.gz'):                     #if shufflereads.pl and trimclean.pl have been run, move to next readset
                        print(f'{isolate}.cleaned.fastq.gz already exists. Moving to next sample.')
                    else:
                        print(f'Running shufflereads.pl on {isolate}')
                        subprocess.run(f'singularity exec {singularity_loc}/lyveset.sif \
                            run_assembly_shuffleReads.pl {isolate_r1} {isolate_r2} > {d.shuffle}/{isolate}.shuffled.fastq', shell=True)
                        print('')

                        if os.path.isfile(f'{d.shuffle}/{isolate}.shuffled.fastq'):
                            print(f'Running trimclean.pl on {isolate}')
                            subprocess.run(f'singularity exec {singularity_loc}/lyveset.sif \
                                run_assembly_trimClean.pl --numcpus {cpu} -i {d.shuffle}/{isolate}.shuffled.fastq \
                                -o {d.lyveset}/reads/{isolate}.cleaned.fastq.gz --nosingletons', shell=True)
                            os.remove(f'{d.shuffle}/{isolate}.shuffled.fastq')
                            print(f'Removed {isolate}.fastq')
                            print('')
                else:
                    pass
            samples.close()
            os.rmdir(d.shuffle)

            print('Now running launch_set.pl')
            subprocess.run(f'singularity exec {singularity_loc}/lyveset.sif \
            launch_set.pl {d.lyveset} -ref {d.lyveset}/reference/{ref_seq} \
            --allowedFlanking --min_alt_frac 0.95 --min_coverage 20 --numcpus {cpu}', shell=True)

output SNP matrix and tree (look to abigails script here)

def run_pipeline(dirs, cfg, cfg_f, VERSION):
    """This function runs the main tree builing pipeline"""

    if os.path.isfile(f'{dirs.base}/logfile.txt'):
        print('Logfiles exist from previous run. Preappending "prev_" to those and generating new logfile\n')
        os.rename(f'{dirs.base}/logfile.txt', f'{dirs.base}/prev_logfile.txt')
        os.rename(f'{dirs.base}/logfile_full.txt', f'{dirs.base}/prev_logfile_full.txt')
    else:
        pass

    file_full = open('logfile_full.txt', 'a')
    logger(f'''\n****************************************************************************************************
You are running Phylogenetic Tree Building {VERSION}. This pipeline is maintained by Andrew Lang
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
#    if os.path.isfile(f'{dirs.base}/BRR_multiqc_config.yaml'):
#        logger('Found MultiQC yaml config file already present in current directory. Will use for MultiQC metric output.\n')
#    else:
#        brr_foundation.generate_multiqcConfig(dirs)
#        logger(f'No .yaml config file for MultiQC found in {dirs.base}. This file has now been created.\n')

    file = open('logfile.txt', 'w')
    brr_foundation.get_memory(dirs, cfg)
    brr_foundation.singularity_check(cfg)
    brr_foundation.database_check(cfg)
    make_directory(dirs.metrics)
    dirs.sample_list = get_sample_list(dirs.base)
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
           VERSION = '28FEB2020'
           run_pipeline(dirs, cfg, config_f, VERSION)
       except Exception:
           traceback.print_exc()
           exit()
