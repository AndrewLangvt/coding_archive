#!/usr/bin/env python3

# isolate_analysis.py
# Andrew S. Lang
# Created: 02OC2019
# Last Modified: 12NOV2019

import sys
import os
import argparse
import shutil
import subprocess
import time
import glob
from Bio import SeqIO
from collections import OrderedDict

# ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ -

# This section of code is the only area that should be adjusted by the user,
# unless you know what you're doing. The names of conda environments should
# not need changing if you established your work environment as outlined in
# the github page. Software and database locations may be adjusted to fit
# your system.

# Names of Conda enviornments
snkmk_env = 'reference_free'    # Snakemake env
roary_env = 'roary'             # Roary evn
iqtree_env = 'iqtree'           # IQ-TREE env

# Locations of Software & Databases
URF_loc = '/home/workflows/UPHL_functional/UPHL_reference_free_docker.smk'    # URF smk file
blast_nt_loc = '/home/workflows/blast_dbs/BLAST_nt'                           # BLAST nt database
singularity_loc = '/home/workflows/singularity_files'                         # Singularity files
figtree_loc = '/home/workflows/programs/FigTree_v1.4.3/lib/figtree.jar'       # Figtree jar file

# ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ - ~ -

parser = argparse.ArgumentParser(description='Run UPHL Reference-Free Pipeline and Tree Generation')
parser = argparse.ArgumentParser(
    epilog='''This script will run the UPHL Reference-Free Pipeline, and can execute
the follow-on tasks of comparing assembled & annotated genomes via Roary & IQ-TREE or
running Lyve-SET to conduct a high-quality SNP analysis.''')
parser._action_groups.pop()
required = parser.add_argument_group('Required Arguments')
#required.add_argument('-d', '--directory', help='Directory name containing sequencing files NOTE: You must execute this script from the directory above.', required=True)
required.add_argument('-a', '--analysis', help='Dictates what analysis to run. "ref-free" will execute the reference-free pipeline. "tree-build" will run Roary & IQ-TREE. "lyveset" will run Lyve-SET. "lyveset-prep" generates the sample_list.txt file needed for Lyve-SET (automatically generated at the end of tree-building stage)', required=True)
optional = parser.add_argument_group('Optional Arugments')
optional.add_argument('-c', '--cores', help='processing cores to use (Default = 32)', default='32')
optional.add_argument('-l', '--lyveset_dir', help='name for lyve-SET directory (Default = lyveset)', default='lyveset')
optional.add_argument('-t', '--treebuilding_dir', help='which treebld directory to use (Default = tree_building)', default='tree_building')


args = parser.parse_args()

class Directories:
    """Generates strings for all directory locations required in this analysis. NOTE- does NOT make the directories, this simply allows for easy future manipulation of I/O"""

    def __init__(self):
        self.base = ''
        self.raw = ''
        self.treebld = ''
        self.roary = ''
        self.iqtree = ''
        self.lyveset = ''
        self.shuffle = ''

def gen_directories():

    dirs = Directories()
    base_dir = os.getcwd()
    dirs.base = base_dir
    dirs.seq = f'{dirs.base}/Sequencing_reads'
    dirs.raw = f'{dirs.base}/Sequencing_reads/Raw'
    dirs.treebld = f'{dirs.base}/tree_building'
    dirs.roary = f'{dirs.treebld}/roary_output'
    dirs.iqtree = f'{dirs.treebld}/iqtree_out'
    dirs.shuffle = f'{dirs.base}/shuffle'
    return dirs

def create_dirs_and_moveFs(d, cpu):
    """Create directories for URF pipeline, and move reads into the necessary folder"""

    seq_files = os.listdir(d.base)
    seq_file_types = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']

    if os.path.isdir(d.raw) == True:                                            # If directories already exist, move new read data into folder
        print('''
You have already run URF Pipeline on this directory.
URF will assemble the new sequencing files you have added.
Now moving newly-added files to Sequencing_reads/Raw.''')
        for file in seq_files:                                                  # Moving sequencing files to the Sequencing_reads/Raw directory
            if file.endswith(tuple(seq_file_types)):
                try:
                    shutil.move(os.path.join(d.base,file), os.path.join(d.raw,file))
                    print(f'Moving {file} to {d.raw}')
                    time.sleep(0.5)
                except BaseException:
                    print(f'\n***ERROR: Issue moving {file} to {d.raw}. Quitting this analysis.')
                    exit()
    else:
        try:
            os.mkdir(d.seq)                                                     # Create /Sequencing_reads subdirectory in base dir
        except FileExistsError:
            print('Subdirectory "Sequencing_reads" already exists.')
        except OSError:
            print(f'\n***ERROR: Failed to create "Sequencing_reads" subdirectories in {d.base}')
            exit()
        else:
            print(f'Created "Sequencing_reads" subdirectory in {d.base}')
            time.sleep(1)

        try:
            os.mkdir(d.raw)                                                     # Create /Sequencing_reads/Raw subdirectory in base dir
        except FileExistsError:
            print('Subdirectory "Raw" already exists.')
        except OSError:
            print(f'\n***ERR)R: Failed to create "Raw" subdirectories in {d.base}.')
            exit()
        else:
            print(f'Created "Raw" subdirectory in {d.seq}.')
            time.sleep(1)

        if os.path.isdir(raw) == True:                                          # Testing to see if directories exist
            print(f'Successfully created subdirectories in {d.base}.\n')
            time.sleep(1)
            for file in seq_files:                                              # Moving sequencing files to the Sequencing_reads/Raw directory
                if file.endswith(tuple(seq_file_types)):
                    print(f'Moving {file} to {d.raw}')
                    try:
                        shutil.move(os.path.join(d.base,file), os.path.join(d.raw,file))
                        time.sleep(0.5)
                    except BaseException:
                        print(f'Issue moving {file} to {d.raw}. Quitting this analysis.')
                        exit()
        else:
            print('\n***ERROR:Issue with directory structure. Confirm you have write privileges for the provided directory.\n')
            exit()

def gen_sampleList(d):
    """Generates a config file for lyveset runs"""

    if os.path.isfile(f'{d.base}/samples_list.txt'):
        print(f'samples_for_lyveset.txt already exists in {d.base}. Moving old txt file to "oldList_samples_for_lyveset.txt" and generating new one.')
        shutil.move(f'{d.base}/samples_list.txt', f'{d.base}/prev_samples_list.txt')
    else:
        pass

    sample_file = open(f'{d.base}/samples_list.txt', 'w')
    avg_snps = open(f'{d.treebld}/average_SNP_counts.tsv', 'r')
    reference = avg_snps.readline().split('\t')[0]
    avg_snps.close()

    sample_file.write(f'reference\t{reference}\n')
    for sample_seq in os.listdir(d.raw):
        if sample_seq.endswith('L001_R1_001.fastq.gz'):
            sample = sample_seq.rstrip('L001_R1_001.fastq.gz')
            sample_file.write(f'Y\t{sample}\n')
    sample_file.close()

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

def run_pipeline(base, cpu):
    """Run the URF pipeline"""

    URF_cmd = f'snakemake --snakefile {URF_loc} --directory {d.base} --use-singularity --singularity-args "--bind {d.base}:/data,{blast_nt_loc}:/blast/blastdb" --singularity-prefix {singularity_loc} --cores {cpu}'

    print('''\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
You have selected to run the UPHL Reference-free (URF) Pipeline.
Attempting to build directory structure.''')

    create_dirs_and_moveFs(d, cpu)
    print(' ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n')
    print(f'Will now run the following command:\n {URF_cmd}')
    print(' ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n')
    time.sleep(3)
    subprocess.run(f'source activate {snkmk_env} && {URF_cmd} && conda deactivate', shell = True)

    print('\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n')
    print('Creating "tree_building" directory for follow-on analyses.')
    try:
        os.mkdir(d.treebld)                                                     # Create /tree_building subdirectory in base dir
    except FileExistsError:
        print('Subdirectory "tree_building" already exists.')

def run_tree_build(d, cpu):
    """RUN ROARY ON GFFs FROM URF PIPELINE.
    THEN GENERATE SNP MATRIX FROM CORE GENOME ALIGNMENT.
    FINALLY, BUILD TREE FROM ALIGNMENT AND VISUALIZE"""

    iqtree_pref = 'iqtree'

    #Expected output files form roary and iqtree
    roary_out_f = f'{d.roary}/core_gene_alignment.aln'
    iqtree_out_f = f'{d.iqtree}/{iqtree_pref}.ckp.gz'

    #Commands issued - Note - **DO NOT** change I/O names here, do so above in "Naming" block
    roary_cmd = f'roary -p {cpu} -f {d.roary} -e -n {d.treebld}/*gff'
    iqtree_cmd = f'iqtree -s {roary_out_f} -t RANDOM -m HKY+I+R -bb 1000 -nt {cpu} -pre {d.iqtree}/{iqtree_pref}'

    print('''\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
You have selected to run the tree-build portion of this wrapper.
This runs your assembled genomes through Roary and IQ-TREE.
Please ensure the GFFs you wish to compare are located in the
"tree_buildling" folder within the folder containing your URF run.

For example, if I ran URF on /pipe_run I would place the
GFFs I want to analyze in a folder called "tree_building" within "pipe_run".
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n''')
    time.sleep(3)

    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    # DIRECTORY CHECK
    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    print(f'''
~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

DIRECTORY STRUCTURE

''')
    print('Checking directory structure and file location for tree buliding.')
    time.sleep(0.5)

    if os.path.isdir(d.treebld) == True:
        if not os.listdir(d.treebld):
            print('\n***ERROR:"tree_building" is empty. Move GFFs to this location for analysis. Quitting.\n')
            exit()
        else:
            print('"tree_building" exists and contains files. Proceeding.\n')
    else:
         print(f'***ERROR:"tree_building" directory does not exist in your specified folder, {d.base}.\n')
         exit()

    time.sleep(1)

    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    # RUNNING ROARY
    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    print(f'''
~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

ROARY ALIGNMENT

''')
    if os.path.isfile(f'{roary_out_f}'):                                        # If roary has already run, skip and exit
        print('Roary has already been run on this directory.')
        print('Move Roary output folder from "tree_buildling" and re-issue command.')
        pass
    else:
        print(f'Directory structure check complete. Will now run the following command:\n{roary_cmd}\n')
        time.sleep(2)
        subprocess.run(f'source activate {roary_env} && {roary_cmd} && conda deactivate', shell = True)

    if os.path.isfile(f'{roary_out_f}'):                                        # If roary output exists, run iqtree
        try:
            os.mkdir(d.iqtree)                                                  # Create iqtree dir, but not if it exists
        except FileExistsError:
            pass
        print('Core gene alignment file exists from Roary run.')

        print(f'''
~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

IQTREE BUILDING

''')
        time.sleep(2)
        if os.path.isfile(f'{iqtree_out_f}'):
            print('You have already run IQ-TREE on this dataset.')
            print('Delete or move "iqtree_out" and re-issue command.\n')
            time.sleep(1)
            pass
        else:
            print(f'Now running IQ-TREE as follows:{iqtree_cmd}')
            subprocess.run(f'source activate {iqtree_env} && {iqtree_cmd} && conda deactivate', shell=True)
    else:
        print('\n***ERROR: core gene alignment file not present. Roary may have failed.\n')
        exit()

    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    # BUILDLING SNP MATRIX
    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    print(f'''
~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

SNP MATRIX

''')

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


    if os.path.isfile(f'{d.iqtree}/{iqtree_pref}.treefile'):
        print('\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
        print('Completion of IQ-TREE stage. Will now open figtree for viewing.\n')
        time.sleep(2)
        subprocess.run(f'java -jar {figtree_loc}', shell = True)
    else:
        print('\n***ERROR: IQ-TREE file not created, or is corrupted. Please consult log file.\n')
        exit()

    gen_sampleList(d)

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

if __name__ == '__main__':

    args = parser.parse_args()
    analysis_to_run = args.analysis
    cpu = args.cores
    lyveset_d = args.lyveset_dir.strip('/')

    d = gen_directories()
    d.lyveset = f'{d.base}/{lyveset_d}'

    if str(analysis_to_run) == 'ref-free':
        run_pipeline(d, cpu)
    elif str(analysis_to_run) == 'tree-build':
        run_tree_build(d, cpu)
        gen_sampleList(d)
    elif str(analysis_to_run) == 'lyveset-prep':
        d.treebld = args.treebuilding_dir
        gen_sampleList(d)
    elif str(analysis_to_run) == 'lyveset':
        run_lyveset(d, cpu)
    else:
        print('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ')
        print("What analysis do you want to run? I'm not a mind reader... yet... heh heh heh...")
        print('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n')
