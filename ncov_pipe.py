#!/usr/bin/env python3

import sys
import os
import json
import shutil
import glob
import datetime
import csv
import subprocess
import gzip
import shutil
from datetime import datetime as dt

config = '/home/alang/Desktop/COVID19/ncov_pipe_MA/genomic_analyses/sars_cov2_config.json'
cromwell_jar = '/home/alang/Desktop/COVID19/ncov_pipe_MA/genomic_analyses/cromwell-65.jar'
run_results_file_loc = "/home/alang/Desktop/COVID19/ncov_pipe_MA/genomic_analyses/cluster3_test/shared_databases/"
metadata_file_loc = "/home/alang/Desktop/COVID19/ncov_pipe_MA/genomic_analyses/cluster3_test/shared_databases/"
vx_consensus_loc = "/home/shared_databases/COVID19/consensus_genomes"
masphl_consensus_loc = "/home/alang/Desktop/COVID19/sequencing_batches/all_consensus_deidentified"
workflow_WDL = "/home/alang/Desktop/COVID19/ncov_pipe_MA/genomic_analyses/workflows/wf_viral_pipeline_local.wdl"

class SampleData:
    """Defines all aspects of sample data as a class, enabling passing to any function."""

    def __init__(self):
        self.samplename = ''
        self.read1 = ''
        self.read2 = ''
        self.PatientFirstName = ''
        self.PatientLastName = ''
        self.PatientDateOfBirth = ''
        self.PatientAddress = ''
        self.City = ''
        self.State = ''
        self.PostalCode = ''
        self.County = ''
        self.Provider = ''
        self.ProviderAddress = ''
        self.Source = ''
        self.SourceDescription = ''
        self.SpecimenType = ''
        self.CollectionDate = ''
        self.RecievedDate = ''
        self.TestName = ''
        self.Result = ''
        self.GISAID_ID = ''
        self.Sequencing_Reason = ''
        self.Sequencing_worksheet__ = ''
        self.collection = ''
        self.deidentified = ''
        self.iso_source = ''
        self.iso_state = ''
        self.iso_country = 'USA'
        self.iso_continent = 'North America'

class Samplestruct:
    """Defines all aspects of sample data as a class, enabling passing to any function."""

    def __init__(self):
        self.samplename = ''
        self.read1 = ''
        self.read2 = ''
        self.iso_state = ''
        self.iso_country = ''
        self.iso_continent = ''
        self.collection = ''
        self.deidentified = ''
        self.iso_source = ''

def abbrev_to_state(abbreviation):
        states = {'AL': 'Alabama',
                        'AK': 'Alaska',
                        'AS': 'American Samoa',
                        'AZ': 'Arizona',
                        'AR': 'Arkansas',
                        'CA': 'California',
                        'CO': 'Colorado',
                        'CT': 'Connecticut',
                        'DE': 'Delaware',
                        'DC': 'District of Columbia',
                        'FL': 'Florida',
                        'GA': 'Georgia',
                        'GU': 'Guam',
                        'HI': 'Hawaii',
                        'ID': 'Idaho',
                        'IL': 'Illinois',
                        'IN': 'Indiana',
                        'IA': 'Iowa',
                        'KS': 'Kansas',
                        'KY': 'Kentucky',
                        'LA': 'Louisiana',
                        'ME': 'Maine',
                        'MD': 'Maryland',
                        'MA': 'Massachusetts',
                        'MI': 'Michigan',
                        'MN': 'Minnesota',
                        'MS': 'Mississippi',
                        'MO': 'Missouri',
                        'MT': 'Montana',
                        'NE': 'Nebraska',
                        'NV': 'Nevada',
                        'NH': 'New Hampshire',
                        'NJ': 'New Jersey',
                        'NM': 'New Mexico',
                        'NY': 'New York',
                        'NC': 'North Carolina',
                        'ND': 'North Dakota',
                        'MP': 'Northern Mariana Islands',
                        'OH': 'Ohio',
                        'OK': 'Oklahoma',
                        'OR': 'Oregon',
                        'PA': 'Pennsylvania',
                        'PR': 'Puerto Rico',
                        'RI': 'Rhode Island',
                        'SC': 'South Carolina',
                        'SD': 'South Dakota',
                        'TN': 'Tennessee',
                        'TX': 'Texas',
                        'UT': 'Utah',
                        'VT': 'Vermont',
                        'VI': 'Virgin Islands',
                        'VA': 'Virginia',
                        'WA': 'Washington',
                        'WV': 'West Virginia',
                        'WI': 'Wisconsin',
                        'WY': 'Wyoming'}
        if abbreviation in states.keys():
                return states[abbreviation]
        else:
                return abbreviation

def struct_construct(sampleclass):
        struct = Samplestruct()
        struct.samplename = sampleclass.samplename
        struct.read1 = sampleclass.read1
        struct.read2 = sampleclass.read2
        struct.collection = sampleclass.collection
        struct.deidentified = sampleclass.deidentified
        struct.iso_state = sampleclass.iso_state
        struct.iso_country = sampleclass.iso_country
        struct.iso_continent = sampleclass.iso_continent
        struct.iso_source = sampleclass.iso_source
        return(struct)

def date():
        stamp = str(datetime.datetime.now().strftime("%Y%m%d"))
        return(stamp)

def timestamp():
        stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        return(stamp)

def gz_size(fname):
    with gzip.open(fname, 'rb') as f:
        return f.seek(0, whence=2)

def run_command_logger(command):
    """ Runs a command and captures STDOUT to logfile as well as printing to screen."""

    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    file_full = open(f'{date()}_logfile_full.txt', 'a')
    file_full.write(f'RUNNING COMMAND: {command}\n')
    file = open(f'{date()}_logfile.txt', 'a')
    file.write(f'RUNNING COMMAND: {command}\n')
    sys.stdout.write(f'RUNNING COMMAND: {command}\n')
    for line in p.stdout:
        sys.stdout.write(line)
        file_full.write(line)
    p.wait()
    file_full.close()

def logger(string):
    """Prints string to screen, as well as appending it to logfile."""

    file = open(f'{date()}_logfile.txt', 'a')
    file_full = open(f'{date()}_logfile_full.txt', 'a')
    print(string)
    file.write(f'{string}\n')
    file_full.write(f'{string}\n')
    file.close()
    file_full.close()

def date_reformat(cur_date):
        while cur_date != '':
                reform_date = dt.strptime(cur_date, '%m/%d/%Y').strftime('%Y-%m-%d')
                return reform_date
        else:
                return cur_date

def make_directory(dirname):
    """Makes a directory if it does not exist"""

    try:
        os.mkdir(dirname)
        logger(f'{dirname} folder created.')
    except FileExistsError:
        pass

def viral_refbased_prep(inputfile):
        # Read in BtB NGS-MAVEN Report CSV file.

        #logger(f'\n{timestamp()}: VIRAL REFBASED CONSENSUS CALL PREPARATION\n')

        logger(f'\n{timestamp()}: VIRAL REFBASED CONSENSUS CALL PREPARATION\n')

        # Make "reads" directory and move all readfiles that contain data into that folder
        make_directory('raw_reads')
        sample_read_IDs = []
        for readfile in glob.glob(f'{os.getcwd()}/*fastq.gz'):
                if gz_size(readfile):                                           # returns 0 if file is empty
                        readfile_base = os.path.basename(readfile)
                        shutil.move(readfile, f'{os.getcwd()}/raw_reads/{readfile_base}')
                elif not gz_size(readfile):
                        readfile_base = os.path.basename(readfile)
                        logger(f'No read data in file {readfile_base}. Moving to "failed_reads"')
                        make_directory('failed_reads')
                        shutil.move(readfile, f'{os.getcwd()}/failed_reads/{readfile_base}')

        for readfile in glob.glob(f'{os.getcwd()}/raw_reads/*fastq.gz'):
                readfile_base = os.path.basename(readfile)
                readname = readfile_base.split('_')[0] # Caputre just the specimen ID, stripping _S*_L001_R1_001.fastq.gz
                if readname not in sample_read_IDs:
                        sample_read_IDs.append(readname)

        # Parse metadata file and assign attributes from each line
        with open(inputfile, 'r', encoding='utf-8-sig') as infile:
                sample_metadata = list(csv.DictReader(infile, delimiter=','))
                sample_dict = {}
                today = datetime.date.today()
                for sample_line in sample_metadata:
                        samplename = sample_line["SampleNo"]
                        if samplename not in sample_dict.keys():
                                sampledata = SampleData()
                                sampledata.samplename = sample_line["SampleNo"]
                                sampledata.PatientFirstName = sample_line["PatientFirstName"]
                                sampledata.PatientLastName = sample_line["PatientLastName"]
                                sampledata.PatientDateOfBirth = date_reformat(sample_line["PatientDateOfBirth"])
                                sampledata.PatientAddress = sample_line["PatientAddress"]
                                sampledata.City = sample_line["City"]
                                sampledata.State = sample_line["State"]
                                sampledata.PostalCode = sample_line["PostalCode"]
                                sampledata.County = sample_line["County"]
                                sampledata.Provider = sample_line["Provider"]
                                sampledata.ProviderAddress = sample_line["ProviderAddress"]
                                if sample_line["Source"] == 'NP/OP Swab':
                                        sampledata.Source = 'Nasopharyngeal/Orophanyngeal swab'
                                        sampledata.iso_source = 'Nasopharyngeal/Orophanyngeal swab'
                                else:
                                        sampledata.Source = sample_line["Source"]
                                        sampledata.iso_source = sample_line["Source"]
                                sampledata.SourceDescription = sample_line["SourceDescription"]
                                sampledata.SpecimenType = sample_line["SpecimenType"]
                                sampledata.CollectionDate = date_reformat(sample_line["CollectionDate"])
                                sampledata.RecievedDate = date_reformat(sample_line["RecievedDate"])
                                sampledata.TestName = sample_line["TestName"]
                                sampledata.Result = sample_line["Result"]
                                sampledata.GISAID_ID = sample_line["GISAID_ID"]
                                sampledata.Sequencing_Reason = sample_line["Sequencing_Reason"]
                                sampledata.Sequencing_worksheet__ = sample_line["Sequencing_worksheet__"]
                                if sampledata.PatientDateOfBirth != '':
                                        sampledata.age = today.year - int(sampledata.PatientDateOfBirth.split('-')[0]) - ((today.month, today.day) < (int(sampledata.PatientDateOfBirth.split('-')[1]), int(sampledata.PatientDateOfBirth.split('-')[2])))
                                else:
                                        sampledata.age = ''
                                # name-specific attributes for pipeline input.
                                sampledata.deidentified = sample_line["GISAID_ID"]
                                full_state = abbrev_to_state(sample_line["State"])
                                sampledata.iso_state = full_state
                                sampledata.collection = date_reformat(sample_line["CollectionDate"])
                                sample_dict[samplename] = sampledata
                        else:
                                logger(f'Duplicate entries for {samplename} in Metadata file.')
                                pass

        def VXsort(VX):
                '''Sorts VX IDs- I.e. puts 21VX20 before 21VX110'''
                if 'VX' in VX:
                        return VX.split('VX')[1]
                else:
                        return VX

        sample_read_IDs.sort(key=VXsort)
        for sampleID in sample_read_IDs:
                read1 = glob.glob(f'{os.getcwd()}/raw_reads/{sampleID}*R1*fastq.gz')[0]
                read2 = glob.glob(f'{os.getcwd()}/raw_reads/{sampleID}*R2*fastq.gz')[0]
                if os.path.isfile(read1) and os.path.isfile(read2):
                        if sampleID in sample_dict.keys():
                                # Appends read data for specimens with data in Metadata file
                                sampledata = sample_dict[sampleID]
                                sampledata.read1 = read1
                                sampledata.read2 = read2
                                sample_dict[sampleID] = sampledata
                                logger(f'Data & readfiles found for {sampleID}')
                        else:
                                # Captures instances where there is no Metadata, but a readfile exists.
                                sampledata = SampleData()
                                sampledata.samplename = sampleID
                                sampledata.deidentified = ''
                                sampledata.collection = ''
                                sampledata.read1 = read1
                                sampledata.read2 = read2
                                sample_dict[sampleID] = sampledata
                                logger(f'WARNING: No metadata entry for {sampleID}. If this is incorrect, please add to the metadata file.')
                else:
                        logger(f'WARNING: Could not find L and R read files for {sampleID}. Excluding from analsis.')
        for sample_meta in sample_dict.keys():
                if sample_meta not in sample_read_IDs:
                        logger(f'WARNING: Metadata file contains entry or {sample_meta} but cannot locate readfiles for this specimen.\n   {sample_meta} will be excluded from this analysis')

        # combine entries into single string, with formatting for WDL workflow input.
        outfile = open('inputs_viral_refbased_local.json', 'w')
        samples = []
        for sample in sample_read_IDs:
                sampleinfo = sample_dict[sample]                        # Only writing to JSON file for those speicmens that have read data
                if os.path.isfile(sampleinfo.read1) and os.path.isfile(sampleinfo.read2):
                        pipe_struct = struct_construct(sampleinfo)   # extracting information to generate a sctruct
                        samples.append(pipe_struct)

        outfile.write('{"viral_pipeline_local.inputsamples": ')
        outfile.write(json.dumps(samples, default=lambda x: x.__dict__, indent = 4))
        outfile.write(',\n')

        # concatenate config file at the end
        config_file = open(config, 'r')
        for line in config_file:
                outfile.write(line)

        outfile.write('}')
        infile.close()
        config_file.close()
        outfile.close()

def run_viral_refbased_wf(inputs, program):

        logger(f'\n{timestamp()}: VIRAL REFERENCE BASED ASSEMBLY\n')
        if program == 'cromwell':
                command = f'java -jar {cromwell_jar} run -i {inputs} {workflow_WDL}'
                run_command_logger(command)
        else:
                command = f'source activate miniwdl && miniwdl run -i {inputs} {workflow_WDL} && conda deactivate'
                run_command_logger(command)
        if len(glob.glob(f'{os.getcwd()}/*_viral_pipeline_local')) != 0: # if this returns a non-empty list
                wdl_output = os.path.basename(max(glob.glob(f'{os.getcwd()}/*_viral_pipeline_local')))    # Getting the folder name of most recent WDL execution
                pipe_start = wdl_output.replace('_viral_pipeline_local', '')
                if os.path.isfile(f'{wdl_output}/workflow.log'):
                        f = open(f'{wdl_output}/workflow.log')
                        lines = f.read().splitlines()
                        last_line = lines[-1]
                        status = last_line.split(' ')[-1]
                        if status == 'done':
                                logger(f'\n{timestamp()}: VIRAL REFERENCE BASED ASSEMBLY CLEANUP\n')
                                logger(f'\n\nSuccessful completion of Viral Refbased Assembly Workflow. Now copying files to {os.getcwd()}/pipeline_outputs')
                                if program == 'cromwell':
                                        pass
                                else:
                                        miniwdl_cleanup(f'{wdl_output}/outputs.json')
                        else:
                                logger('\nERROR: Viral Refbased Assembly Workflow did not complete successuflly\n~ ~ ~ Status is not "done" ~ ~ ~ \nPlease consult logfile & contact Andrew.S.Lang@mass.gov 6179836279 for troubleshooting.\n')
                else:
                        logger('\nERROR: Viral Refbased Assembly Workflow did not complete successuflly. \n~ ~ ~ No workflow.log file~ ~ ~ \nPlease consult logfile & contact Andrew.S.Lang@mass.gov 6179836279 for troubleshooting.\n')
        else:
                logger('\nERROR: Viral Refbased Assembly Workflow did not complete successuflly. \n~ ~ ~ No workflow folder~ ~ ~ \nPlease consult logfile & contact Andrew.S.Lang@mass.gov 6179836279 for troubleshooting.\n')

def miniwdl_cleanup(jsonfile):
        ''' Read in JSON file denoting location of MiniWDL outputs. Create foldernames based upon task output name, and copy files to location'''
        jsonf = open(jsonfile, 'r')
        data = json.load(jsonf)

        base_outputdir = 'pipeline_outputs'

        for task,taskout in data.items():
                # Make directory for output files if DNE
                task_dir = task.split('.')[-1]
                os.makedirs(f'{base_outputdir}/{task_dir}',exist_ok=True)

                if isinstance(taskout, str) and taskout is not None:                                            # workflow.task: "singleOUT"
                        pass
                        filename = os.path.basename(taskout)
                        if not os.path.isfile(f'{base_outputdir}/{task_dir}/{filename}'):
                                shutil.copyfile(taskout, f'{base_outputdir}/{task_dir}/{filename}')
                elif isinstance(taskout, list):                                                                 # workflow.task: ["firstOUT", "secondOUT"]
                        for output in taskout:
                                if isinstance(output, str) and output is not None:
                                        filename = os.path.basename(output)
                                        if not os.path.isfile(f'{base_outputdir}/{task_dir}/{filename}'):
                                                shutil.copyfile(output, f'{base_outputdir}/{task_dir}/{filename}')
                                elif isinstance(output, list):                                                  #  workflow.task: [["firstOUT-1", "secondOUT-1"], ["firstOUT-2", "secondOUT-2"]
                                        for single_out in output:
                                                if isinstance(single_out, str) and single_out is not None:
                                                        filename = os.path.basename(single_out)
                                                        if not os.path.isfile(f'{base_outputdir}/{task_dir}/{filename}'):
                                                                shutil.copyfile(single_out, f'{base_outputdir}/{task_dir}/{filename}')

        merged_metrics = data[f'viral_pipeline_local.merged_metrics']
        if not os.path.isfile(f'{os.getcwd()}/run_results.csv'):
                shutil.copyfile(merged_metrics, f'{os.getcwd()}/run_results.csv')
        else:
                shutil.move(f'{os.getcwd()}/run_results.csv', f'{os.getcwd()}/{timestamp()}.run_results.csv')
                shutil.copyfile(merged_metrics, f'{os.getcwd()}/run_results.csv')
        jsonf.close()

def cromwell_cleanup(jsonfile):
        cwd = os.getcwd()                                                               # create primary output directory
        base = f'{cwd}/covid'
        try:
                os.mkdir(base)
                print(f'{base} directory created')
        except FileExistsError:
                print(f'{base} directory already exists. Will not overwrite')

        with open(jsonfile, 'r') as jsonf:                              # parse JSON file
                data = json.load(jsonf)
                outputs_dict = data["outputs"]
                for task in outputs_dict.keys():
                        task_dir = task.split('.')[1]           # strips the nCoV19_pipeline. prefix from the JSON key
                        try:
                                os.mkdir(f'{base}/{task_dir}')
                                print(f'covid/{task_dir} directory created.')
                        except FileExistsError:
                                print(f'covid/{task_dir} already exists. Will not overwrite')

                        task_items = outputs_dict[task]
                        if type(task_items) == list:                # if multiple task outputs, iterate over
                                for item in task_items:
                                        if item is not None:
                                                try:
                                                        fname = os.path.basename(item)
                                                        shutil.copyfile(item, f'{base}/{task_dir}/{fname}')
                                                        print(f'Moved file {fname} to {task_dir}')
                                                except:
                                                        print(f'Could not copy {fname} to covid/{task_dir}. File still exists here: {item}')
                        elif type(task_items) == str:                   # if single task output, only action that singular item
                                item = task_items
                                try:
                                        fname = os.path.basename(item)
                                        shutil.copyfile(item, f'{base}/{task_dir}/{fname}')
                                        print(f'Moved file {fname} to {task_dir}')
                                except:
                                        print(f'Could not copy {fname} to covid/{task_dir}. File still exists here: {item}')

                        else:
                                print("empty")
                        print('\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n')

def metadata_generate(btb_file, batch_id):
        # Read in BtB NGS-MAVEN Report CSV file.

        # Parse metadata file and assign attributes from each line
        with open(btb_file, 'r', encoding='utf-8-sig') as infile:
                sample_metadata = list(csv.DictReader(infile, delimiter=','))
                sample_dict = {}
                today = datetime.date.today()
                for sample_line in sample_metadata:
                        samplename = sample_line["SampleNo"]
                        if samplename not in sample_dict.keys():
                                sampledata = SampleData()
                                sampledata.samplename = sample_line["SampleNo"]
                                sampledata.PatientFirstName = sample_line["PatientFirstName"]
                                sampledata.PatientLastName = sample_line["PatientLastName"]
                                sampledata.PatientDateOfBirth = date_reformat(sample_line["PatientDateOfBirth"])
                                sampledata.PatientAddress = sample_line["PatientAddress"]
                                sampledata.City = sample_line["City"]
                                sampledata.State = sample_line["State"]
                                sampledata.PostalCode = sample_line["PostalCode"]
                                sampledata.County = sample_line["County"]
                                sampledata.Provider = sample_line["Provider"]
                                sampledata.ProviderAddress = sample_line["ProviderAddress"]
                                if sample_line["Source"] == 'NP/OP Swab':
                                        sampledata.Source = 'Nasopharyngeal/Orophanyngeal swab'
                                        sampledata.iso_source = 'Nasopharyngeal/Orophanyngeal swab'
                                else:
                                        sampledata.Source = sample_line["Source"]
                                        sampledata.iso_source = sample_line["Source"]
                                sampledata.SourceDescription = sample_line["SourceDescription"]
                                sampledata.SpecimenType = sample_line["SpecimenType"]
                                sampledata.CollectionDate = date_reformat(sample_line["CollectionDate"])
                                sampledata.RecievedDate = date_reformat(sample_line["RecievedDate"])
                                sampledata.TestName = sample_line["TestName"]
                                sampledata.Result = sample_line["Result"]
                                sampledata.GISAID_ID = sample_line["GISAID_ID"]
                                sampledata.Sequencing_Reason = sample_line["Sequencing_Reason"]
                                sampledata.Sequencing_worksheet__ = sample_line["Sequencing_worksheet__"]
                                if sampledata.PatientDateOfBirth != '':
                                        sampledata.age = today.year - int(sampledata.PatientDateOfBirth.split('-')[0]) - ((today.month, today.day) < (int(sampledata.PatientDateOfBirth.split('-')[1]), int(sampledata.PatientDateOfBirth.split('-')[2])))
                                else:
                                        sampledata.age = ''
                                # name-specific attributes for pipeline input.
                                sampledata.deidentified = sample_line["GISAID_ID"]
                                sampledata.iso_state = sample_line["State"]
                                sampledata.iso_continent = 'USA'
                                sampledata.collection = date_reformat(sample_line["CollectionDate"])
                                sample_dict[samplename] = sampledata
                        else:
                                print(f'Duplicate entries for {samplename} in Metadata file.')
                                pass

        outfile = open(f'Metadata_{batch_id}.csv', 'w')
        outfile.write('Accession_Number,MASPHL ID,Batch ID,Collection_Date,Reason_For_Sequencing,Provider,PatientFirstName,PatientLastName,PatientDateOfBirth,Age,PatientAddress,City,State,PostalCode,County,Source,Result,Ct Value (N1),Ct Value (N2),Ct Value (RP)\n')
        for sampleID in sample_dict.keys():
                sampledata = sample_dict[sampleID]
                if sampledata.Sequencing_worksheet__ == batch_id:   # if sample is in the desired worksheet
                        outfile.write(f'{sampleID},{sampledata.GISAID_ID},{sampledata.Sequencing_worksheet__},{sampledata.CollectionDate},{sampledata.Sequencing_Reason},{sampledata.Provider},{sampledata.PatientFirstName},{sampledata.PatientLastName},{sampledata.PatientDateOfBirth},{sampledata.age},{sampledata.PatientAddress},{sampledata.City},{sampledata.State},{sampledata.PostalCode},{sampledata.County},{sampledata.Source},{sampledata.Result},,,\n')

                        print(f'Data found for {sampleID}')
                else:
                        pass

if __name__ == '__main__':

        task = sys.argv[1]

        if task == 'prep':
                if len(sys.argv) == 3:
                        print(len(sys.argv))
                        metadata = sys.argv[2]
#                       viral_wdl_prep(metadata)
                        viral_refbased_prep(metadata)
                else:
                        print('ERROR: Please provide input metadata file')
        elif task == 'analyze':
                if os.path.isfile(f'{os.getcwd()}/inputs_viral_refbased_local.json'):
                        if len(sys.argv) == 3:
                                if sys.argv[2] == 'cromwell':
                                        run_viral_wf(f'{os.getcwd()}/inputs_viral_refbased_local.json', 'cromwell')
                        else:
                                run_viral_refbased_wf(f'{os.getcwd()}/inputs_viral_refbased_local.json', 'miniwdl')
                else:
                        print('ERROR: inputs_viral_refbased.json does not exist. Please ensure you have run the "prep" stage first')
        elif task == 'metadata':
                btb_file = sys.argv[2]
                batch_id = sys.argv[3]
                metadata_generate(btb_file, batch_id)
        else:
                print('\nPlease provide a valid task to execute. Options are "prep" and "analyze"\n')
