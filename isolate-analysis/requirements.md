# We use Conda to manage different packages. As such, if you do not have it- I highly suggest you get it!

## Download and install Anaconda/Miniconda
[need help deciding which to install?](https://unidata.github.io/online-python-training/choosing.html)
```
curl -O https://repo.anaconda.com/archive/Anaconda-latest-Linux-x86_64.sh
bash Anaconda-latest-Linux-x86_64.sh
```

## OR

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh
bash Miniconda3-latest-Linux-x86.sh
```

Now log in an log out to update your environment

Regardless of which you chose, you need to update conda
```
conda update conda
conda update anaconda
```

## Good to go. Now create your necessary environments

Two ways to go about this. Either do a fresh install, or install from the `yaml` files in this repo. Doing a fresh install may give you the latest-and-greatest, but sometimes new releases can "break" programs.

### ------ Fresh installs ------ 

#### Snakemake
```
conda create -n snakemake
conda activate snakemake
conda install -c bioconda -c conda-forge snakemake
```
Once complete
```
conda deactivate
```
* NOTE: the first time you activate a conda envt, you likely will need to "specify the language" being used. Choose `bash` 

#### Roary
```
conda create -n roary
conda activate roary
conda install -c bioconda roary
```
Once complete
```
conda deactivate
```

#### IQ-TREE
```
conda create -n iqtree
conda activate iqtree
conda install -c bioconda iqtree
```
Once complete
```
conda deactivate
```

### ------ Install from .yml file ------ 

```
conda env create -f Isolate_assembly_analysis/conda_env_ymls/snakemake.yml
conda env create -f Isolate_assembly_analysis/conda_env_ymls/roary.yml
conda env create -f Isolate_assembly_analysis/conda_env_ymls/iqtree.yml
```

All of your conda environments are now created. To activate them individually, simply call `conda activate <env_name>`
The `URF_pipeline.py script will automatically do this for you, so no need to worry about it for execution of this pipeline.

## Next, we need Singularity, as this is what the URF Pipeline uses in lieu of direct installs of each software package
The official documentation can be found [here](https://singularity.lbl.gov/install-linux) and [here](https://sylabs.io/guides/3.0/user-guide/installation.html).

On Centos
```
yum install singularity
```
On Ubuntu
```
apt-get install -y singularity-container
```
If you don't have root or sudo privileges, ask your administrator to install singularity.

## Finally, let's pull down the UPHL Reference Free (URF) Pipeline Repo
Documentation on this pipeline [here](https://github.com/StaPH-B/UPHL)

```
git clone https://github.com/StaPH-B/UPHL.git
cd UPHL
git init
```
** You can sync with the most-current version of UPHL by executing `git pull` from within this directory at any time. 

** NOTE: The most current URF has a bug in the `multiqc` step. I have pulled a version that is slightly older, but functional. To pull this verison:

MA-PHL does not own nor maintain this pipeline nor version, and will be reverting to the offical URF as soon as the multiqc step is repaired. 

## Very last step is to have a local BLAST DB install

** NOTE: You can run the URF_treebld.py script without a local DB, it will throw an error on the final MultiQC step, but all other processes will successfully complete. 

```
mkdir blast_nt_db
for i in {0..84}; do curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.$i.tar.gz; done
for i in $(ls nt*gz); do tar zxpvf $i; done
```
you now have a local install of the BLAST-nt db. You can also accomplish this by downloading the BLAST suite, and using their `update_blastdb.pl` script. More instructions [here](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
