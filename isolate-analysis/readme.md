# Isolate_assembly_analysis

This repo contains all scripting information pertaining to the dockerized UPHL Reference-Free (URF) pipeline, and subsequent analyses wrapped in a single python script.
This is designed to run on Prometheus, the linux machine at MA-PHL. With some minor adjusments, it can be made to run on any machine. To make these adjustments, change
the `Names of Conda Environments` and `Locations of Software & Databases` at the beginning of the script.

**DISCLAIMER: This workflow has not been validated and MA-PHL does not endorse it's use in a diagnostic setting in any way.**

### Importing data to your account on Prometheus

1. Place all `fastq`/`fq` (gzipped or not) in a single folder for analysis
2. Open terminal in that directory, or if executing from command line- `cd` into that directory.

# The workflows

The python wrapper script `isolate_analysis.py` will execute two workflows, depending on what follows the `-a` flag.

1. UPHL Reference-free pipeline (`-a ref-free`). More info [here](https://github.com/StaPH-B/UPHL).
2. Alignment and Tree building (`-a tree-build`) using [Roary](https://sanger-pathogens.github.io/Roary/) and [IQ-TREE](http://www.iqtree.org/doc/Quickstart#minimal-command-line-examples)
3. High-quality SNP analysis (`-a lyveset`) via Lyve-SET. More info [here](https://github.com/lskatz/lyve-SET)

## Usage
All isolates must first be run through the URF pipeline (`ref-free`). Once complete, this script will allow you to diverge into either `tree-build` or `lyveset`.
  - To carry out `tree-build` analysis, copy `gff`s from `ALL_gff` directory from URF run into `tree_buidling` directory and execute `tree-build` step.
  - To carry out `lyveset` analysis, alter the `samples_list.txt` to reflect which samples you want to include (change `Y` to `N` to exclude), and save it as `samples_lyveset.txt`. `REF` is automatically determined from the "least different" assembly generated in `ref-free` analysis step.

### Code examples
```
isolate_analysis.py -a <workflow>

Required Arguments:
  -a ANALYSIS, --analysis ANALYSIS
                        Dictates what analysis to run. "ref-free" will execute
                        the reference-free pipeline. "tree-build" will run
                        Roary & IQ-TREE. "lyveset" will run Lyve-SET.

Optional Arugments:
  -c CORES, --cores CORES
                        processing cores to use (Default = 32)
  -l LYVESET_DIR, --lyveset_dir LYVESET_DIR
                        name for lyve-SET directory (Default = lyveset)
  -t TREEBUILDING_DIR, --treebuilding_dir TREEBUILDING_DIR
                        which treebld directory to use (Default =
                        tree_building)

This script will run the UPHL Reference-Free Pipeline, and can execute the
follow-on tasks of comparing assembled & annotated genomes via Roary & IQ-
TREE.

```

### URF Pipeline

By selecting this pipeline, all of your isolates will be assembled, and you can view a `multiqc_report.html` that will visualize metrics of your raw and processed data, as well as assemblies. This file is located in `logs`.

- **NOTE**: Generally, it has taken ~15 minutes per read dataset.

*An example command for the URF Pipeline:*

        isolate_analysis.py -a ref-free

### Tree Building

This workflow will take your selected `GFF`s and identify the pan-genome using Roary, from which IQ-TREE then builds a tree for visualization in figtree.

- **NOTE**: Once the URF Pipeline has completed, you must *COPY* `GFF` files for which you want to build a tree into the `tree_building` directory now present within your current working directory.

*An example command for the Tree Build workflow:*

        isolate_analysis.py -a tree-build

### Lyve-SET
This will use the "best" assembly (least different) identified in the URF analysis as a reference, and conduct a high-quality SNP analysis on all isolates with a `Y` in the `samples_list.txt` file.

- **NOTE**: this script will create the lyveset output in `lyveset` as a default. You can adjust this with a `-l` flag, (i.e `-l lyveset_alternate` would put all lyveset output into `lyveset_alternate`.)

*An example command for the Lyveset analysis:*

        isolate_analysis.py -a lyveset

If you wish to run Lyve-SET without running the tree-build step, you must first run `isolate_analysis.py -a lyveset-prep` to generate the `samples_list.txt` file. The `-t` flag dictates which of your Roary/IqTREE directories you'd like to reference to determine the "best" reference for Lyve-SET.

*An example command for Lyveset_prep:*

        isolate_analysis.py -a lyveset-prep -t tree_building_cluster01

---
### [*Software and Package Requirements*](https://gitlab.com/AndrewLang/isolate_assembly_analysis_/blob/master/Requirements.md)

1) Snakemake, roary, and iqtree conda environments
2) Singularity installed
3) Clone of URF pipeline repo
