***DISCLAIMER: This pipeline is not intended for any clinical or diagnostic purposes and currently is not CLIA validated.***

# MA-PHL Organismal Typing Pipeline
This pipeline is designed to take raw short-read sequenes (Paired-End), and generate assemblies to characterize Genus/species and serotype. The primary output file,
`compiled_metrics.tsv` contains the organismal identification, as well as assembly metrics, read quality information. All of the data in `compiled_metrics.tsv` is
easily visualized in the `multiqcReport.html`.

Requirements for this pipeline:
1. Singularity (see [HERE](https://sylabs.io/docs/) for installation instructions)
2. Python v3.6 or newer

## Typing Pipeline:
  1. Read Trimming ([SeqyClean](https://github.com/ibest/seqyclean))
  2. Read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
  3. Genome assembly & QC ([Shovill](https://github.com/tseemann/shovill) & [QUAST](http://quast.sourceforge.net/))
  4. Metagenomic contamination assessment ([MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2))
  5. Read alignment ([BWA](http://bio-bwa.sourceforge.net/bwa.shtml))
  6. Organism Identification ([MASH](https://mash.readthedocs.io/en/latest/))
  7. Serotyping (E. coli:[serotypefinder](https://bitbucket.org/genomicepidemiology/serotypefinder/src/master/), Salmonella:[SeqSero2](https://github.com/denglab/SeqSero2), Strep. pyogenes: [emm_typing](https://github.com/phe-bioinformatics/emm-typing-tool), Strep. pneumoniae: [seroba](https://github.com/sanger-pathogens/seroba), L. pneumophila: [legsta](https://github.com/tseemann/legsta))
  8. Antibiotic Resistance Characterization ([Abricate](https://github.com/tseemann/abricate))
  9. Compile and visualize metrics ([MultiQC](https://multiqc.info/))

This pipeline has a config file that enables the user to adjust many of the thresholds, computational resource limits, and parameters. A full list of commands can be found in the `Commands` section of the [config](https://gitlab.com/ma_ngs/brr_pipeline/blob/master/brr_config.ini).

Some of the adjustable parameters dictate what will "raise flags" in this pipeline. There is a `warning_flags.tsv` file that will be generated, indicating if your data passes or fails, given the metrics in your config. A few examples are:
```
GENOME SIZE - base pair size range of expected genome (for all organisms not in below "SpecificGenomeSizes" section). Default = 1800000:2200000
MAPPING RATE - percent reads mapped to assembled genome. Default = 99.50
AVERAGE COVERAGE - average read coverage of assembled genome. Default = 20.00
AVERAGE READ QUALITY - average quality of trimmed reads. Default = 30.00
SpecificGenomeSizes - These are expected genome sizes for specific organisms. The pipeline will identify the organism with MASH, and compare it to these
    Salmonella_enterica = 4400000:5200000
    Salmonella_bongori = 4200000:4600000
    Streptococcus_pyogenes = 1600000:2000000
    Streptococcus_pneumoniae = 1800000:2200000
    Escherichia_coli = 4400000:4800000
    Legionella_pneumophila = 3200000:3600000
    Campylobacter_jejuni = 1500000:1700000
    Klebsiella_pneumoniae = 5000000:5400000
```

## RUNNING THE PIPELINE
1. Place all reads into a single directory
2. If you would like to incorporate SRR files into your analysis, paste the Run Number (one per line) into a file called "SRR". See example below
```
SRR10294869
SRR10294873
SRR10375928
SRR10375927
```
3. Run `brr_maineffort.py`
    - if this is the first time running `brr_maineffort.py` in this directory, it will generate a new config file. Make sure the parameters for this are correct (primarily location of singularity files and metaphlan2 database)
    - if you already have the config file present (or have used the `-c` flag to indicate another location of the config), the pipeline will begin

## PIPELINE OUTPUTS
1. **compiled_metrics.tsv**: found in the directory where you originally placed your reads, and it contains many of the basic organism/assembly/mapping metrics (all of this can be found visually in the multiqc report)
2. **multiqc_report.html**: opens in your browser and visualizes many of the metrics from your analysis.
    *Sections of MultiQC Report*
    - WARNINGS : Indicates if your data passed the thresholds contained within config file
    - ORGANISMAL INFO : contains information pertaining to organism identified, and serotype
    - ASSEMBLY SPECS : Assembly and read mapping metrics
    - MASH : Organisms identified by Min-Hash method
    - MetaPhlAn : Characterizes metagenomic composition of your samples (good for identifying contamination)
    - Abricate Sections : Each Abricate section contains a heatmap of AR genes identified. Heatmaps of putative drug resistance is identified for CARD, NCBI, and Resfinder
    - **NOTE: For abricate predicted phenotype heatmaps, if identified resistance genes are overlapping, they currently will be counted as separate hits for AR (i.e. may drive up counts)**
    - Drug Resistance Sections : Putative drug resistance profiles based upon abricate outputs. **NOTE: Only for edification purposes, NOT diagnostic**
    - QUAST : QC metrics for your assemblies (most of this is contained within the "Assembly Specs" section)
    - FastQC : Basic metrics for read quality. `PE1` or `PE2` extensions reflect trimmed reads.
3. **output_metrics directory**: contains summary files and more detailed metrics for the various BI tools implemented here
    - `abricate_database_summary` files contain summations of Antibiotic Resistance queries of each database by Abricate
    - `quast_metrics.tsv`: represented in `QUAST` section of multiqc report.
    - `metagenomic_composition.tsv`: represented in `MetaPhlAn` sectoin of multiqc report. Contains metagenomic identification of each sample
    - `top_hits_MASH.tsv`: Contains the top-10 organismal hits for each isolate. If your isolate did not pass threshold for ID with MASH, look in this file to identify other potential hits.
    - *Serotypes*- all of these are represented in the `ORGANISMAL INFO` section of the multiqc report.
      - `e_coli_serotypes.tsv`: serotyping summary for Escherichia coli samples
      - `seqsero2_serotypes.tsv`: serotyping summary for Salmonella samples
      - `emm_serotypes.tsv`: serotyping summary for Streptococcus pyogenes samples
      - `legsta_serotypes.tsv`: serotyping summary for Legionella pneumophila samples
      - `seroba_serotypes.tsv`: serotyping summary for Streptococcus pneumoniae samples
    - `QC-metrics_trimmed.tsv`: qualtiy metrics for trimmed reads. Some of this information is represented in `ASSEMBLY SPECS` of the multiqc report
    - `warning_flags.tsv`: Contains PASS/FAIL evaluations of your data based upon the thresholds in the config file.

# Follow-on Tree Gen Pipeline

[Roary](https://sanger-pathogens.github.io/Roary/)

[Lyve-SET](https://github.com/lskatz/lyve-SET)

### Troubleshooting
Please feel free to contact the developer, Andrew Lang, directly at Andrew.Lang@massmail.state.ma.us, or open an issue on this repo. Additionally, we are continuously looking for ways to improve the pipeline. If you have suggestions, please don't be shy!
