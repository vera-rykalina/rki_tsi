# HiVtime Pipeline
This is a seamless pipeline developed to automate HIV-phyloTSI method. The [method](https://www.medrxiv.org/content/10.1101/2022.05.15.22275117v1) is meant to estimate time since infection (TSI) from HIV deep-sequencing data. Here is a link for the HIV-phyloTSI GitHub [page](https://github.com/BDI-pathogens/HIV-phyloTSI/tree/main).

## Pipeline workflow
![Plot](/images/HIVtime.png)

## Tools
The pipeline includes the following tools:
- KRAKEN2
- FASTP
- ALIENTRIMMER
- FASTQC
- MULTIQC
- SPADES
- CD-HIT
- SHIVER (C. Wymant)
- KALLISTO
- PHYLOSCANNER (C. Wymant & M. Hall)
- IQTREE
- HIV-phyloTSI (T. Golubchik)

## Installation
Make sure that R with the phyloscannerR package installed globally. All other software tools with their dependencies get installed automatically via conda directives. 

Download Kraken2 database, we currently use kraken2_nt_20231129.

Then use the path to kraken DB for the parameter --krakendb.

Make it default by krakendb = "/yourlocation/kraken2_nt_20231129"

The [database](https://benlangmead.github.io/aws-indexes/k2) used by the pipeline is nt Database (29 November, 2023) which is stored within HPC resources. 

Alternativelly, one can install all the tools manually, using .yml recipes (see Environments folder).

## Usage

We recommend to install mamba for quicker installation process (in base env)

```sh
conda install mamba -c conda-forge
```

Create a nextflow environment, using **nextflow.yml** file:

```sh
conda env create -n nextflow -f Environmets/nextflow.yml
```

Activate the *nextflow* environment:
```sh
conda activate nextflow
```

Run the pipeline (example run, for more information use --help): 

```sh
nextflow hivtime.nf \
 -c hivtime_profile.config \
 --outdir Results \
 --fastq "/path_to_your_fastq_files/*R{1,2}*.fastq.gz" \
 -profile rki_slurm,rki_mamba \
--krakendb /path_to_krakendb/kraken2_nt_20231129/ \
--mode paired \
--alignment ../rki_tsi/data/alignments/HIV1_COM_2022_genome_DNA.fasta \
--primers ../rki_tsi/data/primers/primers_GallEtAl2012.fasta \
-resume (optional)
```
