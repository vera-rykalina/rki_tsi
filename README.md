# HIVtime
HIVtime is a seamless pipeline developed to implement the HIV-phyloTSI method. The [method](https://www.medrxiv.org/content/10.1101/2022.05.15.22275117v1) was designed to estimate time since infection (TSI) using HIV deep-sequencing data. The source code and documentation are available on the HIV-phyloTSI GitHub repository [page](https://github.com/BDI-pathogens/HIV-phyloTSI/tree/main).

## Pipeline workflow
![Plot](/images/HIVtime_400_dpi.png)

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
The pipeline is written in Nextflow, which can be easily installed using conda.

- Install conda if not installed:

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

- Install mamba if not installed (recommended):

```sh
conda install mamba -n base -c conda-forge
```

- Clone the repo

```sh
git clone https://github.com/vera-rykalina/rki_tsi.git
```

- Install nextflow, using **nextflow.yml** file

```sh
conda env create -n nextflow -f env/nextflow.yml
```

All other software tools with their dependencies get installed automatically within the pipeline via conda directives. 

## Prerequisite
- Download Kraken2 database

Then use the path to kraken DB for the parameter **--krakendb** or make it default by an assignment: krakendb = "/yourlocation/kraken2_nt_20231129"

I personally use this [database](https://benlangmead.github.io/aws-indexes/k2) (nt Database 29 November, 2023).



## Usage
- Activate the *nextflow* environment:
```sh
conda activate nextflow
```
- Run the pipeline (example run, for more information use **--help**): 

```sh
nextflow hivtime.nf \
 -c hivtime_profile.config \
 --outdir Results \
 --fastq "/path_to_your_fastq_files/*R{1,2}*.fastq.gz" \
 -profile rki_slurm,rki_mamba \
--krakendb /path_to_krakendb/kraken2_nt_20231129/ \
--mode paired \
--primers ../rki_tsi/data/primers/primers_GallEtAl2012.fasta \
-resume (optional)
```
