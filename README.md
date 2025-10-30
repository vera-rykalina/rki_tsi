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



## Usage Examples (for more information use **--help**)
- Activate the *nextflow* environment:
```sh
conda activate nextflow
```

- PE reads, full-length genomes: 

```sh
nextflow hivtime.nf \
 -c hivtime_profile.config \
 --outdir Results \
 --fastq "/path_to_your_fastq_files/*R{1,2}*.fastq.gz" \
 -profile rki_slurm,rki_conda \
--krakendb /path_to_krakendb/kraken2_nt_20231129/ \
--primers ../rki_tsi/data/primers/primers_GallEtAl2012.fasta \
--genome full \
-resume (optional)
```

- SE reads, full-length genomes: 

```sh
nextflow hivtime.nf \
 -c hivtime_profile.config \
 --outdir Results \
 --fastq "/path_to_your_fastq_files/*R{1,2}*.fastq.gz" \
 -profile rki_slurm,rki_conda \
--krakendb /path_to_krakendb/kraken2_nt_20231129/ \
--primers ../rki_tsi/data/primers/primers_GallEtAl2012.fasta \
--genome full \
--mode single \
-resume (optional)
```

- PE reads, partial genomes: 

```sh
nextflow hivtime.nf \
 -c hivtime_profile.config \
 --outdir Results \
 --fastq "/path_to_your_fastq_files/*R{1,2}*.fastq.gz" \
 -profile rki_slurm,rki_conda \
--krakendb /path_to_krakendb/kraken2_nt_20231129/ \
--primers ../rki_tsi/data/primers/primers_GallEtAl2012.fasta \
--genome partial \
--modelname SK \
-resume (optional)
```

## Models
The current version of the HIVtime includes 2 models: the original HIV-phyloTSI (full genomes) and a new model for partial genomes (SK PCR settings). 

![Plot](/images/hiv_genes_coverage_clipped.png)

Below is the PCR settings scheme which visually represents PCR product coverage. The primer sequences can be found here: _/rki_tsi/data/primers/primers_sk_validation_. The model I built for partial genome includes partial _gag_ and _pol_ regions that cover 26.30% (vs 69.99% of HIV-phyloTSI) of the HXB2 genome (K03455). 

I also developed a fully automated Nextflow pipeline for remodeling to be able to generate models for different lab settings (partial or full _gag_ and _pol_ genomic regions). Below is the pipeline metro map.

## Remodeling pipeline workflow
![Plot](/images/retraining_500_dpi.png)

## Notes

- Input FASTQ files

Our FASTQ files follow a specific filename structure: ^([\w]+_)*HIV\d{2}-\d{5}_([\w]+_)*R[12]_\d{3}\.fastq\.gz$

THis structure contains the HIV prefix before the sample ID (HIV\d{2}-\d{5}). I use this HIV prefix to extract the sample ID with an additional line of code. If you have a different filename pattern, simply comment out this line of code in the hivtime.nf script, like this (or use your own method for splitting the filename):

```sh
if (params.mode == 'paired') {
        ch_input_fastq = Channel
        .fromFilePairs( params.fastq, checkIfExists: true )
        //.map { tuple ( it[0].split("HIV")[1].split("_")[0], [it[1][0], it[1][1]]) } <- change is here

        
} else { ch_input_fastq = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .map { file -> [file.simpleName, [file]]}
       //.map {tuple ( it[0].split("HIV")[1].split("_")[0], it[1][0])} <- change is here

}
```

- Configuration

Please be aware that this pipeline is computationally intensive (especially the phylogeny part), which is why it is configured to run on an HPC system. My Slurm profile (rki_slurm) is optimized for a batch size of about 20 samples. If you plan to use different settings, consider creating your own profile.