nextflow.enable.dsl = 2

// Change is required! Specify your projectDir here
projectDir = "/scratch/rykalinav/rki_tsi/"


// Parameters for shiver
params.alientrimmer = "${projectDir}/bin/AlienTrimmer.jar"
params.gal_primers = "${projectDir}/data/primers_GallEtAl2012.fasta"
params.illumina_adapters = "${projectDir}/data/adapters_Illumina.fasta"



log.info """
====================================================
                  TSI PIPELINE
====================================================
             Author: Vera Rykalina
       Affiliation: Robert Koch Institute 
        Acknowledgement: Tanay Golubchik
              Created: 25 June 2024
           Last Updated: 25 July 2024
====================================================
         """

// error codes
params.profile = null
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}


Set modes = ['paired', 'single']
if ( ! (params.mode in modes) ) {
    exit 1, "Unknown mode. Choose from " + modes
}


process RAW_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
  publishDir "${params.outdir}/01_raw_fastqc/${id}", mode: "copy", overwrite: true
 // debug true

  input:
    tuple val(id), path(reads)
  output:
    path "${id}*_fastqc.html", emit: html
    path "${id}*_fastqc.zip",  emit: zip
  script:
    
    """
    [ -f *R1*.fastq.gz ] && mv *R1*.fastq.gz ${id}_raw.R1.fastq.gz
    [ -f *R2*.fastq.gz ] && mv *R2*.fastq.gz ${id}_raw.R2.fastq.gz
    
    fastqc *.fastq.gz
  
    """
  
}


process FASTP {
  label "fastp"
  conda "${projectDir}/env/fastp.yml"
  publishDir "${params.outdir}/02_fastp_trimmed/${id}", mode: "copy", overwrite: true
  //debug true

  input:
    tuple val(id), path(reads)

  output:
    tuple val(id), path("${id}_fastp.R{1,2}.fastq.gz"), emit: reads
    tuple val(id), path("${id}_fastp.json"),            emit: json
    tuple val(id), path("${id}_fastp.html"),            emit: html

 script:
    set_paired_reads = params.mode == 'single' ? '' : "--in2 ${reads[1]} --out2 ${id}_fastp.R2.fastq.gz --unpaired1 ${id}.SE.R1.fastq.gz --unpaired2 ${id}.SE.R2.fastq.gz"
    """
    
    fastp \
        --in1 ${reads[0]} \
        --out1 ${id}_fastp.R1.fastq.gz \
        ${set_paired_reads} \
        --adapter_fasta ${params.illumina_adapters} \
        --json ${id}_fastp.json \
        --html ${id}_fastp.html \
        --low_complexity_filter \
        --overrepresentation_analysis \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread ${task.cpus}

    """
}

process FASTP_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
  publishDir "${params.outdir}/03_trimmed_fastqc/${id}", mode: "copy", overwrite: true
 // debug true

  input:
    tuple val(id), path(reads)
  output:
    path "${id}*_fastqc.html", emit: html
    path "${id}*_fastqc.zip",  emit: zip
 
  script:
    
    """
    fastqc ${reads}
    """
}


process ALIENTRIMMER {
  conda "${projectDir}/env/multiqc.yml"
  publishDir "${params.outdir}/04_primer_trimmed/${id}", mode: "copy", overwrite: true
  debug true
  
  input:
     tuple val(id), path(reads)
  
  output:
    tuple val(id), path("${id}_alientrimmer.R.{1,2}.fastq.gz"), emit: reads
    tuple val(id), path("${id}_alientrimmer.R.S.fastq.gz"), emit: singletons


  script:

  if (params.mode == "paired"){
  """
  java -jar ${params.alientrimmer} \
       -1 ${reads[0]} \
       -2 ${reads[1]} \
       -a ${params.gal_primers} \
       -o ${id}_alientrimmer.R \
       -k 9 \
       -z
  """
  } else if (params.mode == "single") {
    """
      java -jar ${params.alientrimmer} \
           -i ${reads[0]} \
           -a ${params.gal_primers} \
           -o ${id}_alientrimmer.R \
           -k 15 \
           -z
    """
  }

}


process ALIENTRIMMER_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
  publishDir "${params.outdir}/05_alientrimmed_fastqc/${id}", mode: "copy", overwrite: true
 // debug true

  input:
    tuple val(id), path(reads)
  output:
    path "${id}*_fastqc.html", emit: html
    path "${id}*_fastqc.zip",  emit: zip
 
  script:
    
    """
    fastqc ${reads}
    """
}

process MULTIQC {
  conda "${projectDir}/env/multiqc.yml"
  publishDir "${params.outdir}/06_multiqc", mode: "copy", overwrite: true
  debug true
  
  input:
    path report_files

  output:
    path "multiqc_report.html", emit: report
 
  script:
  """
  multiqc .
  """
}



// **************************************INPUT CHANNELS***************************************************
if ( !params.fastq ) {
    exit 1, "input missing, use [--fastq]"
}

if (params.mode == 'paired') {
        ch_input_fastq = Channel
        .fromFilePairs( params.fastq, checkIfExists: true )
        .map{ tuple ( it[0].split("HIV")[1].split("_")[0], [it[1][0], it[1][1]]) }
        //.fromFilePairs( "${projectDir}/RawData/*_R{1,2}*.fastq.gz", checkIfExists: true )

        
} else { ch_input_fastq = Channel
        .fromPath( params.fastq, checkIfExists: true )
        //.fromPath( "${projectDir}/RawData/*.fastq.gz", checkIfExists: true )
        .map { file -> [file.simpleName, [file]]}
        .map {tuple ( it[0].split("HIV")[1].split("_")[0], it[1][0])}.view()
}



workflow {
    ch_raw_fastqc = RAW_FASTQC ( ch_input_fastq )
    ch_fastp_trimmed = FASTP ( ch_input_fastq )
    ch_fastp_fastqc = FASTP_FASTQC ( ch_fastp_trimmed.reads) 
    ch_primer_trimmed = ALIENTRIMMER ( ch_fastp_trimmed.reads)
    ch_alientrimmer_fastqc = ALIENTRIMMER_FASTQC ( ch_primer_trimmed.reads) 
    ch_multiqc = MULTIQC ( ch_raw_fastqc.zip.concat(ch_fastp_fastqc.zip).concat(ch_alientrimmer_fastqc.zip).collect() )

    

}


// fastaq primer trimming
//fastaq sequence_trim 07-00462_fastp.R1.fastq.gz 07-00462_fastp.R2.fastq.gz 07-00462_fastp_trimmed.R1.fastq.gz \
//07-00462_fastp_trimmed.R2.fastq.gz../../../DataShiverInit/primers_GallEtAl2012.fasta --revcomp