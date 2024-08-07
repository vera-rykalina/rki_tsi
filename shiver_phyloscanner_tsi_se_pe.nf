nextflow.enable.dsl = 2

// Run

/* nextflow shiver_phyloscanner_tsi_se_pe.nf \
-c rki_profile.config 
--outdir output 
--fastq "/scratch/rykalinav/rki_recency/Pipeline/RawData/*R{1,2}*.fastq.gz" 
-profile rki_slurm,rki_mamba \
--krakendb /scratch/databases/kraken2_nt_20231129/ \
--mode paired 
-resume
*/

// Change is required! Specify your projectDir here
projectDir = "/scratch/rykalinav/rki_tsi/"


// Parameters for shiver
params.alientrimmer = "${projectDir}/bin/AlienTrimmer.jar"
params.gal_primers = "${projectDir}/data/primers_GallEtAl2012.fasta"
params.illumina_adapters = "${projectDir}/data/adapters_Illumina.fasta"


krakendb = params.krakendb
// taxid of HIV-1 
params.taxid = "11676"


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


// kraken2
process CLASSIFY {

    label "kraken"
    conda "${projectDir}/env/kraken.yml"
    publishDir "${params.outdir}/02_classified_reads/${id}", mode: "copy", overwrite: true, pattern: "*.txt"

    input:
        tuple val(id), path(reads)
        val (krakendb)

    output:
        tuple val(id), path("${id}_classified.R*.fastq"),     emit: classified_fastq
        tuple val(id), path("${id}_unclassified.R*.fastq"),   emit: unclassified_fastq
        tuple val(id), path("${id}_kraken.out.txt"),          emit: kraken_output
        tuple val(id), path("${id}_kraken.report.txt"),       emit: kraken_report


    script:
        set_paired = params.mode == 'paired' ? '--paired' : ''
        set_out_name = params.mode == 'paired' ? '#' : ''

            """
             kraken2 \
              --threads ${task.cpus} \
              --db ${krakendb} \
              ${set_paired} \
              --classified-out ${id}_classified.R${set_out_name}.fastq \
              --unclassified-out ${id}_unclassified.R${set_out_name}.fastq \
              --output ${id}_kraken.out.txt \
              --report ${id}_kraken.report.txt \
              --gzip-compressed \
              ${reads}
          """     
}


// krakentools
process EXTRACT {
    label "krakentools"
    conda "${projectDir}/env/krakentools.yml"
    publishDir "${params.outdir}/03_filtered_reads/${id}", mode: "copy", overwrite: true


    input:
        tuple val(id), path(reads)
        tuple val(id), path(kraken_output)
        tuple val(id), path(kraken_report)
        val (taxid)
 
    
    output:
        tuple val(id), path("${id}_filtered.R*.fastq")
    
    script:
        set_paired_reads = params.mode == 'single' ? '' : "-2 ${reads[1]} -o2 ${id}_filtered.R2.fastq"
           """
            extract_kraken_reads.py \
                -1 ${reads[0]} \
                -o ${id}_filtered.R1.fastq \
                ${set_paired_reads} \
                -k ${kraken_output} \
                --report ${kraken_report} \
                --include-children \
                --taxid ${taxid} \
                --fastq-output
        """
}


// merge unclassified and filtered reads
process MERGE {
    publishDir "${params.outdir}/04_merged_reads", failOnError: true, mode: "copy", overwrite: true

    input:
        tuple val(id), path(unclassified), path(filtered)

    output:
        tuple val("${id}"), path("${id}_filtered.R*.fastq.gz")

    script:
    if (params.mode == "paired") {
        """
        gzip -c ${unclassified[0]} ${filtered[0]} > ${id}_filtered.R1.fastq.gz
        gzip -c ${unclassified[1]} ${filtered[1]} > ${id}_filtered.R2.fastq.gz
        """
    } else if (params.mode == "single") {
       """
       gzip -c ${unclassified[0]} ${filtered[0]} > ${id}_filtered.R1.fastq.gz
       """
    }
}


process KRAKEN_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
  publishDir "${params.outdir}/05_kraken_fastqc/${id}", mode: "copy", overwrite: true
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

process FASTP {
  label "fastp"
  conda "${projectDir}/env/fastp.yml"
  publishDir "${params.outdir}/06_fastp_trimmed/${id}", mode: "copy", overwrite: true
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
        --length_required 80 \
        --thread ${task.cpus}

    """
}

process FASTP_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
  publishDir "${params.outdir}/07_trimmed_fastqc/${id}", mode: "copy", overwrite: true
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
  publishDir "${params.outdir}/08_primer_trimmed/${id}", mode: "copy", overwrite: true
  //debug true
  
  input:
     tuple val(id), path(reads)
  
  output:
    tuple val(id), path("${id}_alientrimmer.R{1,2}.fastq.gz"), emit: reads
    //tuple val(id), path("${id}_alientrimmer.R.S.fastq.gz"), emit: singletons


  script:

  if (params.mode == "paired"){
  """
  java -jar ${params.alientrimmer} \
       -1 ${reads[0]} \
       -2 ${reads[1]} \
       -a ${params.gal_primers} \
       -o ${id}_alientrimmer.R \
       -k 15 \
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
    mv ${id}_alientrimmer.R.fastq.gz ${id}_alientrimmer.R1.fastq.gz
    """
  }

}


process ALIENTRIMMER_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
  publishDir "${params.outdir}/09_alientrimmed_fastqc/${id}", mode: "copy", overwrite: true
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
  publishDir "${params.outdir}/10_multiqc", mode: "copy", overwrite: true
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


process SPADES {
  label "spades"
  conda "${projectDir}/env/spades.yml"
  publishDir "${params.outdir}/11_spades", mode: "copy", overwrite: true
  
  input:
    tuple val(id), path(reads) 

  output:
    path "${id}"
    tuple val (id), path ("${id}/${id}_spades_contigs.fasta"), emit: spadescontigs
 
  script:
    if (params.mode == "paired"){
    """
    spades.py \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${id} \
    --only-assembler \
    --threads ${task.cpus}
  

    mv ${id}/contigs.fasta ${id}/${id}_spades_contigs.fasta
    """

 }   else if (params.mode == "single") {
    """
    spades.py \
    -s ${reads[0]} \
    -o ${id} \
    --only-assembler \
    --threads ${task.cpus}

    mv ${id}/contigs.fasta ${id}/${id}_spades_contigs.fasta
    """

  }

}

process METASPADES {
  label "spades"
  conda "${projectDir}/env/spades.yml"
  publishDir "${params.outdir}/12_metaspades", mode: "copy", overwrite: true
  
  input:
    tuple val(id), path(reads) 

  output:
    path "${id}"
    tuple val(id), path("${id}/${id}_metaspades_contigs.fasta"), emit: spadescontigs
 
  script:
    if (params.mode == "paired") {
    """
    spades.py \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${id} \
    --only-assembler \
    --meta \
    --threads ${task.cpus}
  

    mv ${id}/contigs.fasta ${id}/${id}_metaspades_contigs.fasta
    """
    } else if (params.mode == "single") {
    
    """
    mkdir ${id}
    touch ${id}/${id}_metaspades_contigs.fasta
    """
    }

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
        .map {tuple ( it[0].split("HIV")[1].split("_")[0], it[1][0])}

}



workflow {
    ch_raw_fastqc = RAW_FASTQC ( ch_input_fastq )
    ch_classified_reads = CLASSIFY ( ch_input_fastq, krakendb )
    ch_filtered_reads = EXTRACT ( ch_classified_reads.classified_fastq, ch_classified_reads.kraken_output, ch_classified_reads.kraken_report, params.taxid )
    ch_merged_reads = MERGE ( ch_classified_reads.unclassified_fastq.combine(ch_filtered_reads, by:0) )
    ch_kraken_fastqc = KRAKEN_FASTQC ( ch_merged_reads )
    ch_fastp_trimmed = FASTP (  ch_merged_reads )
    ch_fastp_fastqc = FASTP_FASTQC ( ch_fastp_trimmed.reads) 
    ch_primer_trimmed = ALIENTRIMMER ( ch_fastp_trimmed.reads)
    ch_alientrimmer_fastqc = ALIENTRIMMER_FASTQC ( ch_primer_trimmed.reads) 
    //ch_multiqc = MULTIQC ( ch_raw_fastqc.zip.concat(ch_fastp_fastqc.zip).concat(ch_alientrimmer_fastqc.zip).concat(ch_kraken_fastqc.zip).collect() )
    ch_multiqc = MULTIQC ( ch_kraken_fastqc.zip.concat(ch_alientrimmer_fastqc.zip).concat(ch_fastp_fastqc.zip).concat(ch_raw_fastqc.zip).collect() )
    ch_spades = SPADES ( ch_primer_trimmed.reads )
    ch_metaspades = METASPADES ( ch_primer_trimmed.reads )
    ch_spades_combined = ch_spades.spadescontigs.combine(ch_metaspades.spadescontigs, by:0).view()
    ch_all_contigs = ch_spades_combined.collectFile(name: "contigs.fasta", storeDir: "${projectDir}/${params.outdir}/13_all_contigs")
    

}


// fastaq primer trimming
//fastaq sequence_trim 07-00462_fastp.R1.fastq.gz 07-00462_fastp.R2.fastq.gz 07-00462_fastp_trimmed.R1.fastq.gz \
//07-00462_fastp_trimmed.R2.fastq.gz../../../DataShiverInit/primers_GallEtAl2012.fasta --revcomp

// cd-hit-est -i reads.fa -o output.fa -c 0.95 -n 10 -d 999 -M 0 -T 0