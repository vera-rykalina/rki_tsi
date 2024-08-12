nextflow.enable.dsl = 2

// Run

/* nextflow shiver_phyloscanner_tsi_se_pe.nf \
-c rki_profile.config 
--outdir output 
--fastq "/scratch/rykalinav/rki_recency/Pipeline/RawData/*R{1,2}*.fastq.gz" 
-profile rki_slurm,rki_mamba \
--krakendb /scratch/databases/kraken2_nt_20231129/ \
--mode paired 
--alignment /scratch/rykalinav/rki_tsi/data/HIV1_COM_2021_genome_DNA.fasta \
--primers /scratch/rykalinav/rki_tsi/data/primers_GallEtAl2012.fasta \
-resume
*/

// Change is required! Specify your projectDir here
projectDir = "/scratch/rykalinav/rki_tsi/"


// Parameters for shiver
params.alientrimmer = "${projectDir}/bin/AlienTrimmer.jar"
params.illumina_adapters = "${projectDir}/data/adapters_Illumina.fasta"
params.config = "${projectDir}/bin/config.sh"
params.remove_whitespace = "${projectDir}/bin/tools/RemoveTrailingWhitespace.py"
params.alignment = "${projectDir}/data/HIV1_COM_2021_genome_DNA.fasta"

//alignment = params.alignment
primers = params.primers
//params.gal_primers = "${projectDir}/data/primers_GallEtAl2012.fasta"

// Parameters for kraken
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
           Last Updated: 8 August 2024
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
    path "${id}*_fastqc.html", emit: Html
    path "${id}*_fastqc.zip",  emit: Zip
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
        tuple val(id), path("${id}_classified.R*.fastq"),     emit: ClassifiedFastq
        tuple val(id), path("${id}_unclassified.R*.fastq"),   emit: UnclassifiedFastq
        tuple val(id), path("${id}_kraken.out.txt"),          emit: KrakenOutput
        tuple val(id), path("${id}_kraken.report.txt"),       emit: KrakenReport


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
    publishDir "${params.outdir}/04_merged_reads/${id}", failOnError: true, mode: "copy", overwrite: true

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
    path "${id}*_fastqc.html", emit: Html
    path "${id}*_fastqc.zip",  emit: Zip
 
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
    tuple val(id), path("${id}_fastp.R{1,2}.fastq.gz"), emit: Reads
    tuple val(id), path("${id}_fastp.json"),            emit: Json
    tuple val(id), path("${id}_fastp.html"),            emit: Html

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
    path "${id}*_fastqc.html", emit: Html
    path "${id}*_fastqc.zip",  emit: Zip
 
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
     val (primers)
  
  output:
    tuple val(id), path("${id}_alientrimmer.R*.fastq.gz")



  script:

  if (params.mode == "paired"){
  """
  java -jar ${params.alientrimmer} \
       -1 ${reads[0]} \
       -2 ${reads[1]} \
       -a ${primers} \
       -o ${id}_alientrimmer.R \
       -k 15 \
       -z
  """
  } else if (params.mode == "single") {
    """
      java -jar ${params.alientrimmer} \
           -i ${reads[0]} \
           -a ${primers} \
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
    path "${id}*_fastqc.html", emit: Html
    path "${id}*_fastqc.zip",  emit: Zip
 
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
    path "multiqc_report.html"
 
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
    tuple val (id), path ("${id}/${id}_spades_contigs.fasta"), emit: SpadesContigs
 
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
    tuple val(id), path("${id}/${id}_metaspades_contigs.fasta"), emit: MetaspadesContigs
 
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


process MERGE_CONTIGS {
  publishDir "${params.outdir}/14_merged_contigs", mode: "copy", overwrite: true
  
  input:
    tuple val(id), path(spades_contigs), path(metaspades_contigs) 

  output:
    tuple val (id), path ("${id}_merged_contigs.fasta")
 
  script:

    """
   cat ${spades_contigs} ${metaspades_contigs} > "${id}_merged_contigs.fasta"
    """
}


// cluster contigs
process CD_HIT_EST {
  label "cd_hit_est"
  conda "${projectDir}/env/cd-hit.yml"
  publishDir "${params.outdir}/15_clustered_contigs", mode: "copy", overwrite: true
  
  input:
    tuple val(id), path(contigs) 

  output:
    tuple val (id), path ("${id}_clustered_contigs.fasta")
 
  script:

    """
    cd-hit-est \
       -i ${contigs} \
       -o ${id}_clustered_contigs.fasta \
       -c 0.9 \
       -n 10 \
       -d 999 \
       -l 299 \
       -T ${task.cpus}
    """

}


// SHIVER PART (including KALLISTO)
process SHIVER_INIT {
  conda "${projectDir}/env/shiver.yml"
  publishDir "${projectDir}/${params.outdir}/16_init_dir", mode: "copy", overwrite: true

  input:
     val (alignment)
     val (primers)
     
  output:
     path "InitDir", emit: InitDir
     path "InitDir/ExistingRefsUngapped.fasta", emit: ExistingRefsUngapped
     path "InitDir/IndividualRefs/*.fasta", emit: IndividualRefs
  script:
  
  """
  shiver_init.sh \
    InitDir \
    ${params.config} \
    ${alignment} \
    ${params.illumina_adapters} \
    ${primers}
  """  
}


process FASTQ_RENAME_HEADER {
  //conda "${projectDir}/env/shiver.yml"
  publishDir "${params.outdir}/17_renamed_reads", mode: "copy", overwrite: true
  debug true

  input:
    tuple val(id), path(reads)
  output:
    tuple val("${id}"), path("${id}_renamed_R{1,2}.fastq.gz")
  
  script:
   if (params.mode == "paired") {
   """
   zcat ${reads[0]} |\
      awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' |\
      sed 's/:N:.*//' |\
      gzip -c > ${id}_renamed_R1.fastq.gz
   
   rm ${reads[0]}

   
    zcat ${reads[1]} |\
      awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' |\
      sed 's/:N:.*//' |\
      gzip -c > ${id}_renamed_R2.fastq.gz
  
   rm ${reads[1]}
   """
} else if (params.mode == "single") {
      """
      zcat ${reads[0]} |\
      awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' |\
      sed 's/:N:.*//' |\
      gzip -c > ${id}_renamed_R1.fastq.gz

    rm ${reads[0]}
    """
  }

}

process KALLISTO_INDEX {
  conda "${projectDir}/env/kallisto.yml"
  publishDir "${projectDir}/${params.outdir}/18_kallisto_idx", mode: "copy", overwrite: true

  input:
     path fasta
  
  output:
     path "*.idx"
 
  script:
  """
  kallisto index --index ExistingRefsUngapped.idx ${fasta}
  """
}

process KALLISTO_QUANT {
  conda "${projectDir}/env/kallisto.yml"
  publishDir "${projectDir}/${params.outdir}/19_kallisto_quant", mode: "copy", overwrite: true
  debug true

  input:
     tuple path(index), val(id), path(reads)

  output:
     tuple val(id), path("${id}/${id}_abundance.tsv")
   
  script:
   if (params.mode == "paired") {
   """
   kallisto quant \
    -i ${index} \
    -o ${id} \
    --plaintext ${reads[0]} ${reads[1]}

    mv ${id}/abundance.tsv ${id}/${id}_abundance.tsv
  """
   } else if (params.mode == "single") {
     """
     kallisto quant \
    -i ${index} \
    -o ${id} \
    --plaintext \
    --single \
    --fragment-length 200 \
    --sd 20 \
    ${reads[0]}

    mv ${id}/abundance.tsv ${id}/${id}_abundance.tsv
    """
   }
}

process BEST_ALIGNMENT {
  publishDir "${projectDir}/${params.outdir}/20_best_ref", mode: "copy", overwrite: true
  debug true

  input:
     path (alignments)
     tuple val(id), path(abundancetsv)

  output:
     tuple val(id), path("${id}_bestRef.fasta")
   
  script:
   """
   BestRef=\$(sort -k5 -g ${abundancetsv} | tail -n1 | cut -f1)
   echo "Sample ID: " ${id} "\nBest Kallisto Reference: " \${BestRef} 
   mv \${BestRef}.fasta  ${id}_bestRef.fasta
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
        .map {tuple ( it[0].split("HIV")[1].split("_")[0], it[1][0])}

}



workflow {
    ch_raw_fastqc = RAW_FASTQC ( ch_input_fastq )
    ch_classified_reads = CLASSIFY ( ch_input_fastq, krakendb )
    ch_filtered_reads = EXTRACT ( ch_classified_reads.ClassifiedFastq, ch_classified_reads.KrakenOutput, ch_classified_reads.KrakenReport, params.taxid )
    ch_merged_reads = MERGE ( ch_classified_reads.UnclassifiedFastq.combine( ch_filtered_reads, by:0 ) )
    ch_kraken_fastqc = KRAKEN_FASTQC ( ch_merged_reads )
    ch_fastp_trimmed = FASTP (  ch_merged_reads )
    ch_fastp_fastqc = FASTP_FASTQC ( ch_fastp_trimmed.Reads ) 
    ch_primer_trimmed = ALIENTRIMMER ( ch_fastp_trimmed.Reads, primers )
    ch_alientrimmer_fastqc = ALIENTRIMMER_FASTQC ( ch_primer_trimmed ) 
    ch_multiqc = MULTIQC ( ch_raw_fastqc.Zip.concat(ch_fastp_fastqc.Zip).concat(ch_alientrimmer_fastqc.Zip).concat(ch_kraken_fastqc.Zip).collect() )
    ch_spades = SPADES ( ch_primer_trimmed )
    ch_metaspades = METASPADES ( ch_primer_trimmed )
    ch_spades_combined = ch_spades.SpadesContigs.combine( ch_metaspades.MetaspadesContigs, by:0 )
    ch_merged_contigs = MERGE_CONTIGS ( ch_spades_combined )
    ch_cd_hit_est = CD_HIT_EST ( ch_merged_contigs )
    ch_fastq_renamed_header = FASTQ_RENAME_HEADER ( ch_merged_reads )
    ch_initdir = SHIVER_INIT ( params.alignment, primers )
    ch_kallisto_index = KALLISTO_INDEX ( ch_initdir.ExistingRefsUngapped )
    ch_kallisto_index_reads = ch_kallisto_index.combine( ch_fastq_renamed_header )
    ch_kallisto_quant = KALLISTO_QUANT( ch_kallisto_index_reads )
    ch_bestRef = BEST_ALIGNMENT ( ch_initdir.IndividualRefs, ch_kallisto_quant )
   

}


// fastaq primer trimming
//fastaq sequence_trim 07-00462_fastp.R1.fastq.gz 07-00462_fastp.R2.fastq.gz 07-00462_fastp_trimmed.R1.fastq.gz \
//07-00462_fastp_trimmed.R2.fastq.gz../../../DataShiverInit/primers_GallEtAl2012.fasta --revcomp

// cd-hit-est -i reads.fa -o output.fa -c 0.9 -n 10 -d 999