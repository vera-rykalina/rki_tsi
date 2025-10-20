nextflow.enable.dsl = 2

// Run

/* nextflow hivtime.nf \
-c hivtime_profile.config \
--outdir output_folder_name \
--fastq "/path/rki_tsi/fastq/*R{1,2}*.fastq.gz" \ 
--krakendb /path/databases/kraken2/kraken2_nt_20231129/ \
--primers /path/rki_tsi/data/primers/primers_GallEtAl2012.fasta \
--genome full \
--modelname SK \
-profile rki_slurm,rki_conda \
-resume
*/


//**************************************************PARAMETERS*******************************************************
// Parameters for kraken
// a newer db: kraken2_nt_20240904
//krakendb = params.krakendb

// taxid of HIV-1 
params.taxid = "11676"


// Parameters for shiver
params.alientrimmer = "${projectDir}/bin/AlienTrimmer.jar"
params.adapters = "${projectDir}/data/adapters/adapters.fasta"
params.config_se = "${projectDir}/bin/config_se.sh"
params.config_pe = "${projectDir}/bin/config_pe.sh"
params.alignment = "${projectDir}/data/alignments/HIV1_COM_2022_genome_DNA.fasta"
//primers = params.primers



// Parameters for phyloscanner
params.two_refs = "${projectDir}/data/refs/2refs_HXB2_C.BW.fasta"
params.excision_coordinates = "${projectDir}/data/phyloscanner/DrugResistancePositionsInHXB2.txt"
params.windows_oneline = "${projectDir}/data/phyloscanner/windows250_VR_norms_oneline.txt"
params.hiv_distance_normalisation = "${projectDir}/data/phyloscanner/HIV_DistanceNormalisationOverGenome.csv"
params.k = 15

// Parameters for HIV-PhyloTSI
params.seqprotocol = "amplicons"
params.mode = "paired"


params.normalisation = "${projectDir}/bin/tools/CalculateTreeSizeInGenomeWindows.py"



// help message
params.help = false

if (params.help) { exit 0, helpMSG() }


// error codes
params.profile = null
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory, use [--outdir]"
}

Set modes = ['paired', 'single']
if ( ! (params.mode in modes) ) {
    exit 1, "Unknown mode. Choose from " + modes
}


Set seqprotocols = ['amplicons', 'capture']
if ( ! (params.seqprotocol in seqprotocols) ) {
    exit 1, "Unknown protocol. Choose from " + seqprotocols
}


if ( !params.fastq ) {
    exit 1, "Missing input, use [--fastq]"
}

params.krakendb = null
if ( !params.krakendb ) {
    exit 1, "Missing input, use [--krakendb]"
}

params.primers = null
if ( !params.primers ) {
    exit 1, "Missing input, use [--primers]"
}


Set genomes = ['full', 'partial']
if ( ! (params.genome in genomes) ) {
    exit 1, "Unknown genome. Choose from " + genomes
}

if ( params.genome == 'partial' ) {
  params.modelfolder = "${projectDir}/bin/models/partial"
  modlelname = params.modelname
} else if ( params.genome == 'full' ) {
  params.modelfolder = "${projectDir}/bin/models/full"
}


params.modelname = null
if (params.genome == "partial" && !params.modelname ) {
    exit 1, "Missing model name, use [--modelname]"
}

def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_red = "\u001B[31m";
    c_dim = "\033[2m";
    log.info """
  

    ${c_blue}HIVtime: A Seamless Pipeline for Estimating Time since HIV Infection${c_blue}
    ====================================================================================================================================================
    Author: Vera Rykalina
    ${c_blue}Affiliation: Robert Koch Institute${c_blue}
    Acknowledgement: Tanya Golubchik, Chris Wymant
    Created: 25 June 2024
    ====================================================================================================================================================
  

    ${c_yellow}Usage examples:${c_reset}

    nextflow hivtime.nf -c profile.config --fastq '*R{1,2}.fastq.gz' --krakendb db --primers primers.fasta -profile profile --outdir output --genome full

    nextflow hivtime.nf -c profile.config --fastq '*R{1,2}.fastq.gz' --krakendb db --primers primers.fasta -profile profile --outdir output --genome partial --modelname SK
   
    
    ${c_green}Settings:${c_reset}  
    
    --fastq             Path to a FASTQ files e.g.: '*R{1,2}*.fastq.gz'.

    --krakendb          Path to a Kraken2 database. [recommended: kraken2_nt_20231129].
    
    --primers           Path of a FASTA file containing the primer sequences to be clipped.
    
    --outdir            Name for an output directory e.g. output [string].

    --mode              Choose from [paired, single] [default: paired].

    --adapters          Define the path of a FASTA file containing the adapter sequences to be clipped [default: data/adapters/adapters.fasta].

    --alignment         Define the path of a FASTA file containing the HIV alignment [default: data/alignments/HIV1_COM_2022_genome_DNA.fasta].
    
    --seqprotocol       Choose from [amplicons, capture]  [default: amplicons].

    --genome            Choose from [full, partial].

    --modelname         Only valid for partial genomes. Use model name (e.g. SK, IGS etc.)

    """
}



process RAW_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
 //publishDir "${params.outdir}/raw_fastqc/${id}", mode: "copy", overwrite: true
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


process FASTP {
  label "fastp"
  conda "${projectDir}/env/fastp.yml"
  //publishDir "${params.outdir}/fastp_trimmed/${id}", mode: "copy", overwrite: true
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
        --adapter_fasta ${params.adapters} \
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
 //publishDir "${params.outdir}/trimmed_fastqc/${id}", mode: "copy", overwrite: true
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
  //publishDir "${params.outdir}/primer_trimmed/${id}", mode: "copy", overwrite: true
  //debug true
  
  input:
     tuple val(id), path(reads)
     val (params.primers)
  
  output:
    tuple val(id), path("${id}_alientrimmer.R*.fastq.gz")

  
 // Add options: -l 80 \ -q 20 \ 
  script:

  if (params.mode == "paired"){
  """
  java -jar ${params.alientrimmer} \
       -1 ${reads[0]} \
       -2 ${reads[1]} \
       -a ${params.primers} \
       -o ${id}_alientrimmer.R \
       -k 15 \
       -l 80 \
       -q 20 \
       -z

  rm -f ${id}_alientrimmer.R.S.fastq.gz

  mv ${id}_alientrimmer.R.1.fastq.gz ${id}_alientrimmer.R1.fastq.gz
  mv ${id}_alientrimmer.R.2.fastq.gz ${id}_alientrimmer.R2.fastq.gz
  """
  
  } else if (params.mode == "single") {
    """
      java -jar ${params.alientrimmer} \
           -i ${reads[0]} \
           -a ${params.primers} \
           -o ${id}_alientrimmer.R \
           -k 15 \
           -l 80 \
           -q 20 \
           -z
    mv ${id}_alientrimmer.R.fastq.gz ${id}_alientrimmer.R1.fastq.gz
    """
  }

}


process ALIENTRIMMER_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
 //publishDir "${params.outdir}/alientrimmed_fastqc/${id}", mode: "copy", overwrite: true
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



// kraken2
process CLASSIFY {

    label "kraken"
    conda "${projectDir}/env/kraken.yml"
    publishDir "${params.outdir}/01_classified_reads/${id}", mode: "copy", overwrite: true, pattern: "*.txt"

    input:
        tuple val(id), path(reads)
        val (params.krakendb)

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
              --db ${params.krakendb} \
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
    //publishDir "${params.outdir}/filtered_reads/${id}", mode: "copy", overwrite: true


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
    //publishDir "${params.outdir}/merged_reads/${id}", failOnError: true, mode: "copy", overwrite: true

    input:
        tuple val(id), path(unclassified), path(filtered)

    output:
        tuple val("${id}"), path("${id}_kraken.R*.fastq.gz")

    script:
    if (params.mode == "paired") {
        """
        cat ${unclassified[0]} ${filtered[0]} | gzip -c > ${id}_kraken.R1.fastq.gz
        cat ${unclassified[1]} ${filtered[1]} | gzip -c > ${id}_kraken.R2.fastq.gz
        """
    } else if (params.mode == "single") {
       """
       cat ${unclassified[0]} ${filtered[0]} | gzip -c > ${id}_kraken.R1.fastq.gz
       """
    }
}


process KRAKEN_FASTQC {
  conda "${projectDir}/env/fastqc.yml"
  //publishDir "${params.outdir}/kraken_fastqc/${id}", mode: "copy", overwrite: true
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
  publishDir "${params.outdir}/02_multiqc", mode: "copy", overwrite: true
  debug true
  
  input:
    path report_files

  output:
    path "multiqc_report.html",                        emit: Report
    path "multiqc_data/multiqc_data.json",             emit: Json
    path "multiqc_data/multiqc_general_stats.txt",     emit: Txt
 
  script:
  """
  multiqc .
  """
}


process SPADES {
  label "spades"
  conda "${projectDir}/env/spades.yml"
  //publishDir "${params.outdir}/spades", mode: "copy", overwrite: true
  
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
    --threads ${task.cpus} || mkdir -p ${id} && touch ${id}/contigs.fasta
  

    mv ${id}/contigs.fasta ${id}/${id}_spades_contigs.fasta
    """

 }   else if (params.mode == "single") {
    """
    spades.py \
    -s ${reads[0]} \
    -o ${id} \
    --only-assembler \
    --threads ${task.cpus} || mkdir -p ${id} && touch ${id}/contigs.fasta

    mv ${id}/contigs.fasta ${id}/${id}_spades_contigs.fasta
    """

  }

}

process METASPADES {
  label "spades"
  conda "${projectDir}/env/spades.yml"
  //publishDir "${params.outdir}/metaspades", mode: "copy", overwrite: true
  
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
    --threads ${task.cpus} || mkdir -p ${id} && touch ${id}/contigs.fasta
  

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
  //publishDir "${params.outdir}/merged_contigs", mode: "copy", overwrite: true
  
  input:
    tuple val(id), path(spades_contigs), path(metaspades_contigs) 

  output:
    tuple val (id), path ("${id}_merged_contigs.fasta")
 
  script:

    """
   cat ${spades_contigs} ${metaspades_contigs} > ${id}_merged_contigs.fasta
    """
}


// cluster contigs
process CD_HIT_EST {
  label "cd_hit_est"
  conda "${projectDir}/env/cd-hit.yml"
  publishDir "${params.outdir}/03_clustered_contigs", mode: "copy", overwrite: true
  
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
       -T ${task.cpus} || touch ${id}_clustered_contigs.fasta
    """

}


// SHIVER PART (including KALLISTO)
process SHIVER_INIT {
  conda "${projectDir}/env/shiver.yml"
  //publishDir "${projectDir}/${params.outdir}/init_dir", mode: "copy", overwrite: true

  input:
     val (alignment)
     val (params.primers)
     
  output:
     path "InitDir", emit: InitDir
     path "InitDir/ExistingRefsUngapped.fasta", emit: ExistingRefsUngapped
     path "InitDir/IndividualRefs/*.fasta", emit: IndividualRefs
  script:
  
  if (params.mode == "paired") {
  """
  shiver_init.sh \
    InitDir \
    ${params.config_pe} \
    ${alignment} \
    ${params.adapters} \
    ${params.primers}
  """  
  } else if (params.mode == "single") {

  """
  shiver_init.sh \
    InitDir \
    ${params.config_se} \
    ${alignment} \
    ${params.adapters} \
    ${params.primers}
  """ 
  }
}


process FASTQ_RENAME_HEADER {
  //publishDir "${params.outdir}/renamed_reads", mode: "copy", overwrite: true
  //debug true

  input:
    tuple val(id), path(reads)
  output:
    tuple val("${id}"), path("${id}_renamed_R{1,2}.fastq.gz")
  
  script:
   if (params.mode == "paired") {
   """
   zcat ${reads[0]} |\
     sed 's/kraken:taxid.*//' |\
     sed 's/:N:0:.*//' |\
     awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' |\
     awk '{if (NR%4 == 1)  gsub("/1/","/1",\$1) }1' |\
     gzip -c > ${id}_renamed_R1.fastq.gz
   
   rm ${reads[0]}

   
    zcat ${reads[1]} |\
      sed 's/kraken:taxid.*//' |\
      sed 's/:N:0:.*//' |\
      awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' |\
      awk '{if (NR%4 == 1)  gsub("/2/","/2",\$1) }1' |\
      gzip -c > ${id}_renamed_R2.fastq.gz
  
   rm ${reads[1]}
   """
} else if (params.mode == "single") {
      """
    zcat ${reads[0]} |\
      sed 's/kraken:taxid.*//' |\
      sed 's/:N:0:.*//' |\
      awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' |\
      awk '{if (NR%4 == 1)  gsub("/1/","/1",\$1) }1' |\
      gzip -c > ${id}_renamed_R1.fastq.gz

    rm ${reads[0]}
    """
  }

}

process KALLISTO_INDEX {
  conda "${projectDir}/env/kallisto.yml"
  //publishDir "${projectDir}/${params.outdir}/kallisto_idx", mode: "copy", overwrite: true

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
  //publishDir "${projectDir}/${params.outdir}/kallisto_quant", mode: "copy", overwrite: true
  //debug true

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
  publishDir "${params.outdir}/04_kallisto_best_ref", mode: "copy", overwrite: true
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


process SHIVER_ALIGN {
  label "shiver"
  conda "${projectDir}/env/shiver.yml"
  //publishDir "${params.outdir}/shiver_alignments/${id}", mode: "copy", overwrite: true
  //debug true
  
  input:
    path initdir
    tuple val(id), path(contigs)

  output:
    tuple val("${id}"), path("${id}_wRefs.fasta"), path("${id}.blast")

  script:
    if ( contigs.size() > 0  &&  params.mode == "paired") {
    """
    shiver_align_contigs.sh \
      ${initdir} \
      ${params.config_pe} \
      ${contigs} \
      ${id}

    rm temp_*
    rm *_MergedHits.blast*
    mv ${id}_cut_wRefs.fasta ${id}_wRefs.fasta || mv ${id}_raw_wRefs.fasta ${id}_wRefs.fasta 
    """
  } else if ( contigs.size() > 0  &&  params.mode == "single") {
     """
    shiver_align_contigs.sh \
      ${initdir} \
      ${params.config_se} \
      ${contigs} \
      ${id}

    rm temp_*
    rm *_MergedHits.blast*
    mv ${id}_cut_wRefs.fasta ${id}_wRefs.fasta || mv ${id}_raw_wRefs.fasta ${id}_wRefs.fasta 
    """

   } else {
    """
     printf "There is no contig for sample with ID: ${id}"
     touch ${id}.blast
     touch ${id}_wRefs.fasta
    """
   }
}


process SHIVER_MAP {
  label "shiver"
  conda "${projectDir}/env/shiver.yml"
  publishDir "${params.outdir}/05_shiver_map/${id}", mode: "copy", overwrite: true
  //debug true

  input:
    path initdir
    tuple val(id), path(kallistoRef), path(contigs), path(shiverRef), path(blast)
    tuple val(id), path(reads)
   
  
  output:
    tuple val("${id}"), path("${id}*ref.fasta"), path("${id}*.bam"), path("${id}*.bam.bai"), path("${id}*WithHXB2.csv")
    
  script:
    if ( shiverRef.size() > 0 && params.mode == "paired") {
    """
    shiver_map_reads.sh \
        ${initdir} \
        ${params.config_pe} \
        ${contigs} \
        ${id} \
        ${blast} \
        ${shiverRef} \
        ${reads[0]} \
        ${reads[1]}

    rm temp_* 
    rm *PreDedup.bam
    
    """ 
   } else if ( shiverRef.size() > 0 && params.mode == "single") {
    """
      shiver_map_reads.sh \
        ${initdir} \
        ${params.config_se} \
        ${contigs} \
        ${id} \
        ${blast} \
        ${shiverRef} \
        ${reads[0]} \
    
    rm temp_* 
    rm *PreDedup.bam
    
    """ 
    } else if ( shiverRef.size() <= 0 && params.mode == "paired") {
     """
    touch ${id}.blast
    touch ${id}_merged_contigs.fasta

    shiver_map_reads.sh \
        ${initdir} \
        ${params.config_pe} \
        ${id}_contigs.fasta \
        ${id} \
        ${id}.blast\
        ${kallistoRef} \
        ${reads[0]} \
        ${reads[1]}

    rm temp_* 
    rm *PreDedup.bam
     """
  } else if ( shiverRef.size() <= 0 && params.mode == "single") {
    """
    touch ${id}.blast
    touch ${id}_merged_contigs.fasta

    shiver_map_reads.sh \
        ${initdir} \
        ${params.config_se} \
        ${id}_contigs.fasta \
        ${id} \
        ${id}.blast\
        ${kallistoRef} \
        ${reads[0]} \
        ${reads[1]}

    rm temp_* 
    rm *PreDedup.bam
    """
  }
}


process MAF {
 conda "${projectDir}/env/phylo_tsi.yml"
 publishDir "${params.outdir}/06_maf", mode: "copy", overwrite: true
 //debug true

 input:
  tuple path(hxb2), val(id), path(ref), path(bam), path(bai), path(basefreqs)
    
 output:
  path "${id}_maf.csv"
    
 script:
  if (basefreqs instanceof List) {
  """
  calculate_maf.py -b ${basefreqs[1]} -r ${hxb2} -o ${id}_maf.csv
  """ 
  } else {
  """
  calculate_maf.py -b ${basefreqs} -r ${hxb2} -o ${id}_maf.csv
  """
   }
}

process JOIN_MAFS {
  conda "${projectDir}/env/phylo_tsi.yml"
  publishDir "${params.outdir}/07_joined_maf", mode: "copy", overwrite: true
  //debug true

  input:
    path maf_csvs
    
  output:
    path "*.csv"
    
  script:
    """
    awk 'FNR==1 && NR!=1 { next } { print }' ${maf_csvs} > joined_maf.csv
    
    """
  }


// PHYLOSCANNER PART

process BAM_REF_ID_CSV {
  //publishDir "${params.outdir}/ref_bam_id", mode: "copy", overwrite: true
  //debug true

  input:
    tuple val(id), path(ref), path(bam), path(bai), path(basefreqs)
    
  output:
    path "*_bam_ref_id.csv"
  
  script:
    if (bam instanceof List) {
    """
    for bamfile in *_remap.bam; do
      echo ${id}_remap.bam,${id}_remap_ref.fasta,${id}
    done > ${id}_bam_ref_id.csv
    """ 
  } else {
     """
    for bamfile in *.bam; do
      echo ${id}.bam,${id}_ref.fasta,${id}  
    done > ${id}_bam_ref_id.csv
     """
  }
}


process PHYLOSCANNER_ALIGN_READS {
 label "phyloscanner_align_reads"
 conda "${projectDir}/env/phyloscanner.yml"
 //publishDir "${params.outdir}/phyloscanner_aligned_reads", mode: "copy", overwrite: true 
 //debug true

 input:
  path bam_ref_id_csv, name: "phyloscanner_input.csv"
  path bam_fasta_bai_files

 output:
  path "AlignedReads/*.fasta", emit: AlignedReads
  path "Consensuses/*.fasta", emit: Consensuses
  path "ReadNames/*.csv.gz", emit: ReadsNames
  path "*.csv", emit: WindowCoordinateCorrespondence



 script:
  set_paired = params.mode == 'paired' ? '--merge-paired-reads' : ''
  // remove 9470,9720,9480,9730,9490,9740 from windows
 """
  phyloscanner_make_trees.py \
       ${bam_ref_id_csv} \
       ${set_paired} \
       --quality-trim-ends 25 \
       --alignment-of-other-refs ${params.two_refs} \
       --pairwise-align-to B.FR.83.HXB2_LAI_IIIB_BRU.K03455 \
       --excision-ref B.FR.83.HXB2_LAI_IIIB_BRU.K03455 \
       --excision-coords \$(cat ${params.excision_coordinates}) \
       --no-trees \
       --merging-threshold-a 0 \
       --min-read-count 1 \
       --windows \$(cat ${params.windows_oneline}) 
 """ 
}



process ALIGNED_READS_IQTREE {
  //publishDir "${params.outdir}/reads_all_windows", mode: "copy", overwrite: true
  //debug true

  input:
    tuple val(window), path(aligned_reads)
  
  output:
    tuple val (window), path ("Filtered*.fasta")
    
  script:
  
    if ( aligned_reads[1] ) {

    """
    mv ${aligned_reads[0]} Filtered${aligned_reads[0]}

    mv ${aligned_reads[1]} Filtered${aligned_reads[1]}
    rm Filtered${aligned_reads[0]}
   
    """
    } else {
      """
      mv ${aligned_reads[0]} Filtered${aligned_reads[0]}
      """
    }
  }


process IQTREE {
  label "iqtree"
  conda "${projectDir}/env/phyloscanner.yml"

 input:
  tuple val (window), path (fasta)

 output:
  path "*.treefile", emit: Treefile
  path "*.log", emit: Iqtreelog
 
 script:
 """
  iqtree \
     -s ${fasta} \
     -pre IQTREE_bestTree.InWindow_${window} \
     -m GTR+F+R6 \
     -nt 2 \
     --seed 1
 """ 
 // -nt ${task.cpus} \
 // -nt 2 (for full-length) and -nt 1 (for pol only)
}


process PHYLOSCANNER_R_INSTALL {
  conda "${projectDir}/env/phyloscannerR.yml"
  
  input:

  output:
   path "emptyfile.csv"

  script:
    """
    #!/usr/bin/env Rscript
  
    library(devtools)
    devtools::install_github("BDI-pathogens/phyloscanner/phyloscannerR", ref = "v1.9.2", dependencies=T)
    file.create("emptyfile.csv")
    """ 
}


process PHYLOSCANNER_TREE_ANALYSIS {
  conda "${projectDir}/env/phyloscannerR.yml"
  label "phyloscanner_tree_analysis"
  publishDir "${params.outdir}/09_patStats", mode: "copy", overwrite: true
  //debug true

 input:
  path emptyfile
  path treefile

 output:
   path "*patStats.csv", emit: patstat_csv
   path "*blacklistReport.csv", emit: blacklist_csv
   path "*patStats.pdf", emit: patstat_pdf
   path "nex_trees/*.nex", emit: nex
   //path "*.rda", emit: rda   

 script:
 """
  phyloscanner_analyse_trees.R \
    --skipSummaryGraph \
    --overwrite \
    --outputRDA \
    --outputNexusTree \
    --verbose 1 \
    --windowThreshold 0.5 \
    --allowMultiTrans \
    --directionThreshold 0.33 \
    --readCountsMatterOnZeroLengthBranches \
    --blacklistReport \
    --parsimonyBlacklistK ${params.k} \
    --ratioBlacklistThreshold 0.005 \
    --rawBlacklistThreshold 3 \
    --multifurcationThreshold 1E-5 \
    --outgroupName B.FR.83.HXB2_LAI_IIIB_BRU.K03455 \
    --normRefFileName ${params.hiv_distance_normalisation} \
    --treeFileExtension .treefile IQTREE_bestTree.InWindow "k${params.k}" "s,${params.k}" 

    mkdir nex_trees
    mv *.nex nex_trees
 """ 
}

process PHYLO_TSI {
  conda "${projectDir}/env/phylo_tsi.yml"
  publishDir "${params.outdir}/10_phylo_tsi/reports", mode: "copy", overwrite: true
  debug true

  input:
    path patstat
    path maf
    path regions
    
  output:
    path "phylo_tsi.csv"
  
  script:
    set_protocol = params.seqprotocol == 'amplicons' ? '--amplicons True' : '--amplicons False'

    if ( params.genome == "full" ) {
    """
    HIVPhyloTSI_full_genome.py \
      -d ${params.modelfolder} \
      -p ${patstat} \
      -m ${maf} \
      -o phylo_tsi.csv \
      ${set_protocol}
    """ 
  
  } else if ( params.genome == "partial" ) {
    """
    HIVPhyloTSI_partial_genome.py \
      -d ${params.modelfolder} \
      -p ${patstat} \
      -m ${maf} \
      -r ${regions} \
      -o phylo_tsi.csv \
      --modelname ${params.modelname} \
      ${set_protocol}
    """
  }
}

process MAPPING_NOTES {
  //debug true

  input:
    tuple val(id), path(kallistoRef), path(contigs), path(shiverRef), path(blast)
    
  output:
    path "${id}_mapping_notes.csv"
  
  script:
    if (contigs.size() > 0) {
    """
    echo ${id},"Mapped with SPADES and/or METASPADES contigs" > ${id}_mapping_notes.csv
    """ 
  } else {
     """
    bestref=\$(grep "^>" ${kallistoRef} | sed 's/>//g') 
    echo ${id},"Mapped with reference: \${bestref}" > ${id}_mapping_notes.csv
    """
  }
}

process MULTIQC_READS_REPORT {
  conda "${projectDir}/env/phylo_tsi.yml"
  publishDir "${params.outdir}/10_phylo_tsi/reports", mode: "copy", overwrite: true
  //debug true

  input:
    path multiqc_txt

  output:
    path "multiqc_report.csv"
   
  script:
    """
    multiqc_parser.py -i ${multiqc_txt} -o multiqc_report.csv 
    """ 
}

process REPORT {
  conda "${projectDir}/env/phylo_tsi.yml"
  publishDir "${params.outdir}/10_phylo_tsi", mode: "copy", overwrite: true
  //debug true

  input:
    path phylotsi_csv
    path multiqc_csv
    path mapping_csv

  output:
    path "hivtime.csv"
    path "hivtime.png"
  
  script:
    """
    hivtime_report.py -t ${phylotsi_csv} -q ${multiqc_csv} -m ${mapping_csv} -o hivtime
    """ 
}

// ****************************************************INPUT CHANNELS**********************************************************
ch_ref_hxb2 = Channel.fromPath("${projectDir}/data/refs/HXB2_refdata.csv", checkIfExists: true)
ch_regions = Channel.fromPath("${projectDir}/bin/models/partial/regions.csv", checkIfExists: true)


if (params.mode == 'paired') {
        ch_input_fastq = Channel
        .fromFilePairs( params.fastq, checkIfExists: true )
        .map { tuple ( it[0].split("HIV")[1].split("_")[0], [it[1][0], it[1][1]]) }

        
} else { ch_input_fastq = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .map { file -> [file.simpleName, [file]]}
        .map {tuple ( it[0].split("HIV")[1].split("_")[0], it[1][0])}

}


workflow {
   // ***********************************************************QC*********************************************************************
    ch_raw_fastqc = RAW_FASTQC ( ch_input_fastq )
    ch_fastp_trimmed = FASTP (  ch_input_fastq )
    ch_fastp_fastqc = FASTP_FASTQC ( ch_fastp_trimmed.Reads ) 

    if (params.seqprotocol == "amplicons") {
       ch_primer_trimmed = ALIENTRIMMER ( ch_fastp_trimmed.Reads, params.primers )
       ch_alientrimmer_fastqc = ALIENTRIMMER_FASTQC ( ch_primer_trimmed ) 
       ch_classified_reads = CLASSIFY ( ch_primer_trimmed, params.krakendb )
       ch_filtered_reads = EXTRACT ( ch_classified_reads.ClassifiedFastq, ch_classified_reads.KrakenOutput, ch_classified_reads.KrakenReport, params.taxid )
       ch_merged_reads = MERGE ( ch_classified_reads.UnclassifiedFastq.combine( ch_filtered_reads, by:0 ) )
       ch_kraken_fastqc = KRAKEN_FASTQC ( ch_merged_reads )    
       ch_multiqc = MULTIQC ( ch_raw_fastqc.Zip.concat(ch_fastp_fastqc.Zip).concat(ch_alientrimmer_fastqc.Zip).concat(ch_kraken_fastqc.Zip).collect() )
    } else if (params.seqprotocol == "capture") {
       ch_classified_reads = CLASSIFY ( ch_fastp_trimmed.Reads, params.krakendb )
       ch_filtered_reads = EXTRACT ( ch_classified_reads.ClassifiedFastq, ch_classified_reads.KrakenOutput, ch_classified_reads.KrakenReport, params.taxid )
       ch_merged_reads = MERGE ( ch_classified_reads.UnclassifiedFastq.combine( ch_filtered_reads, by:0 ) )
       ch_kraken_fastqc = KRAKEN_FASTQC ( ch_merged_reads )    
       ch_multiqc = MULTIQC ( ch_raw_fastqc.Zip.concat(ch_fastp_fastqc.Zip).concat(ch_kraken_fastqc.Zip).collect() )
    }
   
    // Creates csv file with rread length and amount of reads 
    ch_multiqc_report = MULTIQC_READS_REPORT ( ch_multiqc.Txt )
    // Contig generation
    ch_spades = SPADES ( ch_merged_reads )
    ch_metaspades = METASPADES ( ch_merged_reads )
    // Combine according to a key that is the first value of every first element, which is a list
    ch_spades_combined = ch_spades.SpadesContigs.combine( ch_metaspades.MetaspadesContigs, by:0 )
    ch_merged_contigs = MERGE_CONTIGS ( ch_spades_combined )
    ch_cd_hit_est = CD_HIT_EST ( ch_merged_contigs )
    // ******************************************************SHIVER*********************************************************************
    ch_fastq_renamed_header = FASTQ_RENAME_HEADER ( ch_merged_reads  )
    ch_initdir = SHIVER_INIT ( params.alignment, params.primers )
    ch_kallisto_index = KALLISTO_INDEX ( ch_initdir.ExistingRefsUngapped )
    ch_kallisto_index_reads = ch_kallisto_index.combine( ch_fastq_renamed_header )
    ch_kallisto_quant = KALLISTO_QUANT( ch_kallisto_index_reads )
    ch_best_ref = BEST_ALIGNMENT ( ch_initdir.IndividualRefs, ch_kallisto_quant )
    ch_wref = SHIVER_ALIGN ( ch_initdir.InitDir, ch_cd_hit_est )
    // Combine according to a key that is the first value of every first element, which is a list
    ch_mapping_args = ch_best_ref.combine(ch_cd_hit_est, by:0).combine(ch_wref, by:0).combine(ch_fastq_renamed_header, by:0)
    ch_mapping_args_non_reads = ch_mapping_args.map {id, bestref, contigs, shiverref, blast, reads  -> tuple (id, bestref, contigs, shiverref, blast)}
    ch_mapping_args_reads = ch_mapping_args.map {id, bestref, contigs, shiverref, blast, reads  -> tuple (id, reads)}
    ch_mapping_out = SHIVER_MAP ( ch_initdir.InitDir, ch_mapping_args_non_reads, ch_mapping_args_reads )
    // Mapping notes
    ch_mapping_notes = MAPPING_NOTES ( ch_mapping_args_non_reads )
    ch_mapping_notes_all = ch_mapping_notes.collectFile( name: "mapping_report.csv", storeDir: "${params.outdir}/10_phylo_tsi/reports" )
    //***********************************************************MAF********************************************************************
    ch_maf_out = MAF (ch_ref_hxb2.combine(ch_mapping_out ))
    ch_joined_maf = JOIN_MAFS ( ch_maf_out.collect() )
    // *******************************************************PHYLOSCANNER******'*****************************************************************
    ch_phyloscanner_csv = BAM_REF_ID_CSV ( ch_mapping_out )
    // An easy way to concatinate bam_ref_id_csv files: use collectFile() operator
    ch_bam_ref_id_all = ch_phyloscanner_csv.collectFile( name: "phloscanner_input.csv", storeDir: "${params.outdir}/08_bam_ref_id_all" )
    ch_mapped_out_no_id = ch_mapping_out.map {id, fasta, bam, bai, csv -> [fasta, bam, bai]}
    ch_aligned_reads = PHYLOSCANNER_ALIGN_READS ( ch_bam_ref_id_all, ch_mapped_out_no_id.flatten().collect() )

    // Combine aligned by phyloscanner reads
    ch_all_aligned_reads = ch_aligned_reads.AlignedReads.flatten().filter(~/^(.(?!(PositionsExcised)))*$/)
    ch_aligned_reads_positions_excised = ch_aligned_reads.AlignedReads.flatten().filter(~/.*PositionsExcised.*/)

    ch_window_all_reads = ch_all_aligned_reads
        .map { it -> tuple (it.getSimpleName().split("InWindow_")[1], it ) }

    ch_window_pos_exc_reads = ch_aligned_reads_positions_excised
        .map { it -> tuple (it.getSimpleName().split("_PositionsExcised_")[1], it) }
    
    ch_grouped_aligned_reads = ch_window_all_reads.concat( ch_window_pos_exc_reads )
        .groupTuple(remainder: true)

    ch_aligned_reads_iqtree = ALIGNED_READS_IQTREE ( ch_grouped_aligned_reads )
    ch_iqtree = IQTREE ( ch_aligned_reads_iqtree  )
    ch_install_phyloscannerR = PHYLOSCANNER_R_INSTALL ( )
    ch_analysed_trees = PHYLOSCANNER_TREE_ANALYSIS ( ch_install_phyloscannerR, ch_iqtree.Treefile.collect())
    // *******************************************************HIVPhyloTSI*****************************************************************
    ch_phylo_tsi = PHYLO_TSI ( ch_analysed_trees.patstat_csv, ch_joined_maf, ch_regions )
    // Report
    ch_hivtime_report = REPORT ( ch_phylo_tsi, ch_multiqc_report, ch_mapping_notes_all )
}
