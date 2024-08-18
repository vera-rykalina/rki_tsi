nextflow.enable.dsl = 2

// Run

/* nextflow shiver_phyloscanner_tsi_se_pe.nf \
-c rki_tsi_profile.config 
--outdir output_pe
--fastq "/scratch/rykalinav/rki_recency/Pipeline/RawData/*R{1,2}*.fastq.gz" 
-profile rki_slurm,rki_mamba \
--krakendb /scratch/databases/kraken2_nt_20231129/ \
--mode paired 
--alignment /scratch/rykalinav/rki_tsi/data/HIV1_COM_2022_genome_DNA.fasta \
--primers /scratch/rykalinav/rki_tsi/data/primers/primers_GallEtAl2012.fasta \
-resume
*/


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
        tuple val("${id}"), path("${id}_kraken.R*.fastq.gz")

    script:
    if (params.mode == "paired") {
        """
        gzip -c ${unclassified[0]} ${filtered[0]} > ${id}_kraken.R1.fastq.gz
        gzip -c ${unclassified[1]} ${filtered[1]} > ${id}_kraken.R2.fastq.gz
        """
    } else if (params.mode == "single") {
       """
       gzip -c ${unclassified[0]} ${filtered[0]} > ${id}_kraken.R1.fastq.gz
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

  rm -f ${id}_alientrimmer.R.S.fastq.gz
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
  publishDir "${params.outdir}/14_merged_contigs", mode: "copy", overwrite: true
  
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
       -T ${task.cpus} || touch ${id}_clustered_contigs.fasta
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
  
  if (params.mode == "paired") {
  """
  shiver_init.sh \
    InitDir \
    ${params.config_pe} \
    ${alignment} \
    ${params.illumina_adapters} \
    ${primers}
  """  
  } else if (params.mode == "single") {

  """
  shiver_init.sh \
    InitDir \
    ${params.config_se} \
    ${alignment} \
    ${params.illumina_adapters} \
    ${primers}
  """ 
  }
}


process FASTQ_RENAME_HEADER {
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


process SHIVER_ALIGN {
  label "shiver"
  conda "${projectDir}/env/shiver.yml"
  publishDir "${params.outdir}/21_alignments/${id}", mode: "copy", overwrite: true
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
  publishDir "${params.outdir}/22_mapped/${id}", mode: "copy", overwrite: true
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
    } else if ( shiverRef.size() >= 0 && params.mode == "paired") {
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
  } else if ( shiverRef.size() >= 0 && params.mode == "single") {
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
 conda "${projectDir}/env/tsi-python.yml"
 publishDir "${params.outdir}/23_maf", mode: "copy", overwrite: true
 //debug true

 input:
  tuple val(id), path(ref), path(bam), path(bai), path(basefreqs)
    
 output:
  path "${id}.csv"
    
 script:
  if (basefreqs instanceof List) {
  """
  produce_maf.py ${basefreqs[1]} ${id}.csv
  """ 
  } else {
  """
  produce_maf.py ${basefreqs} ${id}.csv
  """
   }
}

process JOIN_MAFS {
  conda "${projectDir}/env/tsi-python.yml"
  publishDir "${params.outdir}/24_joined_maf", mode: "copy", overwrite: true
  //debug true

  input:
    path mafcsv
    
  output:
    path "*.csv"
    
  script:
    """
    join_mafs.py ${mafcsv}
    """ 
  }


// PHYLOSCANNER PART

process BAM_REF_ID_CSV {
  publishDir "${params.outdir}/25_ref_bam_id", mode: "copy", overwrite: true
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
 publishDir "${params.outdir}/27_phyloscanner_aligned_reads", mode: "copy", overwrite: true 
 debug true

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
       --dont-check-duplicates \
       --read-names-only \
       --merging-threshold-a 0 \
       --min-read-count 1 \
       --windows \$(cat ${params.windows_oneline})
 
 """ 
}

process IQTREE {
  label "iqtree"
  conda "${projectDir}/env/iqtree.yml"
  publishDir "${params.outdir}/28_iqtree_trees", mode: "copy", overwrite: true
  //debug true

 input:
  path fasta

 output:
  path "*.treefile", emit: treefile
  path "*.log", emit: iqtreelog
 
 script:
 """
  iqtree \
     -s ${fasta} \
     -pre IQTREE_bestTree.InWindow_${fasta.getSimpleName().split("Excised_")[1]} \
     -m GTR+F+R6 \
     -nt ${task.cpus} \
     --seed 1
 """ 
}


process PHYLOSCANNER_TREE_ANALYSIS {
 label "phyloscanner_tree_analysis"
 publishDir "${params.outdir}/29_analysed_trees", mode: "copy", overwrite: true
 debug true

 input:
  path treefile

 output:
   path "*patStats.csv", emit: patstat_csv
   path "*blacklistReport.csv", emit: blacklist_csv
   path "*patStats.pdf", emit: patstat_pdf
   //path "*.rda", emit: rda   
   //path "*.nex", emit: nex

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
 """ 
}

process PHYLO_TSI {
  conda "${projectDir}/env/phylo_tsi.yml"
  publishDir "${params.outdir}/30_phylo_tsi", mode: "copy", overwrite: true
  debug true

  input:
    path patstat
    path maf
    
  output:
    path "phylo_tsi.csv"
  
  script:
    """
    HIVPhyloTSI.py \
      -d ${params.model} \
      -p ${patstat} \
      -m ${maf} \
      -o phylo_tsi.csv \
      --amplicons True
    """ 
}

process PRETTIFY_AND_PLOT {
  conda "${projectDir}/env/tsi-python.yml"
  publishDir "${params.outdir}/31_phylo_tsi", mode: "copy", overwrite: true
  debug true

  input:
    path phylo_tsi_csv

  output:
    path "phylo_tsi_prettified.csv"
    path "tsi_barplot.png"
  
  script:
    """
    tsi_prettify_and_plot.py ${phylo_tsi_csv} 
    """ 
}

process MAPPING_NOTES {
  debug true

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


process PHYLOSCANNER_NORMALISATION {
  label "normalisation"
  conda "${projectDir}/env/phyloscanner.yml"
  debug true

  input:
    path (alignment)
    
  output:
    path "*.csv"
  
  script:
   
    """
    CalculateTreeSizeInGenomeWindows.py \
    ${alignment} \
    B.FR.83.HXB2_LAI_IIIB_BRU.K03455 \
    500 \
    251 \
    HIV_COM_2022_genome_DNA_SizeInWindows \
    --x-iqtree "${params.iqtreeargs}" \
    --end 9460 \
    -T ${task.cpus} \
    -I 10
    """
  
}
//**************************************************PARAMETERS*******************************************************
// Change is required! Specify your projectDir here
projectDir = "/scratch/rykalinav/rki_tsi"

// Parameters for kraken
krakendb = params.krakendb
// taxid of HIV-1 
params.taxid = "11676"


// Parameters for shiver
params.alientrimmer = "${projectDir}/bin/AlienTrimmer.jar"
params.illumina_adapters = "${projectDir}/data/adapters/adapters_Illumina.fasta"
params.config_se = "${projectDir}/bin/config_se.sh"
params.config_pe = "${projectDir}/bin/config_pe.sh"
params.alignment = "${projectDir}/data/alignments/HIV1_COM_2022_genome_DNA.fasta"
primers = params.primers



// Parameters for phyloscanner
params.raxmlargs = "raxmlHPC-SSE3 -m GTRCAT -p 1 --no-seq-check"
params.iqtreeargs = "iqtree -m GTR+F+R6 --seed 1"
params.two_refs = "${projectDir}/data/refs/2refs_HXB2_C.BW.fasta"
params.excision_coordinates = "${projectDir}/data/phyloscanner/DrugResistancePositionsInHXB2.txt"
params.windows_oneline = "${projectDir}/data/phyloscanner/windows250_VR_norms_oneline.txt"
params.hiv_distance_normalisation = "${projectDir}/data/phyloscanner/HIV_DistanceNormalisationOverGenome.csv"
params.k = 15

// Parameters for HIV-PhyloTSI
params.model = "${projectDir}/bin/Model"


// ****************************************************INPUT CHANNELS**********************************************************
ch_ref_hxb2 = Channel.fromPath("${projectDir}/data/refs/HXB2_refdata.csv", checkIfExists: true)


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
   // ***********************************************************QC*********************************************************************
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
    // Combine according to a key that is the first value of every first element, which is a list
    ch_spades_combined = ch_spades.SpadesContigs.combine( ch_metaspades.MetaspadesContigs, by:0 )
    ch_merged_contigs = MERGE_CONTIGS ( ch_spades_combined )
    ch_cd_hit_est = CD_HIT_EST ( ch_merged_contigs )
    ch_fastq_renamed_header = FASTQ_RENAME_HEADER ( ch_primer_trimmed )
     // ******************************************************SHIVER*********************************************************************
    ch_initdir = SHIVER_INIT ( params.alignment, primers )
    ch_kallisto_index = KALLISTO_INDEX ( ch_initdir.ExistingRefsUngapped )
    ch_kallisto_index_reads = ch_kallisto_index.combine( ch_fastq_renamed_header )
    ch_kallisto_quant = KALLISTO_QUANT( ch_kallisto_index_reads )
    ch_best_ref = BEST_ALIGNMENT ( ch_initdir.IndividualRefs, ch_kallisto_quant )
    ch_wref = SHIVER_ALIGN ( ch_initdir.InitDir, ch_merged_contigs )
    // Combine according to a key that is the first value of every first element, which is a list
    ch_mapping_args = ch_best_ref.combine(ch_merged_contigs, by:0).combine(ch_wref, by:0).combine(ch_fastq_renamed_header, by:0)
    ch_mapping_args_non_reads = ch_mapping_args.map {id, bestref, contigs, shiverref, blast, reads  -> tuple (id, bestref, contigs, shiverref, blast)}.view()
    ch_mapping_args_reads = ch_mapping_args.map {id, bestref, contigs, shiverref, blast, reads  -> tuple (id, reads)}.view()
    ch_mapping_out = SHIVER_MAP ( ch_initdir.InitDir, ch_mapping_args_non_reads, ch_mapping_args_reads )
     // *********************************************************MAF*********************************************************************
    ch_maf_out = MAF ( ch_mapping_out )
    ch_hxb2_maf = ch_ref_hxb2.combine(ch_maf_out.collect())
    ch_joined_maf = JOIN_MAFS ( ch_hxb2_maf )
    // *******************************************************PHYLOSCANNER******'*****************************************************************
    ch_phyloscanner_csv = BAM_REF_ID_CSV ( ch_mapping_out )
    // An easy way to concatinate bam_ref_id_csv files: use collectFile() operator
    ch_bam_ref_id_all = ch_phyloscanner_csv.collectFile( name: "phloscanner_input.csv", storeDir: "${projectDir}/${params.outdir}/26_bam_ref_id_all" )
    ch_mapped_out_no_id = ch_mapping_out.map {id, fasta, bam, bai, csv -> [fasta, bam, bai]}
    ch_aligned_reads = PHYLOSCANNER_ALIGN_READS ( ch_bam_ref_id_all, ch_mapped_out_no_id.flatten().collect() )
    ch_aligned_reads_positions_excised = ch_aligned_reads.AlignedReads.flatten().filter(~/.*PositionsExcised.*/)
    ch_iqtree = IQTREE ( ch_aligned_reads_positions_excised )
    ch_analysed_trees = PHYLOSCANNER_TREE_ANALYSIS ( ch_iqtree.treefile.collect() )
    ch_phylo_tsi = PHYLO_TSI( ch_analysed_trees.patstat_csv, ch_joined_maf )
    ch_prettified_tsi = PRETTIFY_AND_PLOT( ch_phylo_tsi )
     // Mapping notes
    ch_mapping_notes = MAPPING_NOTES( ch_mapping_args_non_reads )
    ch_mapping_notes_all = ch_mapping_notes.collectFile(name: "mapping_report.csv", storeDir: "${projectDir}/${params.outdir}/31_phylo_tsi")

    // Normalisation for phyloscanner (for option --normRefFileName)
    ch_ref_normalisation = PHYLOSCANNER_NORMALISATION ( params.alignment ) 
}



// fastaq primer trimming
//fastaq sequence_trim 07-00462_fastp.R1.fastq.gz 07-00462_fastp.R2.fastq.gz 07-00462_fastp_trimmed.R1.fastq.gz \
//07-00462_fastp_trimmed.R2.fastq.gz../../../DataShiverInit/primers_GallEtAl2012.fasta --revcomp

// cd-hit-est -i reads.fa -o output.fa -c 0.9 -n 10 -d 999

/*
-c <infile>
This option allows the alien sequence file to be indicated. Each alien oligonucleotide sequence must be
written in one line, and may not exceed 32,500 nucleotides. Standard degenerate bases are admitted, i.e.
character states M, R, W, S, Y, K, B, D, H, V, N, and X. Lines beginning by the characters ‘#’, ‘%’ or
‘>’ are not considered. Input file name may not be a number (see last page). */

/*  ./../tools/CalculateTreeSizeInGenomeWindows.py \
../../data/alignments/HIV1_COM_2022_genome_DNA.fasta \ 
B.FR.83.HXB2_LAI_IIIB_BRU.K03455 \
500 \
250 \
HIV_COM_2022_genome_DNA_SizeInWindows \
--x-iqtree "iqtree -m GTR+F+R6 -nt 2 --seed 0" \ 
--end 9460 \ 
-T 2 \ 
-I 10
*/