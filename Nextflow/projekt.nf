#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.input               = params.input        ?: '/root/PROJEKT_NEXTFLOW/data/*_{1,2}.fastq.gz'
params.outdir              = params.outdir       ?: '/root/PROJEKT_NEXTFLOW/results'
params.fastqc_threads      = (params.fastqc_threads      ?: 2)  as int
params.trimmomatic_threads = (params.trimmomatic_threads ?: 6)  as int
params.star_threads        = (params.star_threads        ?: 4)  as int
params.adapters            = params.adapters     ?: '/root/PROJEKT_NEXTFLOW/TruSeq3-PE.fa'
params.star_index          = params.star_index   ?: '/root/PROJEKT_NEXTFLOW/STAR_index'
params.multiqc_title       = params.multiqc_title ?: 'Analiza MultiQC'
params.ref                 = params.ref          ?: '/root/PROJEKT_NEXTFLOW/reference/NC_045512.2.fasta'


process FASTQC {
  tag { meta.id }
  publishDir "${params.outdir}/fastqc", mode: 'copy'
  cpus params.fastqc_threads
  conda 'bioconda::fastqc=0.11.9'

  input:
  tuple val(meta), path(fastq_file)

  output:
  tuple val(meta), path('*_fastqc.zip'),  emit: fastqc_zip
  tuple val(meta), path('*_fastqc.html'), emit: fastqc_html

  script:
  """
  set -eo pipefail
  fastqc -t ${task.cpus} -o . ${fastq_file}
  """
}


process MULTIQC {
  tag "multiqc"
  publishDir "${params.outdir}/multiqc", mode: 'copy'
  conda 'bioconda::multiqc=1.14'

  input:
  path fastqc_zips

  output:
  path '*_multiqc_report.html',  emit: report
  path '*_multiqc_report_data',  emit: data

  script:
  """
  set -eo pipefail
  multiqc . --title '${params.multiqc_title}'
  """
}


process TRIMMOMATIC {
  tag { meta.id }
  publishDir "${params.outdir}/trimmed", mode: 'copy'
  cpus params.trimmomatic_threads
  conda 'bioconda::trimmomatic=0.39'

  input:
  tuple val(meta), path(r1), path(r2)

  output:
  tuple val(meta), path("${meta.id}_1_paired.fq.gz"), path("${meta.id}_2_paired.fq.gz")

  script:
  def clip = file(params.adapters).exists() ? "ILLUMINACLIP:${params.adapters}:2:30:10:8:true" : ""
  """
  set -eo pipefail
  trimmomatic PE -threads ${task.cpus} \
    ${r1} ${r2} \
    ${meta.id}_1_paired.fq.gz ${meta.id}_1_unpaired.fq.gz \
    ${meta.id}_2_paired.fq.gz ${meta.id}_2_unpaired.fq.gz \
    ${clip} SLIDINGWINDOW:4:20 MINLEN:36
  """
}


process STAR_ALIGN {
  tag { meta.id }
  publishDir "${params.outdir}/aligned", mode: 'copy'
  cpus params.star_threads
  conda 'bioconda::star=2.7.11b bioconda::samtools=1.17'

  input:
  tuple val(meta), path(t1), path(t2)

  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}_Log.final.out")

  script:
  """
  set -eo pipefail




  
  STAR --runThreadN ${task.cpus} \
       --genomeDir "${params.star_index}" \
       --readFilesIn ${t1} ${t2} \
       --readFilesCommand zcat \
       --outFileNamePrefix ${meta.id}_ \
       --outSAMtype BAM Unsorted


  samtools sort -@ ${task.cpus} -m 256M \
    -o ${meta.id}.bam ${meta.id}_Aligned.out.bam

  rm -f ${meta.id}_Aligned.out.bam
  """
}


process INDEX_BAM {
  tag { meta.id }
  publishDir "${params.outdir}/aligned", mode: 'copy'
  cpus 2
  conda 'bioconda::samtools=1.17'

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")

  script:
  """
  set -eo pipefail
  if [[ "\$(basename "${bam}")" != "${meta.id}.bam" ]]; then
    mv "${bam}" "${meta.id}.bam"
  fi
  samtools index -@ ${task.cpus} "${meta.id}.bam"
  """
}


process COVERAGE {
  tag { meta.id }
  publishDir "${params.outdir}/coverage", mode: 'copy'
  conda 'bioconda::samtools=1.20'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path("${meta.id}.depth.txt"), path("${meta.id}.coverage_summary.txt")

  script:
  """
  set -eo pipefail
  samtools depth -a ${bam} > ${meta.id}.depth.txt

  awk 'BEGIN{cov=0; n=0; ge10=0} {cov+=\$3; n++; if(\$3>=10) ge10++} END{
    if(n>0){
      print "mean_coverage\\t"cov/n;
      print "bases>=10x(%)\\t"100*ge10/n
    } else {
      print "mean_coverage\\t0";
      print "bases>=10x(%)\\t0"
    }}' ${meta.id}.depth.txt > ${meta.id}.coverage_summary.txt
  """
}

process VARIANT_CALL {
  tag { meta.id }
  publishDir "${params.outdir}/variants", mode: 'copy'
  conda 'bioconda::bcftools=1.20 bioconda::samtools=1.20'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi")

  script:
  """
  set -eo pipefail
  if [[ ! -f "${params.ref}" ]]; then
    echo "FATAL: Brak referencji: ${params.ref}" >&2
    exit 1
  fi

  bcftools mpileup -Ou -f ${params.ref} ${bam} \
    | bcftools call -mv -Oz -o ${meta.id}.vcf.gz

  bcftools index -t ${meta.id}.vcf.gz
  """
}


process MERGE_VARIANTS {
  publishDir "${params.outdir}/variants", mode: 'copy'
  conda 'bioconda::bcftools=1.20'

  input:
  path vcfs
  path tbis

  output:
  path "merged_variants.vcf.gz"
  path "merged_variants.vcf.gz.tbi"

  script:

  def vcf_list = (vcfs instanceof List ? vcfs : [vcfs]).collect{ it.getName() }.join(' ')
  """
  set -eo pipefail
  bcftools merge -Oz -o merged_variants.vcf.gz ${vcf_list}
  bcftools index -t merged_variants.vcf.gz
  """
}


workflow {
  log.info """
  ============================================
   FASTQC → MULTIQC   +   TRIMMOMATIC → STAR
   + index BAM → coverage → variants → MERGE
  ============================================
  input      : ${params.input}
  outdir     : ${params.outdir}
  adapters   : ${params.adapters}
  STAR index : ${params.star_index}
  ref fasta  : ${params.ref}
  """

  raw_files_ch = Channel
    .fromPath(params.input, checkIfExists: true)
    .map { f ->
      def id = f.simpleName
               .replaceFirst(/_R?1$/, '')
               .replaceFirst(/_R?2$/, '')
               .replaceFirst(/\.fastq(?:\.gz)?$/, '')
      tuple([id:id], f)
    }

  fastqc_out = FASTQC(raw_files_ch)
  zip_list   = fastqc_out.fastqc_zip.map{ it[1] }.collect()
  MULTIQC(zip_list)


  pairs_ch = Channel
    .fromFilePairs(params.input, checkIfExists: true)
    .map { sid, files -> tuple([id: "${sid}"], files[0], files[1]) }

  trimmed         = TRIMMOMATIC(pairs_ch)
  aligned         = STAR_ALIGN(trimmed)                           // (meta, bam, log)
  aligned_for_idx = aligned.map { meta, bam, log -> tuple(meta, bam) } // zrzucamy log
  indexed         = INDEX_BAM(aligned_for_idx)                    // (meta, bam, bai)

  COVERAGE(indexed)

  variants  = VARIANT_CALL(indexed)                               // (meta, vcf, tbi)
  all_vcfs  = variants.map { meta, vcf, tbi -> vcf }.collect()
  all_tbis  = variants.map { meta, vcf, tbi -> tbi }.collect()
  MERGE_VARIANTS(all_vcfs, all_tbis)
}


workflow.onComplete {
  println ""
  println "✅ Gotowe"
  println "  FASTQC    : ${params.outdir}/fastqc"
  println "  MULTIQC   : ${params.outdir}/multiqc/*_multiqc_report.html"
  println "  TRIM      : ${params.outdir}/trimmed"
  println "  STAR BAM  : ${params.outdir}/aligned/*.bam (+ .bai)"
  println "  COVERAGE  : ${params.outdir}/coverage/*.depth.txt, *.coverage_summary.txt"
  println "  VARIANTS  : ${params.outdir}/variants/*.vcf.gz(.tbi)"
  println "  MERGED    : ${params.outdir}/variants/merged_variants.vcf.gz(.tbi)"
}

  println "  VARIANTS  : ${params.outdir}/variants/*.vcf.gz(.tbi)"
}
