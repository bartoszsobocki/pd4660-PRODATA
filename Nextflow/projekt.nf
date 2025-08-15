#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* ====== PARAMS (Twoje ścieżki) ====== */
params.input               = params.input        ?: '/root/PROJEKT_NEXTFLOW/data/*_{1,2}.fastq.gz'
params.outdir              = params.outdir       ?: '/root/PROJEKT_NEXTFLOW/results'
params.fastqc_threads      = (params.fastqc_threads      ?: 2)  as int
params.trimmomatic_threads = (params.trimmomatic_threads ?: 6)  as int
params.star_threads        = (params.star_threads        ?: 6)  as int
params.adapters            = params.adapters     ?: '/root/PROJEKT_NEXTFLOW/TruSeq3-PE.fa'
params.star_index          = params.star_index   ?: '/root/PROJEKT_NEXTFLOW/STAR_index'
params.multiqc_title       = params.multiqc_title ?: 'Analiza MultiQC'
params.ref                 = params.ref          ?: '/root/PROJEKT_NEXTFLOW/reference/NC_045512.2.fasta'
params.star_sort_ram = (params.star_sort_ram ?: 1200000000) as long


/* ====== FASTQC (RAW) ====== */
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
  fastqc -t ${task.cpus} -o . ${fastq_file}
  """
}

/* ====== MULTIQC (z FASTQC) ====== */
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
  multiqc . --title '${params.multiqc_title}'
  """
}

/* ====== TRIMMOMATIC (niezależnie od FASTQC/MULTIQC) ====== */
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
  set -euo pipefail
  trimmomatic PE -threads ${task.cpus} \
    ${r1} ${r2} \
    ${meta.id}_1_paired.fq.gz ${meta.id}_1_unpaired.fq.gz \
    ${meta.id}_2_paired.fq.gz ${meta.id}_2_unpaired.fq.gz \
    ${clip} SLIDINGWINDOW:4:20 MINLEN:36
  """
}

/* ====== STAR (mapowanie po TRIMMOMATIC) ====== */
process STAR_ALIGN {
  tag { meta.id }
  publishDir "${params.outdir}/aligned", mode: 'copy'
  cpus params.star_threads
  conda 'bioconda::star=2.7.11b'

  input:
  tuple val(meta), path(t1), path(t2)

  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}_Log.final.out")

  script:
  """
  set -euo pipefail

  # sprawdźmy czy index istnieje (genomeParameters.txt)
  if [[ ! -f "${params.star_index}/genomeParameters.txt" ]]; then
    echo "FATAL: Brak pliku ${params.star_index}/genomeParameters.txt (czy na pewno zbudowany STAR_index?)." >&2
    exit 1
  fi

  STAR --runThreadN ${task.cpus} \
       --genomeDir ${params.star_index} \
       --readFilesIn ${t1} ${t2} \
       --readFilesCommand zcat \
       --outFileNamePrefix ${meta.id}_ \
       --outSAMtype BAM SortedByCoordinate

  mv ${meta.id}_Aligned.sortedByCoord.out.bam ${meta.id}.bam
  """
}

/* ====== SAMTOOLS INDEX (indeksowanie BAM) ====== */
process INDEX_BAM {
  tag { meta.id }
  publishDir "${params.outdir}/aligned", mode: 'copy'
  conda 'bioconda::samtools=1.20'

  input:
  tuple val(meta), path(bam), path(log)

  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")

  script:
  """
  set -euo pipefail
  samtools index -@ 2 ${bam}
  ln -sf ${bam} ${meta.id}.bam
  ln -sf ${bam}.bai ${meta.id}.bam.bai
  """
}

/* ====== POKRYCIE (samtools depth) ====== */
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
  set -euo pipefail
  # pełne pokrycie per pozycja
  samtools depth -a ${bam} > ${meta.id}.depth.txt

  # prosta statystyka: średnie pokrycie i % pozycji z >=10x
  awk 'BEGIN{cov=0; n=0; ge10=0} {cov+=\$3; n++; if(\$3>=10) ge10++} END{if(n>0){print "mean_coverage\\t"cov/n; print "bases>=10x(%)\\t"100*ge10/n}else{print "mean_coverage\\t0"; print "bases>=10x(%)\\t0"}}' ${meta.id}.depth.txt > ${meta.id}.coverage_summary.txt
  """
}

/* ====== WARIANTY (bcftools mpileup + call) ====== */
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
  set -euo pipefail

  if [[ ! -f "${params.ref}" ]]; then
    echo "FATAL: Brak referencji: ${params.ref}" >&2
    exit 1
  fi

  bcftools mpileup -Ou -f ${params.ref} ${bam} \
    | bcftools call -mv -Oz -o ${meta.id}.vcf.gz

  bcftools index -t ${meta.id}.vcf.gz
  """
}

/* ================= WORKFLOW ================ */
workflow {

  log.info """
  ============================================
   FASTQC → MULTIQC   +   TRIMMOMATIC → STAR
   + index BAM → coverage → variants
  ============================================
  input      : ${params.input}
  outdir     : ${params.outdir}
  adapters   : ${params.adapters}
  STAR index : ${params.star_index}
  ref fasta  : ${params.ref}
  """

  /* Kanał do FASTQC (pojedyncze pliki) */
  raw_files_ch = Channel
    .fromPath(params.input, checkIfExists: true)
    .map { f ->
      def id = f.simpleName
              .replaceFirst(/_R?1$/, '')
              .replaceFirst(/_R?2$/, '')
              .replaceFirst(/\.fastq(?:\.gz)?$/, '')
      tuple([id: id], f)
    }

  /* FASTQC + MULTIQC */
  fastqc_out = FASTQC(raw_files_ch)
  zip_list   = fastqc_out.fastqc_zip.map{ it[1] }.collect()
  MULTIQC(zip_list)

  /* Pary do TRIMMOMATIC/STAR */
  pairs_ch = Channel
    .fromFilePairs(params.input, checkIfExists: true)
    .map { sid, files -> tuple([id: "${sid}"], files[0], files[1]) }

  trimmed   = TRIMMOMATIC(pairs_ch)
  aligned   = STAR_ALIGN(trimmed)
  indexed   = INDEX_BAM(aligned)
  COVERAGE(indexed)
  VARIANT_CALL(indexed)
}

/* ====== Podsumowanie ====== */
workflow.onComplete {
  println ""
  println "✅ Gotowe"
  println "  FASTQC    : ${params.outdir}/fastqc"
  println "  MULTIQC   : ${params.outdir}/multiqc/*_multiqc_report.html"
  println "  TRIM      : ${params.outdir}/trimmed"
  println "  STAR BAM  : ${params.outdir}/aligned/*.bam (+ .bai)"
  println "  COVERAGE  : ${params.outdir}/coverage/*.depth.txt, *.coverage_summary.txt"
  println "  VARIANTS  : ${params.outdir}/variants/*.vcf.gz(.tbi)"
}
