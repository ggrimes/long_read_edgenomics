nextflow.enable.dsl=2

params.reads = "$baseDir/data/*_{1,2}.fastq."
params.outdir = "results"
params.reference = "/home/training/lr_genomic/Mapping/Av_TH0426.fna"

log.info """\
         LR - N F   P I P E L I N E
         ===================================
         reads           : ${params.reads}
         outdir          : ${params.outdir}
         reference       : ${params.reference}
         """
         .stripIndent()



Channel
  .fromFilePairs(params.reads)
  .set{read_ch}

Channel
  .fromPath(params.reference)
  .set{ref_ch}



read_ch.view()


/*
3.3 Mapping to the reference using NGMLR
We will use another long-read mapper: CoNvex Gap-cost alignMents for Long Reads (ngmlr)
to align reads to the reference genome. This is a long-read mapper designed to sensitively
align PacBio or Oxford Nanopore reads to (large) reference genomes including those spanning
complex structural variations.
*/

process ngmlr {
  cpus 8
  tag "${sampleID} ngmlr mapping"

  input:
    path(read_file)
    path(reference)

  /*
  sort sam file and export to bam
  */
  output:
  path("${sampleID}_sorted.bam*")

  script:
  sampleID=read_file.simpleName
  """
  ngmlr \
  -t ${task.cpus} \
  -r ${reference} \
  -x ont \
  -q ${sampleID}.fastq.gz \
  -o ${sampleID}.sam;

  samtools view \
  -b ${sampleID}.sam | \
  samtools sort - \
  > ${sampleID}_sorted.bam;

  samtools index \
  ${sampleID}_sorted.bam
  """
}

/*
Calling Structural Variants
We will explore two SV callers: Sniffles and CuteSV. Sniffles is quite prominent in the field and
has been used widely to detect SVs from long reads generated from both PacBio and
Nanopore platforms. CuteSV, however is relatively new and is targeting Nanopore data in
particular. It is always a good idea to use more than one SV caller as the efficiency of the latter
depends on the data type along with the SV type and size and using more than one SVcaller
will give us more confidence for a call made.
*/
process sniffles {
  cpus 8
  publishDir "results/SV"
  tag "${sampleID} sniffles sv calling"

  input:
  path(bam)

  output:
  path("${sampleID}_sniffles.vcf")

  script:
  sampleID=read_file.simpleName
  """
  sniffles \
  -m ${sampleID}.bam \
  -v ${sampleID}_sniffles.vcf \
  -t 8
  """

}

/*
SURVIVOR stats sniffles.vcf -1 -1 -1 sniffles_summary_stats >
sniffles_summary.txt
*/
process survivor {
  publishDir "results/SV"
  input:
  path(vcf)

  output:
  path("sniffles_summary.txt")

  script:
  """
  SURVIVOR stats \
  ${vcf} \
  -1 \
  -1 \
  -1 sniffles_summary_stats >\
  sniffles_summary.txt
  """

}


workflow {
  ngmlr(read_ch,ref_ch)
  sniffles(ngmlr.out)
  survivor(sniffles.out)
}
