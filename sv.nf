params.accession = 'SRR7449790'
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






/*
3.3 Mapping to the reference using NGMLR
We will use another long-read mapper: CoNvex Gap-cost alignMents for Long Reads (ngmlr)
to align reads to the reference genome. This is a long-read mapper designed to sensitively
align PacBio or Oxford Nanopore reads to (large) reference genomes including those spanning
complex structural variations.
*/

process ngmlr{
  cpus 8
  tag "${sampleID} ngmlr mapping"

  input:
    tuple(val(sampleID),path(read_file)) from read_ch
    path(reference) from ref_ch

  /*
  sort sam file and export to bam
  */
  output:
  path("${sampleID}_sorted.bam*") into sam_out

  script:
  """
  ngmlr \
  -t ${task.cpus} \
  -r ${reference} \
  -x ont \
  -q ${sampleID}.fastq.gz
  -o ${sampleID}.sam;

  samtools view \
  -b ${sampleID}.sam | \
  samtools sort - \
  > ${sampleID}_sorted.bam;

  samtools index \
  ${sampleID}_sorted.bam
  """

}

sam_out.view()
