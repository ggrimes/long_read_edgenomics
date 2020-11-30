/*
 *  LR pipeline implemented with Nextflow derivied from edgenomic training
 *
 * Authors:
 * - Graeme Grimes <graeme.grimes@igmm.ed.ac.uk>
 */

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.accession = 'SRR7449790'
params.reads = "$baseDir/data/*_{1,2}.fq"
params.outdir = "results"

log.info """\
         LR - N F   P I P E L I N E
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()



 Channel
   .fromFilePairs( params.reads )
   .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
   .into { read_pairs_ch }


/*
Plotting tool for long read sequencing data and alignments.
https://github.com/wdecoster/NanoPlot
*/
process NanoPlot {
  cpus 4
  tag "${sample_id} Nanoplot"
  label "big_mem"
  publishDir "${params.outdir}/${sample_id}/QC",  mode:'copy'

  input:
    tuple sample_id, file(reads_file) from read_pairs_ch

  output:
    path("${sample_id}_nanoplot")

  script:
  """
  NanoPlot \
  -t ${task.cpus} \
  --fastq ${reads_file} \
  --loglength \
  -o ${sample_id}_nanoplot \
  --plots dot
  """
}
