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
map using minimap2 (Kmer minimizers)
*/
process mapping {
  cpus 8
  tag "${sampleID} mapping"

  input:
    tuple(val(sampleID),path(read_file)) from read_ch
    path(reference) from ref_ch

  /*
  sort sam file and export to bam
  */
  output:
  path("${sampleID}_mm2.sam") into sam_out

script:
  """
  minimap2 \
  -t ${task.cpus}  \
  -x map-ont \
  -a ${reference} \
  ${read_file} \
  > ${sampleID}_mm2.sam
  """
}


process sam2bam {

publishDir "${params.outdir}/${sampleID}",  mode:'copy'


input:
path(sam) from sam_out

output:
path(${sampleID}_sorted.bam) into sorted_bam_out


script:
sampleID=sam.simpleName
"""
samtools view \
-b -\
S ${sam} |\
 samtools sort - \
> ${sampleID}_sorted.bam
"""
}

/*
process bam_index {

input:
path(bam) from bam_index_out

output:
path("*.bam")

script:
"""
samtools index ${bam}
"""
}


process flagstat{

publishDir "results/QC/stats"

input:
   tuple(val(sampleID),path(files)) from bam_index_out

output:
  path("stat.txt")

script:
"""
samtools flagstat \
SRR7449790_Av_CIP_mm2_sorted.bam > stats.txt
"""
}


process NanoPlot{
  cpus 8
  tag "${sample_id} Nanoplot"
  label "big_mem"
  publishDir "${params.outdir}/${sample_id}/QC",  mode:'copy'

  input:
    path(bam) from bam_srt_ch

  output:
    path("${sample_id}_nanoplot")
    path("${reads_file}") into read_out_ch


  script:
  sample_id =bam.simpleName
  """
  NanoPlot –t 8 \
  --bam SRR7449790_Av_CIP_mm2_sorted.bam \
  --loglength
  -o nanoplot_bam \
  --plots dot
  """
}
*/
