nextflow.enable.dsl=2

params.reads = "*.{bam,bai}"
params.outdir = "results"
params.reference = "/home/training/lr_genomic/Mapping/Av_TH0426.fna"

log.info """\
         LR - N F   P I P E L I N E
         ===================================
         bam           : ${params.bam}
         outdir          : ${params.outdir}
         reference       : ${params.reference}
         """
         .stripIndent()

include { nglmr;
  sniffles;
  survivor;
  cutesv;
  survivor as s2;
  survivor as s3;
  survivor as s4;
  } from './modules/sv.nf'

 
Channel
  .fromFilePairs(params.reads)
  .set{bam_ch}

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



workflow {
  sniffles(bam_ch)
  survivor(sniffles.out)
  cutesv(ngmlr.out,ref_ch)
  s2(cutesv.out)
  filteringSV(sniffles.out)
  s3(filteringSV.out)
  survmerge(cutesv.out,filteringSV.out)
  s4(survmerge.out)
}
