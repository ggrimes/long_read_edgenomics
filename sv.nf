nextflow.enable.dsl=2

params.bam = "*.{bam,bai}"
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
  filteringSV;
  survmerge;
  } from './modules/sv.nf'


Channel
  .fromFilePairs(params.bam) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
  .set{bam_ch}

Channel
  .fromPath(params.reference)
  .set{ref_ch}







workflow {
  sniffles(bam_ch)
  survivor(sniffles.out)
  cutesv(bam_ch,ref_ch)
  s2(cutesv.out)
  filteringSV(sniffles.out)
  s3(filteringSV.out)
  survmerge(cutesv.out,filteringSV.out)
  s4(survmerge.out)
}
