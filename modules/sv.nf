process ngmlr {
  cpus 8
  tag "${sampleID} ngmlr mapping"

  input:
    path(read_file)
    path(reference)

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
  tuple(val(sampleID),file(bam))

  output:
  path("${sampleID}_sniffles.vcf")

  script:
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
  path("${sampleID}_summary.txt")

  script:
  sampleID=vcf.simpleName
  """
  SURVIVOR stats \
  ${vcf} \
  -1 \
  -1 \
  -1 ${sampleID}_summary_stats > \
  ${sampleID}_summary.txt
  """

}

process cutesv {
  publishDir "results/SV"
  tag "cuteSV"
  cpus 8

  input:
  tuple(val(sampleID),file(bam))
  path(reference)

  output:
  path("cuteSV.vcf")

  script:
  """
  cuteSV \
  -t ${task.cpus} \
  --max_cluster_bias_INS 100 \
  --diff_ratio_merging_INS 0.3 \
  --max_cluster_bias_DEL 100 \
  --diff_ratio_merging_DEL 0.3 \
  ${sampleID}.bam \
  ${reference} \
  cuteSV.vcf \
  .
  """

}

/*
The SVs are classified as PRECISE or IMPRECISE based on the breakpoint locations by sniffles.
We will only retain PRECISE variants called by sniffles:
*/
process filteringSV {
 publishDir "results/SV"

 input:
 path(vcf)

 output:
 path("${sampleID}_precise.vcf")

 script:
 sampleID=vcf.simpleName

 """
 grep -v "IMPRECISE" ${vcf} > ${sampleID}_precise.vcf
 """

}


/*
merge vcf files
*/
process survmerge {
 publishDir "results/SV"

 input:
 path(vcf1)
 path(vcf2)

 output:
 path("Survivor_merged.vcf")

 script:

 """
 ls ${vcf1} ${vcf2} > sample_files;
 SURVIVOR merge \
 sample_files \
 1000 2 1 1 1 30 \
 Survivor_merged.vcf
 """

}
