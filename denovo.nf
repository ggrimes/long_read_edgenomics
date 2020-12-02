
params.reads="R10_q9l500.fastq.gz"
params.genomeSize="5M"
params.conf="shasta.conf"
params.flyout="flye"
params.shastaout="shasta"
params.quastout="quast_flye_shasta"
params.reference="reference/K12_MG1655.fna"
params.busco_dir="busco"

log.info """\
         LR_ASSEMBLY - N F   P I P E L I N E
         ===================================
         reads           : ${params.reads}
         genomeSize      : ${params.genomeSize}
         reference       : ${params.reference}
         flyout          : results/${params.flyout}
         shastaout       : results/${params.shastaout}
}
         """
         .stripIndent()

/*
We will be using the FASTQ data for Escherichia coli K12MG1655 sequenced on the latest
R10.4 flowcell on ONT MinION platform. The raw fast5 data can be downloaded from ENA
under the study PRJEB36648. This is a test dataset generated by Dr. Rasmus H. Kirkegaard
from Albersten lab (https://albertsenlab.org/).
We have already generated the trimmed data using NanoFilt (Average quality:9 and
Length:500) and the file is made available. We have also got E.coli K12 MG1655 reference
genome and ODB database for gammaproteobacteria and other polishing models, which will
be used at a later stage.
*/


/*
Flye is one of the most widely used tools for assembling genomes from long noisy reads.
Shasta is a relatively new assembler but is extremely fast in generating assemblies. Both of
them are designed for a wide range of datasets from small bacteria to large mammals. It is
always a good idea to generate the assemblies with more than one assembler, as each one of
them might be better/worse in different aspects.

*/


Channel
  .fromPath(params.reads)
  .set{reads}

Channel
    .fromPath(params.conf)
    .set{conf}

Channel
     .fromPath(params.reference)
     .set{reference}

Channel
      .fromPath(params.busco_dir)
      .set{busco}


process flye {
  tag "flye assembly ${sampleID}"
  label "highmem"
  cpus 8
  publishDir "results", mode: 'copy'

  input:
  path(reads)

  output:
  path("${params.flyout}") into flyeout,flyeout2

  script:
  sampleID=reads.simpleName
  """
  flye \
  --nano-raw  ${reads}\
  -g ${params.genomeSize} \
  -o ${params.flyout} \
  -t ${task.cpus}
  """

}

/*
The Shasta assembler provides a lot of output files that allows exploring details of many data
structures used during assembly. The assembled contig/scaffold sequences are provided in
the file “Assembly.fasta”.
*/

Channel
  .fromPath(params.reads)
  .set{reads}

process shasta {
  cpus 8
  publishDir "results", mode: 'copy'

  input:
  path(reads)
  path(conf)

  output:
  path("shasta") into shastaout,shastaout2

  script:
  sampleID=reads.simpleName
  """
  #shasta need unzip reads
  gunzip -c ${reads} > ${sampleID}.fastq

  shasta \
  --input ${sampleID}.fastq \
  --assemblyDirectory ${params.shastaout} \
  --threads  ${task.cpus}  \
  --config ${conf}
  """
}


/*
4.4 Assessment of assemblies
4.4.1 Using QUAST
QUAST is a quality assessment tool for evaluating and comparing genome assemblies.
QUAST can evaluate assemblies both with a reference genome, as well as without a reference.
*/


process quast {
cpus 8
publishDir "results", mode: 'copy'

  input:
  path(flyeout)
  path(shastaout)
  path(reference)

  output:
  path("${params.quastout}")

  script:
  """
  quast \
  -o ${params.quastout} \
  -t ${task.cpus} \
  --labels flye,shasta ${flyeout}/assembly.fasta ${shastaout}/Assembly.fasta \
  -r ${reference}
  """
}



/*
BUSCO:
BUSCO provides a rich source of data to assess the quality and completeness of
genome assemblies, gene annotations, and transcriptomes by comparing them to
OrthoDB’s sets of Benchmarking Universal Single-Copy Orthologs.
We will run this on both assemblies generated so far.
*/
flyeout2.concat(shastaout2)

process busco {
  cpus 8
  publishDir "results", mode: 'copy'

  input:
  path(assembly) from flyeout2
  path(busco)

  output:
  path("${sampleID}_busco")

  script:
  sampleID=assembly.simpleName
  """
  busco \
  -i ${assembly}/*.fasta \
  -o ${sampleID}_busco \
  -m
  genome \
  -l gammaproteobacteria_odb10/ \
  -f \
  -c ${task.cpus}
  """
}
