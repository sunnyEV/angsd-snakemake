configfile: "config.yaml"

import os
from sys import stderr

# printing to stderr
def message(*args, **kwargs):
    print(*args, file=stderr, **kwargs)

# make sure the tmp directory exists
os.makedirs(config['tmpDir'], exist_ok=True)


## TODO: Add some informative messages to summarize configuration settings!


################################################################################
## Default Outputs
##  - this defines what running simply `snakemake` would generate
################################################################################
rule all:
 input:
   "qc/raw/fastqc_summary.html"

## Question: Why is the `input` key used? Why not `output`?

   
################################################################################
## Downloading Data
################################################################################

rule download_genome:
  output:
    fasta=config['genome']['fasta']
  params:
    url=config['genome']['url']
  shell:
    """
    wget -O {output.fasta} '{params.url}'
    """

rule download_txome:
  output:
    gtf=config['txome']['gtf']
  params:
    url=config['txome']['url']
  shell:
    """
    wget -O {output.gtf} '{params.url}'
    """

rule download_metadata:
  output:
    tsv=config['metadata']['tsv']
  params:
    url=config['metadata']['url']
  shell:
    """
    wget --content-disposition -O {output.tsv} '{params.url}'
    """

rule download_mapping:
  output:
    tsv=config['sample_map']['tsv']
  params:
    url=config['sample_map']['url']
  shell:
    """
    wget --content-disposition -O {output.tsv} '{params.url}'
    """

rule join_ftp_links:
  input:
    metadata=config['metadata']['tsv'],
    sample_map=config['sample_map']['tsv']
  output:
    "metadata/sample_ftp_map.tsv"
  shell:
    """
    join --header -t $'\t' {input.sample_map} {input.metadata} > {output}
    """

rule download_sample_fastqs:
  input:
    tsv="metadata/sample_ftp_map.tsv"
  output:
    expand("data/fastq/{condition}_{biorep}/{condition}_{biorep}_{lane}.fastq.gz",
           lane=list(config['lanes']), allow_missing=True)
  wildcard_constraints:
    condition="(WT|SNF2)",
    biorep="[0-9]+"
  shell:
    """
    path='data/fastq/{wildcards.condition}_{wildcards.biorep}/{wildcards.condition}_{wildcards.biorep}_'
    mkdir -p $(dirname $path)
    awk '$3 == "{wildcards.condition}" && $4 == "{wildcards.biorep}"' {input.tsv} |
    while read -r f; do
      lane=$(echo $f | cut -d' ' -f2)
      file=$(echo $f | cut -d' ' -f5)
      wget -O ${{path}}${{lane}}.fastq.gz "ftp://${{file}}"
    done
    """

################################################################################
## FastQC on replicates
################################################################################
rule fastqc:
  input:
    fastq="data/fastq/{condition}_{biorep}/{condition}_{biorep}_{lane}.fastq.gz"
  output:
    html="qc/fastqc/{condition}_{biorep}/{condition}_{biorep}_{lane}_fastqc.html",
    archive="qc/fastqc/{condition}_{biorep}/{condition}_{biorep}_{lane}_fastqc.zip"
  wildcard_constraints:
    condition="(WT|SNF2)",
    biorep="[0-9]+",
    lane="[0-9]+"
  conda: "envs/angsd.yaml"
  shell:
    """
    outdir=$(dirname {output.html})
    mkdir -p $outdir
    fastqc -o $outdir {input.fastq}
    """

rule multiqc_fastq:
  input:
    expand("qc/fastqc/{condition}_{biorep}/{condition}_{biorep}_{lane}_fastqc.html",
           condition=list(config['conditions']),
           biorep=list(config['biologicalReplicates']),
           lane=list(config['lanes']))
  output:
    html="qc/fastqc/fastqc_summary.html",
    data="qc/fastqc/fastqc_summary_data.zip"
  conda: "envs/angsd.yaml"
  shell:
    """
    multiqc -z -m fastqc -n fastqc_summary -o qc/raw qc/raw
    """

################################################################################
## Aligning data with STAR
################################################################################
rule star_genome_index:
  input:
    fasta=config['genome']['fasta'],
    gtf=config['txome']['gtf']
  output:
    idx="data/idx/sacCer3_STARindex/SAindex"
  params:
    tmpDir=lambda _: config["tmpDir"] + "/STAR_genome"
  threads: 1
  conda: "envs/angsd.yaml"
  shell:
    """
    mkdir -p $(dirname {params.tmpDir})
    mkdir -p data/idx
    STAR --runMode genomeGenerate \\
      --runThreadN {threads} \\
      --genomeDir $(dirname {output.idx}) \\
      --genomeFastaFiles {input.fasta} \\
      --sjdbGTFfile {input.gtf} \\
      --sjdbOverhang 50 \\
      --outTmpDir {params.tmpDir}
    rm -rf {params.tmpDir}
    """

rule star_align:
  input:
    fastq=expand("data/fastq/{condition}_{biorep}/{condition}_{biorep}_{lane}.fastq.gz",
                 lane=list(config['lanes']), allow_missing=True),
    idx="data/idx/sacCer3_STARindex/SAindex"
  output:
    "data/bam/STAR/{condition}_{biorep}.Aligned.sortedByCoord.out.bam"
  params:
    prefix=lambda wcs: "data/bam/STAR/%s_%s." % (wcs.condition, wcs.biorep),
    inputStr=lambda _, input: ",".join(input[:-1])
  threads: 1
  conda: "envs/angsd.yaml"
  ## TODO: add additional SAM tags
  shell:
    """
    mkdir -p $(dirname {output})
    STAR --runMode alignReads \\
      --runThreadN {threads} \\
      --genomeDir $(dirname {input.idx}) \\
      --readFilesIn {params.inputStr} \\
      --readFilesCommand zcat \\
      --outFileNamePrefix {params.prefix} \\
      --outSAMtype BAM SortedByCoordinate
    """

rule bam_index:
  input:
    bam=lambda wcs: "data/bam/%s/%s.bam" % (wcs.aligner, wcs.sample)
  output:
    "data/bam/{aligner}/{sample}.bam.bai"
  wildcard_constraints:
    aligner="(STAR|BWA)"
  conda: "envs/angsd.yaml"
  shell:
    """
    samtools index {input.bam}
    """

rule filter_chr:
  input:
    bam=lambda wcs: "data/bam/%s/%s.Aligned.sortedByCoord.out.bam" % (wcs.aligner, wcs.sample)
  output:
    bam="data/bam/{aligner}/{sample}.{chrom}.bam",
    bai="data/bam/{aligner}/{sample}.{chrom}.bam.bai"
  wildcard_constraints:
    aligner="(STAR|BWA)",
    chrom="chr(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI|M)"
  conda: "envs/angsd.yaml"
  shell:
    """
    samtools view -b -h {input.bam} {wildcards.chrom} > {output.bam}
    samtools index {output.bam}
    """
