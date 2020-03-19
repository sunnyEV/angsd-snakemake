configfile: "config.yaml"

rule all:
 input:
   "qc/raw/fastqc_summary.html",
   "qc/trimmed/WT_2/fastqc_compare_summary.html"

rule download_metadata:
  output:
    "metadata/PRJEB5348.txt"
  params:
    url=config['metadataURL']
  shell:
    """
    wget --content-disposition -O {output} '{params.url}'
    """

rule download_mapping:
  output:
    "metadata/ERP004763_sample_mapping.tsv"
  params:
    url=config['sampleMappingURL']
  shell:
    """
    wget --content-disposition -O {output} '{params.url}'
    """

rule join_ftp_links:
  input:
    sampleMap="metadata/ERP004763_sample_mapping.tsv",
    metadata="metadata/PRJEB5348.txt"
  output:
    "metadata/sample_ftp_map.tsv"
  shell:
    """
    join --header -t $'\t' {input.sampleMap} {input.metadata} > {output}
    """

rule download_sample_fastqs:
  input:
    "metadata/sample_ftp_map.tsv"
  output:
    expand("data/fastq/raw/{conditionio}_{brep}/{condition}_{biorep}_{lane}.fastq.gz",
           lane=[1,2,3,4,5,6,7], allow_missing=True)
  wildcard_constraints:
    condition="(WT|SNF2)",
    biorep="[0-9]+"
  shell:
    """
    path='data/fastq/raw/{wildcards.condition}_{wildcards.brep}/{wildcards.condition}_{wildcards.brep}_'
    mkdir -p $(dirname $path)
    awk '$3 == "{wildcards.condition}" && $4 == "{wildcards.biorep}"' {input} |
    while read -r f; do
      lane=$(echo $f | cut -d' ' -f2)
      file=$(echo $f | cut -d' ' -f5)
      wget -O ${{path}}${{lane}}.fastq.gz ftp://${{file}}
    done
    """

rule fastqc:
  input:
    "data/fastq/raw/{condition}_{biorep}/{condition}_{biorep}_{lane}.fastq.gz"
  output:
    html="qc/raw/{condition}_{biorep}/{condition}_{biorep}_{lane}_fastqc.html",
    archive="qc/raw/{condition}_{biorep}/{condition}_{biorep}_{lane}_fastqc.zip"
  wildcard_constraints:
    condition="(WT|SNF2)",
    biorep="[0-9]+",
    lane="[0-9]+"
  conda:
    "envs/angsd.yaml"
  shell:
    """
    outdir=$(dirname {output.html})
    mkdir -p $outdir
    fastqc -o $outdir {input}
    """

rule multiqc:
  input:
    expand("qc/raw/{condition}_{biorep}/{condition}_{biorep}_{lane}_fastqc.html",
           condition=list(config['conditions']),
           biorep=list(config['biologicalReplicates']),
           lane=[1,2,3,4,5,6,7])
  output:
    html="qc/raw/fastqc_summary.html",
    data="qc/raw/fastqc_summary_data.zip"
  conda:
    "envs/angsd.yaml"
  shell:
    """
    multiqc -z -m fastqc -n fastqc_summary -o qc/raw qc/raw
    """

rule trim_galore:
  input:
    "data/fastq/raw/{condition}_{biorep}/{condition}_{biorep}_{lane}.fastq.gz"
  output:
    fastq="data/fastq/trimmed/{condition}_{biorep}/{condition}_{biorep}_{lane}_trimmed.fq.gz",
    stats="data/fastq/trimmed/{condition}_{biorep}/{condition}_{biorep}_{lane}.fastq.gz_trimming_report.txt"
  wildcard_constraints:
    condition="(WT|SNF2)",
    biorep="[0-9]+",
    lane="[0-9]+"
  conda:
    "envs/angsd.yaml"
  shell:
    """
    trim_galore -o $(dirname {output.fastq}) {input}
    """

rule fastqc_trimmed:
  input:
    "data/fastq/trimmed/{condition}_{biorep}/{condition}_{biorep}_{lane}_trimmed.fq.gz"
  output:
    html="qc/trimmed/{condition}_{biorep}/{condition}_{biorep}_{lane}_trimmed_fastqc.html",
    archive="qc/trimmed/{condition}_{biorep}/{condition}_{biorep}_{lane}_trimmed_fastqc.zip"
  wildcard_constraints:
    condition="(WT|SNF2)",
    biorep="[0-9]+",
    lane="[0-9]+"
  conda:
    "envs/angsd.yaml"
  shell:
    """
    outdir=$(dirname {output.html})
    mkdir -p $outdir
    fastqc -o $outdir {input}
    """

rule multiqc_compare:
  input:
    expand("qc/raw/{condition}_{biorep}/{condition}_{biorep}_{lane}_fastqc.html",
           lane=[1,2,3,4,5,6,7], allow_missing=True),
    expand("qc/trimmed/{condition}_{biorep}/{condition}_{biorep}_{lane}_trimmed_fastqc.html",
           lane=[1,2,3,4,5,6,7], allow_missing=True)
  output:
    html="qc/trimmed/{condition}_{biorep}/fastqc_compare_summary.html",
    data="qc/trimmed/{condition}_{biorep}/fastqc_compare_summary_data.zip"
  params:
    rawDir=lambda wcs: "qc/raw/%s_%s" % (wcs.condition, wcs.biorep),
    trimDir=lambda wcs: "qc/trimmed/%s_%s" % (wcs.condition, wcs.biorep)
  conda:
    "envs/angsd.yaml"
  shell:
    """
    multiqc -z -m fastqc -n fastqc_compare_summary -o {params.trimDir} {params.rawDir} {params.trimDir}
    """

rule star_align:
  input:
    fastq=expand("data/fastq/raw/{condition}_{biorep}/{condition}_{biorep}_{lane}.fastq.gz",
                 lane=[1,2,3,4,5,6,7], allow_missing=True),
    idx="data/idx/sacCer3_STARindex"
  output:
    "data/bam/STAR/{condition}_{biorep}.Aligned.sortedByCoord.out.bam"
  params:
    prefix=lambda wcs: "data/bam/STAR/%s_%s." % (wcs.condition, wcs.biorep),
    inputStr=lambda _, input: ",".join(input[:-1])
  threads: 1
  conda:
    "envs/angsd.yaml"
  shell:
    """
    mkdir -p $(dirname {output})
    STAR --runMode alignReads \\
      --runThreadN {threads} \\
      --genomeDir {input.idx} \\
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

rule star_genome_index:
  input:
    fasta=config['genomeFASTA'],
    gtf=config['txomeGTF']
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
