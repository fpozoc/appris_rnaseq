import os

configfile: "config/config.yaml"

def output_all():
    for genome in config["genomes"]:
        for ann in config["annotations"][genome]:
            if config["annotations"][genome][ann]["enabled"]:
                for dataset in config["samples"]:
                    for sample in config["samples"][dataset]:
                        for replicate in config["samples"][dataset][sample]:
                            #yield "out/{dataset}/{genome}/bwa/{sid}.{rep}.bam".format(genome=genome,sid=sample,rep=replicate,dataset=dataset)
                            yield [s.format(genome=genome,sid=sample,rep=replicate,dataset=dataset,ann=ann) for s in [
                                    "out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.out.bam",
                                   "out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.toTranscriptome.out.bam",
                                   "out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.out.sorted.bam.bai",
                                   "out/{dataset}/{genome}/featureCounts/{ann}/{sid}.{rep}.out",
                                   "out/{dataset}/{genome}/RSEM/{ann}/{sid}.{rep}",
                                   "out/{dataset}/{genome}/RSEM/{ann}/{sid}.{rep}/genes.results",
                                   "out/{dataset}/{genome}/RSEM/{ann}/{sid}.{rep}/isoforms.results",
                                ]
                            ]

rule all:
    input:
        output_all()

rule get_genomes:
    output:
       "out/common/genomes/{genome}/main.fasta"
    resources:
        mem_mb = 4000
    threads: 1
    params:
        url = lambda wc: config["genomes"][wc.genome]["url"]
    shell:"""
        curl {params.url} | gunzip > {output}
    """

rule get_annotation:
    output:
       main="out/common/genomes/{genome}/{ann}/main.gtf",
       cds="out/common/genomes/{genome}/{ann}/cds.gtf"
    resources:
        mem_mb = 4000
    threads: 1
    params:
        url = lambda wc: config["annotations"][wc.genome][wc.ann]["url"]
    shell:"""
        curl {params.url} | gunzip > {output.main}
        grep -P '\tCDS\t' {output.main} > {output.cds}
    """

rule get_chrsizes:
    input:
       fa= "out/common/genomes/{genome}/main.fasta"
    output:
       chr="out/common/genomes/{genome}/chrsizes"
    resources:
        mem_mb = 4000
    threads: 1
    conda: "envs/environment.yaml"
    shell:"""
        samtools faidx {input.fa} 
        cut -f1,2 {input.fa}.fai > {output}
    """

rule merge_fastqs:
    input:
        r=lambda wc: "seq/{}/{}".format(wc.dataset,config["samples"][wc.dataset][wc.sid][wc.rep][wc.r])
    output:
        merged="out/{dataset}/merge_fastqs/{sid}.{rep}_{r}.fastq.gz",
    log:
        stderr="log/{dataset}/merge_fastqs/{sid}.{rep}_{r}.err",
    threads: 1
    conda: "envs/environment.yaml"
    resources:
        mem_mb = 1000
    shell:
        "cat {input.r} > {output.merged} 2> {log.stderr}"

def input_cutadapt(wc,r):
    files = config["samples"][wc.dataset][wc.sid][wc.rep][r]
    if isinstance(files, str):
        return "seq/{}/{}".format(wc.dataset,files)
    else:
        return "out/{dataset}/merge_fastqs/{sid}.{rep}_{r}.fastq.gz".format(sid=wc.sid,rep=wc.rep,r=r,dataset=wc.dataset)

rule cutadapt_pe:
    input:
        r1=lambda wc: input_cutadapt(wc,"r1"),
        r2=lambda wc: input_cutadapt(wc,"r2")
    output:
        r1="out/{dataset}/cutadapt/{sid}.{rep}_r1.fastq.gz",
        r2="out/{dataset}/cutadapt/{sid}.{rep}_r2.fastq.gz"
    log:
        stdout="log/{dataset}/cutadapt/{sid}.{rep}.out",
        stderr="log/{dataset}/cutadapt/{sid}.{rep}.err"
    threads: 3
    conda: "envs/environment.yaml"
    resources:
        mem_mb = 16000,
        walltime = 1440
    shell:
       "cutadapt -j {threads} {config[params][cutadapt_pe]} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.stdout} 2> {log.stderr}"

idx_cmd = "STAR --runMode genomeGenerate --genomeDir {params.d} --genomeFastaFiles {input} > {log.stdout} 2> {log.stderr}"

rule STAR_idx:
    shadow:"shallow"
    input:
        "out/common/genomes/{genome}/main.fasta"
    output:
        f="out/common/indexes/star/{genome}/genomeParameters.txt"
    params:
        d="out/common/indexes/star/{genome}"
    threads: 1
    conda: "envs/environment.yaml"
    resources:
        mem_mb = 64000,
        threads = 16,
        walltime = 1440
    log:
        stdout="log/star_idx/{genome}.out",
        stderr="log/star_idx/{genome}.err"
    shell: idx_cmd

rule STAR:
    input:
        r1="out/{dataset}/cutadapt/{sid}.{rep}_r1.fastq.gz",
        r2="out/{dataset}/cutadapt/{sid}.{rep}_r2.fastq.gz",
        genome_idx_f="out/common/indexes/star/{genome}/genomeParameters.txt",
        annfile="out/common/genomes/{genome}/{ann}/main.gtf"
    output: 
        bam="out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.out.bam",
        tbam="out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.toTranscriptome.out.bam"
    threads: 12
    params:
        prefix="out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}",
        genome_idx = "out/common/indexes/star/{genome}"
    log:    
        stdout="log/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}.out",
        stderr="log/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}.err"
    conda: "envs/environment.yaml"
    resources:
        mem_mb = 128000,
        walltime = 90
    benchmark:
        "out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/_STAR.bmk"
    shell:         
        """
        STAR --sjdbGTFfile {input.annfile} --genomeDir {params.genome_idx} --outFileNamePrefix {params.prefix}/ --readFilesIn {input.r1} {input.r2} --runThreadN {threads} {config[params][STAR]} > {log.stdout} 2> {log.stderr}
    """

STAR --sjdbGTFfile {input.annfile} --genomeDir {params.genome_idx} --outFileNamePrefix {params.prefix}/ --readFilesIn {input.r1} {input.r2} --quantMode TranscriptomeSAM --runThreadN {threads} {config[params][STAR]} > {log.stdout} 2> {log.stderr}

rule samtools_sort:
   input:
       "out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.out.bam"
   output:
       "out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.out.sorted.bam"
   params:
       "-m 16G"
   threads:8
   resources:
       mem_mb = 32000,
       walltime = 1440
   wrapper:
       "0.30.0/bio/samtools/sort"

rule samtools_index:
   input:
       "out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.out.sorted.bam"
   output:
       "out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.out.sorted.bam.bai"
   params:
       "" # optional params string
   wrapper:
       "0.30.0/bio/samtools/index"

rule featurecounts:
   input:
      bam="out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.out.sorted.bam",
      ann="out/common/genomes/{genome}/{ann}/cds.gtf"
   output:
       "out/{dataset}/{genome}/featureCounts/{ann}/{sid}.{rep}.out"
   threads: 4
   params:
       mem_mb=8000,
       walltime = 1440
   log:    
       stdout="log/{dataset}/{genome}/featureCounts/{ann}/{sid}.{rep}.out",
       stderr="log/{dataset}/{genome}/featureCounts/{ann}/{sid}.{rep}.err"
   conda: "envs/environment.yaml"
   shell:"""
       featureCounts -O -M --fraction -T {threads} -t CDS -g gene_id -a {input.ann} -o {output} {input.bam} > {log.stdout} 2> {log.stderr}
   """

rule RSEM_prepare_reference:
   input:
       fasta="out/common/genomes/{genome}/main.fasta",
       gtf="out/common/genomes/{genome}/{ann}/main.gtf",
   output: 
       seq="out/common/indexes/RSEM/{genome}/reference.seq"
   params:
     genome_dir="out/common/indexes/RSEM/{genome}/{genome}"
   threads: 4
   resources:
       mem_mb=8000,
       walltime=1440
   log:
      stdout="log/RSEM_idx/{genome}.out",
       stderr="log/RSEM_idx/{genome}.out"
   conda: "envs/environment.yaml"
   shell:"""
       rsem-prepare-reference --gtf {input.gtf} --num-threads {threads} {input.fasta} {params.genome_dir} > {log.stdout} 2> {log.stderr}
   """

rule RSEM_calculate_expression:
   input:
       bam="out/{dataset}/{genome}/STAR/{ann}/{sid}.{rep}/Aligned.toTranscriptome.out.bam",
   params:
       prefix="out/{dataset}/{genome}/RSEM/{ann}/{sid}.{rep}/{sid}.{rep}",
       genome_idx="out/common/indexes/RSEM/{genome}/{genome}"
   output:
       genes="out/{dataset}/{genome}/RSEM/{ann}/{sid}.{rep}/genes.results",
       isoforms="out/{dataset}/{genome}/RSEM/{ann}/{sid}.{rep}/isoforms.results"
   threads: 16
   resources: 
       mem_mb=8000,
       walltime=1440
   log:
       stdout="log/{dataset}/{genome}/RSEM/{ann}/{sid}.{rep}.out",
       stderr="log/{dataset}/{genome}/RSEM/{ann}/{sid}.{rep}.err"
   conda: "envs/environment.yaml"
   shell:"""
       rsem-calculate-expression --num-threads {threads} --bam {input.bam} --no-bam-output --paired-end {params.genome_idx} {params.prefix} > {log.stdout} 2> {log.stderr}
   """ 
