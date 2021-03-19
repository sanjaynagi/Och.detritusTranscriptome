# An example collection of Snakemake rules imported in the main Snakefile.


rule FastQC:
    input:
        "resources/reads/{sample}_{n}.fq.gz"
    output:
        html="resources/reads/qc/{sample}_{n}_fastqc.html",
        zip="resources/reads/qc/{sample}_{n}_fastqc.zip"
    log:
        "logs/FastQC/{sample}_{n}_QC.log"
    params:
        outdir="--outdir resources/reads/qc"
    wrapper:
        "0.72.0/bio/fastqc"


rule trinityAssembly:
    input:
        left=expand("resources/reads/{sample}_1.fq.gz", sample=samples),
        right=expand("resources/reads/{sample}_2.fq.gz", sample=samples)
    output:
        "results/trinity_out_dir/transcriptome.fasta"
    log:
        'logs/trinity.log'
    params:
        extra=""
    threads: 6
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_gb=50
    wrapper:
        "0.72.0/bio/trinity"

rule transRate:
    input:
        assembly = "results/trinity_out_dir/transcriptome.fasta",
        left = "resources/reads/volatiles_1_1.fq.gz",
        right = "resources/reads/volatiles_2_2.fq.gz"
    output:
        directory("results/transrate")
    log:
        "logs/transrate.log"
    conda:
        "../envs/rnaseq.yaml"
    threads: 16
    shell:
        """
        transrate --assembly {input.assembly} \
            --left {input.left} --right {input.right} \
            --output {output} 
        """

rule transdecoderLongORFs:
    input:
        fasta = "results/trinity_out_dir/transcriptome.fasta",
        #gene_trans_map="test.gtm" # optional gene-to-transcript identifier mapping file (tab-delimited, gene_id<tab>trans_id<return> )
    output:
        "results/transdecoder_dir/longest_orfs.pep"
    log:
        "logs/transdecoder/longorfs.log"
    params:
        extra=""
    wrapper:
        "0.72.0/bio/transdecoder/longorfs"

rule transdecoderPredict:
    input:
        fasta = "results/trinity_out_dir/transcriptome.fasta",
        longorfs = "results/transdecoder_dir/longest_orfs.pep"
        #pfam_hits = "pfam_hits.txt", # optionally retain ORFs with hits by inputting pfam results here (run separately)
        #blastp_hits = "blastp_hits.txt", # optionally retain ORFs with hits by inputting blastp results here (run separately)
    output:
        "results/transdecoder/transdecoder.bed",
        "results/transdecoder/transdecoder.cds",
        "results/transdecoder/transdecoder.pep",
        "results/transdecoder/transdecoder.gff3"
    log:
        "logs/transdecoder/predict.log"
    params:
        extra=""
    wrapper:
        "0.72.0/bio/transdecoder/predict"

rule Busco:
    input:
        "results/trinity_out_dir/transcriptome.fasta"
    output:
        "results/busco/transcriptome.busco.tsv",
    log:
        "logs/busco.log"
    threads: 8
    params:
        mode="transcriptome",
        lineage_path="diptera_odb10",
        # optional parameters
        extra=""
    wrapper:
        "0.72.0/bio/busco"


rule KallistoIndex:
	input:
		fasta = "results/trinity_out_dir/transcriptome.fasta"
	output:
		index = "resources/reference/kallisto.idx"
	log:
		"logs/kallisto/index.log"
	wrapper:
		"0.72.0/bio/kallisto/index"


rule KallistoQuant:
	input:
		fastq = expand("resources/reads/{{sample}}_{n}.fq.gz", n=[1,2]),
		index = "resources/reference/kallisto.idx"
	output:
		directory("results/quant/{sample}")
	log:
		"logs/kallisto/quant_{sample}.log"
	params:
		extra = "-b 100"
	threads:24
	wrapper:
		"0.72.0/bio/kallisto/quant"

#rule DifferentialGeneExpression:
#	input:
#		samples = config['samples'],
#		gene_names = config['ref']['genenames'],
#		DEcontrasts = "resources/DE.contrast.list",
#		counts = expand("results/quant/{sample}", sample=samples)
#	output:
#		csvs = expand("results/genediff/{comp}.csv", comp=config['contrasts']),
#		xlsx = "results/genediff/RNA-Seq_diff.xlsx",
#		pca = "results/plots/PCA.pdf",
#		countstats = "results/quant/count_statistics.tsv"
#	priority: 10
#	conda:
#		"../envs/diffexp.yaml"
#	log:
#		"logs/DifferentialGeneExpression.log"
#	script:
#		"../scripts/DeseqGeneDE.R"