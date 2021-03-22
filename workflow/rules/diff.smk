
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