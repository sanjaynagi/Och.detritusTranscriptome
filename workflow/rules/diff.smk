
rule SalmonIndex:
	input:
		"results/Ae.detritus.cdhit.transcriptome.fa"
	output:
		directory("resources/reference/transcriptome.idx")
	log:
		"logs/salmon/index.log"
	wrapper:
		"0.72.0/bio/salmon/index"

rule Salmon:
	input:
		r1 = "resources/reads/{sample}_1.fq.gz",
		r2 = "resources/reads/{sample}_2.fq.gz",
		index = "resources/reference/transcriptome.idx"
	output:
		directory("results/quant/{sample}")
	log:
		"logs/salmon/quant_{sample}.log"
	params:
		extra = "--dumpEq"
	threads:24
	wrapper:
		"0.72.0/bio/salmon/quant"


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
