
rule SalmonIndex:
	input:
		"results/trinity_out_dir/Trinity.fasta"
	output:
		directory("resources/reference/transcriptome.idx")
	log:
		"logs/salmon/index.log"
	shell:
		"salmon index -t {input} -i {output} 2> {log}"

rule Salmon:
	input:
		r1 = "resources/reads/{sample}_1.fq.gz",
		r2 = "resources/reads/{sample}_2.fq.gz",
		index = "resources/reference/transcriptome.idx"
	output:
		quant = directory("results/quant/{sample}")
	log:
		"logs/salmon/quant_{sample}.log"
	threads:8
	shell:
		"salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -p {threads} --dumpEq -o {output.quant} 2> {log}"


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
