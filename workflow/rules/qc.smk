

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
