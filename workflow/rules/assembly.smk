
rule trinityAssembly:
    input:
        left=expand("resources/reads/{sample}_1.fq.gz", sample=samples),
        right=expand("resources/reads/{sample}_2.fq.gz", sample=samples)
    output:
        "results/trinity_out_dir/Trinity.fasta"
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
        assembly = "results/trinity_out_dir/Trinity.fasta",
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
        fasta = "results/trinity_out_dir/Trinity.fasta",
        #gene_trans_map="test.gtm" # optional gene-to-transcript identifier mapping file (tab-delimited, gene_id<tab>trans_id<return> )
    output:
        "Trinity.fasta.transdecoder_dir/longest_orfs.pep"
    log:
        "logs/transdecoder/longorfs.log"
    params:
        extra=""
    wrapper:
        "0.72.0/bio/transdecoder/longorfs"

rule transdecoderPredict:
    input:
        fasta = "results/trinity_out_dir/Trinity.fasta",
        longorfs = "Trinity.fasta.transdecoder_dir/longest_orfs.pep"
        #pfam_hits = "pfam_hits.txt", # optionally retain ORFs with hits by inputting pfam results here (run separately)
        #blastp_hits = "blastp_hits.txt", # optionally retain ORFs with hits by inputting blastp results here (run separately)
    output:
        "Trinity.fasta.transdecoder.bed",
        "Trinity.fasta.transdecoder.cds",
        "Trinity.fasta.transdecoder.pep",
        "Trinity.fasta.transdecoder.gff3"
    log:
        "logs/transdecoder/predict.log"
    params:
        extra=""
    wrapper:
        "0.72.0/bio/transdecoder/predict"

rule mv1:
    input:
        longorfs = "Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        bed = "Trinity.fasta.transdecoder.bed"
    output:
        "results/transdecoder/Trinity.fasta.transdecoder.bed"
    log:
        "logs/mv.log"
    shell:
        """
        mv Trinity.fasta.transdecoder_dir results/transdecoder 2> {log}
        mv Trinity.fasta.trans* results/transdecoder/ 2> {log}
        """

rule Busco:
    input:
        fasta = "results/trinity_out_dir/Trinity.fasta"
    output:
        "results/busco/transcriptome.busco.tsv",
    log:
        "logs/busco.log"
    threads: 8
    params:
        mode="transcriptome",
        lineage="diptera_odb10",
    conda:
        "../envs/busco.yaml"
    shell:
        """
        busco -in {input.fasta} --out results/busco --force --cpu {threads} --mode {params.mode} --linea$
        """

rule makeBlastDB:
    input:
        "resources/uniprot_sprot.fasta"
    output:
        touch(".db.complete")
    log:
        "logs/blastdb.log"
    conda:
        "../envs/rnaseq.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

rule BlastNucl:
    input:
        "results/Trinity_out_dir/Trinity.fasta",
        ".db.complete"
    output:
        touch(".blast.complete")
    log:
        "logs/blastnucl.log"
    conda:
        "../envs/rnaseq.yaml"
    threads: 12
    shell:
        """
        blastx -query Trinity.fasta -db resources/uniprot_sprot.fasta -out blastx.outfmt6 \
        -evalue 1e-20 -num_threads {threads} -max_target_seqs 1 -outfmt 6
        """

rule BlastProt:
    input:
        fasta = "results/Trinity_out_dir/Trinity.fasta",
        longorfs = "Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        db = ".db.complete"
    output:
        touch(".blastprot.complete")
    log:
        "logs/blastpro.log"
    conda:
        "../envs/rnaseq.yaml"
    threads: 12
    shell:
        """
        blastp -query {input.longestorfs} \
        -db resources/uniprot_sprot.fasta \
        -num_threads {threads} \
        -max_target_seqs 1 \
        -outfmt 6 > blastp.outfmt6
        """


rule hmmScan:
    input:
        longorfs = "Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        pfam = "resources/Pfam-A.hmm"
    output:
        pfam = "results/TrinotatePFAM.out"
    conda:
        "../envs/trinotate.yaml"
    threads: 12
    log:
        "logs/hmmscan.log"
    shell:
        """
        hmmscan --cpu {threads} --domtblout {output.pfam} {input.pfam} 2> {log}
        """

rule tmhmm:
    input:
        longorfs = "Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        pfam = "resources/Pfam-A.hmm"
    output:
        tmhmm = "results/tmhmm.out"
    conda:
        "../envs/trinotate.yaml"
    threads: 12
    log:
        "logs/tmhmm.log"
    shell:
        """
        tmhmm --short < {input.longorfs} > {output.tmhmm} 2> {log}
        """
