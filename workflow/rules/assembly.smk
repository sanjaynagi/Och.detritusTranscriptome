
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


rule cdhit:
    input:
        assembly = "results/trinity_out_dir/Trinity.fasta",
    output:
        "results/Ae.detritus.cdhit.transcriptome.fa"    
    log:
        "logs/cdhit.log"
    threads: 16
    shell:
        """
        cd-hit-est -i results/trinity_out_dir/Trinity.fasta \
        -o {output} -c 0.95 -n 8 -T {threads} -M 4000
        """

rule transRate:
    input:
        assembly = "results/trinity_out_dir/Trinity.fasta",
        left = "resources/reads/volatiles_1_1.fq.gz",
        right = "resources/reads/volatiles_2_2.fq.gz"
    output:
        directory("results/transrate")
    log:
        "logs/transrate.log"
    threads: 16
    shell:
        """
        workflow/scripts//transrate-1.0.3-linux-x86_64/transrate --assembly {input.assembly} \
            --left {input.left} --right {input.right} \
            --output {output} 
        """

rule transdecoderLongORFs:
    input:
        fasta = "results/trinity_out_dir/Trinity.fasta",
        gene_trans_map="results/trinity_out_dir/Trinity.fasta.gene_trans_map"
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
        #fasta = "results/trinity_out_dir/Trinity.fasta"
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
    output:
        directory("results/busco/Ae.detritus_prot"),
    log:
        "logs/busco.log"
    threads: 8
    params:
        mode="prot",
        lineage="diptera_odb10",
    conda:
        "../envs/busco.yaml"
    shell:
        """
        busco -i {input.longorfs} -l {params.lineage} -o Ae.detritus_prot --out_path results/busco -m {params.mode} --cpu {threads}
        """

rule makeBlastDB:
    input:
        "resources/uniprot_sprot.fasta"
    output:
        "resources/uniprot.sprot.db.dmnd"
    log:
        "logs/blastdb.log"
    shell:
        "diamond makedb --in {input} -d resources/uniprot_sprot.db --threads 24"

rule DiamondNucl:
    input:
        fasta = "results/trinity_out_dir/Trinity.fasta",
        db = "resources/uniprot.sprot.db.dmnd"
    output:
        blast = "results/Ae.det_blastx.outfmt6"
    log:
        "logs/diamondnucl.log"
    threads: 12
    shell:
        """
        diamond blastx -d {input.db} \
        -q {input.fasta} --outfmt 6 --threads {threads} --out {output} \
        -b12 -c1 --max-target-seqs 1 
        """

rule DiamondProt:
    input:
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        db = "resources/uniprot.sprot.db.dmnd"
    output:
        blast = "results/Ae.det_blastp.outfmt6"
    log:
        "logs/blastpro.log"
    threads: 12
    shell:
        """
        diamond blastp -d {input.db} \
        -q {input.longorfs} --outfmt 6 --threads {threads} --out {output} -b12 -c1 --max-target-seqs 1 
        """


rule hmmScan:
    input:
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
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
        hmmsearch --cpu {threads} --domtblout {output.pfam} {input.pfam} {input.longorfs} 2> {log}
        """

rule tmhmm:
    input:
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        pfam = "resources/Pfam-A.hmm"
    output:
        tmhmm = "results/tmhmm.out"
    threads: 12
    log:
        "logs/tmhmm.log"
    shell:
        """
        ./tmhmm-2.0c/bin/tmhmm --short < {input.longorfs} > {output.tmhmm} 2> {log}
        """


rule trinotate:
    input:
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        pfam = "results/TrinotatePFAM.out",
        fasta = "results/trinity_out_dir/Trinity.fasta",
        genemap = "results/trinity_out_dir/Trinity.fasta.gene_trans_map",
        blastp =  "results/Ae.det_blastp.outfmt6",
        blastx =  "results/Ae.det_blastx.outfmt6",
        tmhmm = "results/tmhmm.out"
    output:
        report = "results/trinotate_annotation_report.xls",
        go = "results/go_annotations.txt"
    conda:
        "../envs/trinotate.yaml"
    threads: 12
    log:
        "logs/trinotate.log"
    shell:
        """
        Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate

        Trinotate Trinotate.sqlite init \
        --gene_trans_map {input.genemap} \
        --transcript_fasta {input.fasta} \
        --transdecoder_pep {input.longorfs} 2> {log}

        Trinotate Trinotate.sqlite LOAD_swissprot_blastp {input.blastp} 2>> {log}
        Trinotate Trinotate.sqlite LOAD_swissprot_blastx {input.blastx} 2>> {log}
        Trinotate Trinotate.sqlite LOAD_pfam {input.pfam} 2>> {log}
        Trinotate Trinotate.sqlite LOAD_tmhmm {input.tmhmm} 2>> {log}

        Trinotate Trinotate.sqlite report > {output.report} 2>> {log}

        extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls {output.report} \
        -G --include_ancestral_terms > {output.go} 2>> {log}
        """



