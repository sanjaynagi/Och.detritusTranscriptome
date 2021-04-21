### rule file for de novo transcriptome assembly and annotation

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


rule transdecoderLongORFs:
    input:
        fasta = "results/trinity_out_dir/Trinity.fasta",
        gene_trans_map = "results/trinity_out_dir/Trinity.fasta.gene_trans_map"
    output:
        "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep"
    log:
        "logs/transdecoder/longorfs.log"
    params:
        extra=""
    wrapper:
        "0.72.0/bio/transdecoder/longorfs"

rule transdecoderPredict:
    input:
        fasta = "results/trinity_out_dir/Trinity.fasta",
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep"
        #pfam_hits = "pfam_hits.txt", # optionally retain ORFs with hits by inputting pfam results here (run separately)
        #blastp_hits = "blastp_hits.txt", # optionally retain ORFs with hits by inputting blastp results here (run separately)
    output:
        "Ae.detritus.transdecoder.bed",
        "Ae.detritus.transdecoder.cds",
        "Ae.detritus.transdecoder.pep",
        "Ae.detritus.transdecoder.gff3"
    log:
        "logs/transdecoder/predict.log"
    params:
        extra=""
    wrapper:
        "0.72.0/bio/transdecoder/predict"


"""
The transdecoder wrapper above will throw output files into the current working directory
(the snakemake workflow base directory) and so it may be required to move/rename some files into results
at this point 
"""


rule mv1:
    input:
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        bed = "Ae.detritus.transdecoder.bed"
    output:
        "results/transdecoder/Ae.detritus.transdecoder_dir/longest_orfs.pep"
    log:
        "logs/mv.log"
    shell:
        """
        mv Ae.detritus.transdecoder_dir results/transdecoder 2> {log}
        mv Ae.detritus.trans* results/transdecoder/ 2> {log}
        """

rule Busco_prot:
    input:
        longorfs = "results/transdecoder/Ae.detritus.transdecoder_dir/longest_orfs.pep", 
    output:
        directory("results/busco/Ae.detritus_prot"),
    log:
        "logs/busco_prot.log"
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

rule Busco_nucl:
    input:
        fasta = "results/trinity_out_dir/Trinity.fasta"
    output:
        directory("results/busco/Ae.detritus_nucl"),
    log:
        "logs/busco_nucl.log"
    threads: 8
    params:
        mode="transcriptome",
        lineage="diptera_odb10",
    conda:
        "../envs/busco.yaml"
    shell:
        """
        busco -i {input.fasta} -l {params.lineage} -o Ae.detritus_nucl --out_path results/busco -m {params.mode} --cpu {threads}
        """


rule eggnog_dl:
    output:
        db = directory("resources/eggnog")
    log:
        "logs/eggnog_dl.log"
    shell:
        """
        download_eggnog_data.py --data_dir resources/eggnog/
        """


rule eggnog_nucl:
    input:
        fasta = "results/trinity_out_dir/Trinity.fasta",
        db = "resources/eggnog"
    output:
        annot = "results/eggnog_aedes_nucl.emapper.annotations"
    log:
        "logs/eggnog_nucl.log"
    threads: 12
    shell:
        """
        emapper.py -i {input.fasta} -m diamond -o results/eggnog_aedes_nucl --data_dir {input.db} --cpu 24 --translate
        """


rule eggnog_prot:
    input:
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep",
        db = "resources/eggnog"
    output:
        annot = "results/eggnog_aedes_prot.emapper.annotations"
    log:
        "logs/eggnog_prot.log"
    threads: 12
    shell:
        """
        emapper.py -i {input.longorfs} -m diamond -o results/eggnog_aedes_prot --data_dir {input.db} --cpu 24
        """

rule hmmScan:
    input:
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
        pfam = "resources/Pfam-A.hmm"
    output:
        pfam = "results/Ae.detritusPFAM.out"
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
    output:
        tmhmm = "results/tmhmm.out"
    threads: 12
    log:
        "logs/tmhmm.log"
    shell:
        """
        ./tmhmm-2.0c/bin/tmhmm --short < {input.longorfs} > {output.tmhmm} 2> {log}
        """


rule signalP:
    input:
        longorfs = "results/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep", 
    output:
        signalP = "results/signalp.out"
    log:
        "logs/signalp.log"
    shell:
        """
        cd signalp-5.0b/bin/ && ./signalp -fasta {workflow.basedir}/../{input.longorfs} -format short -prefix {workflow.basedir}/../{output.signalP}
        """



#### OLD stuff - trinotate and manual diamond searches ####

# rule trinotate:
#     input:
#         longorfs = "results/transdecoder/Ae.detritus.transdecoder_dir/longest_orfs.pep", 
#         pfam = "results/Ae.detritusPFAM.out",
#         fasta = "results/Ae.detritus.clust.transcriptome.fa",
#         genemap = "results/trinity_out_dir/Trinity.fasta.gene_trans_map",
#         blastp =  "results/Ae.det_blastp.outfmt6",
#         blastx =  "results/Ae.det_blastx.outfmt6",
#         tmhmm = "results/tmhmm.out",
#         signalp = "results/signalp.out"
#     output:
#         report = "results/trinotate_annotation_report.xls",
#     conda:
#         "../envs/trinotate.yaml"
#     threads: 12
#     log:
#         "logs/trinotate.log"
#     shell:
#         """
#         workflow/scripts/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite init \
#         --gene_trans_map {input.genemap} \
#         --transcript_fasta {input.fasta} \
#         --transdecoder_pep {input.longorfs} 2> {log}

#         workflow/scripts/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp {input.blastp} 2>> {log}
#         workflow/scripts/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx {input.blastx} 2>> {log}
#         workflow/scripts/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_pfam {input.pfam} 2>> {log}
#         workflow/scripts/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_tmhmm {input.tmhmm} 2>> {log}
#         workflow/scripts/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_signalp {input.signalp} 2>> {log}

#         workflow/scripts/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite report > {output.report} 2>> {log}
#         """


# rule trinotateGO:
#     input:
#         report = "results/trinotate_annotation_report.xls",
#     output:
#         go = "results/go_annotations.txt"
#     conda:
#         "../envs/trinotate.yaml"
#     threads: 12
#     log:
#         "logs/GOtrinotate.log"
#     shell:
#         """
#         workflow/scripts/Trinotate-Trinotate-v3.2.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls {input.report} \
#         -G --include_ancestral_terms > {output.go} 2> {log}
#         """



# #        Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate


# rule makeBlastDB:
#     input:
#         "resources/uniprot_sprot.pep"
#     output:
#         "resources/uniprot_sprot.db.dmnd"
#     log:
#         "logs/blastdb.log"
#     shell:
#         "diamond makedb --in {input} -d resources/uniprot_sprot.db --threads 24"

# rule DiamondNucl:
#     input:
#         fasta = "results/Ae.detritus.clust.transcriptome.fa",
#         db = "resources/uniprot_sprot.db.dmnd"
#     output:
#         blast = "results/Ae.det_blastx.outfmt6"
#     log:
#         "logs/diamondnucl.log"
#     threads: 12
#     shell:
#         """
#         diamond blastx -d {input.db} \
#         -q {input.fasta} --outfmt 6 --threads {threads} --out {output} \
#         -b12 -c1 --max-target-seqs 1 
#         """

# rule DiamondProt:
#     input:
#         longorfs = "results/transdecoder/Ae.detritus.transdecoder_dir/longest_orfs.pep", 
#         db = "resources/uniprot_sprot.db.dmnd"
#     output:
#         blast = "results/Ae.det_blastp.outfmt6"
#     log:
#         "logs/blastpro.log"
#     threads: 12
#     shell:
#         """
#         diamond blastp -d {input.db} \
#         -q {input.longorfs} --outfmt 6 --threads {threads} --out {output} -b12 -c1 --max-target-seqs 1 
#         """
