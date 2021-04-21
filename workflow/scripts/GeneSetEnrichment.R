#!/usr/bin/env Rscript
#' This script runs GSEA analysis on each contrast in the workflow
#' on gene lists, ranked using fold-change from DE data 

# GSEA RNA-Seq Ag
library(GO.db)
library(gage)
library(fgsea)
library(data.table)
library(glue)
library(tidyverse)


## functions ##
runEnrich = function(rankedList, GeneSetList, outName){
  
  for (set in c("GO")){
    # pathway enrichment with fgsea for both kegg and  GO terms
    fgseaRes = fgsea(pathways = GeneSetList, 
                     stats = na.omit(rankedList),
                     minSize  = 15,
                     maxSize  = 1000,
                     nperm=10000)
    
    fgseaRes %>% arrange(padj) %>% fwrite(., file = glue("results/gsea/{outName}.{set}.tsv"), sep = "\t")
    
    topPathwaysUp = fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown = fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways = c(topPathwaysUp, rev(topPathwaysDown))
    pdf(glue("results/gsea/{outName}.{set}.fgsea.pdf"))
    plotGseaTable(GeneSetList[topPathways], na.omit(rankedList), fgseaRes)
    null = dev.off()
    print(glue("{set} complete"))
  }
}

###### configuration - metadata and parameters ######
samples = fread("config/samples.tsv") %>% as.data.frame()
contrasts = c("Control_EPNexposed", "Control_volatiles", "EPNexposed_volatiles")
 
GeneSetList = list()
######### GO Terms #########
go = fread("results/eggnog_aedes_nucl.emapper.annotations", sep="\t", skip=3, fill=TRUE) 
go = go %>% as.data.frame() %>% 
  distinct() %>% 
  dplyr::rename("GeneID" = `#query_name`, "GO.id" = GOs) %>% 
  dplyr::select(GeneID, GO.id)

# Aggregate to gene level
go$GeneID = gsub('(.*)_\\w+', '\\1', go$GeneID)
go_counts = go %>% group_by(GeneID, GO.id) %>% tally()

# Remove comments 
go_counts = go_counts[4:nrow(go_counts),]
# Dummy df without gene desc for duplicated() 
dddf = go_counts[,c(1,3)]
# Get the gene descriptions with highest counts (first genes)
# If two descs with equal N, take first desc at random
go_counts = go_counts[!duplicated(dddf$GeneID),][,1:2]

#separate into distint rows 
go = go_counts %>% separate_rows(GO.id, sep=",")

#get descriptions
goterms = Term(GOTERM) %>% enframe() %>% dplyr::rename("GO.id" = "name", "description" = "value")


# join with our go terms
go = inner_join(go, goterms)
go = go %>% dplyr::select(description, GeneID)
# turn into a large list of terms and all the genes with that term
GOlist = go %>% pivot_wider(names_from="description", values_from="GeneID") %>% purrr::transpose()


##### gene DE GSEA ####
for (comp in contrasts[1]){
  print(glue("Running GO enrichment analyses for {comp}"))
  # make ranked list using DE results, rank based on log2foldchange
  rank = fread(glue("results/genediff/{comp}.csv")) %>% 
    arrange("log2FoldChange") %>% 
    dplyr::select(c('GeneID', all_of("log2FoldChange"))) %>% 
    deframe()
  runEnrich(rankedList = sort(rank, decreasing=TRUE), GeneSetList = GOlist[[1]], outName = glue("genediff/{comp}.DE"))
}




sessionInfo()