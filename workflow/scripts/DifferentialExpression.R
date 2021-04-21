 #!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#' This script performs differential expression analysis at the 
#' gene level using DESeq2. It outputs heatmaps, PCAs, gene counts,
#' as well as DE results. Final DE results for all contrasts are 
#' saved in an .xslx report 

## Differential expression
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(data.table)
library(ggrepel)
library(openxlsx)
library(glue)
library(RColorBrewer)
library(EnhancedVolcano)

#read metadata and get contrasts
samples = fread("config/samples.tsv", sep="\t") %>% as.data.frame()
contrasts = c("Control_EPNexposed", "Control_volatiles", "EPNexposed_volatiles")
gene_names = fread("results/eggnog_gene_annots.tsv", sep="\t")

# remove sample 4 EPN as its an outlier on PCAs
samples = samples %>% filter(Name != "EPN_exposed_4")


##### define functions ######
round_df = function(df, digits) {
  
  #' This function rounds all the numeric columns of a data.frame
  nums = vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] = round(df[,nums], digits = digits)
  (df)
}

vst_pca = function(counts, samples, colourvar, name="PCA", st="", comparison=""){
  
  #' This function takes counts and sample metadata as input and builds a DESeq2 dataset
  #' Returning normalised and variance-stabilised counts, and performs PCA on the data
  #' Plotting and saving to pdf

  # make DESeq dataset
  dds = DESeqDataSetFromMatrix(countData = counts, 
                               colData = samples, 
                               design = ~ treatment)
  ###### estimate paramters and normalise 
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  vsd = varianceStabilizingTransformation(dds)
  normcounts = counts(dds, normalized=TRUE)
  vstcounts = assay(vsd)
  vstcounts = vstcounts[order(apply(vstcounts,1,sum),decreasing=TRUE),]
  
  #### write pca of samples to pdf
  pca2=prcomp(t(vstcounts),center=TRUE)
  pc = data.frame(pca2$x) %>% rownames_to_column("Name")
  pc = left_join(pc, samples)

  pdf(glue("results/plots/{name}{st}{comparison}.pdf"))  
  print(ggplot(data=pc,aes(x=PC1, y=PC2, colour=treatment)) + 
    geom_point(size=6, alpha=0.8) + 
    geom_text_repel(aes(label=Name), colour="black") + 
    theme_light() + 
    labs(title=glue("{name} {st} {comparison}"),
         x=glue("PC1 - Variance explained - {round(summary(pca2)$importance[2,][1], 3)}"),
         y=glue("PC2 - Variance explained - {round(summary(pca2)$importance[2,][2], 3)}"))  + 
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background=element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour="black", size=18),
          axis.text.y = element_text(colour="black", size=18)))
  null = dev.off()
  
  return(list(vstcounts, dds, normcounts))
}


#### main ####
cat("\n", "------------- Kallisto - DESeq2 - RNASeq Differential expression ---------", "\n")

#### Read counts for each sample
df = list()
for (sample in samples$Name){
  df[[sample]]= fread(glue("results/quant/{sample}/abundance.tsv"), sep = "\t")
}

counts = data.frame('GeneID' = df[[1]]$target_id)
#get read counts for each gene and fill table
for (sample in samples$Name){
  reads = df[[sample]]$est_counts
  counts = cbind(counts, reads)
}

#rename columns
colnames(counts) = c("GeneID", samples$Name)

## aggregate to gene level
counts$GeneID = gsub('(.*)_\\w+', '\\1', counts$GeneID)
counts = counts %>% group_by(GeneID) %>% summarise_all(sum)

### count total reads counted
counts = counts %>% column_to_rownames('GeneID')
#### get count statistics for each sample and plot ####
count_stats = apply(counts, 2, sum) %>% enframe(name="Sample", value="total_counts") # total counts
ngenes = nrow(counts)
count_stats$genes_zerocounts = apply(counts, 2, function(x){sum(x==0)}) # genes with zero counts 
count_stats$genes_lessthan10counts = apply(counts, 2, function(x){sum(x<10)}) # genes with less than 10 counts
count_stats = count_stats %>% mutate("proportion_zero" = genes_zerocounts/ngenes,
                                     "proportion_low" = genes_lessthan10counts/ngenes)
count_stats %>% fwrite(., "results/quant/count_statistics.tsv",sep="\t")

print("Counting and plotting total reads per sample...")
pdf("results/quant/total_reads_counted.pdf")
ggplot(count_stats, aes(x=Sample, y=total_counts, fill=samples$treatment)) + 
  geom_bar(stat='identity') + 
  theme_light() +
  ggtitle("Total reads counted (mapped to Trinity Ae.detritus transcriptome assembly)") +
  theme(axis.text.x = element_text(angle=45),
        plot.title = element_text(hjust = 0.5))
null = dev.off() 
 
# round numbers to be whole, they are not due averaging across transcripts
counts = counts %>% rownames_to_column('GeneID') %>% 
  mutate_if(is.numeric, round) %>% column_to_rownames('GeneID')

# remove lowly expressed genes
counts = counts[rowMeans(counts) > 5,]


############ Plots PCA with all data, and performs DESeq2 normalisation ########################
res = vst_pca(counts, samples, colourvar = 'treatment', name="PCA")
vstcounts = res[[1]]
dds = res[[2]]
normcounts = res[[3]]

### write out raw and normalised counts 
counts %>% 
  rownames_to_column("GeneID") %>% 
  fwrite(., "results/quant/rawcounts.tsv", sep="\t", row.names = FALSE)

normcounts %>% 
  as.data.frame() %>% 
  rownames_to_column("GeneID") %>% 
  round_df(., 1) %>% 
  fwrite(., "results/quant/normcounts.tsv", sep="\t", row.names = FALSE)

# calculate correlations between samples based on the count data, and plot heatmap
correlations = cor(vstcounts)
pdf("results/plots/heatmap_correlations.pdf")
pheatmap(correlations)
garbage = dev.off()


results_list = list()
names_list = list()
######### subset data and run DESeq for each combo we need, store in xlsx ########
for (cont in contrasts){
  control = str_split(cont, "_")[[1]][1] # get first of string, which is control 
  case = str_split(cont, "_")[[1]][2] # get case 
  controls = which(samples$treatment %in% control)
  cases = which(samples$treatment %in% case)
  ## Perform PCA for each comparison
  subcounts = counts[,c(controls, cases)]
  subsamples = samples[c(controls, cases),]
  # make treatment a factor with the 'susceptible' as reference
  subsamples$treatment = as.factor(as.character(subsamples$treatment))
  subsamples$treatment = relevel(subsamples$treatment, control)
  null = vst_pca(subcounts, subsamples, colourvar='treatment', "PCA_", comparison=cont)[[]]
  
  # perform DE with each comparison, using DESeq dataset (dds) from whole library
  cat("\n", glue("--- Running DESeq2 differential expression analysis on {cont} ---"), "\n")
  cds = nbinomWaldTest(dds)
  results = results(cds, contrast = c("treatment", case, control)) %>% as.data.frame() 
  results = results[order(results$padj),] #order by pvalue 
  results = results %>% rownames_to_column("GeneID") %>% mutate("FC" = (2^log2FoldChange))
    
  ### absolute difference
  #### Get rowsums of counts, grouping by case/control. Then get difference of counts and join with DE results
  readdiff = data.frame(t(rowsum(t(subcounts), group = subsamples$treatment, na.rm = T))) #transpose and get rowsums for each group
  readdiff$absolute_diff = readdiff[,case] - readdiff[,control] #get difference
  readdiff = data.frame(readdiff) %>% rownames_to_column('GeneID')
  results = unique(left_join(results, readdiff[,c('GeneID','absolute_diff')]))
  
  # join DE results with normal gene names
  results = unique(left_join(results, gene_names))
  fwrite(results, glue("results/genediff/{cont}.csv")) #write to csv 
  # volcano plot for each comparison, using EnhancedVolcano. First make vector of labels which is AGAPs unless a gene name exists
  #labels = results %>% mutate("Gene_name" = case_when(Gene_name == "" ~ GeneID,
   #                                  Gene_name == NA ~ GeneID,
    #                                 TRUE ~ Gene_name)) %>% select(Gene_name) %>% deframe()

  # store names of comparisons for xlsx report sheets
  results_list[[cont]] = results
  names_list[[cont]] = cont
  pdf(glue("results/genediff/Volcano_plot_{cont}.pdf"))
  print(EnhancedVolcano(results_list[[cont]],
                        lab=results_list[[cont]]$GeneID,
                        x='log2FoldChange',
                        y='pvalue',
                        title = cont))
  null = dev.off()
  cat("\n", glue("{cont} complete!"), "\n")
}

#### write to excel file on diff sheets #### 
sheets = unlist(names_list)
wb <- createWorkbook("Workbook")

for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = FALSE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file="results/genediff/Ae.detritus_DE_noEPN4.xlsx", overwrite = TRUE)

sessionInfo()
