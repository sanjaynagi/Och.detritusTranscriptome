library(data.table)
library(tidyverse)
library(glue)

# Read Eggnog annotation data
ae_det_annot = fread("results/eggnog_aedes_nucl.emapper.annotations", 
                     sep="\t", skip=3, fill=TRUE)

# Rename and select 
annot_small = ae_det_annot %>%
  rename("GeneID" = `#query_name`, 
         "GeneDescription" = `eggNOG free text desc.`) %>% 
  select(GeneID, GeneDescription)

# Aggregate to gene level
annot_small$GeneID = gsub('(.*)_\\w+', '\\1', annot_small$GeneID)
desc_counts = annot_small %>% group_by(GeneID, GeneDescription) %>% tally()

# Remove comments 
desc_counts = desc_counts[4:nrow(desc_counts),]

# Dummy df without gene desc for duplicated() 
dddf = desc_counts[,c(1,3)]

# Get the gene descriptions with highest counts (first genes)
# If two descs with equal N, take first desc at random
desc_counts = desc_counts[!duplicated(dddf$GeneID),]

# Remove counts column 
gene_descs = desc_counts[,1:2]  
  
# Write to file
fwrite(gene_descs, "results/eggnog_gene_annots.tsv", sep="\t")
