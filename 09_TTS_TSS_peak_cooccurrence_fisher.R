library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tidyverse)

setwd("/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/analysis/TTS_TSS_peak_cooccurrance")
gb_anno_table<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/raw_data/N4_filtered-annotated-peaks_minfoldch5_genebody_plus_500bpStart-to-plus-600bpStop_111924.tsv",
header=T,sep="\t",comment.char="",quote="\"",stringsAsFactors=F)

pro_anno_table<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/raw_data/N4_filtered-annotated-peaks_minfoldch5_minus-2000bp-to-plus-500bp_111623.tsv",
header=T,sep="\t",comment.char="",quote="\"",stringsAsFactors=F)


genes_TTS <- gb_anno_table %>%
    filter(species == "Arabidopsis_thaliana_Col-0",
            n_cons_species_minfrac0 == 4) %>%
    select(tf,gene) %>%
    unique() %>%
    group_by(tf) %>%
    summarise(genes = list(gene)) %>%
    deframe()  # Convert to named list

genes_TSS <- pro_anno_table %>%
    filter(species == "Arabidopsis_thaliana_Col-0",
            n_cons_species_minfrac0 == 4) %>%
    select(tf,gene) %>%
    unique() %>%
    group_by(tf) %>%
    summarise(genes = list(gene)) %>%
    deframe()  # Convert to named list

TF_list <- intersect(names(genes_TSS),names(genes_TTS))

# Initialize results data frame
results <- data.frame(TF = character(),
                      TF_TSS_target_num=numeric(),
                      TF_TTS_target_num=numeric(),
                      cooccur_target_num=numeric(),
                      Fisher_pvalue = numeric(),
                      P_TTS_given_TSS = numeric(),
                      P_TSS_given_TTS = numeric(),
                      Dependency = character(),
                      Logistic_pvalue = numeric(),
                      stringsAsFactors = FALSE)

# Iterate over each TF
for (TF in TF_list) {
  # Get gene lists for this TF
  tts_genes <- genes_TTS[[TF]]
  tss_genes <- genes_TSS[[TF]]
  
  # Create a union of all genes involved
  all_genes <- unique(c(tts_genes, tss_genes))

  # Construct the contingency table
  sample_size <- length(tts_genes)
  pop_size <- 27000
  pop_success <- length(tss_genes)
  sample_success <- length(intersect(tss_genes,tts_genes))

  ## probability of getting more than sample_success intersects between TTS and TSS targets
  phyper_res <- phyper(sample_success-1, pop_success, 27000-pop_success,sample_size,lower.tail=FALSE)



  # Save results
  results <- rbind(results, data.frame(
    TF = TF,
    TF_TSS_target_num=length(tss_genes),
    TF_TTS_target_num=length(tts_genes),
    cooccur_target_num=length(intersect(tss_genes,tts_genes)),
    phyper_pvalue = phyper_res,
    stringsAsFactors = FALSE
  ))
}

# Save results to file
write.table(results, "TF_dependency_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Print first few results
print(head(results))
