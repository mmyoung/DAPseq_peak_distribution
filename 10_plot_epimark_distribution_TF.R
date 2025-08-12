library(GenomicFeatures)
library(ChIPseeker)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/analysis/peak_with_epimark_distribution")

pro_anno_table<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/raw_data/N4_filtered-annotated-peaks_minfoldch5_minus-2000bp-to-plus-500bp_111623.tsv",
header=T,sep="\t",comment.char="",quote="\"",stringsAsFactors=F)

read.csv('/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/raw_data/ath-258-tf-info_simple.csv', header=T) -> tf_info

tf_list <- pro_anno_table %>% filter(species == "Arabidopsis_thaliana_Col-0") %>% pull(tf) %>% unique()

ATG_gr <- read.table("/clusterfs/jgi/groups/gentech/seqtech/plant_multidap_data/genomes/annotations/Arabidopsis_thaliana_Col-0_cds_primary.gff",header = F,stringsAsFactors = F) %>%
makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V4",end.field = "V5",strand.field = "V7")

k27ac_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K27ac_leaf.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
    makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

k27me3_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K27me3_Col0.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
    makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

k36me3_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K36me3_WT.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
    makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

k4me2_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K4me2_WT_12.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
    makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

k4me3_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K4me3_WT_12.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
    makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

for(TF in tf_list){

    TF_c4_peak_gr <- pro_anno_table %>%
        filter(tf==!!TF,
            species == "Arabidopsis_thaliana_Col-0",
            n_cons_species_minfrac0 == 4) %>%
        makeGRangesFromDataFrame(seqnames.field = "peak_chr",start.field = "peak_start",end.field = "peak_end")

    TF_c1_peak_gr <- pro_anno_table %>%
        filter(tf==!!TF,
            species == "Arabidopsis_thaliana_Col-0",
            n_cons_species_minfrac0 == 1) %>%
        makeGRangesFromDataFrame(seqnames.field = "peak_chr",start.field = "peak_start",end.field = "peak_end")

    if(length(TF_c4_peak_gr)>100 & length(TF_c1_peak_gr)>100){
        tagMatrixList <- lapply(list(K27ac=k27ac_gr,K27me3=k27me3_gr,K36me3=k36me3_gr,k4me2=k4me2_gr,k4me3=k4me3_gr,c4=TF_c4_peak_gr,c1=TF_c1_peak_gr), 
                            getTagMatrix, 
                            windows=makeBioRegionFromGranges(ATG_gr,type = "start_site",upstream = 3000, downstream = 600,by="gene"))

        TF_family <- tf_info %>% filter(gene_id==TF) %>% pull(tf_family)
        TF_name <- tf_info %>% filter(gene_id==TF) %>% pull(tf_name)
    
        plotAvgProf(tagMatrixList, xlim=c(-3000, 600),
        xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency",origin_label = "ATG") +
        ggtitle(paste0(TF_name,"(",TF_family,")\n","c1#:",length(TF_c1_peak_gr),";","c4#:",length(TF_c4_peak_gr))) ->plot

        filename <- file.path(paste0(TF_family,"_",TF, "&epimark_ATG_avgplot.pdf"))
        pdf(width = 7, height = 6, filename)
        print(plot)
        dev.off()
        
    }
 

}



