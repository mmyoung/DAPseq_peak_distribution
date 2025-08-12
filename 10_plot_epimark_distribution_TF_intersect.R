library(GenomicFeatures)
library(ChIPseeker)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/analysis/peak_with_epimark_distribution_intersect")

pro_anno_table<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/raw_data/N4_filtered-annotated-peaks_minfoldch5_minus-2000bp-to-plus-500bp_111623.tsv",
header=T,sep="\t",comment.char="",quote="\"",stringsAsFactors=F)

read.csv('/clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/raw_data/ath-258-tf-info_simple.csv', header=T) -> tf_info

tf_list <- pro_anno_table %>% filter(species == "Arabidopsis_thaliana_Col-0") %>% pull(tf) %>% unique()

ATG_gr <- read.table("/clusterfs/jgi/groups/gentech/seqtech/plant_multidap_data/genomes/annotations/Arabidopsis_thaliana_Col-0_cds_primary.gff",header = F,stringsAsFactors = F) %>%
makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V4",end.field = "V5",strand.field = "V7")

k27ac_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K27ac_leaf.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
    makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

# k27me3_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K27me3_Col0.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
#     makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

# k36me3_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K36me3_WT.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
#     makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

# k4me2_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K4me2_WT_12.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
#     makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

# k4me3_gr<- read.table("/clusterfs/jgi/groups/gentech/homes/romalley/reference/epimarks_bed/H3K4me3_WT_12.peak.idr.hammock.gz",header = F,stringsAsFactors = F) %>%
#     makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")

for(TF in tf_list){

    TF_c4_peak_gr <- pro_anno_table %>%
        filter(tf==!!TF,
            species == "Arabidopsis_thaliana_Col-0",
            n_cons_species_minfrac0 == 4) %>%
        makeGRangesFromDataFrame(seqnames.field = "peak_chr",start.field = "peak_start",end.field = "peak_end")

    # TF_c1_peak_gr <- pro_anno_table %>%
    #     filter(tf==!!TF,
    #         species == "Arabidopsis_thaliana_Col-0",
    #         n_cons_species_minfrac0 == 1) %>%
    #     makeGRangesFromDataFrame(seqnames.field = "peak_chr",start.field = "peak_start",end.field = "peak_end")

    k27ac_tf<- subsetByOverlaps(x = TF_c4_peak_gr, ranges = k27ac_gr)

    hits <- findOverlaps(TF_c4_peak_gr, k27ac_gr)
    k27ac_tf <- TF_c4_peak_gr[queryHits(hits)]
    tf_k27ac <- k27ac_gr[subjectHits(hits)]

    k27ac_only <- k27ac_gr[-subjectHits(hits)]
    tf_only <- TF_c4_peak_gr[-queryHits(hits)]

    if(length(k27ac_tf)>100 & length(k27ac_only)>100 & length(tf_only)>100){
        tagMatrixList <- lapply(list(K27ac=k27ac_only,tf_only=tf_only,tf_bind_in_k27ac=k27ac_tf), 
                            getTagMatrix, 
                            windows=makeBioRegionFromGranges(ATG_gr,type = "start_site",upstream = 3000, downstream = 600,by="gene"))

        TF_family <- tf_info %>% filter(gene_id==TF) %>% pull(tf_family)
        TF_name <- tf_info %>% filter(gene_id==TF) %>% pull(tf_name)
    
        plotAvgProf(tagMatrixList, xlim=c(-3000, 600),
        xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency",origin_label = "ATG") +
        ggtitle(paste0(TF_name,"(",TF_family,")\n","peak in k27ac#:",length(k27ac_tf),"tf_only#:",length(tf_only))) ->plot

        filename <- file.path(paste0(TF_family,"_",TF, "_k27ac_ATG_avgplot.pdf"))
        pdf(width = 7, height = 6, filename)
        print(plot)
        dev.off()
        
    }
 

}



