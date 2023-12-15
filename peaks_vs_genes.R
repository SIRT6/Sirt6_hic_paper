setwd('/Users/annaponomareva/Documents/Research/SIRT6_enh/promoters filter/')

library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)
library(soGGi)
library(ChIPseeker)
library(Rsubread)
library(GenomicRanges)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggplot2)

load('/Users/annaponomareva/Documents/Research/SIRT6_enh/dif no old 3/anno_atac_dif.RData')
load('/Users/annaponomareva/Documents/Research/SIRT6_enh/annotation ATAC/anno_atac.RData')

###---Filter peaks: exon + intron + promoter---###
res_exon <- list()
for (i in 1:length(anno_atac)) {
  x <- as.data.frame(anno_atac[[i]]@anno)
  x_p <- x[grepl("Exon", x$annotation) | grepl("Intron", x$annotation) | grepl("Promoter", x$annotation),]
  res_exon[[i]] <- x_p
}
names(res_exon) <- names(anno_atac)

res_exon_dif <- list()
for (i in 1:length(anno_atac_dif)) {
  x <- as.data.frame(anno_atac_dif[[i]])
  x_p <- x[grepl("Exon", x$annotation) | grepl("Intron", x$annotation) | grepl("Promoter", x$annotation),]
  res_exon_dif[[i]] <- x_p
}
names(res_exon_dif) <- names(anno_atac_dif)



###---GO annotation for filtered peaks---###
library(clusterProfiler)
library(Signac)
library(biomaRt)
library(org.Mm.eg.db)

egENSEMBL <- toTable(org.Mm.egENSEMBL)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

load("anno_atac_dif_promoter.RData")
load("/Users/annaponomareva/Documents/Research/SIRT6_enh/genes_universe3.RData")

difpeak_names <- c('KO_old_dif', 'KO_WT_dif', 'KO_young_dif', 'old_WT_dif', 'old_young_dif', 'young_WT_dif')

egos_p <- c()
go_plots_p <- c()

go_atac_data <- res_exon_dif[c(2,3,1,4,5)]
go_names <- difpeak_names[c(2,3,1,4,5)]
go_names_plot <- c("KO-WT", "KO-Young", "KO-Old", "Old-WT", "Old-Young")

for (i in 1:length(go_atac_data)) {
  sample_data <- makeGRangesFromDataFrame(go_atac_data[[i]], keep.extra.columns = T)
  sample_data <- sample_data[!sample_data %over% blacklist_mm10]
  sample_data <- as.data.frame(sample_data)
  sample_data$transcriptId <- gsub("\\..*","",sample_data$transcriptId)
  
  sample_genes <- getBM(attributes = c('ensembl_gene_id', "description",
                                       "chromosome_name", "start_position",
                                       "end_position", "ensembl_transcript_id"),
                        filters = "ensembl_transcript_id",
                        values = sample_data$transcriptId,
                        mart = mart)
  
  # ego_i <- enrichGO(gene = sample_genes$ensembl_gene_id,
  #                   universe = genes_universe3,
  #                   keyType = "ENSEMBL",
  #                   OrgDb = org.Mm.eg.db,
  #                   ont = "all")
  #                   # qvalueCutoff = 1,
  #                   # pvalueCutoff = 1)
                    
  ego_i <- enrichGO(gene = sample_genes$ensembl_gene_id,
                    universe = genes_universe3, 
                    keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db, 
                    ont = "BP")                                 #only Biological Processes (BP)
  ego_i <- clusterProfiler::simplify(ego_i)
  
  if (nrow(ego_i@result) != 0) {
    egos_p[[i]] <- ego_i
    # go_plot_i <- barplot(ego_i, split = "ONTOLOGY", showCategory = 5,
    #                      title = go_names[i],
    #                      font.size = 10,
    #                      # cols.use = c("lightgrey", "blue"),
    # ) +
    #   scale_fill_continuous(name = 'FDR p-value') +
    #   facet_grid(ONTOLOGY~., scale = "free")
    
    go_plot_i <- barplot(ego_i, showCategory = 5,            #only Biological Processes (BP)
                         title = go_names_plot[i],
                         font.size = 10) +
      scale_fill_continuous(name = 'FDR p-value')              
    
    go_plots_p[[i]] <- go_plot_i
    
  } else {
    egos_p[[i]] <- ego_i
    go_plots_p[[i]] <- 'no enriched terms'
  }
}

names(egos_p) <- go_names                               #promoter + exon + intron (only BP)
names(go_plots_p) <- go_names  
go_plots_p_BP <- go_plots_p
save(go_plots_p_BP, file = 'go_plots_p_BP.RData')       #promoter + exon + intron (only BP) - plots












