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

# load("anno_atac_dif_promoter.RData")
# load("/Users/annaponomareva/Documents/Research/SIRT6_enh/genes_universe3.RData")

difpeak_names <- c('KO_old_dif', 'KO_WT_dif', 'KO_young_dif', 'old_WT_dif', 'old_young_dif', 'young_WT_dif')

egos_p <- c()
go_plots_p <- c()

go_atac_data <- res_exon_dif[c(2,3,1,4,5)]
go_names <- difpeak_names[c(2,3,1,4,5)]
go_names_plot <- c("KO-WT", "KO-Young", "KO-Old", "Old-WT", "Old-Young")














