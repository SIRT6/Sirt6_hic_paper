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

anno_atac_dif
anno_atac

res_p <- list()
res_no_p <- list()
for (i in 1:(length(anno_atac_dif)-1)) {
  x <- as.data.frame(anno_atac_dif[[i]]@anno)
  x_p <- x[grepl("Promoter", x$annotation),]
  res_p[[i]] <- x_p
  x_no_p <- x[!grepl("Promoter", x$annotation),]
  res_no_p[[i]] <- x_no_p
}
names(res_p) <- names(anno_atac_dif)[1:5]
names(res_no_p) <- names(anno_atac_dif)[1:5]
save(res_p, file = "anno_atac_dif_promoter.RData")    #resulting filtered files (dif)
save(res_no_p, file = "anno_atac_dif_no_p.RData")


res_p <- list()
res_no_p <- list()
for (i in 1:length(anno_atac)) {
  x <- as.data.frame(anno_atac[[i]]@anno)
  x_p <- x[grepl("Promoter", x$annotation),]
  res_p[[i]] <- x_p
  x_no_p <- x[!grepl("Promoter", x$annotation),]
  res_no_p[[i]] <- x_no_p
}
names(res_p) <- names(anno_atac)
names(res_no_p) <- names(anno_atac)
save(res_p, file = "anno_atac_promoter.RData")     #resulting filtered files (narrow)
save(res_no_p, file = "anno_atac_no_p.RData")
