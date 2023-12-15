setwd('/Users/annaponomareva/Documents/Research/SIRT6_enh/')

######################------------Dif ATAC-seq wuthout old3------------######################
library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)
library(soGGi)


peaks <- dir(path = 'narrowPeak_no_old3', pattern = "*.narrowPeak", full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)

names(myPeaks) <- c("KO_R1", "KO_R2", "KO_R3",
                    "old_R1", "old_R2",
                    "WT_R1", "WT_R2", "WT_R3",
                    "young_R1", "young_R2", "young_R3")
Group <- factor(c("KO", "KO", "KO", 
                  'old', 'old', 
                  'WT', 'WT', 'WT',
                  'young', 'young', 'young'))

runConsensusRegions1 <- function(testRanges,method="majority", overlap="any"){
  if(length(testRanges) > 1){
    
    reduced <- reduce(unlist(testRanges))
    consensusIDs <- paste0("consensus_",seq(1,length(reduced)))
    mcols(reduced) <- 
      do.call(cbind,lapply(testRanges,function(x)(reduced %over% x)+0))
    if(method=="majority"){
      reducedConsensus <- reduced[rowSums(as.data.frame(mcols(reduced))) > length(testRanges)/2,]
    }
    if(method=="none"){
      reducedConsensus <- reduced
    }
    if(is.numeric(method)){
      reducedConsensus <- reduced[rowSums(as.data.frame(mcols(reduced))) > method,]
    }
    consensusIDs <- paste0("consensus_",seq(1,length(reducedConsensus)))
    mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)),consensusIDs)
    return(reducedConsensus)
    
  }
}

consensusToCount <- runConsensusRegions1(GRangesList(myPeaks), "none")
consensusToCount <- as.data.frame(consensusToCount)
consensusToCount$seqnames <- paste("chr", consensusToCount$seqnames, sep = "")
consensusToCount <- makeGRangesFromDataFrame(consensusToCount, keep.extra.columns = T)

library(Signac) #blacklist_mm10 is here
consensusToCount <- consensusToCount[!consensusToCount %over% blacklist_mm10 
                                     & !seqnames(consensusToCount) %in% "chrM"]

save(consensusToCount, file = "dif no old 3/consensusToCount_sirt6.RData")                    ###go to server


################################################################################
library(DESeq2)
load('dif no old 3/countsFromATAC_sirt6_no_old3.RData')
# load("dif no old 3/consensusToCount_sirt6.RData")
metaData <- data.frame(Group, row.names = colnames(myCounts))
atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount)
atacDDS <- DESeq(atacDDS)
save(atacDDS, file = "dif no old 3/atacDDS_all.RData")


################################################################################
load("dif no old 3/atacDDS_all.RData")

#PCA plot
atac_Rlog <- rlog(atacDDS)
# p <- plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog), returnData = F)

library(ggplot2)
library(ggfortify)
library(ggrepel)
# library(ReportingTools)

# atac_Rlog_df <- makeDESeqDF(atacDDS)

plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog), returnData = F) +
  # geom_label(aes(label = name), size = 3) +
  scale_color_manual(values = c('red', 'hotpink4', 'steelblue1', '#00FF7F'),
                     aes(name = 'Group')) +
  ggtitle('PCA plot: without sample "Old3"') +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        panel.background = element_rect(colour = "darkgrey")) +
  geom_text_repel(hjust=-0.25, vjust=0, size = 3.5, color = 'black',
            label = c('KO1', 'KO2', 'KO3', 'old1', 'old2', 
                      'WT1', 'WT2', 'WT3', 'young1', 'young2', 'young3'),
            box.padding = 0.25) +
  theme_bw()

# ggsave('PCA_sirt6_all.png', height = 6, width = 7.5, dpi = 600)
ggsave('PCA_sirt6_no_old_3_new.png', height = 5, width = 7, dpi = 500)

library(BSgenome.Mmusculus.UCSC.mm10)
library(tracktables)

WT_KO_dif <- results(atacDDS, c("Group", "KO", "WT"), format = "GRanges")
WT_KO_dif <- WT_KO_dif[order(WT_KO_dif$pvalue)]
save(WT_KO_dif, file='dif no old 3/WT_KO_dif.RData')

young_old_dif <- results(atacDDS, c("Group", "old", "young"), format = "GRanges")
young_old_dif <- young_old_dif[order(young_old_dif$pvalue)]
save(young_old_dif, file='dif no old 3/young_old_dif.RData')


WT_KO_dif_df <- as.data.frame(WT_KO_dif)
young_old_dif_df <- as.data.frame(young_old_dif)   #NA in padj


# Select significant values (without old3)

nrow(WT_KO_dif_df) #66318
nrow(WT_KO_dif_df[WT_KO_dif_df$padj<0.05,]) #537
nrow(WT_KO_dif_df[WT_KO_dif_df$pvalue<0.05,]) #7306

nrow(young_old_dif_df) #66318
nrow(young_old_dif_df[young_old_dif_df$padj<0.05,]) #18283
nrow(young_old_dif_df[young_old_dif_df$pvalue<0.05,]) #11483

nrow(WT_KO_dif_df[WT_KO_dif_df$padj<0.05 & WT_KO_dif_df$log2FoldChange > log2(1.5),]) #235
nrow(young_old_dif_df[young_old_dif_df$padj<0.05 & young_old_dif_df$log2FoldChange > log2(1.5),]) #3287



#Select values in SIRT6 locus

start <- 81619787 - 1000000
end <- 81629797 + 1000000

# load('dif no old 3/WT_KO_dif.RData')
# load('dif no old 3/young_old_dif.RData')

WT_KO_dif_df <- as.data.frame(WT_KO_dif)
young_old_dif_df <- as.data.frame(young_old_dif)

nrow(WT_KO_dif_df) #66318
nrow(young_old_dif_df) #66318

WT_KO_dif_df <- na.omit(WT_KO_dif_df)  #there were no missing values
young_old_dif_df <- na.omit(young_old_dif_df)  #12858 rows were removed

nrow(WT_KO_dif_df) #66318
nrow(young_old_dif_df) #53460


# WT_KO_locus <- WT_KO_dif_df[WT_KO_dif_df$padj<0.05 & WT_KO_dif_df$log2FoldChange > log2(1.5),]
WT_KO_locus <- filter(WT_KO_locus, padj<0.05)
# WT_KO_locus <- WT_KO_dif_df[WT_KO_dif_df$padj<0.05,] #without lfc
WT_KO_locus$st_sirt <- start
WT_KO_locus$st_max <- ifelse(WT_KO_locus$start > WT_KO_locus$st_sirt, WT_KO_locus$start, WT_KO_locus$st_sirt)
WT_KO_locus$end_sirt <- end
WT_KO_locus$end_min <- ifelse(WT_KO_locus$end > WT_KO_locus$end_sirt, WT_KO_locus$end_sirt, WT_KO_locus$end)
WT_KO_locus$over <- WT_KO_locus$end_min - WT_KO_locus$st_max
WT_KO_peaks_in_locus <- WT_KO_locus[WT_KO_locus$over > 0,]
nrow(WT_KO_peaks_in_locus) #2 - no in chr10

write.csv2(WT_KO_peaks_in_locus, file = 'dif no old 3/WT_KO_peaks_in_locus.csv', quote = F)


# young_old_locus <- young_old_dif_df[young_old_dif_df$padj<0.05 & young_old_dif_df$log2FoldChange > log2(1.5),]
young_old_locus <- filter(young_old_dif_df, padj<0.05)
# young_old_locus <- young_old_dif_df[young_old_dif_df$padj<0.05,]  #without lfc
young_old_locus$st_sirt <- start
young_old_locus$st_max <- ifelse(young_old_locus$start > young_old_locus$st_sirt, young_old_locus$start, young_old_locus$st_sirt)
young_old_locus$end_sirt <- end
young_old_locus$end_min <- ifelse(young_old_locus$end > young_old_locus$end_sirt, young_old_locus$end_sirt, young_old_locus$end)
young_old_locus$over <- young_old_locus$end_min - young_old_locus$st_max
young_old_peaks_in_locus <- young_old_locus[young_old_locus$over > 0,]
nrow(young_old_peaks_in_locus) #60
young_old_peaks_in_locus <- young_old_peaks_in_locus[young_old_peaks_in_locus$seqnames == 'chr10',]
nrow(young_old_peaks_in_locus) #7 - chr10

write.csv2(young_old_peaks_in_locus, file = 'dif no old 3/young_old_peaks_in_locus.csv', quote = F)



######------All dif peaks combo------######
KO_WT_dif <- results(atacDDS, c("Group", "KO", "WT"), format = "GRanges")
KO_WT_dif <- KO_WT_dif[order(KO_WT_dif$pvalue)]
# save(KO_WT_dif, file='dif no old 3/dif peaks/KO_WT_dif.RData')


KO_old_dif <- results(atacDDS, c("Group", "KO", "old"), format = "GRanges")
KO_old_dif <- KO_old_dif[order(KO_old_dif$pvalue)]
# save(KO_old_dif, file='dif no old 3/dif peaks/KO_old_dif.RData')

KO_young_dif <- results(atacDDS, c("Group", "KO", 'young'), format = "GRanges")
KO_young_dif <- KO_young_dif[order(KO_young_dif$pvalue)]
# save(KO_young_dif, file='dif no old 3/dif peaks/KO_young_dif.RData')

old_young_dif <- results(atacDDS, c("Group", "old", 'young'), format = "GRanges")
old_young_dif <- old_young_dif[order(old_young_dif$pvalue)]
# save(old_young_dif, file='dif no old 3/dif peaks/old_young_dif.RData')

old_WT_dif <- results(atacDDS, c("Group", "old", 'WT'), format = "GRanges")
old_WT_dif <- old_WT_dif[order(old_WT_dif$pvalue)]
# save(old_WT_dif, file='dif no old 3/dif peaks/old_WT_dif.RData')

young_WT_dif <- results(atacDDS, c("Group", "young", 'WT'), format = "GRanges")
young_WT_dif <- young_WT_dif[order(young_WT_dif$pvalue)]
# save(young_WT_dif, file='dif no old 3/dif peaks/young_WT_dif.RData')






load('dif no old 3/dif peaks/KO_old_dif.RData')
sample <- as.data.frame(KO_old_dif)
sample <- filter(sample, padj<0.05)
nrow(sample)







