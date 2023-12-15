# setwd("C:/Users/Аня/Desktop/Research/SIRT6_enh/annotation ATAC/")
setwd('/Users/annaponomareva/Documents/Research/SIRT6_enh/')

library(Rsubread)
library(ChIPseeker)
library(GenomicRanges)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)



######------Annotation of narrow peaks------######

anno_atac <- c()

atac_peaks <- dir(path = '/Users/annaponomareva/Documents/Research/SIRT6_enh/narrowPeak_no_old3/', pattern = "*.narrowPeak", full.names = TRUE)

for (i in 1:length(atac_peaks)) {
  sample <- read.csv(atac_peaks[i], header = FALSE, sep = "\t")
  sample$V1 <- as.character(sample$V1)
  sample$V1 <- paste("chr", sample$V1, sep = "")
  
  sample_GRange <- makeGRangesFromDataFrame(sample, start.field = "V2", end.field = "V3", seqnames.field = "V1")
  sample_Anno <- annotatePeak(sample_GRange, tssRegion = c(-3000, 3000), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
  
  anno_atac <- append(anno_atac, sample_Anno)
}

names(anno_atac) <- c('KO_R1', 'KO_R2', 'KO_R3', 'old_R1', 'old_R2', 
                      'WT_R1', 'WT_R2', 'WT_R3', 'young_R1', 'young_R2', 'young_R3')   #numbers 1-11 are saved

save(anno_atac, file = 'anno_atac.RData')





