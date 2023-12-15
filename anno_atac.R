# setwd("C:/Users/Аня/Desktop/Research/SIRT6_enh/annotation ATAC/")
setwd('/Users/annaponomareva/Documents/Research/SIRT6_enh/')

library(Rsubread)
library(ChIPseeker)
library(GenomicRanges)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)



######------Annotation------######

anno_atac <- c()

atac_peaks <- dir(path = 'C:/Users/Аня/Desktop/Research/SIRT6_enh/narrowPeak_no_old3/', pattern = "*.narrowPeak", full.names = TRUE)

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



######------SIRT6 locus overlaps------######
load('anno_atac.RData')

start <- 81619787 - 1000000
end <- 81629797 + 1000000

over_sitr6_Dist_Intron <- c()

for (i in 1:length(anno_atac)) {
  sample <- as.data.frame(anno_atac[[i]]@anno)
  sample <- sample[sample$annotation == "Distal Intergenic" |
                     grepl("Intron", sample$annotation) == TRUE,]
  sample <- sample[sample$start >= start & sample$end <= end,]
  
  sample <- makeGRangesFromDataFrame(sample, keep.extra.columns = T)
  
  over_sitr6_Dist_Intron[[i]] <- sample
}

names(over_sitr6_Dist_Intron) <- c('KO_R1', 'KO_R2', 'KO_R3', 'old_R1', 'old_R2', 
                      'WT_R1', 'WT_R2', 'WT_R3', 'young_R1', 'young_R2', 'young_R3')   #numbers 1-11 are saved

for (i in 1:length(anno_atac)) {
  print(c(names(over_sitr6_Dist_Intron)[i], nrow(as.data.frame(over_sitr6_Dist_Intron[i]))))    #print nrow()
}

save(over_sitr6_Dist_Intron, file = 'over_sitr6_Dist_Intron.RData')



######------Enhancers overlaps------######
library(regioneR)

enh_mm10_cortex <- read.csv('/Users/annaponomareva/Documents/Research/SIRT6_enh/enh_mm10_sirt6_cortex.bed', header = F, sep='\t')
colnames(enh_mm10_cortex) <- c('seqnames', 'start', 'end')
# enh_mm10_cortex <- makeGRangesFromDataFrame(enh_mm10_cortex, keep.extra.columns = T)

for (i in 1:length(over_sitr6_Dist_Intron)){
  x <- as.data.frame(over_sitr6_Dist_Intron[i])
  x <- x[,c(1,2,3)]
  colnames(x) <- c('seqnames', 'start', 'end')
  print(overlapRegions(x, enh_mm10_cortex))   #all are zeros
}





