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


######------Number of peaks in the initial *narrowPeak files------######

narrowPeaks_path <- dir(path = '/Users/annaponomareva/Documents/Research/SIRT6_enh/narrowPeak_no_old3/', 
                   pattern = "*.narrowPeak", full.names = TRUE)
narrowPeaks <- c()
for (i in 1:length(narrowPeaks_path)) {
  sample <- read.csv(narrowPeaks_path[i], header = FALSE, sep = "\t")
  sample$V1 <- as.character(sample$V1)
  sample$V1 <- paste("chr", sample$V1, sep = "")
  narrowPeaks[[i]] <- sample
}

narrow_names <- c('KO_R1', 'KO_R2', 'KO_R3', 'old_R1', 'old_R2',
                  'WT_R1', 'WT_R2', 'WT_R3', 'young_R1', 'young_R2', 'young_R3')

names(narrowPeaks) <- narrow_names


narrow_df <- data.frame(type = c('WT', 'WT', 'WT', 'KO', 'KO', 'KO', 
                                 'young', 'young', 'young', 'old', 'old'),
                        sample = c('WT_R1', 'WT_R2', 'WT_R3', 'KO_R1', 'KO_R2', 'KO_R3',
                                   'young_R1', 'young_R2', 'young_R3', 'old_R1', 'old_R2'),
                        num_peaks = c(nrow(narrowPeaks$WT_R1), nrow(narrowPeaks$WT_R2), nrow(narrowPeaks$WT_R3),
                                      nrow(narrowPeaks$KO_R1), nrow(narrowPeaks$KO_R2), nrow(narrowPeaks$KO_R3),
                                      nrow(narrowPeaks$young_R1), nrow(narrowPeaks$young_R2), nrow(narrowPeaks$young_R3),
                                      nrow(narrowPeaks$old_R1), nrow(narrowPeaks$old_R2)))
narrow_df$type <- factor(narrow_df$type, levels = c('WT', 'KO', 'young', 'old'))

ggplot(narrow_df, aes(x=type, y=num_peaks, fill=type, 
                      label=c("WT_R1", "WT_R2", "WT_R3", "KO_R1", "KO_R2", "KO_R3", 
                              "young_R1", "young_R2", "young_R3", "old_R1", "old_R2"))) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.7) +
  ggtitle('Number of ATAC-seq peaks') +
  ylab('Number of peaks') +
  xlab('Type of samples') +
  ylim(c(7500, 42500)) +
  scale_fill_brewer(palette="RdYlBu", name = 'Type') +
  geom_text(hjust=-0.4, vjust=0.5, size = 4) +
  theme_bw()+
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  # geom_text(label = num_peaks, hjust=-0.4, vjust=-0.8, size = 3)

ggsave(('pictures/Initial peak numbers new.png'), 
       width = 6.9, height = 7.2, dpi = 500)




######------Number of peaks in the promoter, intron and distal intergenic regions------######

load('annotation ATAC/anno_atac.RData')

promoter_nums <- c()

for (i in 1:length(anno_atac)) {
  sample <- as.data.frame(anno_atac[[i]]@anno)
  sample <- sample[grepl("Promoter", sample$annotation) == TRUE,]
  
  promoter_nums <- append(promoter_nums, nrow(sample))
}

# Order of samples: ('KO_R1', 'KO_R2', 'KO_R3', 'old_R1', 'old_R2', 'WT_R1', 'WT_R2', 'WT_R3', 'young_R1', 'young_R2', 'young_R3')

promoter_df <- data.frame(type = c('KO', 'KO', 'KO','old', 'old',
                                   'WT', 'WT', 'WT', 'young', 'young', 'young'),
                        sample = c('KO_R1', 'KO_R2', 'KO_R3', 'old_R1', 'old_R2',
                                   'WT_R1', 'WT_R2', 'WT_R3', 'young_R1', 'young_R2', 'young_R3'),
                        num_peaks = promoter_nums)
promoter_df$type <- factor(promoter_df$type, levels = c('WT', 'KO', 'young', 'old'))

ggplot(promoter_df, aes(x=type, y=num_peaks, fill=type, label=num_peaks)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.4) +
  ggtitle('Number of ATAC-seq peaks in\n promoter regions') +
  ylab('Number of peaks') +
  xlab('Type of samples') +
  # ylim(c(7500, 42500)) +
  scale_fill_brewer(palette="RdYlBu", name = 'Type') +
  geom_text(hjust=-0.5, vjust=0.5, size = 3)

ggsave(('pictures/Promoter peak numbers.png'), 
       width = 6.8, height = 6.43, dpi = 300)



intron_distal_nums <- c()

for (i in 1:length(anno_atac)) {
  sample <- as.data.frame(anno_atac[[i]]@anno)
  sample <- sample[sample$annotation == "Distal Intergenic" |
                     grepl("Intron", sample$annotation) == TRUE,]
  
  intron_distal_nums <- append(intron_distal_nums, nrow(sample))
}

intron_distal_df <- data.frame(type = c('KO', 'KO', 'KO','old', 'old',
                                   'WT', 'WT', 'WT', 'young', 'young', 'young'),
                          sample = c('KO_R1', 'KO_R2', 'KO_R3', 'old_R1', 'old_R2',
                                     'WT_R1', 'WT_R2', 'WT_R3', 'young_R1', 'young_R2', 'young_R3'),
                          num_peaks = intron_distal_nums)
intron_distal_df$type <- factor(intron_distal_df$type, levels = c('WT', 'KO', 'young', 'old'))

ggplot(intron_distal_df, aes(x=type, y=num_peaks, fill=type, label=num_peaks)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.4) +
  ggtitle('Number of ATAC-seq peaks in\n distal intergenic and intron regions') +
  ylab('Number of peaks') +
  xlab('Type of samples') +
  # ylim(c(7500, 42500)) +
  scale_fill_brewer(palette="RdYlBu", name = 'Type') +
  geom_text(hjust=-0.5, vjust=0.5, size = 3)

ggsave(('pictures/Intron and DistInt peak numbers.png'), 
       width = 6.8, height = 6.43, dpi = 300)


