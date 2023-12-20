setwd('/Users/annaponomareva/Documents/Research/SIRT6_enh/promoters filter/')


######------Read data------######
anno_atac_no_p <- get(load('/Users/annaponomareva/Documents/Research/SIRT6_enh/promoters filter/anno_atac_no_p.RData'))

library(dplyr)
library(rtracklayer)
cortex_data <- as.data.frame(import("chip_data/His.Neu.05.AllAg.Cortical_neuron.bed", format="bed"))
# cortex_data <- cortex_data[grepl('H3K27ac', cortex_data$name) | grepl('H3K27me3', cortex_data$name),]
cortex_data_9me2 <- as.data.frame(import("chip_data/His.Neu.05.H3K9me2.AllCell.bed", format="bed"))
cortex_data_9me3 <- as.data.frame(import("chip_data/His.Neu.05.H3K9me3.AllCell.bed", format="bed"))

SRX13085667 <- read.table('chip_data/SRX13085667.05.bed')        #H3K27ac +
SRX13085667 <- SRX13085667 %>% rename('seqnames' = "V1",
                                  'start' = 'V2',
                                  'end' = 'V3')

SRX4546325 <- read.table('chip_data/SRX4546325.05.bed')      #H3K27me3 +
SRX4546325 <- SRX4546325 %>% rename('seqnames' = "V1",
                                    'start' = 'V2',
                                    'end' = 'V3')
SRX4546326 <- read.table('chip_data/SRX4546326.05.bed')      #H3K27me3 +
SRX4546326 <- SRX4546326 %>% rename('seqnames' = "V1",
                                    'start' = 'V2',
                                    'end' = 'V3')

SRX13085666 <- read.table('chip_data/SRX13085666.05.bed')      #H3K4me1 +
SRX13085666 <- SRX13085666 %>% rename('seqnames' = "V1",
                                      'start' = 'V2',
                                      'end' = 'V3')

SRX2357518 <- read.table('chip_data/SRX2357518.05.bed')      #H3K9me2 +
SRX2357518 <- SRX2357518 %>% rename('seqnames' = "V1",
                                    'start' = 'V2',
                                    'end' = 'V3')
SRX2357528 <- read.table('chip_data/SRX2357528.05.bed')      #H3K9me2 +
SRX2357528 <- SRX2357528 %>% rename('seqnames' = "V1",
                                    'start' = 'V2',
                                    'end' = 'V3')

SRX10291312 <- read.table('chip_data/SRX10291312.05.bed')      #H3K9me3 +
SRX10291312 <- SRX10291312 %>% rename('seqnames' = "V1",
                                    'start' = 'V2',
                                    'end' = 'V3')
SRX10291313 <- read.table('chip_data/SRX10291313.05.bed')      #H3K9me3 +
SRX10291313 <- SRX10291313 %>% rename('seqnames' = "V1",
                                      'start' = 'V2',
                                      'end' = 'V3')
SRX10291314 <- read.table('chip_data/SRX10291314.05.bed')      #H3K9me3 +
SRX10291314 <- SRX10291314 %>% rename('seqnames' = "V1",
                                      'start' = 'V2',
                                      'end' = 'V3')
SRX10291315 <- read.table('chip_data/SRX10291315.05.bed')      #H3K9me3 +
SRX10291315 <- SRX10291315 %>% rename('seqnames' = "V1",
                                      'start' = 'V2',
                                      'end' = 'V3')


###---Union rep1+rep2---###             #remove random and sex chromosomes
library(regioneR)
ls('package:regioneR')

mice_chr <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
              'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrM')

h3k27ac <- as.data.frame(SRX13085667)
h3k27ac <- h3k27ac[h3k27ac$seqnames %in% mice_chr,]

h3k27me3 <- as.data.frame(mergeRegions(SRX4546325, SRX4546326))
h3k27me3 <- h3k27me3[h3k27me3$seqnames %in% mice_chr,]

h3k4me1 <- as.data.frame(SRX13085666)
h3k4me1 <- h3k4me1[h3k4me1$seqnames %in% mice_chr,]

h3k9me2 <- as.data.frame(mergeRegions(SRX2357518, SRX2357528))
h3k9me2 <- h3k9me2[h3k9me2$seqnames %in% mice_chr,]

h3k9me3_1 <- mergeRegions(SRX10291312, SRX10291313)
h3k9me3_2 <- mergeRegions(SRX10291314, SRX10291315)
h3k9me3 <- as.data.frame(mergeRegions(h3k9me3_1, h3k9me3_2))
h3k9me3 <- h3k9me3[h3k9me3$seqnames %in% mice_chr,]




#####################------------------- MARKUP -------------------#####################

enhancer_types <- data.frame('seqnames' = '',
                             'start' = '',
                             'end' = '',
                             'enh_type' = '')
enhancer_types <- enhancer_types[-1,]

latent <- h3k27ac[,c(1,2,3)]
latent$enh_type <- 'latent'
colnames(latent) <- colnames(enhancer_types)

active <- overlapRegions(h3k27ac, h3k4me1)
active$start <- apply(active[, c(2,4)], 1, min)
active$end <- apply(active[, c(3,5)], 1, max)
active <- active[,c(1,7,8)]
active$enh_type <- 'active'
colnames(active) <- colnames(enhancer_types)

primed_1 <- overlapRegions(h3k27ac, h3k4me1)
primed_1 <- primed_1[,c(1,4,5)]     #h3k4me1 overlapping with h3k27ac
h3k4me1_1 <- h3k4me1[,1:3]
h3k4me1_1$name <- paste(h3k4me1_1$seqnames, h3k4me1_1$start, h3k4me1_1$end, sep='')
primed_1$name <- paste(primed_1$chr, primed_1$startB, primed_1$endB, sep='')
primed <- h3k4me1_1[!(h3k4me1_1$name %in% primed_1$name),]
primed <- primed[,1:3]
primed$enh_type <- 'primed'
colnames(primed) <- colnames(enhancer_types)

poised <- overlapRegions(h3k27me3, h3k4me1)
poised$start <- apply(poised[, c(2,4)], 1, min)
poised$end <- apply(poised[, c(3,5)], 1, max)
poised <- poised[,c(1,7,8)]
poised$enh_type <- 'poised'
colnames(poised) <- colnames(enhancer_types)

repressed2 <- overlapRegions(h3k27me3, h3k9me2)
repressed2$start <- apply(repressed2[, c(2,4)], 1, min)
repressed2$end <- apply(repressed2[, c(3,5)], 1, max)
repressed2 <- repressed2[,c(1,7,8)]
repressed2$enh_type <- 'repressed'
colnames(repressed2) <- colnames(enhancer_types)

repressed3 <- overlapRegions(h3k27me3, h3k9me3)
repressed3$start <- apply(repressed3[, c(2,4)], 1, min)
repressed3$end <- apply(repressed3[, c(3,5)], 1, max)
repressed3 <- repressed3[,c(1,7,8)]
repressed3$enh_type <- 'repressed'
colnames(repressed3) <- colnames(enhancer_types)

enhancer_types <- rbind(active, primed, latent, poised, repressed2, repressed3)
enhancer_types <- unique(enhancer_types)
save(enhancer_types, file='enhancer_types.RData')                 # resulting file with enhancer annotation


