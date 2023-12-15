setwd('/Users/annaponomareva/Documents/Research/SIRT6_enh/')

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


######------Annotation of dif peaks + Barplot------######

anno_atac_dif <- c()
bar_atac_dif <- c()

atac_peaks <- dir(path = '/Users/annaponomareva/Documents/Research/SIRT6_enh/dif no old 3/dif peaks/', 
                  pattern = "*.RData", full.names = TRUE)
difpeak_names <- c('KO_old_dif', 'KO_WT_dif', 'KO_young_dif', 'old_WT_dif', 'old_young_dif', 'young_WT_dif')

new_colors_11 <- c("#00CD00", "bisque", '#FF7F50', '#FFA07A', 'gold1',
                      'lightgoldenrod1', 'lightskyblue', '#FF83FA', 'seashell4','#CDC5BF','seashell2')
new_colors_10 <- c("#00CD00", '#FF7F50', '#FFA07A', 'gold1',
                      'lightgoldenrod1', 'lightskyblue', '#FF83FA', 'seashell4','#CDC5BF','seashell2')                     

for (i in 1:length(atac_peaks)) {
  sample <- as.data.frame(get(load(atac_peaks[i])))
  sample <- filter(sample, padj<0.05)
  
  if (nrow(sample) > 0) {
    sample <- makeGRangesFromDataFrame(sample, keep.extra.columns = T)
    sample_Anno <- annotatePeak(sample, tssRegion = c(-3000, 3000), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
    
    anno_atac_dif <- append(anno_atac_dif, sample_Anno)
    
    if (nrow(sample_Anno@anno@elementMetadata[grepl("Downstream", sample_Anno@anno@elementMetadata$annotation) == T,])>0) {
      sample_bar <- plotAnnoBar(sample_Anno) +
        ggtitle(difpeak_names[i]) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = new_colors_11, 
                          name = "Features")
      
      bar_atac_dif[[i]] <- sample_bar
    } else {
      sample_bar <- plotAnnoBar(sample_Anno) +
        ggtitle(difpeak_names[i]) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = new_colors_10, 
                          name = "Features")
      
      bar_atac_dif[[i]] <- sample_bar
    }
    
  } else {
    anno_atac_dif <- append(anno_atac_dif, 'no peaks')
    bar_atac_dif[[i]] <- 'no peaks'
  }
}

names(anno_atac_dif) <- difpeak_names   #numbers 1-6 are saved
names(bar_atac_dif) <- difpeak_names

# save(anno_atac_dif, file = 'dif no old 3/anno_atac_dif.RData')
# save(bar_atac_dif, file = 'dif no old 3/bar_atac_dif.RData')

load('dif no old 3/anno_atac_dif.RData')
load('annotation ATAC/anno_atac.RData')

new_colors_11_1 <- c("bisque","#00CD00", '#FF7F50', '#FFA07A', 'gold1',
                            'lightgoldenrod1', 'lightskyblue', '#FF83FA', 'seashell4','#CDC5BF','seashell2')

anno_whole_dif <- plotAnnoBar(anno_atac_dif[c(2,3,1,4,5)]) +           #all samples in one picture
  ggtitle('Feature distribution') +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = new_colors_11_1, 
                    name = "Feature")

ggsave('anno_whole_dif.png', plot = anno_whole_dif, 
       width = 7, height = 6, dpi = 300)



######------GO Annotation------######
#Genome preparing
library(org.Mm.eg.db)
library(annotate)
library(stringr)
library(biomaRt)

egENSEMBL <- toTable(org.Mm.egENSEMBL)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
listAttributes(mart)

mm_genome <- getBM(attributes = c('ensembl_gene_id', "description",
                                  "chromosome_name", "start_position",
                                  "end_position", 'ensembl_transcript_id'), 
                   filters = 'ensembl_gene_id', 
                   values = egENSEMBL$ensembl_id, 
                   mart = mart)                     #Annotation of mouse genes

mm_genome$chromosome_name <- paste("chr", mm_genome$chromosome_name ,sep = "")
mm_genome$width <- mm_genome$end_position - mm_genome$start_position
mm_ensemble <- mm_genome %>% dplyr::select(chromosome_name, start_position, end_position, width,
                                           ensembl_gene_id, ensembl_transcript_id, description)

#Universe preparing 
transcript_universe <- c()
for (i in 1:length(anno_atac)) {
  sample_i <- anno_atac[[i]]@anno
  sample_i <- sample_i[!sample_i %over% blacklist_mm10]
  sample_i <- as.data.frame(sample_i)
  sample_i$transcriptId <- gsub("\\..*","",sample_i$transcriptId)
  sample_tr <- unique(sample_i$transcriptId)
  transcript_universe <- c(transcript_universe, sample_tr)
}

transcript_universe <- unique(transcript_universe)      #30008 transcripts
genes_universe <- getBM(attributes = c('ensembl_gene_id', "ensembl_transcript_id"),
                      filters = "ensembl_transcript_id",
                      values = transcript_universe,
                      mart = mart)
genes_universe <- unique(genes_universe$ensembl_gene_id)     #15454 genes - the nearest genes to all annotated peaks
                                                             #too strict p_adj
library(readxl)
table_genes <- read_excel('41419_2022_5542_MOESM2_ESM.xlsx', 
                          sheet = "Suppl. Table 1", range = cell_rows(2:53702))     #from article
table_genes <- table_genes[table_genes$`WT-1`+table_genes$`WT-2`+table_genes$`WT-3`+table_genes$`WT-4`>4,]
genes_universe2 <- unique(table_genes$`Ensembl ID`)       #21717 genes - from article (Dima)

genes_universe3 <- unique(c(genes_universe, genes_universe2))    #union 1 + 2
save(genes_universe3, file = 'genes_universe3.RData')

#GO plots
library(clusterProfiler)

library(Signac)
head(blacklist_mm10)

go_atac_data <- anno_atac_dif[c(2,3,1,4,5)]
go_names <- difpeak_names[c(2,3,1,4,5)]
names(go_atac_data) <- c()

go_plots <- c()
egos <- c()

for (i in 1:length(go_atac_data)) {
  sample_data <- go_atac_data[[i]]@anno
  sample_data <- sample_data[!sample_data %over% blacklist_mm10]
  sample_data <- as.data.frame(sample_data)
  sample_data$transcriptId <- gsub("\\..*","",sample_data$transcriptId)
  
  sample_genes <- getBM(attributes = c('ensembl_gene_id', "description",
                                      "chromosome_name", "start_position",
                                      "end_position", "ensembl_transcript_id"),
                       filters = "ensembl_transcript_id",
                       values = sample_data$transcriptId,
                       mart = mart)
  
  ego_i <- enrichGO(gene = sample_genes$ensembl_gene_id,
                       universe = genes_universe3, 
                       keyType = "ENSEMBL",
                       OrgDb = org.Mm.eg.db, 
                       ont = "all",
                       # qvalueCutoff = 1,
                       # pvalueCutoff = 1
  )
  
  if (nrow(ego_i@result) != 0) {
    ego_i <- clusterProfiler::simplify(ego_i)
    egos[[i]] <- ego_i
    go_plot_i <- barplot(ego_i, split = "ONTOLOGY", showCategory = 5,
                         title = go_names[i],
                         font.size = 10,
                         # cols.use = c("lightgrey", "blue"),
    ) +
      scale_fill_continuous(name = 'FDR p-value') +
      facet_grid(ONTOLOGY~., scale = "free")
    
    go_plots[[i]] <- go_plot_i
    
  } else {
    egos[[i]] <- ego_i
    go_plots[[i]] <- 'no enriched terms'
  }
}

for (i in 1:length(go_plots)){
  if (i == 2) {
  }
  else {
    ggsave(paste(go_names[i], '_GOplot.png'), plot = go_plots[[i]], 
         width = 5.5, height = 6.7, dpi = 300)
  }
}

names(egos) <- go_names
save(egos, file = 'egos.RData')
save(go_plots, file = 'go_plots.RData')

load('egos.RData')



# ###---Only BP---###
# res_p_e_i <- list()
# for (i in 1:length(anno_atac)) {
#   x <- as.data.frame(anno_atac[[i]])
#   x_p <- x[grepl("Exon", x$annotation) | grepl("Intron", x$annotation) | grepl("Promoter", x$annotation),]
#   res_p_e_i[[i]] <- x_p
# }
# names(res_p_e_i) <- names(anno_atac)    #promoter+exon+intron






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

# narrow_df <- data.frame (R1 = c(nrow(narrowPeaks$KO_R1), nrow(narrowPeaks$WT_R1), nrow(narrowPeaks$young_R1), nrow(narrowPeaks$old_R1)),
#                          R2 = c(nrow(narrowPeaks$KO_R2), nrow(narrowPeaks$WT_R2), nrow(narrowPeaks$young_R2), nrow(narrowPeaks$old_R2)),
#                          R3 = c(nrow(narrowPeaks$KO_R3), nrow(narrowPeaks$WT_R3), nrow(narrowPeaks$young_R3), NA),
#                          row.names = c('KO', 'WT', 'young', 'old'))

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




######------ANOVA: promoter VS intron and distal intergenic regions------######

intron_distal_df$group <- 'DI_intron'
promoter_df$group <- 'promoter'

anova_df <- rbind(intron_distal_df, promoter_df)

anova_res <- aov(num_peaks ~ group + type, data = anova_df)
summary(anova_res) #not significant

anova_res_promoter <- aov(num_peaks ~ type, data = promoter_df)
summary(anova_res_promoter) #not significant

anova_res_DI_intron <- aov(num_peaks ~ type, data = intron_distal_df)
summary(anova_res_DI_intron) #not significant




######------Peaks VS Genes------######
# load('dif no old 3/dif peaks/KO_WT_dif.RData')
# KO_WT_dif_df <- as.data.frame(KO_WT_dif)

KO_WT_dif_anno <- as.data.frame(anno_atac_dif$KO_WT_dif@anno)
KO_WT_dif_anno$transcriptId <- gsub("\\..*","",KO_WT_dif_anno$transcriptId)
rownames(KO_WT_dif_anno) <- 1:nrow(KO_WT_dif_anno)
nrow(KO_WT_dif_anno)     #537

KO_WT_gene_tr <- getBM(attributes = c('ensembl_gene_id', "ensembl_transcript_id"),
                     filters = "ensembl_transcript_id",
                     values = KO_WT_dif_anno$transcriptId,
                     mart = mart,
                     uniqueRows = FALSE)

nrow(KO_WT_gene_tr)    #520 (unique transcripts)

KO_WT_tr <- data.frame(ensembl_gene_id = 1,
                       ensembl_transcript_id = KO_WT_dif_anno$transcriptId)
KO_WT_tr <-KO_WT_tr[order(KO_WT_tr$ensembl_transcript_id),]
KO_WT_tr[duplicated(KO_WT_tr),]

KO_WT_gene_tr <- rbind(KO_WT_gene_tr, KO_WT_tr[duplicated(KO_WT_tr),])
nrow(KO_WT_gene_tr)        #537!!! -> order by transcript name and work on "1"

KO_WT_gene_tr <- KO_WT_gene_tr[order(KO_WT_gene_tr$ensembl_transcript_id),]
for (i in 2:nrow(KO_WT_gene_tr)) {
  if (KO_WT_gene_tr$ensembl_gene_id[i] == 1) {
    KO_WT_gene_tr$ensembl_gene_id[i] <- KO_WT_gene_tr$ensembl_gene_id[i-1]
  }
}
KO_WT_gene_tr[duplicated(KO_WT_gene_tr),]

KO_WT_gene_tr <- KO_WT_gene_tr[order(KO_WT_gene_tr$ensembl_transcript_id),]
KO_WT_dif_anno <- KO_WT_dif_anno[order(KO_WT_dif_anno$transcriptId),]
KO_WT_dif_anno$ensembl_gene_id <- KO_WT_gene_tr$ensembl_gene_id

table_gene_3 <- read_excel('41419_2022_5542_MOESM2_ESM.xlsx', 
                           sheet = "Suppl. Table 3", range = cell_rows(2:18886))     #from article

gene_from_transcript <- function(anno_df) {
  anno_df$transcriptId <- gsub("\\..*","",anno_df$transcriptId)
  rownames(anno_df) <- 1:nrow(anno_df)
  
  df_gene_tr <- getBM(attributes = c('ensembl_gene_id', "ensembl_transcript_id"),
                         filters = "ensembl_transcript_id",
                         values = anno_df$transcriptId,
                         mart = mart)
  
  
  # df_tr <- data.frame(ensembl_gene_id = 1,
  #                        ensembl_transcript_id = anno_df$transcriptId)
  # df_tr <-df_tr[order(df_tr$ensembl_transcript_id),]
  # df_gene_tr <- rbind(df_gene_tr, df_tr[duplicated(df_tr),])
  # 
  # df_gene_tr <- df_gene_tr[order(df_gene_tr$ensembl_transcript_id),]
  # 
  # for (i in 1:nrow(df_gene_tr)) {
  #   if (df_gene_tr$ensembl_gene_id[i] == 1) {
  #     df_gene_tr$ensembl_gene_id[i] <- df_gene_tr$ensembl_gene_id[i-1]
  #   }
  # }
  # 
  # df_gene_tr <- df_gene_tr[order(df_gene_tr$ensembl_transcript_id),]
  # anno_df <- anno_df[order(anno_df$transcriptId),]
  # anno_df$ensembl_gene_id <- df_gene_tr$ensembl_gene_id
  
  anno_df$ensembl_gene_id <- NA
  anno_df$geneLFC <- NA
  
  for (i in 1:nrow(anno_df)) {
    if (anno_df$transcriptId[i] %in% df_gene_tr$ensembl_transcript_id) {
      anno_df$ensembl_gene_id[i] <- df_gene_tr[df_gene_tr$ensembl_transcript_id == anno_df$transcriptId[i],]$ensembl_gene_id
    }
    if (anno_df$ensembl_gene_id[i] %in% table_gene_3$`Ensembl ID`) {
      anno_df$geneLFC[i] <- as.numeric(table_gene_3[table_gene_3$`Ensembl ID` == anno_df[i,]$ensembl_gene_id,]$`log2(Fold Change)`)
    }
  }
  
  return(anno_df)
}

KO_WT_dif_anno_LFC <- gene_from_transcript(anno_df = as.data.frame(anno_atac_dif$KO_WT_dif@anno))
ggplot(KO_WT_dif_anno_LFC, aes(log2FoldChange, geneLFC)) +
  geom_point() +
  xlim(-3.5,4) +
  ggtitle('KO vs WT') +
  xlab('peakLFC')

ggsave(('pictures/LFC_KO_WT_dif.png'), 
       width = 8.3, height = 6.67, dpi = 300)             #table 3 is only for KOvsWT


KO_WT_dif_anno_LFC_abs <- KO_WT_dif_anno_LFC[KO_WT_dif_anno_LFC$log2FoldChange*KO_WT_dif_anno_LFC$geneLFC > 0,]
ggplot(KO_WT_dif_anno_LFC_abs, aes(log2FoldChange, geneLFC)) +
  geom_point() +
  # geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  # geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_segment(aes(x = -2.5, y = 0, xend = 4.5, yend = 0), arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "black") +
  geom_segment(aes(x = 0, y = -1.5, xend = 0, yend = 4), arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "black") +
  xlim(-2.5,4.5) +
  ylim(-1.5,4) +
  ggtitle('KO vs WT') +
  xlab('peakLFC')

ggsave(('pictures/LFC_KO_WT_dif_abs.png'), 
       width = 6.6, height = 6.6, dpi = 300)             #table 3 is only for KOvsWT


# old_young_dif_anno_LFC <- gene_from_transcript(anno_df = as.data.frame(anno_atac_dif$old_young_dif@anno))
# ggplot(old_young_dif_anno_LFC, aes(log2FoldChange, geneLFC)) +
#   geom_point() +
#   xlim(-2.5,2) +
#   ggtitle('old vs young') +
#   xlab('peakLFC')
# 
# ggsave(('pictures/LFC_old_young_dif.png'), 
#        width = 8.3, height = 6.67, dpi = 300)
# 
# old_WT_dif_anno_LFC <- gene_from_transcript(anno_df = as.data.frame(anno_atac_dif$old_WT_dif@anno))
# ggplot(old_WT_dif_anno_LFC, aes(log2FoldChange, geneLFC)) +
#   geom_point() +
#   xlim(-3,2.5) +
#   ggtitle('old vs WT') +
#   xlab('peakLFC')
# 
# ggsave(('pictures/LFC_old_WT_dif.png'), 
#        width = 8.3, height = 6.67, dpi = 300)




######------Number of differential peaks------######

load('dif no old 3/anno_atac_dif.RData')
names(anno_atac_dif)

# narrow_df <- data.frame (R1 = c(nrow(narrowPeaks$KO_R1), nrow(narrowPeaks$WT_R1), nrow(narrowPeaks$young_R1), nrow(narrowPeaks$old_R1)),
#                          R2 = c(nrow(narrowPeaks$KO_R2), nrow(narrowPeaks$WT_R2), nrow(narrowPeaks$young_R2), nrow(narrowPeaks$old_R2)),
#                          R3 = c(nrow(narrowPeaks$KO_R3), nrow(narrowPeaks$WT_R3), nrow(narrowPeaks$young_R3), NA),
#                          row.names = c('KO', 'WT', 'young', 'old'))

narrow_df <- data.frame(sample = names(anno_atac_dif)[1:5],
                        num_peaks = c(length(anno_atac_dif$KO_old_dif@anno),
                                      length(anno_atac_dif$KO_WT_dif@anno),
                                      length(anno_atac_dif$KO_young_dif@anno),
                                      length(anno_atac_dif$old_WT_dif@anno),
                                      length(anno_atac_dif$old_young_dif@anno)))
# narrow_df$type <- factor(narrow_df$type, levels = c('WT', 'KO', 'young', 'old'))

ggplot(narrow_df, aes(x=sample, y=num_peaks, fill = sample)) +
  geom_bar(stat="identity") +
  ggtitle('Number of differential ATAC-seq peaks') +
  ylab('Number of peaks') +
  xlab('Pair of comparison') +
  # ylim(c(7500, 42500)) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="RdYlBu") +
  geom_text(aes(label=num_peaks), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=5.5) +
  theme(panel.background = element_rect(color='grey',fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        plot.title = element_text(size=16)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_x_discrete(labels=c("KO_old_dif" = "KO vs Old", "KO_WT_dif" = "KO vs WT", 'KO_young_dif' = 'KO vs Young',
                            'old_WT_dif' = 'Old vs WT', 'old_young_dif' = 'Old vs Young'))
  

ggsave(('pictures/Differential peak numbers.png'),
       width = 6.5, height = 6.2, dpi = 500)

