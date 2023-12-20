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
  
  ego_i <- enrichGO(gene = sample_genes$ensembl_gene_id,
                    universe = genes_universe3,
                    keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db,
                    ont = "all")
                    # qvalueCutoff = 1,
                    # pvalueCutoff = 1)
                    
  # ego_i <- enrichGO(gene = sample_genes$ensembl_gene_id,
  #                   universe = genes_universe3, 
  #                   keyType = "ENSEMBL",
  #                   OrgDb = org.Mm.eg.db, 
  #                   ont = "BP")
  
  if (nrow(ego_i@result) != 0) {
    ego_i <- clusterProfiler::simplify(ego_i)
    egos_p[[i]] <- ego_i
    go_plot_i <- barplot(ego_i, split = "ONTOLOGY", showCategory = 5,
                         title = go_names[i],
                         font.size = 10,
                         # cols.use = c("lightgrey", "blue"),
    ) +
      scale_fill_continuous(name = 'FDR p-value') +
      facet_grid(ONTOLOGY~., scale = "free")
    
    # go_plot_i <- barplot(ego_i, showCategory = 5,
    #                      title = go_names_plot[i],
    #                      font.size = 10) +
    #   scale_fill_continuous(name = 'FDR p-value')
    
    go_plots_p[[i]] <- go_plot_i
    
  } else {
    egos_p[[i]] <- ego_i
    go_plots_p[[i]] <- 'no enriched terms'
  }
}

names(egos_p) <- go_names
names(go_plots_p) <- go_names



###---Plots: Genes vs Peaks---###

library(ggrepel)

abs_plots <- function(ego_result, categories, anno_result, name_of_dif) {
  
  plots_categories_up_down <- c()
  genes_datasets <- c()
  
  for (i in 1:length(categories)) {
    a <- as.data.frame(ego_result)
    a <- a[a$Description == categories[i],]
    go_genes <- strsplit(a$geneID, split = '/')
    go_genes <- go_genes[[1]]     # list of genes in categories[i]
    
    anno_result_df <- as.data.frame(anno_result)
    anno_result_df$transcriptId <- gsub("\\..*","",anno_result_df$transcriptId)
    rownames(anno_result_df) <- 1:nrow(anno_result_df)
    
    category_df <- getBM(attributes = c('ensembl_gene_id', "ensembl_transcript_id"),
                         filters = "ensembl_gene_id",
                         values = go_genes,
                         mart = mart)
    
    category_trID <- category_df$ensembl_transcript_id
    
    category_df <- filter(anno_result_df, transcriptId %in% category_trID)
    
    df_gene_tr <- getBM(attributes = c('ensembl_gene_id', "ensembl_transcript_id"),
                        filters = "ensembl_transcript_id",
                        values = category_df$transcriptId,
                        mart = mart)
    
    category_df$ensembl_gene_id <- NA
    category_df$geneLFC <- NA
    
    for (ii in 1:nrow(category_df)) {
      if (category_df$transcriptId[ii] %in% df_gene_tr$ensembl_transcript_id) {
        category_df$ensembl_gene_id[ii] <- df_gene_tr[df_gene_tr$ensembl_transcript_id == category_df$transcriptId[ii],]$ensembl_gene_id
      }
      if (category_df$ensembl_gene_id[ii] %in% table_gene_3$`Ensembl ID`) {
        category_df$geneLFC[ii] <- as.numeric(table_gene_3[table_gene_3$`Ensembl ID` == category_df[ii,]$ensembl_gene_id,]$`log2(Fold Change)`)
      }
    }
    
    category_df <- category_df[category_df$log2FoldChange*category_df$geneLFC > 0,]   # if need all genes (not only co-directional) - remove the line
    # pearson_up <- cor(category_df[category_df$log2FoldChange>0,]$log2FoldChange,
    #                   category_df[category_df$log2FoldChange>0,]$geneLFC, method = 'pearson')
    # 
    # pearson_down <- cor(category_df[category_df$log2FoldChange<0,]$log2FoldChange,
    #                     category_df[category_df$log2FoldChange<0,]$geneLFC, method = 'pearson')
    plot_up_down_i <- ggplot(category_df, aes(log2FoldChange, geneLFC, label = ensembl_gene_id)) +
      geom_point(size = 2.2) +
      ggtitle(paste(categories[i])) +
      geom_segment(aes(x = -2.5, y = 0, xend = 2.5, yend = 0), arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "black") +
      geom_segment(aes(x = 0, y = -1.5, xend = 0, yend = 2), arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "black") +
      xlim(-2.5,2.5) +
      ylim(-1.5,2) +
      xlab('peakLFC') +
      theme(plot.title = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(axis.title.x = element_text(size = 10)) +
      # geom_text(hjust=-0.1, vjust=0.5, size = 2.5) +
      geom_text_repel(aes(label = ensembl_gene_id),
                      box.padding   = 0.5, 
                      point.padding = 0.1,
                      segment.size  = 0.1,
                      size = 2.5)
    
    plots_categories_up_down[[i]] <- plot_up_down_i
    
    genes_datasets[[i]] <- category_df
  }
  
  names(genes_datasets) <- categories
  
  result = list(genes = genes_datasets,
                plots = plots_categories_up_down)
  
  return(result)
}

ego_result <- egos_p$KO_WT_dif
x <- ego_result@result

library(readxl)
table_gene_3 <- read_excel('/Users/annaponomareva/Documents/Research/SIRT6_enh/41419_2022_5542_MOESM2_ESM.xlsx', 
                           sheet = "Suppl. Table 3", range = cell_rows(2:18886))     #from article


KO_WT_categories_abs_p <- abs_plots(ego_result = egos_p$KO_WT_dif,      
                                  categories = c(x[grepl("neuro", x$Description) == TRUE |       #x[x$ONTOLOGY == 'CC',][c(2:4),]$Description
                                                     grepl("synap", x$Description) == TRUE |
                                                     grepl("dendrit", x$Description) == TRUE |
                                                     grepl("axon", x$Description) == TRUE,]$Description),
                                  anno_result = res_exon_dif$KO_WT_dif,
                                  name_of_dif = 'KO vs WT')

KO_WT_categories_abs_p      #plots for "neuro" categories

library(ggpubr)
ggarrange(KO_WT_categories_abs_p$plots[[1]], KO_WT_categories_abs_p$plots[[2]],
          ncol = 2, nrow = 1)


###---Plots: Unique Genes vs Peaks (Merging categories)---###

unique_KO_WT_cat <- rbind(KO_WT_categories_abs_p$genes[[1]], 
                          KO_WT_categories_abs_p$genes[[2]])

# for (i in 3:length(KO_WT_categories_abs_p$genes)) {                                #if >=3 categories were identified
#   unique_KO_WT_cat <- rbind(unique_KO_WT_cat, KO_WT_categories_abs_p$genes[[i]])
# }

unique_KO_WT_cat <- unique(unique_KO_WT_cat)
unique_KO_WT_cat$category <- NA
for (i in 1:nrow(unique_KO_WT_cat)) {
  if (unique_KO_WT_cat[i,]$ensembl_gene_id %in% KO_WT_categories_abs_p$genes$`neuron to neuron synapse`$ensembl_gene_id) {
    unique_KO_WT_cat[i,]$category <- 'neuron to neuron synapse'
  }
}

unique_KO_WT_cat_names <- getBM(attributes = c('ensembl_gene_id', "external_gene_name"),
                                filters = "ensembl_gene_id",
                                values = unique_KO_WT_cat$ensembl_gene_id,
                                mart = mart)

unique_KO_WT_cat_names[nrow(unique_KO_WT_cat_names)+1,] <- c('ENSMUSG00000030102', 'Itpr1')            # 2 peaks are overlapping with Itpr1 
unique_KO_WT_cat_names <- unique_KO_WT_cat_names[order(unique_KO_WT_cat_names$ensembl_gene_id),]

unique_KO_WT_cat <- unique_KO_WT_cat[order(unique_KO_WT_cat$ensembl_gene_id),]
unique_KO_WT_cat$gene_name <- unique_KO_WT_cat_names$external_gene_name

ggplot(unique_KO_WT_cat, aes(log2FoldChange, geneLFC, label = gene_name)) +
  geom_point(size = 2.2) +
  ggtitle('KO vs WT: "neuro"') +
  geom_segment(aes(x = -2.5, y = 0, xend = 2.5, yend = 0), arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "black") +
  geom_segment(aes(x = 0, y = -1.5, xend = 0, yend = 2), arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "black") +
  xlim(-2.5,2.5) +
  ylim(-1.5,2) +
  xlab('peakLFC') +
  theme(plot.title = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 10)) +
  # geom_text(hjust=-0.1, vjust=0.5, size = 2.5) +
  geom_text_repel(aes(label = gene_name),
                  box.padding   = 0.5, 
                  point.padding = 0.1,
                  segment.size  = 0.1,
                  size = 2.5)
ggsave(('LFC_KO_WT_neuro_unique_3.png'),      # Unique genes, which are overlapping with differential peaks (plot)
       width = 6, height = 6, dpi = 500)

write.csv(unique_KO_WT_cat, file = 'unique_KO_WT_neuro.csv', sep = '\t')  # Unique genes, which are overlapping with differential peaks (file)

