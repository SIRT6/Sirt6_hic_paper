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




######------GO Annotation of differential peaks------######
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
# save(egos, file = 'egos.RData')
# save(go_plots, file = 'go_plots.RData')
