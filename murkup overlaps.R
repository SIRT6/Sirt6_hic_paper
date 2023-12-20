setwd('/Users/annaponomareva/Documents/Research/SIRT6_enh/promoters filter/')

library(dplyr)
library(rtracklayer)
library(ChIPseeker)

load('enhancer_types.RData')
anno_atac_no_p <- get(load('anno_atac_no_p.RData'))
anno_atac_dif_no_p <- get(load('anno_atac_dif_no_p.RData'))


library(regioneR)
ls('package:regioneR')


##################------ Active, Primed, Latent + Peaks ------##################
enhancer_types_narrow <- enhancer_types[enhancer_types$enh_type == 'active' |
                                          enhancer_types$enh_type == 'primed' |
                                          enhancer_types$enh_type == 'latent',]
narrow_overlap <- c()
for (i in 1:length(anno_atac_no_p)) {
  peak <- anno_atac_no_p[[i]]
  over <- overlapRegions(peak, enhancer_types, colA = c(6,12,13), colB = c(4))
  over <- over %>% dplyr::rename ('seqnames' = 'chr',
                                  'start' = 'startA',
                                  'end' = 'endA',
                                  'start_mark' = 'startB',
                                  'end_mark' = 'endB')
  narrow_overlap[[i]] <- over
}
names(narrow_overlap) <- names(anno_atac_no_p)
save(narrow_overlap, file = 'narrow_overlap.RData')




##################------ Poised, Repressed + Dif.Peaks ------##################
enhancer_types_dif <- enhancer_types[enhancer_types$enh_type == 'poised' |
                                          enhancer_types$enh_type == 'repressed',]
dif_overlap <- c()
for (i in 1:length(anno_atac_dif_no_p)) {
  peak <- anno_atac_dif_no_p[[i]]
  over <- overlapRegions(peak, enhancer_types, colA = c(6,12,13), colB = c(4))
  over <- over %>% dplyr::rename ('seqnames' = 'chr',
                                  'start' = 'startA',
                                  'end' = 'endA',
                                  'start_mark' = 'startB',
                                  'end_mark' = 'endB')
  dif_overlap[[i]] <- over
}
names(dif_overlap) <- names(anno_atac_dif_no_p)
save(dif_overlap, file = 'dif_overlap.RData')


###################--------------Number of potential enhancers-------------###################

number_overlaps <- data.frame('sample' = '',
                              'type' = '',
                              'num_overlaps' = '')

number_overlaps <- number_overlaps[-1,]

for (i in 1:length(narrow_overlap)) {
  x <- narrow_overlap[[i]]
  new_row_active <- c(names(narrow_overlap)[i], 'active', nrow(x[x$enh_type == 'active',]))
  new_row_primed <- c(names(narrow_overlap)[i], 'primed', nrow(x[x$enh_type == 'primed',]))
  new_row_latent <- c(names(narrow_overlap)[i], 'latent', nrow(x[x$enh_type == 'latent',]))
  new_row_poised <- c(names(narrow_overlap)[i], 'poised', nrow(x[x$enh_type == 'poised',]))
  new_row_repressed <- c(names(narrow_overlap)[i], 'repressed', nrow(x[x$enh_type == 'repressed',]))
  number_overlaps <- rbind(number_overlaps, new_row_active, new_row_primed, new_row_latent, new_row_poised, new_row_repressed)
}

for (i in 1:length(dif_overlap)) {
  x <- dif_overlap[[i]]
  new_row_active <- c(names(dif_overlap)[i], 'active', nrow(x[x$enh_type == 'active',]))
  new_row_primed <- c(names(dif_overlap)[i], 'primed', nrow(x[x$enh_type == 'primed',]))
  new_row_latent <- c(names(dif_overlap)[i], 'latent', nrow(x[x$enh_type == 'latent',]))
  new_row_poised <- c(names(dif_overlap)[i], 'poised', nrow(x[x$enh_type == 'poised',]))
  new_row_repressed <- c(names(dif_overlap)[i], 'repressed', nrow(x[x$enh_type == 'repressed',]))
  number_overlaps <- rbind(number_overlaps, new_row_active, new_row_primed, new_row_latent, new_row_poised, new_row_repressed)
}

colnames(number_overlaps) <- c('sample', 'type', 'num_overlaps')


number_overlaps$no_pr_peaks <- c(rep(nrow(anno_atac_no_p$KO_R1), 5), rep(nrow(anno_atac_no_p$KO_R2),5),
                                 rep(nrow(anno_atac_no_p$KO_R3), 5), rep(nrow(anno_atac_no_p$old_R1), 5),
                                 rep(nrow(anno_atac_no_p$old_R2), 5), rep(nrow(anno_atac_no_p$WT_R1), 5),
                                 rep(nrow(anno_atac_no_p$WT_R2), 5), rep(nrow(anno_atac_no_p$WT_R3), 5),
                                 rep(nrow(anno_atac_no_p$young_R1), 5), rep(nrow(anno_atac_no_p$young_R2), 5),
                                 rep(nrow(anno_atac_no_p$young_R3), 5), rep(nrow(anno_atac_dif_no_p$KO_old_dif), 5),
                                 rep(nrow(anno_atac_dif_no_p$KO_WT_dif), 5), rep(nrow(anno_atac_dif_no_p$KO_young_dif), 5),
                                 rep(nrow(anno_atac_dif_no_p$old_WT_dif), 5), rep(nrow(anno_atac_dif_no_p$old_young_dif), 5))

number_overlaps$num_overlaps <- as.numeric(number_overlaps$num_overlaps)
number_overlaps$fract_overlap <- (number_overlaps$num_overlaps) / (number_overlaps$no_pr_peaks)

save(number_overlaps, file = 'number_overlaps.RData')

load('number_overlaps.RData')
load('enhancer_types.RData')
load('narrow_overlap.RData')
load('dif_overlap.RData')
enhancer_types
number_overlaps


ggplot(number_overlaps[!(number_overlaps$type=='poised'| number_overlaps$type=='repressed'),][1:33,], aes(x=sample, y=fract_overlap, fill = type)) +
  geom_bar(position="dodge", stat="identity") +
  coord_flip() +
  ggtitle('Enhancer candidates') +
  xlab('Sample') +
  ylab('Enhancer fraction') +
  scale_fill_brewer(palette="Accent") +
  geom_text(aes(label=num_overlaps), hjust = 1.1, color="black",
            position = position_dodge(width = 0.9), size=4.2) +
  theme(panel.background = element_rect(color='grey',fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        plot.title = element_text(size=16)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  labs(fill = "Enhancer \ntype")+
  scale_x_discrete(labels=c("KO_R1" = "KO-1", "KO_R2" = "KO-2", "KO_R3" = "KO-3",
                            "WT_R1" = "WT-1", "WT_R2" = "WT-2", "WT_R3" = "WT-3",
                            "old_R1" = "Old-1", "old_R2" = "Old-2",
                            "young_R1" = "Young-1", "young_R2" = "Young-2", "young_R3" = "Young-3"))

ggplot(number_overlaps[56:80,], aes(x=sample, y=fract_overlap, fill = type)) +
  geom_bar(position="dodge", stat="identity") +
  coord_flip() +
  ggtitle('Enhancer candidates') +
  xlab('Sample') +
  ylab('Enhancer fraction') +
  ylim(c(0, 0.5))+
  scale_fill_brewer(palette="Accent") +
  geom_text(aes(label=num_overlaps), hjust = 0, color="black",
            position = position_dodge(width = 0.9), size=4.2) +
  theme(panel.background = element_rect(color='grey',fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        plot.title = element_text(size=16)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  labs(fill = "Enhancer \ntype")+
  scale_x_discrete(labels=c("KO_old_dif" = "KO vs Old", "KO_WT_dif" = "KO vs WT", 'KO_young_dif' = 'KO vs Young',
                            'old_WT_dif' = 'Old vs WT', 'old_young_dif' = 'Old vs Young'))




###################--------------UP/DOWN peaks : potential enhancers-------------###################

dif_overlap_up <- c()
for (i in 1:length(anno_atac_dif_no_p)) {
  peak <- anno_atac_dif_no_p[[i]]
  peak <- peak[peak$log2FoldChange>0,]
  over <- overlapRegions(peak, enhancer_types, colA = c(6,12,13), colB = c(4))
  over <- over %>% dplyr::rename ('seqnames' = 'chr',
                                  'start' = 'startA',
                                  'end' = 'endA',
                                  'start_mark' = 'startB',
                                  'end_mark' = 'endB')
  dif_overlap_up[[i]] <- over
}
names(dif_overlap_up) <- names(anno_atac_dif_no_p)
save(dif_overlap_up, file = 'dif_overlap_up.RData')

dif_overlap_down <- c()
for (i in 1:length(anno_atac_dif_no_p)) {
  peak <- anno_atac_dif_no_p[[i]]
  peak <- peak[peak$log2FoldChange<0,]
  over <- overlapRegions(peak, enhancer_types, colA = c(6,12,13), colB = c(4))
  over <- over %>% dplyr::rename ('seqnames' = 'chr',
                                  'start' = 'startA',
                                  'end' = 'endA',
                                  'start_mark' = 'startB',
                                  'end_mark' = 'endB')
  dif_overlap_down[[i]] <- over
}
names(dif_overlap_down) <- names(anno_atac_dif_no_p)
save(dif_overlap_down, file = 'dif_overlap_down.RData')

number_overlaps_up_down <- data.frame('sample' = '',
                              'type' = '',
                              'num_overlaps' = '',
                              'dif_dir' = '')

number_overlaps_up_down <- number_overlaps_up_down[-1,]

for (i in 1:length(dif_overlap_up)) {
  x <- dif_overlap_up[[i]]
  new_row_active <- c(names(dif_overlap_up)[i], 'active', nrow(x[x$enh_type == 'active',]), 'up')
  new_row_primed <- c(names(dif_overlap_up)[i], 'primed', nrow(x[x$enh_type == 'primed',]), 'up')
  new_row_latent <- c(names(dif_overlap_up)[i], 'latent', nrow(x[x$enh_type == 'latent',]), 'up')
  new_row_poised <- c(names(dif_overlap_up)[i], 'poised', nrow(x[x$enh_type == 'poised',]), 'up')
  new_row_repressed <- c(names(dif_overlap_up)[i], 'repressed', nrow(x[x$enh_type == 'repressed',]), 'up')
  number_overlaps_up_down <- rbind(number_overlaps_up_down, new_row_active, new_row_primed, new_row_latent, new_row_poised, new_row_repressed)
}

for (i in 1:length(dif_overlap_down)) {
  x <- dif_overlap_down[[i]]
  new_row_active <- c(names(dif_overlap_down)[i], 'active', nrow(x[x$enh_type == 'active',]), 'down')
  new_row_primed <- c(names(dif_overlap_down)[i], 'primed', nrow(x[x$enh_type == 'primed',]), 'down')
  new_row_latent <- c(names(dif_overlap_down)[i], 'latent', nrow(x[x$enh_type == 'latent',]), 'down')
  new_row_poised <- c(names(dif_overlap_down)[i], 'poised', nrow(x[x$enh_type == 'poised',]), 'down')
  new_row_repressed <- c(names(dif_overlap_down)[i], 'repressed', nrow(x[x$enh_type == 'repressed',]), 'down')
  number_overlaps_up_down <- rbind(number_overlaps_up_down, new_row_active, new_row_primed, new_row_latent, new_row_poised, new_row_repressed)
}

colnames(number_overlaps_up_down) <- c('sample', 'type', 'num_overlaps', 'dif_dir')

number_overlaps_up_down$no_pr_peaks <- c(rep(nrow(anno_atac_dif_no_p$KO_old_dif), 5),
                                 rep(nrow(anno_atac_dif_no_p$KO_WT_dif), 5), rep(nrow(anno_atac_dif_no_p$KO_young_dif), 5),
                                 rep(nrow(anno_atac_dif_no_p$old_WT_dif), 5), rep(nrow(anno_atac_dif_no_p$old_young_dif), 5),
                                 rep(nrow(anno_atac_dif_no_p$KO_old_dif), 5),
                                 rep(nrow(anno_atac_dif_no_p$KO_WT_dif), 5), rep(nrow(anno_atac_dif_no_p$KO_young_dif), 5),
                                 rep(nrow(anno_atac_dif_no_p$old_WT_dif), 5), rep(nrow(anno_atac_dif_no_p$old_young_dif), 5))

number_overlaps_up_down$num_overlaps <- as.numeric(number_overlaps_up_down$num_overlaps)
number_overlaps_up_down$fract_overlap <- (number_overlaps_up_down$num_overlaps) / (number_overlaps_up_down$no_pr_peaks)

save(number_overlaps_up_down, file = 'number_overlaps_up_down.RData')

ggplot(number_overlaps_up_down[1:25,], aes(x=sample, y=fract_overlap, fill = type)) +
  geom_bar(position="dodge", stat="identity") +
  coord_flip() +
  ggtitle('Enhancer candidates: UP') +
  xlab('Sample') +
  ylab('Enhancer fraction') +
  ylim(c(0, 0.5))+
  scale_fill_brewer(palette="Accent") +
  geom_text(aes(label=num_overlaps), hjust = 0, color="black",
            position = position_dodge(width = 0.9), size=4.2) +
  theme(panel.background = element_rect(color='grey',fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        plot.title = element_text(size=16)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  labs(fill = "Enhancer \ntype")+
  scale_x_discrete(labels=c("KO_old_dif" = "KO vs Old", "KO_WT_dif" = "KO vs WT", 'KO_young_dif' = 'KO vs Young',
                            'old_WT_dif' = 'Old vs WT', 'old_young_dif' = 'Old vs Young'))

ggplot(number_overlaps_up_down[26:50,], aes(x=sample, y=fract_overlap, fill = type)) +
  geom_bar(position="dodge", stat="identity") +
  coord_flip() +
  ggtitle('Enhancer candidates: DOWN') +
  xlab('Sample') +
  ylab('Enhancer fraction') +
  ylim(c(0, 0.5))+
  scale_fill_brewer(palette="Accent") +
  geom_text(aes(label=num_overlaps), hjust = 0, color="black",
            position = position_dodge(width = 0.9), size=4.2) +
  theme(panel.background = element_rect(color='grey',fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        plot.title = element_text(size=16)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  labs(fill = "Enhancer \ntype")+
  scale_x_discrete(labels=c("KO_old_dif" = "KO vs Old", "KO_WT_dif" = "KO vs WT", 'KO_young_dif' = 'KO vs Young',
                            'old_WT_dif' = 'Old vs WT', 'old_young_dif' = 'Old vs Young'))









