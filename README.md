# SIRT6 Enhancers


## anno_atac.R
1) Annotation of the *narrowPeak files
2) Pictures with initial number of peaks (whole, promoter, intron + distal intergenic)

## dif_atac_sirt6_no_old3.R
1) Differential peaks obtaining (all combinations: KO-WT, KO-old, KO-young, old-WT, old-young, WT-young[no significant peaks])
2) Annotation of the differential peaks
3) Barplot with annotation all differential peaks
4) Pictures with number of differential peaks

## GO_dif.R
1) GO Annotation of differential peaks
2) GO-plots

## promoter_filter.R
Filtering narrow peaks and differential peaks from promoter regions (4 resulting files: narrow/dif peaks + only promoters/without promoters)

## peaks_vs_genes.R
1) Filtering narrow and differential peaks: get promoter+exon+intron regions
2) Plots for co-directional genes from "neuro"-categories (there are only 2 categories) for KO-WT dif peaks (promoter+exon+intron regions)
3) Plot for co-directional unique genes from "neuro"-categories for KO-WT dif peaks (promoter+exon+intron regions)

Pictures were obtained only for KO-WT, because "41419_2022_5542_MOESM2_ESM.xlsx" (data about gene expression from mitohondial article) contains
info only about KO-WT pair

## enhancer_annotation.R
Creating a table with enhancer types (from ChIP-seq data): active, latent, primed, poised, repressed

## enhancer_annotation_overlaps.R
1) Creating a table with enhancer candidates (overlapping narrow and dif peaks with ChIP-seq annotation)
2) Plots with number of potential enhancers (narrow, differential, UP/DOWN)
