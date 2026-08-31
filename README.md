[![DOI:10.1101/2025.11.10.687599](http://img.shields.io/badge/DOI-10.1101/2025.11.10.687599-B31B1B.svg)](https://doi.org/10.1101/2025.11.10.687599)

# The Divergent 3D Genome Landscapes of Aging and Neurodegenerative mouse models

This repository contains the notebooks to reproduce analyses and panels from the associated preprint.

## Repository structure

The repository is organized by manuscript figure:

| Files | Analysis |
| --- | --- |
| `Figure1_cis_trans_contacts.ipynb` | Cis/trans contact ratios |
| `Figure1_DLR_ICF.ipynb` | Distal-to-local and inter-chromosomal contact fractions |
| `Figure1_nuclear_measurements.ipynb` | Nuclear measurements analysis |
| `Figure1_PCA.ipynb` | PCA of compartment eigenvector values |
| `Figure1_scalings.ipynb` | Contact-probability scaling curves |
| `Figure2_chr_map_saddle_plots.ipynb` | Chromosome maps, saddle plots, and compartment strength |
| `Figure2_chr_map_saddle_plots_comp_strength.ipynb` | Additional compartment-strength analysis |
| `Figure2_diff_compartments.ipynb` | Differential compartment analysis |
| `Figure3_TADs.ipynb` | TAD boundaries, average TADs, and boundary strength |
| `Figure3_average_loop_up_down_folded_loops.ipynb` | Average loops, loop strength, and differential loops |
| `Figure3_loops_TADs.ipynb` | Loop and TAD decile analysis |
| `Figure3_loops_reg_elem.ipynb` | Association of loops with regulatory elements |
| `Figure3_scaling_slopes.ipynb` | Contact-scaling slope analysis |
| `Figure4_compartments.ipynb` | Compartment analysis in aging and neurodegeneration models |
| `Figure4_hilbert_curves.ipynb` | Hilbert-curve and circular visualizations |
| `Figure4_loops.ipynb` | Loop analysis in aging and neurodegeneration models |
| `Figure4_TADs.ipynb` | TAD analysis in aging and neurodegeneration models |
| `Figure5.ipynb` | Transcriptomic, pathway-enrichment, and integrative analyses |
| `FigureS2B.ipynb` | Coverage and mappability analysis for Supplementary Figure S2B |
| `FigureS3G_Loops_and_CTCF.ipynb` | CTCF peak overlap at differential loop anchors, stratified by regulatory-element class |
| `downsampling.py` | Helper script for creating and downsampling Cooler contact matrices |

## Requirements

The notebooks use both Python and R. 

- Python 3.8 or later
- R 4.4 or later, with IRkernel
- Jupyter Notebook or JupyterLab
- `cooltools >= 0.5.4` (the strictest version check used in the notebooks)

Python dependencies used across the repository include:

```text
bioframe cooler coolpuppy cooltools cytoolz h5py matplotlib multiprocess
networkx numpy packaging pandas pysam rpy2 scanpy scikit-learn scipy seaborn
statsmodels tqdm adpbulk
```

The R analyses use packages from CRAN and Bioconductor, principally:

```text
AnnotationHub ChIPpeakAnno ChIPseeker ComplexHeatmap DESeq2 DOSE
EnhancedVolcano GENOVA GeneOverlap GenomicFeatures GenomicRanges HilbertCurve
IRanges ReactomePA TxDb.Mmusculus.UCSC.mm10.ensGene babelgene biomaRt circlize
clusterProfiler cowplot dplyr enrichR enrichplot ensembldb ggVennDiagram
ggplot2 ggprism ggpubr ggrepel ggsci ggtree gplots gprofiler2 limma
org.Mm.eg.db patchwork pheatmap purrr reshape2 reticulate rstatix scales
stringr tidyr tidyverse tximport wesanderson
```

## Citation

If you use this code, please cite:

Kashuk, E., Smirnov, D., Eremenko, E., Tsitrina, A., Kriukov, D., Golova, A.,
Kaluski-Kopatch, S., Khrameeva, E., & Toiber, D. (2026). *The Divergent 3D
Genome Landscapes of Aging and Neurodegenerative mouse models*. Genome Research. [https://doi.org/10.1101/gr.281732.125](https://doi.org/10.1101/gr.281732.125)

```bibtex
@article{kashuk2025divergent,
  author  = {Kashuk, Ekaterina and Smirnov, Dmitrii and Eremenko, Ekaterina and
             Tsitrina, Alexandra and Kriukov, Dmitrii and Golova, Anastasia and
             Kaluski-Kopatch, Shai and Khrameeva, Ekaterina and Toiber, Debra},
  title   = {The Divergent 3D Genome Landscapes of Aging and Neurodegenerative
             mouse models},
  journal = {Genome Research},
  year    = {2026},
  doi     = {10.1101/gr.281732.125},
}
```
