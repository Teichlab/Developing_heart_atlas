# Developing heart atlas
Data processing and analysis for the papers below:
- [High-resolution atlas of the developing human heart and the great vessels](https://www.biorxiv.org/content/10.1101/2024.04.27.591127v1) (Bayraktar et al)
- [Multiomic analysis reveals developmental dynamics of the human heart in health and disease](https://www.biorxiv.org/content/10.1101/2024.04.29.591736v1) (Cranley et al)<be>
<br>
Processed data of sc/snRNAseq and Visium data will be available for browsing and download via [heartcellatlas.org](https://www.heartcellatlas.org/foetal.html) at the time of publication.

## Contents
### General processing (Bayraktar et al & Cranley et al - Figure 1)
- RNA (including cell type annotations)
- ATAC
- Visium
- VisiumHD
- Xenium-5K
### Overview of cellular niches and TissueTypist (Cranley et al - Figure 2)
- Cell type enrichment in cellular niches
- TissueTypist downstream analysis
The TissueTypist package is [here](https://github.com/Teichlab/TissueTypist)
### The developing sinoatrial node (Cranley et al - Figure 3)
- SAN pacemaker cell gene signatures
- SAN niche with Visium HD
- Cell-cell interaction analysis
### Cardiomyocyte development and maturation (Cranley et al - Figure 4)
- Atrial cardiomyocyte maturation
- Ventricular cardiomyoyte development
- Ventricular transmural axis analyiss using [OrganAxis](https://github.com/nadavyayon/TissueTag)
### Trisomy 21 (Cranley et al - Figure 5)
- Differential abundance test (Milo) based on integrated latent space
- Differentially expressed gene ananlysis and GSEA using Milo cellular neighbourhoods
### Supplementary materials (Cranley et al)
- Epigenetic stability analysis
- Cardiac macrophage development
