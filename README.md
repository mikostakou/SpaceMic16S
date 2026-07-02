# SpaceMic16S

### R code for the manuscript

> **Kostakou et al., 2026** *"The illusion of diversity: sampling design drives conflicting estimates of soil bacterial richness"*
> 

This repository contains the R code underlying the analyses, statistics, and figures presented in the manuscript. The scripts allow full reproduction of the main data-processing steps, statistical models, and figures. 


## Project overview

This study evaluates how soil **sampling design** shapes estimates of bacterial diversity across forest and grassland ecosystems in the German Biodiversity Exploratories.

The central question is how estimates of soil bacterial **richness**, **sample coverage**, and **spatial scaling** change when diversity is inferred from *individual soil cores* versus *composite (physically pooled) samples* — and how these patterns shift under **sequencing-depth standardization**. A recurring result is that seemingly minor design choices can reverse ecological conclusions (e.g. which ecosystem appears more diverse), hence the "illusion of diversity."



## Sampling details

| | |
|---|---|
| **Study system** | Biodiversity Exploratories, Germany |
| **Ecosystems** | Forest and grassland soils |
| **Plots** | 57 plots across the three Exploratories |
| **Sampling design** | 14 individual soil cores per plot |
| **Composite samples** | Physical homogenization of soil cores at the plot level |
| **Target marker** | 16S rRNA gene |
| **Sequencing** | Paired-end amplicon sequencing |



## Data processing and analysis

The repository includes scripts for:

- processing paired-end 16S rRNA amplicon sequencing data with **DADA2**;
- assigning taxonomy using the **SILVA** reference database;
- constructing **phyloseq** objects;
- rarefying sequencing depth;
- generating plot-wise spatial subsampling scenarios;
- estimating richness and sample coverage across increasing sampling effort;
- comparing individual-based and composite-sample diversity estimates;
- evaluating richness scaling with spatial extent;
- fitting the generalized linear mixed models (GLMMs) associated with each analysis;
- producing the main and supplementary figures.



## Data availability

 **Raw sequence data:** deposited at NCBI under BioProject accession number `[PRJNA1242586]`.
 **Metadata:** available through the Biodiversity Exploratories Information System (BExIS)(https://doi.org/10.17616/R32P9Q):
 - sample metadata - BExIS ID 32361, https://doi.org/10.71615/bexis.32361; 
 - taxonomic assignments - BExIS ID 32394, https://doi.org/10.71615/bexis.32394; 
 - sequence variant abundance data - BExIS ID 32420, https://doi.org/10.71615/bexis.32420.



## Contact

Maria Kostakou — Helmholtz Centre for Environmental Research (UFZ), Leipzig
mikostakou@gmail.com

## License

This repository is licensed under the MIT License. See the LICENSE file for details.
