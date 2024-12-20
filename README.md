# A global survey of host, aquatic, and soil microbiomes reveals shared abundance and genomic features between bacterial and fungal generalists

This repository provides the code required to reproduce the analyses of the following publication:

Daniel Loos, Ailton Pereira da Costa Filho, Bas E. Dutilh, Amelia E. Barber, Gianni Panagiotou. **A global survey of host, aquatic, and soil microbiomes reveals shared abundance and genomic features between bacterial and fungal generalists**. *Cell Reports*. 2024. [https://doi.org/10.1101/2022.11.15.515575](https://doi.org/10.1016/j.celrep.2024.114046).

Correspondence: amelia.barber@uni-jena.de or gianni.panagiotou@leibniz-hki.de

## Abstract

Environmental change, coupled with alteration in human lifestyles, is profoundly impacting the microbial communities critical to the health of the Earth and its inhabitants. To identify bacteria and fungi that are resistant and susceptible to habitat change, we analyze thousands of genera detected in 1,580 host, soil, and aquatic samples. This large-scale analysis identifies 48 bacterial and 4 fungal genera that are abundant across the three biomes, demonstrating fitness in diverse environmental conditions. Samples containing these generalists have significantly higher alpha diversity. These generalists play a significant role in shaping cross-kingdom community structure, boasting larger genomes with more secondary metabolism and antimicrobial resistance genes. Conversely, 30 bacterial and 19 fungal genera are only found in a single habitat, suggesting a limited ability to adapt to different and changing environments. These findings contribute to our understanding of microbial niche breadth and its consequences for global biodiversity loss.

## Get Started

Source code is stored in the directory [`src`](src).
This workflow starts with a [nextflow](https://www.nextflow.io/) pipeline to retrieve taxon abundance profiles in directory [`src/its-16s-nf`](src/its-16s-nf).
Then, a [drake](https://docs.ropensci.org/drake/) pipeline is used for subsequent statistical analysis and visualization in directory [`src/plans`](src/plans).
The main entry point for this analysis is the [`Makefile`](Makefile) providing shell commands `make build` for building the [docker image](Dockerfile) and `make analyse` to execute the pipeline.
Briefly, R packages were installed using script [`src/install/install.R`](src/install/install.R).
Then script [`main.R`](main.R) is executed that bootstraps the global R environment using script [`src/defaults.R`](src/defaults.R) and executes the individual drake targets defined in directory [`src/plans`](src/plans).
The methods section of the paper provides further details.

### Notable targets of the pipeline

| Target name                            | Description                                                                                                  |
|----------------------------------------|--------------------------------------------------------------------------------------------------------------|
| `bioprojects`                          | Cohorts used in this meta study                                                                              |
| `samples`                              | Table of biological samples with meta data                                                                   |
| `abundances`                           | Abundance profiles of all samples                                                                            |
| `alphadiv`                             | Alpha diversities per sample subset                                                                          |
| `stepwise_constrained_ordinations`     | dbRDA to ordinate the abundance profile of one kingdom using the best explaining taxa from the other kingdom |
| `selected_generalists_specialists`     | Generalists and specialists of the meta study                                                                |
| `generalists_common_coabundance_graph` | Co-abundances common in most environments                                                                    |
| `fig6_heatmap`                         | Abundance heatmap of the paper                                                                               |

## Supplementary analysis

The get started section describe analysis done with R Drake. Here we describe extra analysis performed with R Targets and Snakemake. They are available under [`src/sup`](src/sup). They include:

| Task                                      | Description                                                                               |
|-------------------------------------------|-------------------------------------------------------------------------------------------|
|Primer detection                           | 16S and ITS primers for bioprojects without publication were discovered based on alignment|
|Generalist and specialists BGCs            | Fungal and Bacterial BGCs were identified with antiSMASH                                  |
|Bacterial generalists and specialists AMR  | Antimicrobial Resistance Genes identification                                             |
| Decontamination                           | Perform decontamination of low biomass samples for 16S and ITS                            |

