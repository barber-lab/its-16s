# A Global Survey of Host, Aquatic, and Soil Microbiomes Reveals Ecological Properties Shared between Bacterial and Fungal Generalists

This repository provides the code required to reproduce the analyses of the following publication:

Daniel Loos, Ailton Pereira da Costa Filho, Amelia E. Barber, Gianni Panagiotou: A Global Survey of Host, Aquatic, and Soil Microbiomes Reveals Ecological Properties Shared between Bacterial and Fungal Generalists (2022). In review

Correspondence: gianni.panagiotou@leibniz-hki.de

## Abstract

Microbiome engineering is a fast-evolving area with relevance for human health, agriculture, and climate management solutions. Despite significant efforts in engineering microbiomes to repair dysbiotic communities, new microbes often fail to establish and/or alter ecosystem function. To identify bacterial and fungal genera with the desired ability to manipulate microbial communities, we retrieved paired 16S and ITS rRNA amplicon sequence data from 1,580 host, soil, and aquatic samples and explored the ecological patterns of the 2,977 bacteria and 1,740 fungal genera detected across all samples. Through this large-scale analysis, we revealed that a small number of bacterial and fungal generalists with high prevalence across all environments positively contribute to the taxonomic diversity of their respective kingdom and explain a large percentage of the variation in the cross-kingdom community structure. We also observed that bacterial and fungal generalists have a significantly higher abundance compared to specialists, or genera whose prevalence was strongly associated with a single habitat - possibly due to their ability to avoid competitive interactions and instead elicit positive ones with other highly prevalent genera. These findings can streamline existing strategies to identify bacterial and fungal inoculants with a higher probability to establish in recipient ecosystems and confer noticeable changes in their structure and function.

## Get Started

Source code is stored in the directory [`src`](src).
This workflow starts with a [nextflow](https://www.nextflow.io/ pipeline to retrieve taxon abundance profiles in directory [`src/its-16s-nf`](src/its-16s-nf).
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
