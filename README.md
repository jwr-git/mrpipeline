# Combined MR-Coloc Pipeline 

This package combines Mendelian randomisation (MR) and colocalisation (coloc) with the aim of providing causal evidence between exposure-outcome pairs. The pipeline favours molecular data as exposure, for example, gene expression or protein abundance quantitative trait loci (QTLs), but can be used for non-molecular exposures as well. The pipeline will also annotate and combine results with drug target-related evidence from the literature, and has functionality to output shareable reports between investigators.

The Automated Mendelian randomisation (MR) Pipeline is designed for people with little to no experience of MR to conduct their own MR analyses and receive quick results in an easily interpretable and shareable format. The Pipeline will ingest data from either the [OpenGWAS Database](https://gwas.mrcieu.ac.uk/) or from local GWAS stored in the [.vcf file format](https://github.com/MRCIEU/gwas-vcf-specification) and conduct a fairly standard selection of MR analyses, including sensitivity analyses such as MR-Egger, heterogeneity and colocalisation analyses. The Pipeline outputs results from these analyses in .html format files, which are browsable and interactive and may be shared without the need to send a multitude of images, data, etc.

# How to Install

Within R, the package can be installed for this GitHub like so:

```remotes::install_github("https://github.com/jwr-git/mrpipeline")```

# Documentation

See here: _coming soon_

# Docker

A docker image will also be available soon.
