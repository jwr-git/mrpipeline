# Automated MR Pipeline 

The Automated Mendelian randomisation (MR) Pipeline is designed for people with little to no experience of MR to conduct their own MR analyses and receive quick results in an easily interpretable and shareable format. The Pipeline will ingest data from either the [OpenGWAS Database](https://gwas.mrcieu.ac.uk/) or from local GWAS stored in the [.vcf file format](https://github.com/MRCIEU/gwas-vcf-specification) and conduct a fairly standard selection of MR analyses, including sensitivity analyses such as MR-Egger, heterogeneity and colocalisation analyses. The Pipeline outputs results from these analyses in .html format files, which are browsable and interactive and may be shared without the need to send a multitude of images, data, etc.

# How to Install

Within R, the package can be installed for this GitHub like so:

```remotes::install_github("https://github.com/jwr-git/mrpipeline")```

# How to Use

The Pipeline has one user-facing function, `mr_pipeline()`, which requires at least two arguments to function. These arguments can be either a single OpenGWAS DB ID or a path to a .vcf file (.gz versions are also acceptable), or a vector of these (a mixture of IDs and .vcf files is allowed). The first ID or vector will be treated the exposures for the MR analysis, and the second ID or vector as the outcome. For example:

```mrpipeline::mr_pipeline_(c("../../IEU-a-2.vcf.gz", "ieu-a-1"), c("ieu-b-85", "ieu-a-1082"))```

Will initiate the Pipeline with a local .vcf file (IEU-a-2.vcf.gz) and [ieu-a-1](https://gwas.mrcieu.ac.uk/datasets/ieu-a-1/) (Adiponectin) as exposures, and [ieu-b-85](https://gwas.mrcieu.ac.uk/datasets/ieu-b-85/) (Prostate cancer) and [ieu-a-1082](https://gwas.mrcieu.ac.uk/datasets/ieu-a-1082/) (Thyroid cancer) as outcomes.

# Output

Coming soon.

# Analysis layout

Coming soon.

