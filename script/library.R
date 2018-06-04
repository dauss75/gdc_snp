# packages <- c("data.table", "dplyr", "caret")
# if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#   install.packages(setdiff(packages, rownames(installed.packages())),repos = "http://cran.us.r-project.org")
# }
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))

# bio_packages<-c("maftools","GenomicFeatures","TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db,annotate")
# if (length(setdiff(bio_packages, rownames(installed.packages()))) > 0) {
#   source("https://bioconductor.org/biocLite.R")
#   biocLite(setdiff(bio_packages, rownames(installed.packages())))
# }
suppressPackageStartupMessages(library(maftools))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(annotate))
# rm(packages, bio_packages)
