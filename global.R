library("shiny", lib.loc="D:/R-3.1.0/library")
library("aroma.affymetrix", lib.loc="D:/R-3.1.0/library")
library("matrixStats", lib.loc="D:/R-3.1.0/library")
library("mosaic", lib.loc="D:/R-3.1.0/library")
library("dendextend", lib.loc="D:/R-3.1.0/library")
library("scales", lib.loc="D:/R-3.1.0/library")
library("colorspace", lib.loc="D:/R-3.1.0/library")
library("dplyr", lib.loc="D:/R-3.1.0/library")
library("reshape2", lib.loc="D:/R-3.1.0/library")
library("grid", lib.loc="D:/R-3.1.0/library")
library("gridExtra", lib.loc="D:/R-3.1.0/library")
library("marray", lib.loc="D:/R-3.1.0/library")
library("d3heatmap", lib.loc="D:/R-3.1.0/library")
library("knitr", lib.loc="D:/R-3.1.0/library")
library("plyr", lib.loc="D:/R-3.1.0/library")
library("rstudio", lib.loc="D:/R-3.1.0/library")
library("lattice", lib.loc="D:/R-3.1.0/library")
library("doBy", lib.loc="D:/R-3.1.0/library")
library("latticeExtra", lib.loc="D:/R-3.1.0/library")
library("gplots", lib.loc="D:/R-3.1.0/library")
library("ggplot2", lib.loc="D:/R-3.1.0/library")
library("tidyr", lib.loc="D:/R-3.1.0/library")
library("RColorBrewer", lib.loc="D:/R-3.1.0/library")
library("BSgenome", lib.loc="D:/R-3.1.0/library")
library("BSgenome.Hsapiens.UCSC.hg18", lib.loc="D:/R-3.1.0/library")
library("GenomicFeatures", lib.loc="D:/R-3.1.0/library")
library("GenomicRanges", lib.loc="D:/R-3.1.0/library")
library("Gviz", lib.loc="D:/R-3.1.0/library")
library("TxDb.Hsapiens.UCSC.hg18.knownGene", lib.loc="D:/R-3.1.0/library")
library("biomaRt", lib.loc="D:/R-3.1.0/library")
require(hgu133plus2.db)

## datasets loaded as lists
Probe_Seq_AffyIndices<- readRDS("D:/chromosome_mapping_HGU133plus2_BothStrands.rda")
gExprs<-readRDS("D:/gExprs.rda")
Gene1<-read.csv("D:/REACTOME_CD28_DEPENDENT_VAV1_PATHWAY_PHpos.csv")
Gene2<-read.csv("D:/REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES_PHpos.csv")

## DataList
GeneChoice <- list("NULL" = NULL,
                   "CD28_VAV1_Pathway" = Gene1,
                   "Second_Messengers" = Gene2)

ALL_groups <- read.csv("D:/HR_Other_ALL_Groups.csv") ##This has 5 patient groups (MLL, Normal, PHpos, Other, E2A_PBX)
Group_count<-as.data.frame(tally(~Norm_HR, data=ALL_groups))
Group_count$Var1<-as.character(Group_count$Var1)##Group counts
ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="may2009.archive.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)




