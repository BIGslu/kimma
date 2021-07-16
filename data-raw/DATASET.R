## code to prepare `DATASET` dataset goes here

library(edgeR)
library(limma)

#Data
load("/Users/kim/Documents/_Altman/MS/P259_MS/P259_pDC_public/data_clean/P259_pDC_clean.RData")
count <- read.csv("/Users/kim/Documents/_Altman/MS/P259_MS/P259_pDC_public/data_clean/P259_pDC_counts.csv")

#Select 5 libraries
targets <- dat.pDC.voom$targets[c(1,2,5,6,9,10),c(5,6,17)]
E <- count[,c("geneName", targets$libID)]

#de-identify
targets$libID <- paste("lib",1:6,sep="")
rownames(targets) <- paste("lib",1:6,sep="")

targets$virus <- rep(c("A","B"),3)

colnames(E) <- c("geneName", paste("lib",1:6,sep=""))

#Subset genes
set.seed(546)
genes <- dat.pDC.voom$genes[sample(1:nrow(dat.pDC.voom$E), 1000, replace=FALSE),]

E <- E[E$geneName %in% genes$geneName, ]

#Move rownames
rownames(E) <- E$geneName
E <- E[,-1]

#Format and normalize
dat.example <- edgeR::DGEList(counts=E, samples=targets, genes=genes)

dat.voom.example <- calcNormFactors(dat.example)
dat.voom.example <- limma::voom(dat.voom.example)

#Add to package
usethis::use_data(dat.example, overwrite = TRUE)
usethis::use_data(dat.voom.example, overwrite = TRUE)
