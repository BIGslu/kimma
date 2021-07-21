## code to prepare `DATASET` dataset goes here

library(edgeR)
library(limma)

#Data
load("data-raw/P259_pDC_clean.RData")
count <- read.csv("data-raw/P259_pDC_counts.csv")

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

dat.voom.example <- edgeR::calcNormFactors(dat.example)
dat.voom.example <- limma::voom(dat.voom.example)

#Kinship data

kin.example <- data.frame(rowname = c("donor1","donor2","donor3"),
                          donor1 = c(1, 0.5, 0.1),
                          donor2 = c(0.5, 1, 0.1),
                          donor3 = c(0.1, 0.1, 1)) %>%
  tibble::column_to_rownames() %>% as.matrix()

#Add to package
usethis::use_data(dat.example, overwrite = TRUE)
usethis::use_data(dat.voom.example, overwrite = TRUE)
usethis::use_data(kin.example, overwrite = TRUE)
