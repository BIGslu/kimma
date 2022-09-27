## code to prepare `DATASET` dataset goes here
library(tidyverse)
library(edgeR)
library(limma)

#Data
load("data-raw/P259_pDC_clean.RData")
count <- read.csv("data-raw/P259_pDC_counts.csv")

#Select 12 libraries
targets <- dat.pDC.voom$targets %>%
  filter(virus.detail != "newHRV" & IL5 != "EOS.supp" & donorID != "donor4" &
           !grepl("AT", donorID)) %>%
  select(group, libID, donorID, median_cv_coverage, virus, asthma) %>%
  mutate(virus = fct_relevel(factor(virus), ref="none")) %>%
  mutate(asthma = fct_relevel(factor(asthma), ref="healthy"))

E <- count[,c("geneName", targets$libID)]

count.filter <- count %>%
  select(geneName, all_of(targets$libID))

#de-identify
targets$libID <- paste("lib",1:12,sep="")
rownames(targets) <- paste("lib",1:12,sep="")
colnames(E) <- c("geneName", paste("lib",1:12,sep=""))
colnames(count.filter) <- c("geneName", paste("lib",1:12,sep=""))

set.seed(546)
targets <- targets %>%
  mutate(donorID = recode(donorID, "AC3"="donor4","AC4"="donor5","AC5"="donor6"),
         batch = ifelse(donorID %in% c("donor1","donor2","donor3"), "1","2"))

#Subset genes
set.seed(546)
genes <- dat.pDC.voom$genes[sample(1:nrow(dat.pDC.voom$E), 1000, replace=FALSE),]

E <- E[E$geneName %in% genes$geneName, ]

count.filter <- count.filter %>%
  filter(geneName %in% genes$geneName)

#Move rownames
rownames(E) <- E$geneName
E <- E[,-1]

#Format and normalize
dat.example <- edgeR::DGEList(counts=E, samples=targets, genes=genes)

dat.voom.example <- edgeR::calcNormFactors(dat.example)
dat.voom.example <- limma::voom(dat.voom.example)

#Kinship data

kin.example <- data.frame(rowname = c("donor1","donor2","donor3","donor4","donor5","donor6"),
                          donor1 = c(1, 0.5, 0.1, 0.2, 0.15, 0.1),
                          donor2 = c(0.5, 1, 0.1, 0.11, 0.23, 0.09),
                          donor3 = c(0.1, 0.1, 1, 0.1, 0.2, 0.15),
                          donor4 = c(0.2, 0.11, 0.1, 1, 0.11, 0.13),
                          donor5 = c(0.15, 0.23, 0.2, 0.11, 1, 0.17),
                          donor6 = c(0.1, 0.09, 0.15, 0.13, 0.17, 1)) %>%
  tibble::column_to_rownames() %>% as.matrix()

#Rename donorID to default
dat.example$samples <- dat.example$samples %>%
  rename(ptID=donorID)
dat.voom.example$targets <- dat.voom.example$targets %>%
  rename(ptID=donorID)

#Rename
example.count <- count.filter
example.dat <- dat.example
example.voom <- dat.voom.example
example.kin <- kin.example

#Add to package
usethis::use_data(example.count, overwrite = TRUE)
usethis::use_data(example.dat, overwrite = TRUE)
usethis::use_data(example.voom, overwrite = TRUE)
usethis::use_data(example.kin, overwrite = TRUE)

######
# SEAsnake outputs for tutorials
library(tidyverse)
fcounts <- read_tsv("data-raw/combined_feature_counts.tsv")
flagstat <- read_tsv("data-raw/combined_flagstat.tsv")
picard <- read_tsv("data-raw/combined_picard.tsv")
patient <- read_csv("data-raw/patient_meta.csv")
sample <- read_csv("data-raw/sample_meta.csv")

example.seasnake <- list()
example.seasnake[["fcounts"]] <- fcounts
example.seasnake[["flagstat"]] <- flagstat
example.seasnake[["picard"]] <- picard
example.seasnake[["patient"]] <- patient
example.seasnake[["sample"]] <- sample

usethis::use_data(example.seasnake, overwrite = TRUE)
