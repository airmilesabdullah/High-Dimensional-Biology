Validation\_data\_set
================
Abdullah Farouk
2018-04-01

Load libraries
==============

``` r
library(knitr)
library(GEOquery)
library(biomaRt)
library(tidyverse)
library(data.table)
library(wateRmelon)
library(FDb.InfiniumMethylation.hg19)
```

``` r
GSE73375 <- read.delim("/Users/abdullah/Desktop/ProjectData/GSE73375_signal.txt",fill = TRUE, header = TRUE)
```

Read meta data and cross reactive data in
=========================================

``` r
sample_meta_data <- read.csv("/Users/abdullah/Desktop/ProjectData/GPL13534_HumanMethylation450_15017482_v.1.1.csv", header = TRUE)

colnames(sample_meta_data)[1] <- "TargetID"

#loading hybridization data (you can download from the link that Amy send in the messanger)
Crossreactive_Magda <- read.csv("/Users/abdullah/Desktop/ProjectData/450KMagdaAnnotation.csv",fill = TRUE, header = TRUE)
```

Eliminate probes with large detection p values
==============================================

``` r
#selecting columns only with p-value
pval_test <- dplyr::select(GSE73375, ends_with(".Pval"))

rownames(pval_test) <- GSE73375[,1]

#find the largest p value accross all samples and store it in the orgininal dataset
GSE73375_pval <- GSE73375

GSE73375_pval$max.Pval <- apply(pval_test,1,max)

#filter out the dataset with the largest p-value > 0.05
GSE73375_pfiltered <- filter(GSE73375_pval, max.Pval < 0.05)
```

Calulate Beta Values
====================

``` r
GSE73375_Bval <- GSE73375_pfiltered

#Initialize for loop
x_1 <- 1     #Column to iterate over
y_1 <- 111   #New column to add to GSE446677_pfiltered dataset
j <- 1       #Number of samples in dataset

for(j in 1:36){
  
  d <- GSE73375_pfiltered[ ,x_1+1] + GSE73375_pfiltered[ ,x_1+2]
  
  GSE73375_Bval[ ,y_1] <- GSE73375_pfiltered[ ,x_1+2]/d
  
  x_1 = x_1 + 3
  
  j = j + 1
  
  y_1 = y_1 + 1
  
}

#Keep columns with beta values
GSE73375_Bval <- GSE73375_Bval[ ,111:146]
```

Add column names
================

``` r
#Add column names
col_names <- c(paste0('P', 1:6), paste0('Q', 1:6), paste0('P', 9:12), paste0('Q', 7:12), paste0('P', 19:20), paste0('Q', 19:20), paste0('P', 13:18), paste0('Q', 13:18))

col_names <- col_names[-which(col_names == 'P3')] #Remove P3
col_names <- col_names[-which(col_names == 'Q3')] #Remove Q3

#For Beta values object
colnames(GSE73375_Bval) <- col_names

rownames(GSE73375_Bval) <- GSE73375_pfiltered[ ,1]
```

Filter out data from X and Y chromosome and QC SNP probe
========================================================

``` r
samples_filt <- sample_meta_data[!sample_meta_data$CHR%in%c("X", "Y"), ]
 
samples_filt <- samples_filt[grep("^rs", samples_filt$TargetID, invert = TRUE), ]
```

Filter out cross hybridization probes
=====================================

``` r
Crossreactive_MagdaDatXY <- Crossreactive_Magda[grep("YES",Crossreactive_Magda$XY_Hits), ]

Crossreactive_MagdaDatAuto <- Crossreactive_Magda[grep("YES",Crossreactive_Magda$Autosomal_Hits), ]

Cross_MagdaXY <-samples_filt[which((samples_filt$TargetID) %in% (Crossreactive_MagdaDatXY$IlmnID)), ]

Cross_MagdaAuto <-samples_filt[which((samples_filt$TargetID) %in% (Crossreactive_MagdaDatAuto$IlmnID)), ]

samples_filt <- samples_filt[!(samples_filt$TargetID) %in% (Crossreactive_MagdaDatXY$IlmnID), ]

samples_filt <- samples_filt[!(samples_filt$TargetID) %in% (Crossreactive_MagdaDatAuto$IlmnID), ]
```

Combining filtered annotation data with B-value data
====================================================

``` r
Filtered_GSE73375_Bval <- GSE73375_Bval[which(rownames(GSE73375_Bval) %in% (samples_filt$TargetID)), ]
```

Create Beta matrix
==================

``` r
#Create matrices for beta and m values
beta_matrix_test <- as.matrix(Filtered_GSE73375_Bval)
```

Quantile Normalization
======================

``` r
beta_norm_test_set <- betaqn(beta_matrix_test)
m_norm_test_set <- beta2m(beta_norm_test_set)
```

Load the data and clean the meta data
-------------------------------------

``` r
GSE73375_1 <- getGEO('GSE73375')
```

    ## Found 1 file(s)

    ## GSE73375_series_matrix.txt.gz

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   ID_REF = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## File stored at:

    ## /var/folders/h1/twl48y653xdd70k7gfzmtpm00000gn/T//RtmpXNA2zd/GPL13534.soft

``` r
Meta_data_73375 <- pData(phenoData(GSE73375_1[[1]]))


Meta_data_73375 <- Meta_data_73375[,(colnames(Meta_data_73375) %in% c("geo_accession",
                                                  "characteristics_ch1", 'characteristics_ch1.2'))]


colnames(Meta_data_73375) <- c("Sample", "Disease_Class", "Gestational_Age")

# Meta_data_73375$Sample <- paste("X", gsub(" placenta", "", Meta_data_73375$Sample), sep="")
# Meta_data_73375$Sample <- gsub("Sample name: ", "", Meta_data_73375$Sample)

#Check if columns are the same, and in the same order
sum(colnames(beta_norm_test_set) == Meta_data_73375$Sample)
```

    ## [1] 0

``` r
#Re code disease class and gestational age
Meta_data_73375$Disease_Class <- recode_factor(Meta_data_73375$Disease_Class, `diagnosis: normotensive` = "normal", `diagnosis: preeclamptic` = "preeclampsia")

Meta_data_73375$Gestational_Age <- recode_factor(Meta_data_73375$Gestational_Age, `maternal age: 19` = 19, `maternal age: 20` = 20, `maternal age: 22` = 22, `maternal age: 24` = 24, `maternal age: 25` = 25, `maternal age: 26` = 26, `maternal age: 27` = 27, `maternal age: 28` = 28, `maternal age: 29` = 29,`maternal age: 30` = 30, `maternal age: 31` = 31, `maternal age: 32` = 32, `maternal age: 33` = 33, `maternal age: 34` = 34, `maternal age: 35` = 35, `maternal age: 36` = 36, `maternal age: 37` = 37, `maternal age: 38` = 38)

Meta_data_73375$Gestational_Age <- as.numeric(as.character(Meta_data_73375$Gestational_Age))
  
Meta_data_73375$Gestational_Age <- ifelse(Meta_data_73375$Gestational_Age < 37, "Preterm", "Term")

#Final cleaning of data
betaV_data_norm = as.data.frame(beta_norm_test_set)
rm(beta_norm) #Removes beta_norm
```

    ## Warning in rm(beta_norm): object 'beta_norm' not found

``` r
#mV_data = as.data.frame(m_matrix)
mV_data_norm = as.data.frame(m_norm_test_set)
rm(m_norm_test_set)
```

Clustering CpGs to CGIs
-----------------------

``` r
# Annotation data pulled from FDb database  
IM_data = features(FDb.InfiniumMethylation.hg19)
# Grab annotation data from hg19 genome  
IM_meta = metadata(FDb.InfiniumMethylation.hg19)
genome(IM_data) = IM_meta[which(IM_meta[,'name']=='Genome'),'value']
IM_data = sort(IM_data)

# Filter down to probes within the 450k assay 
probes_450k = as.data.frame(IM_data[IM_data$platform %in% c("BOTH","HM450"),])
probes_450k$Probe_ID = rownames(probes_450k)
#cginame_onlycg = probes_450k[probes_450k$probeType == "cg", ]
#Obtain 450k annotation data and connect each probe to a CGI
hm450 = getPlatform(platform='HM450', genome='hg19')
```

    ## Fetching coordinates for hg19...

``` r
probe_UCSC_name = getNearestGene(hm450)

#Filter beta values and M values to only include probes found within CGIs
betaV_CGI = betaV_data_norm[rownames(betaV_data_norm) %in% rownames(probe_UCSC_name),]
mV_CGI = mV_data_norm[rownames(mV_data_norm) %in% rownames(probe_UCSC_name),]
#Filter down CGI metadata to those probes found within our dataset
cginame = probe_UCSC_name[rownames(mV_CGI),]
cginame$cginame = cginame$nearestGeneSymbol
cginame$Probe_ID = rownames(cginame)

# Aggregation of the CpGs to the CGI level
betaV_CGI = aggregate(betaV_CGI, by = list(cginame$nearestGeneSymbol), mean, na.rm = T)
rownames(betaV_CGI) = betaV_CGI[, "Group.1"]
betaV_CGI = betaV_CGI[!(names(betaV_CGI) %in% c("Group.1"))]

mV_CGI = aggregate(mV_CGI, by = list(cginame$nearestGeneSymbol), mean, na.rm = T)
rownames(mV_CGI) = mV_CGI[, "Group.1"]
mV_CGI = mV_CGI[!(names(mV_CGI) %in% c("Group.1"))]
```

Keep CGI sites that are the same across the two datasets
========================================================

``` r
m_val_CGI <- readRDS('/Users/abdullah/Desktop/ProjectData/CGI_Data/CGI_norm_m_values.RDS')

#Check rownames
row_names <- rownames(m_val_CGI) #rownames of dataset we are working with
row_names_test_set <- rownames(mV_CGI) #rownames of test dataset
rn <- which(row_names_test_set %in% row_names)
mV_CGI_1 <- mV_CGI[rn, ]

#Check
which(rownames(mV_CGI_1) %in% rownames(mV_CGI)==FALSE)
```

    ## integer(0)

Save all the data
=================

``` r
saveRDS(betaV_CGI, file = 'norm_b_values_test_set')
saveRDS(mV_CGI_1, file = 'm_val_test_set_cgi_match')
```
