
``` r
library(GEOquery)
library(devtools)
library(Biobase)
# source('http://bioconductor.org/biocLite.R')
# biocLite('GEOquery')
```

``` r
load("/Users/abdullah/Downloads/eQTL_melanoma.Rdata")

colnames(exprsdat) = exprs.meta[exprs.meta$geo_accession == colnames(exprsdat), c("cell_line")]
colnames(genodat.clean) = geno.meta.clean[geno.meta.clean$geo_accession == colnames(genodat.clean), c("cell_line")]

# For our expression data
# Are the number of rows without any NA value, the same as total number of rows?  
dim(exprsdat) == dim(exprsdat[complete.cases(exprsdat),])
```

    ## [1] TRUE TRUE

``` r
# For our genotype data
# Are the number of rows without any NA value, the same as total number of rows?  
dim(genodat.clean) == dim(genodat.clean[complete.cases(genodat.clean),])
```

    ## [1] TRUE TRUE

``` r
genodat.numeric = genodat.clean
genodat.numeric[genodat.numeric == "NC"] = NA
genodat.numeric = data.frame(sapply(genodat.numeric, function(x) as.numeric(x) - 1))
rownames(genodat.numeric) = rownames(genodat.clean)


exprsdat_log = log2(exprsdat)
probe_means = rowMeans(exprsdat_log)
probe_vars = apply(exprsdat_log, 1 , sd)
#plot(probe_vars~probe_means, xlab="Mean", ylab="Variance")
```

``` r
random_gene_ix = "208644_at"
random_snp_ix = "rs3219090"

#LM

exprs_random = as.numeric(exprsdat_log[random_gene_ix,])
snp_random = as.numeric(genodat.numeric[random_snp_ix,])
lm_random = lm(exprs_random ~ snp_random)
summary(lm_random)
```

    ## 
    ## Call:
    ## lm(formula = exprs_random ~ snp_random)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.27635 -0.04492 -0.00268  0.05974  0.36638 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.115520   0.030723 101.407   <2e-16 ***
    ## snp_random  -0.001252   0.019846  -0.063     0.95    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.101 on 56 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  7.101e-05,  Adjusted R-squared:  -0.01778 
    ## F-statistic: 0.003977 on 1 and 56 DF,  p-value: 0.9499

The p value of the SNP covariate is 0.95. This suggests that effect is statistically insignifant in explaining the variation observed in gene expression values. This implies no statistical association of the SNP "208644\_at" with gene "rs3219090".

Fill in the blanks using 3 of the 5 options listed below: An eQTL analysis usually involves treating SNPs as covariates and fitting a linear model, where the outcome is gene expression values.
