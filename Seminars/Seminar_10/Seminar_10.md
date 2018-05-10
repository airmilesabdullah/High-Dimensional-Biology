Seminar 10
================
Abdullah Farouk
2018-03-25

Load libraries

``` r
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
    ## ✔ tibble  1.4.2     ✔ dplyr   0.7.4
    ## ✔ tidyr   0.8.0     ✔ stringr 1.3.0
    ## ✔ readr   1.1.1     ✔ forcats 0.3.0

    ## ── Conflicts ───────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(dplyr)
library(ggplot2)
library(ermineR)
library(reshape2)
```

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

ErmineJ

``` r
#Get Data
(geneList <- read_csv("ranked_gene_list.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   logFC = col_double()
    ## )

    ## # A tibble: 17,421 x 2
    ##    gene       logFC
    ##    <chr>      <dbl>
    ##  1 A1BG     -4.57  
    ##  2 A1BG-AS1  0.190 
    ##  3 A1CF     -3.15  
    ##  4 A2M      -1.08  
    ##  5 A2M-AS1   0.760 
    ##  6 A4GALT   -0.131 
    ##  7 A4GNT     0.169 
    ##  8 AA06      0.0212
    ##  9 AAAS     -0.192 
    ## 10 AACS      0.345 
    ## # ... with 17,411 more rows

``` r
(mfScores <- read_csv("gene_multifunctionality_scores.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   MF.score = col_double()
    ## )

    ## # A tibble: 18,402 x 2
    ##    gene     MF.score
    ##    <chr>       <dbl>
    ##  1 MEF2C   0.000762 
    ##  2 MEF2A   0.000729 
    ##  3 OPA1    0.0000886
    ##  4 SLC25A4 0.000116 
    ##  5 TYMP    0.000159 
    ##  6 LONP1   0.000357 
    ##  7 MPV17   0.0000162
    ##  8 MGME1   0.0000578
    ##  9 PIF1    0.000347 
    ## 10 MRPL17  0.0000136
    ## # ... with 18,392 more rows

``` r
# download the latest GO.xml file if it doesn't already exist
if (!file.exists("GO.xml")) { goToday("GO.xml") }

ermineInputGeneScores <- geneList %>% 
  mutate(absolute_logFC = abs(logFC)) %>% 
  select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

head(ermineInputGeneScores) # print the first few rows
```

    ##                     absolute_logFC
    ## SAA2-SAA4|SAA1|SAA2       7.790951
    ## C9                        6.691756
    ## AVPR1A                    6.649647
    ## CFHR3|CFHR4               5.994606
    ## F9                        5.942541
    ## CRP                       5.797936

``` r
enrichmentResult <- precRecall(scores = ermineInputGeneScores, 
                               scoreColumn = 1, # column 1 is the scores 
                               bigIsBetter = TRUE, # larger logFC should be ranked higher
                               annotation = "Generic_human", # ask ermineJ to use the Generic_human annotation file (will automatically download)
                               aspects = "B", # look at only biological processes 
                               iterations = 10000, # 10K sampling iterations so that results are stable
                               geneSetDescription = "GO.xml") # use the GO XML file in current directory

#Now, let's take a look at output$result, which's probably what you care about the most
enrichmentResult$results %>% arrange(MFPvalue)
```

    ## # A tibble: 3,144 x 12
    ##    Name         ID    NumProbes NumGenes RawScore     Pval CorrectedPvalue
    ##    <chr>        <chr>     <int>    <int>    <dbl>    <dbl>           <dbl>
    ##  1 lymphocyte … GO:0…       103      103   0.0430 1.00e-12  0.00000000302 
    ##  2 humoral imm… GO:0…        35       35   0.0735 1.00e-12  0.00000000151 
    ##  3 adaptive im… GO:0…       117      117   0.0407 1.00e-12  0.00000000101 
    ##  4 acute infla… GO:0…        58       58   0.0392 1.00e-12  0.000000000755
    ##  5 platelet de… GO:0…       122      122   0.0340 1.00e-12  0.000000000604
    ##  6 regulation … GO:0…        65       65   0.0371 1.00e-12  0.000000000503
    ##  7 regulation … GO:0…        44       44   0.0499 1.00e-12  0.000000000431
    ##  8 complement … GO:0…        46       46   0.0742 1.00e-12  0.000000000302
    ##  9 complement … GO:0…        32       32   0.0769 1.00e-12  0.000000000274
    ## 10 humoral imm… GO:0…       139      139   0.0577 1.00e-12  0.000000000252
    ## # ... with 3,134 more rows, and 5 more variables: MFPvalue <dbl>,
    ## #   CorrectedMFPvalue <dbl>, Multifunctionality <dbl>, `Same as` <chr>,
    ## #   GeneMembers <chr>

Show the TA your list of top 10 GO terms that experience the largest adjustments to get checked off. And good day!

``` r
enrichmentResult$results %>% 
  select(Name, Pval, MFPvalue) %>% 
  mutate(neg_log_pvalue = -log10(Pval),
         neg_log_mfpvalue = -log10(MFPvalue)) %>% 
  mutate(log_pvalue_change = neg_log_mfpvalue - neg_log_pvalue) %>% 
  arrange(desc(abs(log_pvalue_change))) %>% 
  head(10) %>% 
  knitr::kable()
```

| Name                                                |  Pval|  MFPvalue|  neg\_log\_pvalue|  neg\_log\_mfpvalue|  log\_pvalue\_change|
|:----------------------------------------------------|-----:|---------:|-----------------:|-------------------:|--------------------:|
| regulation of mononuclear cell proliferation        |     0|    0.0651|                12|            1.186419|           -10.813581|
| regulation of lymphocyte proliferation              |     0|    0.0651|                12|            1.186419|           -10.813581|
| positive regulation of T cell activation            |     0|    0.0651|                12|            1.186419|           -10.813581|
| regulation of leukocyte proliferation               |     0|    0.0651|                12|            1.186419|           -10.813581|
| positive regulation of leukocyte cell-cell adhesion |     0|    0.0316|                12|            1.500313|           -10.499687|
| organic hydroxy compound biosynthetic process       |     0|    0.0094|                12|            2.026872|            -9.973128|
| regulation of response to wounding                  |     0|    0.0093|                12|            2.031517|            -9.968483|
| monocarboxylic acid biosynthetic process            |     0|    0.0072|                12|            2.142667|            -9.857333|
| regulation of wound healing                         |     0|    0.0071|                12|            2.148742|            -9.851258|
| cellular amino acid catabolic process               |     0|    0.0068|                12|            2.167491|            -9.832509|
