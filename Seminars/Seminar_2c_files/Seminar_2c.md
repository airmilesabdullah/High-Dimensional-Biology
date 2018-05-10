Seminar 2c
================
Abdullah Farouk
2018-01-17

Central Limit Theorem
---------------------

``` r
#Load packages
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyverse)
```

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
set.seed(1) #Set seed for reproducibility

#Set sample size and number of samples to sample over
sampleSize <- 5
numSamples <- 1000

#Draw samples from a chi square distribution
degreeFreedom <- 1
randomChiSqValues <- rchisq(n = numSamples * sampleSize, df = degreeFreedom)

#Plot the distribution of sample points in randomChiSqValues
tibble(x = randomChiSqValues) %>% 
  ggplot() + 
  geom_density(aes(x = x), color = "red")
```

![](Seminar_2c_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-1-1.png)

``` r
# organize the random values into 1000 sample rows of size n = 5 columns
samples <- matrix(randomChiSqValues, nrow = numSamples, ncol = sampleSize)
sampleMeans <- rowMeans(samples) # Calculate sample means 

#Plot of the sample means
tibble(x = sampleMeans) %>% 
  ggplot() + 
  geom_line(aes(x = x), stat = "density", color = "green") +
  geom_point(aes(x = x, y = 0), color = "yellow", shape = 9, size = 3)
```

![](Seminar_2c_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-1-2.png)

``` r
?points
```
