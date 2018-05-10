Seminars
================
Abdullah Farouk
2018-01-17

### Reproducing graph (Seminar 2b)

``` r
library(tidyverse)
```

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr
    ## Loading tidyverse: dplyr

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
library(ggplot2)

ggplot(data = mpg) + 
geom_point(mapping = aes(x = displ, y = hwy, color = drv, size = class)) + scale_color_brewer(palette="Dark2")
```

    ## Warning: Using size for a discrete variable is not advised.

![](Seminar_2b_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-1-1.png)
