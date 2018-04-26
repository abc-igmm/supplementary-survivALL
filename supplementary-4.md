
Continuous Biomarker Assessment by Exhaustive Survival Analysis - Supplementary Code 4
======================================================================================

#### 2018-04-26

Dominic A. Pearce<sup>1</sup>, Ajit J. Nirmal<sup>2</sup>, Tom Freeman<sup>2</sup>, Andrew H. Sims<sup>1</sup>

<sup>1</sup>Applied Bioinformatics of Cancer, University of Edinburgh Cancer Research Centre, Institute of Genetics and Molecular Medicine, Edinburgh, UK <sup>2</sup>Systems Immunology Group, Division of Genetics and Genomics, The Roslin Institute and Royal (Dick) School of Veterinary Studies, University of Edinburgh, Easter Bush, Midlothian, EH25 9RG
\*<andrew.sims@ed.ac.uk>

 

Bootstrapping and hazard ratio thresholds
-----------------------------------------

 

### Libraries

``` r
library(survivALL)
library(Biobase)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(knitr)
library(viridis)
```

 

When performing a large number of statistical test, as is integral to the `survivALL` rationale, it is important to protect against false positive results using some form of multiple testing correction. For `survivALL` this is implemented as a bootstrapping exercise to determine robust thresholds of hazard ratio significance. In short we calculate, for each point-of-separation a upper and lower limit within which we expect to see hazard ratios occur by chance, and beyond which hazard ratios are unlikely (1 in 20) to have occurred by chance.

To achieve this, we randomly sample our survival data with replacement and calculate survival statistics for all points-of-separation, as we would using for a biomarker under investigation. By repeating this procedure 1,000s or 10,000s of times, we produce a distribution of *expected* hazard ratios, of which we use the mean and standard deviation to calculate our per-point-of-separation significance.

 

``` r
disc <- readRDS("discovery-eset.Rds") #see Supplementary 1 for more info on disc

#bootstrapping data should be created in the format of 1 repetition per column
bs_mtx <- matrix(nrow = ncol(disc), ncol = 10000)

system.time(
            for(i in 1:ncol(bs_mtx)){
                bs_mtx[, i] <- allHR(measure = sample(1:ncol(disc), 
                                                      replace = TRUE),
                                     srv = pData(disc),
                                     time = "t.dmfs",
                                     event = "e.dmfs")
            }
)
```

``` r
kable(bs_mtx[1:20, 1:5])
```

|          V1|          V2|         V3|          V4|          V5|
|-----------:|-----------:|----------:|-----------:|-----------:|
|  -4.2985760|          NA|         NA|          NA|          NA|
|  -1.2497716|  -0.5006925|         NA|          NA|  -1.0721970|
|  -1.5452649|   0.2061469|         NA|          NA|  -1.7503348|
|  -2.0365849|   0.4693223|         NA|          NA|  -1.1198149|
|  -1.5134552|   0.4693223|         NA|          NA|  -0.4904476|
|  -1.0492810|   0.5561135|         NA|   0.8959700|  -0.8227830|
|  -0.7008058|  -0.2029949|         NA|  -0.0048401|  -1.1709956|
|  -0.6418548|   0.0974521|         NA|   0.2846912|  -0.9207821|
|  -0.5619627|   0.3314171|         NA|   0.5103646|  -1.1635300|
|  -0.4622613|   0.5601122|  1.7477283|   0.7147431|  -1.0587921|
|  -0.3395972|   0.0230835|  1.9244026|   0.8508627|  -0.8125205|
|  -0.0282329|   0.1291360|  1.0490482|   1.0287279|  -0.9450526|
|  -0.4110332|   0.3221217|  1.1820121|   1.0954043|  -0.7827766|
|  -0.2180065|   0.4929464|  1.2882057|   0.5551933|  -0.7527642|
|  -0.2180065|   0.1171855|  1.4229195|   0.1901938|  -0.9003891|
|  -0.0875949|   0.1292967|  0.8540924|   0.2286172|  -0.8406196|
|   0.0861964|  -0.1945422|  0.4484826|   0.2751000|  -0.9819058|
|   0.1836266|  -0.3395630|  0.5947533|   0.0135470|  -0.9242119|
|   0.3567013|  -0.2270098|  0.7047795|   0.1290378|  -0.9829026|
|   0.0474136|  -0.1397093|  0.7955054|  -0.1085851|  -0.9119966|

 

In essence, this procedure creates our distribution of expected hazard ratios that can themselves be visualised.

 

``` r
bs_dfr <- melt(bs_mtx)

ggplot(bs_dfr, aes(x = Var1, y = value)) + 
        geom_hline(yintercept = 0, linetype = 3) +
        #geom_point(alpha = 0.0065) + 
        geom_hex(bins = 200) + 
        scale_fill_viridis() + 
        theme_pander() + 
        labs(title = "Figures SC4-1", 
             subtitle = "Hazard ratio distribution calculated from 10,000 random sample orderings")
```

<img src="supplementary-4_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

 

Having calculated our bootstrapped data we then simply hand the matrix to either the `survivALL()` or `plotALL()` functions to handle the subsequent thresholding and p-value calculations. It should be noted that thresholding up to 10,000x can be a long process requiring an investment of time and that, if preceded by a ancillary test of significance using the continuous variable as the classifier itself in the coxph model, simply using the coxph-derived p-value (or a corrected equivalent) will not cause a massive increase in false-positive results.

Session Information
-------------------

``` r
sessioninfo::session_info()
```

─ Session info ────────────────────────────────────────────────────────── setting value
version R version 3.4.1 (2017-06-30) os macOS Sierra 10.12.2
system x86\_64, darwin15.6.0
ui X11
language (EN)
collate en\_GB.UTF-8
tz Europe/London
date 2018-04-26

─ Packages ────────────────────────────────────────────────────────────── package \* version date
a4Base \* 1.24.0 2017-04-25 a4Core \* 1.24.0 2017-04-25 a4Preproc \* 1.24.0 2017-04-25 annaffy \* 1.48.0 2017-04-25 annotate 1.54.0 2017-04-25 AnnotationDbi \* 1.38.2 2017-07-27 assertthat 0.2.0 2017-04-11 audio 0.1-5 2013-12-23 backports 1.1.2 2017-12-13 beepr \* 1.2 2015-06-14 Biobase \* 2.36.2 2017-05-04 BiocGenerics \* 0.22.1 2017-10-06 BiocInstaller \* 1.26.1 2017-09-01 bit 1.1-12 2014-04-09 bit64 0.9-7 2017-05-08 bitops 1.0-6 2013-08-17 blob 1.1.0 2017-06-17 bootstrap 2017.2 2017-02-27 caTools 1.17.1 2014-09-10 clisymbols 1.2.0 2017-05-21 codetools 0.2-15 2016-10-05 colorspace 1.3-2 2016-12-14 cowplot \* 0.9.2 2017-12-17 DBI 0.7 2017-06-18 digest 0.6.15 2018-01-28 evaluate 0.10.1 2017-06-24 foreach \* 1.4.4 2017-12-12 gdata 2.18.0 2017-06-06 genefilter \* 1.58.1 2017-05-06 ggplot2 \* 2.2.1 2016-12-30 ggthemes \* 3.4.2 2018-04-03 glmnet \* 2.0-13 2017-09-22 GO.db \* 3.4.1 2017-08-18 gplots \* 3.0.1 2016-03-30 gridExtra 2.3 2017-09-09 gtable 0.2.0 2016-02-26 gtools 3.5.0 2015-05-29 hexbin \* 1.27.2 2018-01-15 highr 0.6 2016-05-09 hms 0.4.1 2018-01-24 htmltools 0.3.6 2017-04-28 IRanges \* 2.10.5 2017-10-08 iterators 1.0.9 2017-12-12 KEGG.db \* 3.2.3 2017-08-18 KernSmooth \* 2.23-15 2015-06-29 knitr \* 1.19 2018-01-29 labeling 0.3 2014-08-23 lattice 0.20-35 2017-03-25 lava 1.6 2018-01-13 lazyeval 0.2.1 2017-10-29 limma \* 3.32.10 2017-10-13 magrittr \* 1.5 2014-11-22 MASS \* 7.3-48 2017-12-25 Matrix \* 1.2-12 2017-11-15 memoise 1.1.0 2017-04-21 mpm \* 1.0-22 2011-11-25 multtest \* 2.32.0 2017-04-25 munsell 0.4.3 2016-02-13 pander 0.6.1 2017-08-06 pillar 1.1.0 2018-01-14 pkgconfig 2.0.1 2017-03-21 plyr 1.8.4 2016-06-08 prodlim 1.6.1 2017-03-06 R6 2.2.2 2017-06-17 Rcpp 0.12.15 2018-01-20 RCurl 1.95-4.10 2018-01-04 readr \* 1.1.1 2017-05-16 reshape2 \* 1.4.3 2017-12-11 rlang 0.1.6 2017-12-21 rmarkdown 1.8 2017-11-17 rmeta 2.16 2012-10-29 rprojroot 1.3-2 2018-01-03 RSQLite 2.0 2017-06-19 S4Vectors \* 0.14.7 2017-10-08 scales 0.5.0 2017-08-24 sessioninfo 1.0.0 2017-06-21 stringi 1.1.6 2017-11-17 stringr 1.3.0 2018-02-19 SuppDists 1.1-9.4 2016-09-23 survcomp 1.26.0 2017-04-25 survival \* 2.41-3 2017-04-04 survivALL \* 0.9.2.1000 2018-04-24 survivalROC 1.0.3 2013-01-13 tibble 1.4.2 2018-01-22 viridis \* 0.5.1 2018-03-29 viridisLite \* 0.3.0 2018-02-01 withr 2.1.2 2018-03-15 XML 3.98-1.9 2017-06-19 xtable 1.8-2 2016-02-05 yaml 2.1.16 2017-12-12 source
Bioconductor
Bioconductor
Bioconductor
Bioconductor
Bioconductor
Bioconductor
CRAN (R 3.4.1)
CRAN (R 3.4.0)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
Bioconductor
Bioconductor
Bioconductor
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
cran (@2017.2)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
Bioconductor
CRAN (R 3.4.1)
cran (@3.4.2)
CRAN (R 3.4.2)
Bioconductor
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.3)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
Bioconductor
CRAN (R 3.4.1)
Bioconductor
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
Bioconductor
CRAN (R 3.4.0)
CRAN (R 3.4.1)
CRAN (R 3.4.2)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
Bioconductor
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
Bioconductor
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.2)
CRAN (R 3.4.3)
CRAN (R 3.4.1)
Bioconductor
CRAN (R 3.4.1)
Github (<pearcedom/survivALL@53029f1>) CRAN (R 3.4.1)
CRAN (R 3.4.1)
cran (@0.5.1)
CRAN (R 3.4.1)
CRAN (R 3.4.4)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
CRAN (R 3.4.1)
