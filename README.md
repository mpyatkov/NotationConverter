# NotationConverter

<!-- badges: start -->
<!-- badges: end -->

NotationConverter is internal tool  is to make conversion between mm9, mm10 and segex notations

## Installation

You can install the development version of NotationConverter like so:

``` r
# install.packages("devtools")
devtools::install_github("mpyatkov/NotationConverter")
```

## Example

Inside this package only two functions. 

notationConverter - allows convert gene name between "mm9", "mm10" and "segex" notations:

``` r
library(NotationConverter)
library(dplyr)

test.df <- tibble(gname = c("lnc50075", "lnc50076", "lnc_fake_name"), id=c(1,2,3))
# fake_name will be NA in the output data frame
output <- notationConverter(test.df, from = "mm10", to = "segex")
```

exportToSegex - allows convert any proper data.frame to Segex upload file

``` r
library(NotationConverter)
library(dplyr)
library(readr)

test.df <- tibble(gname = "lnc50075",
                   i1 = 0.1,
                   i2 = 2.0,
                   log2fc = -4.321928,
                   adj_pvalue = 0.05)

res <- exportToSegex(test.df, from = "mm10")

# Segex can parse filenames, so to facilitate uploading 
# just give the proper name for output file. Example:
# fn <- "1_dataframe_SAMPLEID_Treatment_vs_SAMPLEID_Control_DiffExp_IntronicMonoExonic.tsv"
# write_tsv(res, fn, col_names = T)
```


