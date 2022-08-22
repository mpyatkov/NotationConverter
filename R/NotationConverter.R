# MM10_SEGEX_LOOKUP <- read_csv("./mm10_to_segex.csv", col_names = T) %>%
#     select(mm9 = gene_long_mm10, mm10 = gene_short_mm10, segex = probe_id, in_segex)

#' Allows convert gene name between "mm9", "mm10" and "segex" notations
#'
#' @import purrr dplyr tidyr
#' @importFrom stats setNames
#' @importFrom rlang :=
#' @param input_df data.frame with column you are going to transform from one notation to another
#' @param from character string denote input notations ("mm9", "mm10" or "segex")
#' @param to character string denote output notations ("mm9", "mm10" or "segex")
#' @param column character string which column will be used as input
#' @param replace_column bool, if TRUE replaces "column" with converted values. if FALSE added converted column to the output data.frame
#'
#' @return data.frame with converted notation. NOTE: if replace_column=TRUE and input data.frame has gene name which does not exist in output notation (fake name), then fake names will be replaced with NA
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' test.df <- tibble(gname = c("lnc50075", "lnc50076", "lnc_fake_name"), id=c(1,2,3))
#' # fake_name will be NA in the output data frame
#' output <- notationConverter(test.df, from = "mm10", to = "segex")

notationConverter <- function(input_df,
                               from = "mm9",
                               to = "mm10",
                               column = "gname",
                               replace_column = T) {

    ## checking for correct from and to parameters
    if (!is_empty(setdiff(c("mm10", "mm9"), c("mm9", "mm10", "segex")))) {
        stop("Something wrong with 'from' and 'to' parameters should be 'mm9','mm10' and 'segex' only")
    }

    ## checking for correct 'column' parameter
    if (!(column %in% colnames(input_df))) {
        stop(paste0("Cannot find '",column,"' in input dataset"))
    }

    if (from == to) {
        stop("'from' and 'to' parameters should be different")
    }

    ## mm9 to mm10
    tmp <- left_join(input_df,
                     MM10_SEGEX_LOOKUP %>%
                         select(!!from, !!to),
                     by=setNames(nm=column, from))

    if (replace_column){
        tmp <- tmp %>%
            select(-!!column) %>%
            rename(!!column := !!sym(to)) %>%
            select(!!column, everything())
    }
    tmp
}

#' Allows convert any proper data.frame to Segex upload file
#'
#' @import purrr dplyr tidyr
#' @importFrom rlang .data
#' @param input_df data.frame with 5 columns:
#'
#' - gene_name
#' - intensity 1
#' - intensity 2
#' - log2fc
#' - adj_pvalue
#'
#' @param from character string "mm9" or "mm10" depends on which notation represents the 'gene name' column
#'
#' @return list with two data.frames
#'
#' - to_segex is a data.frame with 75598 rows which ready to upload to Segex database. NOTE: intensity_1 and intensity_2 swaped because Segex database uses Condtion2/Condtion1 notation
#' - appendix is a data.frame which is empty if "from" option equal to "mm9", for "mm10" it represents the genes (403 genes) which are not represented in Segex database
#' @export
#'
#' @examples
#' library(dplyr)
#' library(readr)
#'
#' test.df <- tibble(gname = "lnc50075",
#'                   i1 = 0.1,
#'                   i2 = 2.0,
#'                   log2fc = -4.321928,
#'                   adj_pvalue = 0.05)
#'
#' res <- exportToSegex(test.df, from = "mm10")
#'
#' # Segex can parse filenames, so to facilitate uploading just give the proper name for output file
#' # fn <- "1_scLoupe_SAMPLEID_Control_vs_SAMPLEID_Treatment_DiffExp_IntronicMonoExonic.tsv"
#' # write_tsv(res, fn, col_names = T)
#'
exportToSegex <- function(input_df, from = "mm9") {

    ## input_df should contain 5 columns:
    ## gname, intensity_1, intensity_2, log2fc, adj_pvalue

    if (!(from %in% c("mm9", "mm10"))) {
        stop("Option 'from' is not correct, should be mm9 or mm10")
    }

    ratio <- gname <- id.1.intensity <- id.2.intensity <- log2_fold_change <- adjusted_p_value <- in_segex <- segex <- fc <- NULL
    ## linearize FC, compute ratio

    tmp <- input_df %>%
        select(gname = 1, id.1.intensity = 2, id.2.intensity = 3,log2_fold_change =4 , adjusted_p_value = 5) %>%
        mutate(ratio = 2^log2_fold_change,
               fc = ifelse(ratio < 1, -1/ratio, ratio))

    ## TODO: fix problem with dot(.data) operator because it seems that .data does not replace it dot (.) operator.
    ##       tmp should be replaced to .data in expression below to see the difference
    tmp <- left_join(MM10_SEGEX_LOOKUP %>% select(segex, gname = !!from, in_segex), tmp, by = "gname") %>%
        replace_na(list(ratio = 1, fc = 1, id.1.intensity = 0, id.2.intensity = 0, adjusted_p_value = 1))

    to_segex <- tmp %>%
        filter(in_segex == TRUE) %>%
        ## SWAP ORDER OF INTENSITIES, BECAUSE SEGEX USES INTENSITY1/INTENSITY2 NOTATION
        select(segex, ratio, fc, intensity_1 = id.2.intensity, intensity_2 = id.1.intensity, pvalue_1 = adjusted_p_value)

    ## there is no record about these genes in segex
    appendix <- tmp %>%
        filter(in_segex == FALSE) %>%
        ## SWAP ORDER OF INTENSITIES, BECAUSE SEGEX USES INTENSITY1/INTENSITY2 NOTATION
        select(segex, ratio, fc, intensity_1 = id.2.intensity, intensity_2 = id.1.intensity, pvalue_1 = adjusted_p_value)

    list(segex = to_segex, appendix = appendix)
}

