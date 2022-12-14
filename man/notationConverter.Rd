% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NotationConverter.R
\name{notationConverter}
\alias{notationConverter}
\title{Allows convert gene name between "mm9", "mm10" and "segex" notations}
\usage{
notationConverter(
  input_df,
  from = "mm9",
  to = "mm10",
  column = "gname",
  replace_column = T
)
}
\arguments{
\item{input_df}{data.frame with column you are going to transform from one notation to another}

\item{from}{character string denote input notations ("mm9", "mm10" or "segex")}

\item{to}{character string denote output notations ("mm9", "mm10" or "segex")}

\item{column}{character string which column will be used as input}

\item{replace_column}{bool, if TRUE replaces "column" with converted values. if FALSE added converted column to the output data.frame}
}
\value{
data.frame with converted notation. NOTE: if replace_column=TRUE and input data.frame has gene name which does not exist in output notation (fake name), then fake names will be replaced with NA
}
\description{
Allows convert gene name between "mm9", "mm10" and "segex" notations
}
\examples{
library(dplyr)

test.df <- tibble(gname = c("lnc50075", "lnc50076", "lnc_fake_name"), id=c(1,2,3))
# fake_name will be NA in the output data frame
output <- notationConverter(test.df, from = "mm10", to = "segex")
}
