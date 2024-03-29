% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/non_targeted.R
\name{mass_list}
\alias{mass_list}
\title{Creates list of atomic mass data}
\usage{
mass_list()
}
\value{
A \emph{data.frame} object, with two columns, (elements,
mass).
}
\description{
\code{mass_list} creates a list of atomic mass data
}
\details{
This function creates a data frame containing the data for the
atomic weights for each element in the periodic table
\insertCite{mass}{gcxgclab}.
}
\examples{
file <- system.file("extdata","sample1.cdf",package="gcxgclab")
frame <- extract_data(file,mod_t=.5)
peaks <- top_peaks(frame$TIC_df, 5)
mz <- find_ms(frame, t_peak=peaks$'T'[1])
masslist <- mass_list()
non_targeted(masslist, mz, THR=0.05)

}
\references{
\insertAllCited{}
}
