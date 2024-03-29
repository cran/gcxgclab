% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/non_targeted.R
\name{non_targeted}
\alias{non_targeted}
\title{Compares MS to atomic mass data}
\usage{
non_targeted(masslist, ms, THR = 0.1, ...)
}
\arguments{
\item{masslist}{a \emph{list} object, a list of atomic weights, ideally the
output of mass_list().}

\item{ms}{a \emph{data.frame} object, a data frame of the mass values and the
percent intensity values, ideally the output of find_ms().}

\item{THR}{a \emph{double} object. The threshold of intensity of which to
include peaks for mass comparison. Default is 0.1.}

\item{...}{a \emph{vector} object. Any further optional inputs which
indicate additional elements to consider in the compound, or restrictions on
the number of a certain element in the compound. Should be in the form
c('X', a, b) where X = element symbol, a = minimum number of atoms, b =
maximum number of atoms. a and b are optional. If no minimum, use a=0, if no
maximum, do not include b.}
}
\value{
A \emph{list} object, a list of vectors containing strings of the
matching compounds.
}
\description{
\code{non_targeted} compares the MS data from a peak to atomic mass
data.
}
\details{
This function takes the MS data from an intensity peak in a sample
and compares it to combinations of atomic masses. Then it approximates the
makeup of the compound, giving the best matches to the MS data. Note that
the default matches will contain only H, N, C, O, F, Cl, Br, I, and Si. The
user can input optional parameters to indicate additional elements to be
considered or restrictions on the number of any specific element in the
matching compounds.
}
\examples{
file <- system.file("extdata","sample1.cdf",package="gcxgclab")
frame <- extract_data(file,mod_t=.5)
peaks <- top_peaks(frame$TIC_df, 5)
mz <- find_ms(frame, t_peak=peaks$'T'[1])
masslist <- mass_list()
non_targeted(masslist, mz, THR=0.05)

}
