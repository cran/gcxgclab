% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gauss_fit.R
\name{gauss_fit}
\alias{gauss_fit}
\title{Fitting to Gaussian curve}
\usage{
gauss_fit(TIC_df, peakcoord)
}
\arguments{
\item{TIC_df}{a \emph{data.frame} object. Data frame with 4 columns
(Overall Time Index, RT1, RT2, TIC), ideally the output from create_df(), or
the first data frame returned from extract_data(), $TIC_df.}

\item{peakcoord}{a \emph{vector} object. The two dimensional time retention
coordinates of the peak of interest. c(RT1,RT2).}
}
\value{
A \emph{list} object with three items. The first \emph{data.frame}
object. A data frame with two columns, (time, guassfit), the time values
around the peak, and the intensity values fitted to the optimal Gaussian
curve. Second, a \emph{vector} object of the fitted parameters (a,b,c).
Third, a \emph{double} object, the area under the fitted Gaussian curve.
}
\description{
\code{gauss_fit} fits data around a peak to a Gaussian curve.
}
\details{
This function fits data around the specified peak to a Gaussian
curve, minimized with nonlinear least squares method nls() from "stats"
package.
}
\examples{
file <- system.file("extdata","sample1.cdf",package="gcxgclab")
frame <- extract_data(file,mod_t=.5)
peaks <- top_peaks(frame$TIC_df, 5)
gaussfit <- gauss_fit(frame$TIC_df, peakcoord=c(peaks$'X'[1], peaks$'Y'[1]))
message(paste('Area under curve =',gaussfit[[3]], 'u^2'))
plot_gauss(frame$TIC_df, gaussfit[[1]])

}
