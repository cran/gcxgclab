% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_peaks.R
\name{plot_peakonly}
\alias{plot_peakonly}
\title{Plot only peaks}
\usage{
plot_peakonly(peak_df, title = "Peaks")
}
\arguments{
\item{peak_df}{a \emph{data.frame} object. A data frame with 4 columns
(Time, X, Y, Peak), ideally the output from top_peaks() or thr_peaks().}

\item{title}{a \emph{string} object. Title placed at the top of the plot.
Default title "Peaks".}
}
\value{
A \emph{ggplot} object. A circle plot of peak intensity in 2D
retention time.
}
\description{
\code{plot_peakonly} plots the peaks from a chromatograph.
}
\details{
This function creates a circle plot of the peak intensity vs
the x and y retention times using \code{\link[ggplot2]{ggplot}} from ggplot2
package \insertCite{ggplot2}{gcxgclab}. The size of the circle indicates the
intensity of the peak.
}
\examples{
file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
frame <- extract_data(file1,mod_t=.5)
peaks <- top_peaks(frame$TIC_df, 5)
plot_peakonly(peaks,title="Top 20 Peaks")

}
\references{
\insertAllCited{}
}
