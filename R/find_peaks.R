#' @title   Threshold Peaks
#'
#' @description  `thr_peaks` finds all peaks above the given threshold.
#'
#' @details  This function finds all peaks in the sample above a given intensity
#' threshold.
#'
#' @param TIC_df a \emph{data.frame} object. Data frame with 4 columns
#' (Overall Time Index, RT1, RT2, TIC), ideally the output from create_df(), or
#' the first data frame returned from extract_data(), $TIC_df.
#' @param THR a \emph{float} object. Threshold for peak intensity. Should be a
#' number between the baseline value and the highest peak intensity. Default
#' suggestion is THR = 100000.
#'
#' @return A \emph{data.frame} object. A data frame with 4 columns (Time, X, Y,
#' Peak) with all peaks above the given threshold, with their time coordinates.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file1,mod_t=.5)
#' thrpeaks <- thr_peaks(frame$TIC_df, 100000)
#' plot_peak(thrpeaks, frame, title="Peaks Above 100,000")
#' plot_peakonly(thrpeaks,title="Peaks Above 100,000")
#'
#' @export
thr_peaks <- function(TIC_df, THR=100000){
  frame <- TIC_df
  maxs <- c()
  xcoords <- c()
  ycoords <- c()
  time <- c()

  if (THR <0){
    stop('Threshold must be greater than zero.')
  }

  # finding points which are highest in "ball" around them
  len <- length(frame$'TIC')
  ydim <- length(unique(frame$RT1))
  xdim <- length(frame$TIC)/ydim
  for (i in 11:(len-11)){
    maxball <- c()
    maxmaxball <- 0
    maxball <- frame$'TIC'[(i-10):(i+10)]
    if (i > (xdim+6) && i < (len-xdim-6)){
      maxball <- append(maxball, frame$'TIC'[(i+xdim-5):(i+xdim+5)])
      maxball <- append(maxball, frame$'TIC'[(i-xdim-5):(i-xdim+5)])
    }
    if (i > (2*xdim+5) && i < len-2*xdim-5){
      maxball <- append(maxball, frame$'TIC'[(i+2*xdim-3):(i+2*xdim+3)])
      maxball <- append(maxball, frame$'TIC'[(i-2*xdim-3):(i-2*xdim+3)])
    }
    maxmaxball <- max(maxball)
    if (frame$'TIC'[i]==maxmaxball){
      maxs <- c(maxs,frame$'TIC'[i])
      xcoords <- c(xcoords,frame$'RT1'[i])
      ycoords <- c(ycoords,frame$'RT2'[i])
      time <- c(time, frame$'Overall Time Index'[i])
    }
  }

  # Finding only the peaks with TIC above the given threshold THR
  dec_maxs <- sort(maxs, decreasing=TRUE)
  dec_x <- c()
  dec_y <- c()
  dec_t <- c()

  loc2 <- which(dec_maxs<THR)
  if(length(loc2>0)){
    N  <- loc2[1]-1
  }
  else{
    N <- length(dec_maxs)
  }
  for (i in 1:N){
    loc <- which(maxs==dec_maxs[i])[1]
    dec_x <- c(dec_x,xcoords[loc])
    dec_y <- c(dec_y,ycoords[loc])
    dec_t <- c(dec_t,time[loc])
  }
  max_df <- as.data.frame(cbind(dec_t, dec_x, dec_y, dec_maxs[1:N]))
  colnames(max_df) <- c('T','X','Y','Peak')
  return(max_df)
}

#' @title   Top Peaks
#'
#' @description  `top_peaks` finds the top N highest peaks.
#'
#' @details  This function finds the top N peaks in intensity in the sample.
#'
#' @param TIC_df a \emph{data.frame} object. Data frame with 4 columns
#' (Overall Time Index, RT1, RT2, TIC), ideally the output from create_df(), or
#' the first data frame returned from extract_data(), $TIC_df.
#' @param N \emph{int} object. The number of top peaks to be found in the
#' sample. N should be an integer >=1. Default suggestion is N = 20.
#'
#' @return A \emph{data.frame} object. A data frame with 4 columns (Time, X, Y,
#' Peak) with the top N peaks, with their time coordinates.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file1,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' plot_peak(peaks, frame, title="Top 20 Peaks")
#' plot_peakonly(peaks,title="Top 20 Peaks")
#'
#' @export
top_peaks <- function(TIC_df, N){
  frame <- TIC_df
  maxs <- c()
  xcoords <- c()
  ycoords <- c()
  time <- c()

  # finding points which are highest in "ball" around them
  len <- length(frame$'TIC')
  ydim <- length(unique(frame$RT1))
  xdim <- length(frame$TIC)/ydim
  for (i in 11:(len-11)){
    maxball <- c()
    maxmaxball <- 0
    maxball <- frame$'TIC'[(i-10):(i+10)]
    if (i > (xdim+6) && i < (len-xdim-6)){
      maxball <- append(maxball, frame$'TIC'[(i+xdim-5):(i+xdim+5)])
      maxball <- append(maxball, frame$'TIC'[(i-xdim-5):(i-xdim+5)])
    }
    if (i > (2*xdim+5) && i < len-2*xdim-5){
      maxball <- append(maxball, frame$'TIC'[(i+2*xdim-3):(i+2*xdim+3)])
      maxball <- append(maxball, frame$'TIC'[(i-2*xdim-3):(i-2*xdim+3)])
    }
    maxmaxball <- max(maxball)
    if (frame$'TIC'[i]==maxmaxball){
      maxs <- c(maxs,frame$'TIC'[i])
      xcoords <- c(xcoords,frame$'RT1'[i])
      ycoords <- c(ycoords,frame$'RT2'[i])
      time <- c(time, frame$'Overall Time Index'[i])
    }
  }

  # Finding only the top N peaks and creating final data frame
  dec_maxs <- sort(maxs, decreasing=TRUE)
  dec_x <- c()
  dec_y <- c()
  dec_t <- c()
  for (i in 1:N){
    loc <- which(maxs==dec_maxs[i])
    dec_x <- c(dec_x,xcoords[loc])
    dec_y <- c(dec_y,ycoords[loc])
    dec_t <- c(dec_t,time[loc])
  }
  max_df <- as.data.frame(cbind(dec_t, dec_x, dec_y, dec_maxs[1:N]))
  colnames(max_df) <- c('T','X','Y','Peak')
  return(max_df)
}

#' @title   Peak Plot
#'
#' @description  `plot_peak` plots peaks on a chromatograph plot.
#'
#' @details  This function circles the identified peaks in a sample over a
#' chromatograph plot (ideally smoothed) using \code{\link[ggplot2]{ggplot}}
#' from ggplot2 package \insertCite{ggplot2}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param peaks a \emph{data.frame} object. A data frame with 4 columns (Time,
#' X, Y, Peak), ideally the output from either thr_peaks() or top_peaks().
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data(). Provides the
#' background GCxGC plot, created with plot_chr().
#' @param title a \emph{string} object. Title placed at the top of the plot.
#' Default title "Intensity with Peaks".
#' @param xlab a \emph{string} object. Label for the x axis. Default is
#' "retention time 1".
#' @param ylab a \emph{string} object. Label for the y axis. Default is
#' "retention time 2".
#' @param circlecolor a \emph{string} object. The desired color of the circles
#' which indicate the peaks. Default color red.
#' @param circlesize a \emph{double} object. The size of the circles which
#' indicate the peaks. Default size 5.
#'
#' @return A \emph{ggplot} object. A plot of the chromatogram heatmap, with
#' identified peaks circled in red.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file1,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' plot_peak(peaks, frame, title="Top 20 Peaks")
#'
#' @export
plot_peak <- function(peaks, data, title='Intensity with Peaks',xlab="retention time 1",ylab="retention time 2",
                      circlecolor="red", circlesize=5){
  frame <- data[[1]]
  X <- peaks$'X'
  Y <- peaks$'Y'
  peak_fig <- plot_chr(data,title=title,xlab=xlab,ylab=ylab,scale='log') +
                ggplot2::geom_point(data=peaks, ggplot2::aes(x=X,y=Y),
                color=circlecolor, size=circlesize, shape=1, fill=NA)
  return(peak_fig)
}



#' @title  Plot only peaks
#'
#' @description  `plot_peakonly` plots the peaks from a chromatograph.
#'
#' @details This function creates a circle plot of the peak intensity vs
#' the x and y retention times using \code{\link[ggplot2]{ggplot}} from ggplot2
#' package \insertCite{ggplot2}{gcxgclab}. The size of the circle indicates the
#' intensity of the peak.
#'
#' @references
#' \insertAllCited{}
#'
#' @param peak_df a \emph{data.frame} object. A data frame with 4 columns
#' (Time, X, Y, Peak), ideally the output from top_peaks() or thr_peaks().
#' @param title a \emph{string} object. Title placed at the top of the plot.
#' Default title "Peaks".
#'
#' @return A \emph{ggplot} object. A circle plot of peak intensity in 2D
#' retention time.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file1,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' plot_peakonly(peaks,title="Top 20 Peaks")
#'
#' @export
plot_peakonly <- function(peak_df,title="Peaks"){
  X <- peak_df$'X'
  Y <- peak_df$'Y'
  z <- peak_df$'Peak'
  fig <- ggplot2::ggplot(peak_df, ggplot2::aes(x = X, y = Y, size = z)) +
    ggplot2::geom_point(alpha=0.33, color = "blue") +
    ggplot2::labs(title= title, x='retention time 1',
                  y= 'retention time 2', size= 'Intensity')+
    ggplot2::theme_light() +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  return(fig)
}

# ©2022 Battelle Savannah River Alliance, LLC
# Notice: These data were produced by Battelle Savannah River Alliance, LLC
# under Contract No 89303321CEM000080 with the Department of Energy. During the
# period of commercialization or such other time period specified by DOE, the
# Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this data to
# reproduce, prepare derivative works, and perform publicly and display
# publicly, by or on behalf of the Government. Subsequent to that period, the
# Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this data to
# reproduce, prepare derivative works, distribute copies to the public,
# perform publicly and display publicly, and to permit others to do so. The
# specific term of the license can be identified by inquiry made to Contract or
# DOE. NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR
# ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
# LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
# USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR
# REPRESENTS THAT ITS USE WORLD NOT INFRINGE PRIVATELY OWNED RIGHTS.
