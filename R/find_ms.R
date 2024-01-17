#' @title  Finds MS
#'
#' @description  `find_ms` Finds mass spectra of a peak.
#'
#' @details This function finds the mass spectra values of a peak in the
#' intensity values of a GCxGC sample at a specified overall time index value.
#' Then outputs a data frame of the mass values and percent intensity values
#' which can then be plotted to product the mass spectra plot.
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param t_peak a \emph{float} object. The overall time index value for when
#' the peak occurs in the GCxGC sample (the 1D time value).
#' @param tolerance a \emph{double} object. The tolerance allowed for the time
#' index. Default is 0.0005.
#'
#' @return A \emph{data.frame} object. A data frame of the mass values and the
#' percent intensity values.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' mz <- find_ms(frame, t_peak=peaks$'T'[1])
#' plot_ms(mz)
#' plot_defect(mz,title="Kendrick Mass Defect, CH_2")
#'
#' @export
find_ms <- function(data, t_peak, tolerance=0.0005){
  ms_data <- data[[2]]
  mass <- dplyr::filter(ms_data, ms_data$time_array <= t_peak/60+tolerance &
                          ms_data$time_array>=t_peak/60-tolerance)
  mz <- as.data.frame(cbind(mass$mass_values,
                            mass$intensity_values/max(mass$intensity_values)))
  names(mz) <- c('MZ','Int')
  filtered <- dplyr::filter(mz, mz$MZ>150)$'Int'
  if (length(filtered>0)){
    if (max(filtered)<0.05){
      mz <-dplyr::filter(mz, mz$MZ<150)
    }
  }
  return(mz)
}


#' @title  Plots the mass spectra of a peak.
#'
#' @description  `plot_ms` Plots the mass spectra of a peak.
#'
#' @details This function produces a line plot of the mass spectra data. The
#' mass values vs the percent intensity values as a percent of the highest
#' intensity using \code{\link[ggplot2]{ggplot}} from ggplot2 package
#' \insertCite{ggplot2}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param ms a \emph{data.frame} object. A data frame of the mass values and the
#' percent intensity values, ideally the output of find_ms().
#' @param title a \emph{string} object. Title placed at the top of the plot.
#' Default title "Mass Spectrum".
#'
#' @return A \emph{ggplot} object. A line plot of the mass spectra data. The
#' mass values vs the percent intensity values as a percent of the highest
#' intensity.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' mz <- find_ms(frame, t_peak=peaks$'T'[1])
#' plot_ms(mz)
#'
#' @export
plot_ms <- function(ms,title="Mass Spectrum"){
  mz <- c()
  int <- c()
  for (i in 1:length(ms$'MZ')){
    mz <- append(mz, c(ms$'MZ'[i], ms$'MZ'[i],ms$'MZ'[i]))
    int <- append(int, c(0, ms$'Int'[i],0))
  }
  mss <- as.data.frame(cbind(mz, int))
  colnames(mss) <- c('x','y')
  fig <- ggplot2::ggplot(mss, ggplot2::aes(x = mz, y = int)) +
    ggplot2::geom_line() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title= title, x='m/z value',
                  y= 'percent intensity')+
    ggplot2::theme_light()+
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  return(fig)
}



#' @title  Finds batch of mass spectra
#'
#' @description  `batch_ms` Finds batch of mass spectra of peaks.
#'
#' @details This function uses find_ms() to find the mass spectra values of
#' a batch list of peaks in intensity values of a GCxGC sample at overall time
#' index values specified in a txt or csv file. It outputs a list of data
#' frames, for each peak, of the mass values and percent intensity values which
#' can then be plotted to product the mass spectra plot.
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param t_peaks a \emph{vector} object. A list of times at which the peaks of
#' interest are located in the overall time index for the sample.
#' @param tolerance a \emph{double} object. The tolerance allowed for the time
#' index. Default is 0.0005.
#'
#' @return A \emph{list} object of \emph{data.frame} objects. Each a data frame
#' of the mass values and the percent intensity values.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' mzs <- batch_ms(frame, t_peaks = peaks$'T'[1:5])
#' for (i in 1:length(mzs)){
#'    print(plot_ms(mzs[[i]], title=paste('Mass Spectrum of peak', i)))
#' }
#'
#' @export
batch_ms <- function(data, t_peaks, tolerance=0.0005){
  ms_data <- data[[2]]
  list <- list()
  for (i in t_peaks){
    mass <- dplyr::filter(ms_data, ms_data$time_array <= i/60+tolerance &
                            ms_data$time_array>=i/60-tolerance)
    mz <- as.data.frame(cbind(mass$mass_values,
                              mass$intensity_values/max(mass$intensity_values)))
    names(mz) <- c('MZ','Int')
    filtered <- dplyr::filter(mz, mz$MZ>150)$'Int'
    if (length(filtered)>0){
      if (max(filtered)<0.05){
        mz <-dplyr::filter(mz, mz$MZ<150)
      }
    }
    list <- append(list, list(mz))
  }
  return(list)
}

#' @title  Plots the Kendrick Mass Defect of a peak
#'
#' @description  `plot_defect` Plots Kendrick Mass Defect of a peak.
#'
#' @details This function produces a scatter  plot of the Kendrick mass defects
#' for mass spectrum data. Plotted using \code{\link[ggplot2]{ggplot}} from
#' ggplot2 package \insertCite{ggplot2}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param ms a \emph{data.frame} object. A data frame of the mass values and the
#' percent intensity values, ideally the output of find_ms().
#' @param compound_mass a \emph{float} object. The exact mass, using most common
#' ions, of the desired atom group to base the Kendrick mass on. Default is
#' 14.01565, which is the mass for CH_2.
#' @param title a \emph{string} object. Title placed at the top of the plot.
#' Default title "Kendrick Mass Defect".
#'
#' @return A \emph{ggplot} object. A line plot of the mass spectra data. The
#' mass values vs the percent intensity values as a percent of the highest
#' intensity.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' mz <- find_ms(frame, t_peak=peaks$'T'[1])
#' plot_ms(mz)
#' plot_defect(mz,title="Kendrick Mass Defect, CH_2")
#'
#' @export
plot_defect <- function(ms,compound_mass = 14.01565,title="Kendrick Mass Defect"){
  ms <- dplyr::filter(ms,ms$Int>0.0001)
  K <- round(compound_mass)/compound_mass
  Km <- ms$MZ*K
  msdf <- round(Km)-Km
  df <- as.data.frame(cbind(Km,msdf,ms$Int))
  colnames(df)<- c('Kendrick Mass', 'Kendrick Mass Defect', 'Intensity')

  fig <- ggplot2::ggplot(df, ggplot2::aes(x = Km, y = msdf)) +
    ggplot2::geom_point(size=((log(ms$Int,10)+4)^2),color="blue",alpha=0.25) +
    ggplot2::labs(title= title, x='Kendrick Mass',
                  y= 'Kendrick Mass Defect')+
    ggplot2::theme_light()+
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
