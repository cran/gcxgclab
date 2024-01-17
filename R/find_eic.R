#' @title  Finds EICs
#'
#' @description  `find_eic` calculates the mass defect for each ion, then finds
#' the specific EICs of interest.
#'
#' @details Extracted Ion Chromatogram (EIC) is a plot of intensity at a chosen
#' m/z value, or range of values, as a function of retention time.
#' This function finds intensity values at the given mass-to-charge (m/z)
#' values, MOI, and in a range around MOI given a tolerance. Calculates the mass
#' defect for each ion, then finds the specific EICs of interest. Returns a data
#' frame of time values, mass values, intensity values, and mass defects.
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param MOI a \emph{float} object. The mass (m/z) value of interest.
#' @param tolerance a \emph{double} object. The tolerance allowed for the MOI.
#' Default is 0.0005.
#'
#' @return eic, a \emph{data.frame} object. A data frame of time values,
#' retention time 1, retention time 2, mass values, intensity values, and mass
#' defects.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file1,mod_t=.5)
#' eic <- find_eic(frame, MOI=92.1397,tolerance=0.005)
#' plot_eic(eic,dim=1,title='EIC for MOI 92.1397')
#' plot_eic(eic,dim=2,title='EIC for MOI 92.1397')
#'
#' @export
find_eic <- function(data, MOI, tolerance=0.0005){
  # Finding specific EICs of interest
  ms_data <- data[[2]]
  eic <- dplyr::filter(ms_data, ms_data$mass_values <= MOI +
                       tolerance & ms_data$mass_values >= MOI - tolerance)
  return(eic)
}


#' @title  Plots the EICs
#'
#' @description  `plot_eic` Plots the EICs
#'
#' @details This function produces a scatter plot of the overall time index vs
#' the intensity values at a given mass of interest using
#' \code{\link[ggplot2]{ggplot}} from ggplot2 package
#' \insertCite{ggplot2}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param eic a \emph{data.frame} object. A data frame of the times and
#' intensity values of the EIC of interest, ideally the output of find_eic().
#' @param title a \emph{string} object. Title placed at the top of the plot.
#' Default title "EIC".
#' @param dim a \emph{integer} object. The time dimensions of the plot, either 1
#' or 2. Default is 1.
#'
#' @return A \emph{ggplot} object. A scatter plot of the overall time index vs
#' the intensity values at a given mass of interest.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file1,mod_t=.5)
#' eic <- find_eic(frame, MOI=92.1397,tolerance=0.005)
#' plot_eic(eic,dim=1,title='EIC for MOI 92.1397')
#' plot_eic(eic,dim=2,title='EIC for MOI 92.1397')
#'
#' @export
plot_eic <- function(eic,title="EIC",dim=1){
  if (dim==1){
    time_array <- eic$'time_array'
    intensity_values <- eic$'intensity_values'
    fig <- ggplot2::ggplot(eic, ggplot2::aes(x = time_array,
                    y = intensity_values)) +
           ggplot2::geom_point() +
           ggplot2::scale_fill_viridis_c() +
           ggplot2::labs(title= title, x='time array',
                    y= 'intensity')+
           ggplot2::theme_light()+
           ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  }
  else if (dim==2){
    #names(eic) <- c('time_array','RT1','RT2','mass_values','TIC','ms_defect')
    #fig <- plot_chr(eic, title=title, scale='log')
    RT1 <- eic$'RT1'
    RT2 <- eic$'RT2'
    Int <- log(eic$'intensity_values',10)
    fig <- ggplot2::ggplot(eic, ggplot2::aes(x = RT1, y = RT2,
                                             color=Int)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_viridis_c() +
      ggplot2::labs(title= title, x='retention time 1',
                    y= 'retention time 2', color = 'Log Intensity')+
      ggplot2::theme_light()+
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::guides(color = ggplot2::guide_colorbar(ticks = FALSE))
  }
  else {stop("dim input must be either 1 or 2.")}

  return(fig)
}


#' @title  Finds batch of EICs
#'
#' @description  `batch_eic` calculates the mass defect for each ion, then finds
#' each listed EICs of interest.
#'
#' @details Extracted Ion Chromatogram (EIC) is a plot of intensity at a chosen
#' m/z value, or range of values, as a function of retention time.
#' This function uses find_eic() to find intensity values at the given
#' mass-to-charge (m/z) values, MOIs, and in a range around MOI given a
#' tolerance. Calculates the mass defect for each ion, then finds the specific
#' EICs of interest. Returns a data frame of time values, mass values, intensity
#' values,and mass defects.
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param MOIs a \emph{vector} object. A vector containing a list of all masses
#' of interest to be investigated.
#' @param tolerance a \emph{double} object. The tolerance allowed for the MOI.
#' Default is 0.0005.
#'
#' @return eic_list, \emph{list} object, containing \emph{data.frame} objects.
#' Data frames of time values, mass values, intensity values, and mass defects
#' for each MOI listed in the input csv or txt file.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file1,mod_t=.5)
#' mois <- c(92.1397, 93.07058)
#' eics <- batch_eic(frame, MOIs=mois ,tolerance = 0.005)
#' for (i in 1:length(eics)){
#'    print(plot_eic(eics[[i]], title=paste("EIC for MOI",mois[i])))
#'    print(plot_eic(eics[[i]], title=paste("EIC for MOI",mois[i]), dim=2))
#' }
#'
#' @export
batch_eic <- function(data, MOIs, tolerance=0.0005){
  ms_data <- data[[2]]
  list <- list()
  for (i in (1:length(MOIs))){
    eic <-dplyr::filter(ms_data, ms_data$mass_values <= MOIs[i] + tolerance &
                        ms_data$mass_values >= MOIs[i] - tolerance)
    list <- append(list, list(eic))
  }
  names(list)<-MOIs
  return(list)
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
