#' @title  Extracts data from cdf file.
#'
#' @description `extract_data` Extracts the data from a cdf file.
#'
#' @details This function opens the specified cdf file using the implemented
#' function \code{\link[ncdf4]{nc_open}} from ncdf4 package, then extracts the
#' data and closes the cdf file using the implemented function
#' \code{\link[ncdf4]{nc_close}} from ncdf4 package
#' \insertCite{ncdf4}{gcxgclab}. It then returns a list of two data frames. The
#' first is a dataframe of the TIC data, the output of create_df(). The second
#' is a data frame of the full MS data, the output of mass_data().
#'
#' @references
#' \insertAllCited{}
#'
#' @param filename a \emph{string} object. The path or file name of the cdf
#' file to be opened.
#' @param mod_t a \emph{float} object. The modulation time for the GCxGC sample
#' analysis. Default is 10.
#' @param shift_time a \emph{boolean} object. Determines whether the Overall
#' Time Index should be shifted to 0. Default is TRUE.
#'
#' @return A \emph{list} object. A list of the extracted data: scan acquisition
#' time, total intensity, mass values, intensity values, and point count.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' plot_chr(frame, title='Raw Data', scale="linear")
#' plot_chr(frame, title='Log Intensity')
#'
#' @export
extract_data <- function(filename,mod_t=10,shift_time=TRUE){
  # not used:
  if (1==0){
    Rdpack::append_to_Rd_list(NULL)
  }

  cdf_file <- ncdf4::nc_open(filename)
  scan_acquisition_time <- ncdf4::ncvar_get(cdf_file, "scan_acquisition_time")
  point_count <- ncdf4::ncvar_get(cdf_file, "point_count")
  total_intensity <- ncdf4::ncvar_get(cdf_file, "total_intensity")
  mass_values <- ncdf4::ncvar_get(cdf_file, "mass_values")
  intensity_values <- ncdf4::ncvar_get(cdf_file, "intensity_values")
  ncdf4::nc_close(cdf_file)

  # creating TIC_df
  # creating a 2D time array
  N <- length(scan_acquisition_time)
  if (shift_time){
    shift <- (scan_acquisition_time - scan_acquisition_time[1])
  }
  else {
    shift <- scan_acquisition_time
  }
  rt_1d <- scan_acquisition_time/60
  acq_rate <- mean(diff(rt_1d))
  mod_period <- mod_t/60
  mod_t_exact <- round(mod_period / acq_rate, digits = 0) *acq_rate
  cycles_2D <- mod_t_exact/acq_rate
  cuts_2d <- N/cycles_2D
  rt_2d <- shift[1:cycles_2D]
  rt_2d_array <- rep(rt_2d, times = cuts_2d)
  if (length(rt_2d_array)<N){
    rt_2d_array <- rep(rt_2d, times = (cuts_2d+1))
  }
  rep_vector <- (1:cuts_2d)*cycles_2D
  rep_1d <- shift[rep_vector]
  rep_1d_fi <- rep(rep_1d, each = cycles_2D)
  if (length(rep_1d_fi)<N){
    rep_1d_fi <- rep(rep_1d, each = (cycles_2D+1))
  }
  dimensions <- c(cuts_2d,cycles_2D)

  rt_frame <- as.data.frame(cbind(shift, rep_1d_fi, rt_2d_array,
                                  total_intensity))

  colnames(rt_frame) <- c('Overall Time Index',
                          'RT1', # Retention Time 1, x axis
                          'RT2', # Retention Time 2, y axis
                          'TIC') # TIC stands for total intensity chromatograph

  # creating MS_df
  if (shift_time){
    time_array <- rep((scan_acquisition_time -
                         scan_acquisition_time[1])/60, point_count)
  }
  else{
    time_array <- rep(scan_acquisition_time/60, point_count)
  }
  rt1 <- rep(rep_1d_fi, point_count)
  rt2 <- rep(rt_2d_array, point_count)
  ms_data <- as.data.frame(cbind(time_array, rt1, rt2, mass_values,
                                 intensity_values))
  names(ms_data) <- c('time_array', 'RT1', 'RT2', 'mass_values',
                      'intensity_values')
  ms_data <- dplyr::filter(ms_data, ms_data$intensity_values > 10)

  data3 <- list(rt_frame,ms_data)
  names(data3) <- c("TIC_df", "MS_df")

  return(data3)
}



#' @title  Plot chromatogram
#'
#' @description  `plot_chr` plots TIC data for chromatogram.
#'
#' @details This function creates a contour plot using of TIC data vs the x and
#' y retention times using \code{\link[ggplot2]{ggplot}} from ggplot2 package
#' \insertCite{ggplot2}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param scale a \emph{string} object. Either 'linear' or 'log'. log refers to
#' logarithm base 10. Default is log scale.
#' @param dim a \emph{integer} object. The time dimensions of the plot, either 1
#' or 2. Default is 2.
#' @param floor a \emph{float} object. The floor value for plotting. Values
#' below floor will be scaled up. Default for linear plotting is 0, default for
#' log plotting is 10^3.
#' @param title a \emph{string} object. Title placed at the top of the plot.
#' Default title "Intensity".
#'
#' @return A \emph{ggplot} object. A contour plot of TIC data plotted in two
#' dimensional retention time.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' plot_chr(frame, title='Raw Data', scale="linear")
#' plot_chr(frame, title='Log Intensity')
#'
#' @export
plot_chr <- function(data,scale="log",dim=2,floor=-1,title="Intensity"){
  frame <- data[[1]]
  time_array <- frame$'Overall Time Index'
  TIC <- frame$'TIC'
  if (dim==2){
    RT1 <- frame$'RT1'
    RT2 <- frame$'RT2'
  }
  else if (dim!=1){
    stop("dim input must be either 1 or 2.")
  }
  if (floor==-1){
    if (scale=='log'){
      floor<-10^3
    }
    else if (scale=='linear'){
      floor <- 0
    }
  }
  for (i in 1:length(TIC)){
    if (TIC[i]<floor){
      TIC[i] <- floor}
  }
  if (scale=="linear"){
    z <- TIC
    str <- 'Intensity'
  }
  else if (scale=="log"){
    for (i in 1:length(TIC)){
      if (TIC[i]<1){
        TIC[i]<-1
      }
    }
    z <- log(TIC,10)
    str <- 'Log Intensity'
  }
  else {
    stop("scale input should be either 'linear' or 'log'.")
  }
  if (title=="Intensity"){
    title <- str
  }

  breaks <- c()
  labels<- c()
  if (floor>0){
    if (scale=='log'){

      x <- floor(log(floor,10))
      breaks <- c(x)
      labels <- c(paste0('<',x))
    }
    if (scale =='linear'){
      breaks <- c(floor)
      labels <- c(paste0('<',floor))
    }

    maxz <- round(max(z))
    for (i in (floor(log(floor,10)+1)):maxz){
      breaks <- append(breaks,i)
      labels <- append(labels,i)
    }
    A <- ggplot2::scale_fill_viridis_c(breaks=breaks,labels=labels)
  }
  else {
    A <- ggplot2::scale_fill_viridis_c()
  }

  if (dim==2){
    fig <- ggplot2::ggplot(frame, ggplot2::aes(x=RT1, y=RT2, fill=z)) +
      ggplot2::geom_tile() +
      A +
      ggplot2::labs(title= title, x='retention time 1',
                    y= 'retention time 2', fill= str)+
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(ticks = FALSE))
  }
  else if (dim==1){
    fig <- ggplot2::ggplot(frame, ggplot2::aes(x = time_array,
                                             y = z)) +
      ggplot2::geom_point() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(title= title, x='time',
                    y= 'intensity')+
      ggplot2::theme_light()+
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  }
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
