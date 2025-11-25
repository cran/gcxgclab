#' @title  Baseline correction
#'
#' @description  `bl_corr` performs baseline correction of the intensity values.
#'
#' @details This function performs baseline correction and baseline subtraction
#' for TIC values.
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param gamma a \emph{float} object. Correction factor between 0 and 1. 0
#' results in almost no values being subtracted to the baseline, 1 results in
#' almost everything except the peaks to be subtracted to the baseline. Default
#' is 0.5.
#' @param subtract a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#'
#' @return A \emph{data.frame} object. A data frame of the overall time index,
#' the x-axis retention time, the y-axis retention time, and the baseline
#' corrected total intensity values.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' sm_frame <- smooth(frame, lambda=10)
#' blc_frame <- bl_corr(sm_frame, gamma=0.5)
#' plot_chr(blc_frame, title='Baseline Corrected')
#'
#' @export
bl_corr <- function(data, gamma=0.5, subtract=NULL){
  frame <- data[[1]]
  if (gamma>1 | gamma<0){
    stop("gamma must be between 0 and 1.")
  }
  gamma <- 3*gamma-1
  N <- length(frame$TIC)
  if (!(is.null(subtract))){
    subtract <- subtract[[1]]
    if (length(subtract$'TIC')!=length(frame$'TIC')){
      stop('To perform baseline subtractions, samples must be of the same size.')
    }
    else{
      frame$'TIC' <- frame$'TIC'-subtract$'TIC'
    }
  }
  bl <- c(rep(NA,N))
  sig <- c(rep(NA,N))
  win_size <- 100
  if (win_size > N){
    win_size <- 1
  }
  else if (N<10000){
    win_size <- round(N/100)
  }
  points <- c(1:round(N/win_size))
  points <- points[1:(length(points)-1)]*win_size
  for (i in c(1:N)){
    if (i %in% points){
      if (i== points[1]){
        window <- c(1:i+win_size)
      }
      else if (i== points[length(points)]){
        window <- c(i-win_size,N)
      }
      else {
        window <- c((i-win_size):(i+win_size))
      }
      bl[i] <- min(frame$TIC[window])
    }
  }
  bl[1] <- bl[points[1]]
  bl[length(bl)] <- bl[points[length(points)]]
  min <- bl
  bl2 <- bl[is.na(bl)==FALSE]
  bl2 <- ptw::whit1(bl2,1000)
  sdev <- stats::sd(bl,na.rm=TRUE)
  bl3 <- bl2+gamma*sdev
  bl[c(1,points,N)] <- bl3
  bl <- zoo::na.approx(bl)
  frame$TIC <- frame$TIC-bl
  frame$TIC[frame$TIC<0]<- 0
  data[[1]] <- frame

  return(data)
}


#' @title  Smoothing
#'
#' @description  `smooth` performs smoothing of the intensity values.
#'
#' @details This function performs smoothing of the intensity values using
#' Whittaker smoothing algorithm \code{\link[ptw]{whit1}} from the ptw package
#' \insertCite{ptw}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param lambda a \emph{float} object. A number (parameter in Whittaker
#' smoothing), suggested between 0 to 10^4. Small lambda is very little
#' smoothing, large lambda is very smooth. Default is lambda = 20.
#' @param dir a \emph{string} object. Either "X", "Y", or "XY" to indicate
#' direction of smoothing. "XY" indicates smoothing in both X (horizontal) and Y
#' (vertical) directions. Default "XY".
#'
#' @return A \emph{data.frame} object. A list of two data frames. A TIC data
#' frame and an MS data frame.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' sm_frame <- smooth(frame, lambda=10)
#' plot_chr(sm_frame, title='Smoothed')
#'
#' @export
smooth <- function(data, lambda=20, dir='XY'){
  frame <- data[[1]]
  if (lambda <0){
    warning('Lambda value is too low. Choose lambda greater than or equal to 0.')
  }
  else if (lambda ==0 ){
    return(data)
  }
  else{
    if (lambda >1000){
      warning('Lambda value is too high. Choose lambda less than 1000.')
    }
    if(!(dir=='X' | dir=='Y' | dir=='XY')){
      stop('Smoothing direction, dir, must be either "X", "Y", or "XY".')
    }
    RT1 <- frame$'RT1'
    RT2 <- frame$'RT2'
    if (dir=='X' | dir=='XY'){
      frame <- dplyr::arrange(frame,RT2,RT1)
      frame$TIC <- ptw::whit1(c(frame$TIC),lambda/100)
      frame <- dplyr::arrange(frame,RT1,RT2)}
    if (dir=='Y' | dir== "XY"){
      frame$TIC <- ptw::whit1(c(frame$'TIC'),lambda)
    }
    data[[1]] <- frame
    return(data)
  }
}


#' @title Phase shift
#'
#' @description  `phase_shift` shifts the phase of the chromatogram.
#'
#' @details This function shifts the phase of the chromatogram up or down by
#' the specified number of seconds.
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param shift a \emph{float} object. The number of seconds to shift the phase
#' by.
#'
#' @return A \emph{data.frame} object. A list of two data frames. A TIC data
#' frame and an MS data frame.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' shifted <- phase_shift(frame, -.2)
#' plot_chr(shifted, title='Shifted')
#'
#' @export
phase_shift <- function(data, shift){
  frame <- data[[1]]
  ms_data <- data[[2]]
  maxy <- max(frame$RT2)
  miny <- min(frame$RT2)
  mod <- maxy-miny
  if (abs(shift)>=mod){
    stop(paste0("Shift must be a number between -mod_time and +mod_time, (-",
                          mod, ",", mod,")."))
  }
  ydim <- length(which(frame$RT1==frame$RT1[1]))
  xdim <- length(frame$RT1)/ydim
  minx <- min(frame$RT1)
  maxx <- max(frame$RT1)
  xdiff <- frame$RT1[ydim+1]-frame$RT1[1]
  if (shift<0){
    Tic <- frame$TIC
    len <- length(frame$TIC)
    loc_1 <- which.min(abs(frame$RT2-frame$RT2[len]-shift))
    loc_2 <- which(frame$RT2==frame$RT2[loc_1])
    loc_3 <- loc_2[length(loc_2)]
    sft <- len-loc_3+1
    Tic <- Tic[-c(1:sft)]
    Tic <- append(Tic,rep(frame$TIC[len],sft))
    frame$TIC <- Tic

    ms_data$time_array <- ms_data$time_array+shift*60
  }
  else if (shift >0){
    Tic <- frame$TIC
    len <- length(frame$TIC)
    sft <- which.min(abs(frame$RT2-frame$RT2[1]-shift))[1]-1
    Tic <- append(rep(frame$TIC[1],sft),Tic)
    Tic <- Tic[-c((len-sft+1):len)]
    if (sft>0){
      frame$TIC <- Tic
      ms_data$time_array <- ms_data$time_array+shift*60
    }
  }
  data[[1]] <- frame
  data[[2]] <- ms_data

  return(data)
}



#' @title  Preprocessing
#'
#' @description  `preprocess` performs full preprocessing on a data file.
#'
#' @details This function performs full preprocessing on a data file. Extracts
#' data and performs smoothing and baseline correction.
#'
#' @param filename a \emph{string} object. The file name or path of the cdf
#' file to be opened.
#' @param mod_t a \emph{float} object. The modulation time for the GCxGC sample
#' analysis.Default is 10.
#' @param shift a \emph{float} object. The number of seconds to shift the phase
#' by. Default is 0 to skip shifting.
#' @param lambda a \emph{float} object. A number (parameter in Whittaker
#' smoothing), suggested between 1 to 10^5. Small lambda is very little
#' smoothing, large lambda is very smooth. Default is lambda = 20.
#' @param gamma a \emph{float} object. Correction factor between 0 and 1. 0
#' results in almost no values being subtracted to the baseline, 1 results in
#' almost everything except the peaks to be subtracted to the baseline. Default
#' is 0.5.
#' @param subtract a \emph{data.frame} object. Data frame containing TIC data
#' from a background sample or blank sample to be subtracted from the sample TIC
#' data.
#' @param images a \emph{boolean} object. An optional input. If TRUE, all images
#' of preprocessing steps will be displayed. Default is FALSE, no images will be
#' displayed.
#'
#' @return A \emph{data.frame} object. A list of two data frames. A TIC data
#' frame and an MS data frame.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- preprocess(file,mod_t=.5,lambda=10,gamma=0.5,images=TRUE)
#'
#' @export
preprocess <- function(filename,mod_t=10,shift=0,lambda=20,gamma=0.5,
                       subtract=NULL, images=FALSE){
  message('Extracting data from file...')
  data <- extract_data(filename,mod_t=mod_t)
  fig1 <- plot_chr(data, title='Raw Data',scale='linear', dim=2)
  if (images){ print(fig1)}
  fig2 <- plot_chr(data,title='Log Intensity',scale='log', dim=2)
  if (images){ print(fig2)}
  if (shift!=0){
    message('Performing phase shift...')
    data <- phase_shift(data,shift)
    fig3 <- plot_chr(data,title='Shifted',scale='log', dim=2)
    if (images){ print(fig3)}
  }
  if (lambda>0){
    message('Performing smoothing...')
    data <- smooth(data,lambda)
    fig4 <- plot_chr(data,title='Smoothed')
    if (images){ print(fig4)}
  }
  if (gamma>0){
    message('Performing baseline correction...')
    data <- bl_corr(data,gamma,subtract)
    fig5 <- plot_chr(data,title='Baseline Corrected')
    if (images){ print(fig5)}
  }
  message('Complete.')
  return(data)
}

#' @title  Batch reprocessing
#'
#' @description  `batch_preprocess` performs full preprocessing on a batch of
#' data files.
#'
#' @details This function performs full preprocessing on a batch of data files.
#' Extracts data and performs peak alignment and performs smoothing and baseline
#' correction.
#'
#' @param path a \emph{string} object. The path to the directory containing the
#' cdf files to be batch preprocessed and aligned.
#' @param mod_t a \emph{float} object. The modulation time for the GCxGC sample
#' analysis. Default is 10.
#' @param shift a \emph{float} object. The number of seconds to shift the phase
#' by. Default is 0 to skip shifting.
#' @param lambda a \emph{float} object. A number (parameter in Whittaker
#' smoothing), suggested between 1 to 10^5. Small lambda is very little
#' smoothing, large lambda is very smooth. Default is lambda = 20.
#' @param gamma a \emph{float} object. Correction factor between 0 and 1. 0
#' results in almost no values being subtracted to the baseline, 1 results in
#' almost everything except the peaks to be subtracted to the baseline. Default
#' is 0.5.
#' @param subtract a \emph{data.frame} object. Data frame containing TIC data
#' from a background sample or blank sample to be subtracted from the sample TIC
#' data.
#' @param THR a \emph{float} object. Threshold for peak intensity for peak
#' alignment. Should be a number between the baseline value and the highest peak
#' intensity. Default is THR = 100000.
#' @param do_align a \emph{boolean} object. An optional input allowing the user
#' to skip alignment of the given data files if alignment is not needed. Default
#' is TRUE.
#' @param use_ref_peak a \emph{boolean} object. Determines if an initial shift
#' to a given reference peak, default is toluene, should be done before aligning
#' all other peaks above given threshold THR. Default is TRUE.
#' @param ref_peak a \emph{float} object. The m/z value of the reference peak
#' for optional initial shift. Default is 92.1397 (toluene).
#' @param images a \emph{boolean} object. An optional input. If TRUE, all images
#' of preprocessing steps will be displayed. Default is FALSE, no images will be
#' displayed.
#'
#'
#' @return A \emph{data.frame} object. A list of pairs of data frames. A TIC
#' data frame and an MS data frame for each file.
#'
#' @examples
#' folder <- system.file("extdata",package="gcxgclab")
#' frame_list <- batch_preprocess(folder,mod_t=.5,lambda=10,gamma=0.5,images=TRUE)
#'
#' @export
batch_preprocess <- function(path=".",mod_t=10,shift=0,lambda=20,gamma=0.5,
                             subtract=NULL,THR=10^5,do_align=TRUE,
                             use_ref_peak=TRUE,ref_peak=92.1397,images=FALSE){
  files <- list.files(path=path, pattern='.cdf')
  if (length(files)==0){
    stop('There are no cdf data files in this directory.')
  }
  if (length(files)==1){
    if(do_align){
      message('There is only one cdf data file in this directory. No aligment possible.')
    }
    dfs <- list(gcxgclab::preprocess(files[1],shift=shift,lambda=lambda,gamma=gamma,
                                     subtract=subtract,images=images))
    names(dfs) <- files
    return(dfs)
  }
  else{
    dfs <- list()

    for (i in 1:length(files)){
      filename <- paste0(path, '/', files[i])
      message('Extracting data from file...')
      data <- extract_data(filename,mod_t=mod_t)
      fig1 <- plot_chr(data, title='Raw Data',scale='linear', dim=2)
      if (images){ print(fig1)}
      fig2 <- plot_chr(data,title='Log Intensity',scale='log', dim=2)
      if (images){ print(fig2)}
      dfs <- append(dfs,list(data))
    }
    if(do_align){
      aligned <- align(dfs,THR=THR,use_ref_peak=use_ref_peak,ref_peak=ref_peak)
    }
    else{
      aligned <- dfs
    }
    for (i in 1:length(files)){
      data <- aligned[[i]]
      if (shift!=0){
        message('Performing phase shift...')
        data <- phase_shift(data,shift)
        fig3 <- plot_chr(data,title='Shifted',scale='log', dim=2)
        if (images){ print(fig3)}
      }
      dfs <- append(dfs,list(data))
    }

    for (i in 1:length(files)){
      data <- dfs[[i]]
      if (lambda>0){
        message('Performing smoothing...')
        data <- smooth(data,lambda)
        fig4 <- plot_chr(data,title='Smoothed')
        if (images){ print(fig4)}
      }
      if (gamma>0){
        message('Performing baseline correction...')
        data <- bl_corr(data,gamma,subtract)
        fig5 <- plot_chr(data,title='Baseline Corrected')
        if (images){ print(fig5)}
      }
      dfs[[i]]<-data
    }
    message('Complete.')
  }
  names(dfs) <- files
  return(dfs)
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
