#' @title  1D Gaussian function
#'
#' @description  `gauss` Defines the 1D Gaussian curve function.
#'
#' @details This function defines a 1D Gaussian curve function.
#'
#' @param a,b,c are \emph{float} objects. Parameters in R^1 for the Gaussian
#' function.
#' @param t a \emph{float} object. The independent variable in R^1 for the
#' Gaussian function.
#'
#' @return A \emph{float} object. The value of the Gaussian function at time t,
#' given the parameters input a,b,c.
#'
#' @export
gauss <- function(a,b,c,t){
  G <- a*exp(-(t-b)^2/(2*c^2))
  return(G)
}

#' @title  2D Gaussian function
#'
#' @description  `gauss2` Defines the 2D Gaussian curve function.
#'
#' @details This function defines a 2D Gaussian curve function.
#'
#' @param a,b1,b2,c1,c2  are \emph{float} objects. Parameters in R^1 for the
#' Gaussian function.
#' @param t1,t2  are \emph{float} objects. The independent variables t=(t1.t2)
#' in R^2 for the Gaussian function.
#'
#' @return A \emph{float} object. The value of the Gaussian function at time
#' t=(t1,t2) given the parameters input a,b1,b2,c1,c2.
#'
#' @export
gauss2 <- function(a,b1,b2,c1,c2,t1,t2){
  G <- a*exp(-(t1-b1)^2/(2*c1^2)-(t2-b2)^2/(2*c2^2))
  return(G)
}


#' @title  Fitting to Gaussian curve
#'
#' @description  `gauss_fit` fits data around a peak to a Gaussian curve.
#'
#' @details This function fits data around the specified peak to a Gaussian
#' curve, minimized with nonlinear least squares method nls() from "stats"
#' package.
#'
#' @param TIC_df a \emph{data.frame} object. Data frame with 4 columns
#' (Overall Time Index, RT1, RT2, TIC), ideally the output from create_df(), or
#' the first data frame returned from extract_data(), $TIC_df.
#' @param peakcoord a \emph{vector} object. The two dimensional time retention
#' coordinates of the peak of interest. c(RT1,RT2).
#'
#' @return A \emph{list} object with three items. The first \emph{data.frame}
#' object. A data frame with two columns, (time, guassfit), the time values
#' around the peak, and the intensity values fitted to the optimal Gaussian
#' curve. Second, a \emph{vector} object of the fitted parameters (a,b,c).
#' Third, a \emph{double} object, the area under the fitted Gaussian curve.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' gaussfit <- gauss_fit(frame$TIC_df, peakcoord=c(peaks$'X'[1], peaks$'Y'[1]))
#' message(paste('Area under curve =',gaussfit[[3]], 'u^2'))
#' plot_gauss(frame$TIC_df, gaussfit[[1]])
#'
#'@export
gauss_fit <- function(TIC_df, peakcoord){

  # Identifying closest peak to given coordinates which will be fit
  frame <- TIC_df
  loc <- which.min(abs(frame$RT1-peakcoord[1])+abs(frame$RT2-peakcoord[2]))
  t_star <- frame$'Overall Time Index'[loc]
  height <- frame$'TIC'[loc]

  t_1 <- t_star-2
  t_2 <- t_star+2
  t1 <- which(frame$'Overall Time Index'>=t_1)[1]
  t2 <- which(frame$'Overall Time Index'>=t_2)[1]
  if (is.na(t2)){
    t2 <- length(frame$`Overall Time Index`)
  }
  done <- 0
  i <- 1
  while (done==0){
    if (frame$'TIC'[(loc-i)]<(1/3*height) | frame$'TIC'[(loc+i)]<(1/3*height)){
      done <- 1
      tt1 <- loc-i
      tt2 <- loc+i
    }
    else if (frame$'TIC'[(loc-i)] - frame$'TIC'[(loc-i-3)]<300 |
             frame$'TIC'[(loc+i+3)] - frame$'TIC'[(loc+i)]>-300){
      done <- 1
      tt1 <- loc-i-1
      tt2 <- loc+i+1
    }
    else{ i <- i+1 }
  }
  t <- frame$'Overall Time Index'[tt1:tt2]
  y <- frame$'TIC'[tt1:tt2]
  sig <- max((frame$'Overall Time Index'[tt2]-frame$'Overall Time Index'[tt1])/3,0.001)

  # Fitting only the peak of the data to the Gaussian
  fitting <- nls.multstart::nls_multstart(y~gauss(a,b,c,t), iter=50,
                            start_lower=list(a=height-1, b=t_star-.2, c=sig-.1),
                            start_upper=list(a=height+1, b=t_star+.2, c=sig+.1),
                            control=stats::nls.control(maxiter=100, warnOnly=TRUE),
                            supp_errors = 'Y')
  alpha <- fitting$m$getPars()[1]
  beta <- fitting$m$getPars()[2]
  cat <- abs(fitting$m$getPars()[3])
  params <- c(alpha,beta,cat)
  time = c(t1:t2)
  gaussfit <- gauss(alpha,beta,cat,frame$'Overall Time Index'[time])
  gauss_r <- as.data.frame(cbind(frame$'Overall Time Index'[time],
                                      gaussfit))
  colnames(gauss_r) <- c('time', 'gaussfit')
  area <- abs(sqrt(2*pi)*alpha*cat)[[1]]
  gauss_return <- list(gauss_r,params,area)
  names(gauss_return) <- c("Gauss_df","params","area")
  return(gauss_return)
}


#' @title  Fitting to 2D Gaussian curve
#'
#' @description  `gauss2_fit` fits data around a peak to a 2D Gaussian curve.
#'
#' @details This function fits data around the specified peak to a 2D Gaussian
#' curve, minimized with nonlinear least squares method nls() from "stats"
#' package.
#'
#' @param TIC_df a \emph{data.frame} object. Data frame with 4 columns
#' (Overall Time Index, RT1, RT2, TIC), ideally the output from create_df(), or
#' the first data frame returned from extract_data(), $TIC_df.
#' @param peakcoord a \emph{vector} object. The two dimensional time retention
#' coordinates of the peak of interest. c(RT1,RT2).
#'
#' @return A \emph{list} object with three items. The first \emph{data.frame}
#' object. A data frame with three columns, (time1, time2, guassfit), the time
#' values around the peak, and the intensity values fitted to the optimal
#' Gaussian curve. Second, a \emph{vector} object of the fitted parameters
#' (a,b1,b2,c1,c2). Third, a \emph{double} object, the volume under the fitted
#' Gaussian curve.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' gaussfit2 <- gauss2_fit(frame$TIC_df, peakcoord=c(peaks$'X'[1], peaks$'Y'[1]))
#' message(paste('Volume under curve =',gaussfit2[[3]],'u^3'))
#' plot_gauss2(frame$TIC_df, gaussfit2[[1]])
#'
#'@export
gauss2_fit <- function(TIC_df, peakcoord){
  frame <- TIC_df
  loc <- which.min(abs(frame$RT1-peakcoord[1])+abs(frame$RT2-peakcoord[2]))
  t_star <- frame$'Overall Time Index'[loc]
  height <- frame$'TIC'[loc]

  ydim <- length(which(frame$'RT1'==frame$'RT1'[1]))
  t_1 <- t_star-2
  t_2 <- t_star+2
  t1 <- which(frame$'Overall Time Index'>=t_1)[1]
  t2 <- which(frame$'Overall Time Index'>=t_2)[1]
  t1_star <- frame$'RT1'[loc]
  t2_star <- frame$'RT2'[loc]
  delta_t2 <- frame$'RT2'[loc]-frame$'RT2'[loc-1]
  delta_t1 <- frame$'RT1'[loc]-frame$'RT1'[loc-ydim]
  done <- 0
  i <- 1
  while (done==0){
    if (frame$'TIC'[(loc-i)]<(1/3*height) | frame$'TIC'[(loc+i)]<(1/3*height)){
      done <- 1
      tt1 <- loc-i
      tt2 <- loc+i
    }
    else if (frame$'TIC'[(loc-i)] - frame$'TIC'[(loc-i-3)]<3000 |
             frame$'TIC'[(loc+i+3)] - frame$'TIC'[(loc+i)]>-3000){
      done <- 1
      tt1 <- loc-i
      tt2 <- loc+i
    }
    else{ i <- i+1 }
  }
  ttimes <- c((tt1-2*ydim):(tt2-2*ydim),(tt1-ydim):(tt2-ydim), tt1:tt2,
              (tt1+ydim):(tt2+ydim), (tt1+2*ydim):(tt2+2*ydim))
  ttimes <- sort(unique(ttimes))
  ttime1 <- frame$'RT1'[ttimes]
  ttime2 <- frame$'RT2'[ttimes]
  y <- frame$'TIC'[ttimes]
  sig1 <- (frame$'RT1'[(tt2+2*ydim)]-frame$'RT1'[(tt1-2*ydim)])/2
  sig2 <- (frame$'RT2'[tt2]-frame$'RT2'[tt1])/2
  time = c((t1-5*ydim):(t2-5*ydim),(t1-4*ydim):(t2-4*ydim),
           (t1-3*ydim):(t2-3*ydim),(t1-2*ydim):(t2-2*ydim),
           (t1-ydim):(t2-ydim), t1:t2,(t1+ydim):(t2+ydim),
           (t1+2*ydim):(t2+2*ydim),(t1+3*ydim):(t2+3*ydim),
           (t1+4*ydim):(t2+4*ydim), (t1+5*ydim):(t2+5*ydim))
  time <- time[time>0 & time<length(frame$TIC)]
  time <- sort(unique(time))
  full_t1 <- frame$'RT1'[time]
  full_t2 <- frame$'RT2'[time]

  # Re-wrapping for peaks close to edge
  max2 <- max(full_t2)
  if (max(full_t2)-min(full_t2)>(max(full_t2)-0.5)){
    if (abs(t2_star-max(full_t2))<abs(t2_star-min(full_t2))){
      # Adjusting for high peak
      for (i in 2:length(ttime2)){
        if (abs(t2_star-ttime2[i])>2.1){
          ttime2[i] <- ttime2[i]+max2+delta_t2
          ttime1[i] <- ttime1[i-1]
        }
      }
      for (j in 2:length(full_t2)){
        if (abs(t2_star-full_t2[j])>2.1){
          full_t2[j] <- full_t2[j]+max2+delta_t2
          full_t1[j] <- full_t1[j-1]
        }
      }
    }
    else {
      # Adjusting for low peak
      for (i in (length(ttime2)-1):1){
        if (abs(t2_star-ttime2[i])>2.1){
          ttime2[i] <- ttime2[i]-max2-delta_t2
          ttime1[i] <- ttime1[i+1]
        }
      }
      for (j in (length(full_t2)-1):1){
        if (abs(t2_star-full_t2[j])>2.1){
          full_t2[j] <- full_t2[j]-max2-delta_t2
          full_t1[j] <- full_t1[j+1]
        }
      }
    }
  }

  # Fitting only the peak of the data to the Gaussian
  fitting <- nls.multstart::nls_multstart(y~gauss2(a,b1,b2,c1,c2,ttime1,ttime2),
                                          iter=10,
                                          start_lower=list(a=height-1,b1=peakcoord[1]-.2,
                                              b2=peakcoord[2]-.1, c1=sig1/3-.1, c2=sig2/3-.1),
                                          start_upper=list(a=height+1,b1=peakcoord[1]+.2,
                                              b2=peakcoord[2]+.1, c1=sig1/3+.1, c2=sig2/3+.1),
                                          algorithm= "port",
                                          control=stats::nls.control(maxiter=100,
                                                                     warnOnly=TRUE),
                                          supp_errors='Y')
  alpha <- fitting$m$getPars()[1]
  beta1 <- fitting$m$getPars()[2]
  beta2 <- fitting$m$getPars()[3]
  cat1 <- abs(fitting$m$getPars()[4])
  cat2 <- abs(fitting$m$getPars()[5])
  params <- c(alpha, beta1, beta2, cat1, cat2)
  gaussfit2 <- gauss2(alpha,beta1,beta2,cat1,cat2,
                     full_t1,full_t2)
  gauss2_r <- as.data.frame(cbind(frame$'Overall Time Index'[time], full_t1,
                                  full_t2, gaussfit2))
  colnames(gauss2_r) <- c('t','time1', 'time2', 'gaussfit')
  volume <- abs(2*pi*alpha*cat1*cat2)[[1]]
  gauss2_return <- list(gauss2_r,params,volume)
  names(gauss2_return) <- c("Gauss_df","params","volume")
  return(gauss2_return)
}


#' @title  Plots a peak with the fitted Gaussian curve.
#'
#' @description  `plot_gauss` Plots a peak with the fitted Gaussian curve.
#'
#' @details This function plots the points around the peak in blue dots, with a
#' line plot of the Gaussian curve fit to the peak data in red, using
#' \code{\link[ggplot2]{ggplot}} from ggplot2 package
#' \insertCite{ggplot2}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param TIC_df a \emph{data.frame} object. Data frame with 4 columns
#' (Overall Time Index, RT1, RT2, TIC), ideally the output from create_df(), or
#' the first data frame returned from extract_data(), $TIC_df.
#' @param gauss_return a \emph{data.frame} object. The output from guass_fit().
#' A data frame with two columns, (time, guassfit), the time values around the
#' peak, and the intensity values fitted to the optimal Gaussian curve.
#' @param title a \emph{string} object. Title placed at the top of the plot.
#'
#' @return A \emph{ggplot} object. A plot of points around the peak with a line
#' plot of the Gaussian curve fit to the peak data.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' gaussfit <- gauss_fit(frame$TIC_df, peakcoord=c(peaks$'X'[1], peaks$'Y'[1]))
#' message(paste('Area under curve =',gaussfit[[3]], 'u^2'))
#' plot_gauss(frame$TIC_df, gaussfit[[1]])
#'
#'@export
plot_gauss <- function(TIC_df, gauss_return, title="Peak fit to Gaussian"){
  frame <- TIC_df
  time <- c(gauss_return$'time'[1],
            gauss_return$'time'[length(gauss_return$'time')])
  t1 <- which(frame$'Overall Time Index'== time[1])
  t2 <- which(frame$'Overall Time Index'== time[2])
  data <- as.data.frame(cbind(gauss_return$'time',
                              frame$'TIC'[t1:t2],
                              gauss_return$'gaussfit'))
  colnames(data) <- c('time', 'TIC', 'gaussfit')

  x <- data$'time'
  y1 <- data$'TIC'
  y2 <- data$'gaussfit'
  fig <- ggplot2::ggplot(data) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y1, color = 'data')) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = y2, color='Gaussian fit')) +
    ggplot2::labs(title= title, x='time', y= 'intensity', color='Legend') +
    ggplot2::theme_light() +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
    ggplot2::scale_colour_manual(values = c('data'="blue", 'Gaussian fit'="red"),
                        guide = ggplot2::guide_legend(override.aes = list(
                          linetype = c("blank", "solid"),
                          shape = c(19, NA))))

  return(fig)
}




#' @title  Plots a 3D peak with the fitted Gaussian curve.
#'
#' @description  `plot_gauss2` Plots a 3D peak with the fitted Gaussian curve.
#'
#' @details This function plots the points around the peak with a
#' contour plot of the Gaussian curve fit to the peak data, using
#' \code{\link[ggplot2]{ggplot}} from ggplot2 package
#' \insertCite{ggplot2}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param TIC_df a \emph{data.frame} object. Data frame with 4 columns
#' (Overall Time Index, RT1, RT2, TIC), ideally the output from create_df(), or
#' the first data frame returned from extract_data(), $TIC_df.
#' @param gauss2_return a \emph{data.frame} object. The output from guass_fit().
#' A data frame with two columns, (time, guassfit), the time values around the
#' peak, and the intensity values fitted to the optimal Gaussian curve.
#' @param title a \emph{string} object. Title placed at the top of the plot.
#'
#' @return A \emph{ggplot} object. A contour plot of the Gaussian curve fit to
#' the peak data.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' gaussfit2 <- gauss2_fit(frame$TIC_df, peakcoord=c(peaks$'X'[1], peaks$'Y'[1]))
#' message(paste('Volume under curve =',gaussfit2[[3]],'u^3'))
#' plot_gauss2(frame$TIC_df, gaussfit2[[1]])
#'
#'@export
plot_gauss2 <- function(TIC_df, gauss2_return, title="Peak fit to Gaussian"){
  x <- gauss2_return$'time1'
  y <- gauss2_return$'time2'
  z <- log(gauss2_return$'gaussfit',base=10)

  fig <- ggplot2::ggplot(gauss2_return) +
    ggplot2::geom_tile(ggplot2::aes(x=x, y=y, fill=z)) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title= title, x="retention time 1", y= "retention time 2",
                  fill = "Log Intensity") +
    ggplot2::theme_light() +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(ticks = FALSE))
  fig
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
