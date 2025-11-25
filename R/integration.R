#' @title  TIC integration
#'
#' @description  `TIC_integrate` calculates the total area under the TIC curve.
#'
#' @details This function calculates the total area under the TIC curve using
#' Simpson's Rule area approximation. Area given in 1 flattened time dimension.
#'
#' @param data a \emph{list} object. Data extracted from a cdf file,
#' ideally the output from extract_data().
#' @param start_t a \emph{float} or \emph{string} object. Value of starting time
#' for integration range. Default is 'first' which will start at the initial
#' time in the data.
#' @param end_t a \emph{float} or \emph{string} object. Value of ending time for
#' integration range. Default is 'last' which will start at the final time in
#' the data.
#'
#' @return A \emph{float} object. The calculated approximation of the area under
#' the TIC curve.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' sm_frame <- smooth(frame, lambda=10)
#' blc_frame <- bl_corr(sm_frame, gamma=0.5)
#' TIC_integrate(blc_frame)
#'
#' @export
TIC_integrate <- function(data,start_t='first',end_t='last'){
  TIC_df <- data[[1]]
  len <- length(TIC_df$TIC)
  if (start_t!='first' & length(start_t)==1 & typeof(start_t)=='double'){
    t0 <- which.min(abs(TIC_df$'Overall Time Index'-start_t))
    TIC_df <- TIC_df[-c(1:t0),]
    len <- length(TIC_df$TIC)
  }
  else if (start_t=='first'){
    t0 <- TIC_df$'Overall Time Index'[1]
  }
  else {
    stop("Start time cannot be understood.")
  }
  if (end_t!='last' & length(start_t)==1 & typeof(start_t)=='double'){
    t1 <- which.min(abs(TIC_df$'Overall Time Index'-end_t))
    TIC_df <- TIC_df[-c(t1+1:len),]
    len <- length(TIC_df$TIC)
  }
  else if (end_t=='last'){
    t1 <- TIC_df$'Overall Time Index'[len]
  }
  else {
    stop("End time cannot be understood.")
  }

  z <- len%%2
  dx <- (TIC_df$`Overall Time Index`[len]-TIC_df$`Overall Time Index`[1])/(len-1)
  if (z==0){
    add <- TIC_df$TIC[1]*dx
    TIC_df <- TIC_df[-1,]
    len <- len-1
    dx <- (TIC_df$`Overall Time Index`[len]-TIC_df$`Overall Time Index`[1])/(len-1)
  }
  else{
    add <- 0
  }
  v <- c(1,4,rep(c(2,4),(len-3)/2),1)
  S <- sum(TIC_df$TIC*v*dx/3)+add
  return(S) #units intensity*sec
}
