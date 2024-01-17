#' @title   Targeted Analysis
#'
#' @description  `targeted` performs targeted analysis for a batch of data
#' files, for a list of masses of interest.
#'
#' @details This function performs targeted analysis for a batch of data
#' files, for a list of masses of interest.
#'
#' @param data_list a \emph{list} object. Data extracted from each cdf file,
#' ideally the output from extract_data().
#' @param MOIs a \emph{vector} object. A vector containing a list of all masses
#' of interest to be investigated.
#' @param RTs a \emph{vector} object. An optional vector containing a list of
#' retention times of interest for the listed masses of interest. Default values
#' if left empty will be at the retention time of the highest intensity for the
#' corresponding mass.
#' @param window_size a \emph{vector} object. An optional vector containing a
#' list of window sizes corresponding to the retention times.  Window will be
#' defined by (RT-window_size, RT+window_size). Default if left
#' empty will be 0.1.
#' @param tolerance a \emph{float} object. The tolerance allowed for the MOI.
#' Default is 0.005.
#' @param images a \emph{boolean} object. An optional input. If TRUE, all images
#' of the found peaks will be displayed. Default is FALSE, no images will be
#' displayed.
#'
#' @return a \emph{data.frame} object. A data frame containing the areas of the
#' peaks for the indicated MOIs and list of files.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' file2 <- system.file("extdata","sample2.cdf",package="gcxgclab")
#' file3 <- system.file("extdata","sample3.cdf",package="gcxgclab")
#' frame1 <- extract_data(file1,mod_t=.5)
#' frame2 <- extract_data(file2,mod_t=.5)
#' frame3 <- extract_data(file3,mod_t=.5)
#' targeted(list(frame1,frame2,frame3),MOIs = c(92.1397, 93.07058),
#' RTs = c(6.930, 48.594), images=TRUE)
#'
#' @export
targeted <- function(data_list, MOIs, RTs=c(), window_size=c(),
                     tolerance=0.005,images=FALSE){
  ms_list <- list()
  for (i in 1:length(data_list)){
    ms_list <- append(ms_list,data_list[[i]][2])
  }
  N <- length(ms_list)
  M <- length(MOIs)
  ms_data <- ms_list
  names <- c()
  for (x in 1:length(ms_list)){
    names <- append(names,paste("Sample",x))
  }
  names(ms_data)<- names
  rnames <- c("RT","----------",names)

  if (length(RTs)>M){
    stop("The number of retention times listed is greater than the number of masses of interest. If listing multiple retention times for one mass, repeat mass value in MOIs input.")
  }
  if (length(RTs)>M){
    stop("The number of window sizes is greater than the number of masses of interest.")
  }

  ##
  eic <- list()
  results1 <- data.frame(matrix(nrow=(N+2),ncol=(M+1)))
  colnames(results1) <- c("|",MOIs)
  rownames(results1) <- rnames
  results1[2,]<- rep("----------",(M+1))
  results1[,1] <- rep("|",(N+2))

  for (z in 1:(N+2)){

  }

  #results2 <- results1
  for (i in 1:N){
    eics <- batch_eic(data_list[[i]],MOIs,tolerance)
    for (j in 1:M){
      forblc <- as.data.frame(cbind(eics[[j]]$time_array,eics[[j]]$RT1,eics[[j]]$RT2,eics[[j]]$intensity_values))
      colnames(forblc) <- c('Overall Time Index','RT1','RT2','TIC')
      if (length(forblc$`Overall Time Index`)>length(unique(forblc$`Overall Time Index`))){
        sums <- c()
        t1 <- c()
        t2 <- c()
        for (k in 1:length(unique(forblc$'Overall Time Index'))){
          sums <- append(sums, sum(forblc$TIC[forblc$`Overall Time Index`==unique(forblc$`Overall Time Index`)[k]]))
          loc2 <- which(forblc$`Overall Time Index`==unique(forblc$`Overall Time Index`[i]))[1]
          t1 <- append(t1,forblc$RT1[loc2])
          t2 <- append(t2,forblc$RT2[loc2])
        }
        forblc <- as.data.frame(cbind(unique(forblc$`Overall Time Index`),t1,t2,sums))
        colnames(forblc) <- c('Overall Time Index','RT1','RT2','TIC')
      }

      if (length(forblc$TIC)>500){
        forblc <- list(forblc)
        names(forblc) <- c('TIC_df')
        forblc <- bl_corr(forblc,.25)$TIC_df
      }
      else{
        forblc$TIC <- forblc$TIC-min(forblc$TIC)
      }

      if (is.null(window_size[j])){
        ws <- .1
      }
      else if (is.na(window_size[j])){
        ws <- .1
      }
      else{
        ws <- window_size[j]
      }
      if (is.null(RTs[j])){
        loc <- which(forblc$TIC==max(forblc$TIC))
      }
      else if (is.na(RTs[j])){
        loc <- which(forblc$TIC==max(forblc$TIC))
      }
      else{
        loc1 <- which.min(abs(forblc$'Overall Time Index'-RTs[j]+ws))
        loc2 <- which.min(abs(forblc$'Overall Time Index'-RTs[j]-ws))
        loc <- which(forblc$TIC==max(forblc$TIC[c(loc1:loc2)]))
        loc <- intersect(loc,c(loc1:loc2))[1]
      }
      pt1 <- which.min(abs(forblc$`Overall Time Index`[loc]-forblc$`Overall Time Index`-ws))
      pt2 <- which.min(abs(forblc$`Overall Time Index`[loc]-forblc$`Overall Time Index`+ws))
      if (pt1==loc){
        pt1 <- loc-1
      }
      if (pt2==loc){
        pt2 <- loc+1
      }
      new <- forblc[pt1:pt2,]
      # L <- length(pt1:pt2)
      # LL <- round(L/2)
      # if (L>50){
      #   newnew <- new[c(1:(LL-20),(LL+20):L),]
      # }
      # else if (L>10){
      #   newnew <- new[c(1:(LL-3),(LL+3):L),]
      # }
      # else {
      #   newnew <- new[c(1:(LL-1),(LL+1):L),]
      # }
      flr <- max(1,(pt1-50))
      clg <- min(length(forblc$TIC),(pt2+50))
      newnew <- forblc$TIC[c(flr:pt1,pt2:clg)]
      a <- max(new$TIC,na.rm=TRUE)
      b1 <- forblc$`Overall Time Index`[loc]
      ##
      thr <- mean(newnew,na.rm=TRUE)+2*stats::sd(newnew,na.rm=TRUE)
      nd <- FALSE
      eicavg <-mean(forblc$TIC,na.rm=TRUE)
      if (a<thr | a<100 | a<2*eicavg){
        a <- 0
        cat <- 0.01
        beta <- b1
        nd <- TRUE
        message(paste("Compound for mass", MOIs[j], "at retention time", RTs[j],
                      "was not detected in", rnames[i+2],
                      "due to low intensity values."))
      }
      else{
        fitting <- NULL
        y_var <- new$TIC
        t_var <- new$`Overall Time Index`
        fitting_list <- list()
        suppressWarnings({
          fitting <- nls.multstart::nls_multstart(y_var~gauss(a,b,c,t_var),
                                                  iter=100,
                                            start_lower=list(b=b1-0.02,c=0.001),
                                            start_upper=list(b=b1+0.02,c=0.01),
                                                  supp_errors = 'Y',
                                            upper = c(Inf,0.01))
        })
        if (!(is.null(fitting))){
          beta <- fitting$m$getPars()[1]
          cat <- abs(fitting$m$getPars()[2])
          if (cat>100){
            nd <- TRUE
            message(paste("Compound for mass", MOIs[j], "at retention time", RTs[j],
                          "was not detected in", rnames[i+2],
                          "due to non-fitting Gaussian peak."))
            a <- 0
            cat <- 0.01
          }
        }
        else{
          nd <- TRUE
          message(paste("Compound for mass", MOIs[j], "at retention time", RTs[j],
                        "was not detected in", rnames[i+2],
                        "due to non-fitting Gaussian peak."))
          a <- 0
          cat <- 0.01
        }
      }
      if (abs(sqrt(2*pi)*a*cat)<.5 & !nd){
        a <- 0
        cat <- 0.01
        beta <- b1
        nd <- TRUE
        message(paste("Compound for mass", MOIs[j], "at retention time", RTs[j],
                      "was not detected in", rnames[i+2],
                      "due to low response.")) #A
      }
      if (a!=0){
        window <- cat*4
        new2 <- dplyr::filter(forblc, forblc$`Overall Time Index`>= beta-window)
        new2 <- dplyr::filter(new2, new2$`Overall Time Index`<= beta+window)
        # if (length(new2$TIC)<2){
        #   a <- 0
        #   cat <- 0.01
        #   beta <- b1
        #   nd <- TRUE
        #   message(paste("Compound for mass", MOIs[j], "at retention time", RTs[j],
        #                 "was not detected in", rnames[i+1],
        #                 "due low response B."))
        # }
        # else{
        if(TRUE){
          gaussfit <- gauss(a,beta,cat,new$`Overall Time Index`)
          gauss_r <- as.data.frame(cbind(new$`Overall Time Index`, gaussfit))
          colnames(gauss_r) <- c('time', 'gaussfit')
          area <- abs(sqrt(2*pi)*a*cat)
          x1 <- new2$`Overall Time Index`
          y1 <- new2$TIC
          x2 <- gauss_r$time
          y2 <- gaussfit

          title <- paste("Full EIC", MOIs[j],names(ms_data)[i])
          x3 <- forblc$'Overall Time Index'
          y3 <- forblc$TIC
          if (images){
          suppressMessages({
            fig <- ggplot2::ggplot(forblc) +
              ggplot2::geom_point(ggplot2::aes(x = x3, y = y3, color = 'data')) +
              ggplot2::geom_line(data=gauss_r,ggplot2::aes(x = x2, y = y2, color='Gaussian fit')) +
              ggplot2::labs(title= title, x='time', y= 'intensity', color='Legend') +
              ggplot2::theme_light() +
              ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
              ggplot2::scale_colour_manual(values = c('data'="blue", 'Gaussian fit'="red"),
                                           guide = ggplot2::guide_legend(override.aes = list(
                                             linetype = c("blank", "solid"),
                                             shape = c(19,  NA))))
            print(fig)
          })
          }

          x1 <- new2$`Overall Time Index`
          y1 <- new2$TIC
          newt <- seq(new2$'Overall Time Index'[1],
                      new2$'Overall Time Index'[length(new2$'Overall Time Index')],
                      (new2$'Overall Time Index'[length(new2$'Overall Time Index')]
                       -new2$'Overall Time Index'[1])/100)

          gaussfit2 <- gauss(a,beta,cat,newt)
          gauss_r2 <- as.data.frame(cbind(newt, gaussfit2))
          colnames(gauss_r2) <- c('time', 'gaussfit')
          x2 <- newt
          y2 <- gaussfit2

          title <- paste("EIC",MOIs[j], names(ms_data)[i])
          if (images){
          suppressMessages({
            fig2 <- ggplot2::ggplot(new2) +
            ggplot2::geom_point(ggplot2::aes(x = x1, y = y1, color = 'data')) +
            ggplot2::geom_line(data=gauss_r2,ggplot2::aes(x = x2, y = y2, color='Gaussian fit')) +
            ggplot2::labs(title= title, x='time', y= 'intensity', color='Legend') +
            ggplot2::theme_light() +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
            ggplot2::scale_colour_manual(values = c('data'="blue", 'Gaussian fit'="red"),
                                         guide = ggplot2::guide_legend(override.aes = list(
                                           linetype = c("blank", "solid"),
                                           shape = c(19, NA))))
          print(fig2)
          })
          }
        }
      }
      if (a==0){
        title <- paste("Full EIC", MOIs[j],names(ms_data)[i])
        x3 <- forblc$'Overall Time Index'
        y3 <- forblc$TIC

        if (images){
        suppressMessages({
          fig <- ggplot2::ggplot(forblc) +
            ggplot2::geom_point(ggplot2::aes(x = x3, y = y3, color = 'data')) +
            ggplot2::labs(title= title, x='time', y= 'intensity', color='Legend') +
            ggplot2::theme_light() +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+
            ggplot2::scale_colour_manual(values = c('data'="blue"),
                                         guide = ggplot2::guide_legend(override.aes = list(
                                           linetype = c("blank"),
                                           shape = c(19))))
          print(fig)
        })
        }
      }
      if (nd){
        area <- 'ND'
      }
      if (is.null(RTs[j])){
        results1[1,(j+1)]<-b1
      }
      else if (is.na(RTs[j])){
        results1[1,(j+1)]<-b1
      }
      else{
        results1[1,(j+1)]<-RTs[j]
      }
      if(area !='ND'){
        area <- round(as.numeric(area),6)
      }
      results1[(i+2),(j+1)] <- area
    }
  }
  rn <- c("MOI",rownames(results1))
  results1 <- as.data.frame(rbind(colnames(results1),results1))
  rownames(results1)<- rn
  colnames(results1)<- rep(" ", length(results1[1,]))
  #results <- list(results1,results2)
  #names(results)<-c('1D Areas','2D Volumes')
  return(results1)
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
