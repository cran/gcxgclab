#' @title   Reference Batch Align
#'
#' @description  `align` aligns peaks from samples to a reference sample's
#' peaks.
#'
#' @details This function aligns the peaks from any number of samples. Peaks are
#' aligned to the retention times of the first peak. If aligning to a reference
#' or standard sample, this should be the first in the lists for data frames and
#' for the mass data. The function comp_peaks() is used to find the
#' corresponding peaks. This function  will return a new list of TIC data frames
#' and a list of mass data. The first sample's data is unchanged, used as the
#' reference. Then a TIC data frame and mass data for each of the given samples
#' containing the peaks and time coordinates of the aligned peaks. The time
#' coordinates are aligned to the first sample's peaks, the peak height and MS
#' is unchanged.
#'
#' @param data_list a \emph{list} object. Data extracted from each cdf file,
#' ideally the output from extract_data().
#' @param THR a \emph{float} object. Threshold for peak intensity. Should be a
#' number between the baseline value and the highest peak intensity. Default
#' is THR = 100000.
#' @param use_ref_peak a \emph{boolean} object. Determines if an initial shift
#' to a given reference peak, default is toluene, should be done before aligning
#' all other peaks above given threshold THR. Default is TRUE.
#' @param ref_peak a \emph{float} object. The m/z value of the reference peak
#' for optional initial shift. Default is 92.1397 (toluene).
#'
#' @return A \emph{list} object. List of aligned data from each cdf file and a
#' list of peaks that were aligned for each file.
#'
#' @examples
#' file1 <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' file2 <- system.file("extdata","sample2.cdf",package="gcxgclab")
#' file3 <- system.file("extdata","sample3.cdf",package="gcxgclab")
#' frame1 <- extract_data(file1,mod_t=.5)
#' frame2 <- extract_data(file2,mod_t=.5)
#' frame3 <- extract_data(file3,mod_t=.5)
#' aligned <- align(list(frame1,frame2,frame3))
#' plot_peak(aligned$Peaks$S1,aligned$S1,title="Reference Sample 1")
#' plot_peak(aligned$Peaks$S2,aligned$S2,title="Aligned Sample 2")
#' plot_peak(aligned$Peaks$S3,aligned$S3,title="Aligned Sample 3")
#'
#' @export
align <- function(data_list, THR=100000, use_ref_peak=TRUE, ref_peak=92.1397){

  ref_df <- data_list[[1]]$TIC_df
  ref_df$'Overall Time Index' <- ref_df$'Overall Time Index'/60
  data_list[[1]]$TIC_df$`Overall Time Index` <- data_list[[1]]$TIC_df$`Overall Time Index`/60
  ref_ms <- data_list[[1]]$MS_df

  df_list <- vector(mode="list", length=length(data_list)-1)
  ms_list <- vector(mode="list", length=length(data_list)-1)

  for (i in 1:(length(data_list)-1)){

    df_list[[i]] <- data_list[[i+1]]$TIC_df


    df_list[[i]]$`Overall Time Index` <- df_list[[i]]$`Overall Time Index`/60
    data_list[[i+1]]$TIC_df$`Overall Time Index` <- data_list[[i+1]]$TIC_df$`Overall Time Index`/60
    ms_list[[i]] <- data_list[[i+1]]$MS_df
  }

  ref_peaks <- thr_peaks(ref_df,THR)
  ref_y <- length(which(ref_df$'RT1'==ref_df$'RT1'[1]))
  ref_x <- length(ref_df$'RT1')/ref_y
  if (length(df_list)==0 | length(ms_list)==0){
    stop('Must provide at least one sample to be aligned to the reference.')
  }
  if (length(df_list)!=length(ms_list)){
    stop('Number of TIC data frames does not match the number of MS data frames.')
  }
  aligned <- vector("list", (length(df_list)+2))

  naming <- paste0('S', seq(1,length(df_list)+1, 1))

  aligned[[1]] <- list(ref_df,ref_ms)
  names(aligned[[1]])<- c('TIC_df','MS_df')
  names(aligned) <- c(naming, 'Peaks')
  aligned$Peaks <- rep(list(ref_peaks), length(df_list)+1)

  if (use_ref_peak){
    tol1 <- find_eic(data_list[[i]],ref_peak,tolerance=10^(-3))
    lt1 <- which(tol1$intensity_values==max(tol1$intensity_values))
    ltt1 <- tol1$time_array[lt1]
  }
  else {
    ref_peak <- 0
  }

  for (i in 1:length(df_list)){

    RT1 <- df_list[[i]]$RT1
    RT2 <- df_list[[i]]$RT2
    new_ms <- ms_list[[i]]
    message(paste("Aligning peaks of Sample",(i+1),"to the reference, Sample 1."))
    pb <- utils::txtProgressBar(min = 0, max = 100, style = 3, char = "=")

    df_y <- length(which(df_list[[i]]$'RT1'==df_list[[i]]$'RT1'[1]))
    df_x <- length(df_list[[i]]$'RT1')/df_y
    if (abs(ref_x-df_x)>1 | abs(ref_y-df_y)>1){
      stop(paste('Peak',i, 'cannot be aligned to the reference.
                           Dimensions do not match.'))
    }

    if (use_ref_peak){
      tol2 <- find_eic(data_list[[i+1]],ref_peak,tolerance=10^(-3))
      lt2 <- which(tol2$intensity_values==max(tol2$intensity_values))
      ltt2 <- tol2$time_array[lt2]
      if (ltt1>8.8 & ltt2>8.8 & ltt1<9.4 & ltt2<9.4){
        df_list[[i]]<- phase_shift(data_list[[i+1]],(ltt1-ltt2))[[1]]
      }
    }
    else{
      df_list[[i]] <- data_list[[i+1]][[1]]
    }


    peaks <- thr_peaks(df_list[[i]],THR)
    cp <- comp_peaks(ref_peaks,peaks)
    todel <- c()

    for (j in 1:length(cp$Tref)){

      utils::setTxtProgressBar(pb, round(10/length(cp$Tref)*j))
      ms_r <- find_ms(data_list[[1]],(cp$Tref[j]*60),tolerance=10^-3)
      ms_a <- find_ms(data_list[[i+1]],(cp$Tal[j]*60),tolerance =10^-3)
      ms_r2 <- ms_r
      ms_a2 <- ms_a
      ms_r2$MZ <- round(ms_r$MZ)
      ms_a2$MZ <- round(ms_a$MZ)
      ms_r2 <- dplyr::arrange(ms_r2, ms_r2$MZ)
      ms_a2 <- dplyr::arrange(ms_a2, ms_a2$MZ)
      k <- 2
      while (k <= length(ms_r2$MZ)){
        if (ms_r2$MZ[k]==ms_r2$MZ[k-1]){
          ms_r2$Int[k-1]<- ms_r2$Int[k-1]+ms_r2$Int[k]
          ms_r2<-ms_r2[-k,]
        }
        else {
          k <- k+1
        }
      }
      k <- 2
      while (k <= length(ms_a2$MZ)){
        if (ms_a2$MZ[k]==ms_a2$MZ[k-1]){
          ms_a2$Int[k-1]<- ms_a2$Int[k-1]+ms_a2$Int[k]
          ms_a2<-ms_a2[-k,]
        }
        else {
          k <- k+1
        }
      }
      ms_r2$Int <- ms_r2$Int/max(ms_r2$Int)
      ms_a2$Int <- ms_a2$Int/max(ms_a2$Int)
      ms_r2 <- ms_r2[ms_r2$Int>.1,]
      ms_a2 <- ms_a2[ms_a2$Int>.1,]

      sum <- 0

      for (x in 1:length(ms_a2$MZ)){
        loc <- which(ms_r2$MZ==ms_a2$MZ[x])
        if (length(loc)>0){
          sum <- sum+ abs(ms_r2$Int[loc]-ms_a2$Int[x])
        }
        else {
          sum <- sum+ms_a2$Int[x]
        }
      }

      for (x in 1:length(ms_r2$MZ)){
        loc <- which(ms_a2$MZ==ms_r2$MZ[x])
        if (length(loc)>0){
          sum <- sum+ abs(ms_a2$Int[loc]-ms_r2$Int[x])
        }
        else {
          sum <- sum+ms_r2$Int[x]
        }
      }

      match_per <- 1-sum/(sum(ms_r2$Int)+sum(ms_a2$Int))
      if (match_per<.7){
        todel <- append(todel,j)
      }
    }

    # remove peaks which do not have matching MS
    if (length(todel)>0){
      cp2 <- cp[-todel,]
    } else{
      cp2 <- cp
    }
    if (length(cp2)==0){
      stop("No peaks match with sufficient spacial and MS confidence.")
    }
    frame <- df_list[[i]]
    for (x in 1:length(cp2[,1])){
      utils::setTxtProgressBar(pb,round(10+89/length(cp2[,1])*x))
      j <- length(cp2[,1])-x+1
      if (abs(cp2$Tref[j]-cp2$Tal[j])>0.0001){
        # Shift TIC
        suppressWarnings(
          g <- gauss_fit(df_list[[i]], c(cp2$Xal[j],cp2$Yal[j]))
        )
        cat <- abs(g[[2]][3])
        l1 <- which.min(abs(cp2$Tref[j]-4*cat-df_list[[i]]$`Overall Time Index`))
        l2 <- which.min(abs(cp2$Tref[j]+4*cat-df_list[[i]]$`Overall Time Index`))
        if (is.na(l2)){
          l2<- length(df_list[[i]]$TIC)
        }
        cut <- df_list[[i]]$TIC[c(l1:l2)]
        ldiff1 <- which.min(abs(cp2$Tref[j]-df_list[[i]]$`Overall Time Index`))
        ldiff2 <- which.min(abs(cp2$Tal[j]-df_list[[i]]$`Overall Time Index`))
        ldiff <- ldiff1-ldiff2

        xdim <- length(unique(frame$RT1))
        ydim <- length(frame$RT2)/xdim
        sgn <- ldiff/abs(ldiff)*abs(round(ldiff/ydim))
        if (abs(ldiff)>(ydim-10)){
          if ((l2+sgn*ydim)>length(frame$TIC)){
            len1 <- length(c(l1+sgn*ydim):length(frame$TIC))
            df_list[[i]]$TIC[c((l1+sgn*ydim):length(frame$TIC))]<- cut[1:len1]
            ldiff <- ldiff-sgn*ydim
          } else{
            df_list[[i]]$TIC[c((l1+sgn*ydim):(l2+sgn*ydim))]<- cut
            ldiff <- ldiff-sgn*ydim
          }
        }
        if (ldiff>0){
          frame$TIC[c(l1:(l1+ldiff-1))]<-df_list[[i]]$TIC[c((l2+1):(l2+ldiff))]
          if ((l2+ldiff)>length(frame$TIC)){
            len1 <- length(c(l1+ldiff):length(frame$TIC))
            frame$TIC[c((l1+ldiff):length(frame$TIC))]<-cut[1:len1]
          } else{
            frame$TIC[c((l1+ldiff):(l2+ldiff))]<-cut
          }
        } else if (ldiff<0){
          frame$TIC[c((l2+ldiff+1):l2)]<-df_list[[i]]$TIC[c((l1+ldiff):(l1-1))]
          frame$TIC[c((l1+ldiff):(l2+ldiff))]<-cut
        }

        # Shift MS
        t1 <- df_list[[i]]$`Overall Time Index`[l1]
        t2 <- df_list[[i]]$`Overall Time Index`[l2]
        loc1 <- which.min(abs(ms_list[[i]]$time_array-t1))
        loc2 <- which.min(abs(ms_list[[i]]$time_array-t2))
        L1 <- which(ms_list[[i]]$time_array==ms_list[[i]]$time_array[loc1])[1]
        L2 <- which(ms_list[[i]]$time_array==ms_list[[i]]$time_array[loc2])
        L3 <- L2[length(L2)]

        tdiff <- cp2$Tref[j]-cp2$Tal[j]
        times <- unique(ms_list[[i]][L1:L3,c(1,2,3)])
        new_t <- c()
        ldiff <- ldiff1-ldiff2

        suppressWarnings({
        temp <- data.frame(xx = which.min(abs(ms_list[[i]]$time_array-times$time_array-tdiff)))
        })
        temp$tloc <- rownames(temp)
        tloc <- temp[which(!is.na(temp$xx)),]
        new_t <- rbind(new_t,ms_list[[i]][tloc$tloc,c(1,2,3)])

        if (abs(ldiff)>(ydim-10)){
          new_t2 <- c()
          rt1cut <- ms_list[[i]]$RT1[L1]
          rt1s <- unique(ms_list[[i]]$RT1[sgn*ms_list[[i]]$RT1>=sgn*rt1cut])
          if (sgn>0){
            rt1s <- rt1s[1]
          } else if (sgn<0){
            rt1s <- rt1s[length(rt1s)]
          }
          LL1 <- which(ms_list[[i]]$RT1==rt1s & ms_list[[i]]$RT2==ms_list[[i]]$RT2[L1])[1]
          LL2 <- which(ms_list[[i]]$RT1==rt1s & ms_list[[i]]$RT2==ms_list[[i]]$RT2[L3])
          LL3 <- LL2[length(LL2)]
          new_t2 <- ms_list[[i]][LL1:LL3,c(1,2,3)]
          new_t2 <- unique(new_t2)


          for (z in 1:length(times$time_array)){
            w1 <- which(times$time_array[z]==ms_list[[i]]$time_array)
            w2 <- which(new_t2$time_array[z]==ms_list[[i]]$time_array)
            ms_list[[i]][w1,c(1,2,3)] <- new_t2[z,]
            ms_list[[i]][w2,c(1,2,3)] <- times[z,]
          }

          ldiff <- ldiff-sgn*ydim
        }

        if (ldiff>0){
          # shift high times down
          for (z in c((length(times$time_array)+1-ldiff):length(times$time_array))){
            lloc <- which(new_t$time_array[z]==ms_list[[i]]$time_array)
            zz <- z-length(times$time_array)+ldiff
            new_ms[lloc,c(1,2,3)] <- times[zz,]
          }
        } else if (ldiff<0){
          # shift low times up
          for (z in c(1:(-ldiff))){
            lloc <- which(new_t$time_array[z]==ms_list[[i]]$time_array)
            zz <- z+length(times$time_array)+ldiff
            new_ms[lloc,c(1,2,3)] <- times[zz,]
          }
        }

        hold <- data.frame(y = match(ms_list[[i]]$time_array, times$time_array))
        hold$loca <- rownames(hold)
        loca <- hold[which(!is.na(hold$y)),]

        new_ms[loca$loca, c(1,2,3)] <- new_t[loca$y,]
      }
    }
    utils::setTxtProgressBar(pb,100)
    cp3 <- cp2[,c(1,2,3,8)]
    colnames(cp3) <- c('T','X','Y','Peak')
    frame$`Overall Time Index` <- frame$`Overall Time Index`*60
    frame <- dplyr::arrange(frame, frame$'Overall Time Index')
    frame$RT1 <- RT1
    frame$RT2 <- RT2
    aligned[[i+1]] <- list(frame,new_ms)
    names(aligned[[i+1]]) <- c('TIC_df','MS_df')
    cp3$T <- cp3$T*60
    aligned$Peaks[[i+1]] <- cp3
    close(pb)
  } # end loop over i, samples being aligned

  names(aligned$Peaks) <- naming

  return(aligned)
}


#' @title   Compare Peaks
#'
#' @description  `comp_peaks` compares peaks of two samples.
#'
#' @details This function find compares the peaks from two samples and
#' correlates the peaks by determining the peaks closest to each other in the
#' two samples, within a certain reasonable distance. Then returns a data frame
#' with a list of the correlated peaks including each of their time coordinates.
#'
#' @param ref_peaks a \emph{data.frame} object. A data frame with 4 columns
#' (Time, X, Y, Peak), ideally the output from either top_peaks() or
#' thr_peaks().
#' @param al_peaks a \emph{data.frame} object. A data frame with 4 columns
#' (Time, X, Y, Peak), ideally the output from either top_peaks() or
#' thr_peaks().
#'
#' @return A \emph{data.frame} object. A data frame with 8 columns containing
#' the matched peaks from the two samples, with the time, x, y, and peak values
#' for each.
#'
#' @export
comp_peaks <- function(ref_peaks, al_peaks){
  N <- length(ref_peaks$'T')
  M <- length(al_peaks$'T')
  # Finding peak in second list with min dist to first list peak
  keep_2 <- c()
  keep_1 <- c()
  for (i in 1:N){
    min_d <- 1
    coord <- 0
    for (j in 1:M){
      distance <- 0.01*(ref_peaks$'X'[i]-al_peaks$'X'[j])^2+
        3*(ref_peaks$'Y'[i] -al_peaks$'Y'[j])^2
      min_d <-  min(min_d,distance)
      if (min_d == distance){coord <- j}
    }
    if (coord != 0){keep_2 <- append(keep_2, coord)
    keep_1 <- append(keep_1, i)
    }
  }

  # Finding peak in first list with min dist to second list peak
  keep_1b <- c()
  keep_2b <- c()
  for (i in 1:M){
    min_d <- 1
    coord <- 0
    for (j in 1:N){
      distance <- 0.01*(ref_peaks$'X'[j]-al_peaks$'X'[i])^2+3*(ref_peaks$'Y'[j]
                                                               -al_peaks$'Y'[i])^2
      min_d <-  min(min_d,distance)
      if (min_d == distance){coord <- j}
    }
    if (coord != 0){
      keep_1b <- append(keep_1b, coord)
      keep_2b <- append(keep_2b, i)
    }
  }

  # Keeping only matched peaks
  ordered1b <- ref_peaks[keep_1b, ]
  ordered2b <- al_peaks[keep_2b, ]

  all_peaks_data <- as.data.frame(cbind(ordered1b$'T', ordered1b$'X',
                                        ordered1b$'Y', ordered1b$'Peak',
                                        ordered2b$'T', ordered2b$'X',
                                        ordered2b$'Y', ordered2b$'Peak'))
  colnames(all_peaks_data) <- c('Tref', 'Xref', 'Yref', 'Peakref', 'Tal', 'Xal',
                                'Yal', 'Peakal')

  # Removing double(+) matched peaks, keeping highest matched peak
  all_peaks_data <- dplyr::arrange(all_peaks_data, 'Tal', dplyr::desc("Peakref"))
  Num <- length(all_peaks_data$'Tal')
  if (Num>2){
    for (j in 2:Num){
      i <- Num-j+2
      if (all_peaks_data$'Tal'[i-1]==all_peaks_data$'Tal'[i]){
        all_peaks_data <- all_peaks_data[-i, ]
      }
    }
    all_peaks_data <- dplyr::arrange(all_peaks_data, 'Tref', dplyr::desc("Peakal"))
    Num <- length(all_peaks_data$'Tref')
    for (j in 2:Num){
      i <- Num -j+2
      if (all_peaks_data$'Tref'[i-1]==all_peaks_data$'Tref'[i]){
        all_peaks_data <- all_peaks_data[-i, ]
      }
    }
  }
  return(all_peaks_data)
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
