#' @title  Creates list of NIST data
#'
#' @description  `nist_list` creates a list of the data from the NIST MS
#' database.
#'
#' @details This function takes the MSP file containing the data from the NIST
#' MS Library database and creates a list of string vectors for each compound in
#' the database.
#'
#' @param nistfile a \emph{string} object, the file name or path of the MSP file
#' for the NIST MS Library database.
#' @param ... additional optional \emph{string} objects, the file names or paths
#' of the MSP file for the NIST MS Library if the data base is broken into
#' multiple files.
#'
#' @return nistlist, a \emph{list} object, a list of string vectors for each
#' compound in the database.
#'
#' @export
nist_list <- function(nistfile,...){
  files <- c(nistfile)
  files2 <- c(...)
  files <- append(files,files2)
  N <- length(files)
  nistlist <- list()
  prog <- 0
  for (i in 1:N){
    nistlist_a <- list()
    nist <- readLines(files[1])
    skip <- which(nist=="")
    if (i==1){
      M <- length(skip)
      pb <- utils::txtProgressBar(min = 0, max = N*M, style = 3, char = "=")
    }
    nistlist_a <- append(nistlist_a,list(nist[1:(skip[1]-1)]))
    for (i in 1:(length(skip)-1)){
      nistlist_a <- append(nistlist_a, list(nist[(skip[i]+1):(skip[i+1]-1)]))
      prog <- prog + round(M/length(skip))
      utils::setTxtProgressBar(pb, prog)
    }
    nistlist <- append(nistlist,nistlist_a)
  }
  close(pb)
  return(nistlist)
}


#' @title  Compares MS to NIST MS database
#'
#' @description  `comp_nist` compares the MS data from a peak to the NIST MS
#' database.
#'
#' @details This function takes the MS data from an intensity peak in a sample
#' and compares it to the NIST MS Library database and determines the compound
#' which is the best match to the MS data.
#'
#' @param nistlist a \emph{list} object, a list of compound MS data from the
#' NIST MS Library database, ideally the output of nist_list().
#' @param ms a \emph{data.frame} object, a data frame of the mass values and the
#' percent intensity values, ideally the output of find_ms().
#' @param cutoff a \emph{float} object, the low end cutoff for the MS data,
#' determined based on the MS devices used for analysis. Default is 50.
#' @param title a \emph{string} object. Title placed at the top of the
#' head-to-tail plot of best NIST Library match. Default title "Best NIST match".
#'
#'
#' @return a \emph{data.frame} object, a list of the top 10 best matching
#' compounds from the NIST database, with their compounds, the index in the
#' nistlist, and match percent.
#'
#' @export
comp_nist <- function(nistlist, ms, cutoff=50, title='Best NIST match'){
  # Organizing the sample ms data to match nist
  mz1 <- as.data.frame(cbind(ms$MZ,ms$Int))
  mz1 <- dplyr::filter(ms, ms$Int>=0.001)
  mz1$Int <- mz1$Int*999
  mz1$MZ <- round(mz1$MZ)
  mz1 <- dplyr::arrange(mz1, dplyr::desc(mz1$MZ), dplyr::desc(mz1$Int))
  drop <- c()
  for (i in 1:(length(mz1$MZ)-1)){
    if (mz1$MZ[i]==mz1$MZ[i+1]){
      drop <- append(drop,i+1)
    }
  }
  mz1 <- mz1[-drop,]
  mz1 <- dplyr::arrange(mz1, dplyr::desc(mz1$Int))
  sum <- sum(mz1$Int)

  # Finding entries with appropriate highest peak
  match_list <- list()
  match_no <- c()
  for (i in 1:length(nistlist)){
    mystr <- toString(mz1$MZ[1])
    if (mystr %in% substring(nistlist[[i]],1,nchar(mystr))){
      loc <- which(substring(nistlist[[i]],1,nchar(mystr))==mystr)
      split <- strsplit(nistlist[[i]][loc], ' ')
      myint <- as.numeric(split[[1]][2])
      if (myint >750){
        match_list <- append(match_list, list(nistlist[[i]]))
        match_no <- append(match_no, i)
      }
    }
  }

  if (length(match_list)==0){
    message('No viable matches in the NIST Library.')
    return(NULL)
  }

  # Calculating match percentages
  names <- c()
  formula <- c()
  ids <- c()
  match_per <- c()
  for (i in 1:length(match_list)){
    len <- length(match_list[[i]])
    per <- 0
    checked <- c()
    for (j in 8:len){
      k <- strsplit(match_list[[i]][j]," ")
      if (k[[1]][1] != "Num"){
        k <- c(as.integer(k[[1]][1]), as.double(k[[1]][2]))
        if (k[1] %in% mz1$MZ & k[1]>cutoff){
          loc <- which(mz1$MZ== k[1])
          per <- per+abs(mz1$Int[loc]-k[2])
          checked <- append(checked, loc)
        }
        else if (k[1]>cutoff){
          per <- per+k[2]
        }
      }
    }
    for (q in 1:length(mz1$MZ)){
      if (!(q %in% checked)){
        per <- per+mz1$Int[q]
      }
    }
    per <- (1-per/sum)*100
    match_per <- append(match_per,per)
    names <- append(names, strsplit(match_list[[i]][1],": ")[[1]][2])
    formula <- append(formula,
                      strsplit(match_list[[i]][2],": ")[[1]][2])
    ids <- append(ids, strsplit(match_list[[i]][5],": ")[[1]][2])
  }

  match_df <- as.data.frame(cbind(names,formula,ids))
  match_df$'Per' <- match_per
  match_df$'NIST number' <- match_no
  colnames(match_df) <- c('Name', 'Formula', 'ID', 'Per', 'NIST Number')
  match_df <- dplyr::arrange(match_df,dplyr::desc(match_df$Per))

  # Removing duplicates
  todel <- c()
  for (i in 2:length(match_df$Name)){
    if (match_df$Name[i]==match_df$Name[i-1] &
        match_df$Per[i]==match_df$Per[i-1]){
      todel <- append(todel,i)
    }
  }
  if (length(todel)>0){
    match_df <- match_df[-todel,]
  }
  # finding final top matches
  top_matches <- match_df[1:10,]
  j<-1
  while (j<=10){
    if (as.double(top_matches$'Per'[j])<0){
      top_matches <- top_matches[1:(j-1),]
      j <- 11
    }
    j<-j+1
  }
  message('The best NIST match is ', top_matches$Name[1], ', Formula: ',
          top_matches$Formula[1],'. ', round(top_matches$Per[1],2), '% match.')
  colnames(top_matches)<-c('Name','Formula','ID','Percent Match','Index')
  print(plot_nist(nistlist, top_matches$'Index'[1], mz1,
                  title))
  top_matches <- top_matches[,c(1,2,5,4)]
  rownames(top_matches) <- c(1:length(top_matches[,1]))
  return(top_matches)
}


#' @title  Plots the mass spectra of a NIST compound.
#'
#' @description  `plot_nist` Plots the mass spectra of a NIST compound.
#'
#' @details This function produces line plot of the mass spectra data from the
#' sample on top, and the mass spectrum from a NIST compound entry on the
#' bottom. The mass values vs the percent intensity values as a percent of the
#' highest intensity using \code{\link[ggplot2]{ggplot}} from ggplot2 package
#' \insertCite{ggplot2}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @param nistlist a \emph{list} object, a list of compound MS data from the
#' NIST MS Library database, ideally the output of nist_list().
#' @param k a \emph{integer} object, the index of the NIST compound in the
#' nistlist input.
#' @param ms a \emph{data.frame} object, a data frame of the mass values and the
#' percent intensity values, ideally the output of find_ms().
#' @param title a \emph{string} object. Title placed at the top of the plot.
#' Default title "Mass Spectrum".
#'
#' @return A \emph{ggplot} object. A line plot of the mass spectra data. The
#' mass values vs the percent intensity values as a percent of the highest
#' intensity.
#'
#' @export
plot_nist <- function(nistlist, k, ms, title = "NIST Mass Spectrum"){
  len <- length(nistlist[[k]])
  if (substring(nistlist[[k]][8],1,3)=="Num"){
    N <- 9
  }
  else if (substring(nistlist[[k]][7],1,3)=="Num"){
    N <- 8
  }
  else {
    for (line in 1:12){
      if (substring(nistlist[[k]][i],1,3)=="Num"){
        N <- i+1
      }
      else {
        stop("Cannot identify start point, send error to Stephanie.")
      }
    }
  }
  mz_nist <- c()
  int_nist <- c()
  for (i in N:len){
    mystr <- strsplit(nistlist[[k]][i]," ")
    mz_nist <- append(mz_nist, c(as.numeric(mystr[[1]][1]),as.numeric(mystr[[1]][1]),as.numeric(mystr[[1]][1])))
    int_nist <- append(int_nist, c(0,as.numeric(mystr[[1]][2])/999,0))
  }
  data <- as.data.frame(cbind(mz_nist,-1*int_nist))
  names(data) <- c("MZ", "Int")
  ms$Int <- ms$Int/999
  x <- c(mz_nist[1],mz_nist[length(mz_nist)])
  y <- c(0,0)
  data2 <- as.data.frame(cbind(x,y))
  fig <- plot_ms(ms, title= title) +
    ggplot2::geom_line(data=data, ggplot2::aes(x=mz_nist, y=-1*int_nist, color="NIST compound")) +
    ggplot2:: geom_line(data=data2, ggplot2::aes(x=x, y=y, color= "Sample data")) +
    ggplot2::labs(title= title, x='m/z value', y='percent intensity', color=' ')+
    ggplot2::scale_colour_manual(values = c('Sample data'="black", 'NIST compound'="blue"),
                                 guide = ggplot2::guide_legend(override.aes = list(
                                   linetype = c("solid", "solid"))))
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
