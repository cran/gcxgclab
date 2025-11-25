#' @title  Creates list of atomic mass data
#'
#' @description  `mass_list` creates a list of atomic mass data
#'
#' @details This function creates a data frame containing the data for the
#' atomic weights for each element in the periodic table
#' \insertCite{mass}{gcxgclab}.
#'
#' @references
#' \insertAllCited{}
#'
#' @return A \emph{data.frame} object, with two columns, (elements,
#' mass).
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' mz <- find_ms(frame, t_peak=peaks$'T'[1])
#' masslist <- mass_list()
#' non_targeted(masslist, mz, THR=0.05)
#'
#' @export
mass_list <- function(){
  element <- c('H','D','T','He','He','Li','Li','Be','B','B','C','C','C','N',
               'N','O','O','O','F','Ne','Ne','Ne','Na','Mg','Mg','Mg','Al','Si',
               'Si','Si','P','S','S','S','S','Cl','Cl','Ar','Ar','Ar','K','K',
               'K','Ca','Ca','Ca','Ca','Ca','Ca','Sc','Ti','Ti','Ti','Ti','Ti',
               'V','V','Cr','Cr','Cr','Cr','Mn','Fe','Fe','Fe','Fe','Co','Ni',
               'Ni','Ni','Ni','Ni','Cu','Cu','Zn','Zn','Zn','Zn','Zn','Ga','Ga',
               'Ge','Ge','Ge','Ge','Ge','As','Se','Se','Se','Se','Se','Se','Br',
               'Br','Kr','Kr','Kr','Kr','Kr','Kr','Rb','Rb','Sr','Sr','Sr','Sr',
               'Y','Zr','Zr','Zr','Zr','Zr','Nb','Mo','Mo','Mo','Mo','Mo','Mo',
               'Tc','Tc','Tc','Ru','Ru','Ru','Ru','Ru','Ru','Ru','Rh','Pd','Pd',
               'Pd','Pd','Pd','Pd','Ag','Ag','Cd','Cd','Cd','Cd','Cd','Cd','Cd',
               'Cd','In','In','Sn','Sn','Sn','Sn','Sn','Sn','Sn','Sn','Sn','Sn',
               'Sb','Sb','Te','Te','Te','Te','Te','Te','Te','Te','I','I','Xe','Xe',
               'Xe','Xe','Xe','Xe','Xe','Xe','Xe','Cs','Ba','Ba','Ba','Ba','Ba',
               'Ba','Ba','La','La','Ce','Ce','Ce','Ce','Pr','Nd','Nd','Nd','Nd',
               'Nd','Nd','Nd','Pm','Pm','Sm','Sm','Sm','Sm','Sm','Sm','Sm','Eu',
               'Eu','Gd','Gd','Gd','Gd','Gd','Gd','Gd','Tb','Dy','Dy','Dy','Dy',
               'Dy','Dy','Dy','Ho','Er','Er','Er','Er','Er','Er','Tm','Yb','Yb',
               'Yb','Yb','Yb','Yb','Yb','Yb','Lu','Lu','Hf','Hf','Hf','Hf','Hf',
               'Hf','Ta','Ta','W','W','W','W','W','Re','Re','Os','Os','Os','Os',
               'Os','Os','Os','Ir','Ir','Pt','Pt','Pt','Pt','Pt','Pt','Au','Hg',
               'Hg','Hg','Hg','Hg','Hg','Hg','Tl','Tl','Pb','Pb','Pb','Pb','Bi',
               'Pb','Pb','At','At','Rn','Rn','Rn','Fr','Ra','Ra','Ra','Ra','Ac',
               'Th','Th','Pa','U','U','U','U','U','Np','Np','Pu','Pu','Pu','Pu',
               'Pu','Pu','Am','Am','Cm','Cm','Cm','Cm','Cm','Cm','Bk','Bk','Cf',
               'Cf','Cf','Cf','Es','Fm','Md','Md','No','Lr','Rf','Db','Sg','Bh',
               'Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og')
  isotope <- c(1,2,3,3,4,6,7,9,10,11,12,13,14,14,15,16,17,18,19,20,21,22,23,24,
               25,26,27,28,29,30,31,32,33,34,36,35,37,36,38,40,39,40,41,40,42,
               43,44,46,48,45,46,47,48,49,50,50,51,50,52,53,54,55,54,56,57,58,
               59,58,60,61,62,64,63,65,64,66,67,68,70,69,71,70,72,73,74,76,75,
               74,76,77,78,80,82,79,81,78,80,82,83,84,86,85,87,84,86,87,88,89,
               90,91,92,94,96,93,92,94,95,96,97,98,100,97,98,99,96,98,99,100,
               101,102,104,103,102,104,105,106,108,110,107,109,106,108,110,111,
               112,113,114,116,113,115,112,114,115,116,117,118,119,120,122,124,
               121,123,120,122,123,124,125,126,128,129,130,127,124,126,128,129,130,
               131,132,134,136,133,130,132,134,135,136,137,138,138,139,136,138,
               140,142,141,142,143,144,145,146,148,150,145,147,144,147,148,149,
               150,152,154,151,153,152,154,155,156,157,158,160,159,156,158,160,
               161,162,163,164,165,162,164,166,167,168,170,169,168,170,171,172,
               173,174,176,175,176,174,176,177,178,179,180,180,181,180,182,183,
               184,186,185,187,184,186,187,188,189,190,192,191,193,190,192,194,
               195,196,198,197,196,198,199,200,201,202,204,203,205,204,206,207,
               208,209,209,210,210,211,211,220,222,223,223,224,226,228,227,230,
               232,231,233,234,235,236,238,236,237,238,239,240,241,242,244,241,
               243,243,244,245,246,247,248,247,249,249,250,251,252,252,257,258,
               260,259,262,267,268,271,272,270,276,281,280,285,284,289,288,293,
               292,294)
  mass <- c(1.00782503223,2.01410177812,3.0160492779,3.0160293201,4.00260325413,
            6.0151228874,7.0160034366,9.012183065,10.01293695,11.00930536,
            12.0000000,13.00335483507,14.0032419884,14.00307400443,
            15.00010889888,15.99491461957,16.99913175650,17.99915961286,
            18.99840316273,19.9924401762,20.993846685,21.991385114,
            22.9897692820,23.985041697,24.985836976,25.982592968,26.98153853,
            27.97692653465,28.97649466490,29.973770136,30.97376199842,
            31.9720711744,32.9714589098,33.967867004,35.96708071,34.968852682,
            36.965902602,35.967545105,37.96273211,39.9623831237,38.9637064864,
            39.963998166,40.9618252579,39.962590863,41.95861783,42.95876644,
            43.95548156,45.9536890,47.95252276,44.95590828,45.95262772,
            46.95175879,47.94794198,48.94786568,49.94478689,49.94715601,
            50.94395704,49.94604183,51.94050623,52.94064815,53.93887916,
            54.93804391,53.93960899,55.93493633,56.93539284,57.93327443,
            58.93319429,57.93534241,59.93078588,60.93105557,61.92834537,
            63.92796682,62.92959772,64.92778970,63.92914201,65.92603381,
            66.92712775,67.92484455,69.9253192,68.9255735,70.92470258,
            69.92424875,71.922075826,72.923458956,73.921177761,75.921402726,
            74.92159457,73.922475934,75.919213704,76.919914154,77.91730928,
            79.9165218,81.9166995,78.9183376,80.9162897,77.92036494,79.91637808,
            81.91348273,82.91412716,83.9114977282,85.9106106269,84.9117897379,
            86.9091805310,83.9134191,85.9092606,86.9088775,87.9056125,
            88.9058403,89.9046977,90.9056396,91.9050347,93.9063108,95.9082714,
            92.9063730,91.90680796,93.90508490,94.90583877,95.90467612,
            96.90601812,97.90540482,99.9074718,96.9063667,97.9072124,98.9062508,
            95.90759025,97.9052868,98.9059341,99.9042143,100.9055769,
            101.9043441,103.9054275,102.9054980,101.9056022,103.9040305,
            104.9050796,105.9034804,107.9038916,109.90517220,106.9050916,
            108.9047553,105.9064599,107.9041834,109.90300661,110.90418287,
            111.90276287,112.90440813,113.90336509,115.90476315,112.90406184,
            114.903878776,111.90482387,113.9027827,114.903344699,115.90174280,
            116.90295398,117.90160657,118.90331117,119.90220163,121.9034438,
            123.9052766,120.9038120,122.9042132,119.9040593,121.9030435,
            122.9042698,123.9028171,124.9044299,125.9033109,127.90446128,
            128.9049836,129.906222748,126.9044719,123.9058920,125.9042983,127.9035310,
            128.9047808611,129.903509349,130.90508406,131.9041550856,
            133.90539466,135.907214484,132.9054519610,129.9063207,131.9050611,
            133.90450818,134.90568838,135.90457573,136.90582714,137.90524700,
            137.9071149,138.9063563,135.90712921,137.905991,139.9054431,
            141.9092504,140.9076576,141.9077290,142.9098200,143.9100930,
            144.9125793,145.9131226,147.9168993,149.9209022,144.9127559,
            146.9151450,143.9120065,146.9149044,147.9148292,148.9171921,
            149.9172829,151.9197397,153.9222169,150.9198578,152.9212380,
            151.9197995,153.9208741,154.9226305,155.9221312,156.9239686,
            157.9241123,159.9270624,158.9253547,155.9242847,157.9244159,
            159.9252046,160.9269405,161.9268056,162.9287383,163.9291819,
            164.9303288,161.9287884,163.9292088,165.9302995,166.9320546,
            167.9323767,169.9354702,168.9342179,167.9338896,169.9347664,
            170.9363302,171.9363859,172.9382151,173.9388664,175.9425764,
            174.9407752,175.9426897,173.9400461,175.9414076,176.9432277,
            177.9437058,178.9458232,179.9465570,179.9474648,180.9479958,
            179.9467108,181.94820394,182.95022275,183.95093092,185.9543628,
            184.9529545,186.9557501,183.9524885,185.9538350,186.9557474,
            187.9558352,188.9581442,189.9584437,191.9614770,190.9605893,
            192.9629216,189.9599297,191.9610387,193.9626809,194.9647917,
            195.96495209,197.9678949,196.96656879,195.9658326,197.96676860,
            198.96828064,199.96832659,200.97030284,201.97064340,203.97349398,
            202.9723446,204.9744278,203.9730440,205.9744657,206.9758973,
            207.9766525,208.9803991,208.9824308,209.9828741,209.9871479,
            210.9874966,210.9906011,220.0113941,222.0175782,223.0197360,
            223.0185023,224.0202120,226.0254103,228.0310707,227.0277523,
            230.0331341,232.0380558,231.0358842,233.0396355,234.0409523,
            235.0439301,236.0455682,238.0507884,236.046570,237.0481736,
            238.0495601,239.0521636,240.0538138,241.0568517,242.0587428,
            244.0642053,241.0568293,243.0613813,243.0613893,244.0627528,
            245.0654915,246.0672238,247.0703541,248.0723499,247.0703073,
            249.0749877,249.0748539,250.0764062,251.0795886,252.0816272,
            252.082980,257.0951061,258.0984315,260.10365,259.10103,262.10961,
            267.12179,268.12567,271.13393,272.13826,270.13429,276.15159,
            281.16451,280.16514,285.17712,284.17873,289.19042,288.19274,
            293.20449,292.20746,294.21392)

  massdf <- as.data.frame(cbind(element,isotope,mass))
  return(massdf)
}

#' @title  Compares MS to atomic mass data
#'
#' @description  `non_targeted` compares the MS data from a peak to atomic mass
#' data.
#'
#' @details This function takes the MS data from an intensity peak in a sample
#' and compares it to combinations of atomic masses. Then it approximates the
#' makeup of the compound, giving the best matches to the MS data. Note that
#' the default matches will contain only H, N, C, O, F, Cl, Br, I, and Si. The
#' user can input optional parameters to indicate additional elements to be
#' considered or restrictions on the number of any specific element in the
#' matching compounds.
#'
#' @param masslist a \emph{list} object, a list of atomic weights, ideally the
#' output of mass_list().
#' @param ms a \emph{data.frame} object, a data frame of the mass values and the
#' percent intensity values, ideally the output of find_ms().
#' @param THR a \emph{double} object. The threshold of intensity of which to
#' include peaks for mass comparison. Default is 0.1.
#' @param ... a \emph{vector} object. Any further optional inputs which
#' indicate additional elements to consider in the compound, or restrictions on
#' the number of a certain element in the compound. Should be in the form
#' c('X', a, b) where X = element symbol, a = minimum number of atoms, b =
#' maximum number of atoms. a and b are optional. If no minimum, use a=0, if no
#' maximum, do not include b.
#'
#' @return A \emph{list} object, a list of vectors containing strings of the
#' matching compounds.
#'
#' @examples
#' file <- system.file("extdata","sample1.cdf",package="gcxgclab")
#' frame <- extract_data(file,mod_t=.5)
#' peaks <- top_peaks(frame$TIC_df, 5)
#' mz <- find_ms(frame, t_peak=peaks$'T'[1])
#' masslist <- mass_list()
#' non_targeted(masslist, mz, THR=0.05)
#'
#' @export
non_targeted <- function(masslist, ms, THR=0.1, ...){
  # Simplifying mz for identification
  tolerance <- 0.0
  mz1 <- dplyr::filter(ms, ms$Int>THR)
  mz1 <- dplyr::arrange(mz1,mz1$MZ,mz1$Int)
  i <- 1
  while (i<length(mz1$'MZ')){
    if (abs(mz1$'MZ'[i]-mz1$'MZ'[i+1])<0.01){
      if (mz1$'Int'[i]< mz1$'Int'[i+1]){
        mz1<-mz1[-i,]
      }
      else{
        mz1<-mz1[-(i+1),]
      }
    }
    else {i<-i+1}
  }

  # Creating mass list, adding additional indicated elements
  max_mz <- mz1$MZ[length(mz1$MZ)]
  def <- max_mz-round(max_mz)
  if (def>0 & def<0.33){
    halogen <- FALSE
    nec <- c(11,1,14,16,28,170)
  }
  else {
    halogen <- TRUE
    nec <- c(11,1,14,16,19,28,36,94,95,170)
  }
  res <- list(...)
  if (length(res)>0){
    for (i in 1:length(res)){
      el <- res[[i]][1]
      iso <- FALSE
      if (grepl(" ",el)){
        iso <- TRUE
        is <- strsplit(el," ")[[1]][2]
        el <- strsplit(el," ")[[1]][1]
        loc2 <- which(masslist$isotope==is)
      }
      loc <- which(masslist$element==el)
      if (iso){
        loc1 <- intersect(loc,loc2)[1]
      }
      else {
        loc1 <- loc[1]
      }
      if (length(loc1)==0){
        stop(paste(res[[i]][1], 'is not an atomic symbol.'))
      }
      else{
        if (!(loc1 %in% nec)){
          nec <- append(nec,loc1)
        }
      }
    }
  }
  masslist2 <- masslist[nec,]
  nel <- length(nec)


  # Separate peak grouping
  len <- length(mz1$'MZ')
  dfs <- list()
  i<-2
  j<-1
  while (i<=len){
    if (mz1$'MZ'[i]-mz1$'MZ'[i-1]>3){
      dfs <- append(dfs, list(mz1[j:(i-1),]))
      #if (i==(len-1)){
       # dfs <- append(dfs, list(mz1[i,]))
      #}
      j<-i
    }
    if (i==len){
      if (mz1$'MZ'[i]-mz1$'MZ'[i-1]<3){
        dfs <- append(dfs, list(mz1[j:(i-1),]))
      }
      dfs <- append(dfs, list(mz1[i,]))
    }
    i <- i+1
  }
  if (len==1){
    dfs <- append(dfs, list(mz1[1,]))
  }
  matches2 <- list()
  mznames <- c()
  mzlist <- c()

  finished <- FALSE
  while (!finished){
    if (tolerance>1){
      message("Could not find matching compounds within reasonable tolerance.")
      return(NULL)
    }

  # calculate average and maxes of each group
  masses <- as.double(masslist2$mass)
  for (i in 1:length(dfs)){
    sum <- 0
    tot <- 0
    for (j in 1:length(dfs[[i]]$'Int')){
      sum <- sum + dfs[[i]]$'Int'[j]
      tot <- tot + dfs[[i]]$'Int'[j]*dfs[[i]]$'MZ'[j]
    }
    maxx <- dfs[[i]]$'MZ'[which.max(dfs[[i]]$'Int')]
    if (tolerance==0){
      mzlist <- append(mzlist,as.character(maxx))
    }
    tol <- round(tolerance*100)
    N <- round(maxx*100)
    mass <- round((tot/sum)*100)
    masses2 <- round(masses*100)

    # solve linear Diophantine problem to find compound makeup
    MM <- floor(maxx/2)
    if (tolerance ==0){
      range <- c(N)
    }
    else {
      range <- c(N-tol,N+tol)
    }
    for (num in range){
      sol <- nilde::nlde(a=masses2, n=num, M=MM)
      if (sol$p.n>0){
        for (n in 1:sol$p.n){
          matches2 <- append(matches2, list(as.vector(sol$solution[,n])))
          mznames <- append(mznames, as.character(maxx))
        }
      }
    }
  }

  # Removing results outside of reasonable element ratios and outside element
  # restrictions from optional input
  todel <- c()
  if (length(res)>0){
    for (i in 1:length(res)){
      if (length(res[[i]])==2){
        a <- as.integer(res[[i]][2])
        loc <- which(masslist2$element==res[[i]][1])
        massi <- which(masses2==round(as.double(masslist2$mass[loc])*100))
        for (j in 1:length(matches2)){
          if (matches2[[j]][massi]<a){
            todel <- append(todel,j)
          }
        }
      }
      else if (length(res[[i]])==3){
        a <- as.integer(res[[i]][2])
        b <- as.integer(res[[i]][3])
        loc <- which(masslist2$element==res[[i]][1])
        massi <- which(masses2==round(as.double(masslist2$mass[loc])*100))
        for (j in 1:length(matches2)){
          if (matches2[[j]][massi]<a | matches2[[j]][massi]>b){
            todel <- append(todel,j)
          }
        }
      }
      else if (length(res[[i]])>3){
        stop('Optional element input', res[[i]], 'has too many
                       entries.')
      }
    }
  }
  for (j in 1:length(matches2)){
    if (matches2[[j]][1]>0){
      HC <- matches2[[j]][2]/matches2[[j]][1]
      NC <- matches2[[j]][3]/matches2[[j]][1]
      OC <- matches2[[j]][4]/matches2[[j]][1]
      if (HC>7 | NC>6 | OC>3){
        todel <- append(todel, j)
      }
    }
    if (halogen){
      if (matches2[[j]][8]>0 & matches2[[j]][9]>0){
        todel <- append(todel, j)
      }
    }
  }
  if (length(todel)>0){
    todel <- unique(todel)
    matches2a <- matches2[-todel]
    mznames2 <- mznames[-todel]
  }
  else {
    matches2a <- matches2
    mznames2 <- mznames
  }

  # converting to list of matches, grouping for each peak
  matches2b <- vector("list",length(mzlist))
  names(matches2b) <- mzlist
  if (length(matches2a)>0){
    for (j in 1:length(matches2a)){
      if (mznames2[j] %in% mzlist){
        loc <- which(mzlist==mznames2[j])
        matches2b[[loc]]<- append(matches2b[[loc]],list(matches2a[[j]]))
      }
      else {
        stop("mass not in mzlist")
      }
    }
  }

  matches3 <- list()
  for (i in 1:length(matches2b)){
    m1 <- cbind(matches2b[[i]][[1]])
    if (length(matches2b[[i]])>1){
      for (j in 2:length(matches2b[[i]])){
        m1 <- rbind(cbind(m1, matches2b[[i]][[j]]))
      }
    }
    matches3 <- append(matches3, list(m1))
  }
  names(matches3) <- mzlist

  matches2c <- matches2b
  # Removing solutions which are incorrect fragments
  todel2 <- c()
  if (length(dfs)>1){
  for (j in 1:length(matches2c[[length(matches2c)]])){
    keep <- FALSE
    if (length(matches2c)>1){
      for (jj in 1:(length(matches2c)-1)){
        if (jj==1 | keep){
          keep_jj <- c()
          for (jjj in 1:length(matches2c[[jj]])){
            if (all(matches2c[[length(matches2c)]][[j]]>=matches2c[[jj]][[jjj]])){
                keep_jj <- append(keep_jj,TRUE)
            }
            else {
              keep_jj <- append(keep_jj,FALSE)
            }
          }
          if (any(keep_jj)){
            keep <- TRUE
          }
          else {
            keep <- FALSE
          }
        }
      }
    }
    if (!keep){
      todel2 <- append(todel2,j)
    }
  }
  }
  # # N rule
  # mzval <- as.numeric(names(matches2c)[length(matches2c)])
  # for (i in 1:length(matches2c[[length(matches2c)]])){
  #   if (round(mzval)%%2==0 & abs(round(mzval)-mzval)<.33){
  #     if (matches2c[[length(matches2c)]][[i]][3]%%2==1){
  #       todel2 <- append(todel2,i)
  #     }
  #   }
  #   else if (round(mzval)%%2==1 & abs(round(mzval)-mzval)<.33){
  #     if (matches2c[[length(matches2c)]][[i]][3]%%2==0){
  #       todel2 <- append(todel2,i)
  #     }
  #   }
  # }

  if (length(matches2c[[length(matches2c)]])-length(todel2)>=10){
    finished <- TRUE
    if (length(todel2)>0){
      matches2c[[length(matches2c)]] <- matches2c[[length(matches2c)]][-todel2]
    }
  }
  if (!finished){
    tolerance <- round(tolerance +0.01,2)
  }

  } # end if finished

  # Convert solutions to strings with elements
  matches2d <- list()
  for (i in 1:length(matches2c)){
    matches2e <- c()
    for (j in 1:length(matches2c[[i]])){
    sent <-''
      for (k in 1:nel){
        if (length(matches2c[[i]])>0){
          if (matches2c[[i]][[j]][k]!= 0){
            if (matches2c[[i]][[j]][k]==1){
              sent <- paste0(sent, masslist2$element[k])
            }
            else {
              sent <- paste0(sent, masslist2$element[k], matches2c[[i]][[j]][k])
            }
          }
        }
        else {
          sent <- c('none')
        }
      }
      matches2e <- append(matches2e, c(sent))
    }
    matches2d <- append(matches2d, list(matches2e))
  }
  names(matches2d)<- mzlist
  final_match <- matches2d[[length(matches2d)]][1:10]

  # Calculating masses of the final results
  final_masses <- c()
  for (i in 1:length(matches2c[[length(matches2c)]])){
    sum <- sum(matches2c[[length(matches2c)]][[i]]*as.numeric(masslist2$mass))
    final_masses <- append(final_masses,sum)
  }
  mass_diff <- abs(final_masses-as.numeric(mzlist[length(mzlist)]))
  p_diff <- 100-100*mass_diff/as.numeric(mzlist[length(mzlist)])
  final_df <- as.data.frame(cbind(matches2c[[length(matches2c)]],
                                  matches2d[[length(matches2d)]],final_masses,
                                  mass_diff,p_diff))
  names(final_df)<- c('Solution','Compound', 'Mass', 'Difference','Percent Match')
  final_df <- dplyr::arrange(final_df, mass_diff)

  message(paste0('The best mass solution found is ', final_df[1,2],'.'))
  return(final_df[1:10,c(2,3)])
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
