# -----------------------------------------------------------------------------#
#' Design plate layout
#'
#' Design layouts of upto four 96-well plates by randomization with balances 
#' across those plates in terms of potential confounders. This can handle more
#' than 384 samples, then it produces layouts for plural 386-well plates.
#' 
#' - Balancing in the randomization are as below,\cr
#'   1) The number of samples of each group devided by 'categorical variables'
#'   are same or similar on each plate
#'     e.g. Male patients on each plate are same across all plates. It lets the
#'     ratios of male:female and patient:control are kept.
#'   2) The means of the values of 'quantitative variables' are homogeneous
#'   across all plates. 
#'     e.g. The mean of ages of a plate is not too much deviated from the mean
#'     of other plates.
#' 
#' @param sinfo a data frame that includes information about samples
#' @param cat_var categorical variables into consideration. The number of 
#'   samples in each category will be balanced (or kept constant) across the 
#'   plates. It should be a value or a vector of some of column names in 
#'   \var{sinfo}. When the number of samples are not a multiple of number of
#'   plates, some variation in numbers cannot be avoided. In such case, the
#'   variables listed to the left are prioritized during the balancing.
#' @param qty_var quantitative variables into consideration. The values of the 
#'   variables will be tested by ANOVA to check the balance and displayed in box
#'   plots. Like \var{cat_var}, it should be in \var{sinfo}
#' @param nNull_per_96 the number of null wells on each 96-well plate
#' @param null_var the name of the column in which 'NULL' will be shown to 
#'   indicate null well in \var{sinfo}. If it is given to NULL as default, the
#'   rows for null well will not be appended.
#' @param plate_nr the plate number(s) to which samples will be allocated. The 
#'   default is a sequential number.
#'   
#' @return The extended data frame of \var{sinfo} is generated with the columns
#' \code{'plate'}, \code{'pos'}, \code{'row'}, and \code{'col'}. \code{'plate'}
#' is designated plate ID. \code{'pos'} is the number indicating position within
#' plate (from top left to bottom right, move column-wise first). \code{'row'}
#' and \code{'col'} are the row letter and column number, respectively. The
#' plots to assist to compare the distribution of quantitative variables in each
#' plate are accompanied.
#' 
#' @examples
#' 
#' n <- 368
#' tg <- data.frame(
#'   "id" = c(paste("a", sample(n*2, n), sep="_"), rep("repl", 8)), 
#'   "age" = c(abs(rnorm(n, 50, 20)), rep(NA, 8)), 
#'   "diag_yrs" = c(abs(rnorm(n, 2)), rep(NA, 8)), 
#'   "gender" = c(rep(c("female", "male"), n/2), rep(NA, 8)), 
#'   "dis" = c(rep(c("case", "control"), each = n/2), rep("repl", 8))
#' )
#' X <- plateLayout.SBA(tg, c("dis", "gender"), c("age", "diag_yrs"), 
#'                     null_var= "id", plate_nr= 1:4)
#' head(tg)
#' head(X)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @include copies_from_Useful4me.R
#' @importFrom Useful2me shuffle
#' @export
# -----------------------------------------------------------------------------#
# created  : 2012-01-24 by Mun-Gwan
# modified : 
#   2012-02-22 by Mun-Gwan
#     1. Fix the problem to allocate multiple samples to same position
#     2. 'null_var'
#   2012-03-21 by Mun-Gwan : allow single 96 plate design
#   2012-03-22 by Mun-Gwan
#     1. add stripchart
#     2. fix adding Nulls for less than 4 plates
#     3. allow to select plate number when only one plate to be designed
#     4. fix the same iNull for all plates when there are more than 2 nulls per
#       plate
#   2012-11-09 by Mun-Gwan
#     1. fix the error that appears when there is too small number of samples in
#        a subgroup
#        - sample(with single value)
#        - solved by using shuffle instead
#   2013-02-04 by Mun-Gwan 
#     : allow the "plate" variable in the output to have larger than 4
#       Previously, the number > 4 has been changed to 5 -> 1, 6 -> 2, ...
#   2013-05-31 by Mun-Gwan 
#     : fix the problem arrisen by (sum(tX$res) %% nPlate) can be 2.4e-14
#     instead of 0
#   2014-10-21 by Mun-Gwan 
#     : Try more balanced distribution when there is remaining samples from
#     equal distribution
#   2015-09-15 by Mun-Gwan
#     : 'iNull' can accept 'plate_nr' when 'plate_nr' is not starting from 1
#   2015-12-08 by Mun-Gwan
#     : allow no catagorical variable
#   2016-05-31 by Mun-Gwan
#     : allow to limit to fixed columns
#   2017-10-26 by Mun-Gwan
#     : nPlate is computed using plate_nr when it was given
# -----------------------------------------------------------------------------#

plateLayout.SBA <- function(sinfo,
                            cat_var = NULL,
                            qty_var = NULL,
                            nNull_per_96 = 2,
                            null_var = NULL,
                            plate_nr = NULL) {
  
  stopifnot(is.null(plate_nr) || all(is.integer(plate_nr) & plate_nr >= 1 ),
            is.null(null_var) || null_var %in% colnames(sinfo))
  
  sinfo <- as.data.frame(sinfo)
  
  nSample <- nrow(sinfo)

  if(is.null(plate_nr)) {
    # calculate the number of required plates
    nPlate <- ceiling(nSample / (96 - nNull_per_96))
    plate_nr <- seq_len(nPlate)
  } else {
    nPlate <- length(plate_nr)
  }

  # compute the number of required columns
  nColumn <- ceiling((nSample + nNull_per_96 * nPlate) / 8)
  # the number of required columns per plate
  last_column <- ceiling(nColumn / nPlate)
  
  # number of sample per plate
  nSample__plate <- last_column * 8
  
  # Warn if the plate cannot be filled with samples and given number of null
  # wells
  if(nSample != nPlate * (nSample__plate - nNull_per_96) ) {
    warning("The number of samples is not enough to fill all wells!\n",
            'The assay will have', 
            (nSample__plate - nNull_per_96) - nSample/nPlate, 
            "more 'null' well(s) per plate.")
    nNull_per_96 <- nSample__plate - floor(nSample/nPlate)
  }
  
  
  #>> Null well position >>>>>------------------------------------------
  
  # the pre-determined position of the first two null wells
  #  plate 1 = B02, D10   plate 2 = C02, E10   plate 3 = D02, F10
  #  plate 4 = E02, G10   plate 5 = F02, D11   plate 6 = G02, E11
  #  plate 7 = B03, F11   plate 8 = C03, G11   plate 9 = D03, B11
  #  plate 10 = E03, C11  plate 11 = F03, D11  plate 12 = G03, E11
  iNull <-  lapply(plate_nr, function(j) {
    pos1 <- 9 + j + floor((j - 1) / 6) * 2
    pos2 <- 76 + ((j - 1) %% 4) + (
      if((j - 1) %% 12 > 7) {
        6 
      } else if((j - 1) %% 12 > 3) {
        8 
      } else 0 
    )
    ea <- c(pos1, pos2)
    
    ##  Move the Null well into the columns allowed for this
    lc8 <- ea > last_column * 8    # the nulls after the last column
    ea[lc8] <- ea[lc8] - (12 - last_column) * 8
    stopifnot(anyDuplicated(ea) == 0)
    return(ea)
  }) 

  # From the third null well, their position will be determined at random.
  if(nNull_per_96 > 2) {
    iNull <- lapply(iNull, function(ea) {
      c(ea, shuffle(c(1:nSample__plate)[ -ea ], nNull_per_96 - 2))
    })
  }
  # if the number of null wells in each 96-well plate is only 1,
  if(nNull_per_96 < 2) iNull <- lapply(iNull, function(j) j[1])
  
  #>> Number of 96-well plate >>>>>-------------------------------------
  
  
  ######################################################################
  #>>>>>  Samples of EACH KIND  <<<<<-----------------------------------
  ######################################################################
  ## - compute the number of samples of each kind to be included in each plate
  ## - total number of sample = fixed number per plate * number of plates +
  ##   residual samples 
  ## - (e.g. if 5 samples to 4 plates, 1 sample per plate and 1 sample on one
  ##   randomly chosen plate)
  null_cat_var <- FALSE
  if(is.null(cat_var)) {
    cat_var <- "tEmPoRaRy"      # artificial cat_var
    stopifnot(all(colnames(sinfo) != "tEmPoRaRy"))
    sinfo[[cat_var]] <- "ALL"
    null_cat_var <- TRUE
  }
  
  tX <- table(sinfo[cat_var], useNA = "ifany")/nPlate		# landscape table
  cat(">> The table of frequencies\n")
  print(tX * nPlate)
  
  #>> Turn to longitudinal table <<<<<<<<<
  tX <- as.data.frame(tX)
  if(length(cat_var) == 1) colnames(tX)[1] <- cat_var
  
  tX <- tX %>% 
    dplyr::filter(Freq != 0) %>%      	# remove 0 frequency
    # attempt to have better balance for the variables appeared earlier
    arrange(rev(desc(parse(text= cat_var) %>% eval()))) # sort by cat_var
  
  # the number of kinds dividing the samples based on the categorical variables
  nr <- nrow(tX)
  
  #>> Random number generation <<<<<<<<<<<
  # the generation of vectors including random numbers for the positions
  # excluding the null wells
  rnd <- seq_along(plate_nr) %>% 
    sapply(. %>% {c(1:nSample__plate)[ -iNull[[.]] ]} %>% 
             shuffle(., nSample__plate - nNull_per_96)
    )
  
  tX <- tX %>% 
    # the minimum number of samples of each kind - rounded down integer
    mutate(int = floor(Freq),
           # cumulative sum of above to be used as an index in 'rnd'
           cumS = cumsum(int)
    )
  
  #>> First position assignment <<<<<<<<<<
  # the position of the samples of each kind - excluding the samples in variable number
  tX$pos <- lapply(1:nr, function(j) {
    if(tX$int[j] == 0) {NULL} 		# if no sample to be assigned
    else {rnd[(if(j == 1) 1 else tX$cumS[j-1]+1):tX$cumS[j], ]}
  })
  
  ######################################################################
  #>>>>>  RESIDUAL samples  <<<<<---------------------------------------
  ######################################################################
  ## - because total number of samples of each kind is not multiple of 'nPlate'
  tX <- tX %>% 
    # the number of residual samples on all plates
    mutate(res = round((Freq - int) * nPlate))
  
  # if the total number of resituals is not multiple of 'nPlate', then stop
  if(sum(tX$res) %% nPlate != 0) 
    stop("The total number of samples should be a multiple of", 
         " the number of 96-well plates!")
  
  if(any(tX$res != 0)) {				# if there is residual samples,
    ### plate ID allocation for residual samples
    
    # randomize plate IDs for all residual samples
    rndPlate <- rep(shuffle(seq_along(plate_nr)), sum(tX$res) / nPlate)
    
    tX$resPlate <- list(NA)
    resCum <- 0
    for(j in 1:nr) {
      if(tX$res[j] != 0) {
        tX$resPlate[[j]] <- sort(rndPlate[(resCum + 1):(resCum + tX$res[j])])
      }
      resCum <- resCum + tX$res[j]
    }
    # confirm the same frequency of plate ID appearance in tX$resPlate
    stopifnot(all(table(unlist(tX$resPlate)) == round(sum(tX$res) / nPlate)))
    
    #> rPlateTable : plate ID table - column = plate ID, row = kind of samples
    rPlateTable <- matrix(FALSE, nr, nPlate)
    # convert the information in 'tX$resPlate' to the table 'rPlateTable' with
    # logical values
    for(j in 1:nr) {rPlateTable[j, tX$resPlate[[j]]] <- TRUE}
    # convert the table to contain randomized position information instead of
    # logical values
    rPlateTable[rPlateTable] <- 
      rnd[(tX$cumS[nr] + 1):(tX$cumS[nr] + (resCum / nPlate)), ]
    # store the information in 'tX$resPos'
    tX$resPos <- lapply(1:nr, function(j) {
      rPlateTable[j,][rPlateTable[j,] != 0]
    })
  }
  
  ## the columns of interesting categorical variables 
  subX <- sinfo[cat_var]
  for(eacat_var in cat_var) {
    # avoid invalid level
    levels(subX[,eacat_var]) <- c(levels(subX[,eacat_var]), "<NA>")
  }
  subX[is.na(subX)] <- "<NA>"		# replace missing with <NA>
  
  ## replace missing with <NA>
  for(eacat_var in cat_var) {
    levels(tX[,eacat_var]) <- c(levels(tX[,eacat_var]), "<NA>")
  }
  tX[,cat_var][is.na(tX[,cat_var])] <- "<NA>"
  
  # logical vector : the rows in 'sinfo' of the samples corresponding to each row in 'tX'
  tX$rowX <- lapply(1:nr, function(j) apply(subX, 1, function(k) {
    all(k == tX[j,cat_var])
  }) )
  
  ### Allocate samples #######################################################
  
  ## initialization
  
  ##  Give warnings if newly generated variables exist in input table
  warnings_if_new_variables_exist <- function(x, chq) {
    if(chq %in% colnames(x)) {
      # replace with "...old"
      colnames(x)[colnames(x) == chq] <- paste0(chq, ".old")
      warning("The variable '", chq, "' has been replaced with '", 
              chq, ".old'")
    }
    return(x)
  }
  sinfo <- sinfo %>% 
    warnings_if_new_variables_exist("plate") %>% 
    warnings_if_new_variables_exist("pos")
  
  sinfo$plate <- NA
  sinfo$pos <- NA
  
  for(j in 1:nr) {
    iX <- shuffle(which(tX$rowX[[j]]))
    # plate number
    sinfo$plate[iX] <- if(tX$res[j] == 0) {
      rep(plate_nr, each = tX$int[j]) 
    } else {
      c(rep(plate_nr, each = tX$int[j]), plate_nr[tX$resPlate[[j]]])
    }
    # position number in the plate given above
    sinfo$pos[iX] <- 
      if(tX$res[j] == 0) tX$pos[[j]] else c(tX$pos[[j]], tX$resPos[[j]])
  }
  
  ### >>>>> 'sinfo' : The END of Plate layout design <<<<< ###
  
  ## check if any multiple samples were allocated to same well
  if(anyDuplicated(paste(sinfo$plate, sinfo$pos)) != 0) {
    stop("By unknown reason, multiple samples were allocated to same well!")
  }
  
  ### Visualize the quality of layout ########################################
  
  ##### Quantitative variables
  ### - compare the values of the samples in each plate
  
  if(! is.null(qty_var) && nPlate > 1) {
    
    cat("\n>> The table of mean values\n")
    ### in tables
    lapply(qty_var, function(qV) {
      sapply(1:nr, function(j) {
        sapply(plate_nr, function(iPlate) {
          sinfo[tX$rowX[[j]] & (sinfo$plate == iPlate), ] %>% 
            `[[`(qV) %>% 
            mean()
        })
      }) %>% 
        `colnames<-`(
          sapply(1:nr, . %>% {
            tX[., cat_var]
          } %>% 
            unlist() %>% 
            paste(collapse = "-")
          )
        ) %>% 
        `rownames<-`(paste("plate", plate_nr)) %>% 
        t()
    }) %>% 
      `names<-`(qty_var) %>% 
      print(.)
    
    ### show a summary table for whole samples in a plate
    cat("\n>> The summary table\n")
    lapply(qty_var, function(qV) {
      tab <- by(sinfo[[qV]], sinfo$plate, summary)
      if(any(sapply(tab, length) > 6)) {
        for(ii in seq(tab)) {
          if(!"NA\'s" %in% names(tab[[ii]])) tab[[ii]][["NA\'s"]] <- NA
        }
      }
      tab %>% 
        do.call("rbind", .) %>% 
        `rownames<-`(paste("plate", plate_nr))
    }) %>% 
      `names<-`(qty_var) %>% 
      print(.)
    
    ### Graphical visualization of the value distribution
    # design proper number of columns and rows to show box plots
    nRowPlot <- min((nr - 1) %/% c(5:8)) + 1
    nColPlot <- min(nr, c(5:8)[which.min((nr - 1) %/% c(5:8))])
    for(qV in qty_var) {
      dev.new(title= qV, width= min(20, nColPlot * 2), height= min(12, nRowPlot * 3))
      yr <- range(sinfo[,qV], na.rm= TRUE)
      # adjust parameters to obtain better looking plots
      opar <- par(
        mfrow = c(nRowPlot, nColPlot),
        mar = 2.5 + c(0, 0, 0.5, -1),
        oma = c(2, 0, 3, 0)
      )
      
      # show in boxplots
      for(j in 1:nr) {
        Xs <- sinfo[tX$rowX[[j]], c(qV, "plate")]
        if(all(is.na(Xs[[qV]]))) {		# no data for the variable is available
          boxplot(rep(1, nrow(Xs)) ~ Xs$plate, ylim = yr)
        } else {
          # if the number of plate is only 1
          if(length(unique(Xs$plate)) == 1) {
            boxplot(Xs[[qV]] ~ Xs$plate, ylim = yr)
          } else {
            if(length(Xs[[qV]]) <= length(unique(Xs$plate)) * 2) {
              boxplot_nEle(Xs[[qV]] ~ Xs$plate, 
                                      ylim = yr, 
                                      test = "")
            } else {
              boxplot_nEle(Xs[[qV]] ~ Xs$plate, 
                                      ylim = yr, 
                                      test = "aov")
            }
            # add dots
            stripchart(Xs[[qV]] ~ Xs$plate, method= "jitter", add= TRUE, 
                       vertical= TRUE, pch= 19, bg= "bisque", cex= 0.4, 
                       ylim= yr, col= "seagreen")
          }
        }
        title(main = paste(unlist(tX[j,cat_var]), collapse = "-"))
      }
      title(main = qV, xlab = "Plate", outer = TRUE, col.main = "blue", 
            cex.main = 2, col.lab = "blue", cex.lab = 1.5, mgp = c(0, 1, 0))
      
      par(opar)		# return to original parameter set for graphics
    }
  }
  ##  Remove the temporarily added column when "cat_var" is missing
  if(null_cat_var) {
    sinfo <- sinfo[, colnames(sinfo) != cat_var]
  }
  
  ### Final tyding up the output #############################################
  
  #> Add nulls at the end of table <<<<<<<
  if(!is.null(null_var)) {		# if user want to add nulls in output
    nullT <- data.frame(matrix(NA, nNull_per_96 * nPlate, ncol(sinfo))) %>% 
      `colnames<-`(colnames(sinfo))
    # Fill the column chosen by user with 'Null' to indicate Null well
    if(null_var %in% colnames(sinfo)) nullT[, null_var] <- "NULL"
    nullT$plate <- rep(plate_nr, each = nNull_per_96)
    nullT$pos <- do.call("c", iNull[seq_along(plate_nr)])
    sinfo <- rbind(sinfo, nullT)
  }
  
  ##  Give warnings if newly generated variables exist in input table
  sinfo <- sinfo %>% 
    warnings_if_new_variables_exist("col") %>% 
    warnings_if_new_variables_exist("row") %>% 
    ## convert position number to column and row ids
    mutate(row = factor(LETTERS[(pos - 1) %% 8 + 1]),
           col = factor((pos - 1) %/% 8 + 1))
  
  # replace plate number with given one
  #? if(!missing(plate_nr)) sinfo$plate <- plate_nr[sinfo$plate]
  
  if(nPlate > 4) {
    sinfo$plate384 <- factor(((sinfo$plate - 1) %/% 4) + 1)
  }
  #>> 384 well plate position <<<<<<<<<<<<
  ##  Give warnings if newly generated variables exist in input table
  sinfo <- sinfo %>% 
    warnings_if_new_variables_exist("row384") %>% 
    warnings_if_new_variables_exist("col384") %>% 
    
    mutate(
      row384 = LETTERS[((pos - 1) %% 8 + 1)*2 - plate %% 2],
      col384 = (((plate - 1) %% 4) - 2) %/% 2 + 2 * ((pos - 1) %/% 8 + 1)
    ) %>% 
    mutate_at(c("row384", "col384", "plate"), as.factor)
  
  return(invisible(sinfo))		# return input data.frame with positional information
}
