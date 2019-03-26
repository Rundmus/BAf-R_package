# -----------------------------------------------------------------------------#
#' parse raw file from FLEXMAP3D
#' 
#' The function to parse raw file attained directly from the FLEXMAP3D
#' instrument
#'
#' @param fn name of a file from FLEXMAP3D
#' 
#' @return 
#' a \code{list}, the elements of which are
#'   \item{head}{\code{'data.frame'} that contains information about the
#'   experiment, such as date, instrument id, etc}
#'   \item{pos}{\code{'data.frame'} indicates the position on plate}
#'   \item{other}{The corresponding data in \code{'data.frame'}}
#'   
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @seealso
#' \code{\link{read_LIMS_SBA_files}}
#' \code{\link{read_FlexMAP3D_csv}}
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-05-27 by Mun-Gwan
# modified : 
#   2011-06-15 by Mun-Gwan : can be applied also to 384 well plate 
#   2013-06-27 by Mun-Gwan : adapt to "SBAe" class
#   2015-05-29 by Mun-Gwan : allow gzip file
#   2019-01-16 by Mun-Gwan : fix a bug
# -----------------------------------------------------------------------------#

parse_FlexMAP3D_csv_2_list <- function(fn) {
  
  out <- list(head = data.frame(), pos = data.frame())
  # to compare if the location data in different section is identical 
  location <- c()		
  
  ## read the header section that includes various information about the current
  ## experiment such as date, time, volumns, etc
  strptime_mdylms <- function(x) strptime(x, "%m/%d/%Y %l:%M:%S %p")
  out$head <- gzfile(fn) %>%
    read.table(., sep = ",", stringsAsFactors = FALSE, nrows = 16, 
               skip = 4, row.names = 1, comment.char = "#") %>%
    t() %>%
    data.frame(stringsAsFactors = FALSE)
  out$head$BatchStartTime <- strptime_mdylms(out$head$BatchStartTime)
  out$head$BatchStopTime  <- strptime_mdylms(out$head$BatchStopTime)
  
  aLines <- readLines(conFn <- gzfile(fn), warn = FALSE)
  close(conFn)
  
  ## retrieve all lines to find the lines to indicate which data follows
  dataTypeLines <- which( substr(aLines, 2, 9) == "DataType" )
  
  ### read every section
  for(j in seq_along(dataTypeLines)) {
    start <- dataTypeLines[j]
    
    # find the data type
    dataType <- strsplit(aLines[start], ",")[[1]][2] %>%
      gsub(" ", "_", . ) %>%
      gsub("\"", "", . ) %>%
      tolower()
    
    # find the last line of the section that is seperated by a blank line
    sistaLine <- which(aLines == "") %>% 
      .[. > start] %>% 
      min() - 1
    
    # read the Section
    tg <- gzfile(fn) %>% 
      read.csv(skip= start, nrows= sistaLine - start - 1, as.is= TRUE) %>% {
        # Replace rownames with sample ID provided in the 2nd column which often
        # contain all "Unknown"s only
        if(names(.)[2] == "Sample") {
          rownames(.) <- make.unique(.[,2]) # 2019-01-16 : make.unique
          .[-2]
        } else .
      } %>% 
      # remove the column "Total.Events"
      .[names(.) != "Total.Events"]
    
    if(colnames(tg)[1] == "Location" && dataType != "warnings/errors") {
      if(identical(tg$Location, location)) {
        tg <- tg[-1]		# without "Location"
      } else {
        ## disassemble "Location" into "plate", "row", "col"
        dLoc <- strsplit(tg$Location, split="[(,)]") %>% 
          do.call("rbind", .) %>% 
          `[`(, 2:3) %>% 
          data.frame() %>% 
          `colnames<-`(c("plate", "pos")) %>% 	# plate, position eg. "A1", "A2"
          mutate(row = factor(substr(pos, 1, 1)),
                 col = factor(as.integer(substr(pos, 2, 3)))) %>% 
          dplyr::select(-2)          # remove dLoc$pos
        
        # if there was not previous section that contains "location"
        if(length(location) == 0) {
          location <- tg[,1]
          tg <- tg[-1]	# without "Location"
          out$pos <- dLoc
        } else {
          # bind all data
          tg <- cbind(dLoc, tg[-1])		# without "Location"
        }
      }
    }
    out[[dataType]] <- tg		# push the data in a list
  }
  return(out)
}



# -----------------------------------------------------------------------------#
#' Generate BAf object parsing raw file from FLEXMAP3D
#' 
#' This generates an BAf object parsing raw file attained directly from the
#' FLEXMAP3D instrument. The samples with median bead counts lower than 35 are
#' marked as "failed" in \code{@assay$fail_flag}
#' 
#' @param file name of a file from FLEXMAP3D
#' @param assayid the assay ID of this sample set and antibody set.
#' @param sinfo a \code{tbl_df} or \code{data.frame} that contains information
#'   about samples, in which each row is for each sample
#' @param bead a \code{tbl_df} or \code{data.frame} having bead information. 
#'   Often it can be obtained from other BAf object. If this is missing, a
#'   simple \code{tbl_df} that has only 'assay' column will be used as default.
#' @param sample_batch_c the name of sample batch column in the \code{sinfo}.
#'   e.g. \code{plate}
#' @param sbaid optionally, SBA ID of this antibody set can be given here.
#' @param bead_count_as_failed,bead_count_alert the cutoff values for bead
#'   count. Labeled in \code{'fail_flag'} as \code{'failed'} if the bead count
#'   is lower than \code{bead_count_as_failed}. If it is below
#'   \code{bead_count_alert}, then marked as \code{'lowBeadCount'}. Not in both
#'   cases, it is marked as \code{'ok'}.
#' @param reading_order the integers in which order the sample has read by
#'   FlexMAP3D. This parameter allows to accept the experimental data which is
#'   conducted in reverse or random order.
#' @param QC_plot if TRUE, some plots for quality check-up are generated like
#'   that of \code{\link{read_LIMS_SBA_files}}
#' @param path_QC where the QC plots will be saved. This is valid only if
#'   \code{QC_plot} is TRUE.
#' @param i_neg_ctrl_sample the index of negative control samples such as
#'   \code{EMPTY}. A logical vector of same length as samples can be given.
#' @param row_var,col_var When these are given, the variable in \code{sinfo} is
#'   compared with the values of row or column in \code{file} to confirm the
#'   order of samples in \code{sinfo}.
#'
#' @return an object of \link{BAf-class}
#'
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @seealso
#' \code{\link{read_LIMS_SBA_files}}
#' \code{\link{parse_FlexMAP3D_csv_2_list}}
#' 
#' @include plot_QC.SBA.R
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-08-31 by Mun-Gwan
# modified : 
#   2012-12-12 by Mun-Gwan 
#     add "fail_beadCount", "bead_count_alert" to 'read_FlexMAP3D_csv'
#   2013-06-27 by Mun-Gwan : adapt to "SBAe" class
#   2015-05-29 by Mun-Gwan : allow gzip file
#   2016-06-30 by Mun-Gwan : fail_flag for beads on the basis of bead count
#   2017-08-17 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

read_FlexMAP3D_csv <- 
  function(file, assayid, sinfo, bead, 
           sample_batch_c = "plate", sbaid, 
           bead_count_as_failed = 20, bead_count_alert = 35, 
           reading_order = 1:nrow(sinfo), 
           QC_plot = FALSE, path_QC = "../output/QC", 
           i_neg_ctrl_sample, row_var, col_var) {
    # a primitive check up of the arguments
    stopifnot(
      is.character(assayid),
      length(assayid) == 1L,
      !missing(file), !missing(sinfo), 
      sample_batch_c %in% names(sinfo), 
      bead_count_as_failed <= bead_count_alert, 
      length(reading_order) == nrow(sinfo)
    )

    # parse the file in the format used in FLEXMAP3D 
    parsed_file <- parse_FlexMAP3D_csv_2_list(file)

    ##  Sample  ----------------------------------------------------------------
    
    ###  Check Length of Samples  ###
    stopifnot(dim(parsed_file[["net_mfi"]])[1L] == dim(sinfo)[1L])
    
    ###  Check Position using 'row_var' or 'col_var' ###
    stop_if_unequal_pos <- function(parsed, x, var) {
      par_p <- parsed$pos[[x]] %>% as.character()
      si_p  <- sinfo[[var]] %>% as.character()
      if(!identical(par_p, si_p)) 
        stop("The ", x, "s in 'file' and 'sinfo' are NOT identical.")
    }
    if(!missing(row_var)) stop_if_unequal_pos(parsed_file, "row", row_var)
    if(!missing(col_var)) stop_if_unequal_pos(parsed_file, "col", col_var)
    
    ##  MFI values  ------------------------------------------------------------
    tg <- as.matrix(parsed_file[["net_mfi"]])[reading_order, ]

    ##  Bead  ------------------------------------------------------------------
    bead <- if(missing(bead)) {
      tibble(id = colnames(tg), 
             assay = rep(assayid, ncol(tg)) %>% factor()
      )
    } else {
      if(!"id" %in% names(bead)) bead$id <- colnames(tg)
      colnames(tg) <- bead$id
      bead %>% 
        as_tibble() %>% 
        mutate(assay = factor( rep(assayid, n()) ))
    }
    
    if(!missing(sbaid)) bead$sba <- sbaid %>% factor()
    
    ##  Bead count  ------------------------------------------------------------
    
    # bead count
    beadCount <- as.matrix(parsed_file$count)[reading_order, ]
    
    ##  Assy sample  -----------------------------------------------------------
    
    #  identify the failed or too low bead count sample / bead
    #  NOTE: 'bead_count_as_failed', 'bead_count_alert' are required.
    get_fail_mark <- function(med) {
      med %>% 
        cut(c(0, bead_count_as_failed, bead_count_alert, Inf), 
            c("failed", "lowBeadCount", "ok"), right= F)
    }
    #  median bead count per sample
    med_bc_s <- apply(beadCount, 1, median, na.rm = TRUE)

    assy_sample <- list(
      "reading_order"     = data.frame(reading_order), 
      "median_bead_count" = data.frame(med_bc_s),
      "fail_flag"         = data.frame(get_fail_mark(med_bc_s))
    ) %>%
      lapply(`names<-`, assayid)
    
    ##  Assy bead  -------------------------------------------------------------
    
    #  median bead count per bead
    med_bc_b <- apply(beadCount, 2, median, na.rm = TRUE)

    assy_bead <- list(
      "median_bead_count" = data.frame(med_bc_b),
      "fail_flag"         = data.frame(get_fail_mark(med_bc_b))
    ) %>%
      lapply(`names<-`, assayid)
    
    ##  generate BAf extracting only one data  ---------------------------------
    TG <- BAf(tg, 
              sinfo = sinfo[reading_order, ] %>% as_tibble(),
              binder = bead,
              sinfo_batch_i  = sample_batch_c, 
              binder_batch_i = "assay",
              assay_sinfo  = assy_sample,
              assay_binder = assy_bead
    )

    if(QC_plot) {
      ##  sort by 'plate' and 'pos' if they are available
      if(all(c("plate", "pos") %in% names(sI(TG)))) {
        tmp <- with(sI(TG), order(plate, pos))
        beadCount <- beadCount[tmp, ]
        TG <- TG[tmp, ]
      }
      if(!missing(i_neg_ctrl_sample)) {
        sid(TG)[i_neg_ctrl_sample] <- "EMPTY_0001"
      }
      if(! "plate" %in% names(sI(TG))) {
        sI(TG)$plate <- "one_plate"
      }
      ask_save.plot_QC.SBA(
        baf = TG,
        prefix = b_batch,
        path = path_QC,
        bead_count = beadCount,
        sample_order= order(batch(baf, "sample")),
        color_tbl_b = default_color_SBA_bead(T),
        color_tbl_s = default_color_SBA_ctrl_samples(sid(TG)),
        main = paste(unique(batch(TG, "binder")), collapse= ", "),
        ylab_signal = "Raw MFI",
        guide_line = FALSE
      )
    }
    return(invisible(TG))
  }

