#' Read the files from LIMS
#' 
#' @param ... 
#' @param fn file name
#'
#' @return \code{\link{data.frame}}
#' @noRd
read_file_lims <- function(fn, ...) {
  read.delim(
    gzfile(fn),
    ...,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "#"
  )
}


#' Read LIMS binder file
#'
#' @param fn file name
#' @param ver version
#'
#' @return \code{\link{tbl_df}}
#' @noRd
read_LIMS_binder_file <- function(fn, ver) {
  colClasses <- c(rep("factor", 2), "integer", "character",
                  rep("factor", switch(ver, 3, 5)))
  fn %>% 
    read_file_lims(na.strings= "", colClasses = colClasses) %>% 
    as_tibble() %>% 
    rename_all(. %>% tolower() %>% gsub(" ", "_", .)) %>% 
    # > BAf id < #
    rename(sba = "bead_array", id = "binder_name") %>% 
    dplyr::select(id, everything()) %>% {
      if(ver == 2) {
        mutate_at(., vars(multiple_target), 
                  factor, c(0, 1), c("single", "multi"))
      } else .
    }
}



# -----------------------------------------------------------------------------#
#' Parse a set of constituent files exported from LIMS for a single project. 
#' 
#' A set of files for a project can be obtained from "Project Summary" in LIMS.
#' This function reads the set and construct an object of BAf class parsing the
#' files. Because this function is so sensitive to the format of data, the
#' functionality of the current version of this function can limited to
#' present-day form (2014-03-04) only with limited backward compatibility.
#' Please, note that any change in LIMS can affect this function. Optionally, it
#' produces some plots in technical aspects.
#' 
#' @param SLid the sample layout ID used in LIMS, e.g. "SL0001". Currently it is
#'   the first six letters in all constituent files.
#' @param path the path to the files from LIMS
#' @param ... for the future
#' @param QC_plot whether the plots in technical aspects should be drawn, such
#'   as intentisy distribution plots across samples or antibodies.
#' @param QC_report whether a primary QC report to be written invoking
#'   \code{\link{write_primary_QC_report.SBA}} or not
#' @param author the author of the report
#' @param path_QC where the QC plots or QC report to be saved
#' @param lowbound_mfi the lower bound of MFI. If a median of a sample is below
#'   this threshold, the data of the sample is neglected afterwords in QC
#'   report.
#' 
#' @return an object of the \code{\link{BAf-class}}
#' 
#' @seealso
#' \code{\link{read_FlexMAP3D_csv}}
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @include plot_QC.SBA.R
#' @export
# -----------------------------------------------------------------------------#
# created  : 2013-07-03 by Mun-Gwan
# modified : 
#   2013-09-16 by Mun-Gwan
#   2014-03-03 by Mun-Gwan : 
#     1) Accept new format of files in LIMS keeping backward compatibility
#     2) Allow to make plots of signal distribution in terms of samples and
#     binders
#   2014-11-21 by Mun-Gwan : bud out the plot part to "plot_QC"
#   2015-05-27 by Mun-Gwan : fix the problem to results in a vector instead of
#     a data frame when additional information have one column
#   2017-08-17 by Mun-Gwan : 
#     1) to adapt to BAf-class
#     2) Remove "additional_sampleInfo"
#   2017-08-17 by Mun-Gwan : to adapt to BAf-class
#   2018-04-11 by Mun-Gwan : change format of 'plate' from '1' to 'SL0001-1'
# -----------------------------------------------------------------------------#

read_LIMS_SBA_files <- function(SLid, 
                                path, 
                                ..., 
                                QC_plot = FALSE, 
                                QC_report = FALSE,
                                author = "Mun-Gwan",
                                path_QC = "../QC",
                                lowbound_mfi = 100) {
  stopifnot(!missing(path))
  stopifnot(!missing(SLid))
  
  if((QC_plot || QC_report) && ! file.exists(path_QC)) {
    warning("The 'path_QC'= ", path_QC, 
            " doesn't exist. So, it has been changed to current folder.")
    path_QC <- "."
  }

  ##  Find file names from LIMS ------------------------------------------------
  
  return_file_name_if_exists <- function(type, SLid_= SLid, path_= path) {
    fn <- SLid_ %>%
      paste0("_", type, ".tsv") %>%         # fixed file names
      { file.path(path_, .) }               # add path
    f_gz <- paste(fn, "gz", sep=".")
    
    if(file.exists(f_gz)) f_gz
    else if(file.exists(fn)) fn
    else NULL
  }
  
  fns <- sapply(
    c("data_intensity",
      "data_bead_count",
      "sample",
      "binder",
      "assay"
    ), function(type) {
      out <- return_file_name_if_exists(type)
      if(is.null(out)) 
        stop("The file ", type, ".tsv or (file).tsv.gz is required.")
      out
    }, 
    simplify= FALSE
  )
  
  names(fns)[names(fns) == "data_intensity"] <- "data"
  
  
  ##  Read 'data' file  ########################################################
  tgD <- read_file_lims(fns$data)
  
  # > measurement values < #
  X <- tgD[, 2:ncol(tgD)] %>% 
    as.matrix() %>% 
    `colnames<-`(names(tgD)[2:ncol(tgD)])
  
  
  
  ##  Read ‘Sample’ file  ######################################################
  ##  - Get header line to find the version of the sample file from LIMS
  for(ii in 1:100) {
    header <- scan(con <- gzfile(fns$sample), what= character(), nlines= ii, 
                   sep= "\t", comment.char= "#", quiet= TRUE)
    close(con)
    if(length(header) != 0) break;
  }
  lth_h <- length(header)
  
  #  header of sample file versions
  s_vers <- list(ver1 = c("Sample id", "Sample pos", "Plate pos", "Tube label", 
                          "Class", "Sample type", "Gender", "Sampling age")
  ) %>% c(list(ver2 = append(.$ver1, "Subtype", 5)))
  #  colClasses of different versions of sample files
  s_colClasses <- list(
    ver1 = c("character", "integer", 
             rep(c("character", "factor"), c(2, 3)), "numeric")
  ) %>% c(list(ver2 = append(.$ver1, "factor", 5))) %>% 
    # add NA for other variables
    lapply(. %>% append(., rep(NA, lth_h - length(.))))
  
  verS <- if(lth_h > 8 && identical(header[1:9], s_vers$ver2)) 2
  else if(identical(header[1:8], s_vers$ver1)) 1
  else {
    stop(paste("The file", fns$sample, "doesn't seem to be in LIMS format.\n", 
               "The pre-defined names in the first few columns are as belows\n",
               paste("1)", s_vers$ver1, sep= "\t"), "\n",
               paste("2)", s_vers$ver2, sep= "\t")))
  }
  tgS <- read_file_lims(
    fns$sample, 
    na.strings = c("NA", ""),
    colClasses = s_colClasses[[verS]]
  ) %>%
    as_tibble()
  rm(s_vers, s_colClasses)
  
  if(!identical(tgD[,1], tgS[["Sample id"]])) {
    stop("The sample IDs in 'data' and those in 'sample' must be identical.")
  }
  
  
  ##  STANDARDIZE the VARIABLE names and change to be CATEGORICAL --------------
  
  tolower_lvl <- . %>% `levels<-`(., levels(.) %>% tolower)
  
  tgS <- tgS %>% 
    rename(
      id    = "Sample id",
      sample_pos = "Sample pos",
      pos    = "Plate pos",
      label  = "Tube label",
      dis    = "Class",
      dis_sub = "Subtype",
      type   = "Sample type",
      gender = "Gender",
      age    = "Sampling age"
    ) %>% 
    `names<-`(make.unique(names(.))) %>%           # when more than std columns were read
    ## disease status
    mutate(dis = tolower_lvl(dis) %>% 
             (function(x, from, to) {
               levels(x)[levels(x) == from] <- to
               x
             })("disease", "case"),
           # sample type
           type = tolower_lvl(type),
           # gender
           gender = gender %>% 
             as.character() %>% 
             factor(c("F","M"), c("female", "male"))
    ) %>% 
    # "Cohort number" extracted form 'Sample id'
    `names<-`(make.unique(c("cohort", names(.)))[-1]) %>% # for the case 'cohort' already exists
    mutate(cohort = substr(id, 1, 5) %>% factor) %>% 
    dplyr::select(cohort, everything())
  
  ##  break down the plate position into "Plate number" and "well position"
  tmp <- tgS$pos %>% strsplit(":") %>% sapply(unlist) %>% t()
  tgS <- tgS %>% 
    #  for the case 'plate' already exists
    `names<-`(make.unique(c("plate", names(.)))[-1]) %>%
    #  for the case 'slayout' already exists
    `names<-`(make.unique(c("slayout", names(.)))[-1]) %>%
    mutate(plate = tmp[, 1] %>% paste0(SLid, "-", .) %>% factor, 
           pos   = tmp[, 2] %>% as.integer,
           id = make.unique_key_ids(id),
           slayout = factor(SLid)) %>% 
    dplyr::select(
      id, cohort, label, dis:age, plate, pos, slayout, sample_pos, everything()
    )
  
  ##  Read ‘Binder’ file  ######################################################
  
  tgB <- read_LIMS_binder_file(fns$binder, verS)
  stopifnot(all.equal(tgB$id, colnames(X)))
  tgB$id <- make.unique_key_ids(tgB$id)

  ##  Read ‘Bead Count’ file  ##################################################
  
  if("data_bead_count" %in% names(fns)) {
    tgBC <- read_file_lims(fns$data_bead_count) %>% {
      colnames(.)[-1] <- make.unique_key_ids(names(.)[-1])
      .
    } %>% 
      `rownames<-`(make.unique_key_ids(.[, 1])) %>% 
      .[, -1] %>% 
      as.matrix()

    if(!all.equal(rownames(tgBC), tgS$id)) {
      stop("The sample IDs in 'data_bead_count' and 'sample'", 
           " must be identical.")
    }
    stopifnot(all.equal(colnames(tgBC), tgB$id))   # binder IDs
    tgB[["median_bead_count"]] <- apply(tgBC, 2, median, na.rm= T)
  }
  
  
  ##  Read ‘Assay’ file  #######################################################
  
  tgA <- read_file_lims(fns$assay, na.strings= "", 
                        colClasses = c("factor", "integer", "numeric"))
  
  tgA <- sapply(tgA[, -1, drop= F], . %>% 
                  by(tgA[, 1], function(x) x) %>%
                  do.call("data.frame", .)
                , simplify= F)
  
  # identify the failed or too low bead count sample 
  tgA$fail_flag <- tgA$median_bead_count %>% 
    mutate_all(. %>% cut(c(0, BEADCOUNT_FAILED, BEADCOUNT_TOO_LOW, Inf), 
                         c("failed", "lowBeadCount", "ok"), right= F) 
    )
  
  
  ##  Generate Codebook --------------------------------------------------------
  
  
  std_cbook_sinfo <- tribble(
    ~column  , ~description,
    "id"     , "Sample ID",
    "cohort" , "In-house ID for a set of sample",
    "label"  , "Label for the sample. often ID from origianl DB",
    "dis"    , "disease status",
    "dis_sub", "subtype of disease",
    "type"   , "sample preparation type",
    "gender" , "gender of subject", 
    "age"    , "age at sample [yrs]",
    "plate"  , "the plate number",
    "pos"    , "well position number on the plate",
    "sample_pos", "#Not Used#"
  ) %>% 
  { if(verS == 1) .[-5, ] else . }
  
  std_cbook_binder <- tribble(
    ~column  , ~description,
    "id"     , "Binder ID",
    "assay"  , "Assay ID (1 assay = 1 SL and 1 SBA)",
    "sba"    , paste("SBA ID (a set of antibodies coupled with beads", 
                     "at the same time)"),
    "bead_pos" , "Bead color code",
    "prest_id" , "PrEST ID",
    "ensg_id"  , "Ensembl Gene ID",
    "gene_name", "Gene name (or Gene symbol)", 
    "gene_description", "Gene description",
    "multiple_target" , "If the antibody target single or multiple proteins"
  )
  
  ##  Create BAf object  #######################################################
  
  TG <- BAf(X, 
            sinfo = tgS,
            binder = tgB,
            sinfo_batch_i  = "plate", 
            binder_batch_i = "assay",
            assay_sinfo  = tgA,
            codebook_sinfo  = std_cbook_sinfo, 
            codebook_binder = std_cbook_binder 
  )
  
  ##  >> Technical plots << ----------------------------------------------------
  
  sample_order= order(sI(TG)$plate, sI(TG)$pos)
  
  if(! "data_bead_count" %in% names(fns)) tgBC <- NULL
  
  if(QC_plot) {
    ask_save.plot_QC.SBA(
      baf = TG,
      prefix = SLid,
      path = path_QC,
      bead_count = tgBC,
      sample_order= sample_order,
      color_tbl_b = default_color_SBA_bead(T),
      color_tbl_s = default_color_SBA_ctrl_samples(sid(TG)),
      main = paste(unique(batch(TG, "binder")), collapse= ", "),
      ylab_signal = "Raw MFI",
      guide_line = TRUE
    )
  }
  if(QC_report) {
    assy <- batch(TG, "binder")
    for(ii in unique(assy)) {
      write_primary_QC_report.SBA(
        baf= TG[, assy == ii], 
        bead_count = tgBC[, assy == ii],
        sample_order = sample_order,
        file = file.path(path_QC, 
                         paste0(SLid, ".", ii, ".primary_QC_report.html")),
        title= paste("Primary QC report of", SLid, "-", ii), 
        author= author,
        date= Sys.Date(),
        color_tbl_b = default_color_SBA_bead(T),
        color_tbl_s = default_color_SBA_ctrl_samples(sid(TG)),
        lowbound_mfi = lowbound_mfi
      ) 
    }
  }
  
  return(invisible(TG))
}


