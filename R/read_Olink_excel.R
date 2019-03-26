# -----------------------------------------------------------------------------#
#' Read Olink Excel
#' 
#' This generates an BAf object parsing an output file in Excel format from
#' Olink instrument
#'
#' @param file name of a file from Olink
#' 
#' @return an object of \link{BAf-class}
#'   
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @seealso
#' \code{\link{read_excel}}
#' @importFrom readxl read_excel
#' @export
# -----------------------------------------------------------------------------#
# created  : 2018-09-12 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

read_Olink_excel <- function(file) {
  rd0 <- readxl::read_excel(
    path= file,
    col_names= F,
    col_types= "text"
  )
  
  ##  find i,j 
  i_panel <- which(pull(rd0, 1) == "Panel")
  i_olinkid <- which(pull(rd0, 1) == "OlinkID")
  i_lod <- which(pull(rd0, 1) == "LOD")
  
  i_data <- c((i_olinkid + 1):(i_lod - 1)) %>% {
    .[!is.na(pull(rd0, 1)[.])]
  }
  j_data <- slice(rd0, i_olinkid) %>% 
    unlist(use.names= F) %>% {
      which(!is.na(.))
    }
  
  ##  Header
  rd_h <- rd0[c(i_panel:i_olinkid, i_lod:(i_lod + 1)), j_data] %>% 
    t() %>% 
    `colnames<-`(unlist(.[1, ])) %>% 
    `[`(-1, ) %>% 
    as.tibble %>% 
    separate(Panel, c("Panel", "Version"), sep= "\\(") %>% 
    mutate(
      Version = sub("\\)$", "", Version),
      LOD = as.numeric(LOD),
      `Missing Data freq.` = as.numeric(`Missing Data freq.`)
    ) %>% 
    select(id = OlinkID, everything())
  
  stopifnot(anyDuplicated(rd_h$id) == 0)
  
  
  ##  Numeric NPX Data
  rd_X <- rd0[i_data, j_data] %>% 
    `names<-`(unlist(rd0[i_olinkid, j_data])) %>% 
    mutate_at(-1, as.numeric) %>% 
    as.data.frame() %>% 
    `rownames<-`(make.unique_key_ids(.$OlinkID)) %>% 
    select(-OlinkID) %>% 
    as.matrix()
  
  
  ##  Plate ID
  j_plate <- which(unlist(rd0[i_panel + 1, ], use.names= F) == "Plate ID")
  rd_plate <- rd0[i_data, j_plate] %>% 
    `names<-`(unlist(slice(rd0, i_panel)[j_plate])) %>% 
    as.data.frame() %>% 
    `rownames<-`(make.unique_key_ids(unlist(rd0[i_data, 1])))
  
  ##  sample Info
  concat_plate <- rd_plate %>% 
    unite("concat", everything()) %>% 
    pull(concat)
  
  sinfo <- tibble(
    id= rownames(rd_X),
    plate= concat_plate
  )
  
  ##  QC warning
  j_warn <- which(unlist(rd0[i_panel + 1, ], use.names= F) == "QC Warning")
  rd_warning <- rd0[i_data, j_warn] %>% 
    `names<-`(unlist(slice(rd0, i_panel)[j_warn])) %>% 
    as.data.frame() %>% 
    `rownames<-`(make.unique_key_ids(unlist(rd0[i_data, 1])))
  
  BAf(
    rd_X, 
    sinfo = sinfo,
    binder = rd_h,
    sinfo_batch_i  = "plate", 
    binder_batch_i = "Panel",
    assay_sinfo  = list(
      "Plate ID"= rd_plate,
      "QC warning"= rd_warning
    )
  )
}
