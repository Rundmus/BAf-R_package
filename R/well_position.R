# -----------------------------------------------------------------------------#
#' Conversion of well position index
#' 
#' The functions to convert the index of well position in 96 or 384-well plate
#' \code{well_rowcol} - Convert position number to letter format.
#' 
#' @param x a vector of position numbers
#' @param row if the letters for rows should be returned
#' @param col if the column number should be returned
#' @param as.factor if the returned vector should be converted to factor
#' 
#' @return 
#' \code{well_rowcol} produces a character vector that has well position letter
#' 
#' @examples 
#'   well_rowcol(sample(1:96, 10))
#' 
#' @name well_position
#' @rdname well_position
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#
# created  : 2012-01-24 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

well_rowcol <- function(x, row = TRUE, col = TRUE, as.factor = FALSE) {
  if(! (row | col)) 
    stop("At least, one of the \'row\' or the \'col\' should be TRUE.")
  stopifnot(all(x <= 96 & x >= 1, na.rm= T))
  
  lth <- length(x)
  
  if(row) {
    tg <- LETTERS[((x - 1) %% 8) + 1]
    if(col) tg <- paste0(tg, formatC((x - 1) %/% 8 + 1, width= 2, flag= "0"))
  } else {
    tg <- (x - 1) %/% 8 + 1
  }
  
  if( as.factor ) tg <- factor(tg)
  return(tg)
}


# -----------------------------------------------------------------------------#
#' \code{well_pos} - Convert letter format to position number
#' 
#' @param rc a vector of position letters composed of row letter and column
#'   number e.g. B02
#' @param n_well the number of wells per plate. It should be 96 or 384. The
#'   default is 96.
#'   
#' @return 
#' \code{well_pos} an integer vector of position number
#' 
#' @examples 
#'     well_pos("B01")
#' 
#' @rdname well_position
#' @export
# -----------------------------------------------------------------------------#

well_pos <- function(rc, n_well = c(96, 384)) {
  if(length(n_well) > 1) n_well <- n_well[1]
  stopifnot(n_well %in% c(96, 384))
  n_rows <- if(n_well == 96) 8 else 16
  
  rc <- as.character(rc)
  
  row <- rc %>% 
    substr(1, 1) %>% 
    toupper() %>% 
    match(LETTERS)
  
  if(any(row > n_rows, na.rm= T)) 
    stop("An (or more) unappropriate row letter is detected!")
  col <- as.numeric(substring(rc, 2))
  
  (col - 1) * n_rows + row
}


# -----------------------------------------------------------------------------#
#' \code{well_row384} - Find row letter in 384 plate layout of the given
#' position and plate number
#' 
#' @param pos a vector of position number on 96-well plate
#' @param plate a vector of plate number of 96-well plates
#'   
#' @return 
#' \code{well_row384} a factor or a vector of characters
#' 
#' @rdname well_position
#' @export
# -----------------------------------------------------------------------------#

well_row384 <- function(pos, plate, as.factor = TRUE) {
  stopifnot(length(pos) == length(plate))
  as_num <- . %>% as.character() %>% as.numeric()
  row384 <- LETTERS[ ((as_num(pos) - 1) %% 8 + 1) * 2 - as_num(plate) %% 2 ]
  
  if( as.factor ) row384 <- factor(row384)
  return(row384)
}


# -----------------------------------------------------------------------------#
#' \code{well_col384} - Find column number in 384 plate layout of the given
#' position and plate number
#' 
#' @return 
#' \code{well_col384} a factor or a vector of characters
#' 
#' @rdname well_position
#' @export
# -----------------------------------------------------------------------------#

well_col384 <- function(pos, plate, as.factor = TRUE) {
  stopifnot(length(pos) == length(plate))
  as_num <- . %>% as.character() %>% as.numeric()
  col384 <- 
    (((as_num(plate) - 1) %% 4) - 2) %/% 2 + 2 * ((as_num(pos) - 1) %/% 8 + 1) 
  if( as.factor ) col384 <- factor(col384)
  return(col384)
}

