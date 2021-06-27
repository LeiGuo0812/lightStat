#' Detect continuous same choices in the inventory
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select
#' @param data A data frame. Each row represents an observation, and each column represents an item
#' @param k The minimum number of continuous same choices to detect.
#' @param index A vector of index of column. Used when there is only a subset of columns is to detected.
#' @param keep A logical indicating whether keep the information about the continuous same choices. If is set FALSE, the observations with same choices will be discarded and the information will not be kept. Default is TRUE.
#'
#' @return A data frame, a combination of the original data frame and two other columns:
#' \itemize{
#' \item{\code{identical}: An integer, indicates the number of  continuous same choices in this observation's answers.}
#' \item{\code{start_index}: The index of the first index of the detected continuous same choices.}
#' }
#' @export
#'
#' @examples
#' data('mtcars')
#' sameChoiceDetect(mtcars, k = 4)
#'
sameChoiceDetect <- function(data, k = NULL, index = NULL, keep = TRUE) {

  D <- as.data.frame(data)
  n = nrow(D)

  if (is.null(index)) {
    c = ncol(D)
  } else {
    D <- D[index]
    c = ncol(D)
  }

  if (is.null(k)) {
    k = ceiling(ncol(D)/3)
  }

  cat(paste('\n','The minimum number of continuous same choices to detect is',k,'\n\n'))

  D[,'identical'] <- 0
  D[,'start_index'] <- NA_character_

  for (i in seq_len(n)) {
    for (j in seq_len(c-k+1)) {
      if (sum(!duplicated(t(D[i,j:(j+k-1)]))) == 1) {
        D[i,'identical'] = D[i,'identical'] + 1
        D[i,'start_index'] =  ifelse(is.na(D[i,'start_index']),j,paste(D[i,'start_index'],j,sep = ';'))
      }
   }
  }
  D <- cbind(data, D[,c('identical','start_index')])

  if (!keep) {
    D <- D %>%
      dplyr::filter(identical == 0) %>%
      select(-identical, -start_index)
  }

  return(D)
}


