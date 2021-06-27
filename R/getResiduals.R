#' Get Residuals with Linear Model
#' @description get residuals of columns in data regressed cov with \code{\link[stats]{lm}}
#' @importFrom dplyr pull
#' @importFrom stats as.formula lm
#' @param data a data frame or tibble, contains the columns to calculate residuals
#' @param cov a data frame or tibble, contains the columns as covariate. The length of row should be the same as data
#' @param addmean logical, should the mean of original data should be added back to the residual, default is true
#'
#' @return a data frame with same dimension with data, containing the residuals of the original columns
#' @export
#'
#' @examples
#' data("mtcars")
#' getResiduals(mtcars[1:3], mtcars[4:5])
#'
#' #get residual without adding mean back
#' getResiduals(mtcars[1:3], mtcars[4:5],addmean = FALSE)
#'
getResiduals <- function(data,cov,addmean = TRUE){
  lm_data <- cbind(data,cov)
  var <- colnames(data)
  cov <- colnames(cov)
  results <- lapply(var, function(x){
    formula <- as.formula(paste(x,paste(cov,collapse = '+'),sep = '~'))
    model <- lm(formula = formula, data = lm_data)
    if (addmean == TRUE) {
      return(model$residuals + mean(pull(model$model[x]),na.rm = T))
    } else {
      return(model$residuals)
    }
  })
  names(results) <- paste0(var,'.residual')
  results %>% as.data.frame() -> results
  return(results)
}
