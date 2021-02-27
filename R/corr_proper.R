#' corr_proper
#' @description chooses proper correlation method (pearson or spearman) based on normality of the tested variables.
#' @param data  a dataframe or tibble.
#' @param columns can be either index or colnames of data, which enables subsetting data to test.
#' @param p.adjust p.adjust receive a string of holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", which would conduct p value adjust.
#'
#' @return    a list. contains the matrix and long form result of correlation matrix:\itemize{
#' \item\code{r:} the symmetric matrix of correlations
#' \item\code{p:} two tailed probability of t for each correlation
#' \item\code{n:} the number of sample used in this pair of correlation
#' \item\code{method:} the correlation method used in this pair of correlation analysis
#' \item\code{t:} value of t-test for each correlation
#' \item\code{se:} standard error of the correlation
#' \item\code{cor_long:} the long format of correlation analysis results}
#' @export
#'
#' @examples
#' data("mtcars")
#'
#' #Use all columns to compute correlaton
#' corr_proper(mtcars)
#'
#' #Use only the first 3 columns,and use fdr correction
#' corr_proper(mtcars, 1:3, 'fdr')
#'
corr_proper <- function(data, columns = 1:length(data), p.adjust = NULL){
  #Sub set data
  data_test <- data[columns]
  #Normality test
  normal <- lapply(data_test, shapiro.test)
  #Prepare result container
  table <- data.frame(row.names = colnames(data_test))
  results <- list(r = table,
                  p = table,
                  n = table,
                  method = table,
                  t = table,
                  se = table
                  )
  #Loop for paire-wise variable
  for (i in 1:length(data_test)) {
    for (j in 1:length(data_test)) {
      # If either variable does not meet normality
      if (normal[[i]]$p.value < 0.05 | normal[[j]]$p.value < 0.05 ) {
        #Corr test with spearman correlation
        cor <- psych::corr.test(data_test[i], data_test[j], method = 'spearman')
        #Record test method
        results$method[i,colnames(data_test)[j]] <- 'spearman'
      } else {
        #If both variables meet nomorlity, use pearson correlation
        cor <- psych::corr.test(data_test[i],data_test[j])
        results$method[i,colnames(data_test)[j]] <- 'pearson'
      }
      #Record results
      results$r[i,colnames(data_test)[j]] <- cor$r[1,1]
      results$p[i,colnames(data_test)[j]] <- cor$p[1,1]
      results$n[i,colnames(data_test)[j]] <- cor$n
      results$t[i,colnames(data_test)[j]] <- cor$t[1,1]
      results$se[i,colnames(data_test)[j]] <- cor$se[1,1]
    }
  }
  #Combine results into long format
  #Define a function to convert long
  cor_long <- function(x, value){
   x %>%
      dplyr::mutate(rownames = rownames(.)) %>%
      tibble::as_tibble() %>%
      dplyr::select(rownames,everything()) %>%
      tidyr::pivot_longer(-1, names_to = 'var2', values_to = value) %>%
      dplyr::rename(var1 = rownames)
  }
  #Combine corr test results in long format
  results[['cor_long']] <- suppressMessages(cor_long(results$r,'r') %>%
    dplyr::filter(var1 != var2) %>%
    dplyr::filter(!duplicated(r)) %>%
    dplyr::left_join(cor_long(results$p,'p')) %>%
    dplyr::left_join(cor_long(results$n,'n')) %>%
    dplyr::left_join(cor_long(results$method,'method')) %>%
    dplyr::left_join(cor_long(results$t,'t')) %>%
    dplyr::left_join(cor_long(results$se,'se')))
  #Conduct p adjust if p.adjust is provided
  if (!is.null(p.adjust)) {
    results$cor_long <- results$cor_long %>%
      dplyr::mutate(p.adj = p.adjust(p, method = p.adjust))
  }

  return(results)
}

