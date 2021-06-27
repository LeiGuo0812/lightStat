#' Correlation Analysis with Proper Method
#' @description chooses proper correlation method (pearson or spearman) based on normality of the tested variables.
#' @importFrom magrittr %>%
#' @importFrom dplyr pull mutate select across all_of filter rename left_join everything
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom psych corr.test
#' @importFrom stats shapiro.test
#' @param data  a dataframe or tibble.
#' @param y the second matrix or dataframe with the same number of rows as data.If y is provided, the correlation just between variables in data and y will be calculated. If y is not provided, pair-wise correlation between variables in data will be calculated.
#' @param p.adjust p.adjust receive a string of holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", which would conduct p value adjust.
#'
#' @return    a list. contains the matrix and long form result of correlation matrix:\itemize{
#' \item\code{r:} the symmetric matrix of correlations
#' \item\code{p:} two tailed probability of t for each correlation
#' \item\code{n:} the number of sample used in this pair of correlation
#' \item\code{method:} the correlation method used in this pair of correlation analysis. If both variables are ordinal, kendall correlation analysis will be used. If both variables are numeric and meet normality, pearson correlation analysis will be used. If either of the two variables doesn't meet normality, spearman correlation analysis will be used. If one variable is ordinal and the other is numerical ,spearman correlation analysis will be used. If either variable is factor, NA will be returned to all results.
#' \item\code{t:} value of t-test for each correlation
#' \item\code{se:} standard error of the correlation
#' \item\code{cor_long:} a tibble, the long format of correlation analysis results, NA results will not be shown}
#' @export
#'
#' @examples
#' data("mtcars")
#' corrProper(mtcars)
#'
#' #Use fdr correction
#' corrProper(mtcars, p.adjust = 'fdr')
#'
#' #Calculate correlation between two data frames
#' corrProper(mtcars[1:3],mtcars[4:5], p.adjust = 'fdr')
#'
corrProper <- function(data, y = NULL, p.adjust = NULL){
  if (is.null(y)) {
    #Prepare result container
    table <- data.frame(row.names = colnames(data))
    results <- list(r = table,
                    p = table,
                    n = table,
                    method = table,
                    t = table,
                    se = table
    )
    #Loop for paire-wise variable
    for (i in 1:length(data)) {
      for (j in 1:length(data)) {
        if(is.character(data[i] %>% pull()) |
           is.character(data[j] %>% pull())){
          cor <- list(r = NA,
                      p = NA,
                      n = NA,
                      t = NA,
                      se = NA)
          results$method[i,colnames(data)[j]] <- NA
        } else if (is.ordered(data[i] %>% pull()) &
                   is.ordered(data[j] %>% pull())) {
          cor <- corr.test(as.numeric(data[i] %>% pull()),
                                  as.numeric(data[j] %>% pull()),
                                  method = 'kendall')
          #Record test method
          results$method[i,colnames(data)[j]] <- 'kendall'
        }
        else if(
          (is.ordered(data[i] %>% pull()) &
           !is.factor(data[j] %>% pull())) |
          (!is.factor(data[i] %>% pull()) &
           is.ordered(data[j] %>% pull()))){
          cor <- corr.test(as.numeric(data[i] %>% pull()),
                                  as.numeric(data[j] %>% pull()),
                                  method = 'spearman')
          #Record test method
          results$method[i,colnames(data)[j]] <- 'spearman'
        } else if (is.factor(data[i] %>% pull()) |
                   is.factor(data[j] %>% pull())) {
          cor <- list(r = NA,
                      p = NA,
                      n = NA,
                      t = NA,
                      se = NA)
          results$method[i,colnames(data)[j]] <- NA
        }
        else {
          #return NA if any vector has a 0 sd.
          if (sd(pull(data[i])) == 0 |
              sd(pull(data[j])) == 0) {
            cor <- list(r = NA,
                        p = NA,
                        n = NA,
                        t = NA,
                        se = NA)
            results$method[i,colnames(data)[j]] <- NA
            warning('sd is 0')
          }
          else{
            #Normality test
            normal <- lapply(c(data[i],data[j]), shapiro.test)
            # If either variable does not meet normality
            if (normal[[1]]$p.value < 0.05 | normal[[2]]$p.value < 0.05 ) {
            #Corr test with spearman correlation
            cor <- corr.test(as.numeric(data[i] %>% pull()),
                                    as.numeric(data[j] %>% pull()),
                                    method = 'spearman')
            #Record test method
            results$method[i,colnames(data)[j]] <- 'spearman'
          } else {
            #If both variables meet normality, use pearson correlation
            cor <- corr.test(as.numeric(data[i] %>% pull()),
                                    as.numeric(data[j] %>% pull()))
            results$method[i,colnames(data)[j]] <- 'pearson'
          }
          }
        }
        #Record results
        results$r[i,colnames(data)[j]] <- cor$r
        results$p[i,colnames(data)[j]] <- cor$p
        results$n[i,colnames(data)[j]] <- cor$n
        results$t[i,colnames(data)[j]] <- cor$t
        results$se[i,colnames(data)[j]] <- cor$se
    }
    }
    } else {
    # If the second data frame exists
    #Prepare result container
    table <- data.frame(row.names = colnames(data))
    results <- list(r = table,
                    p = table,
                    n = table,
                    method = table,
                    t = table,
                    se = table
    )
    #Loop for paire-wise variable
    for (i in 1:length(data)) {
      for (j in 1:length(y)) {
        if(is.character(data[i] %>% pull()) |
           is.character(y[j] %>% pull())){
          cor <- list(r = NA,
                      p = NA,
                      n = NA,
                      t = NA,
                      se = NA)
          results$method[i,colnames(y)[j]] <- NA
        } else if (is.ordered(data[i] %>% pull()) &
                   is.ordered(y[j] %>% pull())) {
          cor <- corr.test(as.numeric(data[i] %>% pull()),
                                  as.numeric(y[j] %>% pull()),
                                  method = 'kendall')
          #Record test method
          results$method[i,colnames(y)[j]] <- 'kendall'
        }
        else if(
          (is.ordered(data[i] %>% pull()) &
           !is.factor(y[j] %>% pull())) |
          (!is.factor(data[i] %>% pull()) &
           is.ordered(y[j] %>% pull()))){
          cor <- corr.test(as.numeric(data[i] %>% pull()),
                                  as.numeric(y[j] %>% pull()),
                                  method = 'spearman')
          #Record test method
          results$method[i,colnames(y)[j]] <- 'spearman'
        } else if (is.factor(data[i] %>% pull()) |
                   is.factor(y[j] %>% pull())) {
          cor <- list(r = NA,
                      p = NA,
                      n = NA,
                      t = NA,
                      se = NA)
          results$method[i,colnames(y)[j]] <- NA
        }
        else {
          #return NA if any vector has a 0 sd.
          if (sd(pull(data[i])) == 0 |
              sd(pull(y[j])) == 0) {
            cor <- list(r = NA,
                        p = NA,
                        n = NA,
                        t = NA,
                        se = NA)
            results$method[i,colnames(data)[j]] <- NA
            warning('sd is 0')
          }
          else {
            #Normality test
            normal <- lapply(c(data[i],y[j]), shapiro.test)
            # If either variable does not meet normality
            if (normal[[1]]$p.value < 0.05 | normal[[2]]$p.value < 0.05 ) {
            #Corr test with spearman correlation
            cor <- corr.test(as.numeric(data[i] %>% pull()),
                                    as.numeric(y[j] %>% pull()),
                                    method = 'spearman')
            #Record test method
            results$method[i,colnames(y)[j]] <- 'spearman'
          } else {
            #If both variables meet nomorlity, use pearson correlation
            cor <- corr.test(as.numeric(data[i] %>% pull()),
                                    as.numeric(y[j] %>% pull()))
            results$method[i,colnames(y)[j]] <- 'pearson'
          }
          }
        }
        #Record results
        results$r[i,colnames(y)[j]] <- cor$r
        results$p[i,colnames(y)[j]] <- cor$p
        results$n[i,colnames(y)[j]] <- cor$n
        results$t[i,colnames(y)[j]] <- cor$t
        results$se[i,colnames(y)[j]] <- cor$se
    }
    }
    }
  #Combine results into long format
  #Define a function to convert long
  cor_long <- function(x, value){
    x %>%
      mutate(rownames = rownames(.)) %>%
      as_tibble() %>%
      select(rownames,everything()) %>%
      pivot_longer(-1, names_to = 'var2', values_to = value) %>%
      rename(var1 = rownames)
  }
  #Combine corr test results in long format
  results[['cor_long']] <- suppressMessages(cor_long(results$r,'r') %>%
                                              filter(var1 != var2) %>%
                                              filter(!duplicated(r)) %>%
                                              filter(!is.na(r)) %>%
                                              left_join(cor_long(results$p,'p')) %>%
                                              left_join(cor_long(results$n,'n')) %>%
                                              left_join(cor_long(results$method,'method')) %>%
                                              left_join(cor_long(results$t,'t')) %>%
                                              left_join(cor_long(results$se,'se')))

  #Conduct p adjust if p.adjust is provided
  if (!is.null(p.adjust)) {
    results$cor_long <- results$cor_long %>%
      mutate(p.adj = p.adjust(p, method = p.adjust))
  }

  return(results)
}

### At the day this package is created, 2021.2.27, one of my favorite actor, Menda Wu, passed away. Thanks to his great contribution to Chinese film, which brought so many joys to my childhood. R.I.P
