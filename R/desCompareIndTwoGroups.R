#' Describe and Compare Two Independent Groups
#' @description gives a easy-to-read result of two independent group description and comparison analysis.
#' @importFrom dplyr filter select across all_of pull
#' @param data a data.frame or tibble
#' @param columns dependent variables needed to be test, could be either string or column index
#' @param group the grouping variable, could be either string or column index
#' @param ci  a logical. Determine whether the ci of effect size should be calculated. Default is False
#' @param n The number of replications to use for bootstrap to Calculate ci of effect size, default is 1000 (larger number will be very time-consuming).
#'
#' @details The comparison is based on the normality and variance quality. If the data don't meet normality, Mann-Whitney test will be conducted. Otherwise the student t-test would be conducted.
#'
#' @return a data frame that contains description and comparison results.\itemize{
#'  \item\code{y:} The dependent variable
#'  \item\code{group:} The name of group
#'  \item\code{n:} The number of sample in group used in the test
#'  \item\code{mean:} mean of tested dependent variable
#'  \item\code{sd:} standard deviation of tested dependent variable
#'  \item\code{median:} median of tested dependent variable
#'  \item\code{groupx_q1:} first quantile of tested dependent variable
#'  \item\code{groupx_q3:} third quantile of tested dependent variable
#'  \item\code{mad:} median absolute deviation (see ?MAD)
#'  \item\code{effect.size:} if method is t-test, calculate cohen's d effect size.If the t-test did not make a homogeneity of variance assumption, (the Welch test), the variance term will mirror the Welch test, otherwise a pooled estimate is used. Interpretation: 0.2 (small effect), 0.5 (moderate effect) and 0.8 (large effect). If method is Mann-Whitney,Wilcoxon effect size (r) will be calculated. Interpretation: 0.10 - < 0.3 (small effect), 0.30 - < 0.5 (moderate effect) and >= 0.5 (large effect).
#'  \item\code{effect.conf.low,effect.conf.high:} lower and upper bound of the effect size confidence interval.}
#'
#' @export
#'
#' @examples
#' data("mtcars")
#' desCompareIndTwoGroups(data = mtcars, columns = -8, group = 8)
#'
#'
desCompareIndTwoGroups  <- function(data, columns, group, ci = FALSE, n = 1000){
  #Turn numerical group into character
  if (is.numeric(group)) {
    group = colnames(data)[group]
  }

  # Prepare result container
  results <- data.frame(row.names = 1:length(data[columns]))

  # Define dependent variables
  var <- colnames(data)[columns]
  # Set progress bar
  pb <- progress::progress_bar$new(total = length(var),
                                   format = 'Calculating[:bar]:percent time used: :elapsed eta: :eta')
  # Loop for dependent variables
  for(i in 1:length(var)){
    #Start the bar
    pb$tick(0)
    # Get dependent variable
    dependent = var[i]

    # Record dependent variable
    results[i,'y'] <- dependent

    # Generate test formula
    formula = formula(paste(dependent, ' ~ ', group))

    # Filter NA values in dependent variable
    data_test <- data %>%
      filter(across(all_of(dependent), ~ !is.na(.x)))

    # Get current groups of dependent variable
    two_groups <- pull(unique(data_test[group]))

    #If there are data in both groups
    if(length(two_groups) == 2){

      #Normality test
      normal <- lapply(two_groups,function(x){
        normal_data <- data_test[data_test[group] == x,] %>% pull()
        return(rstatix::shapiro_test(normal_data))
      })

      #Equation of variance test
      vareq <- suppressWarnings(car::leveneTest(y = pull(data[dependent]), group = pull(data[group]), data = data_test)[1,3])

      #Prepare data for describe analysis
      data_test %>%
        select(all_of(c(dependent,group))) -> desc_data
      #Describe table of group1
      desc_data[desc_data[group]==two_groups[1],] %>%
        select(all_of(dependent)) %>% rstatix::get_summary_stats() -> group1_desc
      #Describe table of group2
      desc_data[desc_data[group]==two_groups[2],] %>%
        select(all_of(dependent)) %>% rstatix::get_summary_stats() -> group2_desc

      #Record describe analysis results
      results[i,'group1'] <- two_groups[1]
      results[i,'group2'] <- two_groups[2]
      results[i,'n1'] <- group1_desc$n
      results[i,'n2'] <- group2_desc$n
      results[i,'mean1'] <- group1_desc$mean
      results[i,'mean2'] <- group2_desc$mean
      results[i,'sd1'] <- group1_desc$sd
      results[i,'sd2'] <- group2_desc$sd
      results[i,'median1'] <- group1_desc$median
      results[i,'median2'] <- group2_desc$median
      results[i,'group1_q1'] <- group1_desc$q1
      results[i,'group2_q1'] <- group2_desc$q1
      results[i,'group1_q3'] <- group1_desc$q3
      results[i,'group2_q3'] <- group2_desc$q3
      results[i,'mad1'] <- group1_desc$mad
      results[i,'mad2'] <- group2_desc$mad

      #If group's does not meet normality
      if(normal[[1]]$p.value < 0.05 | normal[[2]]$p.value < 0.05){
        #Record normality result
        results[i,'normal'] <- 'False'
        #Record test method
        results[i,'method'] <- 'Mann-Whitney'
        #Set variance equation as NA
        results[i,'varequ'] <- NA
        #Mann-Whitney test
        model <- rstatix::wilcox_test(data = data_test, formula)
        #Calculate the effect size
        effect_size <- rstatix::wilcox_effsize(data = data_test, formula = formula, ci = ci, nboot = n)

        #Record the test Results
        results[i,'statistic'] <- model$statistic
        results[i,'p'] <- model$p
        results[i,'effect.size'] <- effect_size$effsize
        if (ci == TRUE) {
          results[i,'effect.conf.low'] <- effect_size$conf.low
          results[i,'effect.conf.high'] <- effect_size$conf.high
        }
        #Else means suitable for t-test
      } else {
        #Record normality test result
        results[i,'normal'] <- 'True'
        #Record test method
        results[i,'method'] <- 't-test'
        #Test the variance equation
        if(vareq < 0.05){
          #The variances doesn't qual in two groups
          results[i,'varequ'] <- 'False'
          #Let var.equal as F
          model <- rstatix::t_test(data = data_test, formula, var.equal = F)
          effect_size <- rstatix::cohens_d(data = data_test, formula = formula, var.equal = F, ci = T, nboot = n)

          results[i,'statistic'] <- model$statistic
          results[i,'p'] <- model$p
          results[i,'effect.size'] <- effect_size$effsize
          if (ci == TRUE) {
            results[i,'effect.conf.low']  <- effect_size$conf.low
            results[i,'effect.conf.high'] <- effect_size$conf.high
          }
        } else{
          #Else means the variances are equal
          results[i,'varequ'] <- 'True'
          model <- rstatix::t_test(data = data_test, formula, var.equal = T)
          effect_size <- rstatix::cohens_d(data = data_test, formula = formula, var.equal = T, ci = T, nboot = n)
          results[i,'statistic'] <- model$statistic
          results[i,'p'] <- model$p
          results[i,'effect.size'] <- effect_size$effsize
          if (ci == TRUE) {
            results[i,'effect.conf.low'] <- effect_size$conf.low
            results[i,'effect.conf.high'] <- effect_size$conf.high
          }
        }
      }
      #This branch solves the situation that only one group has data while the other doesn't
    } else if (length(two_groups) == 1) {
      #Find out which group has data
      index <- which(two_groups == pull(unique(data[group])))
      #Prepare data for describe analysis
      data_test %>%
        select(all_of(c(dependent,group))) -> desc_data
      #If group1 has data
      if (index == 1) {
        desc_data[desc_data[group]==two_groups,] %>%
          select(all_of(dependent)) %>% rstatix::get_summary_stats() -> group_desc
        #Record describe analysis results
        results[i,'group1'] <- two_groups
        results[i,'n1'] <- group_desc$n
        results[i,'mean1'] <- group_desc$mean
        results[i,'sd1'] <- group_desc$sd
        results[i,'median1'] <- group_desc$median
        results[i,'group1_q1'] <- group_desc$q1
        results[i,'group1_q3'] <- group_desc$q3
        results[i,'mad1'] <- group_desc$mad
      } else if (index == 2){
        #If group2 has data
        desc_data[desc_data[group]==two_groups,] %>%
          select(all_of(dependent)) %>% rstatix::get_summary_stats() -> group_desc
        #Record describe analysis results
        results[i,'group2'] <- two_groups
        results[i,'n2'] <- group_desc$n
        results[i,'mean2'] <- group_desc$mean
        results[i,'sd2'] <- group_desc$sd
        results[i,'median2'] <- group_desc$median
        results[i,'group2_q1'] <- group_desc$q1
        results[i,'group2_q3'] <- group_desc$q3
        results[i,'mad2'] <- group_desc$mad
      }
    }
    pb$tick()
  }
  return(results)
}




