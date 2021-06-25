# surv miner custom theme
#https://stackoverflow.com/questions/62995652/center-the-plot-title-in-ggsurvplot
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}