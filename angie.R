library(dplyr)
#load("/home/kunal/Downloads/ChemFieldLab2018(2).rdata")
#The path Angie needs to upload her data file
#WHEN CONNECTED TO VPN AND PATH NOT SEEN
load("//pinwebserver01/C$/Shares/Penobscot Indian Nation Environmental Database/Reports/AnnualBaseline/2018/ChemFieldLab2018.RData")


# Get all "Regular" Ecoli samples (ignoring "Split" and "Duplicate")
ecoli_all <- ChemFieldLab2018 %>%
  filter(Constituent =="E. coli (MPN)") %>%
  filter(SampleType == "Regular")


analyze_90_day_window <- function(ecoli_all, start, end, TEN_PERCENT_THRESHOLD, GEOMETRIC_MEAN_THRESHOLD) {
  # Get data in 90-day window
  ecoli_90 <- ecoli_all %>%
    filter(Date > start) %>%
    filter(Date <= end)

  # Summarize and check thresholds
  ecoli_summary <- ecoli_90 %>%
    group_by(Year, WaterBodyReport, SiteCode, SiteClassification)%>%
    summarise(
      NumResults = n(),
      NumResultsTenPercent = 0.1*n(),
      NumResultsAboveTenPercentThreshold = sum(c(ResultValue > TEN_PERCENT_THRESHOLD)),
      RuleTenPercent = NumResultsAboveTenPercentThreshold <= NumResultsTenPercent,
      geometric_mean = exp(mean(log(ResultValue))),
      RuleGeometricMean= geometric_mean <= GEOMETRIC_MEAN_THRESHOLD,
      Passing= RuleTenPercent&RuleGeometricMean,
    )

  return(ecoli_summary)
}


get_ecoli_results <- function(ecoli_all, TEN_PERCENT_THRESHOLD, GEOMETRIC_MEAN_THRESHOLD) {
  # The initial 90-day window ends on the date of the last sample.
  # The window moves back by 1 day each time until it includes the date of the first sample.
  ecoli_full_results <- c() 
  end <- max(ecoli_all$Date)
  start <- end - 90
  while (start >= (min(ecoli_all$Date) - 1)) {
    # Decide if rules pass for this 90-day window, and put it in full results.
    ecoli_summary <- analyze_90_day_window(ecoli_all, start, end, TEN_PERCENT_THRESHOLD, GEOMETRIC_MEAN_THRESHOLD)

    ecoli_result <- ecoli_summary %>%
      select(Year, WaterBodyReport, SiteCode, SiteClassification, Passing)
    ecoli_full_results <- rbind(ecoli_full_results,  ecoli_result)
    
    # Move the 90-day window back by one day.              
    start <-start - 1
    end <- end - 1
  }
  
  # Check if each area passed in all 90-day windows.
  ecoli_output <- ecoli_full_results %>%
    group_by(Year, WaterBodyReport, SiteCode, SiteClassification)%>%
    summarize(success=all(Passing))
  
  return(ecoli_output)
}

out <- get_ecoli_results(ecoli_all, 236, 64)

failing <- out %>% filter(success != TRUE)
passing <- out %>% filter(success == TRUE)

# Plot and inspect failing sites
library(ggplot2)
failing_with_data <- merge(ecoli_all, failing)

# after the image loads, press "Zoom" to see the whole thing
ggplot(failing_with_data ,
       mapping = aes(x=SampleDatetime,y=ResultValue, color=SiteCode)) +
  geom_point() +
  facet_wrap(facets = vars(SiteCode, WaterBodyReport)) +
 theme(legend.position = "none")

# other things to try:
# * facet_grid
# * Grouping by SiteClassification
# * Change thresholds and see what changes
# * Change 90-day window range and see what happens
