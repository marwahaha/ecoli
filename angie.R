library(dplyr)
load("/home/kunal/Downloads/ChemFieldLab2018(1).rdata")

# Get all Ecoli samples
# Question for Angie: Does sampleType matter?u Should we ignore either "Split" or "Duplicate"?
ecoli_all <- ChemFieldLab2018 %>%
  filter(Constituent =="E. coli (MPN)")


get_ecoli_results <- function(ecoli_all, TEN_PERCENT_THRESHOLD, GEOMETRIC_MEAN_THRESHOLD) {
  # Limiting to last 90 days, then moving the 90-day window back by 1 day each time.
  ecoli_full_results <- c() 
  end <- max(ecoli_all$Date)
  start <- end - 90
  while (start >= (min(ecoli_all$Date) - 1)) {
    # Get data in 90-day window
    ecoli_90 <- ecoli_all %>%
      filter(Date > start) %>%
      filter(Date <= end)
    
    # Note for Angie: There is an issue with TMAS data, with some rows having no SiteClassification
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
    
    # Decide if rules pass for this 90-day window, and put it in full results.
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
