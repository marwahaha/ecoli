library(dplyr)
load("/home/kunal/Downloads/ChemFieldLab2018(2).rdata")
#The path Angie needs to upload her data file
#WHEN CONNECTED TO VPN AND PATH NOT SEEN
#load("//pinwebserver01/C$/Shares/Penobscot Indian Nation Environmental Database/Reports/AnnualBaseline/2018/ChemFieldLab2018.RData")
#WHEN CONNECTED TO SERVER
#load("//pinwebserver01/Shares/Penobscot Indian Nation Environmental Database/Reports/AnnualBaseline/2018/ChemFieldLab2018.RData")


# Get all "Regular" Ecoli samples (ignoring "Split" and "Duplicate")
ecoli_all <- ChemFieldLab2018 %>%
  filter(Constituent =="E. coli (MPN)") %>%
  filter(SampleType == "Regular")


# The threshold values come from Ecoli_maxCol and Ecoli_GeoMaxCol
ecoli_without_threshold <- ecoli_all %>% filter(is.na(Ecoli_maxCol) | is.na(Ecoli_GeoMaxCol))
# FOR ANGIE: WHAT DO WE DO WHEN THERE IS NO THRESHOLD?
ecoli_with_threshold <- ecoli_all %>% filter(!is.na(Ecoli_maxCol) & !is.na(Ecoli_GeoMaxCol))


analyze_90_day_window <- function(ecoli_df, start, end) {
  # Get data in 90-day window
  ecoli_90 <- ecoli_df %>%
    filter(Date > start) %>%
    filter(Date <= end)

  # Summarize and check thresholds
  ecoli_summary <- ecoli_90 %>%
    group_by(Year, WaterBodyReport, SiteCode, SiteClassification, Ecoli_maxCol, Ecoli_GeoMaxCol)%>%
    summarise(
      NumResults = n(),
      NumResultsTenPercent = 0.1*n(),
      NumResultsAboveTenPercentThreshold = sum(c(ResultValue > first(Ecoli_maxCol))),
      RuleTenPercent = NumResultsAboveTenPercentThreshold <= NumResultsTenPercent,
      geometric_mean = exp(mean(log(ResultValue))),
      RuleGeometricMean= geometric_mean <= first(Ecoli_GeoMaxCol),
      Passing= RuleTenPercent & RuleGeometricMean,
      PassInst=RuleTenPercent,# Which sites pass the instantaneous rule
      PassGM=RuleGeometricMean,# Which sites pass the geo mean rule
    )

  return(ecoli_summary)
}


#FOR KUNAL: CAN WE SEPARATE OUT THE PASSING INSTANTANEOUS AND GEO MEAN?
get_ecoli_results <- function(ecoli_df) {
  # The initial 90-day window ends on the date of the last sample.
  # The window moves back by 1 day each time until it includes the date of the first sample.
  ecoli_full_results <- c() 
  end <- max(ecoli_df$Date)
  start <- end - 90
  while (start >= (min(ecoli_df$Date) - 1)) {
    # Decide if rules pass for this 90-day window, and put it in full results.
    ecoli_summary <- analyze_90_day_window(ecoli_df, start, end)

    ecoli_result <- ecoli_summary %>%
      select(Year, WaterBodyReport, SiteCode, SiteClassification, Ecoli_maxCol, Ecoli_GeoMaxCol, Passing)
    ecoli_full_results <- rbind(ecoli_full_results,  ecoli_result)
    
    # Move the 90-day window back by one day.              
    start <-start - 1
    end <- end - 1
  }
  
  # Check if each area passed in all 90-day windows.
  ecoli_output <- ecoli_full_results %>%
    group_by(Year, WaterBodyReport, SiteCode, SiteClassification, Ecoli_maxCol, Ecoli_GeoMaxCol)%>%
    summarize(success=all(Passing))
  
  return(ecoli_output)
}

out <- get_ecoli_results(ecoli_all)

failing <- out %>% filter(success != TRUE)
passing <- out %>% filter(success == TRUE)
unknown <- out %>% filter(is.na(success))

# Plot and inspect failing sites
library(ggplot2)
failing_with_data <- merge(ecoli_all, failing)
failing_with_data$title <- paste("site=", failing_with_data$SiteCode, 
                                 " inst=", failing_with_data$Ecoli_maxCol, 
                                 " geo=", failing_with_data$Ecoli_GeoMaxCol,
                                 sep="")

# after the image loads, press "Zoom" to see the whole thing
ggplot(failing_with_data ,
       mapping = aes(x=SampleDatetime,y=ResultValue, color=SiteCode)) +
  geom_point() +
  facet_wrap(facets = vars(title)) +
 theme(legend.position = "none")

# other things to try:
# Class B and C should only compute 90-day windows between April 15th and October 31s
# * facet_grid
# * Change 90-day window range and see what happens
