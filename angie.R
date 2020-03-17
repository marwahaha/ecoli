library(dplyr)
load("/home/kunal/Downloads/ChemFieldLab2018(2).rdata")
#The path Angie needs to upload her data file
#WHEN CONNECTED TO VPN AND PATH NOT SEEN
#load("//pinwebserver01/C$/Shares/Penobscot Indian Nation Environmental Database/Reports/AnnualBaseline/2018/ChemFieldLab2018.RData")
#WHEN CONNECTED TO SERVER
#load("//pinwebserver01/Shares/Penobscot Indian Nation Environmental Database/Reports/AnnualBaseline/2018/ChemFieldLab2018.RData")


# Get all "Regular" Ecoli samples (ignoring "Split" and "Duplicate")
# FOR ANGIE: SHOULD WE INCLUDE SPLIT AND DUPLICATE SAMPLES?
ecoli_all <- ChemFieldLab2018 %>%
  filter(Constituent =="E. coli (MPN)") %>%
  filter(SampleType == "Regular")


analyze_90_day_window <- function(ecoli_df, start) {
  # Get data in 90-day window
  ecoli_90 <- ecoli_df %>%
    filter(Date > start) %>%
    filter(Date <= (start + 90))

  # Summarize and check thresholds
  # RuleTenPercent is TRUE when the site passes the 10% rule
  # RuleGeometricMean is TRUE when the site passes the geo mean rule
  ecoli_summary <- ecoli_90 %>%
    group_by(Year, WaterBodyReport, SiteCode, SiteClassification, Ecoli_maxCol, Ecoli_GeoMaxCol)%>%
    summarise(
      NumSamples = n(),
      # Ten percent rule
      TenPercentOfSamples = 0.1*n(),
      NumSamplesAboveTenPercentThreshold = sum(c(ResultValue > unique(Ecoli_maxCol))),
      FailingTenPercentRule = NumSamplesAboveTenPercentThreshold > TenPercentOfSamples,
      # Geo mean rule
      geometric_mean = exp(mean(log(ResultValue))),
      EnoughSamplesForGeoMean = (n() >= 6), # Only do Geo Mean rule if at least 6 samples
      FailingGeoMeanRule = geometric_mean > unique(Ecoli_GeoMaxCol),
      Failing = FailingTenPercentRule | (EnoughSamplesForGeoMean & FailingGeoMeanRule),
    )

  return(ecoli_summary)
}


get_ecoli_results <- function(ecoli_df) {
  
  # Earliest window starts 90 days before the first measurement.
  # Latest window starts date of the last measurement.
  earliest_window_start <- min(ecoli_df$Date) - 90
  latest_window_start <- max(ecoli_df$Date) 
  # If you change the above, it makes a big difference in the analysis!
  # One way to change it: 
  # * earliest_window_start <- min(ecoli_df$Date)
  # * latest_window_start <- max(ecoli_df$Date) - 90
  
  ecoli_full_results <- c() 
  start <- earliest_window_start
  while (start <= latest_window_start) {
    # Decide if rules pass for this 90-day window, and put it in full results.
    ecoli_summary <- analyze_90_day_window(ecoli_df, start)
    ecoli_full_results <- rbind(ecoli_full_results,  ecoli_summary)
    # Move the 90-day window forward by one day.              
    start <- start + 1
  }
  
  return(ecoli_full_results)
}

full_results <- get_ecoli_results(ecoli_all)

# Check if each area passed in all 90-day windows.
out <- full_results %>%
  group_by(Year, WaterBodyReport, SiteCode, SiteClassification, Ecoli_maxCol, Ecoli_GeoMaxCol)%>%
  summarize(fail=any(Failing))

failing <- out %>% filter(fail == TRUE)
passing <- out %>% filter(fail != TRUE)
unknown <- out %>% filter(is.na(fail))

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

# other things to do:
# * Class B and C should only compute 90-day windows
#    * between April 15th and October 31st
#   * guess: entirety of windows are fully in these months

# To discuss: The earliest and latest samples affect the pass/fail
#  Since many of the 90-day windows only have these samples in them

#FOR KUNAL: CAN WE SEPARATE OUT THE PASSING INSTANTANEOUS AND GEO MEAN?
# If something is failing, you want to know
# Which of INST or GM is failing
# and which 90-day intervals is it failing
# For each site, for each 90-day window:
# * tell me how mnay samples collected
# * tell me if Geo mean pass
# * tell me if 10% rule pass

