library(dplyr)
library(lubridate)
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


# Function to see if a 90-day window is within April 15 to October 31st
within_b_or_c_window <- function(start) {
  year_of_window <- year(start)
  min_date <- paste(year_of_window, 04, 15, sep="-") %>% ymd() %>% as.Date()
  # Note: because our 90-day windows START on start and END "the day prior to" end,
  # Our 90-day window end date can extend to Nov 1
  max_date <- paste(year_of_window, 11, 1, sep="-") %>% ymd() %>% as.Date()
  return((start >= min_date) & ((start + 90) <= max_date))
}


# Function to check ecoli data against thresholds in a 90-day window
analyze_90_day_window <- function(ecoli_df, start) {
  # Get data within 90-day window
  ecoli_90 <- ecoli_df %>%
    filter(Date >= start) %>%
    filter(Date < (start + 90))

  within_bc_window <- within_b_or_c_window(start)
  
  # Summarize and check thresholds
  ecoli_summary <- ecoli_90 %>%
    group_by(Year, WaterBodyReport, SiteCode, SiteClassification, Ecoli_maxCol, Ecoli_GeoMaxCol) %>%
    summarise(
      NumSamples = n(),
      # Class B/C should only use windows within 4/15 - 10/31
      IsClassBOrC = (unique(SiteClassification) == "B") | (unique(SiteClassification) == "C"),
      WithinBOrCWindow = within_bc_window,
      Is_BC_And_Outside_Window = IsClassBOrC & !WithinBOrCWindow,
      # Ten percent rule
      TenPercentOfSamples = 0.1*n(),
      NumSamplesAboveTenPercentThreshold = sum(c(ResultValue > unique(Ecoli_maxCol))),
      FailingTenPercentRule = NumSamplesAboveTenPercentThreshold > TenPercentOfSamples,
      # Geo mean rule
      geometric_mean = exp(mean(log(ResultValue))),
      EnoughSamplesForGeoMean = (n() >= 6), # Only do Geo Mean rule if at least 6 samples
      FailingGeoMeanRule = geometric_mean > unique(Ecoli_GeoMaxCol),
      # Final results
      # if it's BC, and outside the window, PASS it
      # if it's failing the ten percent rule, FAIL it
      # if it's enough samples for geo mean, and failing the geo mean rule, FAIL it
      Failing = !Is_BC_And_Outside_Window & 
        (FailingTenPercentRule | (EnoughSamplesForGeoMean & FailingGeoMeanRule)),
    )

  return(ecoli_summary)
}

# Function to loop through all 90-day windows
get_ecoli_results <- function(ecoli_df) {
  # Earliest window starts 90 days before the first measurement.
  # Latest window starts date of the last measurement.
  earliest_window_start <- min(ecoli_df$Date) - 90
  latest_window_start <- max(ecoli_df$Date) 
  # If you change the above, it makes a big difference in the analysis!
  # One way to change it: 
  # * earliest_window_start <- min(ecoli_df$Date)
  # * latest_window_start <- max(ecoli_df$Date) - 90
  # To discuss: The earliest and latest samples affect the pass/fail
  #  since many of the 90-day windows only have these samples in them
  
  
  ecoli_full_results <- c() 
  start <- earliest_window_start
  while (start <= latest_window_start) {
    # Decide if rules pass for this 90-day window, and put it in full results.
    ecoli_summary <- analyze_90_day_window(ecoli_df, start)
    
    if (nrow(ecoli_summary) > 0) {
      # the below line doesn't work if there are no rows in ecoli_summary
      ecoli_summary$window_start <- start 
      ecoli_full_results <- rbind(ecoli_full_results,  ecoli_summary)
    }
    
    # Move the 90-day window forward by one day.              
    start <- start + 1
  }
  
  return(ecoli_full_results)
}

# Process the data
full_results <- get_ecoli_results(ecoli_all)

# Check if each site failed in any of the 90-day windows.
sites_all <- full_results %>%
  group_by(Year, WaterBodyReport, SiteCode, SiteClassification, Ecoli_maxCol, Ecoli_GeoMaxCol) %>%
  summarize(fail=any(Failing))

sites_failing <- sites_all %>% filter(fail == TRUE)
sites_passing <- sites_all %>% filter(fail != TRUE)
sites_unknown <- sites_all %>% filter(is.na(fail))


# Plot and inspect failing sites
library(ggplot2)
failing_with_data <- merge(ecoli_all, sites_failing)
failing_with_data$title <- paste("site=", failing_with_data$SiteName, 
                                 " inst=", failing_with_data$Ecoli_maxCol, 
                                 " geo=", failing_with_data$Ecoli_GeoMaxCol,
                                 sep="")

# after the image loads, press "Zoom" to see the whole thing
ggplot(failing_with_data ,
       mapping = aes(x=SampleDatetime,y=ResultValue, color=SiteCode)) +
  geom_point() +
  facet_wrap(facets = vars(title)) +
 theme(legend.position = "none")


## For fun  - Sample analysis
# Most of the time, there aren't enough samples in 90 day window for a geo mean):
# Helpful for feedback -- what do we have to change to collect enough samples?
num_sample_analysis <- full_results %>% 
  group_by(SiteCode, EnoughSamplesForGeoMean) %>% 
  summarise(n()) %>%
  arrange(SiteCode)

## Future ideas:
# * Pull data from portal, and then immediately run this analysis
# * Add standards for each site classification 
    # (inspect other states' protocols and see how similar they are to maine's)