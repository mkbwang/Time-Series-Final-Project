# data cleaning of the eagle data
library(lubridate)
library(dplyr)

setwd("~/UM/Biostatistics/Courses/STATS531/Time-Series-Final-Project")

data = read.table("Verreauxs.accel.txt", sep="\t", header=TRUE)

data$date_time = as_datetime(data$date_time)

data$timepoint = as.POSIXct(data$date_time, format="%H:%M:%S")
data$date_time = as.Date(data$date_time)
unique(data$date_time)

data416 = data %>% filter(date_time=="2013-04-16")
data417 = data %>% filter(date_time=="2013-04-17")
data418 = data %>% filter(date_time=="2013-04-18")
data419 = data %>% filter(date_time=="2013-04-19")
data420 = data %>% filter(date_time=="2013-04-20")

# this one has the most number of datapoints
data421 = data %>% filter(date_time=="2013-04-21")
data422 = data %>% filter(date_time=="2013-04-22")
data423 = data %>% filter(date_time=="2013-04-23")
data424 = data %>% filter(date_time=="2013-04-24")

write.csv(data421, file="eagle_421.csv")
