library(data.table)
library(survival)

returnUnixDateTime<-function(date) {
  returnVal<-as.numeric(as.POSIXct(date, format="%d/%m/%Y", tz="GMT"))
  return(returnVal)
}

numberAdmissions <- function(dateplustime1, intervalThresholdDays) {
  # x <- renal_cbgMergeDT[CHI == 101296274]
  # dateplustime1 <- x$dateplustime1
  # intervalThresholdDays = 5
  
  dateplustime1 <- dateplustime1[order(dateplustime1)]
  diff_dateplustime1 <- diff(dateplustime1)
  diff_dateplustime1_days <- diff_dateplustime1 / (60*60*24)
  diff_dateplustime1_days <- c(0, diff_dateplustime1_days)
  flagNewAdmission <- ifelse(diff_dateplustime1_days > intervalThresholdDays, 1, 0)
  admissionNumber <- cumsum(flagNewAdmission)
  
  return(admissionNumber)
}

# Importing the dataset
renalDataset = read.csv("~/R/_workingDirectory/renal_cbg/GlycaemiaBaselineData.csv")
cbgDataset = read.csv('~/R/GlCoSy/source/CHIsetCombined_CORE_allAdult_09-16.csv')

renalDatasetDT <- data.table(renalDataset)
  renalDatasetDT <- unique(renalDatasetDT)
cbgDataDT <- data.table(cbgDataset)

renal_cbgMergeDT <- merge(renalDatasetDT, cbgDataDT, by.x = 'CHI', by.y = 'ID')
renal_cbgMergeDT$unix_dateStartingRRT <- returnUnixDateTime(renal_cbgMergeDT$Date.started.RRT)

# order by ID, then by dateplustime1
renal_cbgMergeDT <- renal_cbgMergeDT[order(renal_cbgMergeDT$CHI, renal_cbgMergeDT$dateplustime1), ]
renal_cbgMergeDT[, c('admissionNumber') := numberAdmissions(dateplustime1, 5) , by=.(CHI)]
renal_cbgMergeDT[, c('admissionDuration') := (max(dateplustime1) - min(dateplustime1)) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT$admissionDurationDays <- renal_cbgMergeDT$admissionDuration / (60*60*24)
renal_cbgMergeDT[, c('n_cbg_duringAdmission') := .N , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('cbg_in_sequence_duringAdmission') := seq(1, .N, 1) , by=.(CHI, admissionNumber)]

renal_cbgMergeDT[, c('min_glucose_during_admission') := min(yyyy) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('median_glucose_during_admission') := quantile(yyyy)[3] , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('IQR_glucose_during_admission') := quantile(yyyy)[4] - quantile(yyyy)[2] , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('mean_glucose_during_admission') := mean(yyyy) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('sd_glucose_during_admission') := sd(yyyy) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT$cv_glucose_during_admission <- renal_cbgMergeDT$sd_glucose_during_admission / renal_cbgMergeDT$mean_glucose_during_admission

singleRowPerAdmission <- renal_cbgMergeDT[cbg_in_sequence_duringAdmission == 1]



