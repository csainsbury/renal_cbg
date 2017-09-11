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


simpleSurvivalPlot<-function(inputFrame,endDateUnix,ylimMin) {
  # inputFrame <- singleRowPerAdmission_lessThan1DayAdmission_firstAdmission
  # endDateUnix <- max(inputFrame$unix_deathDate)
  # sampleDateUnix <- 
  
  SurvivalData<-inputFrame
  
  DaySeconds<-(60*60*24)
  shortCensorPeriodStartDay  <- DaySeconds
  shortCensorPeriodEndDay    <- DaySeconds*10000
  
  lastDOD<-endDateUnix
  SurvivalData$dateOfDischarge<-inputFrame$dateplustime1 + (intervalToDetermineAdmissionDays * (60*60*24))
  SurvivalData$timeToDeath<-ifelse(SurvivalData$isDead==1,(SurvivalData$unix_deathDate-SurvivalData$dateOfDischarge),0)
  #		SurvivalData$timeToDeath<-SurvivalData$timeToDeath/DaySeconds
  SurvivalData$timeToDeathInterval<-ifelse(SurvivalData$isDead==0,(lastDOD-SurvivalData$dateOfDischarge),SurvivalData$timeToDeath)
  SurvivalData$timeToDeathInterval[is.na(SurvivalData$timeToDeathInterval)]<-0; SurvivalData<-subset(SurvivalData,timeToDeathInterval>0)
  # SurvivalData$timeToDeathInterval<-SurvivalData$timeToDeathInterval/(60*60*24*365.25)
  
  SurvivalData$shortDeathEvent <- SurvivalData$isDead
  # SurvivalData$shortDeathEvent <- ifelse(SurvivalData$isDead==1 & SurvivalData$timeToDeath>=(shortCensorPeriodStartDay) & SurvivalData$timeToDeath<(shortCensorPeriodEndDay),1,0)	
  
  #  SurvivalData$sexDigit<-ifelse(nchar(SurvivalData$charID==9),as.numeric(substr(SurvivalData$charID,8,8)),as.numeric(substr(SurvivalData$charID,9,9)))
  # SurvivalData$sexNumber<-ifelse(SurvivalData$sexDigit%%2==0,1,0)
  #  SurvivalData$sex<-factor(1*(SurvivalData$sexNumber <1),levels=0:1,labels=c("F","M"))
  
  
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (min_glucose_during_admission < 3), data = SurvivalData)
  shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds," days\n n= ",nrow(SurvivalData),", threshold: ",quantile(SurvivalData$hba1cIQRinRange)[3],sep="")
  plot(mfitAge50,mark.time=T,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=5,ylim=c(ylimMin,1))
  
  # mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age_atSampleTime+medianHbA1cInRange+nValsPerIDinRange+(hba1cIQRinRange>=quantile(SurvivalData$hba1cIQRinRange)[3]), data = SurvivalData)
  
  mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age_atSampleTime+diabetesDurationYears+medianHbA1cInRange+nValsPerIDinRange+(hba1cIQRinRange>=quantile(SurvivalData$hba1cIQRinRange)[3]), data = SurvivalData)
  pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
  legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
  summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
  legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)
  
  print(mfitAge50.coxph)
  
}

# Importing the dataset
renalDataset = read.csv("~/R/_workingDirectory/renal_cbg/GlycaemiaBaselineData.csv")
cbgDataset = read.csv('~/R/GlCoSy/source/CHIsetCombined_CORE_allAdult_09-16.csv')

# set admission defining gap
intervalToDetermineAdmissionDays <- 1.2

renalDatasetDT <- data.table(renalDataset)
  renalDatasetDT <- unique(renalDatasetDT)
cbgDataDT <- data.table(cbgDataset)

renal_cbgMergeDT <- merge(renalDatasetDT, cbgDataDT, by.x = 'CHI', by.y = 'ID')
renal_cbgMergeDT$unix_dateStartingRRT <- returnUnixDateTime(renal_cbgMergeDT$Date.started.RRT)
renal_cbgMergeDT$unix_deathDate <- returnUnixDateTime(renal_cbgMergeDT$date.death)
  renal_cbgMergeDT$unix_deathDate[is.na(renal_cbgMergeDT$unix_deathDate)] <- 0
  renal_cbgMergeDT$isDead <- ifelse(renal_cbgMergeDT$unix_deathDate > 0, 1, 0)

# order by ID, then by dateplustime1
renal_cbgMergeDT <- renal_cbgMergeDT[order(renal_cbgMergeDT$CHI, renal_cbgMergeDT$dateplustime1), ]
renal_cbgMergeDT[, c('admissionNumber') := numberAdmissions(dateplustime1, intervalToDetermineAdmissionDays) , by=.(CHI)]
renal_cbgMergeDT[, c('admissionDuration') := (max(dateplustime1) - min(dateplustime1)) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT$admissionDurationDays <- renal_cbgMergeDT$admissionDuration / (60*60*24)
renal_cbgMergeDT[, c('n_cbg_duringAdmission') := .N , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('cbg_in_sequence_duringAdmission') := seq(1, .N, 1) , by=.(CHI, admissionNumber)]

renal_cbgMergeDT[, c('min_glucose_during_admission') := min(yyyy) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('max_glucose_during_admission') := max(yyyy) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('median_glucose_during_admission') := quantile(yyyy)[3] , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('IQR_glucose_during_admission') := quantile(yyyy)[4] - quantile(yyyy)[2] , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('mean_glucose_during_admission') := mean(yyyy) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT[, c('sd_glucose_during_admission') := sd(yyyy) , by=.(CHI, admissionNumber)]
renal_cbgMergeDT$cv_glucose_during_admission <- renal_cbgMergeDT$sd_glucose_during_admission / renal_cbgMergeDT$mean_glucose_during_admission

# generate one row per admission dataset
singleRowPerAdmission <- renal_cbgMergeDT[cbg_in_sequence_duringAdmission == 1]

# remove admissions before dialysis started
singleRowPerAdmission <- singleRowPerAdmission[dateplustime1 > unix_dateStartingRRT]

# Q ? how long does a dialysis admission last
singleRowPerAdmission_dialysisAdmission <- singleRowPerAdmission[admissionDurationDays < 1]

########## survival analysis
# flag first admission
# singleRowPerAdmission_lessThan1DayAdmission[, c('flag_first_admission') := ifelse(admissionNumber == min(admissionNumber), 1, 0) , by=.(CHI)]
# 
# singleRowPerAdmission_lessThan1DayAdmission_firstAdmission <- singleRowPerAdmission_lessThan1DayAdmission[flag_first_admission == 1]

# flag first n admissions, for those with at least n admissions
n = 2
singleRowPerAdmission_dialysisAdmission[, c('flag_first_n_admissions') := ifelse(admissionNumber <= n, 1, 0) , by=.(CHI)]
singleRowPerAdmission_dialysisAdmission[, c('indicate_at_least_n_admission') := ifelse(max(admissionNumber) >= n, 1, 0) , by=.(CHI)]
singleRowPerAdmission_dialysisAdmission[, c('total_N_admissions') := max(admissionNumber), by = .(CHI)]

# histogram of numbers of dialysis admissions
hist(singleRowPerAdmission_dialysisAdmission[admissionNumber == 1]$total_N_admissions, breaks = seq(1, 1000, 1), xlim = c(0, 10))

# subset of first n admissions for those with at least n admissions
singleRowPerAdmission_dialysisAdmission_n_admissionData <- singleRowPerAdmission_dialysisAdmission[indicate_at_least_n_admission == 1]

uniqueN(singleRowPerAdmission_dialysisAdmission_n_admissionData$CHI)


# flag if any hypos in first 10 admissions
singleRowPerAdmission_dialysisAdmission_n_admissionData[, c('flag_if_hypo_in_n_admissions') := ifelse(min(min_glucose_during_admission) < 4, 1, 0) , by=.(CHI)]






