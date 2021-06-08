## Function for processing zinc data
## Inputs:
##.  theData: a dataframe of biomarker measurement values with the following columns:
##.     1. groupId: The age-gender grouping of the sample e.g. WRA/PSC
##.     2. regionName: Name of the subregion containing the household for aggregtion purposes
##.     3. ageInMonths: The age of the individual sampled
##.     4. isPregnant:  A true/false flag indicating the pregnancy status of the individual
##.     5. <biomarkerName>:  The name of the biomarker under consideration e.g. zinc
##.
##.  groupId: The age-gender grouping of the sample e.g. WRA/PSC
##.
##.  magic: A number
#  BiomarkerData <- function(theData, groupId, dsgn) {



#### LOCAL TESTING ONLY ####
# this will be replaced by API code
MyGp<-"WRA"  # change this to change the demographic group data used
k<-6 # change this to change the MN brought in.  ##### ALTERED

MyList<-c(1:4,8:15,k)
MyMN<-length(MyList)  # this only works for one MN NB.
MyAgg<-'regionName' # this is the aggregate data field to be used

# Thresholds
ThresholdUpper<-"No value" #Can have text, e.g. "No value" and script will function correctly
ThresholdLower<-"No value" #
PhysLim<-3000 # fixed figure until a csv to import

#Import some data and restrict to the demo gp and MN of interest
DataUse<-read.csv("C:/Users/FannySandalinas/Documents/zinc.csv") # this is the participant and biomarker data

# creat a random number of serum and RBC create #
library(tidyverse)
set.seed(2)

DataUse <- DataUse %>% mutate(serum_folate = runif(3000, min = 2.01, max = 78.10))

DataUse <- DataUse %>% mutate(RBC_folate = runif(3000, min = 59.69, max = 2098.24))

DataUse<-DataUse[(MyList)]   ## here subset the MN to use in this script
DataUse<-DataUse[which(DataUse$groupId==MyGp),] # select demo gp data only

#### Flags to accomodate differences between surveys
Flag_SurvWeightRun<-1 # where 1 = 'run' as not adjusted in the supplied data and 0 = 'do not run' as already adjusted
Flag_SurvWeightSUpplied<-1 # where not supplied set to 0, supplied = 1. Run if above = 1.
Flag_HaemAltAdjust<-0  # Adjustment not run (0) where already adjusted in the supplied data, run if = 1.
Flag_SmokeAdjust<-0 # Adjustment not run (0) where already adjusted in the supplied data, or smoking not recorded so cannot adjust

#Call the stats function
#Stats<-BiomarkerData(DataUse, 'WRA', DHSdesign) {

#### MAPS tool biomarker script 1a ####

#Packages and macros needed by R to run the script below
#macroBMI not currently used so commented out - AKB
#source('macroBMI.R')  #https://rdrr.io/github/KHP-IDEO/rIDEO/man/grp_bmi.html
library(survey) # required to apply svydesign() function
library(psych) # Required for describe() function

#### Eligible data selection ####
DataUse[,MyMN] = as.numeric(as.character(DataUse[,MyMN])) #Make sure the biomarker column values are numeric

DataUse <- DataUse[!(is.na(DataUse[,MyMN])),] #omit row with NA in the col of interest

DataUse<-DataUse[which(DataUse[,MyMN]<=PhysLim),]   #Excluding physiological implausible concentrations for the specific MN

# assign age categories, then exclude the '0' rows
DataUse$AgeCat<-ifelse(DataUse$groupId=="WRA",
                       ifelse(DataUse$ageInMonths>=(15*12) & DataUse$ageInMonths<(20*12),1,ifelse(DataUse$ageInMonths>=(20*12) & DataUse$ageInMonths<(50*12),2,0)),
                       ifelse(DataUse$groupId=="PSC",
                              ifelse(DataUse$ageInMonths>=(0.5*12) & DataUse$ageInMonths<(2*12),1,ifelse(DataUse$ageInMonths>=(2*12) & DataUse$ageInMonths <(5*12),2,0)),
                              ifelse(DataUse$groupId=="SAC",
                                     ifelse(DataUse$ageInMonths>=(5*12) & DataUse$ageInMonths<(11*12),1,ifelse(DataUse$ageInMonths>=(12*12) & DataUse$ageInMonths <(15*12),2,0)),
                                     ifelse(DataUse$groupId=="MEN",
                                            ifelse(DataUse$ageInMonths>=(15*12) & DataUse$ageInMonths<(30*12),1,ifelse(DataUse$ageInMonths>=(30*12) & DataUse$ageInMonths <(54*12),2,0)),
                                            0))))

DataUse<-DataUse[which(DataUse$AgeCat>0),]

DataUse$DemoGpCat<-paste(DataUse$groupId,".",DataUse$AgeCat)

DataUse<-DataUse[which(DataUse$isPregnant==0),]  #exclude pregnant

#### Adjustments ####
#adjust Haem for altitude (all demo gps)
if(Flag_HaemAltAdjust == 1) {
  #code to adjust for altitude, where HEMOGLOBIN_UNADJ is the unadjusted hemoglobin value given by lab. HEMO is the altitude-adjusted value
  DataUse$HEMO_factor= (-0.032*(DataUse$Altitude*0.0032808) + 0.022*(DataUse$Altitude*0.0032808)^2)
  DataUse$HEMO_ADJUSTED<-((DataUse$HEMOGLOBIN_LAB)-(+DataUse$HEMO_factor))
  DataUse$HEMO<-ifelse(DataUse$Altitude>1000,DataUse$HEMO_ADJUSTED,DataUse$HEMOGLOBIN_LAB)
}

# Adjust Haem for smoking
if(Flag_SmokeAdjust == 1) {
  
  ifelse (Number_cigarettes_known=='YES', ifelse (DataUse$NumberCigarettes<10, DataUse$HEMOS<-DataUse$HEMO, ifelse(DataUse$NumberCigarettes>=10 & DataUse$NumberCigarettes<20,
                                                                                                                   DataUse$HEMOS<-(DataUse$HEMO - 0.3), ifelse(DataUse$NumberCigarettes>=20 & DataUse$NumberCigarettes<40, DataUse$HEMOS<-(DataUse$HEMO-0.5), 
                                                                                                                                                               ifelse(DataUse$NumberCigarettes>=40, DataUse$HEMOS<-DataUse$HEMO-0.7))),
                                                  ifelse (Number_cigarettes_known=='NO', DataUse$HEMOS<-DataUse$HEMO-0.3)))
}

#### Survey Weights ####
# Adjust by survey weight where not already done AND weights supplied, *else* run without use weights to get common output format
#  make sure cluster, strat and weight always have the same name and structure (weight to divide by 1000000)
if(Flag_SurvWeightRun == 1 & Flag_SurvWeightSUpplied ==1){
  
  DHSdesign<-svydesign(id=DataUse$cluster, strata=DataUse$strata, weights = DataUse$weights/1000000, data = DataUse, nest = TRUE)
  options("survey.lonely.psu"='adjust')
} else {
  DHSdesign<-svydesign(ids = ~1, strata=NULL , weights = NULL , data = DataUse)
}


#### END MAPS tool biomarker Scrript 1a ####

#### output data ####

# Set up whether a record is < or > thresholds.
# NB works for those where low concs indicate deficiency, not vice versa
#DataUse$ThresholdComp<-ifelse(DataUse[,MyMN]<thresholdLower,"low",
#ifelse(DataUse[,MyMN]>thresholdUpper,"high","not low"))
#ThresholdCount<-as.data.frame(table(DataUse$ThresholdComp)  )

#for zinc we will need to define the threshold according to fasting and time of the day

#if(MyMN=='DataUse$zinc'){  #does that work? we could also run that whatever the micronutrient is

DataUse$BLOOD_MORNING<-ifelse(DataUse$Sampling_time<=12,1,0)

DataUse$ZINC_DEF<-ifelse(DataUse$zinc<70 & 
                           DataUse$FASTING_STATUS==1 & DataUse$BLOOD_MORNING==1,1,
                         ifelse(DataUse$zinc<66 &DataUse$BLOOD_MORNING==1 & DataUse$FASTING_STATUS==0,1,
                                ifelse(DataUse$zinc<59 & DataUse$BLOOD_MORNING==0,1,0)))
#}

# define folate deficiency  

DataUse$ser_fol_def <- ifelse(DataUse$serum_folate<6.8,1,0)# using macrocytic anaemia as a haematological indicator

DataUse$ser_fol_def1 <- ifelse(DataUse$serum_folate<10,1,0)# using homocysteine concentrations as metabolic indicator

DataUse$rbc_fol_def <- ifelse(DataUse$RBC_folate<226.5,1,0)# using macrocytic anaemia as a haematological indicator

DataUse$rbc_fol_def1 <- ifelse(DataUse$RBC_folate<340,1,0)# using homocysteine concentrations as metabolic indicator

DataUse$rbc_fol_def2 <- ifelse(DataUse$RBC_folate<906,1,0)#  for preventing neural tube defect-affected pregnancies in women of reproductive age at the population levela


# Adjust by survey weight where not already done AND weights supplied
if(Flag_SurvWeightRun == 1 & Flag_SurvWeightSUpplied ==1){
  
  DHSdesign<-svydesign(id=DataUse$cluster, strata=DataUse$strata, weights = DataUse$weights/1000000, data = DataUse, nest = TRUE)
  options("survey.lonely.psu"='adjust')
}





# write boxplot output by admin area # here we can replace zinc by any micronutrient

svyboxplot(serum_folate~DataUse$regionName, plot=FALSE, DHSdesign)


####summary stats#### can we delete the unweighted stats?
#  all dataset
Stat1<-data.frame(describe(DataUse[,MyMN]))
Stat2<-data.frame(t(quantile(DataUse[,MyMN],c(.25, .50, .75),na.rm=TRUE)))
StatOutputAllData<-cbind(Stat1,Stat2) #not all content needed front end

# disagg data summary
Stat3<-describeBy(DataUse[,MyMN],DataUse[,MyAgg],mat = TRUE,digits = 2)
Stat3<-as.data.frame(Stat3) # most summary stats

##weighted summary stats## here we can replace zinc by any micronutrient

svyquantile(~zinc, design = DHSdesign, quantiles = c(0.25,0.5,0.75))#Quantile for the total sample#

boxstat1<- as.data.frame(svyby(~zinc, ~regionName, DHSdesign, svyquantile, quantiles=c(0.0,0.25,0.5,0.75,1), keep.var=F))#Quantile by administrative region#   

boxstat2<- as.data.frame(svyby(~zinc, ~AgeCat, DHSdesign, svyquantile, quantiles=c(0.0,0.25,0.5,0.75,1), keep.var=F))#Quantile by different demographic groups#   

DataUse$weights <- DataUse$weights/1000000

library(srvyr)

strat_design_srvyr <- DataUse %>% as_survey_design(id = cluster, strata = strata, weights = weights, nest = TRUE)

strat_design_srvyr %>%
  group_by(regionName) %>%
  summarize(mean = survey_mean(zinc))

stat <- strat_design_srvyr %>%
  group_by(regionName) %>%
  summarise(
    mean = survey_mean(zinc),
    sd = survey_sd(zinc),
    Q = survey_quantile(zinc, c(0.25, 0.5, 0.75))
  ) %>%
  mutate(IQR = Q_q75 - Q_q25) %>%
  mutate(
    out_upp = Q_q75 + 1.5 * IQR,
    out_low = Q_q25 - 1.5 * IQR
  ) 

a <- left_join(stat, DataUse, by = "regionName")

a %>% select(zinc,out_upp,out_low,regionName) %>% filter(zinc< out_low | zinc>out_upp)

# Stat4<-as.data.frame(B1$stats) #take quartiles from boxplot B1 outputs
#Stat4<-Stat4[2:4,]
# rownames(Stat4)<-c("25th","50th","75th")
#Stat5<-as.data.frame(t(Stat4))
# Stat5$group1<-B1$names

#StatsOutput_agg<-merge(Stat3,Stat5,by = "group1") #not all content needed front end

#deficiency prevalence


MyaggZinc<-svyby(~ZINC_DEF, ~DataUse$regionName, DHSdesign, svyciprop, vartype="ci")

Myaggfolate<-svyby(~ser_fol_def, ~DataUse$regionName, DHSdesign, svyciprop, vartype="ci")
Myaggfolate2<-svyby(~ser_fol_def1, ~DataUse$regionName, DHSdesign, svyciprop, vartype="ci")
Myaggfolate3<-svyby(~rbc_fol_def, ~DataUse$regionName, DHSdesign, svyciprop, vartype="ci")
Myaggfolate4<-svyby(~rbc_fol_def1, ~DataUse$regionName, DHSdesign, svyciprop, vartype="ci")
Myaggfolate5<-svyby(~rbc_fol_def2, ~DataUse$regionName, DHSdesign, svyciprop, vartype="ci")

Myaggfolate1<-cbind(Myaggfolate,Myaggfolate2,Myaggfolate3, Myaggfolate4, Myaggfolate5)
Myaggfolate10<-Myaggfolate1[-c(1,5,9,13,17)]
Myaggfolate100<-Myaggfolate10*100





#### end ####

# print(StatsOutput_agg)

# Return aggregated stats table
print(StatOutputAllData)
print(TableZinc)
}


