#Comment out the command below after first run
#install.packages('readr')

source("Functions.R")

#The following code will allow you to calculate the probability of having a HIV-1 infection
#at the time of diagnosis of a subsequent condition.

#The file "Pat_conditions.csv" list for a given patient their INFCareHIV ID, as
#provided in meta files from April 2023. The dates for each condition diagnosed,
#as well as the conditions. An example file has been provided.

#The code below allows you to generate these probabilities for all patients within the 
#given data set or for single patient (Single_pat function).
#The for function (all_pat) is already set to automatically generate for all patients and the 
#only inputs that you may wish to change are "plots" and "csv".
#By default these are set to false for "plots" and "true" for csv, so no plots will be generated 
#only output csv files. IF "csv" is set to FALSE a table will be printed of the results per patient

#Changing "plots=TRUE" will generate a probability distribution plot for each
#condition. The colour green will indicates the probability of being diagnosed with the
#condition and infected with HIV-1. The red colour indicates the probability
#without HIV-1 infection.

#Changing "csv=TRUE" will generate an output .csv file with the probabilities for
#each condition along with the date of each condition. A file will be generated 
#for all patients.

#import .csv with patient ID, Date and conditions.
Pat_conditions <- read_csv("Pat_conditions.csv")


#Use you want to generate plots and data for every patient.
all_pat(Pat_conditions,plots=FALSE,csv=TRUE)

#Use if you would like to examine a single patient.
Single_pat(ID=1,Pat_conditions,plots=FALSE,csv=FALSE)
