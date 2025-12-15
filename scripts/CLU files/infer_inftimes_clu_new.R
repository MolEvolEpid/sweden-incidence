#script to run infection age inference on cluster

library(biophybreak)

source("run.mbm.scratch.R")
source("mbm.predict.scratch.R")

#find individual from command line
i <- commandArgs(trailingOnly = TRUE)
i <- as.integer(i)
i=1
#load dataframe
load("df.incidence.informed.1B.Rdata")

#subset patient of interest
for (i in 1:100){
    df <- df.informed[i,]
    
    if (length(df$CD4_qc[[1]])<4){
      df$CD4_qc[[1]]<-df$CD4_qc[[1]][1:length(df$CD4_qc[[1]])]
      df$CD4_dates[[1]]<-df$CD4_dates[[1]][1:length(df$CD4_qc[[1]])]
    }else{
      df$CD4_qc[[1]]<-tail(df$CD4_qc[[1]],4)
      df$CD4_dates[[1]]<-tail(df$CD4_dates[[1]],4)
    }
    temp_cd4<-df$CD4_qc[[1]]
    temp_cd4_dates<-df$CD4_dates[[1]]

    for (j in 1:length(df$CD4_qc[[1]])){
      df$CD4_qc[[1]]<-temp_cd4[j]
      print(df$CD4_qc[[1]])
      df$CD4_dates[[1]]<-temp_cd4_dates[j]
      print(df$CD4_dates[[1]])


infection.ages <- run.mbm(df, 
                          n.adapt = 1e4,
                          n.burn = 1e5, 
                          n.iter = 1e6,
                          prior.type = 1, 
                          overall.seed = 20230424+i)

saveRDS(infection.ages, file = paste0("~/Documents/New_incidence_work/Training_Data/single_CD4/",j,"M/infection.ages_", i, ".rds"))
  }
}
