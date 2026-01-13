#script to run infection age inference on cluster

library(biophybreak)

source("run.mbm.scratch.R")
source("mbm.predict.scratch.R")

#find individual from command line
i <- commandArgs(trailingOnly = TRUE)
i <- as.integer(i)

#load dataframe
load("df.incidence.informed.1B.Rdata")
df <- df.informed[i,]
#subset patient of interest



infection.ages <- run.mbm(df, 
                          n.adapt = 1e4,
                          n.burn = 1e5, 
                          n.iter = 1e6,
                          prior.type = 1, 
                          overall.seed = 20230424+i)

saveRDS(infection.ages, file = paste0("~/infection.ages_", i, ".rds"))

