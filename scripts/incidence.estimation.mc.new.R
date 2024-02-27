#script for monte carlo estimation of incidence

source("scripts/incidence.mc.functions.R")

#date to use to tag file input/output
date_tag <- "20231024"

#which prior was used?
prior <- "gamma"

#read in inferred distributions
infection.ages <- readRDS(paste0("infection.ages.", prior, ".all.new.", date_tag, ".rds"))

#read in typical distributions
typical.dists <- readRDS(file = paste0("typical.distributions.risk.ee.", prior, ".new.", date_tag, ".rds"))

#convert from decimal to year of diagnosis
infection.ages$diagnosis_year <- decimal2year(as.numeric(infection.ages$diagnosis_date))

#convert from decimal to year of arrival
infection.ages$arrival_year <- decimal2year(as.numeric(infection.ages$arrival_date))

#convert from decimal to year of removal
infection.ages$removal_year <- decimal2year(as.numeric(infection.ages$removal_date))

#study end date: early 2023
tau1 <- 1970
tau3 <- max(infection.ages$diagnosis_date)+1/365.25 #one day after latest diagnosis time

#find diagnosis probabilities based on time from infection to study end
diagnosis.prob.fns <- diagnosis.prob.fn(typical.dists, tau1 = tau1, tau3 = tau3)

#find typical arrival time distributions
typical.arrival.dists <- typical.arrival.dist(infection.ages, risk_groups = c("HET", "MSM", "IDU", "unknown/other"))

#find icdfs of infection ages
inf_icdfs <- lapply(infection.ages$infection.age, FUN = find.icdf)
saveRDS(inf_icdfs, file = paste0("inf_icdfs.", prior, ".", date_tag, ".rds"))

#find monte carlo estimates for incidence
incidence <- estimate.incidence.mc(infection.ages = infection.ages, 
                                   inf_icdfs = inf_icdfs, 
                                   diagnosis.prob.fns = diagnosis.prob.fns, 
                                   typical.arrival.dists = typical.arrival.dists, 
                                   tau1 = tau1,
                                   tau3 = tau3,
                                   n_samples = 10000)
saveRDS(incidence, file = paste0("incidence.mc.", prior, ".", date_tag, ".rds"))

####calculate additional stats
###proportions of incidence
##percentages endo and exo
#endo
incidence$endo_combined_risk$incidence_prop <- incidence$endo_combined_risk$incidence/incidence$total$incidence
incidence$endo_combined_risk$incidence_prop_means <- rowMeans(incidence$endo_combined_risk$incidence_prop)
incidence$endo_combined_risk$incidence_prop_95 <- 
  apply(incidence$endo_combined_risk$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo
incidence$exo_combined_risk$incidence_prop <- incidence$exo_combined_risk$incidence/incidence$total$incidence
incidence$exo_combined_risk$incidence_prop_means <- rowMeans(incidence$exo_combined_risk$incidence_prop)
incidence$exo_combined_risk$incidence_prop_95 <- 
  apply(incidence$exo_combined_risk$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

##percentages HET, IDU, MSM, and UO
#HET
incidence$combined_ee$HET$incidence_prop <- incidence$combined_ee$HET$incidence/incidence$total$incidence
incidence$combined_ee$HET$incidence_prop_means <- rowMeans(incidence$combined_ee$HET$incidence_prop)
incidence$combined_ee$HET$incidence_prop_95 <- 
  apply(incidence$combined_ee$HET$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#IDU
incidence$combined_ee$IDU$incidence_prop <- incidence$combined_ee$IDU$incidence/incidence$total$incidence
incidence$combined_ee$IDU$incidence_prop_means <- rowMeans(incidence$combined_ee$IDU$incidence_prop)
incidence$combined_ee$IDU$incidence_prop_95 <- 
  apply(incidence$combined_ee$IDU$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#MSM
incidence$combined_ee$MSM$incidence_prop <- incidence$combined_ee$MSM$incidence/incidence$total$incidence
incidence$combined_ee$MSM$incidence_prop_means <- rowMeans(incidence$combined_ee$MSM$incidence_prop)
incidence$combined_ee$MSM$incidence_prop_95 <- 
  apply(incidence$combined_ee$MSM$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#UO
incidence$combined_ee$`unknown/other`$incidence_prop <- incidence$combined_ee$`unknown/other`$incidence/incidence$total$incidence
incidence$combined_ee$`unknown/other`$incidence_prop_means <- rowMeans(incidence$combined_ee$`unknown/other`$incidence_prop)
incidence$combined_ee$`unknown/other`$incidence_prop_95 <- 
  apply(incidence$combined_ee$`unknown/other`$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

##percentages combinations of endo/exo and risk
#endo HET
incidence$endo_split$HET$incidence_prop <- incidence$endo_split$HET$incidence/incidence$total$incidence
incidence$endo_split$HET$incidence_prop_means <- rowMeans(incidence$endo_split$HET$incidence_prop)
incidence$endo_split$HET$incidence_prop_95 <- 
  apply(incidence$endo_split$HET$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#endo IDU
incidence$endo_split$IDU$incidence_prop <- incidence$endo_split$IDU$incidence/incidence$total$incidence
incidence$endo_split$IDU$incidence_prop_means <- rowMeans(incidence$endo_split$IDU$incidence_prop)
incidence$endo_split$IDU$incidence_prop_95 <- 
  apply(incidence$endo_split$IDU$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#endo MSM
incidence$endo_split$MSM$incidence_prop <- incidence$endo_split$MSM$incidence/incidence$total$incidence
incidence$endo_split$MSM$incidence_prop_means <- rowMeans(incidence$endo_split$MSM$incidence_prop)
incidence$endo_split$MSM$incidence_prop_95 <- 
  apply(incidence$endo_split$MSM$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#endo UO
incidence$endo_split$`unknown/other`$incidence_prop <- incidence$endo_split$`unknown/other`$incidence/incidence$total$incidence
incidence$endo_split$`unknown/other`$incidence_prop_means <- rowMeans(incidence$endo_split$`unknown/other`$incidence_prop)
incidence$endo_split$`unknown/other`$incidence_prop_95 <- 
  apply(incidence$endo_split$`unknown/other`$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo HET
incidence$exo_split$HET$incidence_prop <- incidence$exo_split$HET$incidence/incidence$total$incidence
incidence$exo_split$HET$incidence_prop_means <- rowMeans(incidence$exo_split$HET$incidence_prop)
incidence$exo_split$HET$incidence_prop_95 <- 
  apply(incidence$exo_split$HET$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo IDU
incidence$exo_split$IDU$incidence_prop <- incidence$exo_split$IDU$incidence/incidence$total$incidence
incidence$exo_split$IDU$incidence_prop_means <- rowMeans(incidence$exo_split$IDU$incidence_prop)
incidence$exo_split$IDU$incidence_prop_95 <- 
  apply(incidence$exo_split$IDU$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo MSM
incidence$exo_split$MSM$incidence_prop <- incidence$exo_split$MSM$incidence/incidence$total$incidence
incidence$exo_split$MSM$incidence_prop_means <- rowMeans(incidence$exo_split$MSM$incidence_prop)
incidence$exo_split$MSM$incidence_prop_95 <- 
  apply(incidence$exo_split$MSM$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo UO
incidence$exo_split$`unknown/other`$incidence_prop <- incidence$exo_split$`unknown/other`$incidence/incidence$total$incidence
incidence$exo_split$`unknown/other`$incidence_prop_means <- rowMeans(incidence$exo_split$`unknown/other`$incidence_prop)
incidence$exo_split$`unknown/other`$incidence_prop_95 <- 
  apply(incidence$exo_split$`unknown/other`$incidence_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

###proportions of current infections

##percentages combinations of endo/exo and risk
#endo HET
incidence$endo_split$HET$infections_current_prop <- incidence$endo_split$HET$infections_current/incidence$total$infections_current
incidence$endo_split$HET$infections_current_prop_means <- rowMeans(incidence$endo_split$HET$infections_current_prop)
incidence$endo_split$HET$infections_current_prop_95 <- 
  apply(incidence$endo_split$HET$infections_current_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo HET
incidence$exo_split$HET$infections_current_prop <- incidence$exo_split$HET$infections_current/incidence$total$infections_current
incidence$exo_split$HET$infections_current_prop_means <- rowMeans(incidence$exo_split$HET$infections_current_prop)
incidence$exo_split$HET$infections_current_prop_95 <- 
  apply(incidence$exo_split$HET$infections_current_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo UO
incidence$exo_split$`unknown/other`$infections_current_prop <- incidence$exo_split$`unknown/other`$infections_current/incidence$total$infections_current
incidence$exo_split$`unknown/other`$infections_current_prop_means <- rowMeans(incidence$exo_split$`unknown/other`$infections_current_prop)
incidence$exo_split$`unknown/other`$infections_current_prop_95 <- 
  apply(incidence$exo_split$`unknown/other`$infections_current_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

###proportions of undiagnosed cases

#endo HET
incidence$endo_split$HET$undiagnosed_current_prop <- incidence$endo_split$HET$undiagnosed_current/incidence$total$undiagnosed_current
incidence$endo_split$HET$undiagnosed_current_prop_means <- rowMeans(incidence$endo_split$HET$undiagnosed_current_prop)
incidence$endo_split$HET$undiagnosed_current_prop_95 <- 
  apply(incidence$endo_split$HET$undiagnosed_current_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#endo MSM
incidence$endo_split$MSM$undiagnosed_current_prop <- incidence$endo_split$MSM$undiagnosed_current/incidence$total$undiagnosed_current
incidence$endo_split$MSM$undiagnosed_current_prop_means <- rowMeans(incidence$endo_split$MSM$undiagnosed_current_prop)
incidence$endo_split$MSM$undiagnosed_current_prop_95 <- 
  apply(incidence$endo_split$MSM$undiagnosed_current_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

###miscellaneous

#IDU+UO current infection counts and proportions of current infections
#counts
IDU_UO_infections_current <- incidence$combined_ee$IDU$infections_current + incidence$combined_ee$`unknown/other`$infections_current
IDU_UO_infections_current_means <- rowMeans(IDU_UO_infections_current)
IDU_UO_infections_current_95 <- apply(IDU_UO_infections_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
#proportions
IDU_UO_infections_current_prop <- IDU_UO_infections_current/incidence$total$infections_current
IDU_UO_infections_current_prop_means <- rowMeans(IDU_UO_infections_current_prop)
IDU_UO_infections_current_prop_95 <- apply(IDU_UO_infections_current_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo IDU+UO current infection counts and proportions of current infections
exo_IDU_UO_infections_current <- incidence$exo_split$IDU$infections_current + incidence$exo_split$`unknown/other`$infections_current
exo_IDU_UO_infections_current_means <- rowMeans(exo_IDU_UO_infections_current)
exo_IDU_UO_infections_current_95 <- apply(exo_IDU_UO_infections_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
#proportions
exo_IDU_UO_infection_current_prop <- exo_IDU_UO_infections_current/incidence$total$infections_current
exo_IDU_UO_infections_current_prop_means <- rowMeans(exo_IDU_UO_infection_current_prop)
exo_IDU_UO_infections_current_prop_95 <- apply(exo_IDU_UO_infection_current_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo MSM and UO undiagnosed cases
exo_MSM_UO_undiagnosed_current <- incidence$exo_split$MSM$undiagnosed_current + incidence$exo_split$`unknown/other`$undiagnosed_current
exo_MSM_UO_undiagnosed_current_means <- rowMeans(exo_MSM_UO_undiagnosed_current)
exo_MSM_UO_undiagnosed_current_95 <- apply(exo_MSM_UO_undiagnosed_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#exo HET and UO undiagnosed cases
exo_HET_UO_undiagnosed_current <- incidence$exo_split$HET$undiagnosed_current + incidence$exo_split$`unknown/other`$undiagnosed_current
exo_HET_UO_undiagnosed_current_means <- rowMeans(exo_HET_UO_undiagnosed_current)
exo_HET_UO_undiagnosed_current_95 <- apply(exo_HET_UO_undiagnosed_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

#save with additional stats
saveRDS(incidence, file = paste0("incidence.mc.", prior, ".", date_tag, ".more.stats.rds"))

###calculate stats of infection ages at diagnosis over time
#read in informed infection age distributions
infection.ages.informed <- readRDS(file = paste0("infection.ages.gamma.informed.new.", date_tag, ".rds"))

#convert from decimal to year of diagnosis
infection.ages.informed$diagnosis_year <- decimal2year(as.numeric(infection.ages.informed$diagnosis_date))

#find icdfs of informed infection ages
inf_icdfs_informed <- lapply(infection.ages.informed$infection.age, FUN = find.icdf)
saveRDS(inf_icdfs_informed, file = paste0("inf_icdfs_informed.", prior, ".", date_tag, ".rds"))

#find means and standard deviations 
inf_age_draws <- lapply(inf_icdfs_informed, FUN = function(icdf, n_samples) icdf(runif(n_samples)), n_samples = 10000)
inf_age_means <- sapply(inf_age_draws, FUN = mean)
inf_age_sd <- sapply(inf_age_draws, FUN = sd)

#data frame for mean TI for known endo patients
known_endo_indices <- which(infection.ages.informed$endoexo == "endo")
endo_TI_means <- data.frame(year = infection.ages.informed$diagnosis_year[known_endo_indices],
                            route = infection.ages.informed$risk_group[known_endo_indices],
                            mean = inf_age_means[known_endo_indices])
saveRDS(endo_TI_means, file = paste0("endo_TI_means.", prior, ".", date_tag, ".rds"))

risks.grouped <- c("HET", "MSM", "IDU", "unknown/other")
intervals <- 2003:2022
indices <- list()
n_indices <- list()
#means_of_means <- list()
#sds_of_means <- list()
#means_of_sds <- list()
#sds_of_sds <- list()
drawn_infection_ages <- list()
mean_infection_ages_samples <- list()
mean_infection_ages <- list()
quantile_infection_ages <- list()
#number of samples per individual
n_samples_per_ind <- 10000
for(i in seq_along(risks.grouped)){
  indices[[i]] <- list()
  n_indices[[i]] <- numeric(length = length(intervals))
  drawn_infection_ages[[i]] <- list()
  mean_infection_ages_samples[[i]] <- list()
  mean_infection_ages[[i]] <- list()
  quantile_infection_ages[[i]] <- list()
  #means_of_means[[i]] <- numeric(length = length(intervals))
  #sds_of_means[[i]] <- numeric(length = length(intervals))
  #means_of_sds[[i]] <- numeric(length = length(intervals))
  #sds_of_sds[[i]] <- numeric(length = length(intervals))
  for(j in seq_along(intervals)){
    if(i != 4) indices[[i]][[j]] <- which(infection.ages.informed$risk_group == risks.grouped[i] & 
                                            infection.ages.informed$diagnosis_year == intervals[j] & infection.ages.informed$endoexo == "endo")
    else indices[[i]][[j]] <- which(!(infection.ages.informed$risk_group %in% risks.grouped[1:3]) & 
                                      infection.ages.informed$diagnosis_year == intervals[j] & infection.ages.informed$endoexo == "endo") 
    n_indices[[i]][[j]] <- length(indices[[i]][[j]])
    #drawn_indices <- indices[sample(n_indices[[i]][[j]], n_samples_per_int, replace = TRUE)]
    #drawn_infection_ages[[i]][[j]] <- mapply(FUN = function(icdf, u) icdf(u),
    #                                         icdf = inf_icdfs_informed[drawn_indices],
    #                                         u = runif(n_samples_per_int))
    #draw one infection time from each individual for each sample 
    drawn_infection_ages[[i]][[j]] <- simplify2array(lapply(inf_icdfs_informed[indices[[i]][[j]]], 
                                                            FUN = function(icdf, n) icdf(runif(n)),
                                                            n = n_samples_per_ind))
    #find the mean infection age (across individuals) for each sample
    if(n_indices[[i]][[j]] > 0){
      mean_infection_ages_samples[[i]][[j]] <- rowMeans(drawn_infection_ages[[i]][[j]])
      mean_infection_ages[[i]][[j]] <- mean(mean_infection_ages_samples[[i]][[j]])
      quantile_infection_ages[[i]][[j]] <- quantile(mean_infection_ages_samples[[i]][[j]], probs = c(0.025, 0.975))
    } else{
      mean_infection_ages_samples[[i]][[j]] <- NA
      mean_infection_ages[[i]][[j]] <- NA
      quantile_infection_ages[[i]][[j]] <- quantile(mean_infection_ages_samples[[i]][[j]], probs = c(0.025, 0.975), na.rm = TRUE)
    }
    #means_of_means[[i]][[j]] <- mean(inf_age_means[indices[[i]][[j]]])
    #sds_of_means[[i]][[j]] <- sd(inf_age_means[indices[[i]][[j]]])
    #means_of_sds[[i]][[j]] <- mean(inf_age_sd[indices[[i]][[j]]])
    #sds_of_sds[[i]][[j]] <- sd(inf_age_sd[indices[[i]][[j]]])
  }
  mean_infection_ages[[i]] <- unlist(mean_infection_ages[[i]])
  quantile_infection_ages[[i]] <- simplify2array(quantile_infection_ages[[i]])
}

mean_df <- function(years, mean, quantile, n){
  #return(rbind(mean, quantile))
  df <- as.data.frame(t(rbind(mean, quantile)))
  names(df) <- c("Mean", "Low2.5", "High97.5")
  df$Year <- years
  df$nInds <- n
  return(df)
}

mean_infection_ages_df <- mapply(FUN = mean_df, 
                                 years = list(intervals),
                                 mean = mean_infection_ages,
                                 quantile = quantile_infection_ages, 
                                 n = n_indices,
                                 SIMPLIFY = FALSE)
names(mean_infection_ages_df) <- risks.grouped

yearly.infection.ages <- list(mean_df = mean_infection_ages_df,
                              samples = mean_infection_ages_samples)

saveRDS(yearly.infection.ages, file = paste0("mean.infection.ages.by.year.", date_tag, ".rds"))

#plot(intervals, means_of_means[[1]], type = 'l', col = "black",
#     ylim = c(0, 3))
#lines(intervals, means_of_means[[2]], col = "blue")
#lines(intervals, means_of_means[[3]], col = "red")
#lines(intervals, means_of_means[[4]], col = "green")

#plot(intervals, means_of_means[[1]], pch = 16, col = "black",
#     ylim = c(0, 5))
#points(intervals, means_of_means[[1]]+means_of_sds[[1]], pch = 16, col = "black")
#points(intervals, means_of_means[[1]]-means_of_sds[[1]], pch = 16, col = "black")

###find distribution for TI for individuals diagnosed during or after 2015
#find indices of known-endo individuals diagnosed 2015 or later
post_2015_indices_endo <- which(infection.ages.informed$diagnosis_year >= 2015 & infection.ages.informed$endoexo == "endo")
#list for indices of known-endo individuals diagnosed 2015 or later for each risk group
post_2015_indices <- vector(mode = "list", length = length(risks.grouped))
#samples of infection ages for those individuals
draws_2015_plus <- vector(mode = "list", length = length(risks.grouped))
for(i in seq_along(risks.grouped)){
  post_2015_indices[[i]] <- intersect(post_2015_indices_endo,
                                      which(infection.ages.informed$risk_group == risks.grouped[i]))
  draws_2015_plus[[i]] <- unlist(inf_age_draws[post_2015_indices[[i]]])
}
names(draws_2015_plus) <- risks.grouped[1:3]
saveRDS(draws_2015_plus, file = "TI_2015_plus.rds")
sapply(draws_2015_plus, FUN = quantile, probs = c(0.5, 0.75, 0.95))

###make table of summary stats for infection distributions

#find median, 2.5th and 97.5th percentiles for each patient's time since infection
TI_quantiles <- simplify2array(lapply(inf_icdfs, 
                                      FUN = function(icdf, p){icdf(p)},
                                      p = c(0.025, 0.5, 0.975)))

TI_means <- sapply(infection.ages$infection.age, FUN = find.pdf.mean)

TI_summary_stats <- data.frame(patient_ID = infection.ages$patient_ID,
                            diagnosis_date = infection.ages$diagnosis_date,
                            mean = TI_means,
                            median = TI_quantiles[2,],
                            quantile2.5 = TI_quantiles[1,],
                            quantile97.5 = TI_quantiles[3,])

saveRDS(TI_summary_stats, file = "TI_summary_stats.rds")
write.csv(TI_summary_stats, file = "TI_summary_stats.csv", row.names = FALSE)
