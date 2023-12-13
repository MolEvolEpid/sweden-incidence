#script to read in individual infection age distributions
#should be run from the same location as prep.incidence.data.new.R

#date to use to tag file output
date_tag <- "20231024"

prior <- "gamma"
df <- readRDS(paste0("df.incidence.new.", date_tag, ".Rdata"))
data.ind <- readRDS(paste0("informed_indices.", date_tag, ".rds"))

#files <- list.files(indir)

infection.ages.informed <- data.frame(patient_ID = df$patient_ID[data.ind],
                                      last_negative = df$last_neg_test_date[data.ind],
                                      earliest_inftime = df$last_neg_test_date[data.ind] - 2/12,
                                      diagnosis_date = df$first_pos_test_date_adj[data.ind],
                                      ART_start_date = df$ART_start_date[data.ind],
                                      arrival_date = df$arrival_date[data.ind],
                                      removal_date = df$removal_date[data.ind],
                                      diag_in_swe_date = df$diag_in_swe_date[data.ind],
                                      endoexo = df$endoexo[data.ind],
                                      risk_group = df$risk_group[data.ind],
                                      infection.age = I(vector(mode = "list", length = length(data.ind))))

#used to read in output if run locally
infection.distributions <- readRDS(file = paste0("infection.distributions.only.gamma.", date_tag, ".rds"))
infection.ages.informed$infection.age <- infection.distributions$infection_age_dists_diag

saveRDS(infection.ages.informed, file = paste0("infection.ages.gamma.informed.new.", date_tag, ".rds"))

pdf(file = paste0("infection.ages.", prior, ".new.", date_tag, ".pdf"), width = 9, height = 5.0625)
for(i in seq_along(data.ind)){
  plot(infection.ages.informed$infection.age[[i]], main = infection.ages.informed$patient_ID[i])
}
dev.off()

#convert float dates to actual dates (to avoid Jan 1st being x.999 issue)
qc.dates <- as.Date((infection.ages.informed$diagnosis_date-1970)*365.25, origin = "1970-01-01")
year <- as.integer(sapply(qc.dates, FUN = substr, start = 1, stop = 4))

###find typical distributions
typical.dist <- function(dists){
  #find maximum 
  max_ages <- sapply(dists, FUN = function(dist) max(dist$x))
  max_age <- max(max_ages)
  #dense-ish timescale
  t <- seq(0, max_age, by = 0.001)
  #fit a linear approximation function to infection age distributions
  inf_age_fns <- lapply(dists, FUN = function(dist) approxfun(dist$x, dist$y, method = "linear", yleft = 0, yright = 0))
  #find typical distribution
  typical <- rep(0, length = length(t))
  integrals <- numeric(length = length(dists))
  for(i in seq_along(inf_age_fns)){
    integrals[i] <- integrate(inf_age_fns[[i]], lower = 0, upper = max_ages[i], stop.on.error = FALSE)$value
    typical <- typical + inf_age_fns[[i]](t)/integrals[i]
  }
  #make sure integrals are okay
  if(any(integrals < 0.95 | integrals > 1.05)){
    stop(paste0("integrals ", which(integrals < 0.95 || integrals > 1.05), "outside of [0.95,1.05]"))
  } 
  #normalize
  integral <- integrate(approxfun(t, typical, method = "linear", yleft = 0, yright = 0), 
                        lower = 0, upper = max_age, stop.on.error = FALSE, subdivisions = 1000)
  #return(integral)
  if(integral$value < 0.95*length(dists) || integral$value > 1.05*length(dists)) stop("typical dist integral outside of [0.95, 1.05")
  typical <- typical/integral$value
  return(list(x = t, y = typical))
}

#only include years under consideration
years.considered <- 2003:2023

#risk groups only for years under consideration
typical.dist.risk.considered <- list()
indices.risk.considered <- list()
risks.grouped <- c("HET", "MSM", "IDU", "unknown/other")
n_risk.considered <- integer(length = 4) #one for each of HET, MSM, and IDU and one other
for(i in seq_along(risks.grouped)){
  if(i != 4) indices.risk.considered[[i]] <- which(infection.ages.informed$risk_group == risks.grouped[i] & year %in% years.considered)
  else indices.risk.considered[[i]] <- which(year %in% years.considered) #use all individuals if unknown or other
  n_risk.considered[i] <- length(indices.risk.considered[[i]])
  typical.dist.risk.considered[[i]] <- typical.dist(infection.ages.informed$infection.age[indices.risk.considered[[i]]])
}

#endo/exo only for years under consideration
typical.dist.ee.considered <- list()
indices.ee.considered <- list()
ee <- c("endo", "exo")
n_ee.considered <- integer(length = 2) 
for(i in seq_along(ee)){
  indices.ee.considered[[i]] <- which(infection.ages.informed$endoexo == ee[i] & year %in% years.considered)
  n_ee.considered[i] <- length(indices.ee.considered[[i]])
  typical.dist.ee.considered[[i]] <- typical.dist(infection.ages.informed$infection.age[indices.ee.considered[[i]]])
}

#typical distributions for risk group and endo/exo for considered years
typical.dist.risk.ee.considered <- vector(mode = "list", length = length(risks.grouped)*2)
dim(typical.dist.risk.ee.considered) <- c(length(risks.grouped), 2)
indices.risk.ee.considered <- vector(mode = "list", length = length(risks.grouped)*2)
dim(indices.risk.ee.considered) <- c(length(risks.grouped), 2)
n_risk.ee.considered <- matrix(nrow = length(risks.grouped), ncol = 2)
for(i in seq_along(risks.grouped)){
  for(j in seq_along(ee)){
    indices.risk.ee.considered[[i,j]] <- intersect(indices.risk.considered[[i]], 
                                                   indices.ee.considered[[j]])
    n_risk.ee.considered[i,j] <- length(indices.risk.ee.considered[[i,j]])
    typical.dist.risk.ee.considered[[i,j]] <- typical.dist(infection.ages.informed$infection.age[indices.risk.ee.considered[[i,j]]])
  }
}
pdf(file = paste0("typical.dists.risk.ee.", prior, ".", date_tag, ".pdf"), width = 9, height = 5.0625)
for(i in seq_along(risks.grouped)){
  for(j in seq_along(ee)){
    plot(typical.dist.risk.ee.considered[[i,j]], type = 'l', xlim = c(0, 13), ylim = c(0,1),
         main = paste0(risks.grouped[i], ", ", ee[j], ", n = ", n_risk.ee.considered[i,j]))
  }
}
dev.off()

#name rows and columns
rownames(typical.dist.risk.ee.considered) <- risks.grouped
colnames(typical.dist.risk.ee.considered) <- ee

saveRDS(typical.dist.risk.ee.considered, file = paste0("typical.distributions.risk.ee.", prior, ".new.", date_tag, ".rds"))

#downsample typical distributions
downsample.typical.dists <- function(typical.dist, length.out = 512){
  fn <- approxfun(x = typical.dist$x, y = typical.dist$y, method = "linear")
  integral <- integrate(fn, lower = 0, upper = max(typical.dist$x))$value
  if(abs(1-integral) > 0.001) stop("downsampled integral differs by more than 0.001")
  times <- seq(0, max(typical.dist$x), length.out = length.out)
  return(list(x = times, y = fn(times)))
}

typical.dists.downsampled <- mapply(FUN = downsample.typical.dists,
                                    typical.dist.risk.ee.considered, 
                                    SIMPLIFY = FALSE)

dim(typical.dists.downsampled) <- c(length(risks.grouped),length(ee))
#name rows and columns
rownames(typical.dists.downsampled) <- risks.grouped
colnames(typical.dists.downsampled) <- ee

#use typical distributions to fill in distributions for individuals with no biomarkers
#full data frame
infection.ages <- data.frame(patient_ID = df$patient_ID,
                             last_negative = df$last_neg_test_date,
                             earliest_inftime = df$last_neg_test_date - 2/12,
                             diagnosis_date = df$first_pos_test_date_adj,
                             ART_start_date = df$ART_start_date,
                             arrival_date = df$arrival_date,
                             removal_date = df$removal_date,
                             diag_in_swe_date = df$diag_in_swe_date,
                             endoexo = df$endoexo,
                             risk_group = df$risk_group,
                             infection.age = I(vector(mode = "list", length = dim(df)[1])))


#group blood product, MTC, and unknown risk groups together
infection.ages$risk_grouped <- infection.ages$risk_group
infection.ages$risk_grouped[infection.ages$risk_grouped == "BPD"] <- "unknown/other"
infection.ages$risk_grouped[infection.ages$risk_grouped == "MTC"] <- "unknown/other"
infection.ages$risk_grouped[infection.ages$risk_grouped == "unknown"] <- "unknown/other"

#counter for informed distributions
k <- 0
for(i in seq_len(dim(infection.ages)[1])){
  if(i %in% data.ind){
    k <- k+1
    infection.ages$infection.age[[i]] <- infection.ages.informed$infection.age[[k]]
  } else{
    infection.ages$infection.age[[i]] <- typical.dists.downsampled[[infection.ages$risk_grouped[i], infection.ages$endoexo[i]]]
  }
}

saveRDS(infection.ages, file = paste0("infection.ages.", prior, ".all.new.", date_tag, ".rds"))

