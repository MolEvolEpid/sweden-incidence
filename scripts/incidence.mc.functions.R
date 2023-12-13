#functions for monte carlo estimation of incidence

#function to calculate the probability of not being diagnosed given infection in a certain year
diagnosis.prob.fn <- function(typical.dists, tau1 = 2003, tau3 = 2021.8){
  #find breakpoints for years of interest (including partial last year)
  interval.points <- seq(tau1, tau3, by = 1)
  #add point for last or partial last year
  if(floor(tau3) - tau3 != 0) interval.points <- c(interval.points, tau3)
  else interval.points <- c(interval.points, tau3+1)
  
  #convert to cdfs
  typical.cdfs <- lapply(typical.dists, FUN = find.cdf)
  dim(typical.cdfs) <- dim(typical.dists)
  rownames(typical.cdfs) <- rownames(typical.dists)
  colnames(typical.cdfs) <- colnames(typical.dists)
  return(typical.cdfs)
}

#function to calculate typical arrival time densities
typical.arrival.dist <- function(infection.ages, risk_groups){
  #list for typical arrival distributions
  typical.arrival <- vector(mode = "list", length = length(risk_groups))
  for(i in seq_along(risk_groups)){
    indices <- which(infection.ages$risk_grouped == risk_groups[i])
    typical.arrival[[i]] <- density(infection.ages$diagnosis_date[indices] - infection.ages$arrival_date[indices], na.rm = TRUE)
  }
  names(typical.arrival) <- risk_groups
  return(typical.arrival)
}

#simple trapezoidal rule integration function for already-discretized functions
integrate.num <- function(f, lower, upper){
  #find indices of x values closest to lower and upper (rounding down)
  #lower
  if(lower <= min(f$x)) ind_l <- 1
  else if(lower >= max(f$x)) return(0) #must be 0 if lower bound if later than time with density
  else ind_l <- max(which(f$x <= lower))
  #upper
  if(upper <= min(f$x)) return(0) ##must be 0 if upper bound is earlier than time with density
  else if(upper >= max(f$x)) ind_u <- length(f$x)
  else ind_u <- max(which(f$x <= upper))
  #check to make sure it makes sense
  if(ind_l > ind_u) stop(paste0(print(ind_l)," ", print(ind_u)," ","lower bound higher than upper bound"))
  else indices <- ind_l:ind_u
  
  #find differences in x values
  diffs <- diff(f$x[indices])
  #find means of y values
  ymids <- numeric(length = length(diffs))
  for(i in seq_along(ymids)){
    ymids[i] <- (f$y[indices[i]]+f$y[indices[i+1]])/2
  }
  return(sum(diffs*ymids))
}

#function to make a cdf from the pdf
find.cdf.num <- function(pdf, tolerance = 0.01){
  #find maximum time
  max_time <- max(pdf$x)
  
  #approximation function of pdf
  pdf_fn <- approxfun(pdf$x, pdf$y, method = "linear", yleft = 0, yright = 0)
  
  integral.step.vals <- rep(0, length = length(pdf$x))
  for(i in seq_along(pdf$x)[-1]){
    integral.step.vals[i] <- integrate(f = pdf_fn, lower = pdf$x[i-1], upper = pdf$x[i])$value
  }
  cdf.num <- cumsum(integral.step.vals)
  if(max(cdf.num) > 1+tolerance || max(cdf.num) < 1-tolerance) stop("CDF error tolerance exceeded")
  #normalize if within tolerance
  cdf.num <- cdf.num/max(cdf.num)
  
  return(list(x = pdf$x, y = cdf.num))
}

#find unique values in numerical cdf (needed to make inverse function)
find.unique.cdf.values <- function(cdf.num){
  uniques <- sapply(unique(cdf.num$y), FUN = function(unique_cdf_y, cdf.num){
    which(unique_cdf_y == cdf.num$y)[1]}, 
    cdf.num) #first of each value
  cdf.num.unique <- list(x = cdf.num$x[uniques], y = cdf.num$y[uniques])
}

#function to make an inverse cdf function
find.icdf <- function(pdf){
  cdf.num <- find.cdf.num(pdf, tolerance = 0.01)
  cdf.num.unique <- find.unique.cdf.values(cdf.num)
  icdf <- approxfun(x = cdf.num.unique$y, y = cdf.num.unique$x, 
                    method = "linear", yleft = min(cdf.num.unique$x), yright = max(cdf.num.unique$x))
  return(icdf)
}

#function to make a (regular)cdf function from the pdf
find.cdf <- function(pdf, tolerance = 0.01){
  #find maximum time
  max_time <- max(pdf$x)
  
  #approximation function of pdf
  pdf_fn <- approxfun(pdf$x, pdf$y, method = "linear", yleft = 0, yright = 0)
  
  integral.step.vals <- rep(0, length = length(pdf$x))
  for(i in seq_along(pdf$x)[-1]){
    integral.step.vals[i] <- integrate(f = pdf_fn, lower = pdf$x[i-1], upper = pdf$x[i])$value
  }
  cdf.num <- cumsum(integral.step.vals)
  if(max(cdf.num) > 1+tolerance || max(cdf.num) < 1-tolerance) stop("CDF error tolerance exceeded")
  #normalize if within tolerance
  cdf.num <- cdf.num/max(cdf.num)
  
  #find approximation function
  cdf_fn <- approxfun(x = pdf$x, y = cdf.num, method = "linear", yleft = 0, yright = 1)
  return(cdf_fn)
}

date2decimal <- function(date){
  as.double(as.Date(date, "%Y-%m-%d"))/365.25 + 1970
}

decimal2year <- function(year.decimal){
  date <- as.Date((year.decimal-1970)*365.25, origin = "1970-01-01")
  year <- as.integer(sapply(date, FUN = substr, start = 1, stop = 4))
}

#monte carlo function to estimate incidence
estimate.incidence.mc <- function(infection.ages, 
                                  inf_icdfs, 
                                  diagnosis.prob.fns, 
                                  typical.arrival.dists, 
                                  tau1,
                                  tau3,
                                  n_samples){
  #find risk groups
  risk_groups <- rownames(diagnosis.prob.fns)
  
  #check to see if tau1 is a whole number
  if(floor(tau1) - tau1 != 0) stop("Only whole year numbers are supported for tau1")
  #find breakpoints for years of interest (including partial last year)
  interval.points <- seq(tau1, tau3, by = 1)
  #add point for last or partial last year
  if(floor(tau3) - tau3 != 0) interval.points <- c(interval.points, tau3)
  else interval.points <- c(interval.points, tau3+1)
  
  #find icdfs for arrival distributions
  arrival.icdfs <- lapply(typical.arrival.dists, FUN = find.icdf)
  
  #find indices for endo group
  endo_indices <- which(infection.ages$endoexo == "endo")
  #find indices for exo group
  exo_indices <- which(infection.ages$endoexo == "exo")
  
  #find indices of exo individuals with a known arrival date
  known_arrival_indices <- which(infection.ages$endoexo == "exo" & !is.na(infection.ages$arrival_date))
  #find indices of exo individuals with a known arrival date
  unknown_arrival_indices <- which(infection.ages$endoexo == "exo" & is.na(infection.ages$arrival_date))
  
  #make matrix for diagnosis years (since effective diagnosis year may be the arrival year)
  diagnosis_years <- matrix(infection.ages$diagnosis_year, nrow = dim(infection.ages)[1], ncol = n_samples)
  
  #function or random draw arrival times (or use actual ones if available)
  draw.arrivals <- function(arrival.icdf, diagnosis.dates, entry.times, entry.years, n_samples, tau3){
    arrival.times <- matrix(nrow = length(entry.years), ncol = n_samples)
    arrival.years <- matrix(nrow = length(entry.years), ncol = n_samples)
    for(i in seq_along(entry.years)){
      if(is.na(entry.years[i])){
        arrival.times[i,] <- diagnosis.dates[i] - arrival.icdf(runif(n_samples))
        while(any(arrival.times[i,] > tau3)){
          impossible_arrivals <- (arrival.times[i,] > tau3)
          arrival.times[i,][impossible_arrivals] <- diagnosis.dates[i] - arrival.icdf(runif(sum(impossible_arrivals)))
        }
        arrival.years[i,] <- floor(arrival.times[i,])
      } else{
        arrival.times[i,] <- rep(entry.times[i], n_samples)
        arrival.years[i,] <- rep(entry.years[i], n_samples)
      }
    }
    return(list(times = arrival.times, years = arrival.years))
  }
  
  #function to tabulate events in years of interest
  tab.years <- function(year, intervals){
    table <- tabulate(floor(year), nbins = max(intervals))[intervals]
    names(table) <- intervals
    return(table)
  }
  
  tab.years.weight <- function(infection.times, intervals, diagnosis.dates, diagnosis.prob.fn, tau3){
    #find weights
    #weights <- 1/matrix(diagnosis.prob.fn(tau3 - infection.times), nrow = nrow(infection.times))
    weights <- 1 + matrix(rgeom(length(infection.times), diagnosis.prob.fn(tau3 - infection.times)), nrow = nrow(infection.times))
    #incidence of diagnosed individuals
    incidence_diag <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(incidence_diag) <- intervals
    #total incidence
    incidence <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(incidence) <- intervals
    for(i in seq_along(intervals)){
      #find individuals infected the ith interval
      inf.this.int <- (floor(infection.times) == intervals[i])
      incidence_diag[i,] <- colSums(inf.this.int)
      incidence[i,] <- colSums(inf.this.int*weights)
    }
    return(list(incidence_diag = incidence_diag, incidence = incidence))
  }
  
  tab.years.weight.exo <- function(infection.times, arrival.times, arrival.years, intervals, 
                                   diagnosis.dates, diagnosis.years, diagnosis.mat,
                                   removal.years, diagnosis.prob.fn, tau3){
    #find weights
    #weights <- 1/matrix(diagnosis.prob.fn(tau3 - infection.times), nrow = nrow(infection.times))
    weights <- 1 + matrix(rgeom(length(infection.times), diagnosis.prob.fn(tau3 - infection.times)), nrow = nrow(infection.times))
    if(any(weights == Inf)) print(which(weights == Inf))
    if(any(infection.times >= tau3)) print(which(infection.times >= tau3))
    #find which infection times are after arrival times for each individual for each sample
    endos <- (infection.times > arrival.times)
    #adjust diagnosis times if they are earlier than the arrival times
    diagnosis.mat[diagnosis.mat < arrival.years] <- 
      arrival.years[diagnosis.mat < arrival.years]
    #make removal times into matrix
    removal.mat <- matrix(removal.years, nrow = nrow(infection.times), ncol = ncol(infection.times))
    #density of "exo" individuals that represent endogenous infection (diagnosed individuals)
    incidence_endo_diag <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(incidence_endo_diag) <- intervals
    #density of "exo" individuals that represent endogenous infection (total incidence)
    incidence_endo <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(incidence_endo) <- intervals
    diagnosed_endo <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(diagnosed_endo) <- intervals
    removed_endo <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(removed_endo) <- intervals
    for(i in seq_along(intervals)){
      #find individuals infected the ith interval (that are also infected after arrival)
      inf.this.int <- (floor(infection.times*endos) == intervals[i])
      incidence_endo_diag[i,] <- colSums(inf.this.int)
      incidence_endo[i,] <- colSums(inf.this.int*weights)
      #find individuals diagnosed in the ith interval
      diagnosed.this.int <- (diagnosis.mat*endos == intervals[i])
      diagnosed_endo[i,] <- colSums(diagnosed.this.int)
      #find individuals removed from "currently infected in swe" the ith interval
      removed.this.int <- (removal.mat*endos == intervals[i])
      removed_endo[i,] <- colSums(removed.this.int, na.rm = TRUE)
    }
    #density of exo individuals inferred to be infected before arrival
    incidence_exo_diag <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(incidence_exo_diag) <- intervals
    #density of "exo" individuals that represent endogenous infection (total incidence)
    incidence_exo <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(incidence_exo) <- intervals
    removed_exo <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(removed_exo) <- intervals
    diagnosed_exo <- matrix(nrow = length(intervals), ncol = ncol(infection.times))
    rownames(removed_exo) <- intervals
    for(i in seq_along(intervals)){
      #find individual that arrived during this interval
      exo.arrived.this.int <- ((arrival.years*!endos) == intervals[i])
      incidence_exo_diag[i,] <- colSums(exo.arrived.this.int)
      incidence_exo[i,] <- colSums(exo.arrived.this.int*weights)
      #find individuals diagnosed in the ith interval
      diagnosed.this.int <- ((diagnosis.mat*!endos) == intervals[i])
      diagnosed_exo[i,] <- colSums(diagnosed.this.int)
      #find individuals removed from "currently infected in swe" the ith interval
      exo.removed.this.int <- ((removal.mat*!endos) == intervals[i])
      removed_exo[i,] <- colSums(exo.removed.this.int, na.rm = TRUE)
    }
    return(list(incidence_endo_diag = incidence_endo_diag, incidence_endo = incidence_endo, 
                removed_endo = removed_endo, diagnosed_endo = diagnosed_endo,
                incidence_exo_diag = incidence_exo_diag, incidence_exo = incidence_exo, 
                removed_exo = removed_exo, diagnosed_exo = diagnosed_exo))
  }
  
  #endo
  endo_incidence_diag <- vector(mode = "list", length = length(risk_groups))
  names(endo_incidence_diag) <- risk_groups
  endo_incidence <- vector(mode = "list", length = length(risk_groups))
  names(endo_incidence) <- risk_groups
  endo_diagnoses_mats <- vector(mode = "list", length = length(risk_groups))
  names(endo_diagnoses_mats) <- risk_groups
  endo_diagnoses <- vector(mode = "list", length = length(risk_groups))
  names(endo_diagnoses) <- risk_groups
  endo_removals <- vector(mode = "list", length = length(risk_groups))
  names(endo_removals) <- risk_groups
  for(i in seq_along(risk_groups)){
    #find indices of subgroup
    indices <- which(infection.ages$endoexo == "endo" & infection.ages$risk_grouped == risk_groups[i])
    endo_diagnoses_mats[[i]] <- diagnosis_years[indices,]
    endo_diagnoses[[i]] <- apply(endo_diagnoses_mats[[i]], 2, FUN = tab.years, intervals = interval.points[-length(interval.points)])
    infection.times <- t(mapply(FUN = function(icdf, diag, n_samples) diag - icdf(runif(n_samples)),
                                inf_icdfs[indices], 
                                infection.ages$diagnosis_date[indices],
                                n_samples))
    endo_inc <- tab.years.weight(infection.times = infection.times, 
                                 intervals = interval.points[-length(interval.points)],
                                 diagnosis.dates = infection.ages$diagnosis_date[indices],
                                 diagnosis.prob.fn = diagnosis.prob.fns[[i, "endo"]],
                                 tau3)
    endo_incidence_diag[[i]] <- endo_inc$incidence_diag
    endo_incidence[[i]] <- endo_inc$incidence
    endo_removals[[i]] <- tab.years(infection.ages$removal_year[indices], intervals = interval.points[-length(interval.points)])
    #make removals into matrices to be consistent with exo
    endo_removals[[i]] <- matrix(rep(endo_removals[[i]], n_samples), ncol = n_samples)
  }
  #exo (with or without known arrival date)
  exo_i_incidence_diag <- vector(mode = "list", length = length(risk_groups))
  names(exo_i_incidence_diag) <- risk_groups
  exo_i_incidence <- vector(mode = "list", length = length(risk_groups))
  names(exo_i_incidence) <- risk_groups
  exo_e_incidence_diag <- vector(mode = "list", length = length(risk_groups))
  names(exo_e_incidence_diag) <- risk_groups
  exo_e_incidence <- vector(mode = "list", length = length(risk_groups))
  names(exo_e_incidence) <- risk_groups
  arrival <- vector(mode = "list", length = length(risk_groups))
  exo_i_diagnoses <- vector(mode = "list", length = length(risk_groups))
  names(exo_i_diagnoses) <- risk_groups
  exo_e_diagnoses <- vector(mode = "list", length = length(risk_groups))
  names(exo_e_diagnoses) <- risk_groups
  exo_i_removals <- vector(mode = "list", length = length(risk_groups))
  names(exo_i_removals) <- risk_groups
  exo_e_removals <- vector(mode = "list", length = length(risk_groups))
  names(exo_e_removals) <- risk_groups
  for(i in seq_along(risk_groups)){
    indices <- which(infection.ages$endoexo == "exo" & infection.ages$risk_grouped == risk_groups[i])
    infection.times <- t(mapply(FUN = function(icdf, diag, n_samples) diag - icdf(runif(n_samples)),
                                inf_icdfs[indices], 
                                infection.ages$diagnosis_date[indices],
                                n_samples))
    exo_diagnoses_mats <- diagnosis_years[indices,]
    arrival[[i]] <- draw.arrivals(arrival.icdfs[[i]],
                                  infection.ages$diagnosis_date[indices],
                                  infection.ages$arrival_date[indices],
                                  infection.ages$arrival_year[indices],
                                  n_samples, 
                                  tau3)
    exo_inc <- tab.years.weight.exo(infection.times = infection.times, 
                                    arrival.times = arrival[[i]]$times,
                                    arrival.years = arrival[[i]]$years,
                                    intervals = interval.points[-length(interval.points)],
                                    diagnosis.dates = infection.ages$diagnosis_date[indices],
                                    diagnosis.years = infection.ages$diagnosis_year[indices],
                                    diagnosis.mat = exo_diagnoses_mats,
                                    removal.years = infection.ages$removal_year[indices],
                                    diagnosis.prob.fn = diagnosis.prob.fns[[i, "exo"]],
                                    tau3)
    exo_i_incidence_diag[[i]] <- exo_inc$incidence_endo_diag
    exo_i_incidence[[i]] <- exo_inc$incidence_endo
    exo_e_incidence_diag[[i]] <- exo_inc$incidence_exo_diag
    exo_e_incidence[[i]] <- exo_inc$incidence_exo
    exo_i_diagnoses[[i]] <- exo_inc$diagnosed_endo
    exo_e_diagnoses[[i]] <- exo_inc$diagnosed_exo
    exo_i_removals[[i]] <- exo_inc$removed_endo
    exo_e_removals[[i]] <- exo_inc$removed_exo
  }
  
  #add together both types of endo individuals
  endo_incidence_diag_both <- mapply(FUN = `+`, endo_incidence_diag, exo_i_incidence_diag, SIMPLIFY = FALSE)
  endo_incidence_both <- mapply(FUN = `+`, endo_incidence, exo_i_incidence, SIMPLIFY = FALSE)
  
  endo_diagnoses_both <- mapply(FUN = `+`, endo_diagnoses, exo_i_diagnoses, SIMPLIFY = FALSE)
  endo_removals_both <- mapply(FUN = `+`, endo_removals, exo_i_removals, SIMPLIFY = FALSE)
  
  #add risk groups together (still split by endo/exo)
  endo_incidence_diag_combined_risk <- Reduce(`+`, endo_incidence_diag_both)
  endo_incidence_combined_risk <- Reduce(`+`, endo_incidence_both)
  
  endo_diagnoses_combined_risk <- Reduce(`+`, endo_diagnoses_both)
  endo_removals_combined_risk <- Reduce(`+`, endo_removals_both)
  
  exo_incidence_diag_combined_risk <- Reduce(`+`, exo_e_incidence_diag)
  exo_incidence_combined_risk <- Reduce(`+`, exo_e_incidence)
  
  exo_diagnoses_combined_risk <- Reduce(`+`, exo_e_diagnoses)
  exo_removals_combined_risk <- Reduce(`+`, exo_e_removals)
  
  #add both endo and exo incidence (still split by risk group)
  incidence_diag_combined_ee <- mapply(FUN = `+`, endo_incidence_diag_both, exo_e_incidence_diag, SIMPLIFY = FALSE)
  incidence_combined_ee <- mapply(FUN = `+`, endo_incidence_both, exo_e_incidence, SIMPLIFY = FALSE)
  
  diagnoses_combined_ee <- mapply(FUN = `+`, endo_diagnoses_both, exo_e_diagnoses, SIMPLIFY = FALSE)
  removals_combined_ee <- mapply(FUN = `+`, endo_removals_both, exo_e_removals, SIMPLIFY = FALSE)
  
  #add everything together
  total_incidence_diag <- endo_incidence_diag_combined_risk + exo_incidence_diag_combined_risk
  total_incidence <- endo_incidence_combined_risk + exo_incidence_combined_risk
  
  total_diagnoses <- endo_diagnoses_combined_risk + exo_diagnoses_combined_risk
  total_removals <- endo_removals_combined_risk + exo_removals_combined_risk
  
  #function to calculate cumulative and current infections
  find_current <- function(incidence, diagnoses, removals){
    #cumulative incidence
    incidence_cum <- apply(incidence, 2, cumsum)
    #cumulative diagnoses
    diagnoses_cum <- apply(diagnoses, 2, cumsum)
    #cumulative removals
    removals_cum <- apply(removals, 2, cumsum)
    #current (for each year) totals of infections (accounting for removed individuals)
    infections_current <- incidence_cum - removals_cum
    #current (for each year) totals of diagnosed cases (accounting for removed individuals)
    diagnoses_current <- diagnoses_cum - removals_cum
    #current (for each year) totals of undiagnosed cases
    undiagnosed_current <- infections_current - diagnoses_current
    #current (for each year) fraction of diagnosed cases
    diagnosed_fraction_current <- diagnoses_current/infections_current
    #find means for certain quantities
    incidence_means <- rowMeans(incidence)
    incidence_cum_means <- rowMeans(incidence_cum)
    infections_current_means <- rowMeans(infections_current)
    undiagnosed_current_means <- rowMeans(undiagnosed_current)
    diagnosed_fraction_current_means <- rowMeans(diagnosed_fraction_current)
    #find 95% intervals for certain quantities
    incidence_95 <- apply(incidence, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    incidence_cum_95 <- apply(incidence_cum, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    infections_current_95 <- apply(infections_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    undiagnosed_current_95 <- apply(undiagnosed_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    diagnosed_fraction_current_95 <- apply(diagnosed_fraction_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    return(list(incidence = incidence, diagnoses = diagnoses, removals = removals,
                incidence_cum = incidence_cum, diagnoses_cum = diagnoses_cum, removals_cum = removals_cum,
                infections_current = infections_current, diagnoses_current = diagnoses_current, 
                undiagnosed_current = undiagnosed_current, diagnosed_fraction_current = diagnosed_fraction_current,
                incidence_means = incidence_means, incidence_cum_means = incidence_cum_means, 
                infections_current_means = infections_current_means,
                undiagnosed_current_means = undiagnosed_current_means, diagnosed_fraction_current_means = diagnosed_fraction_current_means,
                incidence_95 = incidence_95, incidence_cum_95 = incidence_cum_95, infections_current_95 = infections_current_95,
                undiagnosed_current_95 = undiagnosed_current_95, diagnosed_fraction_current_95 = diagnosed_fraction_current_95))
  }
  
  endo_split <- mapply(FUN = find_current, 
                       incidence = endo_incidence_both, 
                       diagnoses = endo_diagnoses_both,
                       removals = endo_removals_both, 
                       SIMPLIFY = FALSE)
  
  exo_split <- mapply(FUN = find_current, 
                      incidence = exo_e_incidence, 
                      diagnoses = exo_e_diagnoses,
                      removals = exo_e_removals, 
                      SIMPLIFY = FALSE)
  
  endo_combined_risk <- find_current(incidence = endo_incidence_combined_risk, 
                                     diagnoses = endo_diagnoses_combined_risk,
                                     removals = endo_removals_combined_risk)
  
  exo_combined_risk <- find_current(incidence = exo_incidence_combined_risk, 
                                    diagnoses = exo_diagnoses_combined_risk,
                                    removals = exo_removals_combined_risk) 
  
  combined_ee <- mapply(FUN = find_current, 
                        incidence = incidence_combined_ee, 
                        diagnoses = diagnoses_combined_ee,
                        removals = removals_combined_ee, 
                        SIMPLIFY = FALSE)
  
  total <- find_current(incidence = total_incidence, 
                        diagnoses = total_diagnoses,
                        removals = total_removals)
  
  return(list(endo_split = endo_split,
              exo_split = exo_split,
              endo_combined_risk = endo_combined_risk,
              exo_combined_risk = exo_combined_risk,
              combined_ee = combined_ee,
              total = total))
}
