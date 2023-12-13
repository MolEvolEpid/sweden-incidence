#prepare data for use in MBM

library(biophybreak)

#to protect patient privacy, the following input files are not publicly available
data <- readRDS("HIV_Biomarker_data_all.rds")
clinic_dates <- read.csv(file = "Thomas_clinics.csv")
leave_dates <- read.csv(file = "Thomas_avskrivna_datum_230516.csv")
diag_in_swe_dates <- read.csv(file = "Thomas_first_diagnosis_Sweden 230516.csv")

#codes of HIV-2 patients
HIV2_codes <- read.csv("HIV2_koder.csv")
HIV2_all_indices <- which(data$Code %in% HIV2_codes$Kod)

#codes of patients with no diagnosed-in-sweden date
no_diag_in_swe <- diag_in_swe_dates$Kod[which(diag_in_swe_dates$First.diagnosis.Sweden == "")]
no_diag_in_swe_indices <- which(data$Code %in% no_diag_in_swe)

#are any of these overlapping?
print(length(intersect(HIV2_all_indices, no_diag_in_swe_indices)))

#number of patients still in care at end of study (April 4th 2023)
print(dim(data)[1] - length(unique(leave_dates$Kod[leave_dates$Kod %in% data$Code])))

#date to use to tag file output
date_tag <- "20231024"

#find patients
pat_ids <- unique(data$Code)
#number of patients
nPat <- length(unique(pat_ids))

#format data to be used in prepare.HIV.data (can ignore warnings about unknown timezone)
last.neg <- as.Date(data$`Last neg HIV test`, "%Y-%m-%d")
first.art <- as.Date(data$`ART start`, "%Y-%m-%d")
first.pos <- as.Date(data$`First pos HIV test`, "%Y-%m-%d")
first.pos.swe <- as.Date(diag_in_swe_dates$First.diagnosis.Sweden, "%Y-%m-%d")
first.VL <- as.Date(data$`First VL date`, "%Y-%m-%d")
#CD4.dates <- as.Date(integer(nPat), origin = "1970-01-01")
CD4.dates <- lapply(data$CD4_samples, FUN = function(x) as.Date(x[[1]], "%Y-%m-%d"))
CD4 <- lapply(data$CD4_samples, FUN = function(x) x[[2]])
seq.dates <- lapply(data$sequences, FUN = function(x) as.Date(x[[1]], "%Y-%m-%d"))
seqs <- lapply(data$sequences, FUN = function(x) x[[2]])

#arrival date
arrival.date <- as.Date(data$`Arrival date`, "%Y-%m-%d")

#removal date
removal.date <- rep(NA, nPat)
removal.date[leave_dates$Kod] <- as.Date(leave_dates$EndDate, "%Y-%m-%d")

#gender
gender.swe <- data$Gender
#translate
translate_gender <- function(gender){
  if(!is.na(gender)){
    gender.translated <- switch(gender, "Man" = "M", "Kvinna" = "F")
  } else{
    gender.translated <- gender #keep NA as NA
  }
}
gender <- unname(sapply(gender.swe, FUN = translate_gender))

#translate risk groups
translate.risks <- function(risks){
  translated <- character(length = length(risks))
  for(i in seq_along(risks)){
    if(is.na(risks[i]) || risks[i] == "") translated[i] <- "unknown"
    else translated[i] <- switch(risks[i],
                                 "Heterosexuell" = "HET",
                                 "Homo/bisexuell" = "MSM",
                                 "I.V. missbruk" = "IDU",
                                 "Blodprodukter" = "BPD",
                                 "Okänd/Övrig" = "unknown",
                                 "Mor-Barn" = "MTC")
  }
  return(translated)
}
risk.groups <- translate.risks(data$Route)

#patients born in sweden or with a previous negative test after arrival are labeled endo; others are labeled exo
#patients with no birth country are endo unless they have an arrival date, in which case they are labeled exo
#(patients labeled exo may become at least partially endo depending on how the infection age distributions compare to the arrival dates)
endoexo <- rep("exo", nPat)
endoexo[data$`Country birth` == "Sverige"] <- "endo"
endoexo[which(last.neg >= arrival.date)] <- "endo"
endoexo[is.na(data$`Country birth`) & is.na(arrival.date)] <- "endo"

#function to adjust first positive dates based on ART start date
adjust_first_pos_art <- function(first.pos, first.art){
  if(!is.na(first.art)){
    if(is.na(first.pos)){
      first.pos.art.adj <- first.art
    } else if(first.art <= first.pos){
      first.pos.art.adj <- first.art
    }  else{
      first.pos.art.adj <- first.pos
    }
  } else{
    first.pos.art.adj <- first.pos
  }
  return(first.pos.art.adj)
}

#function to adjust first positive dates based on first viral load date
adjust_first_pos_VL <- function(first.pos, first.VL){
  if(!is.na(first.VL)){
    if(is.na(first.pos)){
      first.pos.VL.adj <- first.VL
    } else if(first.VL <= first.pos){
      first.pos.VL.adj <- first.VL
    }  else{
      first.pos.VL.adj <- first.pos
    }
  } else{
    first.pos.VL.adj <- first.pos
  }
  return(first.pos.VL.adj)
}

#function to adjust first positive dates based on sequence times
adjust_first_pos_seq <- function(first.pos, seq.dates){
  #find whether we have any sequences
  if(length(seq.dates) == 0){
    return(first.pos) #use first positive date if no sequences
  }
  #find date of first sequence
  first.seq.date <- min(seq.dates)
  #use first seq date if there is no first positive date
  if(is.na(first.pos)){
    first.pos.seq.adj <- first.seq.date
  } else if(first.seq.date < first.pos){
    if(first.seq.date - first.pos >= 366 | first.seq.date - first.pos < 345){
      first.pos.seq.adj <- first.seq.date
    } else{
      #use 365 days before first pos date if first seq is between 345 and 365 days before first pos
      first.pos.seq.adj <- first.pos - 365
    }
  } else{
    first.pos.seq.adj <- first.pos
  }
  return(first.pos.seq.adj)
}

first.pos.seq.adj <- mapply(FUN = adjust_first_pos_seq, 
                            first.pos, 
                            seq.dates)

first.pos.art.adj <- mapply(FUN = adjust_first_pos_art,
                            first.pos.seq.adj, 
                            first.art)

first.pos.VL.adj <- mapply(FUN = adjust_first_pos_VL,
                            first.pos.art.adj, 
                            first.art)

#function to adjust first positive dates based on first CD4 time (only if first pos is missing)
adjust_first_pos_cd4 <- function(first.pos, CD4.dates){
  if(length(CD4.dates) >= 1){
    if(is.na(first.pos)){
      first.pos.cd4.adj <- min(CD4.dates)
    } else{
      first.pos.cd4.adj <- first.pos
    }
  } else{
    first.pos.cd4.adj <- first.pos
  }
  return(first.pos.cd4.adj)
}

first.pos.cd4.adj <- mapply(FUN = adjust_first_pos_cd4,
                            first.pos.VL.adj, 
                            CD4.dates)


#plot(as.double(first.pos.cd4.adj)/365.25+1970,
#     as.double(first.pos.swe)/365.25+1970, col = "#00000020")
#sum(as.double(first.pos.cd4.adj)/365.25+1970 ==
#    as.double(first.pos.swe)/365.25+1970, na.rm = TRUE)
#sum((as.double(first.pos.cd4.adj)/365.25+1970) -
#      (as.double(first.pos.swe)/365.25+1970) < -2, na.rm = TRUE)

#remove last negative tests that occur on the same day or after the first positive date
last.neg[last.neg >= first.pos.cd4.adj] <- NA

#remove CD4 measurements that occur before the first positive test
remove_early_cd4 <- function(first.pos, CD4.dates, CD4){
  CD4_before_pos <- which(CD4.dates >= first.pos)
  CD4.dates <- CD4.dates[CD4_before_pos]
  CD4 <- CD4[CD4_before_pos]
  return(list(CD4.dates, CD4))
}
CD4_pre_qc <- mapply(FUN = remove_early_cd4,
                     first.pos.cd4.adj, 
                     CD4.dates,
                     CD4, 
                     SIMPLIFY = FALSE)
CD4.dates.pqc <- lapply(CD4_pre_qc, FUN = function(x) x[[1]])
CD4.pqc <- lapply(CD4_pre_qc, FUN = function(x) x[[2]])

#function to calculate pol externally
calculate_pol <- function(sequences){
  nseq <- length(sequences)
  length_seq <- integer(length = nseq)
  pol <- double(length = nseq)
  for(i in seq_len(nseq)){
    seq.split <- strsplit(tolower(as.character(sequences[[i]])), split = "")
    frequencies <- table(seq.split)
    polymorphic_count <- sum(frequencies[names(frequencies) != "-" & names(frequencies) != "n" & names(frequencies) != "a" & 
                                           names(frequencies) != "c" & names(frequencies) != "g" & names(frequencies) != "t"])
    length_seq[i] <- sum(frequencies[names(frequencies) != "-" & names(frequencies) != "n"]) #length of sequence excluding gaps and Ns
    pol[i] <- polymorphic_count/length_seq[i]
    #change NaN to NA
    if(is.na(pol[i])) pol[i] <- NA
  }
  return(list(length = length_seq, pol = pol))
}

length_and_pol <- lapply(seqs, FUN = calculate_pol)
#seq_length <- lapply(length_and_pol, FUN = function(x) x$length)
pol <- lapply(length_and_pol, FUN = function(x) x$pol)
seq_lengths <- lapply(length_and_pol, FUN = function(x) x$length)

#additional quality control
#number of seqs per patient
nseqs <- sapply(seq.dates, length)
nseqs_unique <- sapply(seq.dates, FUN = function(x) length(unique(x)))
#copy polymorphisms so all but one value from repeated measurement days can be removed
pol_no_repeat <- pol
#find patients with multiples sequences on the same date
multiseq_pats <- which(nseqs != nseqs_unique)
for(i in multiseq_pats){
  #loop over unique sequence sample dates
  #if there are multiple from the same date, only use the one with the highest pol value
  seq_date_unique <- unique(seq.dates[[i]])
  for(j in seq_len(nseqs_unique[i])){
    same_date_indices <- which(seq.dates[[i]] == seq_date_unique[j])
    #find index of maximum pol value
    pol_max_ind <- which.max(pol_no_repeat[[i]][same_date_indices])
    #make pol values that are not the max NA
    pol_no_repeat[[i]][same_date_indices[-pol_max_ind]] <- NA
  }
}

#find indices of patients with no first positive date
no_first_pos_inds <- which(is.na(first.pos.cd4.adj))

#will give warnings about rbind on rows that are not the same length, but we are using anything from that
df <- prepare.HIV.data(patient_ID = pat_ids[-no_first_pos_inds], 
                       last_neg_test_date = last.neg[-no_first_pos_inds],
                       first_pos_test_date = as.Date(first.pos.cd4.adj[-no_first_pos_inds], origin = "1970-01-01"),
                       ART_start_date = first.art[-no_first_pos_inds],
                       CD4_dates = CD4.dates.pqc[-no_first_pos_inds],
                       CD4 = CD4.pqc[-no_first_pos_inds],
                       seq_dates = seq.dates[-no_first_pos_inds],
                       seqs = seqs[-no_first_pos_inds],
                       pol_override = pol_no_repeat[-no_first_pos_inds],
                       seq_length_override = seq_lengths[-no_first_pos_inds],
                       VL_dates = first.VL[-no_first_pos_inds], 
                       VL = data$`First VL`[-no_first_pos_inds],
                       cluster_ID = rep("1", nPat-length(no_first_pos_inds)), #just a placeholder
                       date_format = "%Y-%m-%d", 
                       gender = gender[-no_first_pos_inds],
                       birth_location = data$`Country birth`[-no_first_pos_inds],
                       suspected_infection_location = data$`Country transmission`[-no_first_pos_inds],
                       risk_group = risk.groups[-no_first_pos_inds],
                       arrival_date = as.double(arrival.date[-no_first_pos_inds])/365.25+1970,
                       removal_date = as.double(removal.date[-no_first_pos_inds])/365.25+1970,
                       diag_in_swe_date = as.double(first.pos.swe[-no_first_pos_inds])/365.25+1970,
                       endoexo = endoexo[-no_first_pos_inds],
                       find_infection_age_distributions = FALSE)

#check to make sure no additional changes to the first positive date have been made
print(all(df$first_pos_test_date_adj == as.double(first.pos.cd4.adj[-no_first_pos_inds])/365.25+1970, na.rm = TRUE))

#change arrival_date and removal date to double
df$arrival_date <- as.double(df$arrival_date)
df$removal_date <- as.double(df$removal_date)

#find patients with no data
#no.data.ind <- which(is.infinite(df$first_pos_test_date_adj))

#add a date and NA value to patients with length 0 CD4 data
empty.CD4 <- which(sapply(df$CD4_dates, FUN = function(x) (length(x) == 0)))
df$CD4_dates[empty.CD4] <- df$first_pos_test_date_adj[empty.CD4]
df$CD4_qc[empty.CD4] <- NA

#add a date and NA value to patients with length 0 pol data
empty.pol <- which(sapply(df$seq_dates, FUN = function(x) (length(x) == 0)))
df$seq_dates[empty.pol] <- df$first_pos_test_date_adj[empty.pol]
df$pol_qc[empty.pol] <- NA

#find indices of HIV-2 patients
HIV2_indices <- which(df$patient_ID %in% HIV2_codes$Kod)
#remove HIV-2 patients
df.no.HIV2 <- df[-HIV2_indices,]

#find indices of patients with no-non NA values for CD4 and pol
find.no.data.ind <- function(df){
  nPat <- dim(df)[1]
  no.data <- logical(nPat)
  for(i in seq_len(nPat)){
    if(all(is.na(df$CD4_qc[[i]])) && all(is.na(df$pol_qc[[i]]))) no.data[i] <- TRUE
    else no.data[i] <- FALSE
  }
  return(which(no.data))
}
no.data.ind <- find.no.data.ind(df.no.HIV2)
#indices with sufficient data
data.ind <- seq_len(dim(df.no.HIV2)[1])[-no.data.ind]
saveRDS(data.ind, file = paste0("informed_indices.", date_tag, ".rds"))

df.informed <- df.no.HIV2[-no.data.ind,]
saveRDS(df.informed, file = paste0("df.incidence.informed.", date_tag, ".rds"))

#df.qc <- df[-no.data.ind,]
#make sure no sequence dates are earlier than QC'd first positive dates
print(length(which(sapply(df$seq_dates, min, na.rm = TRUE) < df$first_pos_test_date_adj)))
#make sure no last negative test dates are later than QC'd first positive dates
print(length(which(df$last_neg_test_date > df$first_pos_test_date_adj)))

#save quality controlled data
saveRDS(df.no.HIV2, file = paste0("df.incidence.new.", date_tag, ".Rdata"))

#find number of patients remaining in care at end of study
print(dim(df.no.HIV2)[1] - length(which(df.no.HIV2$patient_ID %in% leave_dates$Kod)))

#find number of non-SWE-born
n_non_swe_born <- dim(df.no.HIV2)[1] - sum(df.no.HIV2$birth_location == "Sverige", na.rm = TRUE)
print(n_non_swe_born)

#find number of patients with an arrival date
n_arrival <- sum(!is.na(df.no.HIV2$arrival_date))
print(n_arrival)

#percent of non-SWE-born patients with an arrival date
print(n_arrival/n_non_swe_born)

#number diagnosed (non-strictly) after known arrival date
diag_after_arrival <- sum(df.no.HIV2$first_pos_test_date_adj >= df.no.HIV2$arrival_date, na.rm = TRUE)
print(diag_after_arrival)

#percent diagnosed after arrival
print(diag_after_arrival/n_arrival)

#infer TI (time between infection and diagnosis) distributions
infection.distributions <- run.mbm(df.informed, 
                                   n.adapt = 1e4,
                                   n.burn = 1e5, 
                                   n.iter = 1e6, 
                                   prior.type = 1, #1 is gamma distribution
                                   overall.seed = 202303240)
saveRDS(infection.distributions, file = paste0("infection.distributions.only.gamma.", date_tag, ".rds"))
