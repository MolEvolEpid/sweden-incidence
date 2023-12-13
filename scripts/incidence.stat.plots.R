#date to use to tag file input/output
date_tag <- "20231024"

incidence <- readRDS(file = paste0("incidence.mc.gamma.", date_tag, ".more.stats.rds"))

#extract risk groups
risk_groups <- names(incidence$endo_split)
intervals <- 2003:2022

years = as.numeric(rownames(incidence$total$incidence))
considered_indices <- which(years >= intervals[1] & years <= tail(intervals, 1))
#considered_years <- years[considered_indices]

plot_stats <- function(subgroup, years, considered_indices = seq_along(years), title = ""){
  #shorten name for considered indices
  ci <- considered_indices
  plot(years[ci], subgroup$incidence_means[ci], 
       col = "black", pch = 16,
       xlab = "Year", ylab = "Estimated Number of New Infections",
       ylim = c(0, max(subgroup$incidence_95[2,][ci])),
       main = title)
  arrows(x0 = years[ci], y0 = subgroup$incidence_95[1,][ci], 
         x1 = years[ci], y1 = subgroup$incidence_95[2,][ci], 
         code = 3, angle = 90, length = 0.1, col = "black")
  
  plot(years[ci], subgroup$infections_current_means[ci], type = 'l',
       xlab = "Year", ylab = "Number of Individuals", 
       ylim = c(0, max(subgroup$infections_current_95, na.rm = TRUE)),
       main = title)
  points(years[ci], subgroup$infections_current_means[ci], pch = 16)
  lines(years[ci], rowMeans(subgroup$diagnoses_current[ci,]), col = "blue")
  points(years[ci], rowMeans(subgroup$diagnoses_current[ci,]), col = "blue", pch = 16)
  lines(years[ci], subgroup$undiagnosed_current_means[ci], col = "red")
  points(years[ci], subgroup$undiagnosed_current_means[ci], col = "red", pch = 16)
  legend("topleft", legend = c("current infections", "current diagnosed", "current undiagnosed"), 
         lty = c(1,1,1), col = c("black", "blue", "red"))
  
  plot(years[ci], subgroup$undiagnosed_current_means[ci], 
       col = "red", pch = 16,
       xlab = "Year", ylab = "Number of Undiagnosed Individuals",
       ylim = c(0, max(subgroup$undiagnosed_current_95[2,][ci])),
       main = title)
  arrows(x0 = years[ci], y0 = subgroup$undiagnosed_current_95[1,][ci], 
         x1 = years[ci], y1 = subgroup$undiagnosed_current_95[2,][ci], 
         code = 3, angle = 90, length = 0.1, col = "red")
  
  plot(years[ci], subgroup$diagnosed_fraction_current_means[ci], 
       col = "black", pch = 16,
       xlab = "Year", ylab = "Proportion of Current Infections Diagnosed",
       ylim = c(floor(min(subgroup$diagnosed_fraction_current_95[1,][ci])*20)/20, 1),
       main = title)
  arrows(x0 = years[ci], y0 = subgroup$diagnosed_fraction_current_95[1,][ci], 
         x1 = years[ci], y1 = subgroup$diagnosed_fraction_current_95[2,][ci], 
         code = 3, angle = 90, length = 0.1, col = "black")
  abline(h = 0.95, lty = 2) 
}

pdf(file = paste0("incidence_stats_with_subgroups_", date_tag, ".pdf"), width = 9, height = 5.0625)
plot_stats(subgroup = incidence$total, years = years, considered_indices = considered_indices, title = "Total")

plot_stats(subgroup = incidence$endo_split$HET, years = years, considered_indices = considered_indices, title = "Endo Het")
plot_stats(subgroup = incidence$endo_split$MSM, years = years, considered_indices = considered_indices, title = "Endo MSM")
plot_stats(subgroup = incidence$endo_split$IDU, years = years, considered_indices = considered_indices, title = "Endo IDU")
plot_stats(subgroup = incidence$endo_split$`unknown/other`, years = years, considered_indices = considered_indices, "Endo unknown/other")

plot_stats(subgroup = incidence$exo_split$HET, years = years, considered_indices = considered_indices, title = "Exo Het")
plot_stats(subgroup = incidence$exo_split$MSM, years = years, considered_indices = considered_indices, title = "Exo MSM")
plot_stats(subgroup = incidence$exo_split$IDU, years = years, considered_indices = considered_indices, title = "Exo IDU")
plot_stats(subgroup = incidence$exo_split$`unknown/other`, years = years, considered_indices = considered_indices, "Exo unknown/other")

plot_stats(subgroup = incidence$endo_combined_risk, years = years, considered_indices = considered_indices, title = "Endo (all risk groups)")
plot_stats(subgroup = incidence$exo_combined_risk, years = years, considered_indices = considered_indices, title = "Exo (all risk groups)")

plot_stats(subgroup = incidence$combined_ee$HET, years = years, considered_indices = considered_indices, title = "Endo+Exo Het")
plot_stats(subgroup = incidence$combined_ee$MSM, years = years, considered_indices = considered_indices, title = "Endo+Exo MSM")
plot_stats(subgroup = incidence$combined_ee$IDU, years = years, considered_indices = considered_indices, title = "Endo+Exo IDU")
plot_stats(subgroup = incidence$combined_ee$`unknown/other`, years = years, considered_indices = considered_indices, "Endo+Exo unknown/other")
dev.off()

library(ggplot2)
library(gridExtra)
library(RColorBrewer)

#function to put stuff into format for use with ggplot
format_df <- function(incidence, quantity, years = c(2003, 2022)){
  years_all = as.numeric(rownames(incidence$total$incidence))
  considered_indices <- which(years_all >= years[1] & years_all <= years[2])
  ci <- considered_indices #shorten name
  n_years <- length(incidence$endo_split$HET$incidence_means[ci])
  quantity_vector <- c(incidence$endo_split$HET[[quantity]][ci],
                       incidence$endo_split$MSM[[quantity]][ci],
                       incidence$endo_split$IDU[[quantity]][ci],
                       incidence$endo_split$`unknown/other`[[quantity]][ci],
                       incidence$exo_split$HET[[quantity]][ci],
                       incidence$exo_split$MSM[[quantity]][ci],
                       incidence$exo_split$IDU[[quantity]][ci],
                       incidence$exo_split$`unknown/other`[[quantity]][ci])
  year_vector <- rep(as.numeric(names(incidence$endo_split$HET$incidence_means))[ci], 8)
  group_vector_single <- c("endo_HET", "endo_MSM", "endo_IDU", "endo_UO", 
                           "exo_HET", "exo_MSM", "exo_IDU", "exo_UO")
  group_vector <- rep(group_vector_single, each = n_years)
  df <- data.frame(q = quantity_vector,
                   year = year_vector,
                   group = group_vector)
  return(df)
}

#function to put stuff into format for use with ggplot (including error limits)
format_df_error <- function(incidence, quantity, quantity_error, years = c(2003, 2022)){
  years_all = as.numeric(rownames(incidence$total$incidence))
  considered_indices <- which(years_all >= years[1] & years_all <= years[2])
  ci <- considered_indices #shorten name
  n_years <- length(incidence$endo_split$HET$incidence_means[ci])
  quantity_vector <- c(incidence$endo_split$HET[[quantity]][ci],
                       incidence$endo_split$MSM[[quantity]][ci],
                       incidence$endo_split$IDU[[quantity]][ci],
                       incidence$endo_split$`unknown/other`[[quantity]][ci],
                       incidence$exo_split$HET[[quantity]][ci],
                       incidence$exo_split$MSM[[quantity]][ci],
                       incidence$exo_split$IDU[[quantity]][ci],
                       incidence$exo_split$`unknown/other`[[quantity]][ci])
  quantity2.5_vector <- c(incidence$endo_split$HET[[quantity_error]][1,ci],
                          incidence$endo_split$MSM[[quantity_error]][1,ci],
                          incidence$endo_split$IDU[[quantity_error]][1,ci],
                          incidence$endo_split$`unknown/other`[[quantity_error]][1,ci],
                          incidence$exo_split$HET[[quantity_error]][1,ci],
                          incidence$exo_split$MSM[[quantity_error]][1,ci],
                          incidence$exo_split$IDU[[quantity_error]][1,ci],
                          incidence$exo_split$`unknown/other`[[quantity_error]][1,ci])
  quantity97.5_vector <- c(incidence$endo_split$HET[[quantity_error]][2,ci],
                           incidence$endo_split$MSM[[quantity_error]][2,ci],
                           incidence$endo_split$IDU[[quantity_error]][2,ci],
                           incidence$endo_split$`unknown/other`[[quantity_error]][2,ci],
                           incidence$exo_split$HET[[quantity_error]][2,ci],
                           incidence$exo_split$MSM[[quantity_error]][2,ci],
                           incidence$exo_split$IDU[[quantity_error]][2,ci],
                           incidence$exo_split$`unknown/other`[[quantity_error]][2,ci])
  year_vector <- rep(as.numeric(names(incidence$endo_split$HET$incidence_means))[ci], 8)
  group_vector_single <- c("endo_HET", "endo_MSM", "endo_IDU", "endo_UO", 
                           "exo_HET", "exo_MSM", "exo_IDU", "exo_UO")
  group_vector <- rep(group_vector_single, each = n_years)
  df <- data.frame(q = quantity_vector,
                   q2.5 = quantity2.5_vector,
                   q97.5 = quantity97.5_vector,
                   year = year_vector,
                   group = group_vector)
  return(df)
}

#make axis labels
xlabels <- 2003:2022
#xlabels[!(xlabels %in% c(2005, 2010, 2015, 2020))] <- ""
xlabels[seq(1, length(xlabels), by =2)] <- ""

ylabels <- 1000*(0:9)
ylabels[c(2,4,6,8,10)] <- ""

pdf(file = paste0("incidence_area_plots1_", date_tag, ".pdf"), width = 9, height = 5.0625)
ggplot(data = format_df(incidence, "incidence_means"), aes(x = year, y = q, fill = group)) + geom_area() +
  theme_bw() + 
  scale_x_continuous(limits = c(2003,2022), expand = c(0,0), 
                     breaks = 2003:2022, labels = xlabels) +
  scale_y_continuous(limits = c(0,525), expand = c(0,0)) +
  scale_fill_manual(values = brewer.pal(8, "Paired")[c(1,3,5,7,2,4,6,8)]) +
  labs(fill = "") +
  xlab("Year") +
  ylab("Estimated Incidence (number of new infections)")

ggplot(data = format_df(incidence, "infections_current_means"), aes(x = year, y = q, fill = group)) + geom_area() +
  theme_bw() +
  scale_x_continuous(limits = c(2003,2022), expand = c(0,0), 
                     breaks = 2003:2022, labels = xlabels) +
  scale_y_continuous(limits = c(0,9400), expand = c(0,0), 
                     breaks = 1000*(0:9), labels = ylabels) +
  scale_fill_manual(values = brewer.pal(8, "Paired")[c(1,3,5,7,2,4,6,8)]) +
  labs(fill = "") +
  xlab("Year") +
  ylab("Estimated Current Number of Infected Individuals")

ggplot(data = format_df(incidence, "undiagnosed_current_means"), aes(x = year, y = q, fill = group)) + geom_area() +
  theme_bw() +
  scale_x_continuous(limits = c(2003,2022), expand = c(0,0), 
                     breaks = 2003:2022, labels = xlabels) +
  scale_y_continuous(limits = c(0,550), expand = c(0,0)) +
  scale_fill_manual(values = brewer.pal(8, "Paired")[c(1,3,5,7,2,4,6,8)]) +
  labs(fill = "") +
  xlab("Year") +
  ylab("Estimated Number of Undiagnosed Individuals")
dev.off()

#function to make long format data frame of raw samples
format_df_raw <- function(incidence, quantity, years = c(2003, 2022)){
  years_all = as.numeric(rownames(incidence$total$incidence))
  considered_indices <- which(years_all >= years[1] & years_all <= years[2])
  ci <- considered_indices #shorten name
  n_years <- length(ci)
  #Year <- as.character(years_all[ci])
  quantity_vector <- rbind(incidence$endo_split$HET[[quantity]][ci,],
                           incidence$endo_split$MSM[[quantity]][ci,],
                           incidence$endo_split$IDU[[quantity]][ci,],
                           incidence$endo_split$`unknown/other`[[quantity]][ci,],
                           incidence$exo_split$HET[[quantity]][ci,],
                           incidence$exo_split$MSM[[quantity]][ci,],
                           incidence$exo_split$IDU[[quantity]][ci,],
                           incidence$exo_split$`unknown/other`[[quantity]][ci,])
  group_vector_single <- c("endo_HET", "endo_MSM", "endo_IDU", "endo_UO", 
                           "exo_HET", "exo_MSM", "exo_IDU", "exo_UO")
  group_vector <- rep(group_vector_single, each = n_years)
  year_vector <- rep(years_all[ci], 8)
  sample_df <- as.data.frame(quantity_vector)
  sample_df$Group <- as.factor(group_vector)
  sample_df$Year <- as.factor(year_vector)
  sample_df <- reshape2::melt(sample_df, id = c("Year", "Group"))
  sample_df <- sample_df[,-3] #remove column for sample number names
  return(sample_df)
}

diag_prop_raw_df <- format_df_raw(incidence, quantity = "diagnosed_fraction_current", years = c(2003, 2022))
incidence_raw_df <- format_df_raw(incidence, quantity = "incidence", years = c(2003, 2022))

#more general function to format any group into long format
format_df_raw_general <- function(data, years = c(2003, 2022)){
  years_all = as.numeric(rownames(data))
  considered_indices <- which(years_all >= years[1] & years_all <= years[2])
  ci <- considered_indices #shorten name
  n_years <- length(ci)
  #Year <- as.character(years_all[ci])
  quantity_vector <- data[ci,]
  #group_vector <- rep(group_vector_single, each = n_years)
  year_vector <- years_all[ci]
  sample_df <- as.data.frame(quantity_vector)
  #sample_df$Group <- as.factor(group_vector)
  sample_df$Year <- as.factor(year_vector)
  sample_df <- reshape2::melt(sample_df, id = c("Year"))
  sample_df <- sample_df[,-2] #remove column for sample number names
  return(sample_df)
}

diag_prop_raw_df_all_endo <- format_df_raw_general(data = incidence$endo_combined_risk$diagnosed_fraction_current, years = c(2003, 2022))
diag_prop_raw_df_all_exo <- format_df_raw_general(data = incidence$exo_combined_risk$diagnosed_fraction_current, years = c(2003, 2022))
diag_prop_raw_df_all_all <- format_df_raw_general(data = incidence$total$diagnosed_fraction_current, years = c(2003, 2022))
diag_prop_raw_df_both_HET <- format_df_raw_general(data = incidence$combined_ee$HET$diagnosed_fraction_current, years = c(2003, 2022))
diag_prop_raw_df_both_IDU <- format_df_raw_general(data = incidence$combined_ee$IDU$diagnosed_fraction_current, years = c(2003, 2022))
diag_prop_raw_df_both_MSM <- format_df_raw_general(data = incidence$combined_ee$MSM$diagnosed_fraction_current, years = c(2003, 2022))

ggplot(data = diag_prop_raw_df, 
       aes(x = Year, y = value, fill = Group, color = Group)) + geom_violin(adjust = 4, trim = TRUE) + 
  scale_fill_manual(values = brewer.pal(8, "Paired")[c(1,3,5,7,2,4,6,8)]) +
  scale_color_manual(values = brewer.pal(8, "Paired")[c(1,3,5,7,2,4,6,8)]) +
  theme_bw()

#names of groups by endoexo status and risk group
group_names <- levels(diag_prop_raw_df$Group)
#defined colors
colors <- brewer.pal(8, "Paired")[c(1,5,3,7,2,6,4,8)]
#more sparse x axis labels
xlabels_sparse <- rep("", length(2003:2022))
xlabels_sparse[c(3,8,13,18)] <- c(2005, 2010, 2015, 2020)

#plots for old version with proportions instead of percentages
p_prop_v <- list()
for(i in seq_along(group_names)){
  p_prop_v[[i]] <- ggplot(data = diag_prop_raw_df[diag_prop_raw_df$Group == group_names[i],], aes(x = Year, y = value)) +
    geom_violin(color = colors[i], fill = colors[i], adjust = 4) +
    #geom_errorbar(aes(ymin = q2.5, ymax = q97.5), color = colors[i]) + 
    geom_hline(yintercept = 0.95) +
    ylim(c(0.7, 1)) +
    xlab("Year") +
    ylab("Proportion Diagnosed") +
    scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
    theme_bw() +
    ggtitle(group_names[i])
  if(i >= 5) p_prop_v[[i]] <- p_prop_v[[i]] + theme(axis.title.y = element_blank())
  if(!(i %in% c(4,8))) p_prop_v[[i]] <- p_prop_v[[i]] + theme(axis.title.x = element_blank())
}

pdf(file = paste0("diagnosed_proportion_all_violin_", date_tag, ".pdf"), width = 7.5, height = 8)
grid.arrange(p_prop_v[[1]], p_prop_v[[5]], 
             p_prop_v[[3]], p_prop_v[[7]],
             p_prop_v[[2]], p_prop_v[[6]],
             p_prop_v[[4]], p_prop_v[[8]],
             widths = c(1,1), heights = c(1,1,1,1.08),
             nrow = 4)
dev.off()

#change values to percentages instead of proportions
#diag_prop_raw_df_pc <- diag_prop_raw_df
#diag_prop_raw_df_pc$value <- 100*diag_prop_raw_df_pc$value

#diag_prop_raw_df_all_endo_pc <- diag_prop_raw_df_all_endo
#diag_prop_raw_df_all_endo_pc$value <- 100*diag_prop_raw_df_all_endo_pc$value

#diag_prop_raw_df_all_exo_pc <- diag_prop_raw_df_all_exo
#diag_prop_raw_df_all_exo_pc$value <- 100*diag_prop_raw_df_all_exo_pc$value

#diag_prop_raw_df_all_all_pc <- diag_prop_raw_df_all_all
#diag_prop_raw_df_all_all_pc$value <- 100*diag_prop_raw_df_all_all_pc$value

#defined colors
colors <- brewer.pal(8, "Paired")[c(1,3,5,7,2,4,6,8)]

p_prop_v <- list()
for(i in seq_along(group_names)){
  p_prop_v[[i]] <- ggplot(data = diag_prop_raw_df[diag_prop_raw_df$Group == group_names[i],], aes(x = Year, y = 100*value)) +
    geom_violin(color = colors[i], fill = colors[i], adjust = 4) +
    geom_hline(yintercept = 95) +
    xlab("Year") +
    scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
    scale_y_continuous(expand = c(0,0), limits = c(70, 100), labels = c("70%", "80%", "90%", "100%")) +
    theme_bw() +
    ggtitle(group_names[i]) +
    theme(axis.title.y = element_blank())
  if(!(i %in% c(3,7))) p_prop_v[[i]] <- p_prop_v[[i]] + theme(axis.title.x = element_blank())
}

#all endos
p_prop_v[[9]] <- ggplot(data = diag_prop_raw_df_all_endo, aes(x = Year, y = 100*value)) +
  geom_violin(color = "#888888", fill = "#888888", adjust = 4) +
  geom_hline(yintercept = 95) +
  xlab("Year") +
  scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
  scale_y_continuous(expand = c(0,0), limits = c(70, 100), labels = c("70%", "80%", "90%", "100%")) +
  theme_bw() +
  ggtitle("endo_all") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

p_prop_v[[10]] <- ggplot(data = diag_prop_raw_df_all_exo, aes(x = Year, y = 100*value)) +
  geom_violin(color = "#000000", fill = "#000000", adjust = 4) +
  geom_hline(yintercept = 95) +
  xlab("Year") +
  scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
  scale_y_continuous(expand = c(0,0), limits = c(70, 100), labels = c("70%", "80%", "90%", "100%")) +
  theme_bw() +
  ggtitle("exo_all") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

pdf(file = paste0("diagnosed_proportion_all_violin_new_", date_tag, ".pdf"), width = 7.5, height = 8)
grid.arrange(p_prop_v[[9]], p_prop_v[[10]],
             p_prop_v[[1]], p_prop_v[[5]], 
             p_prop_v[[2]], p_prop_v[[6]],
             p_prop_v[[3]], p_prop_v[[7]],
             widths = c(1,1), heights = c(1,1,1,1.08),
             nrow = 4)
dev.off()

art_vs <- read.csv(file = "95_goals.csv")

UNAIDS_goals <- c("Diagnosed" = "#000000", "On ART" = "#000000", "Viral Suppression" = "#000000")

#everything together
p_prop_v[[11]] <- ggplot(data = diag_prop_raw_df_all_all, aes(x = Year, y = 100*value, col = "Diagnosed")) +
  geom_violin(adjust = 4) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = art_all, col = "On ART"), shape = 1) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = vs_all, col = "Viral Suppression"), shape = 4) +
  geom_hline(yintercept = 95) +
  xlab("Year") +
  labs(color = "Goal") +
  scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
  scale_y_continuous(expand = c(0,0), limits = c(70, 100), labels = c("70%", "80%", "90%", "100%")) +
  scale_color_manual(values = UNAIDS_goals)+
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

plot(p_prop_v[[11]]+theme(legend.position = c(0.8, 0.25)))
pdf(file = paste0("95-95-95_goals_hollow_", date_tag, ".pdf"), width = 7.5, height = 5)
plot(p_prop_v[[11]]+theme(legend.position = "none"))
dev.off()

#UNAIDS_goals <- c("Diagnosed" = "#648FFF", "On ART" = "#DC267F", "Viral Suppression" = "#000000")
#UNAIDS_goals <- c("Diagnosed" = "#80BFFF", "On ART" = "#BC165F", "Viral Suppression" = "#000000")
#UNAIDS_goals <- c("Diagnosed" = "#80BFFF", "On ART" = "red", "Viral Suppression" = "green4")
#UNAIDS_goals <- c("Diagnosed" = "#80E5FF", "On ART" = "#EC8E13", "Viral Suppression" = "#A607FF")
#UNAIDS_goals <- c("Diagnosed" = "#80E5FF", "On ART" = "#EC8E13", "Viral Suppression" = "#8B02D8")
#use colorblind-friendly colors
UNAIDS_goals <- c("Diagnosed" = "#67C9E2", "On ART" = "#EC8E13", "Viral Suppression" = "#8B02D8")

#blank data frame for legend
blank_df <- data.frame(group = c("Diagnosed", "On ART", "Viral Suppression"),
                       value = rep(2000, 3))

#data frame for means of years of interest
diag_prop_mean_df <- data.frame(Year = intervals, means = incidence$total$diagnosed_fraction_current_means[considered_indices])

#everything together (more color)
p_prop_v_all <- ggplot(data = diag_prop_raw_df_all_all, aes(x = as.integer(as.character(Year)), y = 100*value)) +
  geom_violin(aes(group = cut_width(as.integer(as.character(Year)),1)), col = "#80E5FF", fill = "#80E5FF") +
  geom_point(data = diag_prop_mean_df, aes(x = Year, y = 100*means), shape = 16,  col = "#67C9E2") +
  geom_line(data = diag_prop_mean_df, aes(x = Year, y = 100*means), col = "#67C9E2") +
  geom_point(data = art_vs, aes(x = year, y = art_all), shape = 16, col = "#EC8E13") +
  geom_line(data = art_vs, aes(x = year, y = art_all), col = "#EC8E13") +
  geom_point(data = art_vs, aes(x = year, y = vs_all), shape = 16, col = "#8B02D8") +
  geom_line(data = art_vs, aes(x = year, y = vs_all), col = "#8B02D8") +
  geom_hline(yintercept = 95) +
  geom_point(data = blank_df, aes(x = value, col = group)) +
  xlab("Year") +
  labs(color = "UNAIDS Goal") +
  scale_x_continuous(expand = c(0,0), limits = c(2002.5, 2022.5), breaks = 2003:2022, labels = xlabels_sparse, 
                     minor_breaks = 2003:2002) +
  scale_y_continuous(expand = c(0,0), limits = c(75, 100), labels = c("75%", "80%", "85%", "90%", "95%","100%")) +
  scale_color_manual(values = UNAIDS_goals)+
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(text = element_text(size = 14))
  #theme(legend.position = "none")
plot(p_prop_v_all+theme(legend.position = c(0.8, 0.25)))
pdf(file = paste0("95-95-95_goals_colors_lines_points_", date_tag, ".pdf"), width = 7.5, height = 5)
plot(p_prop_v_all+theme(legend.position = c(0.8, 0.25)))
dev.off()

p_prop_v[[12]] <- ggplot(data = diag_prop_raw_df_both_HET, aes(x = Year, y = 100*value)) +
  geom_violin(adjust = 4, col = colors[5]) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = art_het), col = colors[5], shape = 1) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = vs_het), col = colors[5], shape = 4) +
  geom_hline(yintercept = 95) +
  xlab("Year") +
  labs(color = "Goal") +
  scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
  scale_y_continuous(expand = c(0,0), limits = c(65, 100), labels = c(NA, "70%", "80%", "90%", "100%")) +
  scale_color_manual(values = UNAIDS_goals)+
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ggtitle("HET")
plot(p_prop_v[[12]])

p_prop_v[[13]] <- ggplot(data = diag_prop_raw_df_both_IDU, aes(x = Year, y = 100*value)) +
  geom_violin(adjust = 4, col = colors[6]) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = art_idu), col = colors[6], shape = 1) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = vs_idu), col = colors[6], shape = 4) +
  geom_hline(yintercept = 95) +
  xlab("Year") +
  labs(color = "Goal") +
  scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
  scale_y_continuous(expand = c(0,0), limits = c(65, 100), labels = c(NA,"70%", "80%", "90%", "100%")) +
  scale_color_manual(values = UNAIDS_goals)+
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ggtitle("IDU")
plot(p_prop_v[[13]])

p_prop_v[[14]] <- ggplot(data = diag_prop_raw_df_both_MSM, aes(x = Year, y = 100*value)) +
  geom_violin(adjust = 4, col = colors[7]) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = art_msm), col = colors[7], shape = 1) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = vs_msm), col = colors[7], shape = 4) +
  geom_hline(yintercept = 95) +
  xlab("Year") +
  labs(color = "Goal") +
  scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
  scale_y_continuous(expand = c(0,0), limits = c(65, 100), labels = c(NA, "70%", "80%", "90%", "100%")) +
  scale_color_manual(values = UNAIDS_goals)+
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  ggtitle("MSM")
plot(p_prop_v[[14]])

#all again, but with title and different limits
p_prop_v[[15]] <- ggplot(data = diag_prop_raw_df_all_all, aes(x = Year, y = 100*value)) +
  geom_violin(adjust = 4, col = "#000000") +
  geom_point(data = art_vs, aes(x = as.factor(year), y = art_all), col = "#000000", shape = 1) +
  geom_point(data = art_vs, aes(x = as.factor(year), y = vs_all), col = "#000000", shape = 4) +
  geom_hline(yintercept = 95) +
  xlab("Year") +
  labs(color = "Goal") +
  scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
  scale_y_continuous(expand = c(0,0), limits = c(65, 100), labels = c(NA, "70%", "80%", "90%", "100%")) +
  scale_color_manual(values = UNAIDS_goals)+
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ggtitle("All")
plot(p_prop_v[[15]])

pdf(file = paste0("95-95-95_goals_combined_endoexo_", date_tag, ".pdf"), width = 7.5, height = 9)
grid.arrange(p_prop_v[[15]],
             p_prop_v[[12]],
             p_prop_v[[13]],
             p_prop_v[[14]],
             heights = c(1,1,1,1.08),
             nrow = 4)
dev.off()

#set widths manually
widths <- rep(0.7, 8)
#widths[3] <- 1
p_inc_v <- list()
for(i in seq_along(group_names)){
  p_inc_v[[i]] <- ggplot(data = incidence_raw_df[incidence_raw_df$Group == group_names[i],], aes(x = Year, y = value)) +
    geom_violin(color = colors[i], fill = colors[i], adjust = 2, width = widths[i], scale = "width") +
    #geom_errorbar(aes(ymin = q2.5, ymax = q97.5), color = colors[i]) + 
    ylim(c(0, 230)) +
    xlab("Year") +
    ylab("Incidence") +
    scale_x_discrete(breaks = 2003:2022, labels = xlabels_sparse) +
    theme_bw() +
    ggtitle(group_names[i])
  if(i >= 5) p_inc_v[[i]] <- p_inc_v[[i]] + theme(axis.title.y = element_blank())
  if(!(i %in% c(4,8))) p_inc_v[[i]] <- p_inc_v[[i]] + theme(axis.title.x = element_blank())
}

pdf(file = paste0("incidence_all_violin_", date_tag, ".pdf"), width = 7.5, height = 8)
grid.arrange(p_inc_v[[1]], p_inc_v[[5]], 
             p_inc_v[[3]], p_inc_v[[7]],
             p_inc_v[[2]], p_inc_v[[6]],
             p_inc_v[[4]], p_inc_v[[8]],
             widths = c(1,1), heights = c(1,1,1,1.08),
             nrow = 4)
dev.off()

###Plot diagnosed proportion for all groups separately
diag_prop_df <- format_df_error(incidence = incidence, 
                                quantity = "diagnosed_fraction_current_means",
                                quantity_error = "diagnosed_fraction_current_95",
                                years = c(2003, 2022))
#names of groups by endoexo status and risk group
group_names <- unique(diag_prop_df$group)
#defined colors
colors <- brewer.pal(8, "Paired")[c(1,5,3,7,2,6,4,8)]

p_prop <- list()
for(i in seq_along(group_names)){
  p_prop[[i]] <- ggplot(data = diag_prop_df[diag_prop_df$group == group_names[i],], aes(x = year, y = q)) +
    geom_point(color = colors[i]) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), color = colors[i]) + 
    geom_hline(yintercept = 0.95) +
    ylim(c(0.75, 1)) +
    xlab("Year") +
    ylab("Proportion Diagnosed") +
    theme_bw() +
    ggtitle(group_names[i])
  if(i >= 5) p_prop[[i]] <- p_prop[[i]] + theme(axis.title.y = element_blank())
  if(!(i %in% c(4,8))) p_prop[[i]] <- p_prop[[i]] + theme(axis.title.x = element_blank())
}

pdf(file = paste0("diagnosed_proportion_all_", date_tag, ".pdf"), width = 7.5, height = 8)
grid.arrange(p_prop[[1]], p_prop[[5]], 
             p_prop[[3]], p_prop[[7]],
             p_prop[[2]], p_prop[[6]],
             p_prop[[4]], p_prop[[8]],
             widths = c(1,1), heights = c(1,1,1,1.08),
             nrow = 4)
dev.off()

###plot average infection ages for patients diagnosed in different years
yearly.infection.ages <- readRDS(file = paste0("mean.infection.ages.by.year.", date_tag, ".rds"))
plots <- list()
for(i in seq_along(yearly.infection.ages$mean_df)){
  plots[[i]] <- ggplot(data = yearly.infection.ages$mean_df[[i]], aes(x = Year, y = Mean)) +
    geom_point() +
    theme_bw() +
    ylim(0,5.5) +
    xlab("Diagnosis Year") +
    ylab("Mean Infection Age at Diagnosis (Years)") +
    geom_errorbar(aes(ymin = Low2.5, ymax = High97.5)) +
    ggtitle(names(yearly.infection.ages$mean_df)[i])
}

pdf(file = paste0("mean_infection_age_at_diag_", date_tag, ".pdf"), width = 9, height = 6)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2)
dev.off()


#make plots for number of individuals
p_nInd <- list()
for(i in seq_along(yearly.infection.ages$mean_df)){
  p_nInd[[i]] <- ggplot(data = yearly.infection.ages$mean_df[[i]], aes(x = Year, y = nInds)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    ylim(0,90) +
    xlab("Diagnosis Year") +
    ylab("Number of Diagnoses") +
    ggtitle(names(yearly.infection.ages$mean_df)[i])
}

pdf(file = paste0("diagnoses_over_time_", date_tag, ".pdf"), width = 9, height = 6)
grid.arrange(p_nInd[[1]], p_nInd[[2]], p_nInd[[3]], p_nInd[[4]], nrow = 2)
dev.off()

colors <- brewer.pal(8, "Paired")[c(1,5,3)]

#plots all in one without the 'other' category
pts <- list()
for(i in 1:3){
  pts[[i]] <- ggplot(data = yearly.infection.ages$mean_df[[i]], aes(x = Year, y = Mean)) +
    geom_point(color = colors[i]) +
    theme_bw() +
    xlab("Diagnosis Year") + 
    scale_y_continuous(limits = c(0,4.5), expand = c(0,0)) +
    geom_errorbar(aes(ymin = Low2.5, ymax = High97.5), color = colors[i])
  if(i == 1){
    #labels.raw <- ggplot_build(pts[[i]])$layout$panel_params[[1]]$y$breaks
    #labels.padded <- stringr::str_pad(labels.raw[!is.na(labels.raw)], width = 4, pad = " ")
    pts[[i]] <- pts[[i]] + ylab("TI at Diagnosis (Years)") +
      theme(axis.title.y = element_text(margin = margin(r = 5)))
      #scale_y_continuous(labels = labels.padded)
  } 
  else pts[[i]] <- pts[[i]] + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
for(i in 1:3){
  pts[[i+3]] <- ggplot(data = yearly.infection.ages$mean_df[[i]], aes(x = Year, y = nInds)) +
    geom_bar(stat = "identity", color = colors[i], fill = colors[i]) +
    theme_bw() +
    #ylim(0,90) +
    ylab("") +
    scale_x_continuous(limits = c(2002.5,2022.5), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
    theme(axis.title.x = element_blank()) +
    ggtitle(names(yearly.infection.ages$mean_df)[i])
  if(i == 1) pts[[i+3]] <- pts[[i+3]] + ylab("Number of Diagnoses") + theme(axis.title.y = element_text(margin = margin(r = -5)))
  else pts[[i+3]] <- pts[[i+3]] + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

pdf(file = paste0("mean_infection_ages_combined_", date_tag, ".pdf"), width = 9, height = 6)
grid.arrange(pts[[4]], pts[[6]], pts[[5]],
             pts[[1]], pts[[3]], pts[[2]],
             widths = c(1.1,1,1),
             heights = c(1,1),
             nrow = 2)
dev.off()

#print stat values used in paper

#percent of exogenous cases in 2009
print(round(100*incidence$exo_combined_risk$incidence_prop_means["2009"], digits = 1))
print(round(100*incidence$exo_combined_risk$incidence_prop_95[,"2009"], digits = 1))

#percent of exo-HET cases in 2009
print(round(100*incidence$exo_split$HET$incidence_prop_means["2009"], digits = 1))
print(round(100*incidence$exo_split$HET$incidence_prop_95[,"2009"], digits = 1))

#total incidence between 2009 and 2014
print(round(incidence$total$incidence_means[as.character(2009:2014)]))
print(round(incidence$total$incidence_95[,as.character(2009:2014)]))

#range of mean percent of HET incidence (endo and exo) between 2003 and 2021
print(round(100*min(incidence$combined_ee$HET$incidence_prop_means[as.character(2003:2021)]), digits = 1))
print(round(100*max(incidence$combined_ee$HET$incidence_prop_means[as.character(2003:2021)]), digits = 1))

#incidence percent of exogenous MSM in 2021
print(round(100*incidence$exo_split$MSM$incidence_prop_means["2021"], digits = 1))
print(round(100*incidence$exo_split$MSM$incidence_prop_95[,"2021"], digits = 1))

#percent of endogenous incidence in 2022
print(round(100*incidence$endo_combined_risk$incidence_prop_means["2022"], digits = 1))
print(round(100*incidence$endo_combined_risk$incidence_prop_95[,"2022"], digits = 1))

#percent of endogenous incidence in 2003
print(round(100*incidence$endo_combined_risk$incidence_prop_means["2003"], digits = 1))
print(round(100*incidence$endo_combined_risk$incidence_prop_95[,"2003"], digits = 1))

#percent of endogenous incidence in 2022
print(round(100*incidence$exo_combined_risk$incidence_prop_means["2022"], digits = 1))
print(round(100*incidence$exo_combined_risk$incidence_prop_95[,"2022"], digits = 1))

#number of exo HET cases in 2018
print(round(incidence$exo_split$HET$incidence_means["2018"]))
print(round(incidence$exo_split$HET$incidence_95[,"2018"]))

#number of exo HET cases in 2021
print(round(incidence$exo_split$HET$incidence_means["2021"]))
print(round(incidence$exo_split$HET$incidence_95[,"2021"]))

#number of exo HET cases in 2022
print(round(incidence$exo_split$HET$incidence_means["2022"]))
print(round(incidence$exo_split$HET$incidence_95[,"2022"]))

#existing endogenous cases in 2003
print(round(incidence$endo_combined_risk$infections_current_means["2003"]))
print(round(incidence$endo_combined_risk$infections_current_95[,"2003"]))
#existing endogenous cases in 2022
print(round(incidence$endo_combined_risk$infections_current_means["2022"]))
print(round(incidence$endo_combined_risk$infections_current_95[,"2022"]))
#percent increase 2003 to 2022
print(round(100*incidence$endo_combined_risk$infections_current_means["2022"]/incidence$endo_combined_risk$infections_current_means["2003"]-100, digits = 1))

#existing exo cases in 2003
print(round(incidence$exo_combined_risk$infections_current_means["2003"]))
print(round(incidence$exo_combined_risk$infections_current_95[,"2003"]))
#existing exo cases in 2022
print(round(incidence$exo_combined_risk$infections_current_means["2022"]))
print(round(incidence$exo_combined_risk$infections_current_95[,"2022"]))
#percent increase 2003 to 2022
print(round(100*incidence$exo_combined_risk$infections_current_means["2022"]/incidence$exo_combined_risk$infections_current_means["2003"]-100, digits = 1))

#percent of existing cases that are exo HET in 2022
print(round(100*incidence$exo_split$HET$infections_current_prop_means["2022"], digits = 1))
print(round(100*incidence$exo_split$HET$infections_current_prop_95[,"2022"], digits = 1))

#IDU+UO current infection counts and proportions of current infections
#counts
IDU_UO_infections_current <- incidence$combined_ee$IDU$infections_current + incidence$combined_ee$`unknown/other`$infections_current
IDU_UO_infections_current_means <- rowMeans(IDU_UO_infections_current)
IDU_UO_infections_current_95 <- apply(IDU_UO_infections_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
#proportions
IDU_UO_infections_current_prop <- IDU_UO_infections_current/incidence$total$infections_current
IDU_UO_infections_current_prop_means <- rowMeans(IDU_UO_infections_current_prop)
IDU_UO_infections_current_prop_95 <- apply(IDU_UO_infections_current_prop, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
#number of existing cases that are IDU and unknown/other in 2022
print(round(IDU_UO_infections_current_means["2022"]))
print(round(IDU_UO_infections_current_95[,"2022"]))
#percent of existing cases that are IDU and unknown/other in 2022
print(round(100*IDU_UO_infections_current_prop_means["2022"], digits = 1))
print(round(100*IDU_UO_infections_current_prop_95[,"2022"], digits = 1))

#number of existing cases that are exo and unknown/other in 2022
print(round(incidence$exo_split$`unknown/other`$infections_current_means["2022"]))
print(round(incidence$exo_split$`unknown/other`$infections_current_95[,"2022"]))
#percent of existing cases that are exo and unknown/other in 2022
print(round(100*incidence$exo_split$`unknown/other`$infections_current_prop_means["2022"], digits = 1))
print(round(100*incidence$exo_split$`unknown/other`$infections_current_prop_95[,"2022"], digits = 1))

#number of undiagnosed endo HET cases in 2003
print(round(100*incidence$endo_split$HET$undiagnosed_current_prop_means["2003"], digits = 1))
print(round(100*incidence$endo_split$HET$undiagnosed_current_prop_95[,"2003"], digits = 1))
#number of undiagnosed endo HET cases in 2022
print(round(100*incidence$endo_split$HET$undiagnosed_current_prop_means["2022"], digits = 1))
print(round(100*incidence$endo_split$HET$undiagnosed_current_prop_95[,"2022"], digits = 1))

#number of undiagnosed endo MSM cases in 2003
print(round(100*incidence$endo_split$MSM$undiagnosed_current_prop_means["2003"], digits = 1))
print(round(100*incidence$endo_split$MSM$undiagnosed_current_prop_95[,"2003"], digits = 1))
#number of undiagnosed endo MSM cases in 2022
print(round(100*incidence$endo_split$MSM$undiagnosed_current_prop_means["2022"], digits = 1))
print(round(100*incidence$endo_split$MSM$undiagnosed_current_prop_95[,"2022"], digits = 1))

#exo HET and UO undiagnosed cases
exo_HET_UO_undiagnosed_current <- incidence$exo_split$HET$undiagnosed_current + incidence$exo_split$`unknown/other`$undiagnosed_current
exo_HET_UO_undiagnosed_current_means <- rowMeans(exo_HET_UO_undiagnosed_current)
exo_HET_UO_undiagnosed_current_95 <- apply(exo_HET_UO_undiagnosed_current, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
#number of undiagnosed exo HET and exo unknown/other in 2021
print(round(exo_HET_UO_undiagnosed_current_means["2021"], digits = 1))
print(exo_HET_UO_undiagnosed_current_95[,"2021"])
#number of undiagnosed exo HET and exo unknown/other in 2022
print(round(exo_HET_UO_undiagnosed_current_means["2022"], digits = 1))
print(exo_HET_UO_undiagnosed_current_95[,"2022"])
#percent increase
print(round(100*exo_HET_UO_undiagnosed_current_means["2022"]/exo_HET_UO_undiagnosed_current_means["2021"]-100, digits = 1))

#range of mean TI for IDU group
print(round(min(yearly.infection.ages$mean_df$IDU$Mean, na.rm = TRUE), digits = 2))
print(round(max(yearly.infection.ages$mean_df$IDU$Mean, na.rm = TRUE), digits = 2))

#range of mean TI for MSM group
print(round(min(yearly.infection.ages$mean_df$MSM$Mean, na.rm = TRUE), digits = 2))
print(round(max(yearly.infection.ages$mean_df$MSM$Mean, na.rm = TRUE), digits = 2))    

#percent diagnosed for endogenous infections 2022
print(round(100*incidence$endo_combined_risk$diagnosed_fraction_current_means["2022"], digits = 1))
print(round(100*incidence$endo_combined_risk$diagnosed_fraction_current_95[,"2022"], digits = 1))
#percent diagnosed for exo infections 2022
print(round(100*incidence$exo_combined_risk$diagnosed_fraction_current_means["2022"], digits = 1))
print(round(100*incidence$exo_combined_risk$diagnosed_fraction_current_95[,"2022"], digits = 1))

#find means of TI for known endo HET patients in 2015 and 2022
endo_TI_means <- readRDS(file = paste0("endo_TI_means.", prior, ".", date_tag, ".rds"))
endo_HET_2015_means <- with(endo_TI_means, mean[year == 2015 & route == "HET"])
endo_HET_2022_means <- with(endo_TI_means, mean[year == 2022 & route == "HET"])
#Mann-Whitney-Wilcoxon test
print(wilcox.test(endo_HET_2015_means, endo_HET_2022_means))

#print 75th percentile TI for known endo patients diagnosed
TI_2015_plus <- readRDS("TI_2015_plus.rds")
print(round(sapply(TI_2015_plus, FUN = quantile, probs = c(0.5, 0.75, 0.95)), digits = 1))
