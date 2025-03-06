
######################################################################
# B-Cell and Survival Analysis for ALL RD Study
#
# This script loads and cleans two datasets – one containing long-format
# B-cell data and another (Paul’s Dataset) with clinical outcomes. It then 
# merges the datasets, creates new time-to-event variables, defines the final 
# analytic sample, performs descriptive statistics, survival analysis, and 
# multivariable logistic regression analyses.
#
# Author: Paul George
# Date: 6 July 2024
######################################################################

###############################
# 1. Setup Environment
###############################

# Load required libraries
library(readxl)       # Read Excel files
library(tidyverse)    # Data manipulation and visualization
library(lubridate)    # Date manipulation
library(MASS)         # Statistical functions

# Set working directory (adjust the path as needed)

# Clear workspace
rm(list = ls())

###############################
# 2. Load and Clean Bcell_long Data
###############################

# Load the long-format B-cell dataset
load('Bcell_long.rda')
ls(Bcell_long)

# For patient with corporate_id 'C880828', update cytogenetics status to "Yes"
Bcell_long[Bcell_long$corporate_id == 'C880828', 'ABL-class or Ph'] <- 'Yes'
Bcell_long[Bcell_long$corporate_id == 'C880828', 'HR_cytogenetics'] <- 'Yes'

###############################
# 3. Load Paul’s Dataset and Select Variables
###############################

# Read Paul’s Dataset v5 from Excel
Pauls_Dataset_v5 <- read_xlsx('temp/Pauls_Dataset v5_1106023.xlsx')

# Select only the relevant columns: study_id, date of relapse/refractory, last contact, and date of death (dod)
Pauls_Dataset_v5_short <- Pauls_Dataset_v5 %>% 
  dplyr::select(study_id, `date_relapse-refractory`, last_contact, dod)

# Check the first few rows
head(Pauls_Dataset_v5_short)

###############################
# 4. Convert Date Columns to Date Format
###############################

# Define a function to clean and convert Excel date columns to Date objects
convert_to_date <- function(date_column) {
  # Remove non-numeric characters and convert to numeric
  numeric_dates <- as.numeric(gsub("[^0-9]", "", date_column))
  # Replace zero values with NA
  numeric_dates[numeric_dates == 0] <- NA
  # Convert from Excel date format (origin = "1899-12-30")
  dates <- as.Date(numeric_dates, origin = "1899-12-30")
  return(dates)
}

# Apply the conversion function to the three date columns using mutate and across
Pauls_Dataset_v5_short <- Pauls_Dataset_v5_short %>%
  mutate(
    `date_relapse-refractory` = convert_to_date(`date_relapse-refractory`),
    last_contact = convert_to_date(last_contact),
    dod = convert_to_date(dod)
  )

# Confirm changes with a peek at the data
head(Pauls_Dataset_v5_short)

###############################
# 5. Merge Datasets
###############################

# Left join Pauls_Dataset_v5_short to Bcell_long by common column (assumed to be study_id)
Bcell_long <- left_join(Bcell_long, Pauls_Dataset_v5_short)

# Inspect a subset of variables to ensure correct merging
checkdf <- Bcell_long %>% 
  dplyr::select(study_id, dod, relapse_refractory, `date_relapse-refractory`, last_contact, poor_outcome)

###############################
# 6. Select Variables of Interest
###############################

Bcell.df <- Bcell_long %>% 
  dplyr::select(
    study_id, corporate_id, age_at_diagnosis, age_at_diagnosis_dichot, 
    Gender, Race, Ethnicity, `Date of Birth`,
    `Diagnosis Date`, diagnosis_year, HR_cytogenetics, wbc_value, leukemia_type, 
    combined_protocol, induction_risk,
    refractory_disease, relapse_disease, patient_died, poor_outcome, 
    induction_mrd_result, consol_mrd_yn, consol_mrd_result, 
    PO_rr_alive, PO_rr_dead, PO_dead_no_disease, Deceased, 
    relapse_refractory, dod, `date_relapse-refractory`, last_contact
  )

###############################
# 7. Create Time-to-Event Variables
###############################

# Create a new variable: date_of_event_or_last_contact is the earliest non-missing event among dod, relapse/refractory date, and last contact
Bcell.df <- Bcell.df %>%
  mutate(date_of_event_or_last_contact = pmin(dod, `date_relapse-refractory`, last_contact, na.rm = TRUE))

# Calculate time to event or last contact (in years) using Diagnosis Date as baseline
Bcell.df <- Bcell.df %>%
  mutate(time_to_event_or_last_contact = as.numeric(date_of_event_or_last_contact - `Diagnosis Date`)/365.25)

# Similarly, create date_of_last_contact as the maximum of the three dates
Bcell.df <- Bcell.df %>%
  mutate(date_of_last_contact = pmax(dod, `date_relapse-refractory`, last_contact, na.rm = TRUE))

# Calculate time to last contact (in years)
Bcell.df <- Bcell.df %>%
  mutate(time_to_last_contact = as.numeric(date_of_last_contact - `Diagnosis Date`)/365.25)

# Check key variables
checkdf <- Bcell.df %>% 
  dplyr::select(study_id, `Diagnosis Date`, dod, relapse_refractory, `date_relapse-refractory`, 
                last_contact, poor_outcome, date_of_event_or_last_contact, time_to_event_or_last_contact, 
                date_of_last_contact, time_to_last_contact)

###############################
# 8. Define Final Analytic Sample
###############################

# Exclude patients with HR cytogenetics (keep those with 'No')
not_HR.df <- Bcell.df %>% filter(HR_cytogenetics == 'No')

# Remove patients with WBC values > 50 (check counts by age group first)
not_HR.df %>%
  group_by(age_at_diagnosis_dichot) %>%
  summarize(
    wbc_gt_50 = sum(wbc_value > 50),
    wbc_lt_50 = sum(wbc_value < 50)
  )

not_HR.df <- not_HR.df %>% filter(wbc_value < 50)

# Additional removals based on clinical criteria:
# - For patients <10 years old, keep only those with Standard Risk induction
final_df <- not_HR.df %>% 
  filter((age_at_diagnosis_dichot == '<10' & induction_risk == "Standard Risk") | 
           age_at_diagnosis_dichot == '>=10')

# Remove a specific patient with corporate_id 'C2637677' (e.g., steroid pre-treatment exclusion)
final_df <- final_df %>% filter(corporate_id != 'C2637677')

###############################
# 9. Initial Descriptive Analyses
###############################

# Summary of time to event or last contact
summary(final_df$time_to_event_or_last_contact)

# Descriptive tables for age at diagnosis
tab <- table(final_df$age_at_diagnosis_dichot)
cbind(Counts = tab, Percentage = prop.table(tab) * 100)

# Gender distribution
tab <- table(final_df$Gender)
cbind(Counts = tab, Percentage = prop.table(tab) * 100)

# Race distribution
tab <- table(final_df$Race)
cbind(Counts = tab, Percentage = prop.table(tab) * 100)

# Ethnicity distribution
tab <- table(final_df$Ethnicity)
cbind(Counts = tab, Percentage = prop.table(tab) * 100)

# Summary of diagnosis year
summary(final_df$diagnosis_year)

# Crosstabs for induction MRD results and other outcomes
table(final_df$age_at_diagnosis_dichot, final_df$induction_mrd_result, final_df$PO_rr_alive, useNA = 'always')
table(final_df$induction_mrd_result, final_df$consol_mrd_result, useNA = 'always')

# Subset: Patients with positive induction MRD
induction_MRD_pos_df <- final_df %>% filter(induction_mrd_result == 'Positive')

# For older patients (>=10 years)
olderkids_df <- final_df %>% filter(age_at_diagnosis_dichot == '>=10')
table(olderkids_df$poor_outcome, olderkids_df$induction_mrd_result)
table(induction_MRD_pos_df$consol_mrd_result, useNA = 'always')
table(induction_MRD_pos_df$age_at_diagnosis_dichot)
table(induction_MRD_pos_df$Deceased, induction_MRD_pos_df$age_at_diagnosis_dichot)
table(induction_MRD_pos_df$Deceased, induction_MRD_pos_df$consol_mrd_result)

table(final_df$poor_outcome)

# Age summary statistics by poor outcome
final_df %>%
  group_by(poor_outcome, age_at_diagnosis_dichot) %>%
  summarize(
    count = n(),
    mean = mean(age_at_diagnosis),
    sd = sd(age_at_diagnosis),
    min = min(age_at_diagnosis),
    median = median(age_at_diagnosis),
    max = max(age_at_diagnosis)
  )

# Fischer exact test for age groups vs. poor outcome
table_data <- table(final_df$poor_outcome, final_df$age_at_diagnosis_dichot)
fisher.test(table_data)

# Gender summaries by poor outcome
final_df$Gender <- as.factor(final_df$Gender)
final_df %>%
  group_by(poor_outcome, Gender) %>%
  summarize(n = n()) %>%
  group_by(poor_outcome) %>%
  mutate(percent = n / sum(n) * 100)

# Chi-squared test for gender and poor outcome
chisq.test(final_df$Gender, final_df$poor_outcome)

# Ethnicity summaries by poor outcome
final_df$Ethnicity <- as.factor(final_df$Ethnicity)
final_df %>%
  group_by(poor_outcome, Ethnicity) %>%
  summarize(n = n()) %>%
  group_by(poor_outcome) %>%
  mutate(percent = n / sum(n) * 100)
chisq.test(final_df$Ethnicity, final_df$poor_outcome)

# Race summaries by poor outcome
final_df %>%
  group_by(poor_outcome, Race) %>%
  summarize(n = n()) %>%
  group_by(poor_outcome) %>%
  mutate(percent = n / sum(n) * 100)

# Fisher test for race and poor outcome
final_df %>%
  group_by(poor_outcome, Race) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = Race, values_from = n, values_fill = 0) %>%
  dplyr::select(-poor_outcome) %>%
  fisher.test()

# Diagnosis year summary by poor outcome
final_df %>%
  group_by(poor_outcome) %>%
  summarize(
    count = n(),
    mean = mean(diagnosis_year),
    sd = sd(diagnosis_year),
    min = min(diagnosis_year),
    median = median(diagnosis_year),
    max = max(diagnosis_year),
    IQR = IQR(diagnosis_year)
  )
t.test(diagnosis_year ~ poor_outcome, data = final_df)

# Time to last contact summary by age group
final_df %>%
  group_by(age_at_diagnosis_dichot) %>%
  summarize(
    count = n(),
    mean = mean(time_to_last_contact),
    sd = sd(time_to_last_contact),
    min = min(time_to_last_contact),
    median = median(time_to_last_contact),
    max = max(time_to_last_contact),
    IQR = IQR(time_to_last_contact)
  )
t.test(time_to_last_contact ~ age_at_diagnosis_dichot, data = final_df)

###############################
# 10. Survival Analysis: Kaplan-Meier Curves
###############################

library(survival)
library(survminer)
library(RColorBrewer)

# Rename age_at_diagnosis_dichot factor levels for clarity
final_df$age_at_diagnosis_dichot <- factor(final_df$age_at_diagnosis_dichot,
                                           levels = c("<10", ">=10"),
                                           labels = c("Age <10 years", "Age >=10 years"))

# Create a survival object using time to event or last contact and poor outcome as event indicator
surv_object <- Surv(time = final_df$time_to_event_or_last_contact, event = final_df$poor_outcome)

# Fit Kaplan-Meier survival curves stratified by age group
fit <- survfit(surv_object ~ age_at_diagnosis_dichot, data = final_df)

# Define a colorblind-friendly palette
color_blind_friendly_palette <- brewer.pal(n = 2, name = "Dark2")

# Plot the survival curves with risk table and confidence intervals
g <- ggsurvplot(fit,
                data = final_df,
                xlab = "Time to event or last contact (years)",
                ylab = "Event-free survival probability",
                ylim = c(0.6, 1), 
                title = "Supplemental Figure S1. Kaplan-Meier Survival Curve",
                subtitle = "by age at diagnosis", 
                legend.labs = c("Age < 10 years", "Age >= 10 years"), 
                surv.median.line = "hv",
                risk.table = TRUE,
                conf.int = TRUE,
                pval = TRUE,
                pval.coord = c(1, 0.65),
                palette = color_blind_friendly_palette,
                linetype = c("solid", "longdash"))
print(g)

# Save the Kaplan-Meier plot as a high-resolution PNG
ggsave("kaplan_meier_plot.png", plot = g$plot, 
       width = 8, height = 5, dpi = 600)

###############################
# 11. Multivariable Cox Proportional Hazards Model
###############################

# Ensure categorical variables are properly formatted and set reference levels
final_df$age_at_diagnosis_dichot <- as.factor(final_df$age_at_diagnosis_dichot)
final_df$Race <- as.factor(final_df$Race)
final_df$Race <- relevel(final_df$Race, ref = "White")
final_df$Ethnicity <- as.factor(final_df$Ethnicity)
final_df$Ethnicity <- relevel(final_df$Ethnicity, ref = "Non-Hispanic/Non-Latino")

# Ensure diagnosis_year is numeric
final_df$diagnosis_year <- as.numeric(as.character(final_df$diagnosis_year))

# Create survival object
surv_obj <- Surv(time = final_df$time_to_event_or_last_contact, event = final_df$poor_outcome)

# Fit the Cox model with multiple covariates
cox_model <- coxph(surv_obj ~ age_at_diagnosis_dichot + Race + Ethnicity + Gender + diagnosis_year, data = final_df)
summary(cox_model)

###############################
# 12. Additional Outcome Analyses
###############################

# Summary statistics of poor outcome by age group
final_df %>%
  group_by(age_at_diagnosis_dichot, poor_outcome) %>%
  summarize(n = n()) %>%
  group_by(age_at_diagnosis_dichot) %>%
  mutate(percent = n / sum(n) * 100)

# Chi-squared test for poor outcome and age group
chisq.test(final_df$poor_outcome, final_df$age_at_diagnosis_dichot)

# Summary statistics of PO_rr_alive by age group
final_df %>%
  group_by(age_at_diagnosis_dichot, PO_rr_alive) %>%
  summarize(n = n()) %>%
  group_by(age_at_diagnosis_dichot) %>%
  mutate(percent = n / sum(n) * 100)
fisher.test(final_df$PO_rr_alive, final_df$age_at_diagnosis_dichot)

# Summary statistics of PO_rr_dead by age group
final_df %>%
  group_by(age_at_diagnosis_dichot, PO_rr_dead) %>%
  summarize(n = n()) %>%
  group_by(age_at_diagnosis_dichot) %>%
  mutate(percent = n / sum(n) * 100)
fisher.test(final_df$PO_rr_dead, final_df$age_at_diagnosis_dichot)

# Summary statistics of PO_dead_no_disease by age group
final_df %>%
  group_by(age_at_diagnosis_dichot, PO_dead_no_disease) %>%
  summarize(n = n()) %>%
  group_by(age_at_diagnosis_dichot) %>%
  mutate(percent = n / sum(n) * 100)
fisher.test(final_df$PO_dead_no_disease, final_df$age_at_diagnosis_dichot)

# Diagnosis year summary by poor outcome
final_df %>%
  group_by(poor_outcome) %>%
  summarize(
    count = n(),
    mean = mean(diagnosis_year),
    sd = sd(diagnosis_year),
    min = min(diagnosis_year),
    median = median(diagnosis_year),
    max = max(diagnosis_year),
    IQR = IQR(diagnosis_year)
  )
t.test(diagnosis_year ~ poor_outcome, data = final_df)

# Time to last contact summary by age group
final_df %>%
  group_by(age_at_diagnosis_dichot) %>%
  summarize(
    count = n(),
    mean = mean(time_to_last_contact),
    sd = sd(time_to_last_contact),
    min = min(time_to_last_contact),
    median = median(time_to_last_contact),
    max = max(time_to_last_contact),
    IQR = IQR(time_to_last_contact)
  )
t.test(time_to_last_contact ~ age_at_diagnosis_dichot, data = final_df)

###############################
# 13. Multivariable Logistic Regression Analyses
###############################

# Ensure variables are correctly formatted
final_df$Race <- as.factor(final_df$Race)
final_df$Race <- relevel(final_df$Race, ref = "White")
final_df$Ethnicity <- as.factor(final_df$Ethnicity)
final_df$Ethnicity <- relevel(final_df$Ethnicity, ref = 'Non-Hispanic/Non-Latino')
Bcell.df$Race <- as.factor(Bcell.df$Race)
Bcell.df$Race <- relevel(Bcell.df$Race, ref = "White")
Bcell.df$Ethnicity <- as.factor(Bcell.df$Ethnicity)
Bcell.df$Ethnicity <- relevel(Bcell.df$Ethnicity, ref = 'Non-Hispanic/Non-Latino')
final_df$age_at_diagnosis_dichot <- as.factor(final_df$age_at_diagnosis_dichot)
Bcell.df$age_at_diagnosis_dichot <- as.factor(Bcell.df$age_at_diagnosis_dichot)

# Model 1: Using Bcell.df with additional clinical variables
model1_age_dichot <- glm(poor_outcome ~ age_at_diagnosis_dichot + Gender + Race + Ethnicity + diagnosis_year + HR_cytogenetics + wbc_value, 
                         data = Bcell.df, family = binomial(link = "logit"))
mod1_sum <- summary(model1_age_dichot)
mod1_pvalues <- mod1_sum$coefficients[, 4]
mod1_OR <- exp(mod1_sum$coefficients[, 1])
mod1_pvalues_df <- data.frame(variable = names(mod1_pvalues), pvalue = mod1_pvalues)
mod1_OR_df <- data.frame(variable = names(mod1_OR), estimate = mod1_OR)
mod1_results <- merge(mod1_pvalues_df, mod1_OR_df, by = "variable")
mod1_results[, 2] <- round(mod1_results[, 2], 2)
mod1_results$OR <- formatC(mod1_results$estimate, format = "f", digits = 2)
mod1_results <- mod1_results %>% dplyr::select(-estimate)
head(mod1_results)
summary(model1_age_dichot)
exp(coef(model1_age_dichot))

# Model 2: Using final_df with fewer variables
model2_age_dichot <- glm(poor_outcome ~ Gender + Race + Ethnicity + diagnosis_year + age_at_diagnosis_dichot, 
                         data = final_df, family = binomial(link = "logit"))
mod2_sum <- summary(model2_age_dichot)
mod2_pvalues <- mod2_sum$coefficients[, 4]
mod2_OR <- exp(mod2_sum$coefficients[, 1])
mod2_pvalues_df <- data.frame(variable = names(mod2_pvalues), pvalue = mod2_pvalues)
mod2_OR_df <- data.frame(variable = names(mod2_OR), estimate = mod2_OR)
mod2_results <- merge(mod2_pvalues_df, mod2_OR_df, by = "variable")
mod2_results[, 2] <- round(mod2_results[, 2], 2)
mod2_results$OR <- formatC(mod2_results$estimate, format = "f", digits = 2)
mod2_results <- mod2_results %>% dplyr::select(-estimate)
mod2_results
summary(model2_age_dichot)

# Sensitivity Analysis: Only patients diagnosed prior to 2019
early_group <- final_df %>% filter(diagnosis_year < 2019)
model2_early <- glm(poor_outcome ~ Gender + Race + Ethnicity + diagnosis_year + age_at_diagnosis_dichot, 
                    data = early_group, family = binomial(link = "logit"))
summary(model2_early)

###############################
# 14. Treatment Protocol Analysis
###############################

# Examine treatment protocols among patients
table(final_df$age_at_diagnosis_dichot)
protocol_table <- final_df %>% 
  filter(!is.na(combined_protocol) & combined_protocol != "Other") %>% 
  group_by(age_at_diagnosis_dichot, combined_protocol) %>% 
  summarize(n = n())
protocol_table <- as.data.frame(protocol_table)
protocol_table
sum(protocol_table$n)

# Check for potential inconsistencies (e.g., children <10 on older protocols)
checkdf <- final_df %>% filter(combined_protocol == 'AALL1131-like (For children > 10 years)' & age_at_diagnosis < 10) %>% 
  dplyr::select(age_at_diagnosis, combined_protocol, age_at_diagnosis_dichot)
checkdf2 <- Bcell.df %>% filter(HR_cytogenetics == 'Yes')

###############################
# 15. Additional Outcome Analysis Using Poor Outcome 2
###############################

# Define poor_outcome2 as poor_outcome == 1 or Positive induction MRD
final_df <- final_df %>%
  mutate(poor_outcome2 = if_else((poor_outcome == 1) | (induction_mrd_result == 'Positive'), 1, 0))
table(final_df$poor_outcome2)

# Age summary by poor_outcome2
final_df %>%
  group_by(poor_outcome2, age_at_diagnosis_dichot) %>%
  summarize(
    count = n(),
    mean = mean(age_at_diagnosis),
    sd = sd(age_at_diagnosis),
    min = min(age_at_diagnosis),
    median = median(age_at_diagnosis),
    max = max(age_at_diagnosis)
  )
table_data <- table(final_df$poor_outcome2, final_df$age_at_diagnosis_dichot)
fisher.test(table_data)

# Gender summary by poor_outcome2
final_df$Gender <- as.factor(final_df$Gender)
final_df %>%
  group_by(poor_outcome2, Gender) %>%
  summarize(n = n()) %>%
  group_by(poor_outcome2) %>%
  mutate(percent = n / sum(n) * 100)
chisq.test(final_df$Gender, final_df$poor_outcome2)

# Ethnicity summary by poor_outcome2
final_df$Ethnicity <- as.factor(final_df$Ethnicity)
final_df %>%
  group_by(poor_outcome2, Ethnicity) %>%
  summarize(n = n()) %>%
  group_by(poor_outcome2) %>%
  mutate(percent = n / sum(n) * 100)
chisq.test(final_df$Ethnicity, final_df$poor_outcome2)

# Race summary by poor_outcome2
final_df %>%
  group_by(poor_outcome2, Race) %>%
  summarize(n = n()) %>%
  group_by(poor_outcome2) %>%
  mutate(percent = n / sum(n) * 100)
final_df %>%
  group_by(poor_outcome2, Race) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = Race, values_from = n, values_fill = 0) %>%
  dplyr::select(-poor_outcome2) %>%
  fisher.test()

# Diagnosis year summary by poor_outcome2
final_df %>%
  group_by(poor_outcome2) %>%
  summarize(
    count = n(),
    mean = mean(diagnosis_year),
    sd = sd(diagnosis_year),
    min = min(diagnosis_year),
    median = median(diagnosis_year),
    max = max(diagnosis_year),
    IQR = IQR(diagnosis_year)
  )
t.test(diagnosis_year ~ poor_outcome2, data = final_df)

###############################
# 16. Logistic Regression with Poor Outcome 2 (Sensitivity Analysis)
###############################

model3 <- glm(poor_outcome2 ~ Gender + Race + Ethnicity + diagnosis_year + age_at_diagnosis_dichot, 
              data = final_df, family = binomial(link = "logit"))
mod3_sum <- summary(model3)
mod3_pvalues <- mod3_sum$coefficients[, 4]
mod3_OR <- exp(mod3_sum$coefficients[, 1])
mod3_pvalues_df <- data.frame(variable = names(mod3_pvalues), pvalue = mod3_pvalues)
mod3_OR_df <- data.frame(variable = names(mod3_OR), estimate = mod3_OR)
mod3_results <- merge(mod3_pvalues_df, mod3_OR_df, by = "variable")
mod3_results[, 2] <- round(mod3_results[, 2], 2)
mod3_results$OR <- formatC(mod3_results$estimate, format = "f", digits = 2)
mod3_results <- mod3_results %>% dplyr::select(-estimate)
mod3_results

# End of Script
