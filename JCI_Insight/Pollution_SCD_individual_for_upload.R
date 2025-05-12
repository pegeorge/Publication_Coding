

############ Analyses ---------------------------------------------------


library(tidyverse)
library(mgcv)
library(splines)
library(fixest)
library(geepack)
library(scales) 


rm(list=ls())

load('Saved R Dataframes/annual_df.rda')


##  let's do some cleaning of the dataframe ------------------------------------------------------------------------------------

table(annual_df$genotype_cat)

#1220
annual_df %>%
  filter(genotype_dichot == 'SCA') %>%
  summarise(count = n_distinct(corp_id)) %>%
  pull(count)

summary(annual_df$years_in_df)

n_distinct(annual_df$corp_id)
# first, remove ED visits per year if the year is less than 0.2

pollution.df <- annual_df %>%
  group_by(corp_id) %>%
  filter(
    !(age_whole_number == min(age_whole_number) & first_year_duration < 0.2),
    !(age_whole_number == max(age_whole_number) & last_year_duration < 0.2)
  ) %>%
  ungroup()


#1188
pollution.df %>%
  filter(genotype_dichot == 'SCA') %>%
  summarise(count = n_distinct(corp_id)) %>%
  pull(count)


###  get the variables to be properly coded

# Convert 'TCD_abnormal_this_age' to a binary variable
pollution.df <- pollution.df %>%
  mutate(TCD_abnormal_binary = case_when(
    TCD_abnormal_this_age == "abnormal" ~ 1,
    TCD_abnormal_this_age == "not abnormal" ~ 0,
    TRUE ~ NA_real_  # Keeps NA as NA
  ))

# Convert 'TCD_abnormal_conditional_this_year' to a binary variable
pollution.df <- pollution.df %>%
  mutate(TCD_abnormal_conditional_binary = case_when(
    TCD_abnormal_conditional_this_age == "abnormal/conditional" ~ 1,
    TCD_abnormal_conditional_this_age == "not abnormal" ~ 0,
    TRUE ~ NA_real_  # Keeps NA as NA
  ))

sum(is.na(pollution.df$TCD_abnormal_this_age))
sum(is.na(pollution.df$TCD_abnormal_conditional_binary))

# change HU to factor variable
pollution.df$HU_this_year = as.factor(pollution.df$HU_this_year)

ls(pollution.df)

#remove rows with annualPM2.5 = 0
pollution.df = pollution.df %>% filter(is.na(pm25) == FALSE)

#1089
pollution.df %>%
  filter(genotype_dichot == 'SCA') %>%
  summarise(count = n_distinct(corp_id)) %>%
  pull(count)

n_distinct(pollution.df$corp_id)

ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = pm25)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Figure 1. Density Plot of Annual PM2.5") + 
  scale_x_continuous(name = "Annual PM2.5 Concentration (ug/m3)") +
  theme_classic()


#standardize the variables 
pollution.df$ec <- (pollution.df$ec - mean(pollution.df$ec)) / (2 * sd(pollution.df$ec))
pollution.df$oc <- (pollution.df$oc - mean(pollution.df$oc)) / (2 * sd(pollution.df$oc))
pollution.df$nh4 <- (pollution.df$nh4 - mean(pollution.df$nh4)) / (2 * sd(pollution.df$nh4))
pollution.df$so4 <- (pollution.df$so4 - mean(pollution.df$so4)) / (2 * sd(pollution.df$so4))
pollution.df$no3 <- (pollution.df$no3 - mean(pollution.df$no3)) / (2 * sd(pollution.df$no3))



##  let's see what this dataframe has
summary(pollution.df$distance_to_hosp_km)
table(pollution.df$genotype_cat)
summary(pollution.df$ED_hospital_days_per_year)
summary(pollution.df$hospital_days_per_year)
summary(pollution.df$ED_days_per_year)

table(pollution.df$TCD_abnormal_binary)
table(pollution.df$TCD_abnormal_conditional_binary)


##  add burden, for 1 year prior, 2 years prior, 3 years prior.  ------------------------------
###    note this is averaged burden (so 1 year / 1, 2 years 2, 3 years / 3)

ls(pollution.df)

pollution.df <- pollution.df %>%
  arrange(corp_id, age_whole_number) %>%  # Ensure data is in chronological order within each corp_id
  group_by(corp_id) %>%
  mutate(
    # For pm25
    pm25_1yr_burden = dplyr::lag(pm25, 1),
    pm25_2yr_burden = (dplyr::lag(pm25, 1) + dplyr::lag(pm25, 2))/2,
    pm25_3yr_burden = (dplyr::lag(pm25, 1) + dplyr::lag(pm25, 2) + dplyr::lag(pm25, 3))/3,
    
    # For ec
    ec_1yr_burden = dplyr::lag(ec, 1),
    ec_2yr_burden = (dplyr::lag(ec, 1) + dplyr::lag(ec, 2))/2,
    ec_3yr_burden = (dplyr::lag(ec, 1) + dplyr::lag(ec, 2) + dplyr::lag(ec, 3))/3,
    
    # For so4
    oc_1yr_burden = dplyr::lag(so4, 1),
    oc_2yr_burden = (dplyr::lag(so4, 1) + dplyr::lag(so4, 2))/2,
    oc_3yr_burden = (dplyr::lag(so4, 1) + dplyr::lag(so4, 2) + dplyr::lag(so4, 3))/3,
    
    # For nh4
    nh4_1yr_burden = dplyr::lag(nh4, 1),
    nh4_2yr_burden = (dplyr::lag(nh4, 1) + dplyr::lag(nh4, 2))/2,
    nh4_3yr_burden = (dplyr::lag(nh4, 1) + dplyr::lag(nh4, 2) + dplyr::lag(nh4, 3))/3,
    
    # For so4
    so4_1yr_burden = dplyr::lag(so4, 1),
    so4_2yr_burden = (dplyr::lag(so4, 1) + dplyr::lag(so4, 2))/2,
    so4_3yr_burden = (dplyr::lag(so4, 1) + dplyr::lag(so4, 2) + dplyr::lag(so4, 3))/3,
    
    # For no3
    no3_1yr_burden = dplyr::lag(no3, 1),
    no3_2yr_burden = (dplyr::lag(no3, 1) + dplyr::lag(no3, 2))/2,
    no3_3yr_burden = (dplyr::lag(no3, 1) + dplyr::lag(no3, 2) + dplyr::lag(no3, 3))/3,
    
    # For SVI_SES
    SES_1yr_burden = dplyr::lag(SVI_SES, 1),
    SES_2yr_burden = (dplyr::lag(SVI_SES, 1) + dplyr::lag(SVI_SES, 2))/2,
    SES_3yr_burden = (dplyr::lag(SVI_SES, 1) + dplyr::lag(SVI_SES, 2) + dplyr::lag(SVI_SES, 3))/3,
    
    # For SVI_household
    housecomp_1yr_burden = dplyr::lag( SVI_household, 1),
    housecomp_2yr_burden = (dplyr::lag( SVI_household, 1) + dplyr::lag( SVI_household, 2))/2,
    housecomp_3yr_burden = (dplyr::lag( SVI_household, 1) + dplyr::lag( SVI_household, 2) + dplyr::lag( SVI_household, 3))/3,
    
    # For  SVI_minority_language
    minorlang_1yr_burden = dplyr::lag( SVI_minority_language, 1),
    minorlang_2yr_burden = (dplyr::lag( SVI_minority_language, 1) + dplyr::lag( SVI_minority_language, 2))/2,
    minorlang_3yr_burden = (dplyr::lag( SVI_minority_language, 1) + dplyr::lag( SVI_minority_language, 2) + dplyr::lag( SVI_minority_language, 3))/3,
    
    # For  SVI_housing_transport
    housetransp_1yr_burden = dplyr::lag( SVI_housing_transport, 1),
    housetransp_2yr_burden = (dplyr::lag( SVI_housing_transport, 1) + dplyr::lag( SVI_housing_transport, 2))/2,
    housetransp_3yr_burden = (dplyr::lag( SVI_housing_transport, 1) + dplyr::lag( SVI_housing_transport, 2) + dplyr::lag( SVI_housing_transport, 3))/3,
    
    # For  SVI
    SVI_1yr_burden = dplyr::lag( SVI, 1),
    SVI_2yr_burden = (dplyr::lag( SVI, 1) + dplyr::lag( SVI, 2))/2,
    SVI_3yr_burden = (dplyr::lag( SVI, 1) + dplyr::lag( SVI, 2) + dplyr::lag( SVI, 3))/3,
    
    # For  ACS_persons_pov
    persons_pov_ACS1yr_burden = dplyr::lag( ACS_persons_pov, 1),
    persons_pov_ACS2yr_burden = (dplyr::lag( ACS_persons_pov, 1) + dplyr::lag( ACS_persons_pov, 2))/2,
    persons_pov_ACS3yr_burden = (dplyr::lag( ACS_persons_pov, 1) + dplyr::lag( ACS_persons_pov, 2) + dplyr::lag( ACS_persons_pov, 3))/3,
    
    # For SVI_SES
    percap_inc1yr_burden = dplyr::lag(ACS_per_cap_inc, 1),
    percap_inc2yr_burden = (dplyr::lag(ACS_per_cap_inc, 1) + dplyr::lag(ACS_per_cap_inc, 2))/2,
    percap_inc3yr_burden = (dplyr::lag(ACS_per_cap_inc, 1) + dplyr::lag(ACS_per_cap_inc, 2) + dplyr::lag(ACS_per_cap_inc, 3))/3,
    
  ) %>%
  ungroup()




## initial table ------------------------------------------------------------------------------
ls(pollution.df)

summary(pollution.df$calendar_year)

pollution.df %>%
  filter(genotype_dichot == 'SCA') %>%
  summarise(n_distinct_corp_id = n_distinct(corp_id))

# Continuous variables
continuous_vars <- c("hospital_days_per_year", 'ED_days_per_year',
                     'WBC_avg', "ANC_avg", 'Hgb_avg', 'creatinine_avg', 
                     "pm25", 
                     "age_whole_number",
                     "SVI", "SVI_SES", 'ADI_natrank',
                     "ACS_per_cap_inc", "distance_to_hosp_km")

continuous_summary <- pollution.df %>% 
  filter(genotype_dichot == 'SCA') %>% 
  dplyr::select(all_of(continuous_vars)) %>%
  summarise(across(.cols = everything(),
                   .fns = list(
                     n = ~sum(!is.na(.)), 
                     min = ~min(., na.rm = TRUE), 
                     `25th` = ~quantile(., probs = 0.25, na.rm = TRUE),
                     mean = ~mean(., na.rm = TRUE), 
                     median = ~median(., na.rm = TRUE), 
                     `75th` = ~quantile(., probs = 0.75, na.rm = TRUE),
                     max = ~max(., na.rm = TRUE), 
                     sd = ~sd(., na.rm = TRUE)
                   )))


reshaped_summary <- continuous_summary %>%
  pivot_longer(
    cols = everything(), 
    names_to = c("variable", "stat"), 
    names_pattern = "(.+)_(.+)"
  ) %>%
  pivot_wider(
    names_from = "stat", 
    values_from = "value"
  )


table_1 = reshaped_summary %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))

#rename the variables
table_1$variable <- c(
  "Annual inpatient hospital days",
  "Annual ED visits",
  "Average annual WBC",
  "Average annual ANC",
  "Average annual Hgb",
  "Average annual creatinine",
  "Annual PM2.5 exposure",
  "Age",
  "Social Vulnerability Index",
  "SVI SES",
  "ADI national rank",
  "ACS per capita income",
  "Distance to hospital (km)"
)

print(table_1)

write.csv(table_1, file = '../Project 2 Analyses and Drafts/Table1.csv')

# Categorical variables
unique_corp_id_df <- pollution.df%>% filter(genotype_dichot == 'SCA') %>% 
  arrange(corp_id, age_whole_number) %>%
  group_by(corp_id) %>%
  slice(1) %>%
  ungroup()


Race_df = as.data.frame(table(unique_corp_id_df$Race))
Sex_df = as.data.frame(table(unique_corp_id_df$Sex))
Sex_df
insurance_df = as.data.frame(table(unique_corp_id_df$insurance))
insurance_df


table1.df = pollution.df %>% filter(genotype_dichot == 'SCA')

table(table1.df$TCD_abnormal_binary, useNA = 'always')

m = mean(unique_corp_id_df$years_in_df)
m 
m*nrow(unique_corp_id_df)

tab_moved = table(unique_corp_id_df$has_moved_ever)
tab_moved
tab_moved[2]/(tab_moved[1]+tab_moved[2])*100

SCA_df = pollution.df %>% filter(genotype_cat == 'SCA') 

rm(SCA_df)


# check that data with visualizations ----------------------------------------------------------------------------------


# 1.  Is there enough variability in our independent variable of interest, PM2.5
# Step 1: Calculate the standard deviation of annual_PM_25 for each person
individual_variability <- pollution.df %>% 
  filter(genotype_dichot == 'SCA') %>% 
  group_by(corp_id) %>%
  summarize(individual_sd = sd(pm25, na.rm = TRUE),
            max_pm = max(pm25, na.rm = TRUE),
            min_pm = min(pm25, na.rm = TRUE),
            range_pm = max_pm - min_pm) %>%
  ungroup()


# Step 2: Calculate the average of these standard deviations
average_within_person_variability <- mean(individual_variability$individual_sd, na.rm = TRUE)
average_within_person_variability

# Join them to the original dataframe
pollution.df = left_join(pollution.df, individual_variability, by = "corp_id")


# Density plot of within-person standard deviation
#            not sure if this is enough individual_variability in terms of range - maybe? 


Supp_Fig_2 = ggplot(individual_variability, aes(x = range_pm)) +
  geom_histogram(fill = "blue", alpha = 0.5, binwidth = 0.1) +  
  labs(title = 'Figure 2. Histogram of PM2.5 Ranges per Individual',
       x = expression("Range of PM"[2.5]*" concentration ("*mu*"g/m"^3*")")) + # Update x axis label with subscripts and greek symbols
  theme_classic()


ggsave(filename = "../Project_2_Figures/Supp_Fig_2.png",
       plot = Supp_Fig_2,
       dpi = 600,         # This is the dots per inch; 300 is recommended for high quality
       width = 6,        # Width in inches; adjust based on your needs
       height = 4,        # Height in inches
       units = "in")      # Units for width and height (inches, cm, mm)

#2.  Check the outcome(s) of interest

str(pollution.df$age_whole_number)

Supp_Fig_4a = ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = hospital_days_per_year)) + 
  geom_histogram(fill = 'blue', alpha = 0.5, binwidth = 1) +
  theme_classic() +
  scale_y_sqrt() +
  labs(
    title = "Annual Inpatient Hospital Days",
    x = "Annual inpatient hospital days",
    y = "Count (Square root scale)"
  )

ggsave(filename = "../Project_2_Figures/Supp_Fig_4a.png",
       plot = Supp_Fig_4a,
       dpi = 600,         # This is the dots per inch; 300 is recommended for high quality
       width = 4,        # Width in inches; adjust based on your needs
       height = 3,        # Height in inches
       units = "in")      # Units for width and height (inches, cm, mm)

no_hosp_days_df = pollution.df %>% 
  filter(hospital_days_per_year == 0)

summary(no_hosp_days_df$ED_days_per_year)
summary(pollution.df$ED_days_per_year)

# Add the age_category variable to pollution.df
pollution.df <- pollution.df %>%
  mutate(
    age_category = case_when(
      age_whole_number >= 0 & age_whole_number <= 5  ~ "0-5",
      age_whole_number >= 6 & age_whole_number <= 10 ~ "6-10",
      age_whole_number >= 11 & age_whole_number <= 18 ~ "11-18",
      TRUE ~ "Unknown"  # for data outside the specified ranges
    )
  )

# Plot the data with the new age_category variable
ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = hospital_days_per_year)) +
  geom_histogram(fill = 'blue', alpha = 0.5, binwidth = 1) +
  facet_wrap(~age_category) +  # Creates a separate plot for each age category
  theme_classic() +
  scale_y_sqrt() +
  labs(
    title = "Histogram of Annual Inpatient Hospital Days by Age Category",
    x = "Annual inpatient hospital days",
    y = "Count (Square root scale)"
  )

Supp_Fig_4b = ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = ED_days_per_year)) + 
  geom_histogram(fill = 'blue', alpha = 0.5, binwidth = 1) +
  theme_classic() +
  scale_y_sqrt() +
  labs(
    title = "Annual ED Visits",
    x = "Annual ED visits",
    y = "Count (Square root scale)"
  )

ggsave(filename = "../Project_2_Figures/Supp_Fig_4b.png",
       plot = Supp_Fig_4b,
       dpi = 600,         # This is the dots per inch; 300 is recommended for high quality
       width = 4,        # Width in inches; adjust based on your needs
       height = 3,        # Height in inches
       units = "in")      # Units for width and height (inches, cm, mm)

pollution.df$ED_hospital_days_per_year
ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = ED_hospital_days_per_year)) + 
  geom_histogram(fill = 'blue', alpha = 0.5, binwidth = 1) +
  theme_classic() +
  scale_y_sqrt() +
  labs(
    title = "Histogram of Annual ED Visits 
+ Hospital Days Per Year",
    x = "Annual ED visits + Inpatient Hospital Days",
    y = "Count (Square root scale)"
  )


ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = hospital_days_per_year)) + 
  geom_density(fill = 'blue', alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Density Plot of Annual Inpatient Hospital Days",
    x = "Annual inpatient hospital days",
    y = "Density"
  )


# Winsorize the data and create binary variable
#note, hospital_days_per_year has several outliers, so I will Winsorize them
# such that all the values about the 99th percentile are set to the 95th percentile value

percentile_95 <- quantile(pollution.df$hospital_days_per_year, 0.95)

pollution.df <- pollution.df %>%
  mutate(HospDays_winsor = ifelse(hospital_days_per_year > percentile_95, percentile_95, hospital_days_per_year), 
         ED_HospDays_winsor = ifelse(ED_hospital_days_per_year > percentile_95, percentile_95, ED_hospital_days_per_year), 
         HospDays_binary = ifelse(hospital_days_per_year == 0, 0, 1))

SCA_df = pollution.df %>% filter(genotype_dichot == 'SCA')

tab_hbd = table(SCA_df$HospDays_binary, useNA = 'always')
tab_hbd
tab_hbd[2]/(tab_hbd[1]+tab_hbd[2])*100

tab_tcd = table(SCA_df$TCD_abnormal_binary)
tab_tcd
tab_tcd[2]/(tab_tcd[1]+tab_tcd[2])*100

tab_insur = table(SCA_df$insurance)
tab_insur
tab_insur[1]/(tab_insur[1]+tab_insur[2] + tab_insur[3])*100
tab_insur[2]/(tab_insur[1]+tab_insur[2] + tab_insur[3])*100
tab_insur[3]/(tab_insur[1]+tab_insur[2] + tab_insur[3])*100

summary(SCA_df$years_in_df)
summary(SCA_df$calendar_year)



checkdf = pollution.df %>% select(corp_id, hospital_days_per_year, HospDays_binary, HospDays_winsor)


ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = HospDays_winsor)) + 
  geom_histogram(fill = 'blue', alpha = 0.5, binwidth = 1) +
  theme_classic()

ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = ED_HospDays_winsor)) + 
  geom_histogram(fill = 'blue', alpha = 0.5, binwidth = 1) +
  theme_classic()

ggplot(pollution.df %>% filter(genotype_dichot == 'SCA'), aes(x = ED_days_per_year)) + 
  geom_histogram(fill = 'blue', alpha = 0.5, binwidth = 1) +
  theme_classic()



# other outcome variables
ggplot(pollution.df, aes(x = HospDays_winsor, fill = genotype_dichot)) + 
  geom_density(alpha = 0.5) +
  theme_classic()

ggplot(pollution.df, aes(x = Hgb_avg, fill = genotype_dichot)) + 
  geom_density(alpha = 0.5) +
  theme_classic()

ggplot(pollution.df, aes(x = WBC_avg, fill = genotype_dichot)) + 
  geom_density(alpha = 0.5) +
  theme_classic()

ggplot(pollution.df, aes(x = ANC_avg, fill = genotype_dichot)) + 
  geom_density(alpha = 0.5) +
  theme_classic()

ggplot(pollution.df, aes(x = creatinine_avg, fill = genotype_dichot)) + 
  geom_density(alpha = 0.5) +
  theme_classic()

# correlations between pollutants


pollutants_df <- pollution.df[,c("pm25", "ec", "nh4", "so4", "no3", "oc")]

# Compute the correlation matrix
cor_matrix <- cor(pollutants_df, use = "complete.obs") # Using complete.obs to handle missing values

library(reshape2) # For melting the correlation matrix

# Transform the correlation matrix for ggplot
melted_cor_matrix <- melt(cor_matrix)

str(melted_cor_matrix)

# Create the heatmap
ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + # Use geom_tile for the heatmap squares
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) + # Add correlation coefficients as text
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + # Use a minimal theme for a nicer look
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x axis labels for better readability
  labs(x = "", y = "", title = "Correlation Matrix of Pollutants", 
       caption = 'note that all pollutants are positively correlated (none are negative correlated)') +
  coord_fixed() # Ensure the aspect ratio is 1 to make the tiles square




# correlation matrix 2 
# Load necessary libraries
library(reshape2)

# Prepare the data
pollutants_df <- pollution.df[, c("pm25", "ec", "nh4", "so4", "no3", "oc")]

# Compute the correlation matrix
cor_matrix <- cor(pollutants_df, use = "complete.obs")

# Transform the correlation matrix for ggplot
melted_cor_matrix <- melt(cor_matrix)

# Create a named vector of labels with proper formatting
labels <- c(
  "pm25" = "PM[2.5]",
  "ec" = "EC",
  "nh4" = "NH[4]^'+'",
  "so4" = "SO[4]^'2-'",
  "no3" = "NO[3]^'-'",
  "oc" = "OC"
)

# Create the heatmap
Supp_Fig_7 = ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + # Use geom_tile for the heatmap squares
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) + # Add correlation coefficients as text
  scale_fill_gradient2(low = "white", high = "red", mid = "orange",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name = "Correlation") +
  scale_x_discrete(labels = function(x) parse(text = labels[x])) +
  scale_y_discrete(labels = function(x) parse(text = labels[x])) +
  theme_minimal() + # Use a minimal theme for a nicer look
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) + # Rotate x-axis labels and remove axis titles
  labs(
    title = "Correlation Matrix of Pollutants",
    caption = "Note that all pollutants are positively correlated (none are negatively correlated)"
  ) +
  coord_fixed() # Ensure the aspect ratio is 1 to make the tiles square


Supp_Fig_7 = ggsave(filename = "../Project_2_Figures/Supp_Fig_7.png",
                    plot = Supp_Fig_7,
                    dpi = 600,         # This is the dots per inch; 300 is recommended for high quality
                    width = 6,        # Width in inches; adjust based on your needs
                    height = 6,        # Height in inches
                    units = "in")      # Units for width and height (inches, cm, mm)


# modeling the results (variable names) ------------------------------------------------
# independent variables: pm25, ec, so4, nh4, so4, no3
# lagged IVs: nh4_1yr_burden (etc), pm25_1yr_burden
# other IVs: creatinine_avg, WBC_avg, Hgb_avg, ANC_avg, HgF_avg, TCD_abnormal_binary, TCD_abnormal_conditional_binary
# other independent vars:  SVI, SVI_1yr_burden, SVI_SES, SES_1yr_burden, insurance, 
# other Independent vars:  HU_this_year, distance_to_hosp_km, yearly_avg_temp
# dependent variables:Hgb_avg, HospDays_winsor, log_hospital_days, TCD_abnormal_binary

library(plm)
library(lmtest)
library(sandwich)
library(broom)

## dataframe to use (pollution.pdata) is a pdata.frame -------------------------
pollution.pdata = pollution.df %>% filter(distance_to_hosp_km < 48.2) 

pollution.pdata %>%
  filter(genotype_dichot == 'SCA') %>%
  summarise(count = n_distinct(corp_id)) %>%
  pull(count)

#1089


save(pollution.pdata, file =  'Saved R Dataframes/pollution.pdata.rda')



pollution.pdata = pdata.frame(pollution.pdata, index = c('corp_id', 'age_whole_number'))

#check how many missing, etc
table(pollution.pdata$age_whole_number, useNA = 'always')
sum(is.na(pollution.pdata$HU_this_year))





###########################################################################################################################
# Fit the fixed effects models, using hospital days per year and PM2.5  --------------------------------------------------------------------------------------------
###########################################################################################################################



# Start here for regression models 

# Annual values -----------------------------------------------------------------------------------------------------------





#  hospital_days_winsor
HospDays_mod <- plm(HospDays_winsor ~ pm25 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)


# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_pm25_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

HospDays_pm25_df

Hosp_pm25_df_IRR <- HospDays_pm25_df %>%
  mutate(
    IRR = exp(Estimate),
    `95% CI lower` = exp(Estimate - 1.96 * Std.Error),
    `95% CI upper` = exp(Estimate + 1.96 * Std.Error)
  ) %>%
  # Optionally round the new columns to 3 decimal places:
  mutate(across(c(IRR, `95% CI lower`, `95% CI upper`), round, 3))

# View the new dataframe
Hosp_pm25_df_IRR








#  hospital_days_binary

HospDays_mod <- plm(HospDays_binary ~ pm25 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDaysBinary_pm25_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  ED_hosp_days_per_year
HospDays_mod <- plm(ED_HospDays_winsor ~ pm25 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_pm25_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  ED_days_per_year
HospDays_mod <- plm(ED_days_per_year ~ pm25 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_pm25_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ pm25 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_pm25_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ pm25 + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_pm25_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ pm25 + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_pm25_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ pm25 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_pm25_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ pm25 + SVI + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_pm25_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# Annual values, interaction with HU -------------------------

#  hospital_days_per_year
HospDays_mod <- plm(HospDays_winsor ~ pm25 + SVI + 
                      as.factor(HU_this_year)*pm25 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_pm25_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year
HospDays_mod <- plm(HospDays_binary ~ pm25 + SVI + 
                      as.factor(HU_this_year)*pm25 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDaysBinary_pm25_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  ED_hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ pm25 + SVI + 
                      as.factor(HU_this_year)*pm25 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_pm25_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  ED_hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ pm25 + SVI + 
                      as.factor(HU_this_year)*pm25 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_pm25_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ pm25 + SVI + 
                      as.factor(HU_this_year)*pm25 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_pm25_interactHU_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ pm25 + SVI + 
                     as.factor(HU_this_year)*pm25  + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_pm25_interactHU_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ pm25 + SVI + 
                     as.factor(HU_this_year)*pm25 + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_pm25_interactHU_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ pm25 + SVI + 
                      as.factor(HU_this_year)*pm25 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_pm25_interactHU_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ pm25 + SVI + 
                   as.factor(HU_this_year)*pm25 + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_pm25_interactHU_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Annual values, subsetting by SES factors  -----------------------------------------------------------------------------------------------------------

# create a new dataframe that only keeps individuals with at least one observation of SVI > 0.50 (higher SVI = higher vulnerability)

# Identify the individuals with any observation of SVI >= 0.50 
individuals_with_high_svi <- pollution.df %>% filter(genotype_dichot == 'SCA') %>% 
  group_by(corp_id) %>%
  filter(any(SVI >= 0.50)) %>%
  ungroup() %>%
  select(corp_id) %>%
  distinct()

# Filter the dataset to keep all observations for these individuals
filtered_pollution_df <- pollution.df %>%
  semi_join(individuals_with_high_svi, by = "corp_id")

# Convert it back to a pdata.frame if needed for panel data analysis
high_SVI_pdata <- pdata.frame(filtered_pollution_df, index = c("corp_id", "age_whole_number"))

n_distinct(high_SVI_pdata$corp_id)


# create a new variable that is SVI_dichot, that == 1 if SVI >= 0.5, and == 0 if SVI <0.5

pollution.pdata = pollution.pdata %>% 
  mutate(SVI_dichot = if_else(SVI >= 0.5, 1, 0))

table(pollution.pdata$SVI_dichot)

#  hospital_days_per_year
HospDays_mod <- plm(HospDays_winsor ~ pm25*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)


# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_pm25_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_binary ~ pm25*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDaysBinary_pm25_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

HospDaysBinary_pm25_interactSVI_df

#  ED_hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ pm25*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_pm25_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  ED_hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ pm25*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_pm25_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ pm25*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_pm25_interactSVI_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ pm25*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_pm25_interactSVI_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ pm25*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_pm25_interactSVI_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ pm25*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_pm25_interactSVI_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ pm25*SVI_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_pm25_interactSVI_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# create a new variable that is annual_per_cap_inc_dichot, that == 1 per cap inc > national median 

#per census, median per capita income was $37,000 (5 year estimate) https://www.census.gov/quickfacts/fact/table/US/INC910221#INC910221 
summary(pollution.pdata$ACS_per_cap_inc)

pollution.pdata = pollution.pdata %>% 
  mutate(annual_per_cap_inc_dichot = if_else(ACS_per_cap_inc >= 37000, 1, 0))

table(pollution.pdata$annual_per_cap_inc_dichot)


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ pm25*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_pm25_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# ED_hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ pm25*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_pm25_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ pm25*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_pm25_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ pm25*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_pm25_interactinc_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ pm25*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_pm25_interactinc_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ pm25*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_pm25_interactinc_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ pm25*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_pm25_interactinc_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ pm25*annual_per_cap_inc_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_pm25_interactinc_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# create a new variable that is insurance dichot, that == 1 if income == private

table(pollution.pdata$insurance)

pollution.pdata = pollution.pdata %>% 
  mutate(insurance_dichot = if_else(insurance == 'Commercial', 1, 0))

table(pollution.pdata$insurance_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ pm25*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_pm25_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

#  ED_hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ pm25*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_pm25_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ pm25*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_pm25_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ pm25*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_pm25_interactinsur_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ pm25*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_pm25_interactinsur_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ pm25*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_pm25_interactinsur_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ pm25*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_pm25_interactinsur_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ pm25*insurance_dichot + 
                   as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_pm25_interactinsur_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





















###########################################################################################################################
# Fit the fixed effects models, using hospital days per year and annual ec  --------------------------------------------------------------------------------------------
###########################################################################################################################


# Annual values -----------------------------------------------------------------------------------------------------------

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ ec + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_ec_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_binary ~ ec + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDaysBinary_ec_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ ec + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_ec_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ ec + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_ec_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ ec + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_ec_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ ec + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_ec_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ ec + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_ec_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ ec + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_ec_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ ec + SVI + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_ec_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# Annual values, interaction with HU -------------------------


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ ec + SVI + 
                      as.factor(HU_this_year)*ec + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_ec_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  ED_hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ ec + SVI + 
                      as.factor(HU_this_year)*ec + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_ec_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ ec + SVI + 
                      as.factor(HU_this_year)*ec + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_ec_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ ec + SVI + 
                      as.factor(HU_this_year)*ec + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_ec_interactHU_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ ec + SVI + 
                     as.factor(HU_this_year)*ec  + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_ec_interactHU_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ ec + SVI + 
                     as.factor(HU_this_year)*ec + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_ec_interactHU_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ ec + SVI + 
                      as.factor(HU_this_year)*ec + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_ec_interactHU_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ ec + SVI + 
                   as.factor(HU_this_year)*ec + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_ec_interactHU_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Annual values, subsetting by SES factors  -----------------------------------------------------------------------------------------------------------

# create a new dataframe that only keeps individuals with at least one observation of SVI > 0.50 (higher SVI = higher vulnerability)

# Identify the individuals with any observation of SVI >= 0.50 
individuals_with_high_svi <- pollution.df %>% filter(genotype_dichot == 'SCA') %>% 
  group_by(corp_id) %>%
  filter(any(SVI >= 0.50)) %>%
  ungroup() %>%
  select(corp_id) %>%
  distinct()

# Filter the dataset to keep all observations for these individuals
filtered_pollution_df <- pollution.df %>%
  semi_join(individuals_with_high_svi, by = "corp_id")

# Convert it back to a pdata.frame if needed for panel data analysis
high_SVI_pdata <- pdata.frame(filtered_pollution_df, index = c("corp_id", "age_whole_number"))

n_distinct(high_SVI_pdata$corp_id)


# create a new variable that is SVI_dichot, that == 1 if SVI >= 0.5, and == 0 if SVI <0.5

pollution.pdata = pollution.pdata %>% 
  mutate(SVI_dichot = if_else(SVI >= 0.5, 1, 0))

table(pollution.pdata$SVI_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ ec*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_ec_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ ec*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_ec_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ ec*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_ec_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ ec*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_ec_interactSVI_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ ec*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_ec_interactSVI_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ ec*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_ec_interactSVI_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ ec*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_ec_interactSVI_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ ec*SVI_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_ec_interactSVI_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# create a new variable that is annual_per_cap_inc_dichot, that == 1 per cap inc > national median

#per census, median per capita income was $37,000 (5 year estimate) https://www.census.gov/quickfacts/fact/table/US/INC910221#INC910221 
summary(pollution.pdata$ACS_per_cap_inc)

pollution.pdata = pollution.pdata %>% 
  mutate(annual_per_cap_inc_dichot = if_else(ACS_per_cap_inc >= 37000, 1, 0))

table(pollution.pdata$annual_per_cap_inc_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ ec*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_ec_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ ec*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_ec_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ ec*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_ec_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ ec*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_ec_interactinc_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ ec*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_ec_interactinc_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ ec*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_ec_interactinc_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ ec*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_ec_interactinc_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ ec*annual_per_cap_inc_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_ec_interactinc_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# create a new variable that is insurance dichot, that == 1 if income == private

table(pollution.pdata$insurance)

pollution.pdata = pollution.pdata %>% 
  mutate(insurance_dichot = if_else(insurance == 'Commercial', 1, 0))

table(pollution.pdata$insurance_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ ec*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_ec_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ ec*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_ec_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ ec*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_ec_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ ec*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_ec_interactinsur_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ ec*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_ec_interactinsur_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ ec*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_ec_interactinsur_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ ec*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_ec_interactinsur_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ ec*insurance_dichot + 
                   as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_ec_interactinsur_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




###########################################################################################################################
# Fit the fixed effects models, using hospital days per year and annual OC  --------------------------------------------------------------------------------------------
###########################################################################################################################


# Annual values -----------------------------------------------------------------------------------------------------------

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ oc + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_oc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_binary ~ oc + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDaysBinary_oc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ oc + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_oc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ oc + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_oc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ oc + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_oc_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ oc + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_oc_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ oc + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_oc_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ oc + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_oc_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ oc + SVI + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_oc_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# Annual values, interaction with HU -------------------------


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ oc + SVI + 
                      as.factor(HU_this_year)*oc + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_oc_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ oc + SVI + 
                      as.factor(HU_this_year)*oc + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_oc_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ oc + SVI + 
                      as.factor(HU_this_year)*oc + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_oc_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ oc + SVI + 
                      as.factor(HU_this_year)*oc + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_oc_interactHU_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ oc + SVI + 
                     as.factor(HU_this_year)*oc  + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_oc_interactHU_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ oc + SVI + 
                     as.factor(HU_this_year)*oc + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_oc_interactHU_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ oc + SVI + 
                      as.factor(HU_this_year)*oc + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_oc_interactHU_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ oc + SVI + 
                   as.factor(HU_this_year)*oc + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_oc_interactHU_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Annual values, subsetting by SES factors  -----------------------------------------------------------------------------------------------------------

# create a new dataframe that only keeps individuals with at least one observation of SVI > 0.50 (higher SVI = higher vulnerability)

# Identify the individuals with any observation of SVI >= 0.50 
individuals_with_high_svi <- pollution.df %>% filter(genotype_dichot == 'SCA') %>% 
  group_by(corp_id) %>%
  filter(any(SVI >= 0.50)) %>%
  ungroup() %>%
  select(corp_id) %>%
  distinct()

# Filter the dataset to keep all observations for these individuals
filtered_pollution_df <- pollution.df %>%
  semi_join(individuals_with_high_svi, by = "corp_id")

# Convert it back to a pdata.frame if needed for panel data analysis
high_SVI_pdata <- pdata.frame(filtered_pollution_df, index = c("corp_id", "age_whole_number"))

n_distinct(high_SVI_pdata$corp_id)


# create a new variable that is SVI_dichot, that == 1 if SVI >= 0.5, and == 0 if SVI <0.5

pollution.pdata = pollution.pdata %>% 
  mutate(SVI_dichot = if_else(SVI >= 0.5, 1, 0))

table(pollution.pdata$SVI_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ oc*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_oc_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ oc*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_oc_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ oc*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_oc_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ oc*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_oc_interactSVI_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ oc*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_oc_interactSVI_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ oc*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_oc_interactSVI_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ oc*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_oc_interactSVI_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ oc*SVI_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_oc_interactSVI_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# create a new variable that is annual_per_cap_inc_dichot, that == 1 per cap inc > national median

#per census, median per capita income was $37,000 (5 year estimate) https://www.census.gov/quickfacts/fact/table/US/INC910221#INC910221 
summary(pollution.pdata$ACS_per_cap_inc)

pollution.pdata = pollution.pdata %>% 
  mutate(annual_per_cap_inc_dichot = if_else(ACS_per_cap_inc >= 37000, 1, 0))

table(pollution.pdata$annual_per_cap_inc_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ oc*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_oc_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ oc*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_oc_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ oc*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_oc_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ oc*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_oc_interactinc_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ oc*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_oc_interactinc_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ oc*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_oc_interactinc_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ oc*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_oc_interactinc_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ oc*annual_per_cap_inc_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_oc_interactinc_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# create a new variable that is insurance dichot, that == 1 if income == private

table(pollution.pdata$insurance)

pollution.pdata = pollution.pdata %>% 
  mutate(insurance_dichot = if_else(insurance == 'Commercial', 1, 0))

table(pollution.pdata$insurance_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ oc*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_oc_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ oc*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_oc_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ oc*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_oc_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ oc*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_oc_interactinsur_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ oc*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_oc_interactinsur_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ oc*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_oc_interactinsur_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ oc*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_oc_interactinsur_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ oc*insurance_dichot + 
                   as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_oc_interactinsur_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))






###########################################################################################################################
# Fit the fixed effects models, using hospital days per year and annual nh4  --------------------------------------------------------------------------------------------
###########################################################################################################################


# Annual values -----------------------------------------------------------------------------------------------------------

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ nh4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_nh4_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_binary ~ nh4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDaysBinary_nh4_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))






#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ nh4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_nh4_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ nh4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_nh4_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ nh4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_nh4_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ nh4 + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_nh4_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ nh4 + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_nh4_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ nh4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_nh4_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ nh4 + SVI + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_nh4_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# Annual values, interaction with HU -------------------------


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ nh4 + SVI + 
                      as.factor(HU_this_year)*nh4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_nh4_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ nh4 + SVI + 
                      as.factor(HU_this_year)*nh4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_nh4_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ nh4 + SVI + 
                      as.factor(HU_this_year)*nh4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_nh4_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ nh4 + SVI + 
                      as.factor(HU_this_year)*nh4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_nh4_interactHU_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ nh4 + SVI + 
                     as.factor(HU_this_year)*nh4  + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_nh4_interactHU_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ nh4 + SVI + 
                     as.factor(HU_this_year)*nh4 + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_nh4_interactHU_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ nh4 + SVI + 
                      as.factor(HU_this_year)*nh4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_nh4_interactHU_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ nh4 + SVI + 
                   as.factor(HU_this_year)*nh4 + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_nh4_interactHU_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Annual values, subsetting by SES factors  -----------------------------------------------------------------------------------------------------------

# create a new dataframe that only keeps individuals with at least one observation of SVI > 0.50 (higher SVI = higher vulnerability)

# Identify the individuals with any observation of SVI >= 0.50 
individuals_with_high_svi <- pollution.df %>% filter(genotype_dichot == 'SCA') %>% 
  group_by(corp_id) %>%
  filter(any(SVI >= 0.50)) %>%
  ungroup() %>%
  select(corp_id) %>%
  distinct()

# Filter the dataset to keep all observations for these individuals
filtered_pollution_df <- pollution.df %>%
  semi_join(individuals_with_high_svi, by = "corp_id")

# Convert it back to a pdata.frame if needed for panel data analysis
high_SVI_pdata <- pdata.frame(filtered_pollution_df, index = c("corp_id", "age_whole_number"))

n_distinct(high_SVI_pdata$corp_id)


# create a new variable that is SVI_dichot, that == 1 if SVI >= 0.5, and == 0 if SVI <0.5

pollution.pdata = pollution.pdata %>% 
  mutate(SVI_dichot = if_else(SVI >= 0.5, 1, 0))

table(pollution.pdata$SVI_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ nh4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_nh4_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ nh4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_nh4_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ nh4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_nh4_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ nh4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_nh4_interactSVI_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ nh4*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_nh4_interactSVI_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ nh4*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_nh4_interactSVI_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ nh4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_nh4_interactSVI_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ nh4*SVI_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_nh4_interactSVI_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# create a new variable that is annual_per_cap_inc_dichot, that == 1 per cap inc > national median

#per census, median per capita income was $37,000 (5 year estimate) https://www.census.gov/quickfacts/fact/table/US/INC910221#INC910221 
summary(pollution.pdata$ACS_per_cap_inc)

pollution.pdata = pollution.pdata %>% 
  mutate(annual_per_cap_inc_dichot = if_else(ACS_per_cap_inc >= 37000, 1, 0))

table(pollution.pdata$annual_per_cap_inc_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ nh4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_nh4_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ nh4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_nh4_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ nh4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_nh4_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ nh4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_nh4_interactinc_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ nh4*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_nh4_interactinc_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ nh4*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_nh4_interactinc_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ nh4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_nh4_interactinc_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ nh4*annual_per_cap_inc_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_nh4_interactinc_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# create a new variable that is insurance dichot, that == 1 if income == private

table(pollution.pdata$insurance)

pollution.pdata = pollution.pdata %>% 
  mutate(insurance_dichot = if_else(insurance == 'Commercial', 1, 0))

table(pollution.pdata$insurance_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ nh4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_nh4_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ nh4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_nh4_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ nh4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_nh4_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ nh4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_nh4_interactinsur_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ nh4*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_nh4_interactinsur_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ nh4*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_nh4_interactinsur_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ nh4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_nh4_interactinsur_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ nh4*insurance_dichot + 
                   as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_ec_interactinsur_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))










###########################################################################################################################
# Fit the fixed effects models, using hospital days per year and annual no3  --------------------------------------------------------------------------------------------
###########################################################################################################################


# Annual values -----------------------------------------------------------------------------------------------------------

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ no3 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_no3_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_binary ~ no3 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDaysBinary_no3_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ no3 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_no3_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ no3 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_no3_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ no3 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_no3_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ no3 + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_no3_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ no3 + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_no3_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ no3 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_no3_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ no3 + SVI + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_no3_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# Annual values, interaction with HU -------------------------


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ no3 + SVI + 
                      as.factor(HU_this_year)*no3 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_no3_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ no3 + SVI + 
                      as.factor(HU_this_year)*no3 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_no3_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ no3 + SVI + 
                      as.factor(HU_this_year)*no3 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_no3_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ no3 + SVI + 
                      as.factor(HU_this_year)*no3 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_no3_interactHU_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ no3 + SVI + 
                     as.factor(HU_this_year)*no3  + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_no3_interactHU_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ no3 + SVI + 
                     as.factor(HU_this_year)*no3 + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_no3_interactHU_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ no3 + SVI + 
                      as.factor(HU_this_year)*no3 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_no3_interactHU_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ no3 + SVI + 
                   as.factor(HU_this_year)*no3 + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_no3_interactHU_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Annual values, subsetting by SES factors  -----------------------------------------------------------------------------------------------------------

# create a new dataframe that only keeps individuals with at least one observation of SVI > 0.50 (higher SVI = higher vulnerability)

# Identify the individuals with any observation of SVI >= 0.50 
individuals_with_high_svi <- pollution.df %>% filter(genotype_dichot == 'SCA') %>% 
  group_by(corp_id) %>%
  filter(any(SVI >= 0.50)) %>%
  ungroup() %>%
  select(corp_id) %>%
  distinct()

# Filter the dataset to keep all observations for these individuals
filtered_pollution_df <- pollution.df %>%
  semi_join(individuals_with_high_svi, by = "corp_id")

# Convert it back to a pdata.frame if needed for panel data analysis
high_SVI_pdata <- pdata.frame(filtered_pollution_df, index = c("corp_id", "age_whole_number"))

n_distinct(high_SVI_pdata$corp_id)


# create a new variable that is SVI_dichot, that == 1 if SVI >= 0.5, and == 0 if SVI <0.5

pollution.pdata = pollution.pdata %>% 
  mutate(SVI_dichot = if_else(SVI >= 0.5, 1, 0))

table(pollution.pdata$SVI_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ no3*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_no3_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ no3*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_no3_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ no3*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_no3_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ no3*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_no3_interactSVI_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ no3*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_no3_interactSVI_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ no3*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_no3_interactSVI_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ no3*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_no3_interactSVI_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ no3*SVI_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_no3_interactSVI_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# create a new variable that is annual_per_cap_inc_dichot, that == 1 per cap inc > national median

#per census, median per capita income was $37,000 (5 year estimate) https://www.census.gov/quickfacts/fact/table/US/INC910221#INC910221 
summary(pollution.pdata$ACS_per_cap_inc)

pollution.pdata = pollution.pdata %>% 
  mutate(annual_per_cap_inc_dichot = if_else(ACS_per_cap_inc >= 37000, 1, 0))

table(pollution.pdata$annual_per_cap_inc_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ no3*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_no3_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ no3*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_no3_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ no3*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_no3_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ no3*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_no3_interactinc_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ no3*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_no3_interactinc_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ no3*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_no3_interactinc_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ no3*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_no3_interactinc_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ no3*annual_per_cap_inc_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_no3_interactinc_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# create a new variable that is insurance dichot, that == 1 if income == private

table(pollution.pdata$insurance)

pollution.pdata = pollution.pdata %>% 
  mutate(insurance_dichot = if_else(insurance == 'Commercial', 1, 0))

table(pollution.pdata$insurance_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ no3*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_no3_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ no3*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_no3_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ no3*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_no3_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ no3*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_no3_interactinsur_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ no3*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_no3_interactinsur_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ no3*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_no3_interactinsur_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ no3*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_no3_interactinsur_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ no3*insurance_dichot + 
                   as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_ec_interactinsur_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

































###########################################################################################################################
# Fit the fixed effects models, using hospital days per year and annual so4  --------------------------------------------------------------------------------------------
###########################################################################################################################


# Annual values -----------------------------------------------------------------------------------------------------------

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ so4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_so4_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_binary ~ so4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDaysBinary_so4_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ so4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_so4_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ so4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_so4_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ so4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_so4_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ so4 + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_so4_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ so4 + SVI + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_so4_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ so4 + SVI + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_so4_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ so4 + SVI + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_so4_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# Annual values, interaction with HU -------------------------


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ so4 + SVI + 
                      as.factor(HU_this_year)*so4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_so4_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ so4 + SVI + 
                      as.factor(HU_this_year)*so4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_so4_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ so4 + SVI + 
                      as.factor(HU_this_year)*so4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_so4_interactHU_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ so4 + SVI + 
                      as.factor(HU_this_year)*so4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_so4_interactHU_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ so4 + SVI + 
                     as.factor(HU_this_year)*so4  + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_so4_interactHU_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ so4 + SVI + 
                     as.factor(HU_this_year)*so4 + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_so4_interactHU_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ so4 + SVI + 
                      as.factor(HU_this_year)*so4 + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_so4_interactHU_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ so4 + SVI + 
                   as.factor(HU_this_year)*so4 + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_so4_interactHU_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Annual values, subsetting by SES factors  -----------------------------------------------------------------------------------------------------------

# create a new dataframe that only keeps individuals with at least one observation of SVI > 0.50 (higher SVI = higher vulnerability)

# Identify the individuals with any observation of SVI >= 0.50 
individuals_with_high_svi <- pollution.df %>% filter(genotype_dichot == 'SCA') %>% 
  group_by(corp_id) %>%
  filter(any(SVI >= 0.50)) %>%
  ungroup() %>%
  select(corp_id) %>%
  distinct()

# Filter the dataset to keep all observations for these individuals
filtered_pollution_df <- pollution.df %>%
  semi_join(individuals_with_high_svi, by = "corp_id")

# Convert it back to a pdata.frame if needed for panel data analysis
high_SVI_pdata <- pdata.frame(filtered_pollution_df, index = c("corp_id", "age_whole_number"))

n_distinct(high_SVI_pdata$corp_id)


# create a new variable that is SVI_dichot, that == 1 if SVI >= 0.5, and == 0 if SVI <0.5

pollution.pdata = pollution.pdata %>% 
  mutate(SVI_dichot = if_else(SVI >= 0.5, 1, 0))

table(pollution.pdata$SVI_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ so4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_so4_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ so4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_so4_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ so4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_so4_interactSVI_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ so4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_so4_interactSVI_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ so4*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_so4_interactSVI_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ so4*SVI_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_so4_interactSVI_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ so4*SVI_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_so4_interactSVI_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ so4*SVI_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_so4_interactSVI_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# create a new variable that is annual_per_cap_inc_dichot, that == 1 per cap inc > national median

#per census, median per capita income was $37,000 (5 year estimate) https://www.census.gov/quickfacts/fact/table/US/INC910221#INC910221 
summary(pollution.pdata$ACS_per_cap_inc)

pollution.pdata = pollution.pdata %>% 
  mutate(annual_per_cap_inc_dichot = if_else(ACS_per_cap_inc >= 37000, 1, 0))

table(pollution.pdata$annual_per_cap_inc_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ so4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_so4_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ so4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_so4_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ so4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_so4_interactinc_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ so4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_so4_interactinc_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ so4*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_so4_interactinc_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ so4*annual_per_cap_inc_dichot + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_so4_interactinc_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ so4*annual_per_cap_inc_dichot + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_so4_interactinc_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ so4*annual_per_cap_inc_dichot + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_so4_interactinc_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))




# create a new variable that is insurance dichot, that == 1 if income == private

table(pollution.pdata$insurance)

pollution.pdata = pollution.pdata %>% 
  mutate(insurance_dichot = if_else(insurance == 'Commercial', 1, 0))

table(pollution.pdata$insurance_dichot)

#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ so4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_so4_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_HospDays_winsor ~ so4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

EDHosp_so4_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


#  hospital_days_per_year

HospDays_mod <- plm(ED_days_per_year ~ so4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ED_days_per_year_so4_interactinsur_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))





# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ so4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_so4_interactinsur_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ so4*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_so4_interactinsur_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ so4*insurance_dichot + 
                     as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_so4_interactinsur_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ so4*insurance_dichot + 
                      as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_abnl_so4_interactinsur_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

# creatinine_avg

creat_mod <- plm(creatinine_avg ~ so4*insurance_dichot + 
                   as.factor(HU_this_year) + SVI + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_ec_interactinsur_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))

















































































# 1 year burden values -----------------------------------------------------------------------------------------------------------

# independent variables: pm25, ec, oc, nh4, so4, no3
# lagged IVs: nh4_1yr_burden (etc), pm25_1yr_burden
# other IVs: creatinine_avg, WBC_avg, Hgb_avg, ANC_avg, HgF_avg, TCD_abnormal_binary, TCD_abnormal_conditional_binary
# other independent vars:  SVI, SVI_1yr_burden, SVI_SES, SES_1yr_burden, insurance, 
# other Independent vars:  HU_this_year, distance_to_hosp_km, yearly_avg_temp
# dependent variables:Hgb_avg, HospDays_winsor, log_hospital_days, TCD_abnormal_binary


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ pm25_1yr_burden + SVI_1yr_burden + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_pm25_1yr_burden_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ pm25_1yr_burden + SVI_1yr_burden + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_pm25_1yr_burden_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ pm25_1yr_burden + SVI_1yr_burden + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_pm25_1yr_burden_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ pm25_1yr_burden + SVI_1yr_burden + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_pm25_1yr_burden_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ pm25_1yr_burden + SVI_1yr_burden + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_pm25_1yr_burden_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ pm25_1yr_burden + SVI_1yr_burden + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_pm25_1yr_burden_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))







# 2 year burden values -----------------------------------------------------------------------------------------------------------

# independent variables: pm25, ec, oc, nh4, so4, no3
# lagged IVs: nh4_1yr_burden (etc), pm25_1yr_burden
# other IVs: creatinine_avg, WBC_avg, Hgb_avg, ANC_avg, HgF_avg, TCD_abnormal_binary, TCD_abnormal_conditional_binary
# other independent vars:  SVI, SVI_1yr_burden, SVI_SES, SES_1yr_burden, insurance, 
# other Independent vars:  HU_this_year, distance_to_hosp_km, yearly_avg_temp
# dependent variables:Hgb_avg, HospDays_winsor, log_hospital_days, TCD_abnormal_binary


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ pm25_2yr_burden + SVI_2yr_burden + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_pm25_2yr_burden_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ pm25_2yr_burden + SVI_2yr_burden +  
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_pm25_2yr_burden_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ pm25_2yr_burden + SVI_2yr_burden + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_pm25_2yr_burden_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ pm25_2yr_burden + SVI_2yr_burden + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_pm25_2yr_burden_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ pm25_2yr_burden + SVI_2yr_burden +  
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_pm25_2yr_burden_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ pm25_2yr_burden + SVI_2yr_burden + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_pm25_2yr_burden_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# 3 year burden values -----------------------------------------------------------------------------------------------------------

# independent variables: pm25, ec, oc, nh4, so4, no3
# lagged IVs: nh4_1yr_burden (etc), pm25_1yr_burden
# other IVs: creatinine_avg, WBC_avg, Hgb_avg, ANC_avg, HgF_avg, TCD_abnormal_binary, TCD_abnormal_conditional_binary
# other independent vars:  SVI, SVI_1yr_burden, SVI_SES, SES_1yr_burden, insurance, 
# other Independent vars:  HU_this_year, distance_to_hosp_km, yearly_avg_temp
# dependent variables:Hgb_avg, HospDays_winsor, log_hospital_days, TCD_abnormal_binary


#  hospital_days_per_year

HospDays_mod <- plm(HospDays_winsor ~ pm25_3yr_burden + SVI_3yr_burden + 
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = quasipoisson(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary( HospDays_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_hospital_days_per_year = coeftest(HospDays_mod, vcov = vcovHC(HospDays_mod, type = "HC1", cluster = "group"))

estimates <- robust_hospital_days_per_year[, "Estimate"]
std_errors <- robust_hospital_days_per_year[, "Std. Error"]
statistics <- robust_hospital_days_per_year[, "t value"]
p_values <- robust_hospital_days_per_year[, "Pr(>|t|)"]
robust_HospDays_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

HospDays_pm25_3yr_burden_df <- robust_HospDays_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# WBC_avg

WBC_pm25_mod <- plm(WBC_avg ~ pm25_3yr_burden + SVI_3yr_burden +  
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = gaussian(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(WBC_pm25_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_WBC_avg = coeftest(WBC_pm25_mod, vcov = vcovHC(WBC_pm25_mod, type = "HC1", cluster = "group"))

estimates <- robust_WBC_avg[, "Estimate"]
std_errors <- robust_WBC_avg[, "Std. Error"]
statistics <- robust_WBC_avg[, "t value"]
p_values <- robust_WBC_avg[, "Pr(>|t|)"]
robust_WBC_pm25_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

WBC_pm25__3yr_burden_df <- robust_WBC_pm25_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# ANC_avg

ANC_avg_mod <- plm(ANC_avg ~ pm25_3yr_burden + SVI_3yr_burden + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(ANC_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_ANC_avg = coeftest(ANC_avg_mod, vcov = vcovHC(ANC_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_ANC_avg[, "Estimate"]
std_errors <- robust_ANC_avg[, "Std. Error"]
statistics <- robust_ANC_avg[, "t value"]
p_values <- robust_ANC_avg[, "Pr(>|t|)"]
robust_ANC_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

ANC_pm25_3yr_burden_df <- robust_ANC_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



# Hgb_avg

Hgb_avg_mod <- plm(Hgb_avg ~ pm25_3yr_burden + SVI_3yr_burden + 
                     as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                     as.numeric(age_whole_number),
                   model = 'within', 
                   family = gaussian(),
                   data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(Hgb_avg_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_Hgb_avg = coeftest(Hgb_avg_mod, vcov = vcovHC(Hgb_avg_mod, type = "HC1", cluster = "group"))

estimates <- robust_Hgb_avg[, "Estimate"]
std_errors <- robust_Hgb_avg[, "Std. Error"]
statistics <- robust_Hgb_avg[, "t value"]
p_values <- robust_Hgb_avg[, "Pr(>|t|)"]
robust_Hgb_avg_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

Hgb_pm25_3yr_burden_df <- robust_Hgb_avg_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# TCD_abnormal_binary

TCD_abnl_mod <- plm(TCD_abnormal_binary ~ pm25_3yr_burden + SVI_3yr_burden +  
                      as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                      as.numeric(age_whole_number),
                    model = 'within', 
                    family = binomial(),
                    data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(TCD_abnl_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_TCD_abnormal_binary = coeftest(TCD_abnl_mod, vcov = vcovHC(TCD_abnl_mod, type = "HC1", cluster = "group"))

estimates <- robust_TCD_abnormal_binary[, "Estimate"]
std_errors <- robust_TCD_abnormal_binary[, "Std. Error"]
statistics <- robust_TCD_abnormal_binary[, "t value"]
p_values <- robust_TCD_abnormal_binary[, "Pr(>|t|)"]
robust_TCD_abnl_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

TCD_pm25_3yr_burden_df <- robust_TCD_abnl_df %>%
  mutate(across(where(is.numeric), round, digits = 3))


# creatinine_avg

creat_mod <- plm(creatinine_avg ~ pm25_3yr_burden + SVI_3yr_burden + 
                   as.factor(HU_this_year) + insurance + distance_to_hosp_km + yearly_min_temp +
                   as.numeric(age_whole_number),
                 model = 'within', 
                 family = gaussian(),
                 data = pollution.pdata %>% filter(genotype_dichot == 'SCA'))

summary(creat_mod)

# Calculate robust standard errors clustered by 'corp_id'
robust_creatinine_avg = coeftest(creat_mod, vcov = vcovHC(creat_mod, type = "HC1", cluster = "group"))

estimates <- robust_creatinine_avg[, "Estimate"]
std_errors <- robust_creatinine_avg[, "Std. Error"]
statistics <- robust_creatinine_avg[, "t value"]
p_values <- robust_creatinine_avg[, "Pr(>|t|)"]
robust_creat_df <- data.frame(Estimate = estimates, Std.Error = std_errors, t.value = statistics, P.value = p_values)

creat_pm25_3yr_burden_df <- robust_creat_df %>%
  mutate(across(where(is.numeric), round, digits = 3))



########## lets visualize these! ---------------------------------

# List names of all dataframes in the environment
dataframe_names <- Filter(function(x) is.data.frame(get(x)), ls())

# Define the prefixes to filter by
prefixes <- c("^HospDays", '^EDHosp', '^ED_days_per_year', "^TCD", "^WBC", "^ANC", "^creat", "^Hgb")

# Create a pattern string from the prefixes for use in grep
pattern <- paste(prefixes, collapse="|")

# Filter dataframe names by the pattern
filtered_names <- grep(pattern, dataframe_names, value=TRUE)

# Print the filtered list of dataframe names
print(filtered_names)

# Initialize an empty list to store the modified dataframes
modified_dfs <- list()

# Loop through each dataframe name
for (df_name in filtered_names) {
  # Retrieve the dataframe by its name
  df <- get(df_name)
  
  # Extract the outcome variable name from the dataframe name (before '_')
  outcome_variable <- sub("_.*", "", df_name)
  
  # Add a new column with the outcome variable name
  df$outcome_variable <- outcome_variable
  
  # Append the modified dataframe to the list
  modified_dfs[[df_name]] <- df
}

# Combine all the modified dataframes into one large dataframe
combined_df <- do.call(rbind, modified_dfs)

# Calculate the 95% confidence interval bounds
combined_df$CI_lower <- combined_df$Estimate - (1.96 * combined_df$Std.Error)
combined_df$CI_upper <- combined_df$Estimate + (1.96 * combined_df$Std.Error)

# If the model names indicating interaction terms are in the row names
combined_df$has_interaction <- grepl("interact", rownames(combined_df))

head(combined_df)

library(dplyr)

#create main_exposure_variable
combined_df <- combined_df %>%
  mutate(main_exposure_variable = case_when(
    grepl("_pm25", rownames(.)) ~ "pm25",
    grepl("_ec", rownames(.)) ~ "ec",
    grepl("_oc", rownames(.)) ~ "oc",
    grepl("_nh4", rownames(.)) ~ "nh4",
    grepl("_no3", rownames(.)) ~ "no3",
    grepl("_so4", rownames(.)) ~ "so4",
    TRUE ~ NA_character_  # default case if none of the above matches
  ))

table(combined_df$main_exposure_variable)

# Add a new variable 'model_name' by extracting the name before the first dot in the row names
combined_df$model_name <- sub("\\..*$", "", rownames(combined_df))

# Add a new variable 'variable_name' by extracting the name after the first dot in the row names
combined_df$variable_name <- sub("^[^.]*\\.", "", rownames(combined_df))

combined_df$variable_name <- sub("^(annual_|Annual_|as.numeric|as.factor)", "", combined_df$variable_name)


# save
regression_results_df_v2_SCA_only = combined_df

save(regression_results_df_v2_SCA_only, file =  'Saved R Dataframes/regression_results_df_v2_SCA_only.rda')






############ Vizualizations ----------------------------------------------



# 1. Load libraries -------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(patchwork)
library(stringr)

# 2. Initialize workspace ------------------------------------------------------
rm(list = ls())
here()  # confirm project root

# 3. Load data -----------------------------------------------------------------
load(here("..", "DB", "Saved R Dataframes",
          "regression_results_df_v2_SCA_only.rda"))

# 4. Export raw table ----------------------------------------------------------
write.csv(
  regression_results_df_v2_SCA_only,
  file = here("..", "R code project 2 individual",
              "regression_results_df_v2_SCA_only.csv"),
  row.names = FALSE
)

# 5. Pre‑processing ------------------------------------------------------------
df1 <- regression_results_df_v2_SCA_only %>%
  filter(!grepl("burden", model_name)) %>%
  group_by(outcome_variable) %>%
  mutate(variable_order = case_when(
    variable_name == "pm25" ~ 1,
    variable_name == "oc"   ~ 2,
    variable_name == "ec"   ~ 3,
    variable_name == "no3"  ~ 4,
    variable_name == "nh4"  ~ 5,
    variable_name == "so4"  ~ 6,
    TRUE                    ~ 7
  )) %>%
  ungroup()

# 6. Define PM2.5 components order ---------------------------------------------
pm25_components <- c("ec", "oc", "no3", "nh4", "so4")

# 7. Subsets for PM2.5, no interactions ---------------------------------------

## 7a. Inflammatory outcomes (betas)
primary_WBC_ANC_df <- df1 %>%
  filter(variable_name == "pm25",
         has_interaction == FALSE,
         outcome_variable %in% c("ANC", "WBC"))

## 7b. Count outcomes (IRRs)
primary_ED_Hosp_df <- df1 %>%
  filter(variable_name == "pm25",
         has_interaction == FALSE,
         outcome_variable %in% c("EDHosp", "ED", "HospDays")) %>%
  mutate(
    Transformed_Estimate = exp(Estimate),
    CI_lower_transformed = exp(Estimate - qnorm(0.975) * Std.Error),
    CI_upper_transformed = exp(Estimate + qnorm(0.975) * Std.Error),
    outcome_variable = recode(outcome_variable,
                              "EDHosp"   = "ED + Hospital Days",
                              "ED"       = "ED Visits",
                              "HospDays" = "Hospital Days"
    )
  )

# drop combined category for Fig_1a
primary_ED_Hosp_df2 <- primary_ED_Hosp_df %>%
  filter(outcome_variable != "ED + Hospital Days")

## 7c. Binary outcomes (ORs)
primary_Binary_df <- df1 %>%
  filter(variable_name == "pm25",
         has_interaction == FALSE,
         outcome_variable %in% c("HospDaysBinary", "TCD")) %>%
  mutate(
    Transformed_Estimate = exp(Estimate),
    CI_lower_transformed = exp(Estimate - qnorm(0.975) * Std.Error),
    CI_upper_transformed = exp(Estimate + qnorm(0.975) * Std.Error),
    outcome_variable = recode(outcome_variable,
                              "HospDaysBinary" = "Hospitalization (Yes)",
                              "TCD"            = "Transcranial Doppler (abnormal)"
    )
  )

# 8. Figure 1 panels -----------------------------------------------------------

## 1a: Count outcomes (IRRs)
Fig_1a <- ggplot(primary_ED_Hosp_df2,
                 aes(x = outcome_variable, y = Transformed_Estimate,
                     color = outcome_variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower_transformed,
                    ymax = CI_upper_transformed),
                width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.05, nudge_x = 0.2) +
  labs(
    title   = "Effect of PM2.5 on Clinical Outcome Variables",
    x       = "Outcome Variable",
    y       = "Incident Rate Ratio",
    caption = "Quasi‑Poisson fixed effects regressions"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

ggsave(
  here("..", "Project_2_Figures", "Fig_1a.tiff"),
  plot  = Fig_1a, dpi = 600,
  width = 5.5, height = 3, units = "in"
)

## 1b: Inflammatory outcomes (betas)
Fig_1b <- ggplot(primary_WBC_ANC_df,
                 aes(x = outcome_variable, y = Estimate,
                     color = outcome_variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.005, nudge_x = 0.2) +
  labs(
    title   = "Effect of PM2.5 on Inflammatory Outcome Variables",
    x       = "Outcome Variable",
    y       = "Beta Coefficient",
    caption = "Fixed effects regressions"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

ggsave(
  here("..", "Project_2_Figures", "Fig_1b.tiff"),
  plot  = Fig_1b, dpi = 600,
  width = 5.5, height = 3, units = "in"
)

## 1c: Binary outcomes (ORs)
Fig_1c <- ggplot(primary_Binary_df,
                 aes(x = outcome_variable, y = Transformed_Estimate,
                     color = outcome_variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower_transformed,
                    ymax = CI_upper_transformed),
                width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.005, nudge_x = 0.2) +
  labs(
    title   = "Effect of PM2.5 on Binary Outcome Variables",
    x       = "Outcome Variable",
    y       = "Odds Ratio",
    caption = "Logistic fixed effects regressions"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

ggsave(
  here("..", "Project_2_Figures", "Fig_1c.tiff"),
  plot  = Fig_1c, dpi = 600,
  width = 6, height = 3, units = "in"
)

## Combine Figure 1
Figure_1_combined <- Fig_1a / Fig_1b / Fig_1c
ggsave(
  here("..", "Project_2_Figures", "Figure_1_combined.tiff"),
  plot  = Figure_1_combined, dpi = 600,
  width = 9, height = 9, units = "in"
)

# 9. Additional plots (no interactions) ----------------------------------------

## All outcomes, beta estimates for pm25
ggplot(df1 %>%
         filter(variable_name == "pm25",
                has_interaction == FALSE,
                main_exposure_variable == "pm25"),
       aes(x = outcome_variable, y = Estimate, color = outcome_variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.04, nudge_x = 0.4) +
  labs(
    title   = "Effect of PM2.5 on Different Outcome Variables",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = "Beta coefficient estimates and p-values"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

## All outcomes, beta for SVI
ggplot(df1 %>%
         filter(variable_name == "SVI",
                has_interaction == FALSE,
                main_exposure_variable == "pm25"),
       aes(x = outcome_variable, y = Estimate, color = outcome_variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.2, nudge_x = 0.3) +
  labs(
    title   = "Effect of SVI on Different Outcome Variables",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = "Beta coefficient estimates and p-values"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

# 10. Supplemental Fig 5: HU effect check --------------------------------------
Supp_Fig_5 <- ggplot(df1 %>%
                       filter(variable_name == "(HU_this_year)1",
                              has_interaction == FALSE,
                              main_exposure_variable == "pm25",
                              outcome_variable %in% c("WBC", "ANC", "Hgb")),
                     aes(x = outcome_variable, y = Estimate,
                         color = outcome_variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.05, nudge_x = 0.3) +
  labs(
    title   = "Effect of Hydroxyurea on Different Outcome Variables",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = "Beta coefficient estimates and p-values"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

ggsave(
  here("..", "Project_2_Figures", "Supp_Fig_5.png"),
  plot  = Supp_Fig_5, dpi = 600,
  width = 5.5, height = 3, units = "in"
)

# 11. Multi‑component plot (no interactions) -----------------------------------

# Prepare outcome order
ordered_outcome_vars <- c("ANC", "WBC", "TCD", "ED", "HospDays")

# Plot all components
ggplot(df1 %>%
         filter(variable_name %in% c("pm25", pm25_components),
                has_interaction == FALSE),
       aes(x = outcome_variable, y = Estimate,
           color = variable_name,
           group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title   = "Effect of PM2.5 Components on Different Outcome Variables",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = paste(
      "95% CI by error bars. Models adjust for insurance, distance to hospital,",
      "average annual temperature, patient age, patient insurance, annual SVI"
    )
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.title = element_text(text = "Component"),
        axis.text.x  = element_text(angle = 0, hjust = 1))

# 12. Fig 4: Ordered components plot -------------------------------------------

# Define desired plotting order
desired_order <- c("pm25", "so4", "no3", "nh4", "ec", "oc")

df_plot <- df1 %>%
  filter(variable_name %in% desired_order,
         has_interaction == FALSE) %>%
  mutate(
    outcome_variable     = factor(outcome_variable,
                                  levels = ordered_outcome_vars),
    variable_name_plot   = factor(variable_name,
                                  levels = rev(desired_order)),
    variable_name        = factor(variable_name,
                                  levels = desired_order)
  ) %>%
  filter(!is.na(outcome_variable))

dodge_width <- 0.6
dodge <- position_dodge(width = dodge_width)

Fig_4 <- ggplot(df_plot,
                aes(x = outcome_variable, y = Estimate,
                    color = variable_name,
                    group = interaction(outcome_variable, variable_name_plot))) +
  geom_point(size = 4, position = dodge) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = dodge) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Effect of PM2.5 Components on Different Outcome Variables",
    x     = "Outcome Variable",
    y     = "Estimate",
    color = "Component"
  ) +
  scale_color_discrete(labels = c(
    expression(PM[2.5]),
    expression(SO[4]^{"2-"}),
    expression(NO[3]^{"-"}),
    expression(NH[4]^{"+"}),
    "EC", "OC"
  )) +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

ggsave(
  here("..", "Project_2_Figures", "Fig_4.tiff"),
  plot  = Fig_4, dpi = 600,
  width = 6, height = 6, units = "in"
)

# 13. Sub‑plots by outcome group ------------------------------------------------

# WBC & ANC
ggplot(df1 %>%
         filter(variable_name %in% c("pm25", pm25_components),
                has_interaction == FALSE,
                outcome_variable %in% c("WBC", "ANC")),
       aes(x = outcome_variable, y = Estimate,
           color = variable_name,
           group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title   = "Effect of PM2.5 Components on Inflammatory Outcomes",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = paste(
      "95% CI by error bars. Models adjust for insurance, distance to hospital,",
      "average annual temperature, patient age, patient insurance, annual SVI"
    )
  ) +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# ED, HospDays, EDHosp
ggplot(df1 %>%
         filter(variable_name %in% c("pm25", pm25_components),
                has_interaction == FALSE,
                outcome_variable %in% c("ED", "HospDays", "EDHosp")),
       aes(x = outcome_variable, y = Estimate,
           color = variable_name,
           group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title   = "Effect of PM2.5 Components on ED & Hospital Outcomes",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = paste(
      "95% CI by error bars. Models adjust for insurance, distance to hospital,",
      "average annual temperature, patient age, patient insurance, annual SVI"
    )
  ) +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# HospDaysBinary & TCD
ggplot(df1 %>%
         filter(variable_name %in% c("pm25", pm25_components),
                has_interaction == FALSE,
                outcome_variable %in% c("HospDaysBinary", "TCD")),
       aes(x = outcome_variable, y = Estimate,
           color = variable_name,
           group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title   = "Effect of PM2.5 Components on Binary Outcomes",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = paste(
      "95% CI by error bars. Models adjust for insurance, distance to hospital,",
      "average annual temperature, patient age, patient insurance, annual SVI"
    )
  ) +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# 14. Interaction terms plots --------------------------------------------------

## 14a. HU interaction
interact_HU_df <- df1 %>%
  filter(grepl(":", variable_name),
         grepl("HU_this_year", variable_name),
         main_exposure_variable == "pm25",
         outcome_variable %in% c("WBC", "ANC", "HospDays",
                                 "TCD", "ED", "EDHosp")) %>%
  mutate(
    outcome_variable = factor(outcome_variable,
                              levels = c("ANC", "WBC", "TCD", "ED", "HospDays")
    ),
    outcome_variable = recode(outcome_variable,
                              "HospDays" = "Hospital Days",
                              "TCD"      = "Transcranial Doppler"
    )
  )

Fig_2a <- ggplot(interact_HU_df,
                 aes(x = outcome_variable, y = Estimate,
                     color = outcome_variable,
                     group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.055, nudge_x = 0.3) +
  labs(
    title   = "HU Interaction Term Estimates and 95% CI",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = "Beta coefficient of HU interaction term"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

ggsave(
  here("..", "Project_2_Figures", "Fig_2a.tiff"),
  plot  = Fig_2a, dpi = 600,
  width = 5.5, height = 3.5, units = "in"
)

## 14b. SVI interaction
interact_SVI_df <- df1 %>%
  filter(grepl(":", variable_name),
         grepl("SVI", variable_name),
         main_exposure_variable == "pm25",
         outcome_variable %in% c("WBC", "ANC", "HospDays",
                                 "TCD", "ED", "EDHosp")) %>%
  mutate(
    outcome_variable = factor(outcome_variable,
                              levels = c("ANC", "WBC", "TCD", "ED", "HospDays")
    ),
    outcome_variable = recode(outcome_variable,
                              "HospDays" = "Hospital Days",
                              "TCD"      = "Transcranial Doppler"
    )
  )

Fig_2b <- ggplot(interact_SVI_df,
                 aes(x = outcome_variable, y = Estimate,
                     color = outcome_variable,
                     group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.055, nudge_x = 0.3) +
  labs(
    title   = "SVI Interaction Term Estimates and 95% CI",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = "Beta coefficient of SVI interaction term"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

ggsave(
  here("..", "Project_2_Figures", "Fig_2b.tiff"),
  plot  = Fig_2b, dpi = 600,
  width = 5.5, height = 3.5, units = "in"
)

## 14c. Insurance interaction
insurance_SVI_df <- df1 %>%
  filter(grepl(":", variable_name),
         grepl("insurance", variable_name),
         main_exposure_variable == "pm25",
         outcome_variable %in% c("WBC", "ANC", "HospDays",
                                 "TCD", "ED", "EDHosp")) %>%
  mutate(
    outcome_variable = factor(outcome_variable,
                              levels = c("ANC", "WBC", "TCD", "ED", "HospDays")
    ),
    outcome_variable = recode(outcome_variable,
                              "HospDays" = "Hospital Days",
                              "TCD"      = "Transcranial Doppler"
    )
  )

Fig_2c <- ggplot(insurance_SVI_df,
                 aes(x = outcome_variable, y = Estimate,
                     color = outcome_variable,
                     group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = paste0("p=", round(P.value, 3))),
            size = 3.5, nudge_y = 0.055, nudge_x = 0.3) +
  labs(
    title   = "Insurance Interaction Term Estimates and 95% CI",
    x       = "Outcome Variable",
    y       = "Estimate",
    caption = "Beta coefficient of insurance interaction term"
  ) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

ggsave(
  here("..", "Project_2_Figures", "Fig_2c.tiff"),
  plot  = Fig_2c, dpi = 600,
  width = 5.5, height = 3.5, units = "in"
)

## Combine Figure 2
Figure_2_combined <- Fig_2a / Fig_2b / Fig_2c
ggsave(
  here("..", "Project_2_Figures", "Figure_2_combined.tiff"),
  plot  = Figure_2_combined, dpi = 600,
  width = 9, height = 9, units = "in"
)

# 15. Burden plots --------------------------------------------------------------

burden_variable_names <- c("pm25",
                           "pm25_1yr_burden",
                           "pm25_2yr_burden",
                           "pm25_3yr_burden")

burden_df <- regression_results_df_v2_SCA_only %>%
  filter(variable_name %in% burden_variable_names,
         has_interaction == FALSE) %>%
  mutate(variable_order = case_when(
    variable_name == "pm25"            ~ 1,
    variable_name == "pm25_1yr_burden" ~ 2,
    variable_name == "pm25_2yr_burden" ~ 3,
    variable_name == "pm25_3yr_burden" ~ 4
  ))

Supp_Fig_6 <- ggplot(burden_df %>%
                       filter(outcome_variable %in% c("WBC", "ANC", "HospDays", "TCD")),
                     aes(x = outcome_variable, y = Estimate,
                         color = variable_name,
                         group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title    = "Effect of Cumulative PM2.5 on Outcomes",
    subtitle = "Burden = previous year(s) average value",
    x        = "Outcome Variable",
    y        = "Estimate",
    color    = "Burden Years",
    caption  = "95% CI by error bars"
  ) +
  scale_color_manual(
    labels = c(
      expression(PM[2.5]~"current year"),
      expression(PM[2.5]~"1-year burden"),
      expression(PM[2.5]~"2-year burden"),
      expression(PM[2.5]~"3-year burden")
    ),
    values = c("purple", "green", "blue", "red")
  ) +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

ggsave(
  here("..", "Project_2_Figures", "Supp_Fig_6.png"),
  plot  = Supp_Fig_6, dpi = 600,
  width = 6, height = 4.5, units = "in"
)

# 16. TCD burden (binary ORs) --------------------------------------------------

TCD_burden_df <- burden_df %>%
  filter(outcome_variable == "TCD") %>%
  mutate(
    Transformed_Estimate = exp(Estimate),
    CI_lower_transformed = exp(Estimate - qnorm(0.975) * Std.Error),
    CI_upper_transformed = exp(Estimate + qnorm(0.975) * Std.Error),
    outcome_variable = recode(outcome_variable,
                              "TCD" = "Transcranial Doppler (abnormal)"
    )
  )

Fig_3 <- ggplot(TCD_burden_df,
                aes(x = outcome_variable,
                    y = Transformed_Estimate,
                    color = variable_name,
                    group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lower_transformed,
                    ymax = CI_upper_transformed),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    title    = "Effect of Cumulative PM2.5 on TCD (Burden)",
    subtitle = "Burden = previous year(s) average value",
    x        = "Outcome Variable",
    y        = "Odds Ratio",
    color    = "Burden Years",
    caption  = "95% CI by error bars"
  ) +
  scale_color_manual(
    labels = c(
      expression(PM[2.5]~"current year"),
      expression(PM[2.5]~"1-year burden"),
      expression(PM[2.5]~"2-year burden"),
      expression(PM[2.5]~"3-year burden")
    ),
    values = c("purple", "green", "blue", "red")
  ) +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

ggsave(
  here("..", "Project_2_Figures", "Fig_3.tiff"),
  plot  = Fig_3, dpi = 600,
  width = 6, height = 4, units = "in"
)

# 17. Interaction summary table -----------------------------------------------

interaction_df <- df1 %>%
  filter(variable_name %in% c("pm25", pm25_components),
         grepl("interact", model_name)) %>%
  group_by(model_name) %>%
  summarise(
    pollutant           = first(variable_name),
    pollutant_estimate  = first(Estimate),
    interaction_term    = nth(variable_name, 2),
    interaction_estimate = nth(Estimate, 2),
    combined_estimate   = pollutant_estimate + interaction_estimate,
    CI_lower_pollutant  = first(CI_lower),
    CI_upper_pollutant  = first(CI_upper),
    CI_lower_combined   = first(CI_lower) + interaction_estimate,
    CI_upper_combined   = first(CI_upper) + interaction_estimate,
    .groups = "drop"
  )

HU_interact_df <- interaction_df %>%
  filter(grepl("HU_this_year", interaction_term))

transformed_df <- HU_interact_df %>%
  mutate(variable_type = list(c("pollutant", "pollutant_plus_interaction"))) %>%
  unnest_longer(variable_type) %>%
  mutate(
    estimate = if_else(variable_type == "pollutant",
                       pollutant_estimate, combined_estimate),
    CI_lower = if_else(variable_type == "pollutant",
                       CI_lower_pollutant, CI_lower_combined),
    CI_upper = if_else(variable_type == "pollutant",
                       CI_upper_pollutant, CI_upper_combined),
    outcome_variable = str_extract(model_name, "^[^_]+")
  ) %>%
  mutate(
    y_order = if_else(row_number() %% 2 == 1,
                      row_number(),
                      row_number() - 0.5),
    variable_type = recode(variable_type,
                           "pollutant"                 = "pollutant_noHU",
                           "pollutant_plus_interaction" = "pollutant_HU"
    )
  )

ggplot(transformed_df %>%
         filter(outcome_variable %in% c("WBC", "ANC"),
                pollutant == "pm25"),
       aes(x = estimate, y = y_order,
           color = variable_type,
           shape = variable_type)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ outcome_variable, scales = "free") +
  labs(
    title   = "Visualizing HU Interaction Terms",
    x       = "Estimate and 95% CI",
    y       = NULL,
    color   = "Variable Type",
    shape   = "Variable Type",
    caption = paste(
      "pollutant_noHU = effect in non-HU group; pollutant_HU = combined effect"
    )
  ) +
  theme_classic() +
  theme(axis.text.y = element_blank())

# 18. SVI interaction overview ------------------------------------------------

interaction_SVI_df <- df1 %>%
  filter(grepl(":", variable_name),
         grepl("SVI", variable_name)) %>%
  mutate(variable_order = case_when(
    main_exposure_variable == "pm25" ~ 1,
    main_exposure_variable == "oc"   ~ 2,
    main_exposure_variable == "ec"   ~ 3,
    main_exposure_variable == "no3"  ~ 4,
    main_exposure_variable == "nh4"  ~ 5,
    main_exposure_variable == "so4"  ~ 6
  ))

ggplot(interaction_SVI_df,
       aes(x = outcome_variable, y = Estimate,
           color = variable_name,
           group = interaction(outcome_variable, variable_order))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title    = "SVI Interaction Term Estimates and 95% CI",
    subtitle = "Higher SVI = higher vulnerability",
    x        = "Outcome Variable",
    y        = "Estimate",
    caption  = paste(
      "95% CI by error bars. Models adjust for insurance,",
      "distance, temperature, age, insurance, annual SVI"
    )
  ) +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
