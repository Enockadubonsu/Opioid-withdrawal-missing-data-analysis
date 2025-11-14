
setwd("C:/Users/adubo/Desktop/PhD BIOSTATS/BIOS 576 C/Project 2/Results")


################################################################################
# Lofexidine for Opioid Withdrawal: Rigorous Longitudinal Analysis
# with Comprehensive Missing Data Sensitivity Analyses
#
# Project: BIOS 576C Project 2
# Analysis: Reanalysis of Yu et al. (2008) Phase III RCT
# Date: October 2024
################################################################################

# ==============================================================================
# 0. SETUP AND PACKAGE LOADING
# ==============================================================================

# Clear environment
rm(list = ls())

# Set seed for reproducibility
set.seed(12345)

# Load required packages
library(tidyverse)
library(nlme)
library(lme4)
library(emmeans)
library(mice)
library(broom.mixed)
library(gt)
library(gtsummary)
library(ggpubr)
library(naniar)
library(VIM)
library(gridExtra)
library(cowplot)

cat("✓ Packages loaded successfully\n\n")

# ==============================================================================
# 1. DATA IMPORT
# ==============================================================================

cat("=== STEP 1: DATA IMPORT ===\n")

# UPDATE THESE PATHS TO YOUR DATA LOCATION
base_end <- read.csv("C:/Users/adubo/Desktop/PhD BIOSTATS/BIOS 576 C/Project 2/base_end.csv", 
                     stringsAsFactors = FALSE)
withdraw <- read.csv("C:/Users/adubo/Desktop/PhD BIOSTATS/BIOS 576 C/Project 2/withdraw.csv", 
                     stringsAsFactors = FALSE)
emesis <- read.csv("C:/Users/adubo/Desktop/PhD BIOSTATS/BIOS 576 C/Project 2/emesis.csv", 
                   stringsAsFactors = FALSE)

cat("✓ Data loaded successfully\n")
cat("  Base/End records:", nrow(base_end), "\n")
cat("  Withdraw records:", nrow(withdraw), "\n")
cat("  Emesis records:", nrow(emesis), "\n\n")

# ==============================================================================
# 2. DATA CLEANING AND MHOWS CALCULATION
# ==============================================================================

cat("=== STEP 2: DATA CLEANING ===\n")

# ------------------------------------------------------------------------------
# 2.1 Clean base_end dataset
# ------------------------------------------------------------------------------

base_clean <- base_end %>%
  mutate(
    PATIENT = as.character(PATIENT),
    treat = factor(treat, levels = c("PLACEBO", "LOFEXIDINE")),
    gender = factor(gender, levels = c(0, 1), labels = c("Male", "Female")),
    depression = as.logical(depression),
    anxiety = as.logical(anxiety),
    IV = as.logical(IV),
    oral_nasal_smoke = as.logical(oral_nasal_smoke),
    smoke = as.logical(smoke),
    dropped = ifelse(!is.na(study_day_drop) & study_day_drop < 11, 1, 0),
    drop_reason_label = factor(drop_reason, 
                               levels = 1:9,
                               labels = c("Completed protocol", "Became ineligible", 
                                          "Toxicity/side effects", "Medical reason (unrelated)",
                                          "Non-compliance", "Request - med not working",
                                          "Request - other reason", "Death", "Other reason")),
    educ_yr = as.numeric(educ_yr),
    age = as.numeric(age)
  ) %>%
  select(PATIENT, treat, age, gender, educ_yr, depression, anxiety, 
         IV, oral_nasal_smoke, smoke, days_use_30d, 
         study_day_drop, drop_reason, drop_reason_label, dropped)

cat("✓ Base dataset cleaned. N =", nrow(base_clean), "\n")

# ------------------------------------------------------------------------------
# 2.2 Clean withdraw dataset - CORRECTED SCORING per Protocol Table 5
# ------------------------------------------------------------------------------

withdraw_clean <- withdraw %>%
  mutate(
    PATIENT = as.character(PATIENT),
    treat = factor(treat, levels = c("PLACEBO", "LOFEXIDINE")),
    day = as.numeric(day),
    
    # Discrete symptoms scoring (Protocol Table 5)
    q1_score = ifelse(Q1 == 2, 1, 0),   # Yawning: 1 point
    q2_score = ifelse(Q2 == 2, 1, 0),   # Lacrimation: 1 point
    q3_score = ifelse(Q3 == 2, 1, 0),   # Rhinorrhea: 1 point
    q4_score = ifelse(Q4 == 2, 1, 0),   # Perspiration: 1 point
    q5_score = ifelse(Q5 == 2, 3, 0),   # Tremor: 3 points
    q6_score = ifelse(Q6 == 2, 3, 0),   # Goose-flesh: 3 points
    q7_score = ifelse(Q7 == 2, 5, 0),   # Restlessness: 5 points
    anorexia_score = ifelse(Q8A == 1 | Q8B == 1 | Q8C == 1, 3, 0),  # Anorexia: 3 points
    
    discrete_symptoms = q1_score + q2_score + q3_score + q4_score + 
      q5_score + q6_score + q7_score + anorexia_score,
    
    continuous_signs = pts_pupil + pts_wt + pts_temp + pts_resp + pts_sysBP
  ) %>%
  select(PATIENT, treat, day, 
         Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8A, Q8B, Q8C, pupil,
         q1_score, q2_score, q3_score, q4_score, q5_score, q6_score, q7_score,
         anorexia_score, discrete_symptoms, 
         pts_pupil, pts_wt, pts_temp, pts_resp, pts_sysBP, continuous_signs,
         Q12DIA, Q12HR, Q14TEMP)

cat("✓ Withdraw dataset cleaned. N =", nrow(withdraw_clean), "\n")

# ------------------------------------------------------------------------------
# 2.3 Clean emesis dataset
# ------------------------------------------------------------------------------

emesis_clean <- emesis %>%
  mutate(
    PATIENT = as.character(PATIENT),
    treat = factor(treat, levels = c("PLACEBO", "LOFEXIDINE")),
    day = as.numeric(STDYDAY)
  ) %>%
  rowwise() %>%
  mutate(
    emesis_episodes = sum(!is.na(c_across(Q1:Q18)) & c_across(Q1:Q18) != "", na.rm = TRUE),
    emesis_points = case_when(
      is.na(emesis_episodes) | emesis_episodes == 0 ~ 0,
      emesis_episodes == 1 ~ 5,
      emesis_episodes == 2 ~ 10,
      emesis_episodes >= 3 ~ 15,
      TRUE ~ 0
    )
  ) %>%
  ungroup() %>%
  select(PATIENT, treat, day, emesis_episodes, emesis_points)

cat("✓ Emesis dataset cleaned. N =", nrow(emesis_clean), "\n")

# ------------------------------------------------------------------------------
# 2.4 Merge and calculate total MHOWS
# ------------------------------------------------------------------------------

mhows_data <- withdraw_clean %>%
  left_join(emesis_clean, by = c("PATIENT", "treat", "day")) %>%
  mutate(
    emesis_points = replace_na(emesis_points, 0),
    MHOWS = discrete_symptoms + continuous_signs + emesis_points
  )

analysis_data <- mhows_data %>%
  left_join(base_clean, by = c("PATIENT", "treat"))

cat("✓ Datasets merged. N observations =", nrow(analysis_data), "\n\n")
#dim(analysis_data)
#head(analysis_data)
# ------------------------------------------------------------------------------
# 2.5 Create analysis dataset with Day 3 as baseline
# ------------------------------------------------------------------------------

baseline_mhows <- analysis_data %>%
  filter(day == 3) %>%
  select(PATIENT, baseline_MHOWS = MHOWS)

analysis_long <- analysis_data %>%
  filter(day >= 3 & day <= 8) %>%
  left_join(baseline_mhows, by = "PATIENT") %>%
  mutate(
    time = day - 3,
    time_factor = factor(day),
    change_from_baseline = MHOWS - baseline_MHOWS,
    treat_numeric = ifelse(treat == "LOFEXIDINE", 1, 0),
    time_numeric = as.numeric(time),
    baseline_MHOWS_centered = baseline_MHOWS - mean(baseline_MHOWS, na.rm = TRUE)
  )

# Keep only participants with baseline
participants_with_baseline <- analysis_long %>%
  filter(day == 3) %>%
  pull(PATIENT) %>%
  unique()

analysis_long <- analysis_long %>%
  filter(PATIENT %in% participants_with_baseline)

cat("✓ Analysis dataset created\n")
cat("  Participants:", length(unique(analysis_long$PATIENT)), "\n")
cat("  Observations:", nrow(analysis_long), "\n\n")

# Check baseline MHOWS
day3_summary <- analysis_long %>%
  filter(day == 3) %>%
  summarise(
    N = n(),
    Mean = mean(MHOWS, na.rm = TRUE),
    SD = sd(MHOWS, na.rm = TRUE),
    Min = min(MHOWS, na.rm = TRUE),
    Max = max(MHOWS, na.rm = TRUE)
  )

cat("Day 3 Baseline MHOWS:\n")
cat("  Mean:", round(day3_summary$Mean, 2), "\n")
cat("  SD:", round(day3_summary$SD, 2), "\n")
cat("  Range:", round(day3_summary$Min, 2), "-", round(day3_summary$Max, 2), "\n\n")

if (day3_summary$Mean < 25) {
  cat("NOTE: Baseline MHOWS lower than Yu et al. (~30-35)\n")
  cat("      This is acceptable - your subset has milder withdrawal\n")
  cat("      Document this as a limitation in your report\n\n")
}
#
analysis_long <- analysis_long %>%
  mutate(treat = factor(treat, levels = c("PLACEBO", "LOFEXIDINE")))

# Save cleaned data
save(analysis_long, base_clean, baseline_mhows, file = "cleaned_analysis_data.RData")

levels(analysis_long$treat)
#dim(analysis_long)
#head(analysis_long)
# ==============================================================================
# 3. DESCRIPTIVE STATISTICS
# ==============================================================================

cat("=== STEP 3: DESCRIPTIVE STATISTICS ===\n")

# ------------------------------------------------------------------------------
# 3.1 Baseline characteristics (Table 1)
# ------------------------------------------------------------------------------

baseline_data <- analysis_long %>%
  filter(day == 3) %>%
  select(PATIENT, treat, age, gender, educ_yr, depression, anxiety,
         IV, oral_nasal_smoke, smoke, days_use_30d, baseline_MHOWS) %>%
  distinct()

table1 <- baseline_data %>%
  select(-PATIENT) %>%
  tbl_summary(
    by = treat,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age, years",
      gender ~ "Gender",
      educ_yr ~ "Education, years",
      depression ~ "Depression",
      anxiety ~ "Anxiety",
      IV ~ "Intravenous use",
      oral_nasal_smoke ~ "Oral/nasal/smoke use",
      smoke ~ "Smoking",
      days_use_30d ~ "Days of use in past 30",
      baseline_MHOWS ~ "MHOWS Day 3 (baseline)"
    )
  ) %>%
  add_overall() %>%
  modify_caption("**Table 1. Baseline Characteristics by Treatment Arm**")

print(table1)

table1 %>%
  as_gt() %>%
  gtsave("Table1_Baseline_Characteristics.html")

cat("✓ Table 1 created\n")

# ------------------------------------------------------------------------------
# 3.2 Missing data summary (Table 2)
# ------------------------------------------------------------------------------

missing_summary <- analysis_long %>%
  group_by(day, treat) %>%
  summarise(n_observed = sum(!is.na(MHOWS)), .groups = "drop") %>%
  pivot_wider(names_from = treat, values_from = n_observed, names_prefix = "n_") %>%
  mutate(
    Total = n_PLACEBO + n_LOFEXIDINE,
    Percent_Placebo = round(100 * n_PLACEBO / max(n_PLACEBO), 1),
    Percent_Lofexidine = round(100 * n_LOFEXIDINE / max(n_LOFEXIDINE), 1)
  )

cat("\n--- Table 2: Missing Data Summary ---\n")
print(missing_summary)
write.csv(missing_summary, "Table2_Missing_Data_Summary.csv", row.names = FALSE)

complete_cases <- analysis_long %>%
  group_by(PATIENT, treat) %>%
  summarise(n_obs = n(), has_all_days = n_obs == 6, .groups = "drop") %>%
  group_by(treat) %>%
  summarise(total = n(), complete = sum(has_all_days), pct_complete = round(100 * complete / total, 1))

cat("\n--- Complete Cases (Days 3-8) ---\n")
print(complete_cases)

# ------------------------------------------------------------------------------
# 3.3 Predictors of missingness
# ------------------------------------------------------------------------------

missingness_analysis <- analysis_long %>%
  filter(day %in% 4:8) %>%
  mutate(missing = is.na(MHOWS)) %>%
  filter(!is.na(baseline_MHOWS))

miss_model <- glm(
  missing ~ treat + baseline_MHOWS + age + gender + IV + days_use_30d,
  data = missingness_analysis,
  family = binomial(link = "logit")
)

miss_results <- broom::tidy(miss_model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(across(where(is.numeric), round, 3))

cat("\n--- Predictors of Missingness (OR) ---\n")
print(miss_results)
write.csv(miss_results, "Missingness_Predictors.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 3.4 Unadjusted means (Table 3)
# ------------------------------------------------------------------------------

table3 <- analysis_long %>%
  group_by(day, treat) %>%
  summarise(
    N = sum(!is.na(MHOWS)),
    Mean = mean(MHOWS, na.rm = TRUE),
    SD = sd(MHOWS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Mean_SD = sprintf("%.2f (%.2f)", Mean, SD)) %>%
  select(day, treat, N, Mean_SD) %>%
  pivot_wider(names_from = treat, values_from = c(N, Mean_SD), names_glue = "{treat}_{.value}")

cat("\n--- Table 3: Unadjusted Means ---\n")
print(table3)
write.csv(table3, "Table3_Unadjusted_Means.csv", row.names = FALSE)

cat("\n✓ Descriptive statistics completed\n\n")

# ==============================================================================
# 4. VISUALIZATIONS
# ==============================================================================

cat("=== STEP 4: VISUALIZATIONS ===\n")

# CONSORT numbers
n_randomized <- length(unique(baseline_data$PATIENT))
n_by_arm <- baseline_data %>% count(treat)

cat("CONSORT Numbers:\n")
cat("  Total randomized:", n_randomized, "\n")
print(n_by_arm)

# Spaghetti plot
fig2 <- ggplot(analysis_long %>% filter(!is.na(MHOWS)),
               aes(x = day, y = MHOWS, group = PATIENT)) +
  geom_line(aes(color = treat), alpha = 0.3) +
  facet_wrap(~treat) +
  scale_color_manual(values = c("PLACEBO" = "#E74C3C", "LOFEXIDINE" = "#3498DB")) +
  labs(title = "Individual MHOWS Trajectories by Treatment Arm",
       subtitle = "Days 3-8 (Day 3 = Baseline)",
       x = "Study Day", y = "MHOWS Score", color = "Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("Figure2_Individual_Trajectories.png", fig2, width = 10, height = 6, dpi = 300)

# Mean trajectories
mean_traj <- analysis_long %>%
  group_by(day, treat) %>%
  summarise(N = sum(!is.na(MHOWS)), Mean = mean(MHOWS, na.rm = TRUE),
            SD = sd(MHOWS, na.rm = TRUE), SE = SD / sqrt(N),
            CI_lower = Mean - 1.96 * SE, CI_upper = Mean + 1.96 * SE, .groups = "drop")

fig3 <- ggplot(mean_traj, aes(x = day, y = Mean, color = treat, fill = treat)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("PLACEBO" = "#E74C3C", "LOFEXIDINE" = "#3498DB")) +
  scale_fill_manual(values = c("PLACEBO" = "#E74C3C", "LOFEXIDINE" = "#3498DB")) +
  labs(title = "Mean MHOWS Trajectories by Treatment Arm",
       subtitle = "Error bars represent 95% confidence intervals",
       x = "Study Day", y = "Mean MHOWS Score", color = "Treatment", fill = "Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("Figure3_Mean_Trajectories.png", fig3, width = 10, height = 6, dpi = 300)

cat("✓ Figures created\n\n")



#NEW 

# ==============================================================================
# 5. PRIMARY ANALYSIS: LINEAR MIXED MODEL (CORRECTED)
# ==============================================================================

cat("=== STEP 5: PRIMARY ANALYSIS ===\n")

# Verify treatment coding before fitting model
cat("Treatment group verification:\n")
cat("  Levels:", levels(analysis_long$treat), "\n")
cat("  Counts:", table(analysis_long$treat), "\n")
cat("  Expected: PLACEBO first, LOFEXIDINE second\n\n")

# Fit primary model
primary_model <- lme(
  fixed = MHOWS ~ treat * time_factor,
  random = ~ 1 + time_numeric | PATIENT,
  data = analysis_long %>% filter(day >= 3),
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
)

cat("✓ Primary model fitted\n")

# Estimated marginal means
emm_primary <- emmeans(primary_model, ~ treat | time_factor)
contrasts_daily <- contrast(emm_primary, method = "pairwise", by = "time_factor")

daily_results <- summary(contrasts_daily, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(Day = as.numeric(as.character(time_factor)), across(where(is.numeric), round, 3))

write.csv(daily_results, "Primary_Daily_Treatment_Effects.csv", row.names = FALSE)

# Average treatment effect (Days 4-8)
emm_avg <- emmeans(primary_model, ~ treat, at = list(time_factor = as.character(4:8)))

# CRITICAL: Check what contrast we're computing
contrast_check <- contrast(emm_avg, method = "pairwise")
cat("\nChecking contrast direction:\n")
print(summary(contrast_check))

# Get the contrast result
contrast_result <- summary(contrast_check, infer = TRUE)
contrast_label <- as.character(contrast_result$contrast[1])

cat("\nContrast computed:", contrast_label, "\n")

# If contrast is "PLACEBO - LOFEXIDINE", we need to flip it
# We want "LOFEXIDINE - PLACEBO" (negative = beneficial)
if (grepl("PLACEBO.*LOFEXIDINE", contrast_label)) {
  cat("Flipping sign (PLACEBO - LOFEXIDINE → LOFEXIDINE - PLACEBO)\n")
  primary_result <- data.frame(
    estimate = -contrast_result$estimate,
    SE = contrast_result$SE,
    df = contrast_result$df,
    lower.CL = -contrast_result$upper.CL,  # Note: swap and negate
    upper.CL = -contrast_result$lower.CL,
    t.ratio = -contrast_result$t.ratio,
    p.value = contrast_result$p.value
  )
} else {
  # Already in correct direction
  primary_result <- as.data.frame(contrast_result)
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║           PRIMARY ANALYSIS RESULT (CORRECTED)                ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Average Treatment Effect (Days 4-8):                    ║\n"))
cat(sprintf("║    Estimate: %7.2f MHOWS points                        ║\n", primary_result$estimate))
cat(sprintf("║    95%% CI: [%6.2f, %6.2f]                              ║\n", primary_result$lower.CL, primary_result$upper.CL))
cat(sprintf("║    p-value: %s                                    ║\n", 
            ifelse(primary_result$p.value < 0.001, "<0.001", sprintf("%.4f", primary_result$p.value))))
cat("║                                                              ║\n")
cat(sprintf("║  Interpretation: Lofexidine reduces MHOWS by %.2f points  ║\n", abs(primary_result$estimate)))
cat("║  compared to placebo (negative = better outcome)             ║\n")
cat("║                                                              ║\n")
if (primary_result$estimate < 0) {
  cat("║  ✓ CORRECT: Negative estimate (Lofexidine improves outcome) ║\n")
} else {
  cat("║  ✗ WARNING: Positive estimate (check contrast direction!)   ║\n")
}
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

write.csv(as.data.frame(primary_result), "Primary_Average_Treatment_Effect.csv", row.names = FALSE)

# Table 4
table4_data <- bind_rows(
  daily_results %>% mutate(Timepoint = paste("Day", Day)) %>%
    select(Timepoint, estimate, lower.CL, upper.CL, p.value),
  primary_result %>% as.data.frame() %>%
    mutate(Timepoint = "Average Days 4-8") %>%
    select(Timepoint, estimate, lower.CL, upper.CL, p.value)
) %>%
  mutate(
    `Treatment Effect (95% CI)` = sprintf("%.2f (%.2f, %.2f)", estimate, lower.CL, upper.CL),
    `p-value` = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  ) %>%
  select(Timepoint, `Treatment Effect (95% CI)`, `p-value`)

write.csv(table4_data, "Table4_Primary_Results.csv", row.names = FALSE)

cat("✓ Primary analysis completed\n\n")

# ==============================================================================
# VERIFICATION: Check that primary and MI have same direction
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║              CONTRAST DIRECTION VERIFICATION                 ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Primary Analysis:                                       ║\n"))
cat(sprintf("║    Estimate: %7.2f                                      ║\n", primary_result$estimate))
cat(sprintf("║    Sign: %s                                              ║\n", 
            ifelse(primary_result$estimate < 0, "Negative (correct)", "Positive (CHECK!)")))
cat("║                                                              ║\n")
cat("║  Expected: NEGATIVE value (Lofexidine reduces MHOWS)         ║\n")
cat("║                                                              ║\n")
cat("║  Both primary and MI analyses should have:                   ║\n")
cat("║    • Same sign (both negative)                               ║\n")
cat("║    • Similar magnitude (within 20-30% is acceptable)         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")


# ==============================================================================
# 6. SENSITIVITY ANALYSES
# ==============================================================================

cat("=== STEP 6: SENSITIVITY ANALYSES ===\n")


sensitivity_results <- list()

# SA1: Alternative covariance structures
cat("Running SA1: Alternative covariance structures...\n")

model_cs <- lme(
  fixed = MHOWS ~ treat * time_factor,
  random = ~ 1 | PATIENT,
  correlation = corCompSymm(form = ~ 1 | PATIENT),
  data = analysis_long %>% filter(day >= 3),
  method = "REML",
  na.action = na.omit
)

model_ar1 <- lme(
  fixed = MHOWS ~ treat * time_factor,
  random = ~ 1 | PATIENT,
  correlation = corAR1(form = ~ time_numeric | PATIENT),
  data = analysis_long %>% filter(day >= 3),
  method = "REML",
  na.action = na.omit
)

extract_avg_effect <- function(model) {
  emm <- emmeans(model, ~ treat, at = list(time_factor = as.character(4:8)))
  contrast_result <- summary(contrast(emm, method = "pairwise"), infer = TRUE)
  return(contrast_result)
}

sensitivity_results$SA1_CS <- extract_avg_effect(model_cs)
sensitivity_results$SA1_AR1 <- extract_avg_effect(model_ar1)

cat("✓ SA1 completed\n")





# SA2: Comprehensive Multiple imputation 

# ==============================================================================
# COMPREHENSIVE MULTIPLE IMPUTATION ANALYSIS
# ==============================================================================

cat("=== RUNNING COMPREHENSIVE MULTIPLE IMPUTATION ===\n\n")

# ------------------------------------------------------------------------------
# 1. PREPARE DATA FOR IMPUTATION
# ------------------------------------------------------------------------------

# Create wide format with additional auxiliary variables
mi_data_wide <- analysis_long %>%
  select(PATIENT, treat, day, MHOWS, age, gender, baseline_MHOWS, 
         IV, days_use_30d, depression, anxiety, smoke) %>%
  distinct(PATIENT, day, .keep_all = TRUE) %>%
  pivot_wider(
    names_from = day,
    values_from = MHOWS,
    names_prefix = "MHOWS_Day",
    values_fn = mean
  ) %>%
  # Keep only one row per patient with patient-level variables
  group_by(PATIENT) %>%
  summarise(
    treat = first(treat),
    age = first(age),
    gender = first(gender),
    baseline_MHOWS = first(baseline_MHOWS),
    IV = first(IV),
    days_use_30d = first(days_use_30d),
    depression = first(depression),
    anxiety = first(anxiety),
    smoke = first(smoke),
    across(starts_with("MHOWS_Day"), ~first(.x))
  ) %>%
  ungroup()

# Check structure
cat("MI data structure:\n")
cat("  N patients:", nrow(mi_data_wide), "\n")
cat("  N variables:", ncol(mi_data_wide), "\n\n")

# Visualize missing data pattern BEFORE imputation
cat("Creating missing data visualizations...\n")

# 1. Missing data heatmap
png("Missing_Data_Heatmap.png", width = 12, height = 8, units = "in", res = 300)
aggr(mi_data_wide %>% select(starts_with("MHOWS_Day")), 
     col = c("lightblue", "red"),
     numbers = TRUE,
     sortVars = TRUE,
     labels = names(mi_data_wide %>% select(starts_with("MHOWS_Day"))),
     cex.axis = 0.7,
     gap = 3,
     ylab = c("Proportion Missing", "Missing Data Pattern"))
dev.off()

# 2. Missing data pattern matrix
png("Missing_Data_Pattern.png", width = 10, height = 8, units = "in", res = 300)
missing_pattern_plot <- md.pattern(
  mi_data_wide %>% select(starts_with("MHOWS_Day"), baseline_MHOWS, age, gender),
  rotate.names = TRUE
)
dev.off()

# Save pattern summary
write.csv(as.data.frame(missing_pattern_plot), 
          "Missing_Data_Pattern_Summary.csv", 
          row.names = TRUE)

cat("✓ Missing data visualizations created\n\n")

# ------------------------------------------------------------------------------
# 2. CONFIGURE IMPUTATION
# ------------------------------------------------------------------------------

# Set up predictor matrix (controls what predicts what)
pred_matrix <- make.predictorMatrix(mi_data_wide)

# Don't use PATIENT as predictor
pred_matrix[, "PATIENT"] <- 0

# Don't impute PATIENT or treatment (these are always observed)
pred_matrix["PATIENT", ] <- 0
pred_matrix["treat", ] <- 0

# Use all auxiliaries to predict MHOWS
mhows_cols <- grep("MHOWS_Day", colnames(mi_data_wide))
for (col in mhows_cols) {
  pred_matrix[col, ] <- 1
  pred_matrix[col, "PATIENT"] <- 0
  pred_matrix[col, "treat"] <- 1  # Treatment predicts outcomes
}

# Set up methods
imp_methods <- make.method(mi_data_wide)
imp_methods["PATIENT"] <- ""  # Don't impute
imp_methods["treat"] <- ""     # Don't impute

# Use PMM for continuous, logistic for binary
imp_methods[mhows_cols] <- "pmm"
if ("gender" %in% names(imp_methods)) imp_methods["gender"] <- "logreg"
if ("IV" %in% names(imp_methods)) imp_methods["IV"] <- "logreg"
if ("depression" %in% names(imp_methods)) imp_methods["depression"] <- "logreg"
if ("anxiety" %in% names(imp_methods)) imp_methods["anxiety"] <- "logreg"
if ("smoke" %in% names(imp_methods)) imp_methods["smoke"] <- "logreg"

cat("Imputation configuration:\n")
cat("  Method: Predictive Mean Matching (PMM) for MHOWS\n")
cat("  Number of imputations: 50\n")
cat("  Maximum iterations: 20\n")
cat("  Including auxiliary variables: Yes\n\n")

# ------------------------------------------------------------------------------
# 3. RUN MICE
# ------------------------------------------------------------------------------

cat("Running MICE imputation (this may take 2-3 minutes)...\n")

set.seed(12345)
imp <- mice(
  mi_data_wide,
  m = 50,                    # 50 imputations
  method = imp_methods,
  predictorMatrix = pred_matrix,
  maxit = 20,               # 20 iterations per imputation
  seed = 12345,
  printFlag = FALSE
)

cat("✓ MICE imputation completed\n\n")

# Check for logged events
if (nrow(imp$loggedEvents) > 0) {
  cat("Logged events during imputation:\n")
  print(imp$loggedEvents)
  write.csv(imp$loggedEvents, "MICE_Logged_Events.csv", row.names = FALSE)
  cat("\nNote: A few logged events is normal and usually not concerning\n\n")
} else {
  cat("✓ No problematic events during imputation\n\n")
}

# ------------------------------------------------------------------------------
# 4. CONVERGENCE DIAGNOSTICS
# ------------------------------------------------------------------------------

cat("Creating convergence diagnostics...\n")

# Overall convergence plot
png("MI_Convergence_Overall.png", width = 12, height = 10, units = "in", res = 300)
plot(imp, layout = c(2, 3))
dev.off()

# Detailed convergence for MHOWS variables only
mhows_vars <- colnames(mi_data_wide)[grep("MHOWS_Day", colnames(mi_data_wide))]

png("MI_Convergence_MHOWS.png", width = 14, height = 10, units = "in", res = 300)
plot(imp, c(mhows_vars[1:min(6, length(mhows_vars))]))
dev.off()

cat("✓ Convergence plots created\n")
cat("  Check 'MI_Convergence_Overall.png' - streams should mix well\n")
cat("  Look for: variance and mean stabilizing across iterations\n\n")

# ------------------------------------------------------------------------------
# 5. DISTRIBUTION DIAGNOSTICS
# ------------------------------------------------------------------------------

cat("Creating distribution diagnostics...\n")

# Strip plots for MHOWS variables (shows observed vs imputed)
png("MI_Stripplot_MHOWS.png", width = 14, height = 10, units = "in", res = 300)
stripplot(imp, MHOWS_Day4 + MHOWS_Day5 + MHOWS_Day6 ~ .imp, 
          pch = 20, cex = 1.2,
          main = "Distribution of Observed vs Imputed Values")
dev.off()

# Density plots with error handling
safe_densityplot <- function(imp_obj, vars) {
  tryCatch({
    # Only plot variables with enough non-missing data
    vars_to_plot <- character(0)
    for (var in vars) {
      n_observed <- sum(!is.na(mi_data_wide[[var]]))
      n_missing <- sum(is.na(mi_data_wide[[var]]))
      if (n_observed >= 3 & n_missing >= 2) {  # Need at least 3 observed and 2 missing
        vars_to_plot <- c(vars_to_plot, var)
      }
    }
    
    if (length(vars_to_plot) > 0) {
      # Create formula
      formula_str <- paste(vars_to_plot[1:min(6, length(vars_to_plot))], collapse = " + ")
      formula_obj <- as.formula(paste("~", formula_str))
      
      densityplot(imp_obj, formula_obj)
    } else {
      plot.new()
      text(0.5, 0.5, "Insufficient data for density plots", cex = 1.5)
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Error creating density plot:\n", e$message), cex = 1.2)
  })
}

png("MI_Densityplot_MHOWS.png", width = 14, height = 10, units = "in", res = 300)
safe_densityplot(imp, mhows_vars)
dev.off()

cat("✓ Distribution diagnostics created\n")
cat("  Check 'MI_Stripplot_MHOWS.png' - imputed should overlap observed\n")
cat("  Blue = observed, Red = imputed\n\n")

# ------------------------------------------------------------------------------
# 6. CORRELATION STRUCTURE EXAMINATION
# ------------------------------------------------------------------------------

cat("Examining correlation structure across imputations...\n")

# Calculate correlations for complete data and each imputation
observed_data <- mi_data_wide %>% 
  select(all_of(mhows_vars)) %>%
  na.omit()

if (nrow(observed_data) >= 3) {
  cor_observed <- cor(observed_data, use = "complete.obs")
  
  # Get correlations from imputed datasets
  cor_imputed_list <- lapply(1:min(10, imp$m), function(i) {
    imputed_complete <- complete(imp, i) %>%
      select(all_of(mhows_vars))
    cor(imputed_complete)
  })
  
  # Average correlation across imputations
  cor_imputed_avg <- Reduce("+", cor_imputed_list) / length(cor_imputed_list)
  
  # Create comparison heatmap
  library(reshape2)
  library(viridis)
  
  cor_observed_long <- melt(cor_observed)
  cor_imputed_long <- melt(cor_imputed_avg)
  
  p1 <- ggplot(cor_observed_long, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    scale_fill_viridis(limits = c(-1, 1)) +
    labs(title = "Correlation: Complete Cases Only",
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(cor_imputed_long, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    scale_fill_viridis(limits = c(-1, 1)) +
    labs(title = "Correlation: After Multiple Imputation (Average)",
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png("MI_Correlation_Comparison.png", width = 14, height = 6, units = "in", res = 300)
  print(gridExtra::grid.arrange(p1, p2, ncol = 2))
  dev.off()
  
  cat("✓ Correlation heatmaps created\n")
  cat("  Compare observed vs imputed correlation structure\n\n")
} else {
  cat("! Insufficient complete cases for correlation comparison\n\n")
}

# ------------------------------------------------------------------------------
# 7. FIT MODELS TO IMPUTED DATA
# ------------------------------------------------------------------------------

cat("Fitting linear mixed models to imputed datasets...\n")

# Convert to long format for modeling

imputed_long_list <- lapply(1:imp$m, function(i) {
  if (i %% 10 == 0) cat("  Processing imputation", i, "of", imp$m, "\n")
  
  complete(imp, i) %>%
    pivot_longer(
      cols = starts_with("MHOWS_Day"),
      names_to = "day",
      values_to = "MHOWS",
      names_prefix = "MHOWS_Day"
    ) %>%
    mutate(
      day = as.numeric(day),
      time_factor = factor(day),
      time_numeric = day - 3,
      .imp = i
    ) %>%
    filter(day >= 3 & day <= 8)
})

# Fit model to each imputation
mi_results_list <- lapply(seq_along(imputed_long_list), function(i) {
  if (i %% 10 == 0) cat("  Fitting model", i, "of", length(imputed_long_list), "\n")
  
  df <- imputed_long_list[[i]]
  
  tryCatch({
    # Fit the model
    model <- lme(
      fixed = MHOWS ~ treat * time_factor,
      random = ~ 1 + time_numeric | PATIENT,
      data = df,
      method = "REML",
      na.action = na.omit,
      control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
    )
    
    # Get treatment effects at each timepoint
    emm <- emmeans(model, ~ treat | time_factor, 
                   at = list(time_factor = as.character(4:8)))
    
    # Average across time (Days 4-8)
    emm_avg <- emmeans(model, ~ treat, 
                       at = list(time_factor = as.character(4:8)))
    
    # Check what emmeans is actually giving you
    emm_check <- emmeans(model, ~ treat, at = list(time_factor = "4"))
    contrast(emm_check, method = "pairwise")
    
    #contrast_avg <- summary(contrast(emm_avg, method = "pairwise"), infer = TRUE)
    #contrast_avg <- summary(contrast(emm_avg, method = "pairwise", reverse = TRUE), infer = TRUE)
    # Explicitly specify the contrast direction
    contrast_avg <- summary(contrast(emm_avg, method = "pairwise"), infer = TRUE)
    #contrast_avg <- summary(contrast(emm_avg, method = "pairwise", reverse = TRUE), infer = TRUE)
    #contrast_avg <- summary(contrast(emm_avg, method = list("Lofexidine - Placebo" = c(1, -1))), infer = TRUE)
    
    return(data.frame(
      imputation = i,
      estimate = contrast_avg$estimate,
      SE = contrast_avg$SE,
      df = contrast_avg$df,
      lower.CL = contrast_avg$lower.CL,
      upper.CL = contrast_avg$upper.CL
    ))
  }, error = function(e) {
    cat("    Error in imputation", i, ":", e$message, "\n")
    return(NULL)
  })
})

# Remove failed imputations
mi_results_list <- mi_results_list[!sapply(mi_results_list, is.null)]
mi_results_df <- do.call(rbind, mi_results_list)

cat("\n✓ Models fitted to", nrow(mi_results_df), "imputed datasets\n\n")

# Save individual imputation results
write.csv(mi_results_df, "MI_Individual_Imputation_Results.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. POOL RESULTS USING RUBIN'S RULES
# ------------------------------------------------------------------------------

cat("Pooling results using Rubin's Rules...\n")

if (nrow(mi_results_df) > 0) {
  mi_estimates <- mi_results_df$estimate
  mi_ses <- mi_results_df$SE
  
  # Rubin's rules
  m <- length(mi_estimates)
  Q_bar <- mean(mi_estimates)              # Pooled estimate
  U_bar <- mean(mi_ses^2)                  # Within-imputation variance
  B <- var(mi_estimates)                    # Between-imputation variance
  T <- U_bar + (1 + 1/m) * B               # Total variance
  SE_pooled <- sqrt(T)
  
  # Relative increase in variance due to nonresponse
  r <- (1 + 1/m) * B / U_bar
  
  # Fraction of missing information
  lambda <- (B + B/m) / T
  
  # Degrees of freedom (Barnard-Rubin adjustment)
  df_old <- (m - 1) / lambda^2
  
  # Observed data df (conservative estimate)
  n_total <- nrow(analysis_long)
  p <- 4  # number of parameters in interaction model
  df_obs <- n_total - p
  df_obs_adj <- df_obs * (1 - lambda)
  
  # Adjusted df
  df_adj <- (df_old * df_obs_adj) / (df_old + df_obs_adj)
  
  # Confidence interval
  t_crit <- qt(0.975, df_adj)
  lower_CL <- Q_bar - t_crit * SE_pooled
  upper_CL <- Q_bar + t_crit * SE_pooled
  
  # P-value
  t_stat <- Q_bar / SE_pooled
  p_value <- 2 * pt(abs(t_stat), df_adj, lower.tail = FALSE)
  
  # Store results
  sensitivity_results$SA2_MI <- data.frame(
    estimate = Q_bar,
    SE = SE_pooled,
    lower.CL = lower_CL,
    upper.CL = upper_CL,
    df = df_adj,
    p.value = p_value,
    FMI = lambda,
    RIV = r
  )
  
  # Create detailed summary
  mi_summary <- data.frame(
    Parameter = c("Pooled Estimate", "Standard Error", "95% CI Lower", "95% CI Upper",
                  "Degrees of Freedom", "t-statistic", "p-value",
                  "Between-imputation variance (B)", "Within-imputation variance (U)",
                  "Total variance (T)", "Relative increase in variance (r)",
                  "Fraction of missing information (λ)"),
    Value = c(Q_bar, SE_pooled, lower_CL, upper_CL, df_adj, t_stat, p_value,
              B, U_bar, T, r, lambda)
  )
  
  write.csv(mi_summary, "MI_Pooled_Results_Detailed.csv", row.names = FALSE)
  
  # Print summary
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║        MULTIPLE IMPUTATION RESULTS (Rubin's Rules)          ║\n")
  cat("╠══════════════════════════════════════════════════════════════╣\n")
  cat(sprintf("║  Pooled Estimate: %7.2f MHOWS points                    ║\n", Q_bar))
  cat(sprintf("║  Pooled SE:       %7.2f                                 ║\n", SE_pooled))
  cat(sprintf("║  95%% CI: [%6.2f, %6.2f]                              ║\n", lower_CL, upper_CL))
  cat(sprintf("║  p-value: %s                                    ║\n", 
              ifelse(p_value < 0.001, "<0.001", sprintf("%.4f", p_value))))
  cat("║                                                              ║\n")
  cat(sprintf("║  Degrees of freedom: %.1f                               ║\n", df_adj))
  cat(sprintf("║  Between-imputation variance (B): %.2f                  ║\n", B))
  cat(sprintf("║  Within-imputation variance (U): %.2f                   ║\n", U_bar))
  cat(sprintf("║  Relative increase in variance: %.2f%%                  ║\n", r * 100))
  cat(sprintf("║  Fraction missing information: %.2f%%                   ║\n", lambda * 100))
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")
  
  cat("Interpretation:\n")
  cat("  • FMI (λ) =", round(lambda * 100, 1), "% - proportion of variation due to missing data\n")
  cat("  • RIV (r) =", round(r, 2), "- variance inflated by", round(r * 100, 1), "% due to missingness\n")
  if (lambda < 0.10) {
    cat("  • Low FMI suggests MAR assumption is likely reasonable\n")
  } else if (lambda < 0.30) {
    cat("  • Moderate FMI suggests some uncertainty about missing data mechanism\n")
  } else {
    cat("  • High FMI suggests substantial impact of missing data - sensitivity analysis crucial\n")
  }
  cat("\n")
  
  cat("✓ SA2 (Multiple Imputation) completed successfully\n\n")
  
} else {
  cat("✗ SA2 (Multiple Imputation) failed - insufficient successful imputations\n\n")
}

# ------------------------------------------------------------------------------
# 9. VISUALIZE POOLED RESULTS
# ------------------------------------------------------------------------------

cat("Creating visualization of imputation variability...\n")

# Forest plot of individual imputation estimates
png("MI_Forest_Plot_Imputations.png", width = 10, height = 12, units = "in", res = 300)

mi_results_df %>%
  mutate(imputation = factor(imputation)) %>%
  ggplot(aes(x = estimate, y = imputation)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = Q_bar, linetype = "solid", color = "blue", size = 1) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), 
                 height = 0.3, alpha = 0.6, color = "darkgray") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = nrow(mi_results_df) + 1, linetype = "solid", size = 1) +
  annotate("point", x = Q_bar, y = nrow(mi_results_df) + 1, 
           size = 4, color = "blue", shape = 18) +
  annotate("errorbarh", xmin = lower_CL, xmax = upper_CL, 
           y = nrow(mi_results_df) + 1, height = 0.5, color = "blue", size = 1.2) +
  annotate("text", x = Q_bar, y = nrow(mi_results_df) + 2, 
           label = sprintf("Pooled: %.2f [%.2f, %.2f]", Q_bar, lower_CL, upper_CL),
           hjust = 0.5, color = "blue", fontface = "bold") +
  labs(title = "Treatment Effect Estimates Across Imputations",
       subtitle = "Individual imputation estimates (gray) and pooled result (blue)",
       x = "Treatment Effect (MHOWS points)",
       y = "Imputation Number") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6))

dev.off()

# Distribution of estimates
png("MI_Estimate_Distribution.png", width = 10, height = 6, units = "in", res = 300)

ggplot(mi_results_df, aes(x = estimate)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, 
                 fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1.2) +
  geom_vline(xintercept = Q_bar, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = c(lower_CL, upper_CL), 
             color = "red", linetype = "dotted", size = 0.8) +
  annotate("text", x = Q_bar, y = Inf, 
           label = sprintf("Pooled: %.2f", Q_bar),
           vjust = 2, color = "red", fontface = "bold") +
  labs(title = "Distribution of Treatment Effect Estimates Across Imputations",
       subtitle = sprintf("Mean = %.2f, SD = %.2f, Range = [%.2f, %.2f]",
                          mean(mi_results_df$estimate),
                          sd(mi_results_df$estimate),
                          min(mi_results_df$estimate),
                          max(mi_results_df$estimate)),
       x = "Treatment Effect Estimate (MHOWS points)",
       y = "Density") +
  theme_classic()

dev.off()

cat("✓ Imputation variability visualizations created\n\n")

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║           MULTIPLE IMPUTATION ANALYSIS COMPLETE              ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat("║  Files created:                                              ║\n")
cat("║    • Missing_Data_Heatmap.png                                ║\n")
cat("║    • Missing_Data_Pattern.png                                ║\n")
cat("║    • MI_Convergence_Overall.png (CHECK THIS!)                ║\n")
cat("║    • MI_Convergence_MHOWS.png                                ║\n")
cat("║    • MI_Stripplot_MHOWS.png                                  ║\n")
cat("║    • MI_Densityplot_MHOWS.png                                ║\n")
cat("║    • MI_Correlation_Comparison.png                           ║\n")
cat("║    • MI_Forest_Plot_Imputations.png                          ║\n")
cat("║    • MI_Estimate_Distribution.png                            ║\n")
cat("║    • MI_Individual_Imputation_Results.csv                    ║\n")
cat("║    • MI_Pooled_Results_Detailed.csv                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")


# Convergence diagnostics
png("MI_Convergence_Diagnostic.png", width = 10, height = 8, units = "in", res = 300)
plot(imp)
dev.off()

png("MI_Density_Diagnostic.png", width = 10, height = 8, units = "in", res = 300)
densityplot(imp)
dev.off()

# SA3: Complete case analysis
cat("Running SA3: Complete case analysis...\n")

complete_patients <- analysis_long %>%
  group_by(PATIENT) %>%
  summarise(n_obs = sum(!is.na(MHOWS))) %>%
  filter(n_obs == 6) %>%
  pull(PATIENT)

analysis_complete <- analysis_long %>%
  filter(PATIENT %in% complete_patients, day >= 3)

model_complete <- lme(
  fixed = MHOWS ~ treat * time_factor,
  random = ~ 1 + time_numeric | PATIENT,
  data = analysis_complete,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

sensitivity_results$SA3_Complete <- extract_avg_effect(model_complete)

primary_ci_width <- primary_result$upper.CL - primary_result$lower.CL
complete_ci_width <- sensitivity_results$SA3_Complete$upper.CL - sensitivity_results$SA3_Complete$lower.CL
precision_loss <- round(100 * (complete_ci_width - primary_ci_width) / primary_ci_width, 1)

cat("✓ SA3 completed (Precision loss:", precision_loss, "%)\n")

# SA4: Covariate adjustment
cat("Running SA4: Covariate adjustment...\n")

model_adjusted <- lme(
  fixed = MHOWS ~ treat * time_factor + baseline_MHOWS + age + gender + IV,
  random = ~ 1 + time_numeric | PATIENT,
  data = analysis_long %>% filter(day >= 3),
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

sensitivity_results$SA4_Adjusted <- extract_avg_effect(model_adjusted)

adjusted_ci_width <- sensitivity_results$SA4_Adjusted$upper.CL - sensitivity_results$SA4_Adjusted$lower.CL
precision_gain <- round(100 * (primary_ci_width - adjusted_ci_width) / primary_ci_width, 1)

cat("✓ SA4 completed (Precision gain:", precision_gain, "%)\n")





# Test if observed outcomes predict future dropout (MAR plausibility)
mar_test <- analysis_long %>%
  arrange(PATIENT, day) %>%
  group_by(PATIENT) %>%
  mutate(
    lagged_MHOWS = lag(MHOWS),
    future_missing = lead(is.na(MHOWS))
  ) %>%
  filter(!is.na(lagged_MHOWS), !is.na(future_missing))

mar_model <- glm(future_missing ~ treat + lagged_MHOWS + baseline_MHOWS + age, 
                 data = mar_test, family = binomial)
summary(mar_model)



#####Considering Partial imputation 

# ==============================================================================
# MULTIPLE IMPUTATION WITH PASSIVE DERIVATION
# ==============================================================================

# Prepare data with COMPONENTS, not total MHOWS
mi_data_components <- analysis_long %>%
  select(PATIENT, treat, day, 
         discrete_symptoms, continuous_signs, emesis_points,
         age, gender, baseline_MHOWS, IV, days_use_30d, 
         depression, anxiety, smoke) %>%
  distinct(PATIENT, day, .keep_all = TRUE) %>%
  pivot_wider(
    names_from = day,
    values_from = c(discrete_symptoms, continuous_signs, emesis_points),
    names_sep = "_Day"
  ) %>%
  group_by(PATIENT) %>%
  summarise(
    treat = first(treat),
    age = first(age),
    gender = first(gender),
    baseline_MHOWS = first(baseline_MHOWS),
    IV = first(IV),
    days_use_30d = first(days_use_30d),
    depression = first(depression),
    anxiety = first(anxiety),
    smoke = first(smoke),
    across(starts_with("discrete_symptoms"), ~first(.x)),
    across(starts_with("continuous_signs"), ~first(.x)),
    across(starts_with("emesis_points"), ~first(.x))
  ) %>%
  ungroup()

# Set up imputation methods
init <- mice(mi_data_components, maxit = 0)
meth <- init$method
pred <- init$predictorMatrix

# Don't impute ID or treatment
meth["PATIENT"] <- ""
meth["treat"] <- ""
pred[, "PATIENT"] <- 0
pred["PATIENT", ] <- 0
pred["treat", ] <- 0

# Impute components using PMM
component_vars <- c(
  grep("discrete_symptoms_Day", names(mi_data_components), value = TRUE),
  grep("continuous_signs_Day", names(mi_data_components), value = TRUE),
  grep("emesis_points_Day", names(mi_data_components), value = TRUE)
)

for (var in component_vars) {
  meth[var] <- "pmm"
}

# Run MICE on COMPONENTS
cat("Running MICE on component variables...\n")
imp_components <- mice(
  mi_data_components,
  m = 50,
  method = meth,
  predictorMatrix = pred,
  maxit = 20,
  seed = 12345,
  printFlag = FALSE
)

cat("✓ Component imputation completed\n\n")

# ==============================================================================
# Step 3: PASSIVELY DERIVE MHOWS from imputed components
# ==============================================================================

# For each imputed dataset, calculate MHOWS from components
imputed_long_list <- lapply(1:imp_components$m, function(i) {
  if (i %% 10 == 0) cat("  Processing imputation", i, "of", imp_components$m, "\n")
  
  # Get completed data with imputed components
  completed <- complete(imp_components, i)
  
  # Reshape to long format
  long_data <- completed %>%
    pivot_longer(
      cols = starts_with(c("discrete_symptoms_Day", "continuous_signs_Day", "emesis_points_Day")),
      names_to = c(".value", "day"),
      names_pattern = "(.+)_Day([0-9]+)"
    ) %>%
    mutate(
      day = as.numeric(day),
      # PASSIVELY DERIVE MHOWS from components
      MHOWS = discrete_symptoms + continuous_signs + emesis_points,
      time_factor = factor(day),
      time_numeric = day - 3,
      .imp = i
    ) %>%
    filter(day >= 3 & day <= 8)
  
  return(long_data)
})

cat("✓ MHOWS passively derived from imputed components\n\n")


# ==============================================================================
# SA6: FIT LINEAR MIXED MODELS TO PASSIVELY-DERIVED MHOWS
# ==============================================================================

cat("=== SA6: Multiple Imputation with Passive Derivation ===\n\n")

# Fit model to each imputation
cat("Fitting linear mixed models to passively-derived MHOWS...\n")

mi_passive_results <- lapply(seq_along(imputed_long_list), function(i) {
  if (i %% 10 == 0) cat("  Fitting model", i, "of", length(imputed_long_list), "\n")
  
  df <- imputed_long_list[[i]]
  
  tryCatch({
    # Fit the same model as primary analysis
    model <- lme(
      fixed = MHOWS ~ treat * time_factor,
      random = ~ 1 + time_numeric | PATIENT,
      data = df,
      method = "REML",
      na.action = na.omit,
      control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
    )
    
    # Get average treatment effect (Days 4-8)
    emm_avg <- emmeans(model, ~ treat, 
                       at = list(time_factor = as.character(4:8)))
    contrast_avg <- summary(contrast(emm_avg, method = "pairwise"), infer = TRUE)
    
    return(data.frame(
      imputation = i,
      estimate = contrast_avg$estimate,
      SE = contrast_avg$SE,
      df = contrast_avg$df,
      lower.CL = contrast_avg$lower.CL,
      upper.CL = contrast_avg$upper.CL
    ))
  }, error = function(e) {
    cat("    Error in imputation", i, ":", e$message, "\n")
    return(NULL)
  })
})

# Remove failed imputations
mi_passive_results <- mi_passive_results[!sapply(mi_passive_results, is.null)]
mi_passive_df <- do.call(rbind, mi_passive_results)

cat("\n✓ Models fitted to", nrow(mi_passive_df), "passively-derived datasets\n\n")

# Save individual results
write.csv(mi_passive_df, "MI_Passive_Individual_Results.csv", row.names = FALSE)


# ==============================================================================
# POOL RESULTS USING RUBIN'S RULES
# ==============================================================================

cat("Pooling results using Rubin's Rules...\n")

if (nrow(mi_passive_df) > 0) {
  mi_estimates <- mi_passive_df$estimate
  mi_ses <- mi_passive_df$SE
  
  # Rubin's rules
  m <- length(mi_estimates)
  Q_bar <- mean(mi_estimates)              # Pooled estimate
  U_bar <- mean(mi_ses^2)                  # Within-imputation variance
  B <- var(mi_estimates)                    # Between-imputation variance
  T <- U_bar + (1 + 1/m) * B               # Total variance
  SE_pooled <- sqrt(T)
  
  # Relative increase in variance due to nonresponse
  r <- (1 + 1/m) * B / U_bar
  
  # Fraction of missing information
  lambda <- (B + B/m) / T
  
  # Degrees of freedom (Barnard-Rubin adjustment)
  df_old <- (m - 1) / lambda^2
  n_total <- nrow(analysis_long)
  p <- 4
  df_obs <- n_total - p
  df_obs_adj <- df_obs * (1 - lambda)
  df_adj <- (df_old * df_obs_adj) / (df_old + df_obs_adj)
  
  # Confidence interval
  t_crit <- qt(0.975, df_adj)
  lower_CL <- Q_bar - t_crit * SE_pooled
  upper_CL <- Q_bar + t_crit * SE_pooled
  
  # P-value
  t_stat <- Q_bar / SE_pooled
  p_value <- 2 * pt(abs(t_stat), df_adj, lower.tail = FALSE)
  
  # Store results
  sensitivity_results$SA6_MI_Passive <- data.frame(
    estimate = Q_bar,
    SE = SE_pooled,
    lower.CL = lower_CL,
    upper.CL = upper_CL,
    df = df_adj,
    p.value = p_value,
    FMI = lambda,
    RIV = r
  )
  
  # Create detailed summary
  mi_passive_summary <- data.frame(
    Parameter = c("Pooled Estimate", "Standard Error", "95% CI Lower", "95% CI Upper",
                  "Degrees of Freedom", "t-statistic", "p-value",
                  "Between-imputation variance (B)", "Within-imputation variance (U)",
                  "Total variance (T)", "Relative increase in variance (r)",
                  "Fraction of missing information (λ)"),
    Value = c(Q_bar, SE_pooled, lower_CL, upper_CL, df_adj, t_stat, p_value,
              B, U_bar, T, r, lambda)
  )
  
  write.csv(mi_passive_summary, "MI_Passive_Pooled_Results.csv", row.names = FALSE)
  
  # Print summary
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║   SA6: MULTIPLE IMPUTATION (PASSIVE DERIVATION) RESULTS     ║\n")
  cat("╠══════════════════════════════════════════════════════════════╣\n")
  cat(sprintf("║  Pooled Estimate: %7.2f MHOWS points                    ║\n", Q_bar))
  cat(sprintf("║  Pooled SE:       %7.2f                                 ║\n", SE_pooled))
  cat(sprintf("║  95%% CI: [%6.2f, %6.2f]                              ║\n", lower_CL, upper_CL))
  cat(sprintf("║  p-value: %s                                    ║\n", 
              ifelse(p_value < 0.001, "<0.001", sprintf("%.4f", p_value))))
  cat("║                                                              ║\n")
  cat(sprintf("║  Degrees of freedom: %.1f                               ║\n", df_adj))
  cat(sprintf("║  Between-imputation variance (B): %.2f                  ║\n", B))
  cat(sprintf("║  Within-imputation variance (U): %.2f                   ║\n", U_bar))
  cat(sprintf("║  Relative increase in variance: %.2f%%                  ║\n", r * 100))
  cat(sprintf("║  Fraction missing information: %.2f%%                   ║\n", lambda * 100))
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")
  
  cat("✓ SA6 (Multiple Imputation - Passive Derivation) completed successfully\n\n")
  
} else {
  cat("✗ SA6 failed - insufficient successful imputations\n\n")
}


# ==============================================================================
# COMPARE DIRECT VS PASSIVE IMPUTATION
# ==============================================================================

cat("=== COMPARING MI APPROACHES ===\n\n")

# Extract results from both MI approaches
mi_direct_est <- sensitivity_results$SA2_MI$estimate
mi_direct_lower <- sensitivity_results$SA2_MI$lower.CL
mi_direct_upper <- sensitivity_results$SA2_MI$upper.CL
mi_direct_fmi <- sensitivity_results$SA2_MI$FMI

mi_passive_est <- sensitivity_results$SA6_MI_Passive$estimate
mi_passive_lower <- sensitivity_results$SA6_MI_Passive$lower.CL
mi_passive_upper <- sensitivity_results$SA6_MI_Passive$upper.CL
mi_passive_fmi <- sensitivity_results$SA6_MI_Passive$FMI

# Create comparison table
mi_comparison <- data.frame(
  Approach = c("MI: Direct Imputation (SA2)", 
               "MI: Passive Derivation (SA6)",
               "Difference"),
  Estimate = c(mi_direct_est, mi_passive_est, mi_passive_est - mi_direct_est),
  CI_Lower = c(mi_direct_lower, mi_passive_lower, NA),
  CI_Upper = c(mi_direct_upper, mi_passive_upper, NA),
  CI_Width = c(mi_direct_upper - mi_direct_lower, 
               mi_passive_upper - mi_passive_lower, NA),
  FMI = c(mi_direct_fmi * 100, mi_passive_fmi * 100, NA),
  p_value = c(sensitivity_results$SA2_MI$p.value,
              sensitivity_results$SA6_MI_Passive$p.value, NA)
)

print(mi_comparison)
write.csv(mi_comparison, "MI_Comparison_Direct_vs_Passive.csv", row.names = FALSE)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║           MI APPROACH COMPARISON                             ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Direct Imputation:    %.2f [%.2f, %.2f]              ║\n", 
            mi_direct_est, mi_direct_lower, mi_direct_upper))
cat(sprintf("║  Passive Derivation:   %.2f [%.2f, %.2f]              ║\n", 
            mi_passive_est, mi_passive_lower, mi_passive_upper))
cat(sprintf("║  Difference:           %.2f points                        ║\n", 
            mi_passive_est - mi_direct_est))
cat(sprintf("║  FMI Direct:           %.1f%%                              ║\n", 
            mi_direct_fmi * 100))
cat(sprintf("║  FMI Passive:          %.1f%%                              ║\n", 
            mi_passive_fmi * 100))
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Interpretation
if (abs(mi_passive_est - mi_direct_est) < 1) {
  cat("Interpretation: Minimal difference between approaches - results robust\n")
} else if (abs(mi_passive_est - mi_direct_est) < 2) {
  cat("Interpretation: Small difference between approaches - generally consistent\n")
} else {
  cat("Interpretation: Notable difference between approaches - warrants discussion\n")
}

# ==============================================================================
# VISUALIZE MI COMPARISON
# ==============================================================================

cat("\nCreating MI comparison visualizations...\n")

# Side-by-side forest plot
mi_forest_data <- data.frame(
  Method = factor(c("Direct Imputation", "Passive Derivation"),
                  levels = c("Passive Derivation", "Direct Imputation")),
  Estimate = c(mi_direct_est, mi_passive_est),
  Lower = c(mi_direct_lower, mi_passive_lower),
  Upper = c(mi_direct_upper, mi_passive_upper)
)

fig_mi_compare <- ggplot(mi_forest_data, aes(x = Estimate, y = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = -11.84, linetype = "dotted", color = "blue", alpha = 0.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, size = 1.2, color = "#9B59B6") +
  geom_point(size = 4, color = "#9B59B6") +
  labs(
    title = "Multiple Imputation Comparison",
    subtitle = "Direct Imputation vs Passive Derivation of MHOWS",
    x = "Treatment Effect (MHOWS points)",
    y = "",
    caption = "Blue dotted line = Primary LMM estimate (-11.84)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave("Figure_MI_Comparison.png", fig_mi_compare, width = 8, height = 4, dpi = 300)

# Distribution comparison
mi_dist_data <- rbind(
  data.frame(Method = "Direct", Estimate = mi_results_df$estimate),
  data.frame(Method = "Passive", Estimate = mi_passive_df$estimate)
)

fig_mi_dist <- ggplot(mi_dist_data, aes(x = Estimate, fill = Method)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = mi_direct_est, color = "#E91E63", linetype = "dashed") +
  geom_vline(xintercept = mi_passive_est, color = "#9C27B0", linetype = "dashed") +
  scale_fill_manual(values = c("Direct" = "#E91E63", "Passive" = "#9C27B0")) +
  labs(
    title = "Distribution of Estimates Across Imputations",
    subtitle = "Comparison of Direct vs Passive Derivation",
    x = "Treatment Effect Estimate",
    y = "Density",
    fill = "MI Method"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("Figure_MI_Distribution_Comparison.png", fig_mi_dist, width = 8, height = 5, dpi = 300)

cat("✓ Visualizations created\n\n")

sa6_row <- data.frame(
  Analysis = "SA6: MI Passive Derivation",
  Estimate = mi_passive_est,
  CI_Lower = mi_passive_lower,
  CI_Upper = mi_passive_upper,
  CI_Width = mi_passive_upper - mi_passive_lower,
  Change_from_Primary = mi_passive_est - (-11.84),
  Percent_Change = round(100 * (mi_passive_est - (-11.84)) / (-11.84), 1),
  p_value = sensitivity_results$SA6_MI_Passive$p.value
)

# Format for display
sa6_display <- sprintf(
  "MI: Passive Derivation | %.2f (%.2f, %.2f) | %.2f | %.2f (%s%%) | %.3f",
  sa6_row$Estimate,
  sa6_row$CI_Lower,
  sa6_row$CI_Upper,
  sa6_row$CI_Width,
  sa6_row$Change_from_Primary,
  ifelse(sa6_row$Change_from_Primary > 0, "+", ""),
  sa6_row$Percent_Change,
  sa6_row$p_value
)

cat("\nAdd this row to Table 6 (under Alternative MAR Approaches):\n")
cat(sa6_display, "\n\n")




# ==============================================================================
# SA5: MNAR Tipping Point Analysis
# ==============================================================================

cat("Running SA5: MNAR tipping point analysis...\n")

# Need zoo package for na.locf
if (!require("zoo", quietly = TRUE)) {
  install.packages("zoo")
  library(zoo)
}

mhows_sd <- sd(analysis_long$MHOWS, na.rm = TRUE)
delta_values <- seq(0, 1.5, by = 0.25)
tipping_results <- data.frame()

for (delta_sd in delta_values) {
  delta_points <- delta_sd * mhows_sd
  
  data_adjusted <- analysis_long %>%
    arrange(PATIENT, day) %>%
    group_by(PATIENT) %>%
    mutate(
      is_missing = is.na(MHOWS),
      MHOWS_imputed = zoo::na.locf(MHOWS, na.rm = FALSE),
      MHOWS_mnar = ifelse(is_missing, MHOWS_imputed + delta_points, MHOWS)
    ) %>%
    ungroup() %>%
    filter(day >= 3)
  
  tryCatch({
    model_mnar <- lme(
      fixed = MHOWS_mnar ~ treat * time_factor,
      random = ~ 1 + time_numeric | PATIENT,
      data = data_adjusted,
      method = "REML",
      na.action = na.omit,
      control = lmeControl(opt = "optim", maxIter = 50)
    )
    
    emm_mnar <- emmeans(model_mnar, ~ treat, at = list(time_factor = as.character(4:8)))
    contrast_mnar <- summary(contrast(emm_mnar, method = "pairwise"), infer = TRUE)
    
    tipping_results <- rbind(
      tipping_results,
      data.frame(
        delta_sd = delta_sd,
        delta_points = delta_points,
        estimate = contrast_mnar$estimate,
        SE = contrast_mnar$SE,
        lower.CL = contrast_mnar$lower.CL,
        upper.CL = contrast_mnar$upper.CL,
        p.value = contrast_mnar$p.value,
        significant = contrast_mnar$p.value < 0.05,
        ci_includes_zero = (contrast_mnar$lower.CL < 0 & contrast_mnar$upper.CL > 0)
      )
    )
  }, error = function(e) {
    cat("  Error at delta =", delta_sd, "\n")
  })
}

tipping_point <- tipping_results %>%
  filter(ci_includes_zero) %>%
  slice(1)

if (nrow(tipping_point) > 0) {
  cat("✓ SA5 completed (Tipping point:", tipping_point$delta_sd, "SD)\n")
} else {
  cat("✓ SA5 completed (No tipping point within 0-1.5 SD)\n")
}

write.csv(tipping_results, "Tipping_Point_Results.csv", row.names = FALSE)

sensitivity_results$SA5_MNAR_0.25SD <- tipping_results %>% filter(delta_sd == 0.25)
sensitivity_results$SA5_MNAR_0.50SD <- tipping_results %>% filter(delta_sd == 0.50)
sensitivity_results$SA5_MNAR_1.00SD <- tipping_results %>% filter(delta_sd == 1.00)

cat("\n✓ All sensitivity analyses completed\n\n")

# ==============================================================================
# 7. CREATE SUMMARY TABLE (TABLE 5)
# ==============================================================================

cat("=== STEP 7: SUMMARY TABLE ===\n")

compile_sa_result <- function(result, name) {
  if (is.data.frame(result) && nrow(result) > 0 && !is.na(result$estimate[1])) {
    data.frame(
      Analysis = name,
      Estimate = result$estimate[1],
      Lower_CL = result$lower.CL[1],
      Upper_CL = result$upper.CL[1],
      p_value = result$p.value[1]
    )
  } else {
    NULL
  }
}

table5 <- bind_rows(
  data.frame(
    Analysis = "Primary: LMM (all available data)",
    Estimate = primary_result$estimate,
    Lower_CL = primary_result$lower.CL,
    Upper_CL = primary_result$upper.CL,
    p_value = primary_result$p.value
  ),
  compile_sa_result(sensitivity_results$SA1_CS, "Alternative cov: Compound Symmetry"),
  compile_sa_result(sensitivity_results$SA1_AR1, "Alternative cov: AR(1)"),
  compile_sa_result(sensitivity_results$SA3_Complete, "Complete Case Analysis"),
  compile_sa_result(sensitivity_results$SA4_Adjusted, "Covariate Adjusted"),
  compile_sa_result(sensitivity_results$SA5_MNAR_0.25SD, "MNAR: δ=0.25 SD"),
  compile_sa_result(sensitivity_results$SA5_MNAR_0.50SD, "MNAR: δ=0.50 SD"),
  compile_sa_result(sensitivity_results$SA5_MNAR_1.00SD, "MNAR: δ=1.00 SD")
) %>%
  mutate(
    CI_Width = Upper_CL - Lower_CL,
    Estimate_CI = sprintf("%.2f (%.2f, %.2f)", Estimate, Lower_CL, Upper_CL),
    p_formatted = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  ) %>%
  select(Analysis, Estimate_CI, CI_Width, p_formatted)

cat("\n--- Table 5: Sensitivity Analysis Summary ---\n")
print(table5)
write.csv(table5, "Table5_Sensitivity_Analysis_Summary.csv", row.names = FALSE)

cat("\n✓ Summary table created\n\n")

# ==============================================================================
# 8. FOREST PLOT AND TIPPING POINT VISUALIZATION
# ==============================================================================

cat("=== STEP 8: FINAL VISUALIZATIONS ===\n")

# Prepare forest plot data
forest_data <- data.frame()

# Add primary
forest_data <- rbind(forest_data, 
                     data.frame(Analysis = "Primary", Estimate = primary_result$estimate,
                                Lower = primary_result$lower.CL, Upper = primary_result$upper.CL, 
                                Type = "Primary"))

# Add other results if they exist
if (!is.na(sensitivity_results$SA1_CS$estimate)) {
  forest_data <- rbind(forest_data,
                       data.frame(Analysis = "CS", Estimate = sensitivity_results$SA1_CS$estimate,
                                  Lower = sensitivity_results$SA1_CS$lower.CL, 
                                  Upper = sensitivity_results$SA1_CS$upper.CL, Type = "Covariance"))
}

if (!is.na(sensitivity_results$SA3_Complete$estimate)) {
  forest_data <- rbind(forest_data,
                       data.frame(Analysis = "Complete Case", 
                                  Estimate = sensitivity_results$SA3_Complete$estimate,
                                  Lower = sensitivity_results$SA3_Complete$lower.CL, 
                                  Upper = sensitivity_results$SA3_Complete$upper.CL, 
                                  Type = "Efficiency"))
}

if (nrow(tipping_results) > 0) {
  forest_data <- rbind(forest_data,
                       data.frame(Analysis = "MNAR δ=0.5", 
                                  Estimate = tipping_results$estimate[tipping_results$delta_sd == 0.5],
                                  Lower = tipping_results$lower.CL[tipping_results$delta_sd == 0.5], 
                                  Upper = tipping_results$upper.CL[tipping_results$delta_sd == 0.5], 
                                  Type = "MNAR"))
}

forest_data <- forest_data %>%
  mutate(Analysis = factor(Analysis, levels = rev(Analysis)))

fig4 <- ggplot(forest_data, aes(x = Estimate, y = Analysis)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = primary_result$estimate, 
             linetype = "dotted", color = "blue", alpha = 0.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, color = Type), 
                 height = 0.2, size = 1) +
  geom_point(aes(color = Type), size = 3) +
  scale_color_manual(values = c("Primary" = "#2C3E50", 
                                "Covariance" = "#3498DB", 
                                "Efficiency" = "#27AE60",
                                "MNAR" = "#E74C3C")) +
  labs(title = "Forest Plot: Sensitivity Analyses",
       subtitle = "Average Treatment Effect (Days 4-8): Lofexidine - Placebo",
       x = "Treatment Effect (MHOWS points)", 
       y = "", 
       color = "Analysis Type") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_line(color = "gray90"))

ggsave("Figure4_Forest_Plot_Sensitivity.png", fig4, width = 10, height = 6, dpi = 300)

# Tipping point curve
if (nrow(tipping_results) > 0) {
  fig5 <- ggplot(tipping_results, aes(x = delta_sd, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.3, fill = "blue") +
    geom_line(size = 1.2, color = "darkblue") +
    geom_point(size = 3, color = "darkblue") +
    labs(title = "Tipping Point Analysis: MNAR Sensitivity",
         subtitle = "Treatment effect as a function of MNAR departure (δ)",
         x = "δ (SD units)", 
         y = "Treatment Effect Estimate\n(with 95% CI)",
         caption = "Shaded region = 95% CI | Red line = null effect") +
    theme_classic(base_size = 12) +
    theme(panel.grid.major = element_line(color = "gray90"))
  
  if (exists("tipping_point") && nrow(tipping_point) > 0) {
    fig5 <- fig5 +
      geom_vline(xintercept = tipping_point$delta_sd, 
                 linetype = "dotted", color = "red", size = 1) +
      annotate("text", x = tipping_point$delta_sd, 
               y = max(tipping_results$upper.CL),
               label = sprintf("Tipping point\nδ = %.2f SD", tipping_point$delta_sd),
               hjust = -0.1, color = "red")
  }
  
  ggsave("Figure5_Tipping_Point.png", fig5, width = 10, height = 6, dpi = 300)
}

cat("✓ Visualizations created\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("                    ANALYSIS COMPLETE!                                          \n")
cat("================================================================================\n\n")

cat("PRIMARY RESULTS:\n")
cat("  Treatment Effect:", round(primary_result$estimate, 2), "points\n")
cat("  95% CI: [", round(primary_result$lower.CL, 2), ",", 
    round(primary_result$upper.CL, 2), "]\n")
cat("  p-value:", ifelse(primary_result$p.value < 0.001, "<0.001", 
                         sprintf("%.4f", primary_result$p.value)), "\n\n")

cat("SENSITIVITY ANALYSIS SUMMARY:\n")
cat("  Complete case precision loss:", precision_loss, "%\n")
if (nrow(tipping_point) > 0) {
  cat("  MNAR tipping point:", tipping_point$delta_sd, "SD\n")
} else {
  cat("  MNAR tipping point: Not reached (highly robust)\n")
}

cat("\nKEY FILES GENERATED:\n")
cat("  • Table1_Baseline_Characteristics.html\n")
cat("  • Table2_Missing_Data_Summary.csv\n")
cat("  • Table3_Unadjusted_Means.csv\n")
cat("  • Table4_Primary_Results.csv\n")
cat("  • Table5_Sensitivity_Analysis_Summary.csv\n")
cat("  • Figure2_Individual_Trajectories.png\n")
cat("  • Figure3_Mean_Trajectories.png\n")
cat("  • Figure4_Forest_Plot_Sensitivity.png\n")
if (nrow(tipping_results) > 0) {
  cat("  • Figure5_Tipping_Point.png\n")
}

cat("\n================================================================================\n")

#===============================================================================
# APPENDIX B.8: MODEL DIAGNOSTICS FOR PRIMARY ANALYSIS
#===============================================================================

cat("=== MODEL DIAGNOSTICS FOR PRIMARY LINEAR MIXED MODEL ===\n")

#------------------------------------------------------------------------------
# SUCCESSFUL DIAGNOSTIC PLOTS
#------------------------------------------------------------------------------

cat("→ Creating comprehensive diagnostic plots...\n")

# Extract model data and residuals (this worked)
model_data <- getData(primary_model)
resid_df <- data.frame(
  Fitted = fitted(primary_model),
  Residuals = resid(primary_model, type = "pearson"),
  Patient = model_data$PATIENT,
  Day = model_data$day,
  Treatment = model_data$treat
)

# Merge with original data for additional variables
analysis_long_diagnostics <- analysis_long %>%
  inner_join(resid_df, by = c("PATIENT" = "Patient", "day" = "Day"))

cat("✓ Diagnostic data prepared successfully\n")

#------------------------------------------------------------------------------
# 1. RESIDUAL DIAGNOSTIC PLOTS (ALL WORKED)
#------------------------------------------------------------------------------

cat("→ Creating residual diagnostic plots...\n")

# 1.1 Residual vs Fitted Plot
p1 <- ggplot(resid_df, aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6, color = "#3498DB") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", color = "#E74C3C", se = TRUE) +
  labs(title = "Residuals vs Fitted Values",
       subtitle = "Primary Linear Mixed Model",
       x = "Fitted Values", y = "Pearson Residuals") +
  theme_classic()

# 1.2 QQ Plot of Residuals
p2 <- ggplot(resid_df, aes(sample = Residuals)) +
  stat_qq(alpha = 0.6, color = "#3498DB") +
  stat_qq_line(color = "#E74C3C") +
  labs(title = "Q-Q Plot of Residuals",
       subtitle = "Checking normality assumption",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_classic()

# 1.3 Histogram of Residuals
p3 <- ggplot(resid_df, aes(x = Residuals)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, fill = "#3498DB", alpha = 0.7) +
  geom_density(color = "#E74C3C", linewidth = 1) +
  labs(title = "Distribution of Residuals",
       subtitle = "With density overlay",
       x = "Pearson Residuals", y = "Density") +
  theme_classic()

# 1.4 Residuals by Treatment Group
p4 <- ggplot(resid_df, aes(x = Treatment, y = Residuals, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("PLACEBO" = "#E74C3C", "LOFEXIDINE" = "#3498DB")) +
  labs(title = "Residuals by Treatment Group",
       subtitle = "Checking homoscedasticity",
       x = "Treatment", y = "Pearson Residuals") +
  theme_classic() +
  theme(legend.position = "none")

# Arrange and save residual plots
residual_grid <- grid.arrange(p1, p2, p3, p4, ncol = 2)
ggsave("Figure6_Residual_Diagnostics.png", residual_grid, width = 12, height = 10, dpi = 300)

cat("✓ Residual diagnostic plots created and saved\n")

#------------------------------------------------------------------------------
# 2. RANDOM EFFECTS DIAGNOSTICS (ALL WORKED)
#------------------------------------------------------------------------------

cat("→ Creating random effects diagnostic plots...\n")

# Extract random effects
ranef_df <- ranef(primary_model)
colnames(ranef_df) <- c("Random_Intercept", "Random_Slope")

# 2.1 Random intercepts vs slopes
p5 <- ggplot(ranef_df, aes(x = Random_Intercept, y = Random_Slope)) +
  geom_point(alpha = 0.6, color = "#9B59B6") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", color = "#E74C3C", se = TRUE) +
  labs(title = "Random Intercepts vs Random Slopes",
       subtitle = "BLUPs for each patient",
       x = "Random Intercept", y = "Random Slope") +
  theme_classic()

# 2.2 QQ plots for random effects
p6 <- ggplot(ranef_df, aes(sample = Random_Intercept)) +
  stat_qq(alpha = 0.6, color = "#3498DB") +
  stat_qq_line(color = "#E74C3C") +
  labs(title = "Q-Q Plot: Random Intercepts",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_classic()

p7 <- ggplot(ranef_df, aes(sample = Random_Slope)) +
  stat_qq(alpha = 0.6, color = "#3498DB") +
  stat_qq_line(color = "#E74C3C") +
  labs(title = "Q-Q Plot: Random Slopes", 
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_classic()

# Arrange and save random effects plots
ranef_grid <- grid.arrange(p5, p6, p7, ncol = 2)
ggsave("Figure7_Random_Effects_Diagnostics.png", ranef_grid, width = 12, height = 8, dpi = 300)

cat("✓ Random effects diagnostic plots created and saved\n")

#------------------------------------------------------------------------------
# 3. AUTOCORRELATION AND TIME TREND DIAGNOSTICS (WORKED)
#------------------------------------------------------------------------------

cat("→ Creating autocorrelation and time trend plots...\n")

# 3.1 Autocorrelation function
acf_data <- acf(resid(primary_model), plot = FALSE)
acf_df <- with(acf_data, data.frame(lag, acf))

p8 <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(xend = lag, yend = 0), size = 1, color = "#3498DB") +
  geom_hline(yintercept = c(-1, 1) * 0.05, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Autocorrelation Function of Residuals",
       subtitle = "Checking independence assumption",
       x = "Lag", y = "Autocorrelation") +
  theme_classic()

# 3.2 Residuals over time by treatment
p9 <- ggplot(analysis_long_diagnostics, aes(x = day, y = Residuals, group = PATIENT)) +
  geom_line(alpha = 0.2, color = "gray50") +
  geom_smooth(aes(group = 1), method = "loess", color = "#E74C3C", se = TRUE) +
  facet_wrap(~treat) +
  labs(title = "Residuals Over Time by Treatment",
       subtitle = "Checking for time trends",
       x = "Study Day", y = "Pearson Residuals") +
  theme_classic()

# Arrange and save final diagnostics
final_diag_grid <- grid.arrange(p8, p9, ncol = 2)
ggsave("Figure8_Final_Diagnostics.png", final_diag_grid, width = 12, height = 6, dpi = 300)

cat("✓ Autocorrelation and time trend plots created and saved\n")

#------------------------------------------------------------------------------
# 4. BASIC MODEL ASSUMPTIONS CHECK (SIMPLIFIED)
#------------------------------------------------------------------------------

cat("→ Performing basic model assumptions checks...\n")

# Simple checks that work without complex data manipulation
assumptions_checks <- data.frame(
  Check = c(
    "Normality of Residuals (Shapiro-Wilk)",
    "Mean of Residuals ≈ 0", 
    "Autocorrelation at Lag 1",
    "Random Effects Normality (Intercepts)",
    "Random Effects Normality (Slopes)"
  ),
  Result = c(
    ifelse(shapiro.test(resid(primary_model))$p.value > 0.05, "PASS", "CAUTION"),
    ifelse(abs(mean(resid(primary_model))) < 0.001, "PASS", "CHECK"),
    ifelse(abs(acf_data$acf[2]) < 0.2, "PASS", "CAUTION"),
    ifelse(shapiro.test(ranef_df$Random_Intercept)$p.value > 0.05, "PASS", "CAUTION"), 
    ifelse(shapiro.test(ranef_df$Random_Slope)$p.value > 0.05, "PASS", "CAUTION")
  ),
  Details = c(
    sprintf("W = %.4f, p = %.3f", 
            shapiro.test(resid(primary_model))$statistic,
            shapiro.test(resid(primary_model))$p.value),
    sprintf("Mean = %.6f", mean(resid(primary_model))),
    sprintf("ACF(1) = %.4f", acf_data$acf[2]),
    sprintf("W = %.4f, p = %.3f", 
            shapiro.test(ranef_df$Random_Intercept)$statistic,
            shapiro.test(ranef_df$Random_Intercept)$p.value),
    sprintf("W = %.4f, p = %.3f", 
            shapiro.test(ranef_df$Random_Slope)$statistic,
            shapiro.test(ranef_df$Random_Slope)$p.value)
  )
)

print(assumptions_checks)
write.csv(assumptions_checks, "Model_Assumptions_Summary.csv", row.names = FALSE)

cat("✓ Basic model assumptions checks completed\n")

#------------------------------------------------------------------------------
# 5. FINAL DIAGNOSTICS SUMMARY
#------------------------------------------------------------------------------

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║               MODEL DIAGNOSTICS SUMMARY                     ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat("║  SUCCESSFULLY COMPLETED DIAGNOSTICS:                        ║\n")
cat("║                                                              ║\n")
cat("║  • Residuals vs Fitted values plot                          ║\n")
cat("║  • Q-Q plot of residuals                                    ║\n")
cat("║  • Residual distribution histogram                          ║\n")
cat("║  • Residuals by treatment group                             ║\n")
cat("║  • Random intercepts vs slopes                              ║\n")
cat("║  • Q-Q plots for random effects                             ║\n")
cat("║  • Autocorrelation function plot                            ║\n")
cat("║  • Residuals over time by treatment                         ║\n")
cat("║                                                              ║\n")
cat("║  FILES GENERATED:                                           ║\n")
cat("║    • Figure6_Residual_Diagnostics.png                       ║\n")
cat("║    • Figure7_Random_Effects_Diagnostics.png                 ║\n")
cat("║    • Figure8_Final_Diagnostics.png                          ║\n")
cat("║    • Model_Assumptions_Summary.csv                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("✓ Comprehensive model diagnostics completed successfully\n")
cat("✓ All diagnostic plots saved for manuscript inclusion\n\n")

# Save updated workspace
save.image("lofexidine_analysis_complete.RData")
cat("✓ Updated workspace saved: lofexidine_analysis_complete.RData\n")

