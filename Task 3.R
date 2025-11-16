# ==============================================================================

# CM2018 HT25-1 Statistics for Medical Engineering 7.5 credits
# Task 3 Treatment Effect in Large-Stage Pancreatic Cancer Patients


# === PACKAGES & LIBRARIES =====================================================

# Install Packages (install once if needed)
#install.packages(c("readr","dplyr","ggplot2", "corrplot", "survminer"))
library(readr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(tidyr)
library(survival)
library(survminer)


# === SETUP ====================================================================

# NOTE: Adjust path as needed.

# Load data (long format)
data_3 <- read_csv("Data_T3.csv", show_col_types = FALSE)

# Drop first (index) column
data_3 <- data_3[,-1]

# Get numerical and categorical data
data_3_num <- data_3[, c(1,3,6,7,8)]
data_3_cat <- data_3[, c(2,4,5)]

# Convert categorical variables to numeric (0/1)
data_3_mod <- data_3
data_3_mod$TreatmentGroup <- ifelse(data_3$TreatmentGroup == "Treatment", 1, 0)
data_3_mod$Sex <- ifelse(data_3$Sex == "Female", 1, 0)


# === MISSING DATA =============================================================

# Column-wise count of missing values
colSums(is.na(data_3))

# Proportion of missing per variable
colMeans(is.na(data_3))

# === Basic plots ==============================================================

# --- Histogram plots for continuous variables (Age, GFR, Time) ----------------

data_3 %>%
  pivot_longer(cols = c(Age, GFR, Time),
               names_to = "variable",
               values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, colour = "black") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw()

# --- Bar plots for categorical variables (TreatmentGroup, Sex, ECOG_PS) -------

data_3_mod %>%
  pivot_longer(cols = c(TreatmentGroup, Sex, ECOG_PS),
               names_to = "variable",
               values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_bar() +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw()

# --- Boxplots (Age & GFR by TreatmentGroup) -----------------------------------

# Set 1 row, 2 columns
par(mfrow = c(1, 2))

# Age by treatment group
boxplot(Age ~ TreatmentGroup, 
        data = data_3_mod,
        col = c("lightblue", "lightgreen"),
        xlab = "Treatment Group",
        ylab = "Age (years)",
        main = "Age Distribution by Treatment Group")

# GFR by treatment group
boxplot(GFR ~ TreatmentGroup, 
        data = data_3_mod,
        col = c("lightblue", "lightgreen"),
        xlab = "Treatment Group",
        ylab = "GFR (mL/min/1.73m²)",
        main = "GFR Distribution by Treatment Group")



# Compare baseline covariates by TreatmentGroup (simple plots)
# Age and GFR by treatment
data_3 %>%
  ggplot(aes(x = TreatmentGroup, y = Age)) +
  geom_boxplot() +
  theme_bw()

data_3 %>%
  ggplot(aes(x = TreatmentGroup, y = GFR)) +
  geom_boxplot() +
  theme_bw()

# ECOG distribution by treatment
data_3 %>%
  ggplot(aes(x = ECOG_PS, fill = TreatmentGroup)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  theme_bw()


# === Covariance and Correlation matrices ======================================

# The covariance matrix contains the covariances between all pairs of variables.
# Each entry tells you how two variables vary together in absolute terms, 
# in the original units of the data.

# The diagonal values show large differences in variance across variables, 
# reflecting their different measurement scales. ID has an especially large 
# variance because it is an identifier rather than a true variable. Off-diagonal 
# covariances are small, confirming weak co-variation between variables. This 
# means the variables vary mostly independently, though their scales differ 
# substantially.

# The correlation matrix is a standardized version of the covariance matrix:
# This re-scales everything so that:
#  * Diagonal entries = 1
#  * Off-diagonal entries = correlations between -1 and +1
# So it only tells you about direction and strength of linear relationships, 
#not their magnitude in the data’s units.

# The variables show very weak correlations, mostly between –0.2 and +0.2, 
# indicating that they are largely independent. This suggests there is no 
# multicollinearity among the independent variables, which is good for modeling. 
# The only mildly positive relationships are between GFR and Time (r ≈ 0.22) 
# and between ECOG_PS and Event (r ≈ 0.20), but these are still weak.

# Set layout: 1 row, 2 columns
par(mfrow = c(1, 2))

# === Covariance matrix plot ===================================================
cov_mat <- cov(data_3_num)

corrplot(cov_mat,
         is.corr    = FALSE,
         method      = "circle",
         type        = "upper",
         col         = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = NA,
         tl.col      = "black",
         tl.srt      = 45,
         title       = "Covariance matrix of independent variables",
         mar         = c(0,0,2,0))

# === Correlation matrix plot ==================================================
cor_mat <- cor(data_3_num, use = "pairwise.complete.obs")

corrplot(cor_mat,
         method        = "circle",
         type          = "upper",
         diag          = TRUE,
         order         = "original",
         col           = colorRampPalette(c("red","white","blue"))(200),
         addCoef.col   = NA,
         tl.col        = "red",
         tl.srt        = 45,
         tl.cex        = 1.2,
         addgrid.col   = "grey90",
         cl.ratio      = 0.2,
         mar           = c(0,0,2,0),
         title         = "Correlation of independent variables")



# === Relationship between continuous variables  ===============================

# --- Age vs GFR ---------------------------------------------------------------

# The plot shows a very weak positive relationship between age and GFR in both 
# groups. The treatment and control lines almost overlap, indicating no 
# meaningful difference between groups and only a minimal age effect on GFR.
data_3$TreatmentGroup <- factor(data_3$TreatmentGroup)

# Linear
ggplot(data_3, aes(x = Age, y = GFR, color = TreatmentGroup)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Relationship between Age and GFR by Treatment Group")

# Nonlinear
ggplot(data_3, aes(x = Age, y = GFR, colour = TreatmentGroup)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title  = "Relationship Between Age and GFR by Treatment Group",
    x      = "Age (years)",
    y      = "GFR (mL/min/1.73m²)",
    color  = "Treatment group"
  ) +
  theme_bw()

# --- Time vs GFR --------------------------------------------------------------

# The plot shows that GFR increases slightly with time in both groups. The 
# treatment group has a slightly steeper upward trend, but the two groups 
# overlap heavily, suggesting no meaningful difference.<
data_3$TreatmentGroup <- factor(data_3$TreatmentGroup)

# Linear
ggplot(data_3, aes(x = Time, y = GFR, color = TreatmentGroup)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Linear GFR over Time by Treatment Group")

# Nonlinear
ggplot(data_3, aes(x = Time, y = GFR, colour = TreatmentGroup)) +
  geom_point(alpha = 0.6,
             position = position_jitter(width = 0, height = 0.5)) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title  = "Nonlinear GFR over Time by Treatment Group",
    x      = "Time (weeks)",
    y      = "GFR (mL/min/1.73m²)",
    color  = "Treatment group"
  ) +
  theme_bw()


# --- ECOG Performance Status over Time ----------------------------------------

# Make sure TreatmentGroup is a factor
data_3$TreatmentGroup <- factor(data_3$TreatmentGroup)

ggplot(data_3, aes(x = Time, y = ECOG_PS, colour = TreatmentGroup)) +
  geom_point(alpha = 0.6, 
             position = position_jitter(width = 0, height = 0.05)) +  # spread points a bit
  geom_smooth(method = "loess", se = TRUE) +                          # trend per group
  scale_y_continuous(breaks = 0:3) +                                  # ECOG is discrete
  labs(
    title = "ECOG Performance Status over Time",
    x     = "Time (weeks)",
    y     = "ECOG PS",
    colour = "Treatment group"
  ) +
  theme_bw()



# === Kaplan–Meier curves & log-rank test ======================================


# --- Cumulative number of deceased over time ----------------------------------

# Survival object
surv_object <- Surv(time = data_3$Time, event = data_3$Event)

# Fit KM curves by TreatmentGroup
fit_KM <- survfit(surv_object ~ TreatmentGroup, data = data_3)

p_KM <- ggsurvplot(
  fit_KM,
  data        = data_3,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = TRUE,
  legend.title = "Group",
  legend.labs  = levels(data_3$TreatmentGroup),
  title       = "Survival and Cumulative Deaths by Treatment Group",
  xlab        = "Time (weeks)",
  ylab        = "Survival probability",
  ggtheme     = theme_bw()
)

# --- By treatment group, compute cumulative number of deaths over time. -------
deaths_cum <- data_3 %>%
  filter(Event == 1) %>%                 # only those who died
  arrange(Time) %>%                      # order by time
  group_by(TreatmentGroup) %>%           # per treatment group
  mutate(CumDeaths = row_number())       # cumulative count within group

# Overall maximum number of deaths (for scaling)
max_deaths_all <- max(deaths_cum$CumDeaths)

# Scale cumulative deaths to [0, 1] so they can be overlaid on survival axis
deaths_cum <- deaths_cum %>%
  mutate(CumDeaths_scaled = CumDeaths / max_deaths_all)

# Add cumulative deaths as dashed steps to the KM plot, with a secondary y-axis.
p_KM$plot <- p_KM$plot +
  geom_step(
    data = deaths_cum,
    aes(x = Time, y = CumDeaths_scaled, color = TreatmentGroup),
    linetype = "dashed"
  ) +
  scale_y_continuous(
    name = "Survival probability",
    sec.axis = sec_axis(
      ~ . * max_deaths_all,
      name = "Cumulative number of deaths"
    )
  )

# Explicit log-rank test (same p-value as above)
survdiff(surv_object ~ TreatmentGroup, data = data_3)


# === Cox proportional hazards models ==========================================

cox_crude <- coxph(Surv(Time, Event) ~ TreatmentGroup, data = data_3)
summary(cox_crude)

# Extract hazard ratio and 95% CI
exp(cbind(HR = coef(cox_crude), confint(cox_crude)))

# Adjusted model: Treatment + covariates
cox_adj <- coxph(
  Surv(Time, Event) ~ TreatmentGroup + Age + Sex + ECOG_PS + GFR,
  data = data_3
)

summary(cox_adj)

# Hazard ratios with 95% CI
exp(cbind(HR = coef(cox_adj), confint(cox_adj)))


# === Statistical tests ========================================================

# Statistical test for treatment effect on GFR e.g., independent t-test or 
# ANCOVA (adjusting for Age, ECOG_PS, etc.)

# There was no significant difference in GFR between patients receiving 
# treatment and those in the control group, t(295.87) = –0.44, p = 0.662.
# Mean GFR values were 78.4 mL/min/1.73 m² for the control group and 
# 79.4 mL/min/1.73 m² for the treatment group, suggesting that the treatment 
# did not affect kidney function.
t.test(GFR ~ TreatmentGroup, data = data_3)

# There was no significant difference in GFR between male and female patients, 
# t(237.47) = 0.02, p = 0.98. The mean GFR was 78.9 mL/min/1.73 m² for females 
# and 78.9 mL/min/1.73 m² for males. This indicates that kidney function, as 
# measured by GFR, did not differ by sex in this cohort of pancreatic cancer 
# patients.
t.test(GFR ~ Sex, data = data_3)