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
data_3_num <- data_3[sapply(data_3, is.numeric)]
data_3_cat <- data_3[!sapply(data_3, is.numeric)]

# Convert categorical variables to numeric (0/1)
data_3$TreatmentGroup <- ifelse(data_3$TreatmentGroup == "Treatment", 1, 0)
data_3$Sex <- ifelse(data_3$Sex == "Female", 1, 0)


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

data_3 %>%
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
        data = data_3,
        col = c("lightblue", "lightgreen"),
        xlab = "Treatment Group",
        ylab = "Age (years)",
        main = "Age Distribution by Treatment Group")

# GFR by treatment group
boxplot(GFR ~ TreatmentGroup, 
        data = data_3,
        col = c("lightblue", "lightgreen"),
        xlab = "Treatment Group",
        ylab = "GFR (mL/min/1.73m²)",
        main = "GFR Distribution by Treatment Group")


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

# === Covariance matrix plot =======================================
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

# === Correlation matrix plot =======================================
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





# --- Relationship between continuous variables ------------------------------------------------

# The plot shows a very weak positive relationship between age and GFR in both 
# groups. The treatment and control lines almost overlap, indicating no 
# meaningful difference between groups and only a minimal age effect on GFR.
ggplot(data_3, aes(x = Age, y = GFR, color = TreatmentGroup)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Relationship between Age and GFR by Treatment Group")


# The plot shows that GFR increases slightly with time in both groups. The 
# treatment group has a slightly steeper upward trend, but the two groups 
# overlap heavily, suggesting no meaningful difference.<
ggplot(data_3, aes(x = Time, y = GFR, color = TreatmentGroup)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Relationship between Time and GFR by Treatment Group")


# --- Visualize survival time by treatment and event ---------------------------

surv_object <- Surv(time = data_3$Time, event = data_3$Event)
fit <- survfit(surv_object ~ TreatmentGroup, data = data_3)

ggsurvplot(fit, data = data_3,
           pval = TRUE,
           risk.table = TRUE,
           title = "Survival Curves by Treatment Group",
           xlab = "Time (weeks)", ylab = "Survival Probability")


# --- Statistical tests --------------------------------------------------------

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


# --- Multiple linear regression (GFR ~ predictors) ----------------------------
# To quantify how treatment, age, and performance status jointly affect kidney function.
lm_gfr <- lm(GFR ~ TreatmentGroup + Age + ECOG_PS, data = data_3)
summary(lm_gfr)





# --- Multiple linear regression: GFR as outcome ------------------------------

# 1) Prep & fit an MLR

#install.packages(("ggfortify", "ggeffects")) # if needed
library(broom)       # tidy(), augment(), glance()
library(car)         # avPlots(), vif()
library(ggeffects)   # ggpredict() for marginal effects
library(ggfortify)   # autoplot() diagnostics

data_3 <- data_3 %>%
  mutate(
    TreatmentGroup = factor(TreatmentGroup),
    ECOG_PS = factor(ECOG_PS)   # treat as factor; drop this line to keep numeric
  )

# Fit model (interaction Age*TreatmentGroup)
lm_gfr <- lm(GFR ~ Age * TreatmentGroup + ECOG_PS, data = data_3)

summary(lm_gfr)
vif(lm_gfr)  # quick multicollinearity check


# 2) Coefficient plot (with 95% CI)

coef_df <- broom::tidy(lm_gfr, conf.int = TRUE)

ggplot(coef_df %>% filter(term != "(Intercept)"),
       aes(x = reorder(term, estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Multiple Linear Regression on GFR",
       x = "", y = "Coefficient (with 95% CI)")


# 4) Actual vs. fitted (goodness-of-fit view)

aug <- broom::augment(lm_gfr)

ggplot(aug, aes(.fitted, GFR, color = TreatmentGroup)) +
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme_minimal() +
  labs(title = "Observed vs Fitted GFR", x = "Fitted GFR", y = "Observed GFR")



# Signed correlations including binary 0/1 vars
corr_vars <- data_3_mod %>% dplyr::select(where(is.numeric))  # includes your 0/1 binaries
cor_mat_signed <- cor(corr_vars, use = "pairwise.complete.obs")

corrplot(cor_mat_signed,
         method = "circle",
         type   = "upper",
         col    = colorRampPalette(c("red","white","blue"))(200),
         tl.col = "black",
         tl.srt = 45,
         title  = "Pearson correlations (numeric + binary 0/1)")



install.packages("DescTools") # if needed
library(DescTools)
library(dplyr)

df <- data_3  # use the original with categoricals intact
num_vars <- names(df)[sapply(df, is.numeric)]
cat_vars <- setdiff(names(df), num_vars)

# Helper: correlation ratio (eta^2) for numeric ~ categorical
eta2_num_cat <- function(x, g) {
  ok <- complete.cases(x, g)
  x <- x[ok]; g <- as.factor(g[ok])
  grand_mean <- mean(x)
  ss_total <- sum((x - grand_mean)^2)
  ss_between <- sum(tapply(x, g, function(v) length(v) * (mean(v) - grand_mean)^2))
  if (ss_total == 0) return(NA_real_)
  ss_between / ss_total  # eta^2 in [0,1]
}

vars <- names(df)
M <- matrix(NA_real_, nrow = length(vars), ncol = length(vars),
            dimnames = list(vars, vars))

for (i in seq_along(vars)) {
  for (j in seq_along(vars)) {
    vi <- df[[vars[i]]]
    vj <- df[[vars[j]]]
    
    if (i == j) {
      M[i, j] <- 1
      next
    }
    
    ni <- is.numeric(vi)
    nj <- is.numeric(vj)
    
    if (ni && nj) {
      # numeric-numeric: absolute Pearson (strength)
      M[i, j] <- abs(suppressWarnings(cor(vi, vj, use = "pairwise.complete.obs")))
    } else if (ni && !nj) {
      # numeric ~ categorical: eta^2
      M[i, j] <- eta2_num_cat(vi, vj)
    } else if (!ni && nj) {
      # numeric ~ categorical (symmetric)
      M[i, j] <- eta2_num_cat(vj, vi)
    } else {
      # categorical-categorical: Cramér's V
      ok <- complete.cases(vi, vj)
      tab <- table(vi[ok], vj[ok])
      if (all(dim(tab) >= 2)) {
        M[i, j] <- CramerV(tab, bias.correct = TRUE)
      } else {
        M[i, j] <- NA_real_
      }
    }
  }
}

# Visualize association strengths (0..1)
corrplot(M,
         is.corr = TRUE,            # treat as a 0..1 "correlation-like" measure
         method  = "circle",
         type    = "upper",
         col     = colorRampPalette(c("white","steelblue"))(200),
         tl.col  = "black",
         tl.srt  = 45,
         title   = "Association strengths (|r|, η², Cramér’s V)")




