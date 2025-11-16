## =========================================================
## Seminar 3 — Task 4 (Adverse Effects & Survival)
## =========================================================
## Structure:
## 1) Setup & data
## 2) EDA (boxplots, proportions, correlation matrix)
## 3) Survival exploration (KM by neutropenia)
## 4) Modelling
##    4.1 Logistic model: neutropenia
##    4.2 Cox model: survival
##    4.3 Logistic model: mortality
## =========================================================

## --- 1) Setup & data -------------------------------------
setwd("Your_working_directory")
options(stringsAsFactors = FALSE)

library(ggplot2)
library(GGally)
library(pscl)
library(pROC)
library(broom)
library(effects)
library(dplyr)
library(survival)
library(survminer)
library(ResourceSelection)
library(scales)

# Load data
data <- read.csv("merged_data.csv")

# Basic checks
print(dim(data))
print(summary(data))
print(table(data$Neutropenia))

# Helper theme & factors
theme_pub <- theme_bw(base_size = 12) + theme(panel.grid = element_blank())
data$Neutro_f <- factor(data$Neutropenia, levels = c(0,1), labels = c("No","Yes"))
data$Sex_f    <- factor(data$Sex)

## =========================================================
## 2) EDA
## =========================================================

# Age by neutropenia
p_age <- ggplot(data, aes(x = Neutro_f, y = Age)) +
  geom_boxplot(fill = "white") +
  labs(x = "Neutropenia", y = "Age (years)", title = "Age by Neutropenia") +
  theme_pub
print(p_age)
print(wilcox.test(Age ~ Neutropenia, data = data))

# GFR by neutropenia
p_gfr <- ggplot(data, aes(x = Neutro_f, y = GFR)) +
  geom_boxplot(fill = "white") +
  labs(x = "Neutropenia", y = "GFR (mL/min/1.73m²)", title = "GFR by Neutropenia") +
  theme_pub
print(p_gfr)
print(wilcox.test(GFR ~ Neutropenia, data = data))

# Time by neutropenia
p_time <- ggplot(data, aes(x = Neutro_f, y = Time)) +
  geom_boxplot(fill = "white") +
  labs(x = "Neutropenia", y = "Time to event (weeks)",
       title = "Time to Event by Neutropenia Status") +
  theme_pub
print(p_time)
print(wilcox.test(Time ~ Neutropenia, data = data))

# Sex vs neutropenia
p_sex <- ggplot(data, aes(x = Sex_f, fill = Neutro_f)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent) +
  labs(x = "Sex", y = "Proportion", fill = "Neutropenia",
       title = "Neutropenia by Sex") +
  theme_pub
print(p_sex)

# ECOG vs neutropenia
p_ecog <- ggplot(data, aes(x = factor(ECOG_PS), fill = Neutro_f)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent) +
  labs(x = "ECOG PS", y = "Proportion", fill = "Neutropenia",
       title = "Neutropenia by ECOG PS") +
  theme_pub
print(p_ecog)
print(chisq.test(table(data$Neutropenia, data$ECOG_PS)))

# Mortality event vs neutropenia
p_event <- ggplot(data, aes(x = factor(Event), fill = Neutro_f)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent) +
  labs(x = "Event", y = "Proportion", fill = "Neutropenia",
       title = "Neutropenia by Mortality") +
  theme_pub
print(p_event)
print(chisq.test(table(data$Neutropenia, data$Event)))

## ---- Correlation matrix ---------------------------------

num_vars <- na.omit(data[, c("Age","ECOG_PS","GFR","Time","Event","Neutropenia")])

upper_cor_heat <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  r <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
  ggplot() +
    geom_tile(aes(x = 1, y = 1, fill = r)) +
    scale_fill_gradient2(limits = c(-1,1),
                         low = "firebrick2", mid = "white", high = "steelblue3") +
    geom_text(aes(x = 1, y = 1, label = sprintf("%.2f", r)), size = 4) +
    theme_void() + theme(legend.position = "none")
}

p_pairs <- ggpairs(
  num_vars,
  lower = list(continuous = GGally::wrap("points", alpha = 0.35, size = 0.5)),
  diag  = list(continuous = "densityDiag"),
  upper = list(continuous = upper_cor_heat)
) + theme_pub
print(p_pairs)

## =========================================================
## 3) Survival exploration (Kaplan–Meier curves)
## =========================================================

# Convert to survival-friendly dataset
# (keep Neutropenia numeric for other models, but create a factor label for KM/Cox)
df_surv <- data %>%
  mutate(
    Event       = as.numeric(Event),         # 1 = death, 0 = censored
    Time        = as.numeric(Time),          # follow-up time in weeks
    Neutro_group = factor(Neutropenia,       # survival label
                          levels = c(0,1),
                          labels = c("No neutropenia","Neutropenia"))
  )

# Fit Kaplan–Meier curve by neutropenia status
km_fit <- survfit(Surv(Time, Event) ~ Neutro_group, data = df_surv)

# Plot KM curve with censoring marks, CI bands, and risk table
# → This shows “raw” survival differences before adjusting for confounders.
ggsurvplot(
  km_fit,
  data = df_surv,
  pval = TRUE,                # log-rank test
  conf.int = TRUE,            # confidence bands
  risk.table = TRUE,          # table of patients at risk over time
  censor = TRUE,              # show censoring ticks
  censor.shape = 124,
  censor.size = 3,
  legend.title = "Neutropenia status",
  legend.labs = c("No neutropenia", "Neutropenia"),
  surv.median.line = "hv",    # vertical/horizontal median survival lines
  xlab = "Time (weeks)",
  ylab = "Survival probability",
  palette = c("#1f78b4", "#e31a1c"),
  ggtheme = theme_minimal(base_size = 14)
)

# Extract median survival by group
surv_median(km_fit)
summary(km_fit)$table


## =========================================================
## 4) MODELLING SECTION
## =========================================================


## ---------------------------------------------------------
## 4.1 Logistic regression — Neutropenia model (toxicity)
## ---------------------------------------------------------
# Purpose:
#   Identify which baseline characteristics predict neutropenia.
#   This helps understand *why* neutropenia occurs.

fit_neut <- glm(
  Neutropenia ~ ECOG_PS + GFR + Age + Sex,
  data = data,
  family = binomial
)

summary(fit_neut)

# Odds ratios (easier to interpret than log-odds)
OR_neut <- exp(cbind(OR = coef(fit_neut), confint(fit_neut)))
print(OR_neut)

# Pseudo-R² shows overall model “fit quality”
pR2(fit_neut)

# ROC curve and AUC (discrimination performance)
roc_neut <- roc(response = data$Neutropenia, predictor = fitted(fit_neut))
plot(roc_neut, col = "darkblue", lwd = 2,
     main = "ROC Curve — Neutropenia Model")
abline(a = 0, b = 1, lty = 2, col = "gray")
auc(roc_neut)   # expected ~0.7

# Predicted probability vs GFR (visualizing GFR as the main driver)
data$pred_neut <- fitted(fit_neut)

p_pred_gfr <- ggplot(data, aes(GFR, pred_neut, color = Neutro_f)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 1) +
  labs(x = "GFR", y = "Predicted probability of neutropenia",
       title = "Predicted Probability vs GFR") +
  scale_color_manual(values = c("steelblue","firebrick")) +
  theme_pub
print(p_pred_gfr)

# Forest plot of adjusted odds ratios
fit_neut_tidy <- tidy(fit_neut, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = recode(term,
                       "Age" = "Age (per year)",
                       "SexMale" = "Male (vs Female)",
                       "ECOG_PS" = "ECOG PS (per unit)",
                       "GFR" = "GFR (per mL/min)"))

ggplot(fit_neut_tidy, aes(x = reorder(term, estimate), y = estimate)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  geom_point(size = 3, color = "firebrick") +
  coord_flip() +
  labs(y = "Odds Ratio (95% CI)",
       title = "Factors Associated with Neutropenia") +
  theme_pub


## ---------------------------------------------------------
## 4.2 Cox proportional hazards model — Survival
## ---------------------------------------------------------
# Purpose:
#   Determine whether neutropenia independently predicts survival,
#   after adjusting for baseline factors (GFR, ECOG, age, sex).

df_cox <- df_surv %>%
  mutate(
    Sex     = factor(Sex),
    ECOG_PS = as.factor(ECOG_PS)
  )

cox_fit <- coxph(
  Surv(Time, Event) ~ Neutro_group + GFR + ECOG_PS + Age + Sex,
  data = df_cox
)

summary(cox_fit)      # hazard ratios
cox.zph(cox_fit)      # proportional hazards assumption check
ggforest(cox_fit, data = df_cox)   # forest plot of HRs

# Adjusted survival curves
# These show survival for each neutropenia group *holding all other covariates constant*.
adj_neut_plot <- ggadjustedcurves(
  fit      = cox_fit,
  data     = df_cox,
  variable = "Neutro_group",
  legend.title = "Neutropenia (adjusted)",
  legend.labs  = c("No neutropenia", "Neutropenia")
)

adj_neut_plot +
  labs(x = "Time (weeks)", y = "Adjusted survival probability") +
  theme_minimal(base_size = 14)

# Predicted survival for hypothetical patients with different GFR
# → Shows GFR is a major survival determinant.
gfr_values <- c(40, 70, 100)

new_patients <- data.frame(
  Neutro_group = factor("No neutropenia", levels = levels(df_cox$Neutro_group)),
  GFR          = gfr_values,
  ECOG_PS      = factor("1", levels = levels(df_cox$ECOG_PS)),
  Age          = median(df_cox$Age),
  Sex          = factor("Male", levels = levels(df_cox$Sex))
)

sf_gfr <- survfit(cox_fit, newdata = new_patients)

ggsurvplot(
  sf_gfr,
  data = new_patients,
  conf.int = TRUE,
  legend.title = "GFR (mL/min/1.73m²)",
  legend.labs  = paste("GFR", gfr_values),
  xlab = "Time (weeks)",
  ylab = "Predicted survival probability",
  ggtheme = theme_minimal(base_size = 14)
)


## ---------------------------------------------------------
## 4.3 Logistic regression — Mortality model (binary outcome)
## ---------------------------------------------------------
# Purpose:
#   Predict death (Event = 1) as a binary outcome.
#   Complements Cox model by asking:
#     “Does neutropenia help predict mortality *beyond* baseline factors?”

log_fit <- glm(
  Event ~ Neutro_group + GFR + ECOG_PS + Age + Sex,
  data = df_cox,
  family = binomial
)

summary(log_fit)

# Odds ratios (easy to read)
exp(cbind(OR = coef(log_fit), confint(log_fit)))

# Model fit
pR2(log_fit)

# ROC curve & AUC
roc_mort <- roc(df_cox$Event, fitted(log_fit))
plot(roc_mort, col = "blue", lwd = 3,
     main = "ROC Curve — Mortality Model")
abline(a = 0, b = 1, lty = 2)
auc(roc_mort)
