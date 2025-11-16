library(visdat)
library(ggplot2)
library(GGally)
library(corrplot)

data <- read.csv("data/Data_T1.csv")

# Handle region, age and sex as factors instead of categorical
data$Region <- factor(data$Region)
data$Sex <- factor(data$Sex, levels = c("Male", "Female"))
data$AgeGroup <- factor(data$AgeGroup, levels = c("20-39", "40-59", "60-79"))

# Creates a new variable, incidence rate, since all subgroups (region-sex-age) don't have identical Npopulations (although they are similar)
data$IncidenceRate <- (data$NewCases / data$Npopulation) * 100000

# Overview of data
str(data)
summary(data)

# Explore missingness
vis_miss(data) # No missing data


# ----------- General visualisations -------------------------------------------
# Incidence rate by region
ggplot(data, aes(x = Region, y = IncidenceRate, fill = Region)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Incidence rate across regions")

# Incidence rate by age group and sex, per region
ggplot(data, aes(x = AgeGroup, y = IncidenceRate, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ Region) + 
  theme_minimal() +
  ggtitle("Incidence across age group and sex, across regions")


# ------------ Lifestyle factors & cancer ----------------------------------------------
# Relationship between standardised CLI and cancer incidence for each region
ggplot(data, aes(x = CLIstd, y = IncidenceRate, color = Region)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  ggtitle("Incidence rate and CLI across regions")
# Healthier lifestyle = lower incidence rate for all regions

# Explore whether CLI affects cancer risk differently by age or sex
ggplot(data, aes(x = CLIstd, y = IncidenceRate, color = AgeGroup)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Sex) +
  theme_minimal() +
  ggtitle("Incidence rate and CLI across age group and sex")
# Healthier lifestyle = lower incidence rate for both sexes for all age groups except females 60-79

# Relationship between smoking prevalence and cancer incidence, higher values mean healthier lifestyle
ggplot(data, aes(x = SmokingPrevalence, y = IncidenceRate, color = Region)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()+
  ggtitle("Incidence rate and smoking prevalence across regions")
# Shows higher incidence rate for smokers in region 2, 3 & 5, lower incidence rate for smokers in region 1 & 4
  
# Explore whether smoking affects cancer risk differently by age or sex
ggplot(data, aes(x = SmokingPrevalence, y = IncidenceRate, color = AgeGroup)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Sex) +
  theme_minimal()+
  ggtitle("Incidence rate and smoking prevalence across age group and sex")
# Smoking increases incidence rates for both sexes for all age groups except females 60-79

# Relationship between median BMI and cancer incidence
ggplot(data, aes(x = BMImedian, y = IncidenceRate, color = Region)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
theme_minimal()+
  ggtitle("Incidence rate and median BMI across regions")
# Higher median BMI drastically increases incidence rate in region 2 & 4, slightly increases rate in region 5 
# and decreases rate in region 1 & 3

# Explore whether BMI affects cancer risk differently by age or sex
ggplot(data, aes(x = BMImedian, y = IncidenceRate, color = AgeGroup)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Sex) +
  theme_minimal()+
  ggtitle("Incidence rate and medican BMI across age group and sex")
# Higher median BMI drastically increases incidence rate for both sexes 40-59, females 60-79 and males 20-39,
# slightly decreases rate for females 20-39 and drastically decreases rate for males 60-79


# ----------- Correlation ------------------------------------------------------
# Correlation plot of numeric variables
cor_data <- cor(numeric_data)
corrplot(cor_data, addCoef.col = "black")
# CLIstd is moderately negatively correlated to incidence rate, meaning healthier lifestyle = lower risk for cancer
# Smoking and BMI is slightly positively correlated to incidence rate, meaning smoking and higher BMI = higher risk for cancer


###########################################################################################################################################

# TASK 1 - PART 2 ____________________________________________________________________

# Load additional required packages
library(MASS)
library(lme4)
library(DHARMa)
library(performance)
library(effects)


# __________________ MAKING THE MODEL _______________________________________________
#####################################################################################

# Since we're modeling counts with different population sizes at risk,
# we should use an offset for population size to handle differences 
data$logPopulation <- log(data$Npopulation)

# Initial Poisson model (also since we have count data)
poisson_model <- glm(NewCases ~ Region + AgeGroup + Sex + CLIstd + 
                       SmokingPrevalence + BMImedian + offset(logPopulation),
                     family = poisson(link = "log"),
                     data = data)

summary(poisson_model)

# Check for overdispersion
dispersion_test <- sum(residuals(poisson_model, type = "pearson")^2) / 
  df.residual(poisson_model)
dispersion_test

# RESULT : Dispersion = 1.039291 which is <1,5 so its okay
# (If we were to have a dispersion we could use Negative Binomial) 


# Model diagnostics using DHARMa
sim_res <- simulateResiduals(poisson_model)
plot(sim_res)

# Check indices of model performance 
model_performance <- performance::model_performance(poisson_model)
print(model_performance)


#_________________ MODEL INTERPRETATION _____________________________________________
#####################################################################################

# Exponentiate coefficients to get Incidence Rate Ratios (IRR)
# We do this since Poission regression model has the coeff. in log-scale 
irr <- exp(coef(poisson_model))
irr_ci <- exp(confint(poisson_model))

results_table <- data.frame(
  Variable = names(irr),
  IRR = round(irr, 3),
  Lower_CI = round(irr_ci[,1], 3),
  Upper_CI = round(irr_ci[,2], 3)
)

print(results_table)

### Where IRR>1 indicates increased risk and IRR<1 indicates protective effect 


#________________REGIONAL AND DEMOGRAPHICS ANALYSIS__________________________________
#####################################################################################

# Regional comparisons
regional_effects <- results_table[grep("Region", results_table$Variable),]
print("Regional Effects (compared to Region1):")
print(regional_effects)

# Age group effects
age_effects <- results_table[grep("AgeGroup", results_table$Variable),]
print("Age Group Effects (compared to 20-39):")
print(age_effects)

# Sex effect
sex_effect <- results_table[grep("Sex", results_table$Variable),]
print("Sex Effect (compared to Male):")
print(sex_effect)


#______________________LIFESTYLE FACTOR ANALYSIS____________________________________
####################################################################################

# Lifestyle factors impact
lifestyle_effects <- results_table[c("CLIstd", "SmokingPrevalence", "BMImedian"),]
print("Lifestyle Factor Effects:")
print(lifestyle_effects)

# Interpret lifestyle factors
cat("\n=== LIFESTYLE FACTOR INTERPRETATION ===\n")
cat("CLIstd (Composite Lifestyle Index):", 
    ifelse(lifestyle_effects["CLIstd", "IRR"] < 1, 
           "Protective effect - healthier lifestyle reduces cancer risk",
           "Risk factor - unhealthy lifestyle increases cancer risk"), "\n")

cat("SmokingPrevalence:", 
    ifelse(lifestyle_effects["SmokingPrevalence", "IRR"] > 1, 
           "Risk factor - smoking increases cancer risk",
           "No significant risk or protective effect"), "\n")

cat("BMImedian:", 
    ifelse(lifestyle_effects["BMImedian", "IRR"] > 1, 
           "Risk factor - higher BMI increases cancer risk",
           "No significant risk or protective effect"), "\n")


#________________MODEL GOODNESS-OF-FIT (VALIDATION)________________________________
###################################################################################

# Residual analysis
par(mfrow = c(2, 2))  # Create a 2x2 grid for plots
plot(poisson_model)      # Generate 4 diagnostic plots
par(mfrow = c(1, 1))   # Reset to single plot layout

# Check for influential points
# We use Cook's distance for measuring observations that have unusual  
# influence of the model
influence_plot <- plot(poisson_model, which = 4)  


# Cook's distance > 1 is generally concerning
# Cook's distance > 4/(n-p) is often used as a threshold
n <- nrow(data)
p <- length(coef(poisson_model))
cooks_threshold <- 4/(n - p)

cat("Cook's distance threshold:", round(cooks_threshold, 3), "\n")
cat("Any points above this may be overly influential\n")

# Likelihood ratio test for model significance
# We compare full model (with all predictors) with 
# null model (intercept + offset)
null_model <- update(poisson_model, . ~ 1 + offset(logPopulation))
lr_test <- anova(null_model, poisson_model, test = "Chisq")
print("Likelihood Ratio Test:")
print(lr_test)


# Pseudo R-squared (just like R-squared in linear regression)
pseudo_r2 <- 1 - (poisson_model$deviance / poisson_model$null.deviance)
cat("Pseudo R-squared:", round(pseudo_r2, 3), "\n")


#___________________VISUALISATION OF MODEL RESULTS__________________________________
####################################################################################

# Plot predicted VS observed
predicted <- predict(poisson_model, type = "response")
observed <- data$NewCases

plot(observed, predicted, 
     xlab = "Observed Cases", ylab = "Predicted Cases",
     main = "Predicted vs Observed Cases",
     pch = 16, col = "blue")
abline(0, 1, col = "red", lwd = 2)

# Plot effects
model_effects <- allEffects(poisson_model)
plot(model_effects)

# Create a summary plot of key risk factors
key_factors <- lifestyle_effects
key_factors$Factor <- rownames(key_factors)

ggplot(key_factors, aes(x = Factor, y = IRR, ymin = Lower_CI, ymax = Upper_CI)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Risk Factors for Pancreatic Cancer",
       subtitle = "Incidence Rate Ratios with 95% Confidence Intervals",
       y = "Incidence Rate Ratio (IRR)",
       x = "Risk Factor") +
  theme_minimal() +
  coord_flip()






