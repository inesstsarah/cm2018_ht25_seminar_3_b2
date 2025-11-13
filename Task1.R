library(visdat)
library(ggplot2)
library(GGally)
library(corrplot)
library(MASS)
library(lme4)
library(DHARMa)
library(performance)
library(effects)
library(car)
library(dplyr)

# Load data 
data <- read.csv("data/Data_T1.csv")

# Creates a variable, incidence rate, since all subgroups (region-sex-age) don't have identical Npopulations (although they are similar)
data$IncidenceRate <- (data$NewCases / data$Npopulation) * 10000

# Overview of data
str(data)
summary(data)

# Explore missingness
vis_miss(data) # No missing data


# ----------- General visualisations -------------------------------------------
# Incidence rate by region
ggplot(data, aes(x = Region, y = IncidenceRate, fill = Region)) +
  geom_boxplot() +
  theme_minimal()

# Incidence rate by age group and sex, per region
ggplot(data, aes(x = AgeGroup, y = IncidenceRate, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ Region) + 
  theme_minimal()


# ------------ Lifestyle factors & cancer ----------------------------------------------
# Relationship between standardised CLI and cancer incidence for each region
ggplot(data, aes(x = CLIstd, y = IncidenceRate, color = Region)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()
# Healthier lifestyle = lower incidence rate for all regions

# Explore whether CLI affects cancer risk differently by age or sex
ggplot(data, aes(x = CLIstd, y = IncidenceRate, color = AgeGroup)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Sex) +
  theme_minimal()
# Healthier lifestyle = lower incidence rate for both sexes for all age groups except females 60-79

# Relationship between smoking prevalence and cancer incidence, higher values mean healthier lifestyle
ggplot(data, aes(x = SmokingPrevalence, y = IncidenceRate, color = Region)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
  theme_minimal()
# Shows higher incidence rate for smokers in region 2, 3 & 5, lower incidence rate for smokers in region 1 & 4
  
# Explore whether smoking affects cancer risk differently by age or sex
ggplot(data, aes(x = SmokingPrevalence, y = IncidenceRate, color = AgeGroup)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Sex) +
  theme_minimal()
# Smoking increases incidence rates for both sexes for all age groups except females 60-79

# Relationship between median BMI and cancer incidence
ggplot(data, aes(x = BMImedian, y = IncidenceRate, color = Region)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
theme_minimal()
# Higher median BMI drastically increases incidence rate in region 2 & 4, slightly increases rate in region 5 
# and decreases rate in region 1 & 3

# Explore whether BMI affects cancer risk differently by age or sex
ggplot(data, aes(x = BMImedian, y = IncidenceRate, color = AgeGroup)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Sex) +
  theme_minimal()
# Higher median BMI drastically increases incidence rate for both sexes 40-59, females 60-79 and males 20-39,
# slightly decreases rate for females 20-39 and drastically decreases rate for males 60-79


# ----------- Correlation ------------------------------------------------------
# Correlations and scatterplots between numeric variables
numeric_data <- data[, c("CLIstd", "SmokingPrevalence", "BMImedian", "IncidenceRate")]
ggpairs(numeric_data)

# Correlation plot of numeric variables
cor_data <- cor(numeric_data)
corrplot(cor_data, addCoef.col = "black")
# CLIstd is moderately negatively correlated to incidence rate, meaning healthier lifestyle = lower risk for cancer
# Smoking and BMI is slightly positively correlated to incidence rate, meaning smoking and higher BMI = higher risk for cancer





# __________________ MAKING THE MODEL _______________________________________________
#####################################################################################

# Since we're modeling counts with different population sizes at risk,
# we should use an offset for population size to handle differences 
data$logPopulation <- log(data$Npopulation)


#__________________STRUCTURE OF MODEL MAKING__________________________________________

# 1. Count data makes Poisson model suitable. We use glm and glmer for creating models
#    since data is not normal nor continuous (Generalised Linear Model)
#    We explore what Poisson model fits the best:

          #  Poisson                            (initial model with independent effects)   <-- (Fixed effects)
          #  Poisson with interactions          (combining for ex. CLIstd and sex)         <-- (Fixed effects)
          #  Poisson with random effects        (good with clustered structured data)      <-- (Mixed effects)

# 2. We we check for multicollinearity through VIFs (Variance Inflation Factors) after
#    the first initial model is made to see if model needs adjustment 

#    If : VIF < 3         No problem    
#         VIF = (3-5)     Moderate correlation
#         VIF > (5-10)    Problematic (consider dropping variable)

# 3. Each model performance will be evaluated through:

          # Dispersion 
          # AIC
          # Likelihood ratio 
          # Residual analysis

# 4. Best fitted model will be chosen for final conclusions and visualizations
#________________________________________________________________________________________



# INITIAL POISSON MODEL 
#________________________________________________________________________________________
poisson_model <- glm(NewCases ~ Region + AgeGroup + Sex + CLIstd + 
                       SmokingPrevalence + BMImedian + offset(logPopulation),
                     family = poisson(link = "log"),
                     data = data)

cat("\n MODEL 1: INITIAL POISSON MODEL  \n")
summary(poisson_model)


# Checking dispersion
disp1   <- sum(residuals(poisson_model, type = "pearson")^2)/df.residual(poisson_model)
cat("Dispersion: ", round(disp1,3), "\n")

# CHECKING FOR MULTICOLLINEARITY
#--------------------------------
cat("\n CHECKING FOR MULTICOLLINEARITY IN INITIAL MODEL \n")
vif_initialModel <- vif(poisson_model)
print(vif_initialModel)
#_________________________________________________________________________________________



# POISSON MODEL WITH INTERACTIONS (CLIstd, SEX, AND AGE)
#_________________________________________________________________________________________
poisson_int <- glm(NewCases ~ Region + AgeGroup*CLIstd + Sex*CLIstd + 
                     SmokingPrevalence + BMImedian + offset(logPopulation),
                   family = poisson(link = "log"), data = data)

cat("\n MODEL 2: POISSON MODEL WITH INTERACTIONS  \n")
summary(poisson_int)


# Checking dispersion
disp2   <- sum(residuals(poisson_int, type = "pearson")^2)/df.residual(poisson_int)
cat("Dispersion: ", round(disp2,3), "\n")
#________________________________________________________________________________________




# POISSON MIXED EFFECT MODEL (RANDOM INTERCEPT BY REGION)
#________________________________________________________________________________________
poisson_glmm <- glmer(NewCases ~ AgeGroup + Sex + CLIstd + SmokingPrevalence +
                        BMImedian + offset(logPopulation) + (1 | Region), 
                      family = poisson(link = "log"), data = data)

cat("\n MODEL 3: POISSON MIXED EFFECTS MODEL  \n")
summary(poisson_glmm)

# Dispersion cant be looked at in mixed effect models since variance is composed of 
# two parts. Second variance comes from random effect of Region which adds variability 
# across clusters.
#_________________________________________________________________________________________


#________________MODEL PERFORMANCE COMPARISON__________________________________________
#######################################################################################
cat("\n  MODEL COMPARISON \n")

# Compare AICs
model_compare <- data.frame(             # Using data.frame to get table
  Model = c(" Initial Poisson", "Poisson interactions", "Poisson GLMM"),
  AIC = c(AIC(poisson_model), AIC(poisson_int), AIC(poisson_glmm))
)
print(model_compare)


# Compare Likelihood ratio tests (OBS! Only for fixed effect glms)
cat("\n  Likelihood Ratio Tests \n")
cat("\n  OBS! Only on fixed effect glm´s (poisson_model and poisson_int")
anova(poisson_model, poisson_int, test = "Chisq")


# (R^2 value comparison is skipped due to glms models having pure R^2,
# while glmm has the pseudo R^2.)


# Compare DHARMa residuals
res1 <- simulateResiduals(poisson_model)
res2 <- simulateResiduals(poisson_int)
res3 <- simulateResiduals(poisson_glmm)

# Plotting DHARMa residuals
par(mfrow = c(3, 2))
plot(res1, main = "DHARMa - Initial Poisson")
plot(res2, main = "DHARMa - Poisson with Interaction")
plot(res3, main = "DHARMa - Poisson GLMM")
par(mfrow = c(1, 1))

# Test zero inflation ("when model says there is too many zeros")

#               If p<0.05  --> too many zeros --> zero inflated
#               If p>0.05  --> no zero inflation

testZeroInflation(res1)
testZeroInflation(res2)
testZeroInflation(res3)


# _______________________CHOSEN MODEL_______________________________________________
######################################################################################
# Based on lowest AIC and best residual diagnostics, Poisson GLMM was chosen as final model
finalModel <- poisson_glmm
######################################################################################




#_________________ INTERPRETATION OF CHOSEN MODEL _____________________________
#####################################################################################

# IRR (Incident )
# We do this since Poission regression model has the coeff. in log-scale 

#        If IRR>1    <-- indicates increased risk
#        If IRR<1    <-- indicates protective effect 


# Extract fixed effects only    <-- (Since glmm has both fixed and random effects)
irr <- exp(fixef(poisson_glmm))

# Confidence intervals for fixed effects only
# beta and Wald helps only using the fixed effects
irr_ci <- exp(confint(poisson_glmm, parm = "beta_", method = "Wald"))

results_table <- data.frame(
  Variable = names(irr),
  IRR = round(irr, 3),
  Lower_CI = round(irr_ci[,1], 3),
  Upper_CI = round(irr_ci[,2], 3)
)

print(results_table)





#___________________VISUALISATION OF MODEL RESULTS__________________________________
####################################################################################


# Effect plot 
plot(allEffects(finalModel))


# Predicted VS Observed cases
pred <- predict(finalModel, type = "response")
obs  <- data$NewCases

ggplot(data, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", lwd = 1.2) +
  labs(title = "Observed vs Predicted Counts",
       x = "Observed cases",
       y = "Predicted cases") +
  theme_minimal()

# DHARMa residuals plot again for final model
res_final <- simulateResiduals(finalModel)
plot(res_final)

# Predicted incidence by age group
# Change Sex and Region for all plots OBS!
newdata_age <- data %>% 
  group_by(AgeGroup) %>% 
  summarize(
    CLIstd = mean(CLIstd),
    SmokingPrevalence = mean(SmokingPrevalence),
    BMImedian = mean(BMImedian),
    Sex = "Male",
    Region = "Region1",
    logPopulation = mean(logPopulation)
  )

newdata_age$pred <- predict(finalModel, newdata = newdata_age, type = "response")

ggplot(newdata_age, aes(x = AgeGroup, y = pred)) +
  geom_col(fill = "steelblue") +
  labs(title = "Predicted Pancreatic Cancer Cases by Age Group",
       y = "Predicted Cases (per stratum)") +
  theme_minimal()


# Predicted incidence VS CLIstd
new_cli <- data.frame(
  CLIstd = seq(min(data$CLIstd), max(data$CLIstd), length = 100),
  SmokingPrevalence = mean(data$SmokingPrevalence),
  BMImedian = mean(data$BMImedian),
  Sex = "Male",
  AgeGroup = "40-59",
  Region = "Region1",
  logPopulation = mean(data$logPopulation)
)

new_cli$pred <- predict(finalModel, newdata = new_cli, type="response")

ggplot(new_cli, aes(x = CLIstd, y = pred)) +
  geom_line(size = 1.2, color = "darkgreen") +
  labs(title = "Predicted Cancer Cases vs. Lifestyle Index (CLIstd)",
       y = "Predicted Cases") +
  theme_minimal()

# Forest plot of final model IRR
irr <- exp(fixef(finalModel))
irr_ci <- exp(confint(finalModel, parm = "beta_", method = "Wald"))

irr_df <- data.frame(
  Variable = names(irr),
  IRR = irr,
  Lower = irr_ci[,1],
  Upper = irr_ci[,2]
)

ggplot(irr_df, aes(x = reorder(Variable, IRR), y = IRR)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Incidence Rate Ratios (IRR) of Final Model",
       y = "IRR (log scale)") +
  scale_y_log10() +
  theme_minimal()

