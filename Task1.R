library(visdat)
library(ggplot2)
library(GGally)
library(corrplot)

data <- read.csv("data/Data_T1.csv")

# Handle region and sex as factors instead of categorical
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
  geom_smooth(method = "lm", se = FALSE) +
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
  geom_smooth(method = "lm", se = FALSE) +
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

