
# Libraries
install.packages("RSiena")
install.packages("jtools")
install.packages("lmerTest")
library(effects)
library(ggplot2)
library(lme4)
library(RSiena)
library(jtools)
library(lmerTest)
library(boot)
install.packages("dplyr")
library(dplyr)
# Read data
data <- read.csv("data/Data_T2.csv")
data

str(data)

# ---- Explore the data visually. ----

# 1. Line plots
# Display the change in tumor volume over time for each category of mice
# and the treatment they have

xeno_group <- data[data$Model == 1,]
tumor_group <- data[data$Model ==0,]
control_group <- data[data$Treatment == 0,]
treated_group <- data[data$Treatment == 1,]

control_xeno <- control_group[control_group$Model==1,]
treated_xeno <- treated_group[treated_group$Model==1,]
control_tumor <- control_group[control_group$Model==0,]
treated_tumor <- treated_group[treated_group$Model==0,]

# Get averages for the plots

avg_ctr_xeno <- control_xeno %>%
  group_by(Time) %>%
  summarise(mean_DV = mean(DV, na.rm = TRUE))
avg_tr_xeno <- treated_xeno %>%
  group_by(Time) %>%
  summarise(mean_DV = mean(DV, na.rm = TRUE))
avg_ctr_tumor <- control_tumor %>%
  group_by(Time) %>%
  summarise(mean_DV = mean(DV, na.rm = TRUE))
avg_tr_tumor <- treated_tumor %>%
  group_by(Time) %>%
  summarise(mean_DV = mean(DV, na.rm = TRUE))

# ---- Plot all line plots ----

ggplot(data=control_xeno, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(size = 0.7) +  
  geom_line(data = avg_ctr_xeno, aes(x = Time, y = mean_DV), inherit.aes = FALSE, colour = "black", size = 1) +
  scale_color_discrete(name = "ID") +
  ggtitle("Control Xenograft Mice") + ylim(300,1400) +
  theme(plot.title = element_text(hjust = 0.5)) 

ggplot(data=treated_xeno, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(size = 0.7) + 
  geom_line(data = avg_tr_xeno, aes(x = Time, y = mean_DV), inherit.aes = FALSE, colour = "black", size = 1) +
  scale_color_discrete(name = "ID") + 
  ggtitle("Treated Xenograft Mice") + ylim(300,1400) +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(data=control_tumor, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(size = 0.7)+
  geom_line(data = avg_ctr_tumor, aes(x = Time, y = mean_DV), inherit.aes = FALSE, colour = "black", size = 1) +
  scale_color_discrete(name = "ID") + 
  ggtitle("Control Tumor Mice") + ylim(300,1400) +
  theme(plot.title = element_text(hjust = 0.5))
  
  
ggplot(data=treated_tumor, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(size = 0.7)+
  geom_line(data = avg_tr_tumor, aes(x = Time, y = mean_DV), inherit.aes = FALSE, colour = "black", size = 1) +
  scale_color_discrete(name = "ID") + 
  ggtitle("Treated Tumor Mice")  + ylim(300,1400) +
  theme(plot.title = element_text(hjust = 0.5))


# ---- Boxplots of each individual from each category ----

model.labs <- c("Tumor Mice", "Xenograft Mice")
names(model.labs) <- c("0", "1")


# Plot of boxplots of different Model's DV for each Treatment
# 0 = Tumor Group, 1 = Xenograft Group
p1 <- ggplot(treated_group, aes(x=ID, y=DV, group=ID,  color = as.factor(Model))) + 
  scale_color_discrete(name = "Model", labels=c("Tumor Mice","Xenograft Mice")) + 
  geom_boxplot() +
  facet_wrap(~Model, labeller = labeller(Model = model.labs)) + 
  ggtitle("Boxplots of Treated Group") + 
  theme(plot.title = element_text(hjust = 0.5))

p1

p2 <- ggplot(control_group, aes(x=ID, y=DV, group=ID,  color = as.factor(Model))) + 
  scale_color_discrete(name = "Model", labels=c("Tumor Mice","Xenograft Mice")) + 
  geom_boxplot() +
  facet_wrap(~Model, , labeller = labeller(Model = model.labs)) + 
  ggtitle("Boxplots of Control Group") + 
  theme(plot.title = element_text(hjust = 0.5))

p2

# Plot of histograms of different Treatment's DV for each Model
treat.labs <- c("Control", "Treated")
names(treat.labs) <- c("0", "1")

p3 <- ggplot(xeno_group, aes(x=ID, y=DV, group=ID,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment", labels=c("Control","Treated")) + 
  geom_boxplot() +
  facet_wrap(~Treatment, labeller = labeller(Treatment = treat.labs)) + 
  ggtitle("Boxplots of Xenograft Group") + 
  theme(plot.title = element_text(hjust = 0.5))

p3

p4 <- ggplot(tumor_group, aes(x=ID, y=DV, group=ID,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment", labels=c("Control","Treated")) + 
  geom_boxplot() +
  facet_wrap(~Treatment,labeller = labeller(Treatment = treat.labs)) + 
  ggtitle("Boxplots of Tumor Group") + 
  theme(plot.title = element_text(hjust = 0.5))

p4

# Boxplot of avg
p5 <- ggplot(tumor_group, aes(x=Treatment, y=DV,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment", labels=c("Control","Treated")) + 
  geom_boxplot() + ylim(300,1250) +
  ggtitle("Averaged Boxplots of Tumor Group") + 
  theme(plot.title = element_text(hjust = 0.5))
p5

# Boxplot of avg
p6 <- ggplot(xeno_group, aes(x=Treatment, y=DV,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment", labels=c("Control","Treated")) + 
  geom_boxplot() +ylim(300,1250) +
  ggtitle("Averaged Boxplots of Xenograft Group") + 
  theme(plot.title = element_text(hjust = 0.5))
p6

# Boxplot over time for each category
p7 <- ggplot(xeno_group, aes(x=Time, y=DV, group = Time,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment", labels=c("Control","Treated")) + 
  geom_boxplot() + facet_wrap(~Treatment,labeller = labeller(Treatment = treat.labs)) +
  ggtitle("Boxplots of Xenograft Group over Time") + 
  theme(plot.title = element_text(hjust = 0.5))
p7

# Boxplot over time for each category
p8 <- ggplot(tumor_group, aes(x=Time, y=DV, group = Time,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment", labels=c("Control","Treated")) + 
  geom_boxplot() + facet_wrap(~Treatment, labeller = labeller(Treatment = treat.labs)) +
  ggtitle("Boxplots of Tumor Group over Time") + 
  theme(plot.title = element_text(hjust = 0.5))
p8

# ----- Modelling -----

# Null model to see baseline
m0_lme <- lmer(DV~1 + (1 | ID), tumor_group)
plot(m0_lme)
summary(m0_lme)

# Mixed effects model treatment w/ partial pooling (tumor group)
tumor_model <- lmer(DV~1 + Treatment + (1|ID), tumor_group)
summary(tumor_model)
plot(tumor_model)
summ(tumor_model)
# MODEL FIT:
# AIC = 26636.98, BIC = 26659.73
#Pseudo-R² (fixed effects) = 0.16
# Pseudo-R² (total) = 0.50
# Fixed effects is 0.16, fixed + inter individual is 0.5
# ICC 0.4, varianxce by individual (not so much)

ranova(tumor_model)

# Mixed effects model treatment w/ pooling (xenograft group)
xeno_model <- lmer(DV~1 + Treatment + (1|ID), xeno_group)
summary(xeno_model)
plot(xeno_model)
summ(xeno_model)
#AIC = 41177.46, BIC = 41201.84
#Pseudo-R² (fixed effects) = 0.10
#Pseudo-R² (total) = 0.44

# Likelihood ratio test tumor group
m_null <- lmer(DV~1 + (1|ID), tumor_group, REML = FALSE)
m_treat <- lmer(DV~1 + Treatment + (1|ID), tumor_group, REML = FALSE)
anova(m_null, m_treat)
# pval = 0.001509, treatment performs better

m_null <- lmer(DV~1 + (1|ID), xeno_group, REML = FALSE)
m_treat <- lmer(DV~1 + Treatment + (1|ID), xeno_group, REML = FALSE)
anova(m_null, m_treat)
# pval = 0.002323 , treatment performs better than null

# Effect Plot
plot(allEffects(xeno_model))

# Try predictions

tumor_group$fit <- predict(m2_lme)   #Add model fits to dataframe

p <- ggplot(tumor_group, aes(x = Time, y = DV, colour = ID)) +
  geom_point(size=3) +
  geom_line(aes(y = predict(m2_lme)),size=1) 
p


# Bootstrap
confint(tumor_model, parm = c(3,4), method ="boot", nsim = 1000, boot.type = "perc")
#                2.5 %    97.5 %
#(Intercept)  636.3281 729.48496
#Treatment   -190.9064 -48.40985


# ----- Conditional growth model ----
cond_model_tumor <- lmer(DV ~ 1 + Treatment * Time + (1 + Time | ID), tumor_group)
summary(cond_model_tumor)
summ(cond_model_tumor)
# AIC = 21324.38, BIC = 21369.89
# Pseudo-R² (fixed effects) = 0.59
# Pseudo-R² (total) = 0.96 
# ICC = 0.79

# Visualize model 
# Create predicted values
tumor_group$predicted_cond <- predict(cond_model_tumor)

# Group-level (averaged) predictions
avg_preds <- tumor_group %>%
  group_by(Treatment, Time) %>%
  summarise(mean_pred = mean(predicted_cond), .groups = "drop")

# Plot

# Plot actual DV and predicted values over time by ID, faceted by Treatment
ggplot(tumor_group, aes(x = Time, y = DV, color = as.factor(ID))) +
  geom_point(alpha = 0.6) +  # actual data points
  geom_line(aes(y = predicted_cond), size = 1) +  # prediction lines
  facet_wrap(~ Treatment, labeller = labeller(Treatment = treat.labs)) +  # separate panels by treatment
  labs(
    title = "Tumor Group - DV over Time by ID, Faceted by Treatment",
    x = "Time (days)",
    y = "DV",
    color = "ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


# Plot only treatment-level prediction lines
ggplot(avg_preds, aes(x = Time, y = mean_pred, color = as.factor(Treatment))) +
  geom_line(size = 1.5) +
  labs(
    title = "Average Predicted DV over Time by Treatment",
    x = "Time (days)",
    y = "Predicted DV",
    color = "Treatment"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

fit2.m <- predict(cond_model_tumor, re.form = NA)


# Plot with marginal fit
# Plot actual DV and predicted values over time by ID, faceted by Treatment
ggplot(tumor_group, aes(x = Time, y = DV, color = as.factor(Treatment))) +
  geom_point(alpha = 0.6) +  # actual data points
  geom_line(aes(y = fit2.m), size = 1) +  # prediction lines
  facet_wrap(~ Treatment, labeller = labeller(Treatment = treat.labs)) +  # separate panels by treatment
  labs(
    title = "Xenograft Group - DV over Time by ID, Marginal Fit, Faceted by Treatment",
    x = "Time (days)",
    y = "DV",
    color = "Treatment"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


plot(allEffects(cond_model_tumor))

cond_model_xeno <- lmer(DV ~ 1 + Treatment * Time + (1 + Time | ID), xeno_group)
summary(cond_model_xeno)
summ(cond_model_xeno)
# AIC = 21324.38, BIC = 21369.89
# Pseudo-R² (fixed effects) = 0.59
# Pseudo-R² (total) = 0.96 
# ICC = 0.79

plot(allEffects(cond_model_xeno))


# Visualize model 
# Create predicted values
xeno_group$predicted_cond <- predict(cond_model_xeno)
fit2.m <- predict(cond_model_xeno, re.form = NA)
fit2.c <- predict(cond_model_xeno, re.form = NULL)
# Group-level (averaged) predictions
avg_preds <- xeno_group %>%
  group_by(Treatment, Time) %>%
  summarise(mean_pred = mean(predicted_cond), .groups = "drop")

# Plot

# Plot actual DV and predicted values over time by ID, faceted by Treatment
ggplot(xeno_group, aes(x = Time, y = DV, color = as.factor(ID))) +
  geom_point(alpha = 0.6) +  # actual data points
  geom_line(aes(y = predicted_cond), size = 1) +  # prediction lines
  facet_wrap(~ Treatment, labeller = labeller(Treatment = treat.labs)) +  # separate panels by treatment
  labs(
    title = "Xenograft Group - DV over Time by ID, Faceted by Treatment",
    x = "Time (days)",
    y = "DV",
    color = "ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Plot with marginal fit
# Plot actual DV and predicted values over time by ID, faceted by Treatment
ggplot(xeno_group, aes(x = Time, y = DV, color = as.factor(Treatment))) +
  geom_point(alpha = 0.6) +  # actual data points
  geom_line(aes(y = fit2.m), size = 1) +  # prediction lines
  facet_wrap(~ Treatment, labeller = labeller(Treatment = treat.labs)) +  # separate panels by treatment
  labs(
    title = "Xenograft Group - DV over Time by ID, Marginal Fit, Faceted by Treatment",
    x = "Time (days)",
    y = "DV",
    color = "Treatment"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Plot only treatment-level prediction lines
ggplot(avg_preds, aes(x = Time, y = fit2.m, color = as.factor(Treatment))) +
  geom_line(size = 1.5) +
  labs(
    title = "Average Predicted DV over Time by Treatment",
    x = "Time (days)",
    y = "Predicted DV",
    color = "Treatment"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
