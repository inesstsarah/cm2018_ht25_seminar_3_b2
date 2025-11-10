
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

ggplot(data=control_xeno, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(linetype = "dashed")+ 
  geom_point() +  scale_color_discrete() +
  ggtitle("Control Xenograph Mice") + 
  theme(plot.title = element_text(hjust = 0.5)) 




ggplot(data=treated_xeno, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(linetype = "dashed")+
  geom_point() + scale_color_discrete(name = "ID") + 
  ggtitle("Treated Xenograph Mice") + 
  theme(plot.title = element_text(hjust = 0.5))



ggplot(data=control_tumor, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(linetype = "dashed")+
  geom_point() + scale_color_discrete(name = "ID") + 
  ggtitle("Control Tumor Mice") + 
  theme(plot.title = element_text(hjust = 0.5))
  
  

ggplot(data=treated_tumor, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(linetype = "dashed")+
  geom_point() + scale_color_discrete(name = "ID") + 
  ggtitle("Treated Tumor Mice") + 
  theme(plot.title = element_text(hjust = 0.5))


# 2. Boxplots of each individual from each category

# Plot of histograms of different Model's DV for each Treatment
# 0 = Tumor Group, 1 = Xenograft Group
p1 <- ggplot(treated_group, aes(x=ID, y=DV, group=ID,  color = as.factor(Model))) + 
  scale_color_discrete(name = "Model") + 
  geom_boxplot() +
  facet_wrap(~Model) + ggtitle("Histograms of Treated Group") + 
  theme(plot.title = element_text(hjust = 0.5))

p1

p2 <- ggplot(control_group, aes(x=ID, y=DV, group=ID,  color = as.factor(Model))) + 
  scale_color_discrete(name = "Model") + 
  geom_boxplot() +
  facet_wrap(~Model) + ggtitle("Histograms of Control Group") + 
  theme(plot.title = element_text(hjust = 0.5))

p2

# Plot of histograms of different Treatment's DV for each Model
p3 <- ggplot(xeno_group, aes(x=ID, y=DV, group=ID,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment") + 
  geom_boxplot() +
  facet_wrap(~Treatment) + ggtitle("Histograms of Xenograft Group") + 
  theme(plot.title = element_text(hjust = 0.5))

p3

p4 <- ggplot(tumor_group, aes(x=ID, y=DV, group=ID,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment") + 
  geom_boxplot() +
  facet_wrap(~Treatment) + ggtitle("Histograms of Tumor Group") + 
  theme(plot.title = element_text(hjust = 0.5))

p4

# Boxplot of avg
p5 <- ggplot(tumor_group, aes(x=ID, y=DV,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment") + 
  geom_boxplot() +
  ggtitle("Histograms of Tumor Group") + 
  theme(plot.title = element_text(hjust = 0.5))
p5

# Boxplot of avg
p6 <- ggplot(xeno_group, aes(x=ID, y=DV,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment") + 
  geom_boxplot() +
  ggtitle("Histograms of Xenograft Group") + 
  theme(plot.title = element_text(hjust = 0.5))
p6

# Boxplot over time for each category
p6 <- ggplot(xeno_group, aes(x=Time, y=DV, group = Time,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment") + 
  geom_boxplot() + facet_wrap(~Treatment) +
  ggtitle("Histograms of Xenograft Group over Time") + 
  theme(plot.title = element_text(hjust = 0.5))
p6

# Boxplot over time for each category
p7 <- ggplot(tumor_group, aes(x=Time, y=DV, group = Time,  color = as.factor(Treatment))) + 
  scale_color_discrete(name = "Treatment") + 
  geom_boxplot() + facet_wrap(~Treatment) +
  ggtitle("Histograms of Tumor Group over Time") + 
  theme(plot.title = element_text(hjust = 0.5))
p7


# Create mixed-effect model 
# Try with linear mixed-effect model
# Fixed effects fit: Complete Pooling

m1_fe1 <- lm(formula = DV ~ Treatment + factor(ID)-1, data = tumor_group)
plot(m1_fe1)
summary(m1_fe1)

# Null model 
m0_lme <- lmer(DV~1 + (1 | ID), tumor_group)
plot(m0_lme)
summary(m0_lme)

# Mixed effects model treatment w/ pooling (tumor group)
m2_lme <- lmer(DV~1 + Treatment + (1|ID), tumor_group)
summary(m2_lme)
plot(m2_lme)
summ(m2_lme)
# MODEL FIT:
# AIC = 26636.98, BIC = 26659.73
#Pseudo-R² (fixed effects) = 0.16
# Pseudo-R² (total) = 0.50
# Fixed effects is 0.16, fixed + inter individual is 0.5
# ICC 0.4, varianxce by individual (not so much)

ranova(m2_lme)

# Mixed effects model treatment w/ pooling (xenograft group)
m3_lme <- lmer(DV~1 + Treatment + (1|ID), xeno_group)
summary(m3_lme)
plot(m3_lme)
summ(m3_lme)
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
plot(allEffects(m2_lme))

# Bootstrap
confint(m2_lme, parm = c(3,4), method ="boot", nsim = 1000, boot.type = "perc")
#                2.5 %    97.5 %
#(Intercept)  636.3281 729.48496
#Treatment   -190.9064 -48.40985


# ----- Conditional growth model ----
model_tumor <- lmer(DV ~ 1 + Treatment * Time + (1 + Time | ID), tumor_group)
summary(model_tumor)
summ(model_tumor)
# AIC = 21324.38, BIC = 21369.89
# Pseudo-R² (fixed effects) = 0.59
# Pseudo-R² (total) = 0.96 
# ICC = 0.79
plot(allEffects(model_tumor))

model_xeno <- lmer(DV ~ 1 + Treatment * Time + (1 + Time | ID), xeno_group)
summ(model_xeno)
# AIC = 21324.38, BIC = 21369.89
# Pseudo-R² (fixed effects) = 0.59
# Pseudo-R² (total) = 0.96 
# ICC = 0.79
plot(allEffects(model_xeno))
summary(model_xeno)
summ(model_xeno)
