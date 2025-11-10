
# Libraries
install.packages("RSiena")
library(ggplot2)
library(lme4)
library(RSiena)

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

#
# Get average of each of the mice at each time point and plot them all into 
# one plot 

lm(formula = DV ~ Time, data = treated_tumor)

# Create mixed-effect model 
# Try with linear mixed-effect model
# Fixed effects fit: Complete Pooling


