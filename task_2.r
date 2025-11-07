
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

# Display the change in tumor volume over time for each category of mice
# and the treatment they have
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
  theme(plot.title = element_text(hjust = 0.5)) +




ggplot(data=treated_xeno, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(linetype = "dashed")+
  geom_point() + scale_color_discrete(name = "ID") + 
  ggtitle("Treated Xenograph Mice") + 
  theme(plot.title = element_text(hjust = 0.5))



ggplot(data=control_tumor, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(linetype = "dashed")+
  geom_point() + scale_color_discrete(name = "ID") + 
  ggtitle("Control Xenograph Mice") + 
  theme(plot.title = element_text(hjust = 0.5))
  
  

ggplot(data=treated_tumor, aes(x=Time, y=DV, group=as.factor(ID), colour = as.factor(ID))) +
  geom_line(linetype = "dashed")+
  geom_point() + scale_color_discrete(name = "ID") + 
  ggtitle("Control Xenograph Mice") + 
  theme(plot.title = element_text(hjust = 0.5))

# Get average of each of the mice at each time point and plot them all into 
# one plot 


# Create mixed-effect model 

# Try with linear mixed-effect model
data_xeno <- data[data$Model==1,]

str(data_xeno)



