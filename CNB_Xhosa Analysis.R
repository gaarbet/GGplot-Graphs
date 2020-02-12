#CNB Xhosa
#meet with ruben to go over
#Get mixed models correctly identified and put onto graphs
#fix difficult to see graphs
#show N's, and significance as well, etc before each slide.


library(psych)
library(dplyr)
library(tidyverse)
library(tidyr)
library(Rmisc)
library(visreg)
library(ggplot2)
library(lme4)
library(lmerTest)
rm(list=ls())
df <- read.csv("xhosa_cnb_merged_clinical.csv")
df <- df[which(df$valid_code == "V" | df$valid_code == "VC" | df$valid_code == "V1"),]

for (i in 1:length(df$PCET_A_valid_code)){
  if(df[i,"PCET_A_valid_code"] != "V"){#V is the only valid code for this file
    df[i,"PCET_A_PCET_ACC2"] <- NA
  }
  if(df[i,"CPF_A_valid_code"] != "V"){
    df[i,"CPF_A_CPF_CR"] <- NA
  }
  if(df[i,"ER40_A_valid_code"] != "V"){
    df[i,"ER40_A_ER40_CR"] <- NA
  }
  if(df[i,"PMAT24_A_valid_code"] != "V"){
    df[i,"PMAT24_A_PMAT24_A_CR"] <- NA
  }
  if(df[i,"SPCPTNL_valid_code"] != "V"){
    df[i,"SPCPTNL.SCPT_TP"] <- NA
  }
  if(df[i,"VSPLOT24_valid_code"] != "V"){
    df[i,"VSPLOT24_VSPLOT24_CR"] <- NA
  }
  if(df[i,"SVOLT_A_valid_code"] != "V"){
    df[i,"SVOLT_A_SVOLT_CR"] <- NA
  }
  if(df[i,"SFNB2_valid_code"] != "V"){
    df[i,"SFNB2_SFNB_MCR"] <- NA
  }
}
df<- df[which(df$CASE.CONTROL == "case" | df$CASE.CONTROL =="control"),]

#add in tests not already converted to z scores
pc_sd <- sd(df$PCET_A_PCET_ACC2,na.rm = T)
pc_mean <- mean(df$PCET_A_PCET_ACC2,na.rm = T)
for (i in 1:length(df$PCET_A_PCET_ACC2)){
  df[i,322] <- (df[i,"PCET_A_PCET_ACC2"] - pc_mean)/pc_sd
}
colnames(df)[322] <- "PCET"

pc_sd <- sd(df$SPCPTNL_SCPT_TP,na.rm = T)
pc_mean <- mean(df$SPCPTNL_SCPT_TP,na.rm = T)
for (i in 1:length(df$SPCPTNL_SCPT_TP)){
  df[i,323] <- (df[i,"SPCPTNL_SCPT_TP"] - pc_mean)/pc_sd
}
colnames(df)[323] <- "CPT"
colnames(df)[45] <- "CPF"
colnames(df)[97] <- "ER40"
colnames(df)[174] <- "SFNB2"
colnames(df)[216] <- "PMAT24"
colnames(df)[223] <- "SVOLT"
colnames(df)[251] <- "VSPLOT24"


#Melt/Reformat file for plotting
x <- df %>% 
  pivot_longer(cols  = c(45,97,174,216,223,251,322,323), names_to = c('test')) 
#For only Case/Control
temp <- x[(x$CASE.CONTROL == "case" | x$CASE.CONTROL =="control"),]

tgc <- summarySE(temp, measurevar = "value", groupvars = c("test","CASE.CONTROL","gender"), na.rm = T)
pd <- position_dodge(0.30)
ggplot(tgc, aes(x=test, y=value,col = gender, shape = CASE.CONTROL,group = interaction(gender,CASE.CONTROL))) + 
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.4, position=pd) +
  geom_line(position=pd, size = .23) +
  geom_point(position=pd, size = 6) +scale_shape_manual(name = "Case/Control",values = c(8,5))+
  scale_x_discrete(labels=c("CPF", "CPT", "ER40","PCET", "PMAT24","SFNB2","SVOLT","VSPLOT"))+
  scale_colour_manual(name = "Gender", values = c("red","blue"))+ggtitle("CNB Xhosa - Accuracy by Case Type/Gender")+
  xlab("Test Name") + ylab("Z score (Accuracy)")+theme(legend.justification = "top")+
  geom_hline(yintercept=0.00, size = .25)+ theme(axis.title.x = element_text(size = 14))+ theme(axis.title.y = element_text(size = 14))

#marginal values
#mixed model
#Dont run before age brackets
temp$test <- as.factor(temp$test)
temp$CASE.CONTROL <- droplevels(temp$CASE.CONTROL)
temp$CASE.NO <- as.factor(temp$CASE.NO)


#me2 <- lmer(value ~ CASE.CONTROL + age + gender + test + (1|CASE.NO), data = temp)
me <- lmer(value ~ (CASE.CONTROL+gender+age+test)^2 + (1 | CASE.NO), data = temp)
summary(me)
#summary(lmer(value ~ CASE.CONTROL+gender+age+test + test:CASE.CONTROL+test:gender+test:age +(1 | CASE.NO), data = temp))
me2 <- lmer(value ~ CASE.CONTROL*test + age + gender + (1|CASE.NO), data = temp)#test/case
summary(me2)
visreg(me2,xvar="test", by = "CASE.CONTROL",points=list(cex=.2, pch=5, col = c("green","orange")),line=list(col=c("green","orange")), breaks = 3,jitter = T, overlay = T, xlab = "Test",ylab = "Zscore (accuracy)")#test/case/zscore





me2 <- lmer(value ~ CASE.CONTROL*gender + age + test + (1|CASE.NO), data = temp)#Case/gender interaction
visreg(me2,xvar="CASE.CONTROL", by = "gender",points=list(cex=.2, pch=5),jitter = T, overlay = T, xlab = "Case/Control",ylab = "Zscore (accuracy)") #case/gender

me2 <- lmer(value ~ CASE.CONTROL*test + age + gender + (1|CASE.NO), data = temp)#test/case
visreg(me2,xvar="test", by = "CASE.CONTROL",points=list(cex=.2, pch=5, col = c("green","orange")),line=list(col=c("green","orange")), breaks = 3,jitter = T, overlay = T, xlab = "Test",ylab = "Zscore (accuracy)")#test/case/zscore

me2 <- lmer(value ~ CASE.CONTROL*age+test + gender + (1|CASE.NO), data = temp)
visreg(me2,xvar="CASE.CONTROL", by = "age",points=list(cex=.2, pch=5), breaks = 3,jitter = T, overlay = T, xlab = "Case/Control",ylab = "Zscore (accuracy)")#case/age

me2 <- lmer(value ~ CASE.CONTROL + test*age + gender + (1|CASE.NO), data = temp)
visreg(me2,xvar="test", by = "age", breaks = 3,points=list(cex=.2, pch=5),jitter = T, overlay = T, xlab = "Test",ylab = "Zscore (accuracy)")#test/age

for(i in 1:length(x$age)){
  if(x[i,"age"] %in% 9:25){
    x[i,"age"] <- "9-25"
  }
  if(x[i,"age"] %in% 26:33){
    x[i,"age"] <- "26-33"
  }
  if(x[i,"age"] %in% 34:41){
    x[i,"age"] <- "34-41"
  }
  if(x[i,"age"] %in% 42:49){
    x[i,"age"] <- "42-49"
  }
  if(x[i,"age"] %in% 50:99){
    x[i,"age"] <- "50+"
  }
  if(x[i,"age"] %in% 0:8){
    x[i,"age"] <- NA
  }
}

temp <- x[(x$CASE.CONTROL == "case" | x$CASE.CONTROL =="control"),]
tgc <- summarySE(temp, measurevar = "value", groupvars = c("age","CASE.CONTROL","gender"), na.rm = T)
#tgc <- tgc[1:24,]
pd <- position_dodge(0.30)
x$age <-as.factor(x$age)
tgc$age <- factor(tgc$age, levels = c("9-25","26-33","34-41","42-49","50+"))
ggplot(tgc, aes(x=age, y=value,col = gender, shape = CASE.CONTROL,group = interaction(gender,CASE.CONTROL))) + 
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.4, position=pd) +
  geom_line(position=pd, size = .23) +
  geom_point(position=pd, size = 6) +scale_shape_manual(name = "Case/Control",values = c(8,5))+
  scale_x_discrete(labels=c("9-25","26-33","34-41","42-49","50+"))+
  scale_colour_manual(name = "Gender", values = c("red","blue"))+ggtitle("CNB Xhosa - Accuracy by Case Type/Gender(No test)")+
  xlab("Age") + ylab("Z score (Accuracy)")+theme(legend.justification = "top")+
  geom_hline(yintercept=0.00, size = .25)+ theme(axis.title.x = element_text(size = 14))+ theme(axis.title.y = element_text(size = 14))
#ggplot no test, first bracket age




#Speed
rm(list=ls())
df <- read.csv("xhosa_cnb_merged_clinical.csv")
df <- df[which(df$valid_code == "V" | df$valid_code == "VC" | df$valid_code == "V1"),]

for (i in 1:length(df$PCET_A_valid_code)){
  if(df[i,"PCET_A_valid_code"] != "V"){
    df[i,"PCET_A_PCET_RTCR"] <- NA
  }
  if(df[i,"CPF_A_valid_code"] != "V"){
    df[i,"CPF_A_CPF_RTCR"] <- NA
  }
  if(df[i,"ER40_A_valid_code"] != "V"){
    df[i,"ER40_A_ER40_RTCR"] <- NA
  }
  if(df[i,"PMAT24_A_valid_code"] != "V"){
    df[i,"PMAT24_A_PMAT24_A_RTCR"] <- NA
  }
  if(df[i,"SPCPTNL_valid_code"] != "V"){
    df[i,"SPCPTNL_SCPT_TPRT"] <- NA
  }
  if(df[i,"VSPLOT24_valid_code"] != "V"){
    df[i,"VSPLOT24_VSPLOT24_RTCR"] <- NA
  }
  if(df[i,"SVOLT_A_valid_code"] != "V"){
    df[i,"SVOLT_A_SVOLT_RTCR"] <- NA
  }
  if(df[i,"SFNB2_valid_code"] != "V"){
    df[i,"SFNB2_SFNB_MRTC"] <- NA
  }
  if(df[i,"MPRACT_valid_code"] != "V"){
    df[i,"MPRACT_MP2RTCR"] <- NA
  }
  if(df[i,"SCTAP_valid_code"] != "V"){
    df[i,"SCTAP.SCTAP_TOT"] <- NA
  }
}
df <- df[which(df$CASE.CONTROL == "case" | df$CASE.CONTROL =="control"),]

pc_sd <- sd(df$PCET_A_PCET_RTCR,na.rm = T)
pc_mean <- mean(df$PCET_A_PCET_RTCR,na.rm = T)
for (i in 1:length(df$PCET_A_PCET_RTCR)){
  df[i,322] <- (df[i,"PCET_A_PCET_RTCR"] - pc_mean)/pc_sd
}
colnames(df)[322] <- "PCET"

pc_sd <- sd(df$SPCPTNL_SCPT_TPRT,na.rm = T)
pc_mean <- mean(df$SPCPTNL_SCPT_TPRT,na.rm = T)
for (i in 1:length(df$SPCPTNL_SCPT_TPRT)){
  df[i,323] <- (df[i,"SPCPTNL_SCPT_TPRT"] - pc_mean)/pc_sd  #
}
colnames(df)[323] <- "CPT"

pc_sd <- sd(df$MPRACT_MP2RTCR,na.rm = T)
pc_mean <- mean(df$MPRACT_MP2RTCR,na.rm = T)
for (i in 1:length(df$MPRACT_MP2RTCR)){
  df[i,324] <- (df[i,"MPRACT_MP2RTCR"] - pc_mean)/pc_sd  
}
colnames(df)[324] <- "MPRACT"

df$SCTAP_SCTAP_TOT <- df$SCTAP_SCTAP_TOT*-1

colnames(df)[46] <- "CPF"
colnames(df)[98] <- "ER40"
colnames(df)[175] <- "SFNB2"
colnames(df)[219] <- "PMAT24"
colnames(df)[224] <- "SVOLT"
colnames(df)[253] <- "VSPLOT24"
colnames(df)[198] <- "SCTAP"

x <- df %>% 
  pivot_longer(cols  = c(46,98,175,198,219,224,253,322,323,324), names_to = c('test')) 
x$value <- x$value *-1

x <- x[which(x$value <= 5),]
x <- x[which(x$value >= -5),]

#For only Case/Control
temp <- x[(x$CASE.CONTROL == "case" | x$CASE.CONTROL =="control"),]

tgc <- summarySE(temp, measurevar = "value", groupvars = c("test","CASE.CONTROL","gender"), na.rm = T)
pd <- position_dodge(0.35)
ggplot(tgc, aes(x=test, y=value,col = gender, shape = CASE.CONTROL,group = interaction(gender,CASE.CONTROL))) + 
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.4, position=pd) +
  geom_line(position=pd, size = .27) +
  geom_point(position=pd, size = 6) +scale_shape_manual(name = "Case/Control",values = c(8,5))+
  scale_x_discrete(labels=c("CPF", "CPT", "ER40","MPRACT","PCET","PMAT24","SCTAP","SFNB2","SVOLT","VSPLOT"))+
  scale_colour_manual(name = "Gender", values = c("red","blue"))+ggtitle("CNB Xhosa - Speed by Case Type/Gender")+
  xlab("Test Name") + ylab("Z score (Speed - Higher is faster)")+theme(legend.justification = "top")+
  geom_hline(yintercept=0.00, size = .25)#higher is faster


temp$test <- as.factor(temp$test)
temp$CASE.CONTROL <- droplevels(temp$CASE.CONTROL)
temp$CASE.NO <- as.factor(temp$CASE.NO)

me2 <- lmer(value ~ CASE.CONTROL + age + gender + test + (1|CASE.NO), data = temp)
me <- lmer(value ~ (CASE.CONTROL+gender+age+test)^2 + (1 | CASE.NO), data = temp)
summary(lmer(value ~ CASE.CONTROL+gender+age+test + test:CASE.CONTROL+test:gender+test:age +(1 | CASE.NO), data = temp))
summary(me)

me2 <- lmer(value ~ CASE.CONTROL*gender + age + test + (1|CASE.NO), data = temp)
visreg(me2,xvar="CASE.CONTROL", by = "gender",points=list(cex=.2, pch=5),jitter = T, overlay = T, xlab = "Case/Control",ylab = "Zscore (speed)")#case/gender/zscore
me2 <- lmer(value ~ CASE.CONTROL*test + age + gender + (1|CASE.NO), data = temp)
visreg(me2,xvar="test", by = "CASE.CONTROL",points=list(cex=.2, pch=5,col = c("green","orange")),line=list(col=c("green","orange")), breaks = 3,jitter = T, overlay = T, xlab = "Test",ylab = "Zscore (speed)")#test/case/zscore
me2 <- lmer(value ~ CASE.CONTROL*age+test + gender + (1|CASE.NO), data = temp)
visreg(me2,xvar="CASE.CONTROL", by = "age",points=list(cex=.2, pch=5), breaks = 3,jitter = T, overlay = T, xlab = "Case/Control",ylab = "Zscore (speed)")#case/age/zscore
me2 <- lmer(value ~ CASE.CONTROL + test*age + gender + (1|CASE.NO), data = temp)
visreg(me2,xvar="test", by = "age", breaks = 3,points=list(cex=.2, pch=5),jitter = T, overlay = T, xlab = "Test",ylab = "Zscore (speed)")#test/age/zscore


for(i in 1:length(x$age)){
  if(x[i,"age"] %in% 9:25){
    x[i,"age"] <- "9-25"
  }
  if(x[i,"age"] %in% 26:33){
    x[i,"age"] <- "26-33"
  }
  if(x[i,"age"] %in% 34:41){
    x[i,"age"] <- "34-41"
  }
  if(x[i,"age"] %in% 42:49){
    x[i,"age"] <- "42-49"
  }
  if(x[i,"age"] %in% 50:99){
    x[i,"age"] <- "50+"
  }
  if(x[i,"age"] %in% 0:8){
    x[i,"age"] <- NA
  }
}
temp <- x[(x$CASE.CONTROL == "case" | x$CASE.CONTROL =="control"),]
tgc <- summarySE(temp, measurevar = "value", groupvars = c("age","CASE.CONTROL","gender"), na.rm = T)
#tgc <- tgc[1:24,]
pd <- position_dodge(0.30)
x$age <-as.factor(x$age)
tgc$age <- factor(tgc$age, levels = c("9-25","26-33","34-41","42-49","50+"))
ggplot(tgc, aes(x=age, y=value,col = gender, shape = CASE.CONTROL,group = interaction(gender,CASE.CONTROL))) + 
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.4, position=pd) +
  geom_line(position=pd, size = .23) +
  geom_point(position=pd, size = 6) +scale_shape_manual(name = "Case/Control",values = c(8,5))+
  scale_x_discrete(labels=c("9-25","26-33","34-41","42-49","50+"))+
  scale_colour_manual(name = "Gender", values = c("red","blue"))+ggtitle("CNB Xhosa - Speed by Case Type/Gender(No test)")+
  xlab("Age") + ylab("Z score (speed)")+theme(legend.justification = "top")+
  geom_hline(yintercept=0.00, size = .25)+coord_cartesian(ylim=c(-2,1))+ scale_y_continuous(breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1))
#ggplot no test, first bracket age


######
######
#EFFICIENCY GRAPH
rm(list=ls())
df <- read.csv("xhosa_cnb_merged_clinical.csv")
df <- df[which(df$valid_code == "V" | df$valid_code == "VC" | df$valid_code == "V1"),]

for (i in 1:length(df$PCET_A_valid_code)){
  if(df[i,"PCET_A_valid_code"] != "V"){
    df[i,"PCET_A_PCET_ACC2"] <- NA
    df[i,"PCET_A_PCET_RTCR"] <- NA
  }
  if(df[i,"CPF_A_valid_code"] != "V"){
    df[i,"CPF_A_CPF_CR"] <- NA
    df[i,"CPF_A_CPF_RTCR"] <- NA
  }
  if(df[i,"ER40_A_valid_code"] != "V"){
    df[i,"ER40_A_ER40_CR"] <- NA
    df[i,"ER40_A_ER40_RTCR"] <- NA
  }
  if(df[i,"PMAT24_A_valid_code"] != "V"){
    df[i,"PMAT24_A_PMAT24_A_CR"] <- NA
    df[i,"PMAT24_A_PMAT24_A_RTCR"] <- NA
  }
  if(df[i,"SPCPTNL_valid_code"] != "V"){
    df[i,"SPCPTNL.SCPT_TP"] <- NA
    df[i,"SPCPTNL_SCPT_TPRT"] <- NA
  }
  if(df[i,"VSPLOT24_valid_code"] != "V"){
    df[i,"VSPLOT24_VSPLOT24_CR"] <- NA
    df[i,"VSPLOT24_VSPLOT24_RTCR"] <- NA
  }
  if(df[i,"SVOLT_A_valid_code"] != "V"){
    df[i,"SVOLT_A_SVOLT_CR"] <- NA
    df[i,"SVOLT_A_SVOLT_RTCR"] <- NA
  }
  if(df[i,"SFNB2_valid_code"] != "V"){
    df[i,"SFNB2_SFNB_MCR"] <- NA
    df[i,"SFNB2_SFNB_MRTC"] <- NA
  }
  
}
df <- df[which(df$CASE.CONTROL == "case" | df$CASE.CONTROL =="control"),]

pc_sd <- sd(df$PCET_A_PCET_RTCR,na.rm = T)
pc_mean <- mean(df$PCET_A_PCET_RTCR,na.rm = T)
for (i in 1:length(df$PCET_A_PCET_RTCR)){
  df[i,322] <- (df[i,"PCET_A_PCET_RTCR"] - pc_mean)/pc_sd
}
colnames(df)[322] <- "PCET_speed"

pc_sd <- sd(df$SPCPTNL_SCPT_TPRT,na.rm = T)
pc_mean <- mean(df$SPCPTNL_SCPT_TPRT,na.rm = T)
for (i in 1:length(df$SPCPTNL_SCPT_TPRT)){
  df[i,323] <- (df[i,"SPCPTNL_SCPT_TPRT"] - pc_mean)/pc_sd 
}
colnames(df)[323] <- "CPT_speed"

pc_sd <- sd(df$PCET_A_PCET_ACC2,na.rm = T)
pc_mean <- mean(df$PCET_A_PCET_ACC2,na.rm = T)
for (i in 1:length(df$PCET_A_PCET_ACC2)){
  df[i,324] <- (df[i,"PCET_A_PCET_ACC2"] - pc_mean)/pc_sd
}
colnames(df)[324] <- "PCET_acc"

pc_sd <- sd(df$SPCPTNL_SCPT_TP,na.rm = T)
pc_mean <- mean(df$SPCPTNL_SCPT_TP,na.rm = T)
for (i in 1:length(df$SPCPTNL_SCPT_TP)){
  df[i,325] <- (df[i,"SPCPTNL_SCPT_TP"] - pc_mean)/pc_sd
}
colnames(df)[325] <- "CPT_acc"



###efficiency columns
for(i in 1:length(df$CPF_A_CPF_CR)){ 
  df[i,326] <- (df[i,"CPF_A_CPF_CR"] + -1*df[i,"CPF_A_CPF_RTCR"])/2
  df[i,327] <- (df[i,"PCET_acc"] + -1*df[i,"PCET_speed"])/2
  df[i,328] <- (df[i,"CPT_acc"] + -1*df[i,"CPT_speed"])/2
  df[i,329] <- (df[i,"ER40_A_ER40_CR"] + -1*df[i,"ER40_A_ER40_RTCR"])/2
  df[i,330] <- (df[i,"PMAT24_A_PMAT24_A_CR"] + -1*df[i,"PMAT24_A_PMAT24_A_RTCR"])/2
  df[i,331] <- (df[i,"SFNB2_SFNB_MCR"] + -1*df[i,"SFNB2_SFNB_MRTC"])/2
  df[i,332] <- (df[i,"SVOLT_A_SVOLT_CR"] + -1*df[i,"SVOLT_A_SVOLT_RTCR"])/2
  df[i,333] <- (df[i,"VSPLOT24_VSPLOT24_CR"] + -1*df[i,"VSPLOT24_VSPLOT24_RTCR"])/2
}
for(i in 1:1){
  colnames(df)[326] <- "CPF_efficiency"
  colnames(df)[327] <- "PCET_efficiency"
  colnames(df)[328] <- "CPT_efficiency"
  colnames(df)[329] <- "ER40_efficiency"
  colnames(df)[330] <- "PMAT24_efficiency"
  colnames(df)[331] <- "SFNB2_efficiency"
  colnames(df)[332] <- "SVOLT_efficiency"
  colnames(df)[333] <- "VSPLOT_efficiency"
}

x <- df %>% 
  pivot_longer(cols  = c(326:333), names_to = c('test')) 
x <- x[which(x$value <= 5),]
x <- x[which(x$value >= -4),]
x$value = as.numeric(as.character(x$value))

temp <- x[(x$CASE.CONTROL == "case" | x$CASE.CONTROL =="control"),]

tgc <- summarySE(temp, measurevar = "value", groupvars = c("test","CASE.CONTROL","gender"),na.rm=T)
pd <- position_dodge(0.30)
ggplot(tgc, aes(x=test, y=value,col = gender, shape = CASE.CONTROL,group = interaction(gender,CASE.CONTROL))) + 
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.4, position=pd) +
  geom_line(position=pd, size = .27) +
  geom_point(position=pd, size = 6) +scale_shape_manual(name = "Case/Control",values = c(8,5))+
  scale_x_discrete(labels=c("CPF", "CPT", "ER40","PCET","PMAT24","SFNB2","SVOLT","VSPLOT"))+
  scale_colour_manual(name = "Gender", values = c("red","blue"))+ggtitle("CNB Xhosa - Efficiency by Case Type/Gender")+
  xlab("Test Name") + ylab("Z score (Efficiency)")+theme(legend.justification = "top")+
  geom_hline(yintercept=0.00, size = .25)

temp$test <- as.factor(temp$test)
temp$CASE.CONTROL <- droplevels(temp$CASE.CONTROL)
temp$CASE.NO <- as.factor(temp$CASE.NO)

me2 <- lmer(value ~ CASE.CONTROL + age + gender + test + (1|CASE.NO), data = temp)
me <- lmer(value ~ (CASE.CONTROL+gender+age+test)^2 + (1 | CASE.NO), data = temp)
summary(lmer(value ~ CASE.CONTROL+gender+age+test + test:CASE.CONTROL+test:gender+test:age +(1 | CASE.NO), data = temp))
summary(me)

me2 <- lmer(value ~ CASE.CONTROL*gender + age + test + (1|CASE.NO), data = temp)
visreg(me2,xvar="CASE.CONTROL", by = "gender",points=list(cex=.2, pch=5),jitter = T, overlay = T, xlab = "Case/Control",ylab = "Zscore (efficiency)")#case/gender/zscore
me2 <- lmer(value ~ CASE.CONTROL*test + age + gender + (1|CASE.NO), data = temp)
visreg(me2,xvar="test", by = "CASE.CONTROL",points=list(cex=.2, pch=5,col = c("green","orange")),line=list(col=c("green","orange")), breaks = 3,jitter = T, overlay = T, xlab = "Test",ylab = "Zscore (efficiency)")#test/case/zscore
me2 <- lmer(value ~ CASE.CONTROL*age+test + gender + (1|CASE.NO), data = temp)
visreg(me2,xvar="CASE.CONTROL", by = "age",points=list(cex=.2, pch=5), breaks = 3,jitter = T, overlay = T, xlab = "Case/Control",ylab = "Zscore (efficiency)")#case/age/zscore
me2 <- lmer(value ~ CASE.CONTROL + test*age + gender + (1|CASE.NO), data = temp)
visreg(me2,xvar="test", by = "age", breaks = 3,points=list(cex=.2, pch=5),jitter = T, overlay = T, xlab = "Test",ylab = "Zscore (efficiency)")#test/age/zscore

for(i in 1:length(x$age)){
  if(x[i,"age"] %in% 9:25){
    x[i,"age"] <- "9-25"
  }
  if(x[i,"age"] %in% 26:33){
    x[i,"age"] <- "26-33"
  }
  if(x[i,"age"] %in% 34:41){
    x[i,"age"] <- "34-41"
  }
  if(x[i,"age"] %in% 42:49){
    x[i,"age"] <- "42-49"
  }
  if(x[i,"age"] %in% 50:99){
    x[i,"age"] <- "50+"
  }
  if(x[i,"age"] %in% 0:8){
    x[i,"age"] <- NA
  }
}
temp <- x[(x$CASE.CONTROL == "case" | x$CASE.CONTROL =="control"),]
tgc <- summarySE(temp, measurevar = "value", groupvars = c("age","CASE.CONTROL","gender"), na.rm = T)
#tgc <- tgc[1:24,]
pd <- position_dodge(0.30)
x$age <-as.factor(x$age)
tgc$age <- factor(tgc$age, levels = c("9-25","26-33","34-41","42-49","50+"))
ggplot(tgc, aes(x=age, y=value,col = gender, shape = CASE.CONTROL,group = interaction(gender,CASE.CONTROL))) + 
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.4, position=pd) +
  geom_line(position=pd, size = .23) +
  geom_point(position=pd, size = 6) +scale_shape_manual(name = "Case/Control",values = c(8,5))+
  scale_x_discrete(labels=c("9-25","26-33","34-41","42-49","50+"))+
  scale_colour_manual(name = "Gender", values = c("red","blue"))+ggtitle("CNB Xhosa - Efficiency by Case Type/Gender(No test)")+
  xlab("Age") + ylab("Z score (efficiency)")+theme(legend.justification = "top")+
  geom_hline(yintercept=0.00, size = .25)+coord_cartesian(ylim=c(-1.5,1))+ scale_y_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1))


#####
#####
#####
#####
#Speed Factor Analysis
rm(list=ls())
df <- read.csv("xhosa_cnb_merged_clinical.csv")
df <- df[which(df$valid_code == "V" | df$valid_code == "VC" | df$valid_code == "V1"),]

for (i in 1:length(df$CPF_A_valid_code)){
  if(df[i,"CPF_A_valid_code"] != "V"){
    df[i,"CPF_A_CPF_RTCR"] <- NA
  }
  if(df[i,"ER40_A_valid_code"] != "V"){
    df[i,"ER40_A_ER40_RTCR"] <- NA
  }
  if(df[i,"PMAT24_A_valid_code"] != "V"){
    df[i,"PMAT24_A_PMAT24_A_RTCR"] <- NA
  }
  if(df[i,"SPCPTNL_valid_code"] != "V"){
    df[i,"SPCPTNL_SCPT_TPRT"] <- NA
  }
  if(df[i,"VSPLOT24_valid_code"] != "V"){
    df[i,"VSPLOT24_VSPLOT24_RTCR"] <- NA
  }
  if(df[i,"SVOLT_A_valid_code"] != "V"){
    df[i,"SVOLT_A_SVOLT_RTCR"] <- NA
  }
  if(df[i,"SFNB2_valid_code"] != "V"){
    df[i,"SFNB2_SFNB_MRTC"] <- NA
  }
  if(df[i,"MPRACT_valid_code"] != "V"){
    df[i,"MPRACT_MP2RTCR"] <- NA
  }
  if(df[i,"SCTAP_valid_code"] != "V"){
    df[i,"SCTAP.SCTAP_TOT"] <- NA
  }
  if(df[i,"PCET_A_valid_code"] != "V"){
    df[i,"PCET_A_PCET_RTCR"] <- NA
  }
}
df <- df[which(df$CASE.CONTROL == "case" | df$CASE.CONTROL =="control"),]

pc_sd <- sd(df$PCET_A_PCET_RTCR,na.rm = T)
pc_mean <- mean(df$PCET_A_PCET_RTCR,na.rm = T)
for (i in 1:length(df$PCET_A_PCET_RTCR)){
  df[i,322] <- (df[i,"PCET_A_PCET_RTCR"] - pc_mean)/pc_sd
}
colnames(df)[322] <- "PCET"

pc_sd <- sd(df$SPCPTNL_SCPT_TPRT,na.rm = T)
pc_mean <- mean(df$SPCPTNL_SCPT_TPRT,na.rm = T)
for (i in 1:length(df$SPCPTNL_SCPT_TPRT)){
  df[i,323] <- (df[i,"SPCPTNL_SCPT_TPRT"] - pc_mean)/pc_sd  
}
colnames(df)[323] <- "CPT"

pc_sd <- sd(df$MPRACT_MP2RTCR,na.rm = T)
pc_mean <- mean(df$MPRACT_MP2RTCR,na.rm = T)
for (i in 1:length(df$MPRACT_MP2RTCR)){
  df[i,324] <- (df[i,"MPRACT_MP2RTCR"] - pc_mean)/pc_sd  
}
colnames(df)[324] <- "MPRACT"

efa_RT <- data.frame(df$CPF_A_CPF_RTCR,df$ER40_A_ER40_RTCR,df$PCET,df$CPT,df$PMAT24_A_PMAT24_A_RTCR,df$SVOLT_A_SVOLT_RTCR,
                     df$VSPLOT24_VSPLOT24_RTCR,df$SFNB2_SFNB_MRTC,df$MPRACT,df$SCTAP_SCTAP_TOT)

colnames(efa_RT) <- c("CPF","ER40","PCET","CPT","PMAT24","SVOLT","VSPLOT","SFNB2","MPRACT","SCTAP")
nfactors(efa_RT)
fa.parallel(efa_RT)
fac <- fa(efa_RT,nfactors = 3, rotate = "oblimin", fm = "pa")
fac
print(fac$loadings, cutoff = 0.3)
fa.diagram(fac, sort = T,main = "CNB Xhosa - Speed Factor Analysis") 


######
######
######
#Accuracy FA

rm(list=ls())
df <- read.csv("xhosa_cnb_merged_clinical.csv")
df <- df[which(df$valid_code == "V" | df$valid_code == "VC" | df$valid_code == "V1"),]

for (i in 1:length(df$PCET_A_valid_code)){
  if(df[i,"PCET_A_valid_code"] != "V"){
    df[i,"PCET_A_PCET_ACC2"] <- NA
  }
  if(df[i,"CPF_A_valid_code"] != "V"){
    df[i,"CPF_A_CPF_CR"] <- NA
  }
  if(df[i,"ER40_A_valid_code"] != "V"){
    df[i,"ER40_A_ER40_CR"] <- NA
  }
  if(df[i,"PMAT24_A_valid_code"] != "V"){
    df[i,"PMAT24_A_PMAT24_A_CR"] <- NA
  }
  if(df[i,"SPCPTNL_valid_code"] != "V"){
    df[i,"SPCPTNL.SCPT_TP"] <- NA
  }
  if(df[i,"VSPLOT24_valid_code"] != "V"){
    df[i,"VSPLOT24_VSPLOT24_CR"] <- NA
  }
  if(df[i,"SVOLT_A_valid_code"] != "V"){
    df[i,"SVOLT_A_SVOLT_CR"] <- NA
  }
  if(df[i,"SFNB2_valid_code"] != "V"){
    df[i,"SFNB2_SFNB_MCR"] <- NA
  }
}
df <- df[which(df$CASE.CONTROL == "case" | df$CASE.CONTROL =="control"),]

pc_sd <- sd(df$PCET_A_PCET_ACC2,na.rm = T)
pc_mean <- mean(df$PCET_A_PCET_ACC2,na.rm = T)
for (i in 1:length(df$PCET_A_PCET_ACC2)){
  df[i,322] <- (df[i,"PCET_A_PCET_ACC2"] - pc_mean)/pc_sd
}
colnames(df)[322] <- "PCET"

pc_sd <- sd(df$SPCPTNL_SCPT_TP,na.rm = T)
pc_mean <- mean(df$SPCPTNL_SCPT_TP,na.rm = T)
for (i in 1:length(df$SPCPTNL_SCPT_TP)){
  df[i,323] <- (df[i,"SPCPTNL_SCPT_TP"] - pc_mean)/pc_sd
}
colnames(df)[323] <- "CPT"

efa_CR <- data.frame(df$CPF_A_CPF_CR, df$ER40_A_ER40_CR, df$PCET,df$PMAT24_A_PMAT24_A_CR, df$SVOLT_A_SVOLT_CR, df$VSPLOT24_VSPLOT24_CR,
                     df$CPT, df$SFNB2_SFNB_MCR)

colnames(efa_CR) <- c("CPF","ER40","PCET","PMAT","SVOLT","VSPLOT","CPT","SFNB2")
nfactors(efa_CR)
fa.parallel(efa_CR)
fac <- fa(efa_CR,nfactors = 3, rotate = "cluster", fm = "pa")
fac
print(fac$loadings, cutoff = 0.25)
fa.diagram(fac, sort = T, main = "CNB Xhosa - Accuracy Factor Analysis") 
##
###
#efficiency

rm(list=ls())
df <- read.csv("xhosa_cnb_merged_clinical.csv")
df <- df[which(df$valid_code == "V" | df$valid_code == "VC" | df$valid_code == "V1"),]

for (i in 1:length(df$PCET_A_valid_code)){
  if(df[i,"PCET_A_valid_code"] != "V"){
    df[i,"PCET_A_PCET_ACC2"] <- NA
    df[i,"PCET_A_PCET_RTCR"] <- NA
  }
  if(df[i,"CPF_A_valid_code"] != "V"){
    df[i,"CPF_A_CPF_CR"] <- NA
    df[i,"CPF_A_CPF_RTCR"] <- NA
  }
  if(df[i,"ER40_A_valid_code"] != "V"){
    df[i,"ER40_A_ER40_CR"] <- NA
    df[i,"ER40_A_ER40_RTCR"] <- NA
  }
  if(df[i,"PMAT24_A_valid_code"] != "V"){
    df[i,"PMAT24_A_PMAT24_A_CR"] <- NA
    df[i,"PMAT24_A_PMAT24_A_RTCR"] <- NA
  }
  if(df[i,"SPCPTNL_valid_code"] != "V"){
    df[i,"SPCPTNL.SCPT_TP"] <- NA
    df[i,"SPCPTNL_SCPT_TPRT"] <- NA
  }
  if(df[i,"VSPLOT24_valid_code"] != "V"){
    df[i,"VSPLOT24_VSPLOT24_CR"] <- NA
    df[i,"VSPLOT24_VSPLOT24_RTCR"] <- NA
  }
  if(df[i,"SVOLT_A_valid_code"] != "V"){
    df[i,"SVOLT_A_SVOLT_CR"] <- NA
    df[i,"SVOLT_A_SVOLT_RTCR"] <- NA
  }
  if(df[i,"SFNB2_valid_code"] != "V"){
    df[i,"SFNB2_SFNB_MCR"] <- NA
    df[i,"SFNB2_SFNB_MRTC"] <- NA
  }
  
}
df <- df[which(df$CASE.CONTROL == "case" | df$CASE.CONTROL =="control"),]

pc_sd <- sd(df$PCET_A_PCET_RTCR,na.rm = T)
pc_mean <- mean(df$PCET_A_PCET_RTCR,na.rm = T)
for (i in 1:length(df$PCET_A_PCET_RTCR)){
  df[i,322] <- (df[i,"PCET_A_PCET_RTCR"] - pc_mean)/pc_sd
}
colnames(df)[322] <- "PCET_speed"

pc_sd <- sd(df$SPCPTNL_SCPT_TPRT,na.rm = T)
pc_mean <- mean(df$SPCPTNL_SCPT_TPRT,na.rm = T)
for (i in 1:length(df$SPCPTNL_SCPT_TPRT)){
  df[i,323] <- (df[i,"SPCPTNL_SCPT_TPRT"] - pc_mean)/pc_sd 
}
colnames(df)[323] <- "CPT_speed"

pc_sd <- sd(df$PCET_A_PCET_ACC2,na.rm = T)
pc_mean <- mean(df$PCET_A_PCET_ACC2,na.rm = T)
for (i in 1:length(df$PCET_A_PCET_ACC2)){
  df[i,324] <- (df[i,"PCET_A_PCET_ACC2"] - pc_mean)/pc_sd
}
colnames(df)[324] <- "PCET_acc"

pc_sd <- sd(df$SPCPTNL_SCPT_TP,na.rm = T)
pc_mean <- mean(df$SPCPTNL_SCPT_TP,na.rm = T)
for (i in 1:length(df$SPCPTNL_SCPT_TP)){
  df[i,325] <- (df[i,"SPCPTNL_SCPT_TP"] - pc_mean)/pc_sd
}
colnames(df)[325] <- "CPT_acc"

###efficiency columns
for(i in 1:length(df$CPF_A_CPF_CR)){ 
  df[i,326] <- (df[i,"CPF_A_CPF_CR"] + -1*df[i,"CPF_A_CPF_RTCR"])/2
  df[i,327] <- (df[i,"PCET_acc"] + -1*df[i,"PCET_speed"])/2
  df[i,328] <- (df[i,"CPT_acc"] + -1*df[i,"CPT_speed"])/2
  df[i,329] <- (df[i,"ER40_A_ER40_CR"] + -1*df[i,"ER40_A_ER40_RTCR"])/2
  df[i,330] <- (df[i,"PMAT24_A_PMAT24_A_CR"] + -1*df[i,"PMAT24_A_PMAT24_A_RTCR"])/2
  df[i,331] <- (df[i,"SFNB2_SFNB_MCR"] + -1*df[i,"SFNB2_SFNB_MRTC"])/2
  df[i,332] <- (df[i,"SVOLT_A_SVOLT_CR"] + -1*df[i,"SVOLT_A_SVOLT_RTCR"])/2
  df[i,333] <- (df[i,"VSPLOT24_VSPLOT24_CR"] + -1*df[i,"VSPLOT24_VSPLOT24_RTCR"])/2
}
for(i in 1:1){
  colnames(df)[326] <- "CPF_efficiency"
  colnames(df)[327] <- "PCET_efficiency"
  colnames(df)[328] <- "CPT_efficiency"
  colnames(df)[329] <- "ER40_efficiency"
  colnames(df)[330] <- "PMAT24_efficiency"
  colnames(df)[331] <- "SFNB2_efficiency"
  colnames(df)[332] <- "SVOLT_efficiency"
  colnames(df)[333] <- "VSPLOT_efficiency"
}

efa_CR <- data.frame(df$CPF_efficiency, df$ER40_efficiency, df$PCET_efficiency,df$PMAT24_efficiency, df$SVOLT_efficiency, df$VSPLOT_efficiency,
                     df$CPT_efficiency, df$SFNB2_efficiency)

nfactors(efa_CR)
fa.parallel(efa_CR)
fac <- fa(efa_CR,nfactors = 4, rotate = "oblimin", fm = "pa")
fac
print(fac$loadings, cutoff = 0.25)
fa.diagram(fac, sort = T, main = "CNB Xhosa - Efficiency Factor Analysis") 
str(x$CASE.CONTROL)
str(x$test)


x$test <- as.factor(x$test)
x$gender <- as.factor(x$gender)
x$CASE.CONTROL<- as.factor(x$CASE.CONTROL)

summary(lme(value~ test +(gender+age+CASE.CONTROL),data=x,random = ~ 1 | valid_code, na.action = na.exclude))

