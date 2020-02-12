plot_20191109

#Plot

library(psych)
library(Amelia)
library(ggplot2)
library(CorrMixed)
library(gtools)
library(mgcv)
library(visreg)
library(lubridate)
library(irr)

df <- read.csv("plot_20191109.csv") 

df[,6] <- as.numeric(mdy(df[,6]))
colnames(df)[6] <- "Time"

#df[,11] <- as.numeric(mdy(df[,11]))
#colnames(df)[11] <- "DOB"
colnames(df)[1] <- "bblid"
colnames(df)[13] <- "Sex"


ids <- unique(na.omit(df$bblid))
for (i in 1:length(ids)) {
  temp <- df[which(df$bblid == ids[i]),]
  if (is.element("N",temp[,c(8)])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
  df[which(df$bblid == ids[i]),] <- temp}

acc <- df$plot_pc
#PTP[PTP<0.3333] <- 0.3333
acc[acc<18] <- NA
df$plot_pc <- acc

spe <- df$plot_rtcr
#PTP[PTP<0.3333] <- 0.3333
spe[spe<200] <- NA
df$plot_rtcr <- spe

#df <- df[which(df$Timepoint > 1),]   # best is > 1
df <- df[which(df$timepoint < 6),]     # best is < 6
df <- df[which(is.na(df$age) == FALSE),]###

Age <- scale(df$age)
Age_Squared <- Age^2
Age_Cubed <- Age^3


set.seed(2)
temp <- amelia(df[,c("plot_pc","plot_rtcr")], m=1)$imputations[[1]]
df[,c("plot_pc","plot_rtcr")] <- temp
#

modx <- lm(plot_pc~Age+Age_Squared+Age_Cubed,data=df,na.action=na.exclude)
acc_r <- scale(residuals(modx,na.action=na.exclude))
#acc_r <- lm(CPF_A.CPF_CR~Age+Age_Squared+Age_Cubed,data=df,na.action=na.exclude)
mody <- lm(plot_rtcr~Age+Age_Squared+Age_Cubed,data=df, na.action = na.exclude)
spe_r <- scale(residuals(mody,na.action = na.exclude))

df <- data.frame(df,acc_r,spe_r)

ids <- unique(na.omit(df$bblid))
for (i in 1:length(ids)) {
  temp <- df[which(df$bblid == ids[i]),]
  temp$Time <- temp$Time - min(temp$Time,na.rm=TRUE)
  df[which(df$bblid == ids[i]),] <- temp}

Age_Median_Split <- as.numeric(quantcut(Age,q=2))

#rescaling age/time to years
df$Time <- df$Time/365.25
Time_Squared <- scale(df$Time)^2

df <- data.frame(df,Age_Median_Split,Time_Squared)

# Run actual models
summary(lme(plot_pc~Age_Squared+(Sex+age)^2,data=df,random = ~ 1 | bblid, na.action = na.exclude))
summary(lme(plot_rtcr~Age_Squared+(Sex+age)^2,data=df,random = ~ 1 | bblid,na.action = na.exclude))

ggplot(df, aes(x = age, y = plot_pc, group = factor(Sex), col = factor(Sex))) + stat_smooth(method = "loess", se = F)+
  xlab("Age (years)")+ylab("Percent Correct") + geom_point(size = .6) +ggtitle("Plot Test - Age vs Percent Correct") + labs(col = "Sex")+
  theme(legend.justification = "top")+theme(plot.title = element_text(hjust=0.5))+ them

#fits a polynomial using local fitting
ggplot(df, aes(x = age, y = plot_rtcr, group = factor(Sex), col = factor(Sex))) + stat_smooth(method = "loess", se = F)+
  xlab("Age (years)")+ylab("Reaction Time (ms)") + geom_point(size = .5) +ggtitle("Plot Test - Age vs Reaction Time") + labs(col = "Sex")+
  coord_cartesian(ylim=c(0,25000))+
  scale_y_continuous(name = "Reaction Time(ms)",breaks = c(5000,10000,15000,20000,25000))+ theme(legend.justification = "top")+
  theme(plot.title = element_text(hjust = 0.5)) 




#visreg(lme(CPF_A_CR~Age_Squared+(Sex+test_sessions_v.age)^2,data=df,random = ~ 1 | bblid,na.action = na.exclude),xvar="test_sessions_v.age",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Score")
#visreg(lme(CPF_A.CPF_RT~Age_Squared+(Sex+test_sessions_v.age)^2,data=df,random = ~ 1 | bblid,na.action = na.exclude),xvar="test_sessions_v.age",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Reaction Time (ms)")

visreg(lme(plot_pc~Age_Squared+(Sex+age)^2,data=df,random = ~ 1 | bblid,na.action = na.exclude),xvar="age",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Score")
visreg(lme(plot_rtcr~Age_Squared+(Sex+age)^2,data=df,random = ~ 1 | bblid,na.action = na.exclude),xvar="age",by="Sex",overlay=TRUE,xlab="Age (Years)",ylab="Reaction Time (ms)")




#################wide
#################wide
#################wide

df <- read.csv("plot_20191109.csv") 

df[,6] <- as.numeric(mdy(df[,6]))
colnames(df)[6] <- "Time"

#df[,11] <- as.numeric(mdy(df[,11]))
#colnames(df)[11] <- "DOB"
colnames(df)[1] <- "bblid"
colnames(df)[13] <- "Sex"

ids <- unique(df$bblid)

ids <- unique(na.omit(df$bblid))
for (i in 1:length(ids)) {
  temp <- df[which(df$bblid == ids[i]),]
  if (is.element("N",temp[,c(8)])) {temp <- matrix(NA,dim(temp)[1],dim(temp)[2])}
  df[which(df$bblid == ids[i]),] <- temp}

df <- df[which(is.na(df$bblid) == FALSE),]

acc <- df$plot_pc
#PTP[PTP<0.3333] <- 0.3333
acc[acc<18] <- NA
df$plot_pc <- acc

spe <- df$plot_rtcr
#PTP[PTP<0.3333] <- 0.3333
spe[spe<200] <- NA
df$plot_rtcr <- spe

#df <- df[which(df$Timepoint > 1),]   # best is > 1
df <- df[which(df$timepoint < 6),]     # best is < 6 , there are no x6 or x7 variables 
df <- df[which(is.na(df$age) == FALSE),]###

Age <- scale(df$age)
Age_Squared <- Age^2
Age_Cubed <- Age^3

TP <- scale(df$timepoint)
TP_Squared <- scale(df$timepoint)^2
TP_Cubed <- scale(df$timepoint)^3

#impute missing

set.seed(2)
temp <- amelia(df[,c("plot_pc","plot_rtcr")], m=1)$imputations[[1]]
df[,c("plot_pc","plot_rtcr")] <- temp
ids <- unique(df$bblid)

for (i in 1:length(ids)) {
  temp <- df[which(df$bblid == ids[i]),]
  temp$Time <- temp$Time - min(temp$Time,na.rm=TRUE)
  df[which(df$bblid == ids[i]),] <- temp}

int <- matrix(NA,dim(df)[1],1)
df <- data.frame(df,int)

#create interval between time points
for (i in 1:length(ids)) {
  temp <- df[which(df$bblid == ids[i]),]
  for (j in 1:dim(temp)[1]) {
    if (temp[j,2] == 8) {try(temp[j,20] <- temp[j,6] - temp[which(temp$timepoint == 7),6])}
    if (temp[j,2] == 7) {try(temp[j,20] <- temp[j,6] - temp[which(temp$timepoint == 6),6])}
    if (temp[j,2] == 6) {try(temp[j,20] <- temp[j,6] - temp[which(temp$timepoint == 5),6])}
    if (temp[j,2] == 5) {try(temp[j,20] <- temp[j,6] - temp[which(temp$timepoint == 4),6])}
    if (temp[j,2] == 4) {try(temp[j,20] <- temp[j,6] - temp[which(temp$timepoint == 3),6])}
    if (temp[j,2] == 3) {try(temp[j,20] <- temp[j,6] - temp[which(temp$timepoint == 2),6])}
    if (temp[j,2] == 2) {try(temp[j,20] <- temp[j,6])}
    if (temp[j,2] == 1) {try(temp[j,20] <- 0)}
    df[which(df$bblid == ids[i]),] <- temp}}

#regressing out age
mod <- lm(plot_pc~Age+Age_Squared+Age_Cubed,data=df,na.action=na.exclude)
acc_r <- scale(residuals(mod,na.action=na.exclude))
mod2 <- lm(plot_rtcr~Age+Age_Squared+Age_Cubed,data=df,na.action=na.exclude)
spe_r <- scale(residuals(mod2,na.action=na.exclude))

#regressing out age and practice...had to remove int
mod3 <- lm(plot_pc~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed,data=df,na.action = na.exclude)
acc_r_p<- scale(residuals(mod3,na.action=na.exclude))
mod4<- lm(plot_rtcr~Age+Age_Squared+Age_Cubed+TP+TP_Squared+TP_Cubed,data=df,na.action = na.exclude)
spe_r_p <- scale(residuals(mod4,na.action=na.exclude))


df<- data.frame(df,acc_r,acc_r_p,spe_r,spe_r_p)


x1 <- df[which(df$timepoint == 1),]
x2 <- df[which(df$timepoint == 2),]
x3 <- df[which(df$timepoint == 3),]
x4 <- df[which(df$timepoint == 4),]
x5 <- df[which(df$timepoint == 5),]
x6 <- df[which(df$timepoint == 6),]
x7 <- df[which(df$timepoint == 7),]

x <- merge(x1,x2,by="bblid",all=TRUE)
x <- merge(x,x3,by="bblid",all=TRUE)
x <- merge(x,x4,by="bblid",all=TRUE)
x <- merge(x,x5,by="bblid",all=TRUE)
x <- merge(x,x6,by="bblid",all=TRUE)
x <- merge(x,x7,by="bblid",all=TRUE)

colnames(x)[21:24] <- paste(colnames(x)[21:24],"_TIME_1",sep="")
colnames(x)[44:47] <- paste(colnames(x)[44:47],"_TIME_2",sep="")
colnames(x)[67:70] <- paste(colnames(x)[67:70],"_TIME_3",sep="")
colnames(x)[90:93] <- paste(colnames(x)[90:93],"_TIME_4",sep="")#
colnames(x)[113:116] <- paste(colnames(x)[113:116],"_TIME_5",sep="")#

acc_age_out <- data.frame(x$acc_r.x_TIME_1,x$acc_r.y_TIME_2,x$acc_r.x_TIME_3,x$acc_r.y_TIME_4,x$acc_r.x_TIME_5)
acc_age_prac_out <- data.frame(x$acc_r_p.x_TIME_1,x$acc_r_p.y_TIME_2,x$acc_r_p.x_TIME_3,x$acc_r_p.y_TIME_4,x$acc_r_p.x_TIME_5)
#speed variable
speed_age_out <- data.frame(x$spe_r.x_TIME_1,x$spe_r.y_TIME_2,x$spe_r.x_TIME_3,x$spe_r.y_TIME_4,x$spe_r.x_TIME_5)
speed_age_prac_out <- data.frame(x$spe_r_p.x_TIME_1,x$spe_r_p.y_TIME_2,x$spe_r_p.x_TIME_3,x$spe_r_p.y_TIME_4,x$spe_r_p.x_TIME_5)
#age and practice regressed out

res <- matrix(NA,4,4)



for (i in 2:5) {
  res[1,(i-1)] <- icc(acc_age_out[,1:i],type="agreement",model="twoway")$value
  res[2,(i-1)] <- icc(acc_age_prac_out[,1:i],type="agreement",model="twoway")$value
  res[3,(i-1)] <- icc(speed_age_out[,1:i],type="agreement",model="twoway")$value
  res[4,(i-1)] <- icc(speed_age_prac_out[,1:i],type="agreement",model="twoway")$value
} 
res






df4 <- df[which (df$Time < 2000),]
#df4 <- df4[-which (df4$Time > 0 & df3$Time <730),]
df4 <- df4[which(df4$ntimepoints >= 2) ,]
df4 <- df4[which(df4$timepoint <= 2),] 
df4 <- subset(df4,duplicated(bblid) | duplicated(bblid, fromLast=TRUE))
i <- 1

for (i in 1:length(df4$bblid)){
  scoref <- df4[i*2,18]
  scorepr <- df4[2*i-1,18]
  days <- df4[i*2, 20]
  df4[i,21] <- days
  df4[i,22] <- scoref - scorepr
}

for (i in 1:length(df4$bblid)){
  if( df4[i,21] <366){
    df4[i,23] <- "0 to 1 Year"
  }
  if( df4[i,21] > 365 & df4[i,21]< 730){
    df4[i,23] <- "1 to 2 Years"
  }
  if( df4[i,21] > 730 & df4[i,21]< 1095){
    df4[i,23] <- "2 to 3 Years"
  }
  if( df4[i,21] > 1095 & df4[i,21]< 1460){
    df4[i,23] <- "3 to 4 Years"
  }
  if( df4[i,21] > 1460 & df4[i,21]< 1825){
    df4[i,23] <- "4 to 5 Years"
  }
  if( df4[i,21] > 1824){
    df4[i,23] <- "5+ Years"
  }
}
#df4 <- df4[which(is.na(df4[,23]) == FALSE),]







x <-1
summarydf <- as.data.frame(df4$bblid)
for (x in 1:length(df4[,23])){
  if (df4[x,23] == "0 to 1 Year"){
    summarydf[x,2] <- df4[x,22]
  }
  
  
}
summarydf <- summarydf[which(is.na(summarydf[,2]) == FALSE),]
finaldf <- as.data.frame(mean(summarydf[,2]))
finaldf[1,2] <- sd(summarydf[,2]) / sqrt(length(summarydf[,2]))
finaldf [1,3] <- "0 to 1 Year"

##
x <-1
summarydf <- as.data.frame(df4$bblid)
for (x in 1:length(df4[,23])){
  if (df4[x,23] == "1 to 2 Years"){
    summarydf[x,2] <- df4[x,22]
    
  }
  if (df4[x,23] == "1 to 2 Years"){
    summarydf[x,2] <- df4[x,22]
  }
  
  
}
summarydf <- summarydf[which(is.na(summarydf[,2]) == FALSE),]
finaldf[2,1] <- mean(summarydf[,2])
finaldf[2,2] <- sd(summarydf[,2]) / sqrt(length(summarydf[,2]))
finaldf [2,3] <- "1 to 2 Years"
##
##
##
x <-1
summarydf <- as.data.frame(df4$bblid)
for (x in 1:length(df4[,23])){
  if (df4[x,23] == "2 to 3 Years"){
    summarydf[x,2] <- df4[x,22]
    
  }
  
  
}
summarydf <- summarydf[which(is.na(summarydf[,2]) == FALSE),]
finaldf[3,1] <- mean(summarydf[,2])
finaldf[3,2] <- sd(summarydf[,2]) / sqrt(length(summarydf[,2]))
finaldf [3,3] <- "2 to 3 Years"
##
x <-1
summarydf <- as.data.frame(df4$bblid)
for (x in 1:length(df4[,23])){
  if (df4[x,23] == "3 to 4 Years"){
    summarydf[x,2] <- df4[x,22]
  }
  
}
summarydf <- summarydf[which(is.na(summarydf[,2]) == FALSE),]
finaldf[4,1] <- mean(summarydf[,2])
finaldf[4,2] <- sd(summarydf[,2]) / sqrt(length(summarydf[,2]))
finaldf [4,3] <- "3 to 4 Years"
##
x <-1
summarydf <- as.data.frame(df4$bblid)
for (x in 1:length(df4[,23])){
  if (df4[x,23] == "5+ Years"){
    summarydf[x,2] <- df4[x,22]
  }
  
}
summarydf <- summarydf[which(is.na(summarydf[,2]) == FALSE),]
finaldf[5,1] <- mean(summarydf[,2])
finaldf[5,2] <- sd(summarydf[,2]) / sqrt(length(summarydf[,2]))
finaldf [5,3] <- "5+ Years"
colnames(finaldf)[3] <- "group"
finaldf$group <- as.factor(finaldf$group)
colnames(finaldf)[1] <- "mean"
colnames(finaldf)[2] <- "se"

Time_Period <- as.factor(df4[,23])
ggplot( data = df4, mapping = aes(Time_Period, df4[,22], fill = Time_Period)) + geom_boxplot()+ 
  xlab("Days Between Time1 & Time 2") + ylab("Score Difference(out of 40)") +  geom_errorbar(aes(ymin=finaldf$mean - finaldf$sd, finaldf$mean + finaldf$sd), width=.2,
                                                                                             position=position_dodge(.9)) 


ggplot(data = finaldf, mapping = aes(finaldf$group, finaldf$mean)) +geom_bar(stat="identity")+
  xlab("Days Between Time1 & Time 2") + ylab("Score Difference(Percentage)")+
  geom_errorbar(aes(ymin=finaldf$mean - 2*finaldf$se, ymax= finaldf$mean + 2*finaldf$se), width=.2,
                position=position_dodge(.9)) 

ggplot(dtata = finaldf,)




























df4 <- df4[1:563,]
Time_Period <- as.factor(df4[,23])
ggplot( data = df4, mapping = aes(Time_Period, df4[,22], fill = Time_Period)) + geom_boxplot()+ xlab("Days Between Time1 & Time 2") + ylab("Score Difference(out of 40)") + title("hi")
#group and set as factor


sum3 <- df[which(df$timepoint %in% "3"),]