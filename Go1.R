#Go1 Data
library(paramap)
library(paran)
library(polycor)
library(psych)
library(mokken)
library(Rmisc)
library(dplyr)
library(tidyr)
#ibrary(devtools)


#library(dplyr)
rm(list=ls())
df <- read.csv("GO1_medical.csv")
demographic <- read.csv("GO1_demographics.csv")
df <- merge(df,demographic)

df <- df[,-17, drop = F]#deletes huntingtons_rtg
df <- df[,-25, drop = F]#deletes med_ms_rt, 
df$smry_med_other_vision_rtg <- as.numeric(df$smry_med_other_vision_rtg)

x <- df %>% 
  pivot_longer(cols  = c(2:42), names_to = c("name")) 
tgc_sex <- summarySE(x, measurevar = "value", groupvars = c("name","sex"), na.rm = T)
tgc_race <- summarySE(x, measurevar = "value", groupvars = c("name","race"), na.rm = T)#combined?
#tgc_race <- tgc[which(tgc$race == 1 | tgc$race == 2 |  tgc$race == 5),]
df<- df[,2:42]
df <- df[,-13]
df <- df[,-27]

cordf <- polychoric(df)$rho
#factors
#nfactors(cordf, n.obs = 9498)
#fa.parallel(df, cor = 'poly')#2-4
#fa.parallel(df, cor = 'poly', quant = 0.97) #3-4 glorfeld, 
#paran(mat = cordf, cfa = T, centile = 99, iterations = 2000, graph = T, n = 9498)#3 glorfeld, 
#VSS.scree(cordf)# 3
#paramap::MAP(cordf)#3 factors

fac <- fa.sort(fa(cordf,4,rotate="oblimin", fm = "ml", n.obs = 9498, cor = "poly"))
fac
print(fac$loadings, cutoff = 0.3)
fa.diagram(fac, sort = T, main = "Go1")



#bifactor, no correlations between latents traits, very low, bifactor structure does not look gooed, the generals are all low as wlel, you can rerun just the omega to see

omega(cordf,4,flip=FALSE, n.obs = 9498, fm = "ml")#
bifac <- omegaSem(cordf,4,flip = F, n.obs = 9498, fm = "ml")#ECV = 0.25, 25% of the common variance is attributable to the general factor
bifac

#
#Male PCA
rm(list=ls())
df <- read.csv("GO1_medical.csv")
demographic <- read.csv("GO1_demographics.csv")
df <- merge(df,demographic)

df <- df[,-17, drop = F]#deletes huntingtons_rtg
df <- df[,-25, drop = F]#deletes med_ms_rt, 
df$smry_med_other_vision_rtg <- as.numeric(df$smry_med_other_vision_rtg)
df_male <- df[which(df$sex == 1),]
df_male<- df_male[,2:42]
df_male <- df_male[,-13, drop = F]#move this before pivot
df_male <- df_male[,-28, drop = F]
cordf <- polychoric(df_male)$rho

nfactors(cordf,n = 4592)
VSS.scree(cordf)
pc <- pca(cordf, nfactors = 3, n.obs = 4592, cor = "poly", rotate = "promax")
print(pc$loadings, cutoff = 0.3)
plot(pc)
fa.diagram(pc)

#Female
rm(list=ls())
df <- read.csv("GO1_medical.csv")
demographic <- read.csv("GO1_demographics.csv")
df <- merge(df,demographic)

df <- df[,-17, drop = F]#deletes huntingtons_rtg
df <- df[,-25, drop = F]#deletes med_ms_rt, 
df$smry_med_other_vision_rtg <- as.numeric(df$smry_med_other_vision_rtg)
df_female <- df[which(df$sex == 2),]
df_female<- df_female[,2:42]
df_female <- df_female[,-13, drop = F]#move this before pivot
df_female <- df_female[,-28, drop = F]
cordf <- polychoric(df_female)$rho

nfactors(cordf,n = 4906)
pc <- pca(cordf, nfactors = 3, n.obs = 4906, cor = "poly", rotate = "promax")
print(pc$loadings, cutoff = 0.3)
plot(pc)
fa.diagram(pc)

#Not white
rm(list=ls())
df <- read.csv("GO1_medical.csv")
demographic <- read.csv("GO1_demographics.csv")
df <- merge(df,demographic)

df <- df[,-17, drop = F]#deletes huntingtons_rtg
df <- df[,-25, drop = F]#deletes med_ms_rt, 
df$smry_med_other_vision_rtg <- as.numeric(df$smry_med_other_vision_rtg)
df <- df[which(df$race != 1),]
df <- df[,2:42]
df <- df[,-13, drop = F]#move this before pivot
df <- df[,-28, drop = F]
cordf <- polychoric(df)$rho

nfactors(cordf, n.obs = 4200)

nfactors(cordf,n = 4200)
pc <- pca(cordf, nfactors = 3, n.obs = 4200, cor = "poly", rotate = "promax")
print(pc$loadings, cutoff = 0.3)
plot(pc)
fa.diagram(pc)
######

#should have one for FA

dffac3 <- df[,c(20,15,14,10,16,13,24,21,29,7,11,18)]#checking monotonicity
sum_temp <- 0
for(i in 1:length(zz$smry_med_metabolic_rtg)){
  for(x in 1:12){
    sum_temp <- sum_temp + zz[i,x]
  }
  zz[i,13] <- sum_temp
  sum_temp <- 0

}

#monotonicity, sep into 4 groups (based on PCA)
df1 <- df[,c(1,2,28)]
df2 <- df[,c(32,33,35,36,38,34)]
df3 <- df[,c(14,20,18,7,15,10,13,11,17,24,16)]#fac1
df4 <- df[,c(4,27,3,6,5,30,8,23,25)]
df5 <- df[,c(12,21,22,26,19,9,29,37,39)]

z <- check.monotonicity(na.omit(df5))
plot(z)
summary(z)


#####
#####PCA
#####
rm(list=ls())
df <- read.csv("GO1_medical.csv")
demographic <- read.csv("GO1_demographics.csv")
df <- merge(df,demographic)

df <- df[,-17, drop = F]#deletes huntingtons_rtg
df <- df[,-25, drop = F]#deletes med_ms_rt, 
df$smry_med_other_vision_rtg <- as.numeric(df$smry_med_other_vision_rtg)
df<- df[,2:42]
df <- df[,-13, drop = F]#deletes
df <- df[,-27, drop = F]#deletes
cordf <- polychoric(df)$rho

#paran(mat = cordf, cfa = F, centile = 99, iterations = 2000, graph = T, n = 9498)#

pc <- pca(cordf, nfactors = 5, n.obs = 9498, cor = "poly", rotate = "promax")
print(pc$loadings, cutoff = 0.3)
plot(pc)
fa.diagram(pc)
