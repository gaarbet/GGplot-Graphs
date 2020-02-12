#dif target vs foils(faces that were not there at the beginning)

library(mirt)
library(psych)
#CPF_B_I.9492.csv
x <- read.csv("CPF_A_I.9492.csv")
file <- read.csv("CPF_A_I.9492.csv")
dfx <- read.csv("CPF_repeated_measures_partial.csv")
z <- read.csv("GO_RACE_3.csv")

x <- na.omit(x[,grepl("_CORR",colnames(x))])

xcor <- polychoric(x)$rho

nfactors(xcor,n.obs=4000)
fa.parallel(xcor,n.obs=4000)
VSS.scree(xcor)

factors <- 2

mod <- fa(xcor,factors)$loadings[,1:factors]
loadingss <- fa(xcor,factors)$loadings[,1:factors]##

bifactor_key <- rep(NA,dim(x)[2])

#assigning items to target or foil factor (can disagree)
for (j in 1:dim(x)[2]) {
  if (mod[j,1] == max(mod[j,])) {bifactor_key[j] <- 1}
  if (mod[j,2] == max(mod[j,])) {bifactor_key[j] <- 2}
  try(if (mod[j,3] == max(mod[j,])) {bifactor_key[j] <- 3})
  try(if (mod[j,4] == max(mod[j,])) {bifactor_key[j] <- 4})
  try(if (mod[j,5] == max(mod[j,])) {bifactor_key[j] <- 5})
  try(if (mod[j,6] == max(mod[j,])) {bifactor_key[j] <- 6})
  try(if (mod[j,7] == max(mod[j,])) {bifactor_key[j] <- 7})
  try(if (mod[j,8] == max(mod[j,])) {bifactor_key[j] <- 8})
  try(if (mod[j,9] == max(mod[j,])) {bifactor_key[j] <- 9})
  try(if (mod[j,10] == max(mod[j,])) {bifactor_key[j] <- 10})
}

mod <- bfactor(x,bifactor_key, optimizer = "NR")

params <- coef(mod,simplify=TRUE)$items
rows <- dim(params)[1]
end1 <- dim(params)[2]-2
params <- params[,1:end1]/1.702
slopes <- params[,1:(end1-1)]
intercepts <- params[,end1]
constant <- sqrt(rowSums(slopes^2)+1)
difficulties <- (intercepts*(-1))/constant

params <- data.frame(slopes[,1],difficulties)  # NORMAL METRIC!!!!

#########**----------------------------Dif----------------------------------------**########


loadings <- bifactor_key # Correspond to screener DEP, other DEP, GAD and SOC items ####

results <- data.frame(matrix(NA, nrow=40, ncol=6)) 
colnames(results) <- c("Item", "Main", "Interaction", "DiffAIC", "ChiSq", "ChiSqP") 

#items <- na.omit(x[,grepl("_CORR",colnames(x))])

df <- merge(file, dfx, by.x = "SESSID", by.y = "test_sessions.datasetid")
df <- merge(df, z, by.x = "test_sessions.bblid", by.y = "bblid")#
#dfx is repeated measures file cause it includes gender/race
#temp <- na.omit(proband_dff[,grepl("_CORR",colnames(proband_dff))])#from the merged file
df$gender_real <- as.factor(df$gender_real)

items <- na.omit(grep("_CORR", colnames(df), value=TRUE))

i=1
for (item in items) {
  # Determine index of item for testing
  item_ind <- which(items %in% item)
  
  #need to get proband df (includes dif factors)
  tmp_cpf_mod <- bfactor(df[,items[!(items %in% item)]], loadings[-item_ind], technical=list(NCYCLES=2000))
  write.csv(summary(tmp_cpf_mod)[[1]], file=paste0("/Users/gaarbet/Documents/Filetest", item, "_ex.csv"), row.names=FALSE)
  
  assign(paste0(item, "_cpf_bifactor"), tmp_cpf_mod)
  # General = Face Memory/Working Memory
  scores <- fscores(tmp_cpf_mod, QMC=TRUE) ### Q: defaults okay?
  df[,paste0("General", item, "_ex")] <- scores[,1]
  df[,paste0("Target", item, "_ex")] <- scores[,2] 
  df[,paste0("Foil", item, "_ex")] <- scores[,3] 

  #General + Target + Foil + Gender + (General*Gender)
  tmp_log_mod <- glm(df[,item] ~ df[,paste0("General", item, "_ex")] + df[,paste0("Target", item, "_ex")] + df[,paste0("Foil", item, "_ex")] + df$gender_real + df[,paste0("General", item, "_ex")]*df$gender_real, family="binomial")
  assign(paste0(item, "_mod"), tmp_log_mod)
  #General + Target + Foil  
  tmp_simp_log_mod <- glm(df[,item] ~ df[,paste0("General", item, "_ex")] + df[,paste0("Target", item, "_ex")] + df[,paste0("Foil", item, "_ex")], family="binomial")
  assign(paste0(item, "_simp_mod"), tmp_simp_log_mod)
  
  # Run a ChiSq Test 
  chisq <- anova(tmp_simp_log_mod, tmp_log_mod, test="Chisq")
  
  # Put relevant statistics in results dataframe
  results[i, "Item"] <- item
  results[i, "Main"] <- summary(tmp_log_mod)$coefficients[4,3]
  results[i, "Interaction"] <- summary(tmp_log_mod)$coefficients[5,3]
  results[i, "DiffAIC"] <- summary(tmp_log_mod)$aic - summary(tmp_simp_log_mod)$aic
  results[i, "ChiSq"] <- chisq$Deviance[2] 
  results[i, "ChiSqP"] <- chisq[[5]][2] 
  i=i+1
}

###


Foil <- df[,c("CPFa_I_TRIAL000001_CORR", "CPFa_I_TRIAL000002_CORR",
                         "CPFa_I_TRIAL000007_CORR",'CPFa_I_TRIAL000011_CORR',
                         "CPFa_I_TRIAL000012_CORR",'CPFa_I_TRIAL000013_CORR',
                         "CPFa_I_TRIAL000015_CORR",'CPFa_I_TRIAL000016_CORR',
                         "CPFa_I_TRIAL000017_CORR",'CPFa_I_TRIAL000018_CORR',
                         "CPFa_I_TRIAL000023_CORR",'CPFa_I_TRIAL000024_CORR',
                         "CPFa_I_TRIAL000026_CORR",'CPFa_I_TRIAL000027_CORR',
                         "CPFa_I_TRIAL000030_CORR",'CPFa_I_TRIAL000032_CORR',
                         "CPFa_I_TRIAL000033_CORR",'CPFa_I_TRIAL000039_CORR',
                         "CPFa_I_TRIAL000038_CORR",'CPFa_I_TRIAL000037_CORR')]

Target <- df[,c("CPFa_I_TRIAL000003_CORR", "CPFa_I_TRIAL000004_CORR",
                                 "CPFa_I_TRIAL000005_CORR",'CPFa_I_TRIAL000006_CORR',
                                 "CPFa_I_TRIAL000008_CORR",'CPFa_I_TRIAL000009_CORR',
                                 "CPFa_I_TRIAL000010_CORR",'CPFa_I_TRIAL000014_CORR',
                                 "CPFa_I_TRIAL000036_CORR",'CPFa_I_TRIAL000040_CORR',
                                 "CPFa_I_TRIAL000019_CORR",'CPFa_I_TRIAL000020_CORR',
                                 "CPFa_I_TRIAL000021_CORR",'CPFa_I_TRIAL000022_CORR',
                                 "CPFa_I_TRIAL000025_CORR",'CPFa_I_TRIAL000028_CORR',
                                 "CPFa_I_TRIAL000029_CORR",'CPFa_I_TRIAL000031_CORR',
                                 "CPFa_I_TRIAL000034_CORR",'CPFa_I_TRIAL000035_CORR')]

dichoDif(Target,group = df$gender_real, focal.name = 1 , method = c("TID", "MH", "Std", #1 is female
                                                                 "Logistic", "Lord"))
dichoDif(Foil,group = df$gender_real, focal.name = 1 , method = c("TID", "MH", "Std",
                                                                  "Logistic","Lord"))
tpm1 <- difLord(Target, group = df$gender_real, focal.name = 1, model = "2PL") #1 is female?
plot(tpm1, plot = "itemCurve", item = 9)

tpm1 <- difLord(Foil, group = df$gender_real, focal.name = 1, model = "2PL") #1 is female?
plot(tpm1, plot = "itemCurve", item = 10)


#if target faces are actually easier for people of the same race, or female faces easier for same sex.
#Do we see dif for targets and foils differently? i.e. black people can tell they saw a black face comparitavely better than 
#   knowing the face was not there?
df <- 

a <- subset(df, gender_real == 1)
b <- subset(df, gender_real == 2)

df <- merge(df, z, by.x = "test_sessions.bblid", by.y = "bblid")
df <- subset(df, White_Black_Other != 3)

dichoDif(Target,group = df$Race_White, focal.name = 1 , method = c("TID", "MH", "Std",
                                                                    "Logistic", "Lord"))
dichoDif(Foil,group = df$Race_White, focal.name = 1 , method = c("TID", "MH", "Std",
                                                                  "Logistic","Lord"))
tpm1 <- difLord(Target, group = df$Race_White, focal.name = 1, model = "2PL") #1 is female?
plot(tpm1, plot = "itemCurve", item = 16)

tpm1 <- difLord(Foil, group = df$Race_White, focal.name = 1, model = "2PL") #1 is whire
plot(tpm1, plot = "itemCurve", item = 19)


#___#_#__#__________
#_#_#__#_#_#_#_#_#_#_
#_#_#_#_#_#_#__#_#_#_#
#_#_#_#___#_#_#_#__##_

#dif target vs foils(faces that were not there at the beginning)

library(mirt)
library(psych)
#CPF_B_I.9492.csv
x <- read.csv("CPF_A_I.9492.csv")
file <- read.csv("CPF_A_I.9492.csv")
dfx <- read.csv("CPF_repeated_measures_partial.csv")
z <- read.csv("GO_RACE_3.csv")

x <- na.omit(x[,grepl("_CORR",colnames(x))])

xcor <- polychoric(x)$rho

nfactors(xcor,n.obs=4000)
fa.parallel(xcor,n.obs=4000)
VSS.scree(xcor)

factors <- 2

mod <- fa(xcor,factors)$loadings[,1:factors]
loadingss <- fa(xcor,factors)$loadings[,1:factors]##

bifactor_key <- rep(NA,dim(x)[2])

#assigning items to target or foil factor (can disagree)
for (j in 1:dim(x)[2]) {
  if (mod[j,1] == max(mod[j,])) {bifactor_key[j] <- 1}
  if (mod[j,2] == max(mod[j,])) {bifactor_key[j] <- 2}
  try(if (mod[j,3] == max(mod[j,])) {bifactor_key[j] <- 3})
  try(if (mod[j,4] == max(mod[j,])) {bifactor_key[j] <- 4})
  try(if (mod[j,5] == max(mod[j,])) {bifactor_key[j] <- 5})
  try(if (mod[j,6] == max(mod[j,])) {bifactor_key[j] <- 6})
  try(if (mod[j,7] == max(mod[j,])) {bifactor_key[j] <- 7})
  try(if (mod[j,8] == max(mod[j,])) {bifactor_key[j] <- 8})
  try(if (mod[j,9] == max(mod[j,])) {bifactor_key[j] <- 9})
  try(if (mod[j,10] == max(mod[j,])) {bifactor_key[j] <- 10})
}

mod <- bfactor(x,bifactor_key, optimizer = "NR")

params <- coef(mod,simplify=TRUE)$items
rows <- dim(params)[1]
end1 <- dim(params)[2]-2
params <- params[,1:end1]/1.702
slopes <- params[,1:(end1-1)]
intercepts <- params[,end1]
constant <- sqrt(rowSums(slopes^2)+1)
difficulties <- (intercepts*(-1))/constant

params <- data.frame(slopes[,1],difficulties)  # NORMAL METRIC!!!!

#########**----------------------------Dif----------------------------------------**########


loadings <- bifactor_key # Correspond to screener DEP, other DEP, GAD and SOC items ####

results <- data.frame(matrix(NA, nrow=40, ncol=6)) 
colnames(results) <- c("Item", "Main", "Interaction", "DiffAIC", "ChiSq", "ChiSqP") 

#items <- na.omit(x[,grepl("_CORR",colnames(x))])

df <- merge(file, dfx, by.x = "SESSID", by.y = "test_sessions.datasetid")
df <- merge(df, z, by.x = "test_sessions.bblid", by.y = "bblid")#
df <- subset(df, White_Black_Other != 3)
#dfx is repeated measures file cause it includes gender/race
#temp <- na.omit(proband_dff[,grepl("_CORR",colnames(proband_dff))])#from the merged file
df$Race_White <- as.factor(df$Race_White)

items <- na.omit(grep("_CORR", colnames(df), value=TRUE))

i=1
for (item in items) {
  # Determine index of item for testing
  item_ind <- which(items %in% item)
  
  #need to get proband df (includes dif factors)
  tmp_cpf_mod <- bfactor(df[,items[!(items %in% item)]], loadings[-item_ind], technical=list(NCYCLES=2000))
  #write.csv(summary(tmp_cpf_mod)[[1]], file=paste0("/Users/gaarbet/Documents/Filetest", item, "_ex.csv"), row.names=FALSE)
  
  assign(paste0(item, "_cpf_bifactor"), tmp_cpf_mod)
  # General = Face Memory/Working Memory
  scores <- fscores(tmp_cpf_mod, QMC=TRUE) ### Q: defaults okay?
  df[,paste0("General", item, "_ex")] <- scores[,1]
  df[,paste0("Target", item, "_ex")] <- scores[,2] 
  df[,paste0("Foil", item, "_ex")] <- scores[,3] 
  
  #General + Target + Foil + Gender + (General*Gender)
  tmp_log_mod <- glm(df[,item] ~ df[,paste0("General", item, "_ex")] + df[,paste0("Target", item, "_ex")] + df[,paste0("Foil", item, "_ex")] + df$Race_White + df[,paste0("General", item, "_ex")]*df$Race_White, family="binomial")
  assign(paste0(item, "_mod"), tmp_log_mod)
  #General + Target + Foil  
  tmp_simp_log_mod <- glm(df[,item] ~ df[,paste0("General", item, "_ex")] + df[,paste0("Target", item, "_ex")] + df[,paste0("Foil", item, "_ex")], family="binomial")
  assign(paste0(item, "_simp_mod"), tmp_simp_log_mod)
  
  # Run a ChiSq Test 
  chisq <- anova(tmp_simp_log_mod, tmp_log_mod, test="Chisq")
  
  # Put relevant statistics in results dataframe
  results[i, "Item"] <- item
  results[i, "Main"] <- summary(tmp_log_mod)$coefficients[4,3]
  results[i, "Interaction"] <- summary(tmp_log_mod)$coefficients[5,3]
  results[i, "DiffAIC"] <- summary(tmp_log_mod)$aic - summary(tmp_simp_log_mod)$aic
  results[i, "ChiSq"] <- chisq$Deviance[2] 
  results[i, "ChiSqP"] <- chisq[[5]][2] 
  i=i+1
}


