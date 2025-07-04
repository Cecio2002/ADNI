# 28/05/2025
#------- MULTINOMIAL MIXED MODELS ----------
library(readxl)
library(brms)

alz <- read_excel("Dataset_Alz_feat_eng.xlsx")
# Convert diagnosis to factor and label levels
alz$DIAGNOSIS <- factor(alz$DIAGNOSIS, levels=c(1,2,3),
                        labels=c("CN","MCI","AD"))
alz$Subject_ID <- factor(alz$Subject_ID)   # ensure subject is a factor
alz$gender <- factor(alz$gender, levels=1:2, labels=c("M","F"))
alz$rs744373_C <-factor(alz$rs744373_C, levels=0:2)
alz$rs11767557_C <-factor(alz$rs11767557_C, levels=0:2)
alz$rs11136000_T <-factor(alz$rs11136000_T, levels=0:2)
alz$rs11771145_A <-factor(alz$rs11771145_A, levels=0:2)
alz$rs3764650_G <-factor(alz$rs3764650_G , levels=0:2)
alz$rs17125944_C <- factor(alz$rs17125944_C , levels=0:2)
alz$rs3851179_A <- factor(alz$rs3851179_A , levels=0:2)
alz$genotype <- factor(alz$genotype, levels= c("2/2", "3/2", "3/3", "4/2","4/3", "4/4"))
alz$Visit.f <- factor(alz$Visit, levels =1:4)


## modello base 1 ( random intercept for Subject_ID) ---- 
alz_fit <- brm(DIAGNOSIS ~ Hippocampus_Total + gender + (1 | Subject_ID),
               data=alz, family=categorical(),
               iter=2000, chains=2, silent=TRUE)
summary(alz_fit)  # population and group-level estimates
report(alz_fit)
# Plot conditional effects for each diagnosis level
plot(conditional_effects(alz_fit, categorical = TRUE), ask = FALSE)


## modello 2 (piu predittori) ----- 
alz_fit.2 <- brm(DIAGNOSIS ~ Hippocampus_Total + Thalamus_Total + Visit + (1 | Subject_ID),
               data=alz, family=categorical(),
               iter=2000, chains=2, silent=TRUE)
summary(alz_fit.2)  # population and group-level estimates

# Plot conditional effects for each diagnosis level
plot(conditional_effects(alz_fit.2, categorical = TRUE), ask = FALSE)


## modello 3 (aggiungamo anche un random slope per la visit) ----- 
alz_fit.3 <- brm(DIAGNOSIS ~ Hippocampus_Total  + Thalamus_Total + Visit + (1 +Visit | Subject_ID),
                 data=alz, family=categorical(),
                 iter=2000, chains=2, silent=TRUE)
summary(alz_fit.3)
plot(conditional_effects(alz_fit.3, categorical = TRUE), ask = FALSE)


## modello 4 (mettiamo dentro qualche variabile categorica preselezionata con analisi esplorativa) ----- 
alz_fit.4 <- brm(DIAGNOSIS ~ Hippocampus_Total  + Thalamus_Total + Visit + rs17125944_C + (1| Subject_ID),
                 data=alz, family=categorical(),
                 iter=2000, chains=2, silent=TRUE)
summary(alz_fit.4)
plot(conditional_effects(alz_fit.4, categorical = TRUE), ask = FALSE)



## modello 5 (lascio dentro solo due numeriche e provo con piu SNP poiche sembrano molto promettenti) ----- 
alz_fit.5 <- brm(DIAGNOSIS ~ Hippocampus_Total  + Thalamus_Total  + rs17125944_C + rs3851179_A+ rs744373_C+ (1| Subject_ID),
                 data=alz, family=categorical(),
                 iter=2000, chains=2, silent=TRUE, control = list(adapt_delta = 0.99))
summary(alz_fit.5)
plot(conditional_effects(alz_fit.5, categorical = TRUE), ask = FALSE)


## modello 6 (levo rs17125944_C perche classe 2 troppo poco numerosa, per adesso provo solo genotype e poi aggiungerò dopo) ----- 
alz_fit.6 <- brm(DIAGNOSIS ~ Hippocampus_Total  + Thalamus_Total + genotype + (1 | Subject_ID),
                 data=alz, family=categorical(),
                 iter = 4000, chains=2, silent=TRUE, control = list(adapt_delta = 0.95,max_treedepth = 15))
summary(alz_fit.6)
plot(conditional_effects(alz_fit.6, categorical = TRUE), ask = FALSE)

