# --- LIBRERIE NECESSARIE ---
library(mvtnorm)
library(rgl)
library(car)
library(MVN)
library(biotools)

# --- LETTURA E PULIZIA DATI ---
alz <- read.csv("alz_finale_unico.csv")
alz$DIAGNOSIS <- factor(alz$DIAGNOSIS)
alz$gender <- factor(alz$gender)

# Selezione variabili neuroanatomiche
brain_vars <- c("Hippocampus_Total", "Thalamus_Total", "SuperiorTemporal_Total", 
                "Insula_Total", "Precuneus_Total", "Mammilare_Total", 
                "CaudalMiddleFrontal_Total", "InfLatVentricle_Total")

alz_sub <- alz[, brain_vars]
alz_sub <- na.omit(alz_sub)  # Rimuove righe con NA

# --- DISTANZA E CLUSTERING ---
dist.matrix <- dist(alz_sub, method = "euclidean")
image(1:nrow(alz_sub), 1:nrow(alz_sub), as.matrix(dist.matrix), main = "Distance Matrix", asp = 1)

# clustering
cluster.dl <- hclust(dist.matrix, method = "ward.D2")

# dendrogramma
plot(cluster.dl, main = "Dendrogram", hang = -0.1, xlab = "", labels = FALSE, cex = 0.6)
rect.hclust(cluster.dl, k = 2)  # numero cluster modificabile

# assegna i cluster
cluster.dl.cut <- cutree(cluster.dl, k = 2)
alz_sub$cluster <- factor(cluster.dl.cut)

# --- COPHENETIC CORRELATION ---
coph.coeff <- cor(dist.matrix, cophenetic(cluster.dl))
cat("Cophenetic correlation coefficient:", coph.coeff, "\n")

# --- PLOT CLUSTER COLORATI ---
pairs(alz_sub[, brain_vars], col = as.numeric(alz_sub$cluster), pch = 19)

# --- MANOVA ---
fit.manova <- manova(as.matrix(alz_sub[, brain_vars]) ~ cluster, data = alz_sub)
summary(fit.manova)

# --- TEST ASSUNZIONI ---
for (k in levels(alz_sub$cluster)) {
  cat("Cluster", k, "\n")
  print(mvn(alz_sub[alz_sub$cluster == k, brain_vars], mvnTest = "hz"))
}

boxM(alz_sub[, brain_vars], alz_sub$cluster)



library(ggplot2)
alz_sub$DIAGNOSIS_GROUPED <- ifelse(alz_sub$DIAGNOSIS %in% c(1, 2), "Non-Alzheimer", "Alzheimer")
alz_sub$DIAGNOSIS_GROUPED <- factor(alz_sub$DIAGNOSIS_GROUPED)
for (var in brain_vars) {
  p1 <- ggplot(alz_sub, aes_string(x = "factor(cluster)", y = var, fill = "factor(cluster)")) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste("Boxplot of", var, "by Cluster"), x = "Cluster", y = var) +
    theme_minimal()
  
  p2 <- ggplot(alz_sub, aes_string(x = "DIAGNOSIS_GROUPED", y = var, fill = "DIAGNOSIS_GROUPED")) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste("Boxplot of", var, "by Diagnosis (Grouped)"), x = "Diagnosis", y = var) +
    theme_minimal()
  
  print(p1)
  print(p2)
}
