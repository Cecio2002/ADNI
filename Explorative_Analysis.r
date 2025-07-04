# ---- Progetto Applied: Analisi Esplorativa ----

# 1. Librerie e dati
# -------------------
library(dplyr)      # manipulazione dati
library(tidyr)      # reshape dati
library(ggplot2)    # grafici
library(viridis)    # palette discrete per ggplot2
library(corrplot)
library(factoextra)
library(ggalluvial)
library(igraph)
library(plotly)

# Carica il dataframe (già importato in R)
df <- df_selected_all_features_grouped_correct

# Converte in factor 'visit' e 'DIAGNOSIS'
df <- df %>% 
  mutate(
    visit     = factor(visit, levels = sort(unique(visit))),
    DIAGNOSIS = factor(DIAGNOSIS,
                       levels = c(1,2,3),
                       labels = c("Safe","Mild","Dementia"))
  )
#----variabili numeriche----

# ---- Section 0a: Scatterplot matrix per tutte le variabili numeriche ----
library(GGally)

numeric_vars <- c(
  "Hippocampus_Total",  "Thalamus_Total",       "SuperiorTemporal_Total",
  "Insula_Total",       "Precuneus_Total",      "Mammilare_Total",
  "CaudalMiddleFrontal_Total", "InfLatVentricle_Total"
)

ggpairs(
  df,
  columns = numeric_vars,
  mapping = aes(color = DIAGNOSIS, shape = visit, alpha = 0.5),
  upper = list(continuous = wrap("cor", size = 2)),     # matrice superiore: corr
  lower = list(continuous = wrap("points", size = 0.7)),# matrice inferiore: punti
  diag  = list(continuous = wrap("densityDiag"))        # diagonale: densità
) +
  theme_minimal(base_size = 12) +
  labs(title = "Scatterplot matrix delle variabili numeriche")


df %>%
  select(all_of(numeric_vars)) %>%             # solo le tue variabili numeric
  pivot_longer(cols       = everything(),
               names_to   = "variable",
               values_to  = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(
    title = "Distribuzioni delle variabili numeriche",
    x     = "Valore",
    y     = "Conteggio"
  ) +
  theme_minimal(base_size = 14)


GGally::ggpairs(
  df[, c(numeric_vars, "DIAGNOSIS", "visit")],
  mapping = aes(color = DIAGNOSIS, shape = visit, alpha = 0.5),
  columns = 1:length(numeric_vars),
  upper   = list(continuous = wrap("density", alpha = 0.3)),
  lower   = list(continuous = wrap("points", size = 0.5))
) +
  theme_minimal(base_size = 10) +
  labs(title="Scatterplot‐matrix: colore=diagnosi, shape=visita")


# ---- Section 1: Traiettorie individuali (Pazienti “changer”) ----

# 1.1 Filtra i pazienti che cambiano diagnosi almeno una volta
df_changer <- df %>%
  group_by(Subject_ID) %>%
  filter(n_distinct(DIAGNOSIS) > 1) %>%
  ungroup()

# 1.2 Definisci le variabili numeriche (colonne 3–10)
numeric_vars <- names(df_changer)[3:10]

# 1.3 Passa in formato long per il plotting
df_long <- df_changer %>%
  select(Subject_ID, visit, all_of(numeric_vars)) %>%
  pivot_longer(
    cols      = all_of(numeric_vars),
    names_to  = "variable",
    values_to = "value"
  )

# 1.4 Plot delle traiettorie individuali
ggplot(df_long, aes(x = visit, y = value, group = Subject_ID, color = Subject_ID)) +
  geom_line(size = 0.8, alpha = 0.8) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title    = "Traiettorie individuali delle variabili numeriche",
    subtitle = "Pazienti con diagnosi variabile – ogni linea un paziente",
    x        = "Visita",
    y        = "Valore"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 12),
    plot.title      = element_text(face = "bold", size = 16)
  )



# ---- Section 2: Conteggio e Distribuzione Diagnosi per Visita ----

# 2.1 Bar-plot facet per diagnosi
df %>%
  group_by(visit, DIAGNOSIS) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = visit, y = n)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ DIAGNOSIS, scales = "free_y") +
  labs(
    x     = "Visita",
    y     = "Numero pazienti",
    title = "Conteggio pazienti per diagnosi e visita"
  ) +
  theme_minimal(base_size = 14)

ggplot(df, aes(x = visit)) +
  geom_bar(
    aes(fill = DIAGNOSIS, group = DIAGNOSIS),
    position = position_dodge(width = 0.7),
    width    = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "Safe"     = "#a1d99b",
      "Mild"     = "#41ab5d",
      "Dementia" = "#6a51a3"  
    )
  ) +
  labs(
    x     = "Visita",
    y     = "Numero pazienti",
    fill  = "Diagnosi",
    title = "Distribuzione delle diagnosi per visita"
  ) +
  theme_minimal(base_size = 14)


# ---- 3. Distribuzioni (Density) per variabili numeriche ----
for (v in numeric_vars) {
  print(
    ggplot(df, aes_string(x=v, color="DIAGNOSIS", fill="DIAGNOSIS")) +
      geom_density(alpha=0.3) +
      facet_wrap(~ visit, scales="free") +
      labs(title=paste("Density di", v, "per Diagnosi e Visita"),
           x=v, y="Densità") +
      theme_minimal()
  )
}

# ---- 4. Heatmap di correlazione tra biomarcatori ----
M <- cor(df[numeric_vars], use="pairwise.complete.obs")
corrplot(M, method="color", type="upper", tl.cex=0.7,
         title="Correlazioni tra variabili numeriche")


# ---- 5. PCA Biplot colorato per diagnosi ----
res.pca <- prcomp(df[numeric_vars], scale.=TRUE, na.action=na.omit)
fviz_pca_biplot(
  res.pca,
  geom.ind="point",
  col.ind=df$DIAGNOSIS,
  addEllipses=TRUE,
  palette="Set2",
  title="PCA Biplot — diagnosi"
)

# ---- 6. Alluvial plot delle transizioni diagnostiche ----
library(ggplot2)
library(dplyr)
library(ggalluvial)

df_alluv <- df %>% mutate(visit = paste0("V", visit))

ggplot(df_alluv,
       aes(x = visit, stratum = DIAGNOSIS, alluvium = Subject_ID,
           y = 1, fill = DIAGNOSIS)) +
  geom_flow(alpha = 0.6) +
  geom_stratum(alpha = 0.9) +
  scale_fill_manual(values = c(
    "Safe"     = "#A1D99B",  # Verde chiaro
    "Mild"     = "#41AB5D",  # Verde medio
    "Dementia" = "#006D2C"   # Verde scuro
  )) +
  labs(title = "Flussi diagnostici V1→V4",
       x = "Visita", y = "Numero pazienti") +
  theme_minimal(base_size = 14)

# ---- 7. Delta (V4−V1) e boxplot per diagnosi iniziale ----

# 1) Definiamo le tue variabili numeriche
numeric_vars <- c(
  "Hippocampus_Total",  "Thalamus_Total",       "SuperiorTemporal_Total",
  "Insula_Total",       "Precuneus_Total",      "Mammilare_Total",
  "CaudalMiddleFrontal_Total", "InfLatVentricle_Total"
)

# 2) Calcoliamo i delta restando in long e poi tornando wide solo su V1 e V4
df_delta <- df %>%
  # Se visit è factor, lo trasformiamo in numerico
  mutate(visit = as.numeric(as.character(visit))) %>%
  select(Subject_ID, visit, all_of(numeric_vars), DIAGNOSIS) %>%
  pivot_longer(
    cols      = all_of(numeric_vars),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  # FILTRIAMO sulle visite 1 e 4
  filter(visit %in% c(1, 4)) %>%
  # Wide solo su visit per ciascuna variabile
  pivot_wider(
    id_cols      = c(Subject_ID, DIAGNOSIS, variable),
    names_from   = visit,
    values_from  = value,
    names_prefix = "V"
  ) %>%
  # Calcoliamo il delta e riordiniamo DIAGNOSIS
  mutate(
    delta     = V4 - V1,
    DIAGNOSIS = factor(DIAGNOSIS, levels = c("Safe", "Mild", "Dementia"))
  )

# 3) Creiamo il boxplot
library(ggplot2)
library(RColorBrewer)

ggplot(df_delta, aes(x = DIAGNOSIS, y = delta, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.6, outlier.size = 0.8) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Delta (V4 − V1) per diagnosi iniziale",
    x     = "Diagnosi a V1",
    y     = "V4 − V1"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# ---- 8. Clustering “soft” con k-means sui profili di slope ----
# Calcola slope per ogni paziente e variabile
df_slopes <- df %>%
  mutate(nv = as.integer(visit)) %>%
  pivot_longer(all_of(numeric_vars), names_to="variable", values_to="value") %>%
  group_by(Subject_ID, variable) %>%
  summarise(slope = if(sum(!is.na(value))>=2) coef(lm(value~nv))[2] else NA_real_,
            .groups="drop") %>%
  pivot_wider(names_from="variable", values_from="slope") %>%
  drop_na()

# k-means su slopes
set.seed(42)
km <- kmeans(scale(df_slopes[-1]), centers=3)
cent <- as.data.frame(km$centers) %>% mutate(cluster=1:3) %>%
  pivot_longer(-cluster, names_to="variable", values_to="value")

ggplot(cent, aes(variable, value, group=cluster, color=factor(cluster))) +
  geom_line(size=1) +
  labs(title="Centroidi k-means sui profili di slope") +
  theme_minimal()


# ---- 10. Network di correlazioni tra biomarcatori ----
M2 <- cor(df[numeric_vars], use="pairwise.complete.obs")
edges <- which(abs(M2) > 0.7 & abs(M2) < 1, arr.ind=TRUE)
g <- graph_from_data_frame(
  data.frame(
    from = rownames(M2)[edges[,1]],
    to   = colnames(M2)[edges[,2]],
    weight = M2[edges]
  ), directed=FALSE
)
plot(g, vertex.label.cex=0.8,
     main="Network correlazioni |r|>0.7")


# ---- 11. Dashboard interattivo con plotly ----
p <- ggplot(df, aes(x=visit, y=Hippocampus_Total, color=DIAGNOSIS)) +
  geom_jitter(width=0.2, size=1) +
  labs(title="Hippocampus_Total per visita e diagnosi") +
  theme_minimal()
ggplotly(p)


#----variabili categoriche---------------------------------------------------

library(dplyr)
library(ggplot2)
library(scales)
dfB<-Dataset_Alz_feat_eng
# Definisci quali visite plottare e le loro etichette
visits_to_plot <- c(1, 2, 3,4)
visit_labels  <- c("visit1", "visit2", "visit3","visit4")

# Assicurati che cat_vars contenga i nomi delle colonne 3:11
cat_vars <- names(dfB)[3:11]
# Codici esadecimali consigliati per azzurro, blu, verde acqua
colori <- c("#8ecae6", "#219ebc", "#023047")

for (v in cat_vars) {
  df_tmp <- dfB %>%
    filter(Visit %in% visits_to_plot) %>%
    mutate(
      Visit = factor(Visit,
                     levels = visits_to_plot,
                     labels = visit_labels),
      Diagnosi = factor(DIAGNOSIS)
    ) %>%
    group_by(Visit, !!sym(v), Diagnosi) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Visit, !!sym(v)) %>%
    mutate(prop = n / sum(n))
  
  p <- ggplot(df_tmp,
              aes_string(x = v, y = "prop", fill = "Diagnosi")) +
    geom_col(position = "stack") +
    facet_wrap(~ Visit, nrow = 1) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = colori) +   # aggiungi questa riga
    labs(
      title = paste("Distribuzione di", v, "nei tre time-point"),
      x     = v,
      y     = "Proporzione",
      fill  = "Diagnosi"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
  
  print(p)
}



#----REDUCTION: PCA E MCA ----
# install.packages(c("readxl","dplyr","FactoMineR","factoextra","rgl"))
library(readxl)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(rgl)

# 1) Carica il dataset
df2 <- read_excel("Dataset_Alz_feat_eng.xlsx")

# ───────────────────────────────────────────────
# 2) PCA sulle variabili numeriche (cols 13–20 già standardizzate)
# ───────────────────────────────────────────────
num_vars <- df2[, 13:20]

pca_res <- prcomp(num_vars, center = FALSE, scale. = FALSE)

# Soglia 70% varianza cumulata
var_perc  <- pca_res$sdev^2 / sum(pca_res$sdev^2)
cumvar    <- cumsum(var_perc)
threshold <- 0.70
k_opt     <- which(cumvar >= threshold)[1]
message("PCA → scelti ", k_opt,
        " PC (", round(cumvar[k_opt]*100,1), "% varianza)\n")

pc_scores <- as.data.frame(
  pca_res$x[, seq_len(k_opt), drop = FALSE]
)
colnames(pc_scores) <- paste0("PC", seq_len(k_opt))

# 3) PARTE “SPIEGAZIONE” DELLA PCA
## 3a) Riassunto della varianza spiegata
print(summary(pca_res))

## 3b) Loadings delle variabili sui primi k_opt PC
loadings <- pca_res$rotation[, seq_len(k_opt), drop = FALSE]
print(loadings)

## 3c) Contributi delle variabili ai primi 3 assi
fviz_contrib(pca_res, choice = "var", axes = 1, top = ncol(num_vars)) +
  ggtitle("Contributi variabili a PC1")
fviz_contrib(pca_res, choice = "var", axes = 2, top = ncol(num_vars)) +
  ggtitle("Contributi variabili a PC2")
fviz_contrib(pca_res, choice = "var", axes = 3, top = ncol(num_vars)) +
  ggtitle("Contributi variabili a PC3")

## 3d) Biplot 3D (PC1–PC3)
ind3d   <- pca_res$x[, 1:3]
var3d   <- pca_res$rotation[, 1:3]
scale_f <- max(abs(ind3d)) / max(abs(var3d)) * 0.8
var3d_s <- var3d * scale_f

open3d()
plot3d(ind3d,
       col   = "blue", size = 5,
       xlab  = "PC1", ylab = "PC2", zlab = "PC3",
       main  = "Biplot 3D PCA (PC1-PC3)")
for (i in seq_len(nrow(var3d_s))) {
  v    <- var3d_s[i, ]
  name <- rownames(var3d_s)[i]
  segments3d(rbind(c(0,0,0), v), col="red", lwd=2)
  text3d(v[1],v[2],v[3], texts=name, adj=c(1.2,1.2), col="red")
}

# ───────────────────────────────────────────────
# 4) MCA sulle 9 variabili categoriche
# ───────────────────────────────────────────────
cat_vars <- df2 %>%
  select(
    gender,
    genotype,
    rs744373_C, rs11767557_C, rs11771145_A,
    rs11136000_T, rs3851179_A, rs17125944_C,
    rs3764650_G
  ) %>%
  mutate(across(everything(), as.factor))

mca_res <- MCA(cat_vars, graph = FALSE)

# 5) PARTE “SPIEGAZIONE” DELLA MCA
## 5a) Riassunto dell’inerzia
print(mca_res$eig)
# Col1 = eigenvalue; Col2 = % of variance; Col3 = cumulative %

## 5b) Loadings (coordinate modalità) su Dim1 e Dim2
var_coord <- mca_res$var$coord[, 1:2, drop = FALSE]
print(var_coord)

## 5c) Contributi (%) modalità a Dim1 e Dim2
var_contrib <- mca_res$var$contrib[, 1:2, drop = FALSE]
colnames(var_contrib) <- c("Dim1","Dim2")
print(var_contrib)

# Somma contributi per variabile
cat_contrib <- data.frame(
  variable = sub("=.*", "", rownames(var_contrib)),
  var_contrib,
  check.names = FALSE
)
var_contrib_by_var <- cat_contrib %>%
  group_by(variable) %>%
  summarise(
    contrib_Dim1 = sum(Dim1),
    contrib_Dim2 = sum(Dim2)
  )
print(var_contrib_by_var)



# ───────────────────────────────────────────────
# 6) Unisco PCA + MCA + chiavi in un unico dataframe
# ───────────────────────────────────────────────
# (uso mca_res$ind$coord e m_opt per estrarre tutte MC1…MCm_opt)
eig_vals    <- mca_res$eig[, "percentage of variance"] / 100
cum_inerc   <- cumsum(eig_vals)
# Numero di dimensioni realmente calcolate
available_dims <- ncol(mca_res$ind$coord)

# Assicuriamoci che m_opt non superi available_dims
m_opt <- min(m_opt, available_dims)
if (is.na(m_opt) || m_opt < 1) m_opt <- ncol(mca_res$ind$coord)

mca_scores <- as.data.frame(
  mca_res$ind$coord[, seq_len(m_opt), drop = FALSE]
)
colnames(mca_scores) <- paste0("MC", seq_len(m_opt))

reduced_df <- bind_cols(
  df2 %>% select(Subject_ID, DIAGNOSIS, Visit),
  pc_scores,
  mca_scores
)

# 7) Verifica
glimpse(reduced_df)
head(reduced_df)


#----LMM----
# 1. Caricamento librerie
library(readxl)    # per read_excel
library(dplyr)     # per manipolazioni dati
library(lme4)      # per lmer
library(stats)     # per prcomp

# 2. Importa il dataset
df <- read_excel("Dataset_Alz_feat_eng.xlsx")

# 3. Pre‐processing dei dati
#    - Assicuriamoci che Subject_ID, gender, DIAGNOSIS e visit siano fattori
df <- df %>%
  mutate(
    Subject_ID = factor(Subject_ID),
    gender     = factor(gender),
    DIAGNOSIS  = factor(DIAGNOSIS),
    visit      = as.numeric(Visit)   # se vuoi usarlo come covariata numerica
  )

# 4. Codifica SNP/genotipo e PCA
#    Supponendo che 'genotype' sia una variabile categorica, creiamo dummies:
geno_dummies <- model.matrix(~ genotype - 1, data = df)
#    Eseguiamo la PCA sui dummy SNP
pca_geno <- prcomp(geno_dummies, center = TRUE, scale. = TRUE)
#    Prendiamo i primi 3 componenti principali
df$PC1 <- pca_geno$x[,1]
df$PC2 <- pca_geno$x[,2]
df$PC3 <- pca_geno$x[,3]

# 5. Specifichiamo il modello misto
#    Sostituisci 'hippocampus_volume' con la tua variabile MRI di interesse
modello <- lmer(
  Hippocampus_Total ~ visit + DIAGNOSIS + gender + PC1 + PC2 + PC3 +
    (1 + visit | Subject_ID),
  data = df,
  REML = FALSE
)

# 6. Riepilogo risultati
summary(modello)


#----MANOVA----
df<-Dataset_Alz_feat_eng
library(dplyr)

# Get only the last visit for each subject
df_last <- df %>%
  group_by(Subject_ID) %>%
  filter(Visit == max(Visit)) %>%
  ungroup()

# Check class balance
table(df_last$DIAGNOSIS)


# Run MANOVA with brain volumes as dependent variables
manova_model <- manova(cbind(Hippocampus_Total, Thalamus_Total, SuperiorTemporal_Total,
                             Insula_Total, Precuneus_Total, Mammilare_Total,
                             CaudalMiddleFrontal_Total, InfLatVentricle_Total) 
                       ~ DIAGNOSIS, data = df_last)

# Global MANOVA test
summary(manova_model, test = "Wilks")

# Follow-up univariate ANOVAs
summary.aov(manova_model)
