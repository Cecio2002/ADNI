# Leggi il file TSV
gwas_data <- read.delim("~/Downloads/gwas-association-downloaded_2025-05-07-MONDO_0004975-withChildTraits.tsv", header = TRUE, sep = "\t")

# Visualizza le prime righe
head(gwas_data)

# Ordina per p-value crescente
#gwas_sorted <- gwas_data[order(gwas_data$P.VALUE), ]

# Seleziona i primi 50 SNP
top_snps <- head(gwas_data$SNPS, 500)



# Carica il file .bim di PLINK
bim_data <- read.table("~/Desktop/Project/Plinkk/ADNI_GO_2_Forward_Bin.bim", header = FALSE)

# Visualizza i primi 10 SNP
head(bim_data)

# Estrai la lista degli SNP presenti nel .bim
bim_snps <- bim_data$V2

# Rimuovi APOE e filtra solo quelli presenti nel .bim
valid_snps <- top_snps[top_snps %in% bim_snps & grepl("^rs\\d+$", top_snps)]

# Quanti ne restano?
length(valid_snps)


# Scrivi in file di testo, uno per riga
writeLines(valid_snps, "Alzheimer_snps_valid.txt")


#uso il file fatto da chat con i 50 SNP piu importanti ( quello fatto prima è una paraculata)
system("~/Desktop/Project/Plinkk/plink --bfile ~/Desktop/Project/Plinkk/ADNI_GO2_GWAS_2nd_orig_BIN  --extract ~/Desktop/Project/Plinkk/snp_alzheimer_list.txt --recode A --out ~/Desktop/Project/Plinkk/patient_snps")


# Carica il file raw in R
raw_data <- read.delim("~/Desktop/Project/Plinkk/patient_snps.raw", header = TRUE)

# Visualizza le prime righe dei dati
head(raw_data)

# Separare i nomi delle colonne dalla prima riga
new_colnames <- strsplit(colnames(raw_data), "\\.")[[1]]

# Visualizza i nomi separati
new_colnames

raw_data_clean <- read.delim("~/Desktop/Project/Plinkk/patient_snps.raw", header = FALSE, skip = 1)

# Separare i dati nelle colonne usando lo spazio come delimitatore
raw_data_clean <- as.data.frame(t(sapply(strsplit(as.character(raw_data_clean$V1), " "), identity)))

# Assegna i nomi delle colonne
colnames(raw_data_clean) <- new_colnames

# Verifica il nuovo dataframe
head(raw_data_clean)

write.csv(raw_data_clean, "~/Desktop/raw_data_clean.csv", row.names = FALSE)

