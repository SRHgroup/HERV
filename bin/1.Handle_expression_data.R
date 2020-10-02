# 2018.08.30
# Anne-Mette 
# R version: 3.5.0

# ---------------------------
# Make unified expression data outoput 
# ---------------------------


# ---------------------------
# LOAD Data
# ---------------------------

load('data/raw_data/Expression_data/ERV_expression.RData')


str(erv_expression)
head(erv_expression)
tail(erv_expression)


# ---------------------------
# Handle Data
# ---------------------------

# make expression values nueric 
erv_expression$mean_exp <- as.numeric(erv_expression$mean_exp)

# Generate unique identifier 
erv_expression$sample_id <- paste(erv_expression$patient, erv_expression$cycle, sep = '')

# treatment identifyer 
erv_expression$treatment <- erv_expression$cycle
erv_expression$treatment[erv_expression$cycle == 'C4D28' | erv_expression$cycle == 'C6D28'] <- "After Treatment"
erv_expression$treatment[erv_expression$cycle == 'C1D1' ] <- "Before Treatment"
erv_expression$treatment[erv_expression$cycle == 'healthy' ] <- "Healthy"

# Change end of trestment (EOT) to the correct cyle 
erv_expression$cycle[erv_expression$cycle == 'EOT'] <- 'C6D28'

# study identifyer 
erv_expression$study <- 'SH'
erv_expression$study[erv_expression$cycle == 'healthy'] <- sapply(strsplit(erv_expression$patient[erv_expression$cycle == 'healthy'], '_'), `[[`, 2)

# Patient identifyer 
# Change patient number to correspond to T cell data 
erv_expression$patient <- gsub("PD", "SH", erv_expression$patient)
erv_expression$patient[erv_expression$cycle == 'healthy'] <- gsub("SRR\\d{7}_GSE\\d{5}_", "", erv_expression$patient[erv_expression$cycle == 'healthy'])

# sample identyfier 
erv_expression$sample <- paste(erv_expression$patient, erv_expression$cycle, sep = '') 


# ID 
erv_expression$ID <- NA
erv_expression$ID[erv_expression$study == 'SH'] <- paste(erv_expression$hugo_symbol[erv_expression$study == 'SH'],
                                                         erv_expression$sample[erv_expression$study == 'SH'], 
                                                         sep = '_')

# alter outlier exp levels 
# -----------------------------

# change all expression values below 0.05 to 0.05 to avoid "dilution" of the scale
erv_expression$mean_exp[erv_expression$mean_exp < 0.05 & erv_expression$mean_exp > 0] <- 0.05
# erv_expression[erv_expression$mean_exp > 150,  ] 
erv_expression$mean_exp[erv_expression$mean_exp > 150] <- 150


# plotting order 
erv_expression$treatment <- factor(erv_expression$treatment, levels=c("Healthy","Before Treatment","C1D8", "C1D28", "After Treatment"))


# ---------------------------
# Save new DF
# ---------------------------
erv_exp_dat <- erv_expression


save( erv_exp_dat , file ='data/plot_data/expression_data/expression.RData')












# 2019.03.27
# Annie Borch
# R version: 3.5.0

# ---------------------------
# Make unified expression data outoput 
# ---------------------------


# ---------------------------
# LOAD Data
# ---------------------------

load('data/raw_data/Expression_data/APM_expression.RData')

str(APM_exp)
head(APM_exp)
tail(APM_exp)


# ---------------------------
# Handle Data AMP genes 
# ---------------------------

# make expression values nueric 
APM_exp$mean_exp <- as.numeric(APM_exp$mean_exp)

# Generate unique identifier 
APM_exp$sample_id <- paste(APM_exp$Patient, APM_exp$Cycle, sep = '')

# Change end of trestment (EOT) to the correct cyle 
APM_exp$Cycle[APM_exp$Cycle == 'EOT'] <- 'C6D28'
unique(APM_exp$Cycle)

# treatment identifyer 
table(APM_exp$treatment)
APM_exp$treatment <- APM_exp$Cycle
APM_exp$treatment[APM_exp$Cycle == 'C4D28' | APM_exp$Cycle == 'C6D28'] <- "After Treatment"
APM_exp$treatment[APM_exp$Cycle == 'C1D1' ] <- "Before Treatment"
APM_exp$treatment[APM_exp$Cycle == 'C1D' ] <- "Before Treatment"
APM_exp$treatment[APM_exp$Cycle == 'healthy' ] <- "Healthy"

# # Change end of trestment (EOT) to the correct cyle 
# APM_exp$Cycle[APM_exp$Cycle == 'EOT'] <- 'C6D28'
# unique(APM_exp$Cycle)
# study identifyer 
APM_exp$study <- 'SH'
APM_exp$study[APM_exp$Cycle == 'healthy'] <- sapply(strsplit(APM_exp$Patient[APM_exp$Cycle == 'healthy'], '_'), `[[`, 2)

# Patient identifyer 
# Change patient number to correspond to T cell data 
APM_exp$Patient <- gsub("PD", "SH", APM_exp$Patient)
#APM_exp$Patient[APM_exp$Cycle == 'healthy'] <- gsub("SRR\\d{7}_GSE\\d{5}_", "", APM_exp$Patient[APM_exp$Cycle == 'healthy'])

# sample identyfier 
APM_exp$sample <- paste(APM_exp$Patient, APM_exp$Cycle, sep = '') 


# ID 
APM_exp$ID <- NA
APM_exp$ID[APM_exp$study == 'SH'] <- paste(APM_exp$hugo_symbol[APM_exp$study == 'SH'],
                                           APM_exp$sample[APM_exp$study == 'SH'], 
                                                         sep = '_')

# alter outlier exp levels 
# -----------------------------

# change all expression values below 0.05 to 0.05 to avoid "dilution" of the scale
APM_exp$mean_exp[APM_exp$mean_exp < 0.05 & APM_exp$mean_exp > 0] <- 0.05
# erv_expression[erv_expression$mean_exp > 150,  ] 
APM_exp$mean_exp[APM_exp$mean_exp > 150] <- 150


# plotting order 
APM_exp$treatment <- factor(APM_exp$treatment, levels=c("Healthy","Before Treatment","C1D8", "C1D28", "After Treatment"))



# ---------------------------
# Save new DF
# ---------------------------
APM_expression <- APM_exp
save(APM_expression, file = 'data/plot_data/expression_data/APM_expression_data.RData')




length(unique(APM_expression$treatment))




# ---------------------------
# Handle Data CTA genes 
# ---------------------------

load('data/raw_data/Expression_data/CTA_expression.Rdata')

# make expression values nueric 
CTA_exp$mean_exp <- as.numeric(CTA_exp$mean_exp)

# Generate unique identifier 
CTA_exp$sample_id <- paste(CTA_exp$Patient, CTA_exp$Cycle, sep = '')

# Change end of trestment (EOT) to the correct cyle 
CTA_exp$Cycle[CTA_exp$Cycle == 'EOT'] <- 'C6D28'
unique(CTA_exp$Cycle)

# treatment identifyer 
table(CTA_exp$treatment)
unique(CTA_exp$Cycle)
CTA_exp$treatment <- CTA_exp$Cycle
CTA_exp$treatment[CTA_exp$Cycle == 'C4D28' | CTA_exp$Cycle == 'C6D28'] <- "After Treatment"
CTA_exp$treatment[CTA_exp$Cycle == 'C1D1' ] <- "Before Treatment"
CTA_exp$treatment[CTA_exp$Cycle == 'C1D' ] <- "Before Treatment"
CTA_exp$treatment[CTA_exp$Cycle == 'healthy' ] <- "Healthy"

# # Change end of trestment (EOT) to the correct cyle 
# APM_exp$Cycle[APM_exp$Cycle == 'EOT'] <- 'C6D28'
# unique(APM_exp$Cycle)
# study identifyer 
CTA_exp$study <- 'SH'
CTA_exp$study[CTA_exp$Cycle == 'healthy'] <- sapply(strsplit(CTA_exp$Patient[CTA_exp$Cycle == 'healthy'], '_'), `[[`, 2)

# Patient identifyer 
# Change patient number to correspond to T cell data 
CTA_exp$Patient <- gsub("PD", "SH", CTA_exp$Patient)
#APM_exp$Patient[APM_exp$Cycle == 'healthy'] <- gsub("SRR\\d{7}_GSE\\d{5}_", "", APM_exp$Patient[APM_exp$Cycle == 'healthy'])

# sample identyfier 
CTA_exp$sample <- paste(CTA_exp$Patient, CTA_exp$Cycle, sep = '') 


# ID 
CTA_exp$ID <- NA
CTA_exp$ID[CTA_exp$study == 'SH'] <- paste(CTA_exp$hugo_symbol[CTA_exp$study == 'SH'],
                                           CTA_exp$sample[CTA_exp$study == 'SH'], 
                                           sep = '_')

# alter outlier exp levels 
# -----------------------------

# change all expression values below 0.05 to 0.05 to avoid "dilution" of the scale
CTA_exp$mean_exp[CTA_exp$mean_exp < 0.05 & CTA_exp$mean_exp > 0] <- 0.05
# erv_expression[erv_expression$mean_exp > 150,  ] 
CTA_exp$mean_exp[CTA_exp$mean_exp > 150] <- 150


# plotting order 
CTA_exp$treatment <- factor(CTA_exp$treatment, levels=c("Healthy","Before Treatment","C1D8", "C1D28", "After Treatment"))



# ---------------------------
# Save new DF
# ---------------------------
CTA_expression <- CTA_exp

save(CTA_expression, file = 'data/plot_data/expression_data/CTA_expression_data.RData')





