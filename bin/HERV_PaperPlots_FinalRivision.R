# Anne-Mette 
# 2017.11.20
# R version: 3.6.2
# edited: 2020.03.23

# -------------------------
# Plotting for the ERV manuscript 
# All plot generated for the publication 
# --------------------------


# Load Libraries 
# ----------------------
# install.packages('openxlsx')
# install.packages('ggrepel')
# install.packages('ggbeeswarm')
# install.packages('gridExtra')
# install.packages('cowplot')
# install.packages("ggsignif")
# install.packages('ggplot2')

library(openxlsx)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(ggsignif)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)

# Load data 
# ----------------------

load('data/plot_data/2019.08.25.HERV_barcode.RData')
load( 'data/plot_data/expression.RData')
clin_HHRH <- read.table('data/plot_data/2018.11.08.clinical_HHRH.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
clin_SH <- read.table('data/plot_data/clinical_SH.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# load in HERv expression gene list 
HERV_library_postive <-read.xlsx('data/HERV expression and T cell_Sunil 29042019.xlsx', sheet = 7)
HERV_library_postive$hugo_symbol

# load AMP expression
load('data/plot_data/APM_expression_data.RData')

# load CTA expression 
load('data/plot_data/CTA_expression_data.RData')



# look at data 
# ----------------------
# baracoda data 
head(barcode_dat)
dim(barcode_dat)
tail(barcode_dat)
table(barcode_dat$sample)
table(barcode_dat$patient)
table(barcode_dat$patient, barcode_dat$cycle)
table(barcode_dat$HLAallele)

# RNAseq expression data 
head(erv_exp_dat)
tail(erv_exp_dat)
table(erv_exp_dat$cycle)

# clinical 
head(clin_HHRH)
head(clin_SH)



# ------------------------------------------------------------------------

# ############### #
#                 #
#   HANDLE DATA   #
#                 #
# ############### #

# ------------------------------------------------------------------------

# --------------------------
# Handle expression data for Heatmap
# --------------------------

# Define colour breaks
# my_breaks = c(0,0.05,0.5,5,50,500)

# remove spaces from names 
head(erv_exp_dat)
erv_exp_dat$hugo_symbol <- gsub(' ','', erv_exp_dat$hugo_symbol)
erv_exp_dat$accession_number <- gsub(' ','', erv_exp_dat$accession_number)

# make plotting data frame without unnessesary samples and genes (viral mimecry genes)
plot_exp <- erv_exp_dat[erv_exp_dat$id != 'ViralMimicry' & # no Viral mimery 
                            erv_exp_dat$treatment %in% c('Healthy', 'Before Treatment', 'After Treatment') & # only treatment for healthy befor and after 
                            erv_exp_dat$hugo_symbol %in% HERV_library_postive$hugo_symbol & # only hugo symbols where T cell responce was relevant 
                            !grepl('GMP|MEP|MGP|GSE69905_1', erv_exp_dat$patient), ] # exclude unvanted healthy donor samples 

head(plot_exp)
# select highes ranking transciprt 
# ...................................
# make temporary df of transcript level expression taking the sum of transcripts across all patients (not healthy donors)
plot_exp_tmp <- aggregate(mean_exp ~ accession_number + hugo_symbol, 
                          plot_exp[plot_exp$treatment != 'Healthy', ],
                          sum)
plot_exp_tmp$hugo_symbol <- as.character(plot_exp_tmp$hugo_symbol)
head(plot_exp_tmp)
# find rank of transcripts to be plottet (row with max sum of expression)
transcript_ranks <- unsplit(lapply(split(plot_exp_tmp, plot_exp_tmp$hugo_symbol),
                                   function(x) {
                                       x$rank <- rank(-x$mean_exp, ties.method = "first")
                                       x
                                   }), plot_exp_tmp$hugo_symbol)

# Final expression df only incluring highest ranking transcript 
final_exp_plot <- plot_exp[plot_exp$accession_number %in% transcript_ranks$accession_number[transcript_ranks$rank == 1], ]


# aggregate according to sum of expression for all patients pr erv 
exp_pr_erv <- aggregate(mean_exp ~ hugo_symbol + treatment, 
                                 data = final_exp_plot,
                                 FUN = sum)

# mean expression of HERVs all all patients 
exp_pr_erv_mean <- aggregate(mean_exp ~ hugo_symbol, 
                        data = final_exp_plot[final_exp_plot$treatment != 'Healthy', ],
                        FUN = mean)

# aggregate according to MEAN of expression for all HERVS pr patient 
mean_exp_pr_patient <- final_exp_plot %>% 
                    filter(treatment %in% c('Before Treatment', 'After Treatment')) %>% 
                    aggregate(mean_exp ~ patient + treatment, data = ., FUN = mean)

# aggregate according to ... 
exp_pr_patient <- aggregate(mean_exp ~ patient + hugo_symbol + treatment, data = final_exp_plot, FUN = mean)



# Antigen presenting machienery (AMP) expression
# ...................................
unique(APM_expression$Patient)
APM_expression <- APM_expression[!APM_expression$Patient=="GSE69905_1",]

# Cancer testis antigen (CTA) expression
# ...................................
CTA_expression <- CTA_expression[!CTA_expression$Patient=="GSE69905_1",]


# --------------------------
# COUNT df 
# --------------------------
head(barcode_dat)

# take out samle BC101 and BC96 
barcode_dat <- barcode_dat[!barcode_dat$patient %in%  c("BC101","BC96"), ]

# generate unique long ID
barcode_dat$LongID <- paste(barcode_dat$sample, barcode_dat$Library.type, sep = '_')

# Remove space from petides and gene names 
barcode_dat$Peptide <- gsub(' ','', barcode_dat$Peptide)
barcode_dat$Protein.name <- gsub(' ','', barcode_dat$Protein.name)


# LOG2FC
# calculation for big data frame, MEAN of Log2FC  
barcode_dat$mean_logfold_signif <- 0
barcode_dat$mean_logfold_signif[barcode_dat$signif == TRUE] <- ave(barcode_dat$log_fold_change[barcode_dat$signif == TRUE], 
                                                                           barcode_dat$LongID[barcode_dat$signif == TRUE], 
                                                                           FUN = mean, na.rm = TRUE)
# calculation for big data frame, SUM of Log2FC 
barcode_dat$sum_logfold_signif <- 0
barcode_dat$sum_logfold_signif[barcode_dat$signif == TRUE] <- ave(barcode_dat$log_fold_change[barcode_dat$signif == TRUE], 
                                                                   barcode_dat$LongID[barcode_dat$signif == TRUE], 
                                                                   FUN = sum, na.rm = TRUE)

# get significant responces pr patient for each library type 
lib_signif_patient <- aggregate(signif ~ patient +  Library.type,sum ,data = barcode_dat)

# calculate number of significant responces pr patient for each HERV 
lib_signif_erv <- barcode_dat %>% 
                    filter(Library.type == 'HERV' & study != 'BC') %>% # select only HERV library and samples from patients 
                    separate_rows(Protein.name, sep = ",") %>% # seperating the hervnamees statid for the same entry (separate_rows)
                    aggregate(data = ., signif ~ Protein.name, sum) # calculating the sum of significant entries 

# make HERV response column
lib_signif_erv$HERVRespons <- ifelse(lib_signif_erv$signif, 'Yes', 'No')
table(lib_signif_erv$HERVRespons)

# generate named vector 
lib_signif_erv_vec <- setNames(lib_signif_erv$HERVRespons, lib_signif_erv$Protein.name )

# Add to data expression data  
exp_pr_erv_mean$HERVRespons <- lib_signif_erv_vec[exp_pr_erv_mean$hugo_symbol]



# FOR SUNIL 
# ---------------------------------
HERV_library_size <- barcode_dat %>% 
                        filter(Library.type == 'HERV') %>% # Only look at HERVS
                        separate_rows(Protein.name, sep = ",") %>% # separate individual HERVs symbols where one peptide is found in several HERVs
                        aggregate(data = ., Peptide ~ Protein.name + HLAallele, FUN = n_distinct) %>% #calculate number of unique peptides for each protein and HLA 
                        pivot_wider(names_from = HLAallele, values_from = Peptide) # make wide dataframe 

# calculate sum of total peptide library
HERV_library_size <- HERV_library_size %>% mutate(sum = rowSums(.[2:5], na.rm = TRUE))

# write to excel file 
write.xlsx(HERV_library_size, file=paste(dir,'data/HERV_peptide_library_size.xlsx', sep = ''))



# pre and post aza
# .........................
# 
herv_exp_pr_patient <- barcode_dat %>% 
                        filter(Library.type == 'HERV' & treatment != 'Healthy') %>% 
                        aggregate(data = ., Peptide ~ patient + HLAallele + treatment, FUN = n_distinct) %>% #calculate number of unique peptides for each protein and HLA 
                        pivot_wider(names_from = HLAallele, values_from = Peptide) # make wide dataframe 

# calculate sum of total peptide library
herv_exp_pr_patient <-  mutate_each(herv_exp_pr_patient,funs(as.numeric), 3:6) %>% 
                        mutate(sum = rowSums(.[3:6], na.rm = TRUE)) %>% 
                        arrange(treatment, patient)

# calculate positive peptides for each HERV
patient_response_peptides <- barcode_dat %>% 
                            filter(Library.type == 'HERV' & treatment != 'Healthy') %>%
                            aggregate(data = ., Peptide ~ patient + treatment + signif , FUN = n_distinct) %>% 
                            arrange(treatment, patient) %>% 
                            pivot_wider(names_from = signif, values_from = Peptide) 
# check i dentical row before combining 
identical(herv_exp_pr_patient$patient,patient_response_peptides$patient) 

# add column to df 
herv_exp_pr_patient$ResponsePeptides <- patient_response_peptides$`TRUE`
herv_exp_pr_patient$ResponsePeptides[is.na(herv_exp_pr_patient$ResponsePeptides)] <- 0

# calculate response fraction 
herv_exp_pr_patient$ResponseFraction <- herv_exp_pr_patient$ResponsePeptides / herv_exp_pr_patient$sum * 100

# generate id 
herv_exp_pr_patient$id <- paste(herv_exp_pr_patient$patient, herv_exp_pr_patient$treatment, sep = '_')

# add expression information 
mean_exp_pr_patient$treatment <- ifelse(mean_exp_pr_patient$treatment=='Before Treatment',yes = 'Before', no = 'After') 
mean_exp_pr_patient$id <- paste(mean_exp_pr_patient$patient, mean_exp_pr_patient$treatment, sep = '_')
# create named vector 
mean_exp_pr_patient_vec <- setNames(mean_exp_pr_patient$mean_exp, mean_exp_pr_patient$id)
herv_exp_pr_patient$mean_exp <- mean_exp_pr_patient_vec[herv_exp_pr_patient$id]

# write to excel file 
write.xlsx(herv_exp_pr_patient, file=paste(dir,'data/HERV_expression_pr_patient.xlsx', sep = ''))

# ---------------------------------





# --------------------------
# Handle clinical data 
# --------------------------

head(clin_HHRH)
dim(clin_HHRH)
head(clin_SH)

# make compatible ID
clin_HHRH$study  <- sapply(strsplit(clin_HHRH$Pt.no, '_'), `[[`, 1)
clin_HHRH$patient <- sapply(strsplit(clin_HHRH$Pt.no, '_'), `[[`, - 1)
clin_HHRH$ID <- paste(clin_HHRH$study, clin_HHRH$patient, sep = '')

# convert stabile desies to non responder 
clin_HHRH$Responder[clin_HHRH$ID == 'RH124'] <- FALSE

# SH 
head(clin_SH)
clin_SH$patient <- gsub("PD", "SH", clin_SH$patient)
clin_SH$study <- 'SH'
table(clin_SH$Response)
clin_SH$Responder <- clin_SH$Response == 'CR'

# merge data from the two studies 
# ..................................
# wanted columns 
HHRH_short <- clin_HHRH[ ,c('ID', 'study', 'Responder', 'Risk_binary')]
colnames(HHRH_short) <- c('patient','study', 'Responder', 'Risk_binary')
SH_short <- clin_SH[ , c('patient','study', 'Responder', 'Risk_binary')]

clin <- rbind(HHRH_short,SH_short )
head(clin)
table(clin$patient)






# --------------------------
# Meta data 
# --------------------------

# Meta data faile for annotation column of expression plots 
meta_df <- clin[clin$patient %in% unique(final_exp_plot$patient), c("patient","Responder")]
meta_df$Responder <- ifelse(meta_df$Responder == TRUE,'Yes', 'No')
rownames(meta_df) <- meta_df$patient
meta_df$patient <- NULL 

# Add annotation of detection of any respons to HERV peptide 
head(lib_signif_patient)
lib_signif_patient[lib_signif_patient$Library.type == 'HERV' & lib_signif_patient$patient %in% rownames(meta_df), ]
lib_signif_patient$HERVRespons <- ifelse(lib_signif_patient$signif == 0, 'No','Yes')
# generate named vecter 
herv_signif_vec <- setNames(lib_signif_patient$HERVRespons, lib_signif_patient$patient )
# add to meta df 
meta_df$HERVRespons <- herv_signif_vec[rownames(meta_df)]






# ------------------------------------------------------------------------


# ############### #
#                 #
#      PLOT       #
#                 #
# ############### #


# ------------------------------------------------------------------------

setwd(paste(dir, 'results/PlotsFinalRevision', sep = ''))

# ------------------------------------
#           FIGURE 2
# ------------------------------------
# set breaks, limist and range 
my_breaks <- c(0,1,2,3,4,6,8,10)
my_range <-  c(.05, 1)
my_lim <- c(0,10)

# T cell reactivity heatmap - barcode heatmap  
# HLA specific plots 

hlaplot_signif <- function(df, hla){
    # select samples with the correct hla and significan values and tested for HERV lib (ensure correct hla type)
    samples <-  unique(df$sample[df$HLAallele == hla & barcode_dat$Library.type=="HERV"])
    peptides <- unique(df$Peptide[df$HLAallele == hla & df$signif == TRUE ]) #
    # df_select <- df[df$Peptide %in% peptides, ] 
    df_select <- df[df$sample %in% samples & df$Peptide %in% peptides, ] 
    # determin plorting order of signifficant
    df_select$signif <- factor(df_select$signif, levels = c('TRUE','FALSE'), labels=c("YES", "NO"))
    df_select$treatment <- factor(df_select$treatment, levels = c('Healthy','Before', 'After'), labels=c('Healthy donors', paste('Patients', 'Pre-AZA',sep='\n'), paste('Patients', 'Post-AZA',sep='\n')))
    df_select$Library.type <- factor(df_select$Library.type, levels = c('HERV','CTA','viral'), labels=c('HERV','CTA','Viral'))
    
    p <- ggplot(df_select, aes(sample, Peptide)) +
        geom_tile(aes(fill = signif, colour = signif, alpha = log_fold_change)) +
        scale_alpha_continuous('Log2FC',
                               range = my_range,
                               breaks = my_breaks,
                               limits = my_lim) +
        scale_color_manual(values = c("#F00000", "#bdbdbd")) +
        scale_fill_manual(values = c("#F00000", "#ABABA9"))+ 
        theme_bw(base_size = 12) +
        labs(x = "", y = "", title = hla) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 6),
              axis.text.y = element_text(size = 8),
              axis.ticks = element_blank(),
              panel.spacing = unit(0.2, "lines"),
              strip.text.y = element_text(angle = 0, hjust = 0),
              strip.background = element_rect(fill = 'white', colour = 'white'),
              legend.position = 'none',
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        facet_grid(Library.type ~ treatment, drop = T, space = "free", scales = "free")
    return(p)
}


# make legend plot
legend_plot <- function(x){
    ggplot(barcode_dat[1:10,], aes(sample, Peptide)) + 
        geom_tile(aes(alpha = log_fold_change), color =  x,fill = x) +
        scale_alpha_continuous('Log2FC',
                               range = my_range,
                               breaks = my_breaks,
                               limits = my_lim) +
        guides(alpha = guide_legend(nrow = 1,
                                    title.position = "left",
                                    label.position = "top",
                                    title.vjust = 1))
}


# make plot opbjects
p1 <- hlaplot_signif(barcode_dat, 'HLA-A01:01')
p2 <- hlaplot_signif(barcode_dat, 'HLA-A02:01')
p3 <- hlaplot_signif(barcode_dat, 'HLA-B07:02')
p4 <- hlaplot_signif(barcode_dat, 'HLA-B08:01')
# gelend opjects 
l <- legend_plot("#F00000")
l1 <- get_legend(l) 
l <- legend_plot("#ABABA9")
l2 <- get_legend(l) 

# Make one plot 
pdf('Figure2.pdf', width = 15, height  = 10)
ggdraw() +
    draw_plot(p1, .003, .72, .44, .28) +
    draw_plot(p2, .46, .47, .54, .53) +
    draw_plot(p3, 0, .41, .39, .3) +
    draw_plot(p4, .005, .075, .43, .32) +
    draw_plot(l2, -.037, -.07, .4, .2) +
    draw_plot(l1, -.037, -.034, .4, .2) +
    draw_plot_label(c("A", "B", "C", "D", "E"), c(0,  .46, 0,  0, .46), c(1, 1, .712, .395, .395), size = 15)
dev.off()



# TEST
# Log2FC of grey values 
# head(barcode_dat)
# barcode_dat[barcode_dat$HLAallele == 'HLA-A02:01' &
#                 barcode_dat$Peptide %in% c('YMRTLLDSI', 'YLKACLTVL') &
#                 barcode_dat$sample %in% c('BC100','HH12C1D1','SH7157C1D1', 'SH7163C1D1', 'SH7164C1D1', 'HH12C3D1', 'SH7157C4D28', 'SH7163C6D28','SH7169C4D1', 'HH18C2D5','BC303'),
#             c('sample','HLAallele','Peptide','count.1','log_fold_change','signif')]







# ------------------------------------

#           FIGURE 5 

# ------------------------------------

# Figure 5A
# --------------------------------
# box plot / beeswarm plot of mean erv expression accross all patients
# subset data 
plot_exp_pr_erv <- final_exp_plot[!final_exp_plot$patient  %in% c('SH7154', 'SH7162'), ]

# generate data frame for statistical plotting 
plot_exp_pr_erv_exponential <- plot_exp_pr_erv %>% mutate(mean_exp = 10^(mean_exp + 0.05))

round(p.adjust(0.008, "BH"), 3)
# plot data 
p1 <- ggplot(plot_exp_pr_erv, aes(x = treatment, y = mean_exp + 0.05)) +
    geom_quasirandom(aes(colour = treatment), size = 1, alpha = .6)  +
    geom_boxplot(aes(), alpha = 0) +
    geom_signif(data = plot_exp_pr_erv_exponential, 
                comparisons = list(c("Healthy", "After Treatment"), 
                                   c("Healthy", "Before Treatment")), 
                test = 'wilcox.test',
               # test = 't.test',
                map_signif_level = T, 
                textsize = 5,
               test.args = list(paired = F),
                y_position= c(3, 2.7),
                tip_length = 0) +
    geom_signif(data = plot_exp_pr_erv_exponential, 
                comparisons = list(c("Before Treatment", "After Treatment")), 
                test = 'wilcox.test',
               # test = 't.test',
                test.args = list(paired = TRUE),
                map_signif_level = T, 
                textsize = 5,
                y_position= 2.5,
                tip_length = 0) +
    scale_y_log10(breaks = c(.1, 1, 10, 100,1000), labels = c('0.1', '1', '10', '100', '1000')) + 
    annotation_logticks(sides = 'l', colour = '#969696') +
    scale_color_manual(values = c('#ABABA9', '#d4b9da','#AA3377')) +
    labs(y = '', x = '') + #HERV Expression (TPM)
    theme_classic( base_size = 15) +
    theme(legend.position = 'none') +
   # scale_x_discrete(labels=c("Healthy\ndonors",'Patients\nPre-AZA','Patients\nPost-AZA')) 
theme(axis.text.x = element_blank(),
      axis.text.y = element_blank()
      )
ggsave(p1, file = "fig_5_A.wilcox.test.pdf", , width = 6, height  = 6)

# calculate wilcox test 
for(p in c('SH7151', 'SH7163', 'SH7167', 'SH7161')){
    print(p)
    bee_df %>% subset(cycle %in% c('C1D1', 'C6D28') & patient == p) %>% # subset correct data 
        wilcox.test(mean_exp ~ cycle, data = .,paired = TRUE) %>% # run test 
        print() # print test output 
}

# run wilcoxonstest 
plot_exp_pr_erv %>% subset(treatment %in% c('Healthy', 'After Treatment')) %>% wilcox.test(mean_exp ~ treatment, data = .)
plot_exp_pr_erv %>% subset(treatment %in% c('Healthy', 'Before Treatment')) %>% wilcox.test(mean_exp ~ treatment, data = .)# run test 
plot_exp_pr_erv %>% subset(treatment %in% c('Before Treatment', 'After Treatment')) %>% wilcox.test(mean_exp ~ treatment, data = ., paired = TRUE)# run test 
    

# Figure 5B
# --------------------------------

# VIRKER MEN GIVER FORKERTE TAL ..... PASSER IKKE
# 
# herv_exp_pr_patient[!is.na(herv_exp_pr_patient$mean_exp), ]
# 
# p2 <- ggplot(herv_exp_pr_patient, aes(x = mean_exp, y = ResponseFraction)) +
#     geom_point(aes(colour = treatment), size = 3, alpha = .6) +
#     scale_y_log10()+#breaks = c(.1, 1, 10, 100,1000), labels = c('0.1', '1', '10', '100', '1000')) + 
#     scale_x_log10()+#breaks = c(.1, 1, 10, 100,1000), labels = c('0.1', '1', '10', '100', '1000')) + 
#     annotation_logticks(sides = 'lb', colour = '#969696') +
#     scale_color_manual(values = c('#d4b9da','#AA3377')) +
#     # labs(y = 'HERV Expression (TPM)', x = '') +
#     theme_classic( base_size = 15) +
#     theme(legend.position = 'none')
    


# Figure 5C 
# --------------------------------

# set plotting order 
exp_pr_erv_mean$HERVRespons <- factor(exp_pr_erv_mean$HERVRespons, levels = c('Yes', 'No'))
# exp dat for test
exp_pr_erv_mean_exponential <- exp_pr_erv_mean
exp_pr_erv_mean_exponential$mean_exp <- 10^(exp_pr_erv_mean_exponential$mean_exp + 0.01)

# plot data
p2 <- ggplot(exp_pr_erv_mean, aes(x = HERVRespons, y = mean_exp + 0.01)) +
    geom_quasirandom(aes(shape = HERVRespons, alpha = HERVRespons), fill = '#AA3377', colour = '#AA3377', size = 3,  stroke = 1)  +
    geom_boxplot(aes(), alpha = 0) +
    geom_signif(data = exp_pr_erv_mean_exponential,
                comparisons = list(c("Yes", "No")),
                test = 'wilcox.test',
                map_signif_level = TRUE,
                textsize = 5,
                y_position= 2,
                tip_length = 0) +
    scale_y_log10(breaks = c(.1, 1, 10, 100), labels = c('0.1', '1', '10', '100')) +
    annotation_logticks(sides = 'l', colour = '#969696') +
    scale_alpha_manual(values = c(.9, 1 )) +
    scale_shape_manual(values = c(21,1)) +
    labs(y = 'Mean HERV expression (TPM)', x = '') +
    theme_classic( base_size = 15) +
    theme(legend.position = 'none') +
    scale_x_discrete(labels=c('T cell\npositive','T cell\nnegative'))




# Figure 5D 
# --------------------------------
# logFC expression heatmap comparing HERV gene expression chanes between healthy doners and patients before treatment


# Logfold change 
# ---------------------------
log_fold = function(x,y){
    log2(x) - log2(y) }
# ---------------------------
# list of breaks for colour legends (make one comon legend across plots)
breaksList = seq(-7, 7, by = .2)
#breaksList = seq(-10, 10, by = .2)
# ---------------------------


# calculate mean expression value for healthy individuals 
mean_healthy <- aggregate(mean_exp ~ hugo_symbol, data = subset(exp_pr_patient, exp_pr_patient$treatment=="Healthy"), mean) 
colnames(mean_healthy) <- c("hugo_symbol", "Healthy")

# extract data for expression of HERVs before treatments 
HERV_expression_before <- subset(exp_pr_patient, exp_pr_patient$treatment=="Before Treatment")

# change to wide format 
HERV_expression_before <- spread(HERV_expression_before, patient, mean_exp) %>% left_join(mean_healthy, .)
HERV_expression_before$treatment <- NULL

# calculate logfold change comaparing mean healthy for alle patinets before treatment¨
foldchange_before_healthy <- data.frame()
count = 0
for (col in colnames(HERV_expression_before[,3:length(colnames(HERV_expression_before))])) {
    count=count+1
    name <- paste0(col,"foldchange",sep="_")
    name <- log_fold(HERV_expression_before[col]+0.001, HERV_expression_before$Healthy+0.001)
    if (count==1) {
        foldchange_before_healthy <- cbind(HERV_expression_before$hugo_symbol, name)
    }
    else {
        foldchange_before_healthy <- cbind(foldchange_before_healthy,name)
    }
}

foldchange_before_healthy <- column_to_rownames(foldchange_before_healthy, var = "HERV_expression_before$hugo_symbol")

# set meta df to relevant patients 
meta_df_1 <- meta_df[rownames(meta_df) %in% colnames(foldchange_before_healthy), ]

# make hetmap 
pmap1 <- as.ggplot(pheatmap(foldchange_before_healthy, 
                 fontsize_col=10,
                 fontsize_row=10,
                 annotation_col = meta_df_1,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(length(breaksList)),
                 breaks = breaksList,
                 legend = FALSE, 
                 annotation_legend = FALSE))



# Figure 5E 
# --------------------------------
# logFC expression heatmap comparing HERV gene expression chanes before and after treatment with AZA

# Remove all healthy samples and samples without "after treatment" sample 
plot_heat <- exp_pr_patient[exp_pr_patient$treatment != 'Healthy' &
                                     !exp_pr_patient$patient  %in% c('SH7154', 'SH7162'), ]
# spread data to wide format to be compatible for logfoldchange calculations 
plot_heat  <- spread(plot_heat, treatment, mean_exp)
# calculate logfold change before and after 
plot_heat$fold_change_b_a <- log_fold(plot_heat[ ,"After Treatment"] + 0.001, plot_heat[ ,"Before Treatment"] + 0.001)

# Generate plotting matrtix 
plot_heat_logfold <- spread(plot_heat[ ,c("patient", "hugo_symbol","fold_change_b_a")], 
                            patient, 
                            fold_change_b_a)
plot_heat_logfold <- column_to_rownames(plot_heat_logfold, var = 'hugo_symbol')

# set meta df to relevant patients 
meta_df_2 <- meta_df[rownames(meta_df) %in% colnames(plot_heat_logfold), ]

# Plot Heatmap 
pmap2 <- as.ggplot(pheatmap(plot_heat_logfold, 
                 fontsize_col=10, 
                 fontsize_row=10,
                 annotation_col = meta_df_2,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(length(breaksList)),
                 breaks = breaksList))


# check minimun and max values to know how large the saturation is 
min(plot_heat_logfold)
max(plot_heat_logfold)
min(foldchange_before_healthy)
max(foldchange_before_healthy)

# Make one plot 
# --------------------------------
pdf('Figure5_new.pdf', width = 14, height  = 18)
ggdraw() +
    draw_plot(p1, x=0.03, y=.725, width = .38, height = .26) +
    draw_plot(p2, x=.75, y=.7, width = .2, height = .28) +
    draw_plot(pmap1, x=0, y=0, width = .45, height = .69) +
    draw_plot(pmap2, x=.47, y=0, width = .55, height = .69) +
    draw_plot_label(c("A", "B", "C", "D", "E"), c(0,  .4, .75,  0, .47), c(1, 1, 1, .7, .7), size = 15)
dev.off()
















# ------------------------------------

#       SUPPLEMENTARY FIGURES

# ------------------------------------



# ------------------------------------
# Supplementary figure 1
# ------------------------------------
# Heatmap of T cell peptides screened all for each HLA alllele 

# function for plotting the heatmap 
hlaplot_all <- function(df, hla){
    # make HLA specific df 
    df_HLA <- df[df$HLAallele == hla, ]
    # order according to significance
    df_HLA <- df_HLA[order(df_HLA$signif, df_HLA$log_fold_change), ]
    # patients with HERV data 
    id_erv <- unique(df_HLA$sample[df_HLA$Library.type == 'HERV'])
    
    # determin plorting order of signifficant 
    df_HLA$treatment <- factor(df_HLA$treatment, levels = c('Healthy','Before', 'After'), labels=c('Healthy donor', 'Patient\nPre-AZA', 'Patient\nPost-AZA')) 
    df_HLA$Library.type <- factor(df_HLA$Library.type, levels = c('HERV','CTA','viral'), labels=c('HERV','CTA','Viral')) 
 
    p <- ggplot(df_HLA[which(df_HLA$sample %in% id_erv), ], aes(x = sample, y = pMHC)) + 
        geom_tile(aes(fill = signif, colour = signif, alpha = log_fold_change)) + #colour = signif
        scale_alpha_continuous('Log2FC',
                               range = c(.1,1),
                               breaks = my_breaks,
                               limits = my_lim) +
        scale_fill_manual(values = c("#ABABA9", "#F00000"))+
        scale_colour_manual(values = c("#DDDDDD00", "#F0000050"))+ # zero'es bake colour transparent / 50 percent 
        theme_bw(base_size = 12) +
        labs(x = "", y = "", title = hla) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 6),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              panel.spacing = unit(0.2, "lines"),
              strip.text.y = element_text(angle = 0, hjust = 0),
              strip.background = element_rect(fill = 'white', colour = 'white'),
              legend.position = 'none',
              panel.border = element_rect(colour = "#ABABA9"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        facet_grid(Library.type ~ treatment, drop = T, space = "free", scales = "free",
                   labeller = labeller(Library.type = c("CTA" = "", "HERV" = "","Viral" = "")))
     return(p)
}

# Make individual plotting objects 
p1 <- hlaplot_all(barcode_dat, 'HLA-A01:01')
p2 <- hlaplot_all(barcode_dat, 'HLA-A02:01')
p3 <- hlaplot_all(barcode_dat, 'HLA-B07:02')
p4 <- hlaplot_all(barcode_dat, 'HLA-B08:01')

# Make one plot  
pdf('FigureS1.pdf', width = 10, height  = 13)      
ggdraw() +
    draw_plot(p1, 0, .55, .47, .45) +
    draw_plot(p2, .5, .55, .47, .45) +
    draw_plot(p3, 0, .1, .47, .45) +
    draw_plot(p4, .5, .1, .47, .45) +
    draw_plot(l2, -.037, -.06, .4, .2) +
    draw_plot(l1, -.037, -.034, .4, .2) +
    draw_text(c("HERV","CTA","Viral","HERV","CTA","Viral"), 
              c(.47,  .468, .468,.97,  .968, .968), 
              c(0.8,.63,.62,0.8,.66,.62), 
              size = 10)+
    draw_text(c("HERV","CTA","Viral","HERV","CTA","Viral"), 
              c(.47,  .468, .468,.97,  .968, .968), 
              c(0.36,.178,.168,0.36,.178,.168), 
              size = 10)+
    
    draw_plot_label(c("A", "B", "C", "D"), c(0,  .5, 0,  .5), c(1, 1, .55, .55), size = 15)
dev.off()





 # ------------------------------------
# Supplementary figure 5
# ------------------------------------
# Transcript level expression plots 

# Expression heatmap (TPM) on transcript level 
p <- ggplot(plot_exp, aes(patient, accession_number)) + 
    geom_tile(aes(fill = mean_exp, colour = "0")) + 
    scale_fill_distiller('TPM', 		                     
                         palette = "Spectral", 
                         trans = "log", 
                         na.value="#0570b0", 
                         breaks = c(.1, 1, 10, 100), 
                         labels = c('0.1', '1', '10', '100')) +
    scale_colour_manual(values = "white") +
    theme_grey(base_size = 9) + 
    labs(x = "", y = "") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size=8),
          axis.text.y = element_text(hjust = 1,vjust = .5),
          axis.ticks = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 10), 
          strip.text.y = element_text(angle = 0, hjust = 0,vjust = .5), 
          strip.background = element_rect(fill = 'white'),
          legend.text = element_text(size = 9),
          legend.margin = margin(0,0,0,0),
          legend.spacing = unit(0, 'cm'),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(order = 2,
                                  title = '',
                                  title.position = 'top',
                                  barwidth = 10),
           colour = guide_legend(title = 'TPM',
                                 title.position = 'top',
                                 override.aes = list(fill = "#0570b0"),
                                 order = 1,
                                 label.hjust = .5,
                                 label.vjust = 1,
                                 label.position = 'bottom')) +
    facet_grid(hugo_symbol ~ treatment, drop = T, space = "free", scales = "free",
               labeller = labeller(treatment = c('Healthy' = 'Healthy donor',
                                                 'Before Treatment' = 'Patient Pre-AZA', 
                                                 'After Treatment' = 'Patient Post-AZA')))
pdf('FigureS5.pdf', width = 7, height = 11)
print(p)
dev.off()




# ------------------------------------
# Supplementary figure 6
# ------------------------------------
my_breaks2 = c(0,0.05,0.5,5,50,500)

# Figure S6A 
# ------------------------------------
# Exppression TMPs of selected transcripts 
p1 <- ggplot(final_exp_plot, aes(patient, hugo_symbol)) + 
    geom_tile(aes(fill = mean_exp, colour = "0")) + 
    scale_fill_distiller('TPM', 		                     
                         palette = "Spectral", 
                         trans = "log", 
                         na.value="#0570b0", 
                         breaks = c(.1, 1, 10, 100), 
                         labels = c('0.1', '1', '10', '100')) +
    scale_colour_manual(values = "white") +
    theme_grey(base_size = 9) + 
    labs(x = "", y = "") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size=8),
          axis.text.y = element_text(hjust = 0),
          axis.ticks = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 13), 
          strip.background = element_rect(fill = 'white'),
          legend.text = element_text(size = 9),
          legend.margin = margin(0,0,0,0),
          legend.spacing = unit(0, 'cm'),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(order = 2,
                                  title = '',
                                  title.position = 'top',
                                  barwidth = 10),
           colour = guide_legend(title = 'TPM',
                                 title.position = 'top',
                                 override.aes = list(fill = "#0570b0"),
                                 order = 1,
                                 label.hjust = .5,
                                 label.vjust = 1,
                                 label.position = 'bottom')) +
    facet_grid( ~ treatment, drop = T, space = "free", scales = "free",
                labeller = labeller(treatment = c('Healthy' = 'Healthy donor',
                                                  'Before Treatment' = 'Patient Pre-AZA', 
                                                  'After Treatment' = 'Patient Post-AZA')))


# Figure S6B 
# ------------------------------------
# samples with all treatment cyckles 
patients <- unique(erv_exp_dat$patient[erv_exp_dat$treatment %in% c('C1D8', 'C1D28')])
# extract DF with full cycles 
completeCycles_plot_exp <- erv_exp_dat[erv_exp_dat$id != 'ViralMimicry' & # no Viral mimery 
                                           erv_exp_dat$treatment != 'Healthy' & # exclude unvanted healthy donor samples
                                           erv_exp_dat$accession_number %in% transcript_ranks$accession_number[transcript_ranks$rank == 1] &
                                           erv_exp_dat$hugo_symbol %in% HERV_library_postive$hugo_symbol & # only hugo symbols where T cell responce was relevant 
                                           erv_exp_dat$patient %in% patients, ]  
table(completeCycles_plot_exp$treatment)
head(completeCycles_plot_exp)

# transform data for significante calculations 
completeCycles_exp_plot_exponential <- completeCycles_plot_exp
completeCycles_exp_plot_exponential$mean_exp <- 10^(completeCycles_exp_plot_exponential$mean_exp +.05)

# generate plot 
p4 <- ggplot(completeCycles_plot_exp, aes(x = treatment, y = mean_exp +.05)) +
    geom_quasirandom(aes(colour = treatment), size = 1, alpha = .8)  +
    geom_boxplot(aes(), alpha = 0) +
    geom_signif(data = completeCycles_exp_plot_exponential,
                comparisons = list(c("Before Treatment", "After Treatment")), 
                test = 'wilcox.test', 
                test.args = list(paired = TRUE),
                map_signif_level = T, 
                textsize = 2,
                vjust = -.2, 
                y_position = 2,
                tip_length = 0) +
    scale_y_log10(breaks = c(.1, 1, 10, 100), labels = c('0.1', '1', '10', '100')) + 
    annotation_logticks(sides = 'l', colour = '#969696') +
    scale_color_manual(values = c('#d4b9da', '#c994c7', '#df65b0','#AA3377')) +
    labs(y = 'HERV Expression (TPM)', x = '') +
    theme_classic( base_size = 10) +
    theme(legend.position = 'none') +
    scale_x_discrete(labels=c("Before",'C1D8','C1D28', "After"))

# WiLCOXONS test 
wilcox.test(completeCycles_plot_exp$mean_exp[completeCycles_plot_exp$treatment == 'Before Treatment'],
            completeCycles_plot_exp$mean_exp[completeCycles_plot_exp$treatment == 'After Treatment'],
            paired = TRUE)



# Figure S6C 
# ------------------------------------
# Selected beeswarm df generation 
bee_df <- erv_exp_dat[erv_exp_dat$cycle != 'healthy' & 
                          erv_exp_dat$id != 'ViralMimicry' & 
                          erv_exp_dat$hugo_symbol %in% HERV_library_postive$hugo_symbol &
                          erv_exp_dat$patient %in% c('SH7151','SH7167','SH7163','SH7161'), ]
# Define plotting order of factors 
bee_df$patient <- factor(bee_df$patient, levels=c("SH7151", "SH7163", "SH7167", "SH7161"))
bee_df$cycle <- factor(bee_df$cycle, levels=c('C6D28', 'C1D28','C1D8', 'C1D1'))

# make plotting funktion 
plot_bee <- function(df){
    df_exponential <- df
    df_exponential$mean_exp <- 10^(df_exponential$mean_exp +.01)
    
    p <- ggplot(df, aes(x = cycle, y = mean_exp + .01)) +
        geom_quasirandom(aes(colour = mean_exp), size = 1) + 
        geom_boxplot(aes(), alpha = 0) +
        geom_signif(data = df_exponential,
                    comparisons = list(c("C1D1", "C6D28")),
                    test = 'wilcox.test',
                    test.args = list(paired = TRUE),
                    map_signif_level = TRUE,
                    textsize = 2,
                    #hjust = -1,
                    # hjust = 5,
                    # vjust = -.5,
                    # vjust = 35,
                    y_position = 2.5,
                    tip_length = 0) +
        scale_colour_distiller(name = "TPM",
                               palette = "Spectral", 
                               trans = "log", 
                               na.value = "#0570b0", 
                               breaks = my_breaks2, 
                               labels = my_breaks2, 
                               guide=FALSE) +
        ylab("HERV Expression (TMP)") +
        theme_minimal(base_size = 15)+
        facet_wrap(~ patient, ncol = 2, drop = TRUE) + 
        scale_y_log10(breaks =  c(.1, 1, 10, 100), 
                      labels = c('0.1', '1', '10', '100')) +
        coord_flip() +
        labs(x = "") +
        theme(strip.background = element_rect(fill = 'white'), 
              panel.spacing = unit(0.5, "lines"),
              text = element_text(size = 10),
              axis.text = element_text(size = 9),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    return(p)
}

# Make plot objects 
p3_1 <- plot_bee(subset(bee_df, bee_df$patient %in% c('SH7151','SH7163')))
p3_2 <- plot_bee(subset(bee_df, bee_df$patient %in% c('SH7167','SH7161')))
               

# paired will coxons (Mann-Whitney) tests for each patient 
for(p in c('SH7151', 'SH7163', 'SH7167', 'SH7161')){
    print(p)
    bee_df %>% subset(cycle %in% c('C1D1', 'C6D28') & patient == p) %>% # subset correct data 
                    wilcox.test(mean_exp ~ cycle, data = .,paired = TRUE) %>% # run test 
                    print() # print test output 
}


# Make one plot
# ------------------------------------
pdf('FigureS6.pdf', width = 10, height  = 9)
ggdraw() +
    draw_plot(p1, 0, 0, .6, 1) +
    draw_plot(p3_1, .6, .43, .38, .22) +
    draw_plot(p3_2, .6, .27, .38, .17) +
    draw_plot(p4, .6, .66, .395, .3) +
    draw_plot_label(c("A", "B", "C"), c(0, .62, .62), c(1, 1, .66), size = 15)
dev.off()













# ------------------------------------
# Supplementary figure S7
# ------------------------------------
# Transcript level expression plots 


# Figure S7A
#.................
# plot of sum of expression of ervs acoress all patients (each dot is an HERV) 
exp_pr_erv_exponential <- exp_pr_erv %>% mutate(mean_exp = 10^(mean_exp + 0.01))

p1 <- ggplot(exp_pr_erv, aes(x = treatment, y = mean_exp + 0.01)) +
    geom_quasirandom(aes(colour = treatment), size = 3, alpha = .8)  +
    geom_boxplot(aes(), alpha = 0) +
    geom_signif(data = exp_pr_erv_exponential, 
                comparisons = list(c("Healthy", "After Treatment"),
                                   c("Healthy", "Before Treatment"),
                                   c("Before Treatment", "After Treatment")), 
                test = 'wilcox.test', 
                test.args = list(paired=F, exact = F),
                map_signif_level = F, 
                textsize = 5,
                # vjust = .5, 
                y_position= c(3.5, 3.1, 2.9),
                tip_length = 0) +
    scale_y_log10(breaks = c(.1, 1, 10, 100,1000), labels = c('0.1', '1', '10', '100', '1000')) + 
    annotation_logticks(sides = 'l', colour = '#969696') +
    scale_color_manual(values = c('#e7e1ef', '#d4b9da','#AA3377')) +
    labs(y = 'HERV expression sum across individual (TPM)', x = '') +
    theme_classic( base_size = 15) +
    theme(legend.position = 'none') +
    scale_x_discrete(labels=c("Healthy\ndonors","Patients\nPre-AZA","Pateints\nPost-AZA")) 


# Figure S7B
#.................
# aggregate according to sum of expression for all HERVS pr patient 
exp_pr_patient_2 <- aggregate(mean_exp ~ patient + treatment,  data = exp_pr_patient,FUN = sum) %>% 
                        subset(!patient  %in% c('SH7154', 'SH7162')) # remove patients without after sample 

# generate data frame for statistical plotting 
exp_pr_patient_2_exponential <- exp_pr_patient_2 %>% mutate(mean_exp = ifelse(mean_exp > 300, 
                                                                              yes = 10^(300 + 0.01), # if values are above 300 TPM make then 300 to avoid infinite values 
                                                                              no = 10^(mean_exp + 0.01))) # exponential calculation
# plot data 
p2 <- ggplot(exp_pr_patient_2, aes(x = treatment, y = mean_exp + 0.01)) +
    geom_quasirandom(aes(colour = treatment), size = 3, alpha = .8)  +
    geom_boxplot(aes(), alpha = 0) +
    geom_signif(data = exp_pr_patient_2_exponential, 
                comparisons = list(c("Healthy", "After Treatment"), 
                                   c("Healthy", "Before Treatment")), 
                test = 'wilcox.test', 
                map_signif_level = F, 
                textsize = 5,
                test.args = list(paired=F, exact = F),
                y_position= c(3, 2.8),
                tip_length = 0) +
    geom_signif(data = exp_pr_patient_2_exponential, 
                comparisons = list(c("Before Treatment", "After Treatment")), 
                test = 'wilcox.test',
                test.args = list(paired = TRUE, exact = F),
                map_signif_level = F, 
                textsize = 5,
                y_position= 2.7,
                tip_length = 0) +
    scale_y_log10(breaks = c(.1, 1, 10, 100,1000), labels = c('0.1', '1', '10', '100', '1000'))+
    annotation_logticks(sides = 'l', colour = '#969696') +
    scale_color_manual(values = c('#e7e1ef', '#d4b9da','#AA3377')) +
    labs(y = 'sum of HERV Expression pr individual (TPM)', x = '') +
    theme_classic( base_size = 15) +
    theme(legend.position = 'none') +
    scale_x_discrete(labels=c("Healthy donors",'Patients\nPre-AZA','Patients\nPost-AZA')) 

# calculate wilcox test 
exp_pr_patient_2 %>% subset(treatment %in% c('Healthy', 'After Treatment')) %>% wilcox.test(mean_exp ~ treatment, data = .)
exp_pr_patient_2 %>% subset(treatment %in% c('Healthy', 'Before Treatment')) %>% wilcox.test(mean_exp ~ treatment, data = .)
exp_pr_patient_2 %>% subset(treatment %in% c('Before Treatment', 'After Treatment')) %>% wilcox.test(mean_exp ~ treatment, data = ., paired = TRUE)


# Figure S7C
#.................

# Make line plot 
exp_pr_patient_2_2 <- subset(exp_pr_patient_2, treatment != 'Healthy')
# plot data 
p3 <- ggplot(exp_pr_patient_2_2, aes(x = treatment, y = mean_exp)) + 
    geom_point(size = 3, alpha = .8, aes(colour = patient)) +
    geom_line(aes(group = patient,colour = patient)) +
    geom_boxplot(alpha = 0) +
    scale_colour_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#000000','#bdbdbd','#02818a', '#e7298a')) +
    geom_signif(data = exp_pr_patient_2_2, 
                comparisons = list(c("Before Treatment", "After Treatment")), 
                test = 'wilcox.test',
                test.args = list(paired = TRUE),
                map_signif_level = TRUE, 
                textsize = 5,
                y_position= 450,
                tip_length = 0) +
    theme_classic( base_size = 15) +
    labs(y = 'sum of HERV Expression pr patient (TPM)', x = '') +
    scale_x_discrete(labels=c('Patients\nPre-AZA','Patients\nPost-AZA')) 

# calculate wilcox test 
exp_pr_patient_2_2 %>% subset(treatment %in% c('Before Treatment', 'After Treatment')) %>% wilcox.test(mean_exp ~ treatment, data = ., paired = TRUE)



pdf('FigureS7.pdf', width = 10, height  = 10)
ggdraw() +
    draw_plot(p1, x = 0, y = .5, width = .45, height =  .45) +
    draw_plot(p2, x = .5, y = .5, width = .45, height =  .45) +
    draw_plot(p3, x = 0, y = 0, width = .6, height =  .45) +
    draw_plot_label(c("A", "B", "C"), c(0, .5, 0), c(1, 1, .5), size = 15)
dev.off()
















# ------------------------------------
# Supplementary figure S9
# ------------------------------------
#  APM genes 
# ------------------------------------

# Fig S9A
# ------------------------------------
# heatmap of antigen presentation pathway (AMP) genes 
# subset data 
PLOTTING_data_APM <- APM_expression %>% aggregate(mean_exp ~ treatment + Hugo + Patient, data = ., sum)

plot_APM <- subset(PLOTTING_data_APM, treatment %in% c('Healthy', 'Before Treatment', 'After Treatment'))
# plot heatmap   
p1 <- ggplot(plot_APM, aes(Patient, Hugo)) + 
    geom_tile(aes(fill = mean_exp, colour = "0")) + 
    scale_fill_distiller(name = "TPM", 
                         palette = "Spectral", 
                         trans = "log", 
                         na.value="#0570b0", 
                         breaks = c(0,0.05,0.5,5,50,1000), 
                         labels =c("0","0.05","0.5","5","50","1000") ) +
    scale_colour_manual(values = "white") +
    theme_grey(base_size = 9) + 
    labs(x = "", y = "") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.text.y = element_text(hjust = 0),
          axis.ticks = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 13), 
          strip.background = element_rect(fill = 'white'),
          legend.text = element_text(size = 9),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(order = 1), colour = guide_legend(title = NULL, override.aes = list(fill = "#0570b0"), order = 2)) +
    facet_grid( ~ treatment, drop = T, space = "free", scales = "free",
                labeller = labeller(treatment = c('Healthy' = 'Healthy donor',
                                                  'Before Treatment' = 'Patient Pre-AZA', 
                                                  'After Treatment' = 'Patient Post-AZA')))

# Fig S9B
# ------------------------------------
# plot avarage dot plot 
p2 <- ggplot(plot_APM, aes(x = treatment, y = mean_exp + 0.01)) +
    geom_quasirandom(aes(colour = treatment), size = 1, alpha = .6)  +
    geom_boxplot(aes(), alpha = 0) +
    geom_signif(data = exp_pr_erv_exponential, 
                comparisons = list(c("Healthy", "After Treatment"),
                                   c("Healthy", "Before Treatment"),
                                   c("Before Treatment", "After Treatment")), 
                test = 'wilcox.test', 
                test.args = list(paired=TRUE),
                map_signif_level = TRUE, 
                textsize = 3,
                # vjust = .5, 
                y_position= c(4.4, 3.9, 3.7),
                tip_length = 0) +
    scale_y_log10(breaks = c(.1, 1, 10, 100,1000, 10000), labels = c('0.1', '1', '10', '100', '1000','10000')) + 
    annotation_logticks(sides = 'l', colour = '#969696') +
    scale_color_manual(values = c('#ABABA9', '#d4b9da','#AA3377')) +
    labs(y = 'Sum of APM Expression (TPM)', x = '') +
    theme_classic(base_size = 13) +
    theme(legend.position = 'none') +
    scale_x_discrete(labels=c("Healthy\ndonors","Patients\nPre-AZA","Pateints\nPost-AZA")) 

pdf('FigureS9.pdf', width = 10, height  = 8)
ggdraw() +
    draw_plot(p1, 0,  0, .57, .99) +
    draw_plot(p2, .62, .46, .38, .5) +
    draw_plot_label(c("A", "B"), c(0, .62), c(1, 1), size = 15)
dev.off()
















# -------------------------------------------------------------------------------------------------
# TEST TEST TEST  TEST TEST TEST  TEST TEST TEST  TEST TEST TEST  TEST TEST TEST  TEST TEST TEST 
# -------------------------------------------------------------------------------------------------










































































#----------------------------------------------------

## Annie are adding some stuf ###

#----------------------------------------------------




#--------------------------------------------------------------------------------------------------------
##                                          CTA genes 
#--------------------------------------------------------------------------------------------------------



#CTA_expression$treatment[CTA_expression$Cycle=='C1D' & CTA_expression$Patient=='SH7161']<- as.factor('Before Treatement')
#unique(CTA_expression$Cycle[CTA_expression$Patient=='SH7161'])
# sum expression pr trancript 

#unique(CTA_expression$treatment[CTA_expression$Patient == "SH7161"])
CTA_expression_sum <-  aggregate(mean_exp ~ Patient+Hugo+treatment,  data = CTA_expression, sum) %>% 
                            mutate(mean_exp = as.numeric(mean_exp))
CTA_plot <- subset(CTA_expression_sum, treatment %in% c("Healthy","Before Treatment","After Treatment"))


### plot tPM 
p1 <- ggplot(CTA_plot, aes(Patient, Hugo)) + 
  geom_tile(aes(fill = mean_exp, colour = "0")) + 
  scale_fill_distiller('TPM', 		                     
                       palette = "Spectral", 
                       trans = "log", 
                       na.value="#0570b0", 
                       breaks = c(.1, 1, 10, 100), 
                       labels = c('0.1', '1', '10', '100')) +
  scale_colour_manual(values = "white") +
  theme_grey(base_size = 9) + 
  labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size=11),
        axis.text.y = element_text(hjust = 0),
        axis.ticks = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 13), 
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 12),
        legend.margin = margin(0,0,0,0),
        legend.spacing = unit(0, 'cm'),
        legend.position = "bottom") +
  guides(fill = guide_colourbar(order = 2,
                                title = '',
                                title.position = 'top',
                                barwidth = 10),
         colour = guide_legend(title = 'TPM',
                               title.position = 'top',
                               override.aes = list(fill = "#0570b0"),
                               order = 1,
                               label.hjust = .5,
                               label.vjust = 1,
                               label.position = 'bottom')) +
  facet_grid( ~ treatment, drop = T, space = "free", scales = "free")



# logfold change 
#-----------------------------

CTA_expression_before_after_df <- CTA_expression_sum %>% 
                                        filter(treatment %in% c("Before Treatment", "After Treatment")) %>% 
                                        pivot_wider(values_from = mean_exp, names_from = treatment) %>% 
                                        mutate(fold_change_b_a = log_fold(`After Treatment`+0.001,`Before Treatment`+0.001))
    
####------------------------------------------------------------------------
# fold change before vs after 
#######------------------------------------------------------------------------
# CTA_p2 <- ggplot(CTA_expression_before_after_df,
#                  aes(Patient, Hugo)) + 
#   geom_tile(aes(fill = fold_change_b_a)) + 
#   scale_fill_distiller(name = "log2 fold change", palette = "Spectral", breaks = c(-10,-5,-3,-1,0,1,3,5,10), labels =c("-10","-5","-3","-1","0","1","3","5","10")) +
#   scale_colour_manual(values = "white") +
#   theme_grey(base_size = 9) + 
#   labs(x = "", y = "") + 
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete(expand = c(0, 0)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
#         axis.text.y = element_text(hjust = 0),
#         axis.ticks = element_blank(),
#         panel.spacing = unit(0.1, "lines"),
#         strip.text.x = element_text(size = 13), 
#         strip.background = element_rect(fill = 'white'),
#         legend.text = element_text(size = 9),
#         legend.position = "bottom") 
# ggsave(CTA_p2 ,file = "/Volumes/SUND/Public/T-cells-and-cancer/SRH group/Group members/Sunil+AM/HERV/results/plots/Annie_adding/Paper_plots_CTA/CTA_before_after_logfold_exp.pdf", height = 20, width = 8 )
  #  guides(fill = guide_colourbar(order = 1), colour = guide_legend(title = NULL, override.aes = list(fill = "#0570b0"), order = 2)) +
 # facet_grid( ~ treatment, drop = T, space = "free", scales = "free")
###### for pheatmap 
CTA_expression_before_after_pheat <- CTA_expression_before_after_df[, c("Patient", "Hugo","fold_change_b_a")]
CTA_expression_before_after_pheat$fold_change_b_a <- as.numeric(CTA_expression_before_after_pheat$fold_change_b_a)
CTA_expression_before_after_pheat <- spread(CTA_expression_before_after_pheat,Patient,fold_change_b_a)
rownames(CTA_expression_before_after_pheat) <- CTA_expression_before_after_pheat$Hugo
CTA_expression_before_after_pheat$Hugo <- NULL
CTA_expression_before_after_pheat <- CTA_expression_before_after_pheat[!colnames(CTA_expression_before_after_pheat) %in% c("SH7154","SH7162")]
## calclulate mean and take all with mean above 0 out 
up_down_CTA_before_after <- as.data.frame(rowMeans(CTA_expression_before_after_pheat)[rowMeans(CTA_expression_before_after_pheat) > 1 ])
CTA_expression_before_after_pheat <- CTA_expression_before_after_pheat[rownames(CTA_expression_before_after_pheat) %in% rownames(up_down_CTA_before_after),]
pmap_CTA_before_after <- pheatmap(CTA_expression_before_after_pheat, fontsize_col=14, fontsize_row=14)
ggsave(pmap_CTA_before_after ,
       file = "/Volumes/SUND/Public/T-cells-and-cancer/SRH group/Group members/Sunil+AM/HERV/results/plots/Annie_adding/Paper_plots_CTA/pmap_CTA_fold_change_before_after_subset.pdf", height = 6, width = 8 )


####------------------------------------------------------------------------
######### fold change before vs healthy ########
#######------------------------------------------------------------------------
mean_healthy <- aggregate(mean_exp ~ Hugo, data = subset(CTA_expression_sum, CTA_expression_sum$treatment=="Healthy"), mean) 
colnames(mean_healthy) <- c("Hugo", "Healthy")
CTA_expression_before <- subset(CTA_expression_sum, CTA_expression_sum$treatment=="Before Treatment")
 
CTA_expression_before <- spread(CTA_expression_before, Patient,mean_exp) %>% left_join(mean_healthy, .)
CTA_expression_before$treatment <- NULL

# calculate logfold change comaparing mean healthy for alle patinets before treatment¨
foldchange_before_healthy <- data.frame()
count=0
for (col in colnames(CTA_expression_before[,3:length(colnames(CTA_expression_before))])) {
  count=count+1
  name <- paste0(col,"foldchange",sep="_")
  name <- log_fold(CTA_expression_before[col]+0.001, CTA_expression_before$Healthy+0.001)
  if (count==1) {
    foldchange_before_healthy <- cbind(CTA_expression_before$Hugo, name)
  }
  else {
  foldchange_before_healthy <- cbind(foldchange_before_healthy,name)
  }
  
}
# 
rownames(foldchange_before_healthy) <- foldchange_before_healthy$`CTA_expression_before$Hugo`
foldchange_before_healthy$`CTA_expression_before$Hugo` <- NULL


foldchange_before_healthy[foldchange_before_healthy[,1:length(colnames(foldchange_before_healthy))]>5] <- 5 
foldchange_before_healthy[foldchange_before_healthy[,1:length(colnames(foldchange_before_healthy))] < (-5)]  <- -5

## calclulate mean and take all with mean above 0 out 
up_down_CTA_before_healthy <- as.data.frame(rowMeans(foldchange_before_healthy)[rowMeans(foldchange_before_healthy) > 0 ])
CTA_foldchange_before_healthy_pheat <- foldchange_before_healthy[rownames(foldchange_before_healthy) %in% rownames(up_down_CTA_before_healthy),]

pmap_CTA_before_healthy <- pheatmap(CTA_foldchange_before_healthy_pheat, fontsize_col=14 , fontsize_row=14)

ggsave(pmap_CTA_before_healthy,
       file = "/Volumes/SUND/Public/T-cells-and-cancer/SRH group/Group members/Sunil+AM/HERV/results/plots/Annie_adding/Paper_plots_CTA/pmap_CTA_fold_change_before_healthy_subset.pdf", height = 12, width = 8 )











# Addition to Figure 5A 
# INDIVIDUAL BOX PLOTS ... do not use 
# expression of HERVs in individual patients / doners 
final_exp_plot_exponential <- final_exp_plot %>% mutate(mean_exp = 10^(mean_exp + 0.01))

p <- ggplot(final_exp_plot, aes(x = treatment, y = mean_exp + 0.01)) +
    geom_quasirandom(aes(colour = patient), size = 3, alpha = .8)  +
    geom_boxplot(aes(), alpha = 0) +
    geom_signif(data = final_exp_plot_exponential, 
                comparisons = list(c("Healthy", "After Treatment")), 
                test = 'wilcox.test', 
                map_signif_level = T, 
                textsize = 5,
                vjust = .5, 
                y_position= 3,
                tip_length = .0002) +
    scale_y_log10(breaks = c(.1, 1, 10, 100,1000), labels = c('0.1', '1', '10', '100', '1000')) + 
    annotation_logticks(sides = 'l', colour = '#969696') +
    labs(y = 'HERV Expression (TPM)', x = '') +
    theme_classic( base_size = 15) +
    facet_wrap(~patient, ncol = 4)+
    theme(legend.position = 'none') +
    scale_x_discrete(labels=c("Healthy\nDoners",'Patients\nPre-AZA','Patients\nPost-AZA')) 





# Supplementary figure XXX
# ------------------------------------
# Heatmap of expression - all transcripst all samples 
# Also devided into different patients and cycles 

erv_exp_dat <- erv_exp_dat[!erv_exp_dat$patient=="GSE69905_1",]

p1 <- ggplot(erv_exp_dat[erv_exp_dat$id != 'ViralMimicry', ], aes(patient, accession_number)) + 
    geom_tile(aes(fill = mean_exp, colour = "0")) + 
    scale_fill_distiller(name = "TPM", palette = "Spectral", trans = "log", na.value="#0570b0", breaks = my_breaks, labels = my_breaks) +
    scale_colour_manual(values = "white") +
    theme_grey(base_size = 9) + 
    labs(x = "", y = "") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.ticks = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text.y = element_text(angle = 0, hjust = 0),
          strip.text.x = element_text(size = 13), 
          strip.background = element_rect(fill = 'white'),
          legend.text = element_text(size = 9),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(order = 1), colour = guide_legend(title = NULL, override.aes = list(fill = "#0570b0"), order = 2)) +
    facet_grid( hugo_symbol~ treatment, drop = T, space = "free", scales = "free")

print(p1)


# additional expression plots emphasizing individual patients 

# heatmap 
p2 <- ggplot(erv_exp_dat[erv_exp_dat$cycle != 'healthy' & erv_exp_dat$id != 'ViralMimicry', ], aes(cycle, accession_number)) + 
    geom_tile(aes(fill = mean_exp, colour = "0")) + 
    scale_fill_distiller(name = "TPM", palette = "Spectral", trans = "log", na.value = "#0570b0", breaks = my_breaks, labels = my_breaks) +
    scale_colour_manual(values = "white") +
    theme_grey(base_size = 9) + 
    labs(x = "", y = "") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = .5),
          axis.ticks = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text.y = element_text(angle = 0, hjust = 0),
          strip.text.x = element_text(angle = 90, hjust = .5),
          strip.background = element_rect(fill = 'white'),
          legend.position="none") +
    # guides(fill = guide_colourbar(order = 1), colour = guide_legend(title = NULL, override.aes = list(fill = "#0570b0"), order = 2)) +
    facet_grid(hugo_symbol ~ patient, drop = T, space = "free", scales = "free")
# print(p2)



### calculate average of before and after treatment for each HERv gene from figure 4A (p1) 
average_pr_gene_treatment <- final_exp_plot %>% 
    aggregate(mean_exp ~ treatment + hugo_symbol, data = ., mean) %>% 
    pivot_wider(values_from = mean_exp,names_from = treatment)
# Make heatmap
p1 <- ggplot(average_pr_gene_treatment, aes(treatment, hugo_symbol)) +
  geom_tile(aes(fill = mean_exp, colour = "0")) +
  scale_fill_distiller(name = "TPM",
                       palette = "Spectral",
                       trans = "log",
                       na.value="#0570b0",
                       breaks = c(0,0.05,0.5,5,50,1000),
                       labels =c("0","0.05","0.5","5","50","1000") ) +
  scale_colour_manual(values = "white") +
  theme_grey(base_size = 9) +
  labs(x = "", y = "") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 0),
        axis.ticks = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 9),
        legend.position = "bottom") +
  guides(fill = guide_colourbar(order = 1), colour = guide_legend(title = NULL, override.aes = list(fill = "#0570b0"), order = 2)) +
  facet_grid( ~ treatment, drop = T, space = "free", scales = "free")

# ggsave(meanframe_genes,file = "meanframe_genes.pdf", height = 10, width = 6 )

# write to excel
write.csv(mean_frame, file='/Volumes/vet/Public/Afdeling-for-Immunologi-Vaccinologi/SRH group/Group members/Sunil+AM/HERV/results/plots/Annie_adding/mean_frame.csv')



# Supplementary figure 3
# ------------------------------------
# Viral response pathway expression 
p4 <- ggplot(erv_exp_dat[erv_exp_dat$id == 'ViralMimicry' & 
                             erv_exp_dat$treatment %in% c('Healthy', 'Before Treatment', 'After Treatment') &
                             !grepl('GMP|MEP|MGP', erv_exp_dat$patient), ],
             aes(patient, hugo_symbol)) + 
    geom_tile(aes(fill = mean_exp, colour = "0")) + 
    scale_fill_distiller(name = "TPM", palette = "Spectral", trans = "log", na.value="#0570b0", breaks = my_breaks, labels = my_breaks) +
    scale_colour_manual(values = "white") +
    theme_grey(base_size = 9) + 
    labs(x = "", y = "") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.text.y = element_text(hjust = 0),
          axis.ticks = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 13), 
          strip.background = element_rect(fill = 'white'),
          legend.text = element_text(size = 9),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(order = 1), colour = guide_legend(title = NULL, override.aes = list(fill = "#0570b0"), order = 2)) +
    facet_grid( ~ treatment, drop = T, space = "free", scales = "free")



p5 <- ggplot(erv_exp_dat[erv_exp_dat$id == 'ViralMimicry' & 
                             erv_exp_dat$treatment %in% c('Healthy', 'Before Treatment', 'After Treatment') & 
                             !grepl('GMP|MEP|MGP', erv_exp_dat$patient), ], 
             aes(x = treatment, y = mean_exp + 0.1)) +
    geom_quasirandom(aes(colour = treatment), size = 1)  +
    geom_boxplot(aes(), alpha = 0) +
    scale_y_log10(breaks = my_breaks) + 
    scale_color_manual(values = c('#d9d9d9', '#969696', '#525252')) +
    labs(y = 'VRP Expression (TPM)', x = '') +
    theme(legend.position = 'none') +
    scale_x_discrete(labels=c("Healthy","Before","After"))




pdf('XXX.pdf', width = 10, height  = 5)
ggdraw() +
    draw_plot(p4, .03,  0, .57, .95) +
    draw_plot(p5, .6, .4, .4, .5) +
    draw_plot_label(c("A", "B"), c(0, .6), c(1, 1), size = 15)
dev.off()


