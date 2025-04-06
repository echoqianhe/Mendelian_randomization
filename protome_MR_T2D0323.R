###########################2024-01-24-Protemic-MR(European Ancestry)###########################
rm(list=ls())
graphics.off()
setwd('D:/SH_desktop/Proteome_MR/pQTL-analysis') 
#install.packages("openxlsx")
library("openxlsx")
if (!require(pacman)) install.packages("pacman"); library(pacman)
p_load("xlsx")
p_load("data.table")
p_load("furrr")
p_load("dplyr")
p_load("ggplot2")
p_load("tibble")
p_load("tidyr")
p_load("purrr")
p_load("readr")
p_load("forcats")
p_load("biomaRt")
p_load_gh("MRCIEU/TwoSampleMR")
p_load_gh("MRCIEU/MRInstruments")
p_load_gh("slowkow/proxysnps")
library(gwasglue)
library(gwasvcf)
suppressWarnings(suppressPackageStartupMessages({
  library(gwasvcf)
  library(VariantAnnotation)
  library(dplyr)
  library(magrittr)
}))
library(ieugwasr)
library(phenoscanner)

###
# Perform: primary proteome - MR analysis [ARIC (EA -hg38) pQTL(proteins) on T2D] =========================================
###
##### 2024-09-23 Filter pQTL pass genome-significance

aric_path <- "ARIC/EA/EA_data/EA/"
aric_data <- ".PHENO1.glm.linear"
protein_names <- list.files(path=aric_path,pattern = paste0("*",aric_data))

tmp_aric_pqtl_exp <- NULL

for (i in protein_names){
  tmp <- suppressWarnings(fread(paste0(aric_path,i),
                                stringsAsFactors = F,
                                data.table = F))
  tmp$P <- as.numeric(tmp$P)
  print(paste0('step 1:', i))
  tmp <- tmp[tmp$P <5e-08,]
  
  if (nrow(tmp) != 0) {
    tmp$tissue <- gsub(aric_data,"",i)
    tmp_aric_pqtl_exp <- rbind(tmp_aric_pqtl_exp, tmp)
    print(paste0('step 2: combine the sig pQTL', gsub(aric_data,"",i)))
  }
}
##### 2024-09-23 Map pQTL to proteins
seq <- read.table('ARIC/EA/EA_data/seqid.txt', sep = "\t", header = T)
seq <- seq %>% dplyr::rename(chr_hg38 = chromosome_name, transcription_start_site_hg38=transcription_start_site)

aric_pqtl_exp <- tmp_aric_pqtl_exp %>% left_join(seq, by = c('tissue'='seqid_in_sample'))

##### 2024-09-23 Perform two sample MR (input: proteins; outcome: T2D)
# 1) format input data
aric_pqtl_exp <- format_data(
  aric_pqtl_exp,
  type = "exposure",
  snp_col = "ID",
  beta_col = 'BETA',
  se_col = 'SE',
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P",
  eaf_col = 'A1_FREQ',
  phenotype_col = "entrezgenesymbol"
)
write.csv(aric_pqtl_exp, file = "aric_pqtl_exp.csv", row.names = F)

aric_proteins <- unique(aric_pqtl_exp$exposure)

# 2) clump input data 
sig_aric_pqtl_exp_clump <- NULL

for (protein in aric_proteins) {
  tmp_aric_pqtl_exp_clump <- aric_pqtl_exp %>% 
    filter(exposure == protein)
  # clump exp data
  tmp_aric_pqtl_exp_clump <- clump_data(tmp_aric_pqtl_exp,
                                       clump_kb = 5000, 
                                       clump_r2 = 0.01)
  sig_aric_pqtl_exp_clump <- rbind(sig_aric_pqtl_exp_clump, tmp_aric_pqtl_exp_clump) 
}

sig_aric_pqtl_exp_clump <- read.csv('F:/SH_desktop/pQTL-analysis/sig_aric_pqtl_exp_clump0129.csv')
# 2) format outcome data
outc <- fread('D:/Kali-manuscript0921/glucose_GWAS/DIAMANTE-EUR.sumstat.txt', sep=" ")

aric_outcome_dat <- format_data(
  dat = outc,
  type = "outcome", 
  snps = sig_aric_pqtl_exp_clump$SNP,
  header = TRUE,
  phenotype_col = "phenotype",  
  snp_col = "rsID",
  beta_col = "Fixed.effects_beta",
  se_col = "Fixed.effects_SE",
  effect_allele_col = "effect_allele",                                                                                                                    
  other_allele_col = "other_allele",
  pval_col = "Fixed.effects_p.value",
  eaf_col = 'effect_allele_frequency'
)
aric_outcome_dat$outcome <- "T2D_EUR"

# 3.1) perform MR (with clump)
aric_H_data_clump <- NULL
aric_pqtl_t2d_singlesnp_results_clump <- NULL
aric_protein_name <- NULL

aric_proteins <- unique(sig_aric_pqtl_exp_clump$exposure)
# notes: 不是所有proteins在outcome_dat中都有结果，如没有，则无法添加PROTEINS
for (protein in aric_proteins) {
  tmp_sig_aric_pQTL_exp_clump <- sig_aric_pqtl_exp_clump %>% 
    filter(exposure == protein)
  
  aric_protein_name <- c(aric_protein_name, protein)
  mm <- length(aric_protein_name)
  
  print(paste0("step1: filter proteins======", aric_protein_name[mm]))
  
  tmp_aric_H_data_clump <- harmonise_data(exposure_dat = tmp_sig_aric_pQTL_exp_clump,
                                     outcome_dat = aric_outcome_dat)
  aric_H_data_clump <- rbind(aric_H_data_clump, tmp_aric_H_data_clump)
  print("step 2: Harnomize name")
  
  if (!(is.data.frame(tmp_aric_H_data_clump) && nrow(tmp_aric_H_data_clump) == 0)){
    tmp_aric_pqtl_t2d_singlesnp_results_clump <- tmp_aric_H_data_clump %>% 
      mr_singlesnp(
        parameters = default_parameters(),
        single_method = "mr_wald_ratio",
        all_method = c('mr_ivw_mre', 'mr_egger_regression',"mr_weighted_median", "mr_weighted_mode"))
    tmp_aric_pqtl_t2d_singlesnp_results_clump$proteins <- protein
  }
  print("step 3: perform MR")
  
  aric_pqtl_t2d_singlesnp_results_clump <- rbind(aric_pqtl_t2d_singlesnp_results_clump,
                                            tmp_aric_pqtl_t2d_singlesnp_results_clump)
}

write.csv(aric_H_data_clump, file = 'aric_pqtl_t2d_H_data_clump(all)_0424.csv', row.names = F)

## ARIC(EA) quality control=========================
aric_pqtl_t2d_singlesnp_results_clump$proceed[aric_pqtl_t2d_singlesnp_results_clump$exposure != aric_pqtl_t2d_singlesnp_results_clump$proteins] <- 0
aric_pqtl_t2d_singlesnp_results_clump$proceed[aric_pqtl_t2d_singlesnp_results_clump$exposure == aric_pqtl_t2d_singlesnp_results_clump$proteins] <- 1

write.csv(aric_pqtl_t2d_singlesnp_results_clump, file = 'aric_pqtl_t2d_mr_results_clump(all)_0424.csv', row.names = F)

# ARIC data selecting significant results/最终剩 1527(aric_proteins)个蛋白====
sig_aric_pqtl_t2d_singlesnp_results_clump <- aric_pqtl_t2d_singlesnp_results_clump %>% filter(
  SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
             'All - MR Egger',
             "All - Weighted median",
             "All - Weighted mode") | p < 0.05/1527)

sig_aric_pqtl_t2d_singlesnp_results_clump <- sig_aric_pqtl_t2d_singlesnp_results_clump[!is.na(sig_aric_pqtl_t2d_singlesnp_results_clump$b),]


# ARIC add nsnp to the significant results (968-sig_aric_proteins)
sig_aric_proteins <- unique(sig_aric_pqtl_t2d_singlesnp_results_clump$proteins)

for (protein in sig_aric_proteins) {
  tmp_sig_aric_pqtl_t2d <- aric_pqtl_t2d_singlesnp_results_clump %>% 
    filter(proteins == protein) %>%
    filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
             SNP != "All - MR Egger" &
             SNP != "All - Weighted median" &
             SNP != "All - Weighted mode")
  print(paste0('ARIC step1: filter snp results=====', protein))
  
  sig_aric_pqtl_t2d_singlesnp_results_clump$nSNP[sig_aric_pqtl_t2d_singlesnp_results_clump$proteins == protein]<- dim(tmp_sig_aric_pqtl_t2d)[1]
  print(paste0('ARIC step 2: add snp numbers======', dim(tmp_sig_aric_pqtl_t2d)[1]))
}
# ARIC 根据SNP列的信息，重新赋值是IVW,Wald ratio,MR_egger方法，筛选结果
sig_aric_pqtl_t2d_singlesnp_results_clump$mr_method <- ifelse(sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - MR Egger" & 
                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - Weighted median" & 
                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - Weighted mode" &
                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$nSNP == 1, "Wald_ratio", 
                                                              ifelse(sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - MR Egger" & 
                                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - Weighted median" &
                                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - Weighted mode" &
                                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$nSNP != 1, "Confused", 
                                                                              ifelse(sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - MR Egger" &
                                                                                       sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - Weighted median" &
                                                                                sig_aric_pqtl_t2d_singlesnp_results_clump$SNP != "All - Weighted mode", "IVW(random_models)", "Sensitivity analysis")))

sig_aric_pqtl_t2d_singlesnp_results_clump <- sig_aric_pqtl_t2d_singlesnp_results_clump %>% filter(sig_aric_pqtl_t2d_singlesnp_results_clump$mr_method != "Confused")

# ARIC 筛选SIG的结果 0.05/ 1563 (968)(sig_aric_proteins) - 67-sig_aric_pqtl_t2d_singlesnp_results_clump sig-results (1 duplicated proteins)
sig_aric_proteins <- unique(sig_aric_pqtl_t2d_singlesnp_results_clump$proteins)
sig_aric_pqtl_t2d_singlesnp_results_withevidence <- NULL

for (protein in sig_aric_proteins) {
  tmp_sig_aric_pqtl_t2d_singlesnp_results_clump <- sig_aric_pqtl_t2d_singlesnp_results_clump %>% filter(proteins == protein)
  tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence1 <- ifelse(tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$mr_method %in% c("IVW(random_models)", "Wald_ratio") & 
                                                                     tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05/1527, "robust",
                                                                   ifelse(tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$mr_method %in% c("IVW(random_models)", "Wald_ratio") &
                                                                            tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05, "probable","no_evidence"))
  print(paste0("Step 1:", protein))
  condition1 <- tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$SNP %in% c("All - MR Egger","All - Weighted median","All - Weighted mode") & 
    tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05
  #condition2 <- tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$SNP == "All - Weighted median" & tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05
  #condition3 <- tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$SNP == "All - Weighted mode" & tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05
  tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[condition1] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[condition2] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[condition3] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2)] <- "non_evidence"
  #
  sig_aric_pqtl_t2d_singlesnp_results_withevidence <- rbind(sig_aric_pqtl_t2d_singlesnp_results_withevidence, tmp_sig_aric_pqtl_t2d_singlesnp_results_clump)
}

sig_aric_pqtl_t2d_singlesnp_results_withevidence <- unique(sig_aric_pqtl_t2d_singlesnp_results_withevidence)

sig_aric_pqtl_t2d_singlesnp_results_withevidence_update <- NULL
for (protein in sig_aric_proteins) {
  tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence <- sig_aric_pqtl_t2d_singlesnp_results_withevidence %>% filter(proteins == protein)
  if(sum(is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence[!(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) ==3){
    tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$evidence3 <- "non_evidence"
  }
  if(sum(is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence[!(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) %in% c(1,2)){
    tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$evidence3 <- "robust"
  }
  if(sum(is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence[!(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) ==0){
    tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$evidence3 <- "no_sensitivity_test"
  }
  print(paste0("step1:", protein))
  sig_aric_pqtl_t2d_singlesnp_results_withevidence_update <- rbind(sig_aric_pqtl_t2d_singlesnp_results_withevidence_update,tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence)
}

sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence_combined <- ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence1 == "robust" & 
                                                                                      sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence3 == "robust", "robust",
                                                                                    ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence1 == "robust"&
                                                                                             sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence3 == "non_evidence","suggestive",
                                                                                           ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence1 == "robust"&
                                                                                                    sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence3 == "no_sensitivity_test" &
                                                                                                    sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$mr_method == "Wald_ratio","probable",
                                                                                           ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence1 == "probable"&
                                                                                                    sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence3 == "robust","probable",
                                                                                                  ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence1 == "probable"&
                                                                                                           sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence3 == "non_evidence","suggestive",
                                                                                                         ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence1 == "no_evidence"&
                                                                                                                  sig_aric_pqtl_t2d_singlesnp_results_withevidence_update$evidence3 != "robust","non_evidence", "suggestive"))))))
    
write.csv(sig_aric_pqtl_t2d_singlesnp_results_withevidence_update, file = "Sig_aric_pqtl_t2d_mr_results_clump_withevidence_0424.csv")

sig_aric_pqtl_t2d_singlesnp_results_clump <- sig_aric_pqtl_t2d_singlesnp_results_withevidence_update %>% 
  filter(evidence_combined %in% c("robust","probable"))

## part 1. ARIC (EA) MR results visualization====================================
##======= 1) volcano plot for discovery dataset 针对所有有MR结果的proteins绘制==========
## selecting significant results/最终剩 968 个蛋白 p1_pqtl_singlesnp_results 为volcano plot 用，与之前的sig_aric_pqtl_singlesnp_results_clump的流程一致
p1_pqtl_singlesnp_results <- aric_pqtl_t2d_singlesnp_results_clump %>% filter(
  SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
             'All - MR Egger') | p < 0.05/1527
)

p1_pqtl_singlesnp_results <- p1_pqtl_singlesnp_results[!is.na(p1_pqtl_singlesnp_results$b),]

# ARIC(EA) add nsnp to the significant results (965-p1_sig_proteins)
p1_sig_proteins <- unique(p1_pqtl_singlesnp_results$proteins)

for (protein in p1_sig_proteins) {
  tmp_sig_pqtl_t2d <- aric_pqtl_t2d_singlesnp_results_clump %>% 
    filter(proteins == protein) %>%
    filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
             SNP != "All - MR Egger")
  print(paste0('step1: filter snp results=====', protein))
  
  p1_pqtl_singlesnp_results$nSNP[p1_pqtl_singlesnp_results$proteins == protein]<- dim(tmp_sig_pqtl_t2d)[1]
  print(paste0('step 2: add snp numbers======', dim(tmp_sig_pqtl_t2d)[1]))
}
# 重新赋值方法，筛选结果p1_pqtl_singlesnp_results 包含所有proteins的MR结果
p1_pqtl_singlesnp_results$mr_method <- ifelse(p1_pqtl_singlesnp_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                p1_pqtl_singlesnp_results$SNP != "All - MR Egger" & p1_pqtl_singlesnp_results$nSNP == 1, 
                                                         "Wald_ratio", ifelse(p1_pqtl_singlesnp_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                                p1_pqtl_singlesnp_results$SNP != "All - MR Egger" & p1_pqtl_singlesnp_results$nSNP != 1, "Confused", 
                                                                              ifelse(p1_pqtl_singlesnp_results$SNP == "All - MR Egger", "MR_egger", "IVW(random_models)")))

p1_pqtl_singlesnp_results <- p1_pqtl_singlesnp_results %>% filter(p1_pqtl_singlesnp_results$mr_method != "Confused")
p1_pqtl_singlesnp_results <- p1_pqtl_singlesnp_results %>% generate_odds_ratios()

p1_pqtl_singlesnp_results <- p1_pqtl_singlesnp_results %>% filter(SNP != 'All - MR Egger')

write.csv(p1_pqtl_singlesnp_results, file = "p1_pqtl_singlesnp_results_0424.csv", row.names = F)

## volcano plot======2024-02-28 合并primary MR analysis and sesitivity MR analysis (delete horizontal pleiotropy)的结果
primary_sesitivity_analysis_consistent_list <- intersect(sig_aric_pqtl_t2d_singlesnp_results_clump$exposure, 
                                                         sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$exposure)
options(ggrepel.max.overlaps = Inf)
#CC6666
pdf("pqtl_t2d_volcona_plot(EA)_0428.pdf",
    width = 8,
    height = 8
)
p1_pqtl_singlesnp_results %>% 
  mutate(change=factor(ifelse(p < 0.05 & or !=1, ifelse(or > 1 ,'Up','Down'),'NoSignificant'),levels=c('Up','Down','NoSignificant'))) %>%
  ggplot(aes(x=or, y=-log10(p), color=change))+geom_point()+
  scale_color_manual(values=c("#9b3a74","#156077","#808080"))+
  ggrepel::geom_label_repel(aes(label=ifelse(exposure %in% primary_sesitivity_analysis_consistent_list, exposure,'')))+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'none'
  )+
  ylab('-log10 (Pval)')+xlab('OR of T2D')+
  ggtitle('Volcano Plot of MR analysis in European Ancestry')+
  geom_vline(xintercept=1,lty=3,col="black",lwd=1)  +
  geom_hline(yintercept = -log10(0.05/1563),lty=3,col="black",lwd=1) 
dev.off()

##===== 2) MR forest plot visualization针对48 consistent个SIG结果绘制 ===========================
# forest plot (EA)======= 2024-02-28 合并primary MR analysis and sesitivity MR analysis (delete horizontal pleiotropy)的结果
p2_sig_pqtl_singlesnp_results <- sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy %>% 
  filter(exposure %in% primary_sesitivity_analysis_consistent_list) %>%
  generate_odds_ratios()


p2_sig_pqtl_singlesnp_results$orlab <- paste0("",
                           ifelse(sprintf("%.2f",p2_sig_pqtl_singlesnp_results$or)<0.0051,
                                  format(p2_sig_pqtl_singlesnp_results$or,scientific = TRUE,digits=3),
                                  sprintf("%.2f",p2_sig_pqtl_singlesnp_results$or)),
                           "(", 
                           ifelse(sprintf("%.2f",p2_sig_pqtl_singlesnp_results$or_lci95)<0.0051,
                                  format(p2_sig_pqtl_singlesnp_results$or_lci95,scientific = TRUE,digits=3),
                                  sprintf("%.2f",p2_sig_pqtl_singlesnp_results$or_lci95)),
                           "-",
                           ifelse(sprintf("%.2f",p2_sig_pqtl_singlesnp_results$or_uci95)<0.0051,
                                  format(p2_sig_pqtl_singlesnp_results$or_uci95,scientific = TRUE,digits=3),
                                  sprintf("%.2f",p2_sig_pqtl_singlesnp_results$or_uci95)),
                           ")")


p2_sig_pqtl_singlesnp_results <- p2_sig_pqtl_singlesnp_results[order(p2_sig_pqtl_singlesnp_results[,'or']),]
p2_sig_pqtl_singlesnp_results$exposure <- factor(p2_sig_pqtl_singlesnp_results$exposure, levels = p2_sig_pqtl_singlesnp_results$exposure)

p2_sig_pqtl_singlesnp_results$x <- 5.7
p2_sig_pqtl_singlesnp_results$y <- seq(1,48,1)

# annotation 
annotation<- NULL
annotation0 <- p2_sig_pqtl_singlesnp_results[c('x','y','nSNP')] 
annotation1 <- p2_sig_pqtl_singlesnp_results[c('x','y','orlab')] %>% dplyr::rename(nSNP=orlab) %>% rbind(annotation0)
annotation2 <- p2_sig_pqtl_singlesnp_results[c('x','y','exposure')] %>% dplyr::rename(nSNP=exposure) 
annotation <- annotation1 %>% rbind(annotation2) 

annotation[1:48,1] <- 0.3
annotation[49:96,1] <- -0.1
annotation[97:144,1] <- -0.6


annotation <- rbind(c(0.3,48.4,"OR(95%CI)"),annotation)
annotation <- rbind(c(-0.1,48.4,"nSNP"),annotation)
annotation <- rbind(c(-0.6,48.4,"Proteins"),annotation)

annotation$x <- as.numeric(annotation$x)
annotation$y <- as.numeric(annotation$y)
# forest plot====
pdf("pqtl_t2d_forest_plot(EA)_0428.pdf",
    width = 8,
    height = 8
)
p2_sig_pqtl_singlesnp_results %>% 
  ggplot(aes(or, exposure)) +
  geom_point(aes(col=or,size = nSNP)) +
  scale_color_continuous_sequential(palette = "ag_Sunset", l1 = 20, c2 = 70, p1 = 1) +
  #scale_color_continuous_diverging(palette = "Tofino") +
  geom_errorbarh(aes(xmax=or_lci95, xmin=or_uci95, col = or), height=0.2, size=1.0) +
  scale_x_continuous(limits = c(-0.8,1.5), breaks = seq(0,2.0,0.5)) +
  geom_vline(aes(xintercept = 1)) + theme_bw()+
  theme(axis.text.y = element_blank(), axis.ticks.y =element_blank()) +
  theme(axis.line.x = element_line(color = "#808080")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), plot.background = element_rect(color="black",size = 0.5)) +
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6)) +
  geom_text(data = annotation, aes(x=x, y=y, label=nSNP), size=2.8) +
 #labs(colour = "P value") +
  #coord_fixed(ratio = 1.2) +
  xlab("Odds Ratio") + ylab(" ") +
  ggtitle(expression("MR: Proteins"  %->% "T2D(European Ancestry)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###
# Perform: coloc-lacalization [ARIC(EA-hg38) sig MR proteins & T2D] 针对67个SIG proteins===================
###
##### 2024-09-23 perform pqtl - proteins and T2D colocalization analysis
library("coloc")

aric_path <- "ARIC/EA/EA_data/EA/"
aric_data <- ".PHENO1.glm.linear"
protein_names <- list.files(path=aric_path,pattern = paste0("*",aric_data))
# only select 67 sig protein gwas files
#seq <- read.table('C:/Users/Administrator/Downloads/seqid.txt', sep = "\t", header = T)
seq <- seq %>% dplyr::rename(chr_hg38 = chromosome_name, transcription_start_site_hg38=transcription_start_site)

# 2024-04-22 filter robust and probable proteins (72 = 65+7)
tmp_aric_pqtl_exp <- sig_aric_pqtl_t2d_singlesnp_results_clump %>% dplyr::select(exposure,proteins,evidence_combined)
# 2024-04-22 合并后发现有少数protein有多个探针(CD14,HDGF,HP,IL15RA,LRP11,MSR1,SVEP1)，根据质控原则，删除重复探针的protein
# delete 7 dubplicated proteins
tmp_aric_pqtl_exp <- tmp_aric_pqtl_exp %>% filter(!(exposure %in% c("CD14","D14","HDGF","HP","IL15RA","LRP11","MSR1","SVEP1")))
tmp_aric_pqtl_exp <- tmp_aric_pqtl_exp %>% left_join(seq, by = c('exposure'='entrezgenesymbol'))

# 2024-04-22 only read 54 proteins GWAS in the following step/tmp_aric_pqtl_exp_proteins is the vector include 61 proteins names
tmp_aric_pqtl_exp_proteins <- tmp_aric_pqtl_exp$seqid_in_sample
tmp_aric_pqtl_exp_proteins <- paste(tmp_aric_pqtl_exp_proteins,".PHENO1.glm.linear", sep = "")

## COLOC: this part code is for coloc analysis
tmp_aric_pqtl_exp <- NULL

for (i in tmp_aric_pqtl_exp_proteins){
  tmp <- suppressWarnings(fread(paste0(aric_path,i),
                                stringsAsFactors = F,
                                data.table = F))
  tmp$P <- as.numeric(tmp$P)
  print(paste0('step 1:', i))
  tmp <- tmp[tmp$P <1e-04,]
  
  if (nrow(tmp) != 0) {
    tmp$tissue <- gsub(aric_data,"",i)
    tmp_aric_pqtl_exp <- rbind(tmp_aric_pqtl_exp, tmp)
    print(paste0('step 2: combine the sig pQTL', gsub(aric_data,"",i)))
  }
}

# select significant original file of 67 proteins (hg38) 
# sig_tmp_aric_pqtl_exp 为包含61个SIG protein的PQTL(P<5e-04) 可用于做COLOC分析
tmp_aric_pqtl_exp <- tmp_aric_pqtl_exp %>% left_join(seq, by = c('tissue'='seqid_in_sample'))
sig_tmp_aric_pqtl_exp <- tmp_aric_pqtl_exp

sig_tmp_aric_pqtl_exp <- sig_tmp_aric_pqtl_exp %>% dplyr::rename(CHR="#CHROM",SNP=ID, BP=POS)

sig_tmp_aric_pqtl_exp <- sig_tmp_aric_pqtl_exp %>% dplyr::select(SNP,CHR,BP,REF,ALT,A1,A1_FREQ,BETA,SE,P,tissue,uniprot_id,entrezgenesymbol)


# convert hg38 to hg37
## liftover sumstats_dt为HG37坐标系的ARIC中的SIG(MR-67) PQTL 
# sumstats_dt 为HG37坐标系的 sig_tmp_aric_pqtl_exp
sumstats_dt <- MungeSumstats::liftover(sumstats_dt = sig_tmp_aric_pqtl_exp,
                                                 ref_genome = "hg38",
                                                 convert_ref_genome = "hg19")
knitr::kable(head(sumstats_dt))

# 合并pqtl和GWAS的数据 [sig_tmp_aric_pqtl_exp_outc] 为[protine - T2D]进行COLOC分析的输入数据
sig_tmp_aric_pqtl_exp_outc <- sig_tmp_aric_pqtl_exp %>% left_join(outc, by = c("SNP"="rsID"))
write.csv(sig_tmp_aric_pqtl_exp_outc, file = "sig_tmp_aric_pqtl_exp_outc_0424.csv", row.names = F)

# remove duplicated SNPs
sig_tmp_aric_pqtl_exp_outc <- sig_tmp_aric_pqtl_exp_outc[order(sig_tmp_aric_pqtl_exp_outc[,'SNP'], sig_tmp_aric_pqtl_exp_outc[,'P']),]
sig_tmp_aric_pqtl_exp_outc <- sig_tmp_aric_pqtl_exp_outc[!duplicated(sig_tmp_aric_pqtl_exp_outc$SNP), ]

# run coloc analysis based on sig_tmp_aric_pqtl_exp_outc
coloc_aric_proteins <- unique(sig_tmp_aric_pqtl_exp_outc$entrezgenesymbol)

coloc_aric_pqtl_t2d <- NULL

for (i in coloc_aric_proteins){
  # dataset 1: pqtl summary; dataset 2: t2d
  tmp_coloc_aric_pqtl_t2d <- sig_tmp_aric_pqtl_exp_outc %>% filter(entrezgenesymbol == i)
  
  tmp_coloc_aric_pqtl_t2d <- tmp_coloc_aric_pqtl_t2d %>% na.omit()
  tmp_coloc_aric_pqtl_t2d <- tmp_coloc_aric_pqtl_t2d[tmp_coloc_aric_pqtl_t2d$P != 0,]
  
  tmp <- coloc.abf(dataset2 = list(pvalues=tmp_coloc_aric_pqtl_t2d$Fixed.effects_p.value, 
                                   type = "cc", s = 0.16 , 
                                   N = 429191, snp=tmp_coloc_aric_pqtl_t2d$SNP),
                   dataset1 = list(pvalues=tmp_coloc_aric_pqtl_t2d$P, type ="quant", 
                                   N = 7213, snp = tmp_coloc_aric_pqtl_t2d$SNP), 
                   MAF = tmp_coloc_aric_pqtl_t2d$A1_FREQ)
  print(paste0('step 1: coloc-analysis=====',i))
  tmp <- tmp$results
  if (nrow(tmp) != 0) {
    tmp$protein <- i
    coloc_aric_pqtl_t2d <- rbind(coloc_aric_pqtl_t2d, tmp)
  }
}


write.csv(coloc_aric_pqtl_t2d, file = "coloc_aric_pqtl_t2d_0424.csv", row.names = F)

sig_coloc_aric_pqtl_t2d <- coloc_aric_pqtl_t2d %>% dplyr::filter(SNP.PP.H4>0.8)

# 查看coloc pph4>0.8和sig mr result的重叠
intersect(sig_coloc_aric_pqtl_t2d$protein, sig_aric_pqtl_t2d_singlesnp_results_clump$exposure)
# 查看coloc pph4>0.8和7个robust 的mr result的重叠
intersect(sig_coloc_aric_pqtl_t2d$protein,sig_aric_pqtl_t2d_singlesnp_results_clump_7robust)

# save file [MR info & coloc info] 29 proteins pass the PP.H4 test 其中 [10 proteins shared the same pqtl(IVs) with MR analysis]
# 10/29 proteins share the same pqtl(IVs) with MR anslysis
# sig_aric_pqtl_t2d_ivs_coloc 为 pass PP.H4 & pass unique pqtl(IVs) check 10 proteins(采用p<1e-04)
write.csv(sig_coloc_aric_pqtl_t2d, file = 'sig_coloc_aric_pqtl_t2d_pph40.8_0424.csv', row.names = F)
# sig_aric_pqtl_t2d_ivs_coloc 为 PASS PP.H4 及 PQTL 在 MR analysis 中的
sig_aric_pqtl_t2d_ivs_coloc <- sig_coloc_aric_pqtl_t2d %>% filter(snp %in% aric_pqtl_t2d_singlesnp_results_clump$SNP)

###
# Perform: check association between pqtl & eqtl(EA eqtlGen & GTEx v8)==========================
###
##### 2024-09-23 Perform pleiotropy test

# 1)select IVs(pQTL) information from aric_pqtl_t2d_singlesnp_results_clump data
sig_aric_pqtl_t2d_ivs <- aric_pqtl_t2d_singlesnp_results_clump %>% 
  filter(aric_pqtl_t2d_singlesnp_results_clump$exposure %in% sig_aric_pqtl_t2d_singlesnp_results_clump$exposure) %>%
  filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & SNP != "All - MR Egger")
# rs4680458 & rs16829322 (perfect proxy for RARRES1)RARRES1 7个robust gene之一
# 
# 2)check if coloc snp was colocalized with IVs(pQTL) 14 + 1(LD Proxy) = 15/proteins
intersect(sig_aric_pqtl_t2d_ivs$SNP,sig_coloc_aric_pqtl_t2d$snp)

###
# horizontal pleiotropy test: check pQTL pleiotropy in 61-7 = 54 proteins (robust & probable MR results) tmp_aric_pqtl_exp_proteins
###
# 1)step-1: check if cis-pqtl associated with other proteins or genes
sig_aric_pqtl_t2d_ivs_qc <- sig_aric_pqtl_t2d_ivs %>% 
  filter(!(exposure %in% c("CD14","D14","HDGF","HP","IL15RA","LRP11","MSR1","SVEP1"))) 
sig_aric_pqtl_exp_clump %>% filter(SNP %in% sig_aric_pqtl_t2d_ivs_qc$SNP & !(exposure %in% sig_aric_pqtl_t2d_ivs_qc$exposure))
# 2024-04-22 identify 2 pQTL may have pleiotropy rs78125821(ERAP2)/rs1050541(SAT2) 且该2个proteins与t2d mr的结果为p>0.05

# 1.1)check if cis-pqtl associated with other genes 
#### read eQTL dataset 
eqtl <- fread("/Volumes/Echo_SSD/eQTL_summary_data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
              sep="\t", header = T)
eqtl <- eqtl %>% filter(SNP %in% sig_aric_pqtl_t2d_ivs_qc$SNP & !(GeneSymbol %in% sig_aric_pqtl_t2d_ivs_qc$exposure))
eqtl <- eqtl %>% filter(Pvalue < 5e-08 & FDR < 0.05)

# 2024-04-22 identify 165 pQTL may have pleiotropy 
# 查看165pqtl对应的基因 有多少测了protein abundence
eqtl_seq <- eqtl %>% left_join(seq, by = c("GeneSymbol" ="entrezgenesymbol"))
eqtl_seq <- eqtl_seq[!is.na(eqtl_seq$seqid_in_sample)]
eqtl_seq <- eqtl_seq[order(eqtl_seq$SNP,eqtl_seq$Pvalue),]
eqtl_seq <- eqtl_seq[!duplicated(eqtl_seq$SNP),]
eqtl_seq <- eqtl_seq[order(eqtl_seq$seqid_in_sample,eqtl_seq$Pvalue),]
# identify the potential pleiotropy pqtl from ARIC
tmp_eqtl_seq_proteins <- eqtl_seq$seqid_in_sample
tmp_eqtl_seq_proteins <- paste(tmp_eqtl_seq_proteins,".PHENO1.glm.linear", sep = "")

# 2)step-2:identify the secondary proteins/genes
tmp_eqtl_seq <- NULL

for (i in tmp_eqtl_seq_proteins){
  tmp <- suppressWarnings(fread(paste0(aric_path,i),
                                stringsAsFactors = F,
                                data.table = F))
  #tmp$P <- as.numeric(tmp$P)
  print(paste0('step 1:', i))
  tmp <- tmp[tmp$ID %in% eqtl_seq$SNP,]
  
  if (nrow(tmp) != 0) {
    tmp$tissue <- gsub(aric_data,"",i)
    tmp_eqtl_seq <- rbind(tmp_eqtl_seq, tmp)
    print(paste0('step 2: combine the sig pQTL', gsub(aric_data,"",i)))
  }
}
tmp_eqtl_seq <- tmp_eqtl_seq %>% left_join(seq, by = c('tissue'='seqid_in_sample'))

# 3)step-3:perform two-sample MR analysis for secondary proteins and T2D
# 3.1) format input data
secondary_pqtl_eqtl_seq_exp <- format_data(
  tmp_eqtl_seq,
  type = "exposure",
  snp_col = "ID",
  beta_col = 'BETA',
  se_col = 'SE',
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P",
  eaf_col = 'A1_FREQ',
  phenotype_col = "entrezgenesymbol"
)
write.csv(secondary_pqtl_eqtl_seq_exp, file = "secondary_pqtl_eqtl_seq_exp0423.csv", row.names = F)

secondary_pqtl_eqtl_proteins <- unique(secondary_pqtl_eqtl_seq_exp$exposure)

# 3.2) clump input data 
sig_secondary_pqtl_eqtl_exp_clump <- NULL

for (protein in secondary_pqtl_eqtl_proteins) {
  tmp_secondary_pqtl_eqtl_seq_exp_clump <- secondary_pqtl_eqtl_seq_exp %>% 
    filter(exposure == protein)
  # clump exp data
  tmp_secondary_pqtl_eqtl_seq_exp_clump <- clump_data(tmp_secondary_pqtl_eqtl_seq_exp_clump,
                                        clump_kb = 5000, 
                                        clump_r2 = 0.01)
  sig_secondary_pqtl_eqtl_exp_clump <- rbind(sig_secondary_pqtl_eqtl_exp_clump, tmp_secondary_pqtl_eqtl_seq_exp_clump) 
}

write.csv(sig_secondary_pqtl_eqtl_exp_clump,'sig_secondary_pqtl_eqtl_exp_clump0423.csv', row.names = F)
# 3.3) format outcome data
outc <- fread('/Volumes/Echo_SSD/Kali-manuscript0921/glucose_GWAS/DIAMANTE-EUR.sumstat.txt', sep=" ")

secondary_outcome_dat <- format_data(
  dat = outc,
  type = "outcome", 
  snps = sig_secondary_pqtl_eqtl_exp_clump$SNP,
  header = TRUE,
  phenotype_col = "phenotype",  
  snp_col = "rsID",
  beta_col = "Fixed.effects_beta",
  se_col = "Fixed.effects_SE",
  effect_allele_col = "effect_allele",                                                                                                                    
  other_allele_col = "other_allele",
  pval_col = "Fixed.effects_p.value",
  eaf_col = 'effect_allele_frequency'
)
secondary_outcome_dat$outcome <- "T2D_EUR"

# 3.4) perform MR (with clump)
secondary_H_data_clump <- NULL
secondary_pqtl_t2d_singlesnp_results_clump <- NULL
secondary_protein_name <- NULL

secondary_proteins <- unique(sig_secondary_pqtl_eqtl_exp_clump$exposure)
# notes: 不是所有proteins在outcome_dat中都有结果，如没有，则无法添加PROTEINS
for (protein in secondary_proteins) {
  tmp_sig_secondary_pqtl_eqtl_exp_clump <- sig_secondary_pqtl_eqtl_exp_clump %>% 
    filter(exposure == protein)
  
  secondary_protein_name <- c(secondary_protein_name, protein)
  mm <- length(secondary_protein_name)
  
  print(paste0("step1: filter proteins======", secondary_protein_name[mm]))
  
  tmp_secondary_H_data_clump <- harmonise_data(exposure_dat = tmp_sig_secondary_pqtl_eqtl_exp_clump,
                                          outcome_dat = secondary_outcome_dat)
  secondary_H_data_clump <- rbind(secondary_H_data_clump, tmp_secondary_H_data_clump)
  print("step 2: Harnomize name")
  
  if (!(is.data.frame(tmp_secondary_H_data_clump) && nrow(tmp_secondary_H_data_clump) == 0)){
    tmp_secondary_pqtl_t2d_singlesnp_results_clump <- tmp_secondary_H_data_clump %>% 
      mr_singlesnp(
        parameters = default_parameters(),
        single_method = "mr_wald_ratio",
        all_method = c('mr_ivw_mre', 'mr_egger_regression',"mr_weighted_median", "mr_weighted_mode"))
    tmp_secondary_pqtl_t2d_singlesnp_results_clump$proteins <- protein
  }
  print("step 3: perform MR")
  
  secondary_pqtl_t2d_singlesnp_results_clump <- rbind(secondary_pqtl_t2d_singlesnp_results_clump,
                                                          tmp_secondary_pqtl_t2d_singlesnp_results_clump)
}

write.csv(secondary_H_data_clump, file = 'secondary_pqtl_eqtl_t2d_H_data_clump(all)0423.csv', row.names = F)

## ARIC(EA) quality control=========================
secondary_pqtl_t2d_singlesnp_results_clump$proceed[secondary_pqtl_t2d_singlesnp_results_clump$exposure != secondary_pqtl_t2d_singlesnp_results_clump$proteins] <- 0
secondary_pqtl_t2d_singlesnp_results_clump$proceed[secondary_pqtl_t2d_singlesnp_results_clump$exposure == secondary_pqtl_t2d_singlesnp_results_clump$proteins] <- 1

# 合并secondary pqtl & eqtl mr的结果 pqtl识别2个secondary proteins - rs78125821(ERAP2)/rs1050541(SAT2) 
secondary_pqtl_t2d_singlesnp_results_clump_pqtl <- aric_pqtl_t2d_singlesnp_results_clump %>% filter(exposure %in% c('ERAP2','SAT2'))
# eqtl 识别 27个 secondary proteins 
# mr分析 识别2 个secondary proteins: horizontal pleiotropy proteins (ADK(rs10824064/rs74656412),TNFSF12(rs1050541/rs11651783))
secondary_pqtl_t2d_singlesnp_results_clump <- rbind(secondary_pqtl_t2d_singlesnp_results_clump,secondary_pqtl_t2d_singlesnp_results_clump_pqtl)

# 查看pleiotropy pqtl的primary proteins[PLAU(rs10824064/rs74656412)-coloc.PPH4>0.8,SHBG(rs1050541/rs11651783)]
sig_aric_pqtl_t2d_ivs_qc %>% filter(SNP %in% c("rs11651783","rs74656412","rs1050541","rs10824064"))

write.csv(secondary_pqtl_t2d_singlesnp_results_clump, file = 'secondary_pqtl_eqtl_t2d_singlesnp_results_clump(all)0423.csv')

# step-4: check if PLAU & ADK - SHBG & TNFSF12 are not in the same pathway(GPS-Prot server/Enrichr)-vertical pleiotropy or horizontal pleiotropy
# none of them in the same pathway

# step-5: delete the horizontal pQTL and run the MR again
sig_aric_pqtl_exp_clump <- read.csv('sig_aric_pqtl_exp_clump0129.csv')
secondary_pqtl_eqtl_snp_list <- secondary_pqtl_t2d_singlesnp_results_clump %>% filter(!(SNP %in% c("All - Inverse variance weighted (multiplicative random effects)",
                                                                                                   "All - MR Egger","All - Weighted median",
                                                                                                   "All - Weighted mode")))
secondary_pqtl_eqtl_snp_list <- unique(secondary_pqtl_eqtl_snp_list$SNP)
sig_aric_pqtl_exp_clump_nopleiotropy <- sig_aric_pqtl_exp_clump %>% filter(!(SNP %in% secondary_pqtl_eqtl_snp_list))
sig_aric_pqtl_exp_nopleiotropy <- aric_pqtl_exp %>% filter(!(SNP %in% secondary_pqtl_eqtl_snp_list))
# 2) format outcome data
outc <- fread('/Volumes/Echo_SSD/Kali-manuscript0921/glucose_GWAS/DIAMANTE-EUR.sumstat.txt', sep=" ")

aric_outcome_dat_nopleiotropy <- format_data(
  dat = outc,
  type = "outcome", 
  snps = sig_aric_pqtl_exp_clump_nopleiotropy$SNP,
  header = TRUE,
  phenotype_col = "phenotype",  
  snp_col = "rsID",
  beta_col = "Fixed.effects_beta",
  se_col = "Fixed.effects_SE",
  effect_allele_col = "effect_allele",                                                                                                                    
  other_allele_col = "other_allele",
  pval_col = "Fixed.effects_p.value",
  eaf_col = 'effect_allele_frequency'
)
aric_outcome_dat_nopleiotropy$outcome <- "T2D_EUR"

# 3.1) perform MR (with clump)
aric_H_data_clump_nopleiotropy <- NULL
aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy <- NULL
aric_protein_name_nopleiotropy <- NULL

aric_proteins_nopleiotropy <- unique(sig_aric_pqtl_exp_clump_nopleiotropy$exposure)

# notes: 不是所有proteins在outcome_dat中都有结果，如没有，则无法添加PROTEINS
for (protein in aric_proteins_nopleiotropy) {
  tmp_sig_aric_pQTL_exp_clump_nopleiotropy <- sig_aric_pqtl_exp_clump_nopleiotropy %>% 
    filter(exposure == protein)
  
  aric_protein_name_nopleiotropy <- c(aric_protein_name_nopleiotropy, protein)
  mm <- length(aric_protein_name_nopleiotropy)
  
  print(paste0("step1: filter proteins======", aric_protein_name_nopleiotropy[mm]))
  
  tmp_aric_H_data_clump_nopleiotropy <- harmonise_data(exposure_dat = tmp_sig_aric_pQTL_exp_clump_nopleiotropy,
                                          outcome_dat = aric_outcome_dat_nopleiotropy)
  aric_H_data_clump_nopleiotropy <- rbind(aric_H_data_clump_nopleiotropy, tmp_aric_H_data_clump_nopleiotropy)
  print("step 2: Harnomize name")
  
  if (!(is.data.frame(tmp_aric_H_data_clump_nopleiotropy) && nrow(tmp_aric_H_data_clump_nopleiotropy) == 0)){
    tmp_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy <- tmp_aric_H_data_clump_nopleiotropy %>% 
      mr_singlesnp(
        parameters = default_parameters(),
        single_method = "mr_wald_ratio",
        all_method = c('mr_ivw_mre', 'mr_egger_regression',"mr_weighted_median", "mr_weighted_mode"))
    tmp_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$proteins <- protein
  }
  print("step 3: perform MR")
  
  aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy <- rbind(aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy,
                                                 tmp_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy)
}


write.csv(aric_H_data_clump_nopleiotropy, file = 'aric_pqtl_t2d_H_data_clump_nopleiotropy(all)0423.csv', row.names = F)

## ARIC(EA) quality control=========================
aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$proceed[aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$exposure != aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$proteins] <- 0
aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$proceed[aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$exposure == aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$proteins] <- 1

write.csv(aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy, file = 'aric_pqtl_t2d_mr_results_clump_nopleiotropy(all)0421.csv', row.names = F)

# ARIC data selecting significant results/最终剩 1563(aric_proteins)个蛋白====
sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy <- aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy %>% filter(
  SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
             'All - MR Egger',
             "All - Weighted median",
             "All - Weighted mode") | p < 0.05/1563)

sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy <- sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy[!is.na(sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$b),]


# ARIC add nsnp to the significant results (968-sig_aric_proteins)
sig_aric_proteins <- unique(sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$proteins)

for (protein in sig_aric_proteins) {
  tmp_sig_aric_pqtl_t2d <- aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy %>% 
    filter(proteins == protein) %>%
    filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
             SNP != "All - MR Egger" &
             SNP != "All - Weighted median" &
             SNP != "All - Weighted mode")
  print(paste0('ARIC step1: filter snp results=====', protein))
  
  sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$nSNP[sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$proteins == protein]<- dim(tmp_sig_aric_pqtl_t2d)[1]
  print(paste0('ARIC step 2: add snp numbers======', dim(tmp_sig_aric_pqtl_t2d)[1]))
}
# ARIC 根据SNP列的信息，重新赋值是IVW,Wald ratio,MR_egger方法，筛选结果
sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$mr_method <- ifelse(sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                             sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - MR Egger" & 
                                                                             sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - Weighted median" & 
                                                                             sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - Weighted mode" &
                                                                             sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$nSNP == 1, "Wald_ratio", 
                                                              ifelse(sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                       sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - MR Egger" & 
                                                                       sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - Weighted median" &
                                                                       sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - Weighted mode" &
                                                                       sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$nSNP != 1, "Confused", 
                                                                     ifelse(sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - MR Egger" &
                                                                              sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - Weighted median" &
                                                                              sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP != "All - Weighted mode", "IVW(random_models)", "Sensitivity analysis")))

sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy <- sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy %>% 
  filter(sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$mr_method != "Confused")

# ARIC 筛选SIG的结果 0.05/ 1563 (968)(sig_aric_proteins) - 67-sig_aric_pqtl_t2d_singlesnp_results_clump sig-results (1 duplicated proteins)
sig_aric_proteins <- unique(sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$proteins)
sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy <- NULL

for (protein in sig_aric_proteins) {
  tmp_sig_aric_pqtl_t2d_singlesnp_results_clump <- sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy %>% filter(proteins == protein)
  tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence1 <- ifelse(tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$mr_method == "IVW(random_models)" & 
                                                                      tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05/1563, "robust",
                                                                    ifelse(tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$mr_method == "IVW(random_models)" &
                                                                             tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05, "probable","no_evidence"))
  print(paste0("Step 1:", protein))
  condition1 <- tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$SNP %in% c("All - MR Egger","All - Weighted median","All - Weighted mode") & 
    tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05
  #condition2 <- tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$SNP == "All - Weighted median" & tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05
  #condition3 <- tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$SNP == "All - Weighted mode" & tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05
  tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[condition1] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[condition2] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[condition3] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2)] <- "non_evidence"
  #
  sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy <- rbind(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy, tmp_sig_aric_pqtl_t2d_singlesnp_results_clump)
}

sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy <- unique(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy)

sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update <- NULL
for (protein in sig_aric_proteins) {
  tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence <- sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy %>% filter(proteins == protein)
  if(sum(is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence[tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$mr_method != "IVW(random_models)",]$evidence2)) ==3){
    tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$evidence3 <- "non_evidence"
  }
  if(sum(is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence[tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$mr_method != "IVW(random_models)",]$evidence2)) %in% c(1,2)){
    tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$evidence3 <- "robust"
  }
  if(sum(is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence[tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$mr_method != "IVW(random_models)",]$evidence2)) ==0){
    tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence$evidence3 <- "no_sensitivity_test"
  }
  print(paste0("step1:", protein))
  sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update <- rbind(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update,tmp_sig_aric_pqtl_t2d_singlesnp_results_withevidence)
}

sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence_combined <- ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence1 == "robust" & 
                                                                                                   sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence3 == "robust", "robust",
                                                                                    ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence1 == "robust"&
                                                                                             sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence3 == "non_evidence","suggestive",
                                                                                           ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence1 == "probable"&
                                                                                                    sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence3 == "robust","probable",
                                                                                                  ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence1 == "probable"&
                                                                                                           sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence3 == "non_evidence","suggestive",
                                                                                                         ifelse(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence1 == "no_evidence"&
                                                                                                                  sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update$evidence3 != "robust","non_evidence", "suggestive")))))

write.csv(sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update, file = "Sig_aric_pqtl_t2d_mr_results_clump_withevidence_nopleiotropy_0424.csv")

sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy <- sig_aric_pqtl_t2d_singlesnp_results_withevidence_nopleiotropy_update %>% 
  filter(evidence_combined %in% c("robust","probable"))

# 查看coloc pph4>0.8和sig mr result的重叠-34 proteins
intersect(sig_coloc_aric_pqtl_t2d$protein, sig_aric_pqtl_t2d_singlesnp_results_clump$exposure)

# 查看coloc pph4>0.8和7个robust 的mr result的重叠 7 robust proteins
intersect(sig_coloc_aric_pqtl_t2d$protein,sig_aric_pqtl_t2d_singlesnp_results_clump_7robust)

# 查看coloc pph4>0.8和sig mr result的重叠(20个重叠) 
intersect(sig_coloc_aric_pqtl_t2d$protein, sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$exposure)

# 下一步：查看 PASS PPH4 的SNP 中有多少在 MR ANALYSIS IVS 
# sig_aric_pqtl_t2d_ivs_coloc 为 PASS PP.H4 及 PQTL 在 MR analysis 中的 (9个蛋白质pass pph4 test & mr analysis)
sig_aric_pqtl_t2d_ivs_coloc <- sig_coloc_aric_pqtl_t2d %>% 
  filter(protein %in% intersect(sig_coloc_aric_pqtl_t2d$protein, sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$exposure)) %>%
  filter(snp %in% aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$SNP)

# 查看coloc pph4>0.8和7个robust 的mr result的重叠
intersect(sig_coloc_aric_pqtl_t2d$protein,sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy_robust)

consistent_proteins_robust_and_probable <- intersect(sig_aric_pqtl_t2d_singlesnp_results_clump$exposure,sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$exposure)
consistent_proteins_robust <- intersect(sig_aric_pqtl_t2d_singlesnp_results_clump_7robust,sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy_robust)

sig_aric_pqtl_t2d_singlesnp_results_clump %>% filter(exposure %in% consistent_proteins_robust_and_probable)

intersect(sig_coloc_aric_pqtl_t2d$protein,consistent_proteins_robust_and_probable)

## LOCUSZOOM 2024-04-24 this part code is for locuszoom plot (only select sig 20 proteins)
tmp_aric_pqtl_exp_locuszoom <- NULL

for (i in protein_names){
  tmp <- suppressWarnings(fread(paste0(aric_path,i),
                                stringsAsFactors = F,
                                data.table = F))
  tmp$P <- as.numeric(tmp$P)
  print(paste0('step 1:', i))
  tmp <- tmp[tmp$P <5e-01,]
  
  if (nrow(tmp) != 0) {
    tmp$tissue <- gsub(aric_data,"",i)
    tmp_aric_pqtl_exp_locuszoom <- rbind(tmp_aric_pqtl_exp_locuszoom, tmp)
    print(paste0('step 2: combine the sig pQTL', gsub(aric_data,"",i)))
  }
}
#seq <- read.table('C:/Users/Administrator/Downloads/seqid.txt', sep = "\t", header = T)
seq <- seq %>% dplyr::rename(chr_hg38 = chromosome_name, transcription_start_site_hg38=transcription_start_site)

tmp_aric_pqtl_exp_locuszoom <- tmp_aric_pqtl_exp_locuszoom %>% left_join(seq, by = c('tissue'='seqid_in_sample'))

# select significant original file of 67 proteins (hg38) 
# sig_tmp_aric_pqtl_exp_locuszoom 为包含67个SIG protein的PQTL(P<5e-01) 可用于做locuszoom plot分析
sig_tmp_aric_pqtl_exp_locuszoom <- tmp_aric_pqtl_exp_locuszoom %>% filter(entrezgenesymbol %in% sig_aric_pqtl_t2d_singlesnp_results_clump$exposure)
sig_tmp_aric_pqtl_exp_locuszoom <- sig_tmp_aric_pqtl_exp_locuszoom %>% dplyr::rename(CHR="#CHROM",SNP=ID, BP=POS)

sig_tmp_aric_pqtl_exp_locuszoom <- sig_tmp_aric_pqtl_exp_locuszoom %>% dplyr::select(SNP,CHR,BP,REF,ALT,A1,A1_FREQ,BETA,SE,P,tissue,uniprot_id,entrezgenesymbol)


# convert hg38 to hg37
## liftover sumstats_dt为HG37坐标系的ARIC中的SIG(MR-67) PQTL 
# sumstats_dt 为HG37坐标系的 sig_tmp_aric_pqtl_exp
sumstats_dt_locuszoom <- MungeSumstats::liftover(sumstats_dt = sig_tmp_aric_pqtl_exp_locuszoom,
                                        ref_genome = "hg38",
                                        convert_ref_genome = "hg19")
knitr::kable(head(sumstats_dt_locuszoom))

# 合并pqtl和GWAS的数据 [sig_tmp_aric_pqtl_exp_outc] 为进行COLOC分析的输入数据
# 2024-04-20 sig_tmp_aric_pqtl_exp_outc与sig_tmp_aric_pqtl_exp_locuszoom_outc的区别为p阈值，绘制locuszoom的数据集纳入p<0.1的variants
# sig_tmp_aric_pqtl_exp_outc 做coloc分析的数据，纳入p<5e-08的阈值
sig_tmp_aric_pqtl_exp_locuszoom_outc <- sumstats_dt_locuszoom %>% left_join(outc, by = c("SNP"="rsID"))

# remove duplicated SNPs
sig_tmp_aric_pqtl_exp_locuszoom_outc <- sig_tmp_aric_pqtl_exp_locuszoom_outc[order(sig_tmp_aric_pqtl_exp_locuszoom_outc[,'SNP'], sig_tmp_aric_pqtl_exp_locuszoom_outc[,'P']),]
sig_tmp_aric_pqtl_exp_locuszoom_outc <- sig_tmp_aric_pqtl_exp_locuszoom_outc[!duplicated(sig_tmp_aric_pqtl_exp_locuszoom_outc$SNP), ]
# remove T2D_EA missing SNPs
sig_tmp_aric_pqtl_exp_locuszoom_outc <- sig_tmp_aric_pqtl_exp_locuszoom_outc[!is.na(sig_tmp_aric_pqtl_exp_locuszoom_outc$Fixed.effects_beta),]

write.table(sig_tmp_aric_pqtl_exp_locuszoom_outc, file = 'F:/SH_desktop/pQTL-analysis/sig_tmp_aric_pqtl_exp_locuszoom_outc(protiens_t2d_coloc_input).csv', sep="\t", row.names = F, col.names = T, quote = F)


## Part 4. ====
##===== 4) locus zoom plot for variants with PPH4>0.8================

library(locuscomparer)
library(gwasglue) 

coloc_aric_proteins <- sig_coloc_aric_pqtl_t2d$protein

gwas_t2d <- NULL
pqtl_t2d <- NULL

for (i in coloc_aric_proteins){
  # dataset 1: pqtl summary; dataset 2: t2d
  print(paste0('step 1: select GWAS and pqtl data=====',i))
  
  
  tmp_coloc_aric_pqtl_t2d <- sig_tmp_aric_pqtl_exp_locuszoom_outc %>% filter(entrezgenesymbol == i)
  
  tmp_gwas_t2d <- tmp_coloc_aric_pqtl_t2d[,c('SNP','Fixed.effects_p.value')]
  tmp_gwas_t2d$outcome <- "T2D_EA"
  
  tmp_pqtl_t2d <- tmp_coloc_aric_pqtl_t2d[,c('SNP','P')]
  tmp_pqtl_t2d$outcome <- i
  
  gwas_t2d <- rbind(gwas_t2d, tmp_gwas_t2d)
  pqtl_t2d <- rbind(pqtl_t2d, tmp_pqtl_t2d)
}

gwas_t2d <- gwas_t2d %>% dplyr::rename(rsid = SNP, pval = Fixed.effects_p.value)
pqtl_t2d <- pqtl_t2d %>% dplyr::rename(rsid = SNP, pval = P)

#variants <- unique(sig_coloc_aric_pqtl_riskfactor_AA$snp)

#options(ggrepel.max.overlaps = Inf)

## locuszoom plot - proteins - T2D coloc PPH2>0.8 lucoszoom plot==========

# 1) GCKR & T2D_EA (GWAS pval = 1.621e-22)
tmp_pqtl_t2d <- pqtl_t2d %>% filter(outcome == "GCKR") %>% dplyr::select(-outcome)
tmp_gwas_t2d <- gwas_t2d %>% dplyr::select(-outcome)

pdf(file = "GCKR_T2D_zoomlocus_plot_EA_0218.pdf",
    width = 6,
    height = 6
)
locuscompare(in_fn1 = tmp_gwas_t2d,
             in_fn2 = tmp_pqtl_t2d,
             title1 = 'GWAS',
             title2 = 'pQTL',
             snp = "rs1260326")
dev.off() 



# 2) HHIP(PROK2) & T2D_EA (GWAS pval = 1.031e-06)
tmp_pqtl_t2d <- pqtl_t2d %>% filter(outcome == "HHIP") %>% dplyr::select(-outcome)

# check lead pqtl
sig_coloc_aric_pqtl_t2d[sig_coloc_aric_pqtl_t2d$protein == 'HHIP',]

pdf(file = "HHIP_T2D_zoomlocus_plot_EA_0218.pdf",
    width = 6,
    height = 6
)
locuscompare(in_fn1 = tmp_gwas_t2d,
             in_fn2 = tmp_pqtl_t2d,
             title1 = 'GWAS',
             title2 = 'pQTL',
             snp = "rs11727676")
dev.off()

# 3) NCAN(SORD,NRP2) & T2D_EA  (GWAS pval = 5.507e-08)
tmp_pqtl_t2d <- pqtl_t2d %>% filter(outcome == "NCAN") %>% dplyr::select(-outcome)

# check lead pqtl
#sig_coloc_aric_pqtl_t2d[sig_coloc_aric_pqtl_t2d$protein == 'NCAN',]

pdf(file = "NCAN_T2D_zoomlocus_plot_EA_0218.pdf",
    width = 6,
    height = 6
)
locuscompare(in_fn1 = tmp_gwas_t2d,
             in_fn2 = tmp_pqtl_t2d,
             title1 = 'GWAS',
             title2 = 'pQTL',
             snp = "rs2228603")
dev.off()


# 4) ANGPTL4(CFD,C5,GMPR2,GRHPR,PGP,OVCA2,PSME2,TBCB) & T2D_EA (GWAS pval = 2.826e-05)
tmp_pqtl_t2d <- pqtl_t2d %>% filter(outcome == "ANGPTL4") %>% dplyr::select(-outcome)

# check lead pqtl
#sig_coloc_aric_pqtl_t2d[sig_coloc_aric_pqtl_t2d$protein == 'PSME2',]

pdf(file = "ANGPTL4_T2D_zoomlocus_plot_EA_0218.pdf",
    width = 6,
    height = 6
)
locuscompare(in_fn1 = tmp_gwas_t2d,
             in_fn2 = tmp_pqtl_t2d,
             title1 = 'GWAS',
             title2 = 'pQTL',
             snp = "rs116843064")
dev.off()


# 5) NELL1(FMOD,CHRDL2,CLEC3B,CCL16) & T2D_EA  (GWAS pval = 0.0006519)
tmp_pqtl_t2d <- pqtl_t2d %>% filter(outcome == "NELL1") %>% dplyr::select(-outcome)

# check lead pqtl
#sig_coloc_aric_pqtl_t2d[sig_coloc_aric_pqtl_t2d$protein == 'PSME2',]

pdf(file = "NELL1_T2D_zoomlocus_plot_EA_0218.pdf",
    width = 6,
    height = 6
)
locuscompare(in_fn1 = tmp_gwas_t2d,
             in_fn2 = tmp_pqtl_t2d,
             title1 = 'GWAS',
             title2 = 'pQTL',
             snp = "rs8176786")
dev.off()


# 6) RARRES1(HYAL1,PLEK,WARS,B3GNT8) & T2D_EA (GWAS pval =0.0006234)
tmp_pqtl_t2d <- pqtl_t2d %>% filter(outcome == "RARRES1") %>% dplyr::select(-outcome)

# check lead pqtl
#sig_coloc_aric_pqtl_t2d[sig_coloc_aric_pqtl_t2d$protein == 'PSME2',]
pdf(file = "RARRES1_T2D_zoomlocus_plot_EA_0218.pdf",
    width = 6,
    height = 6
)
locuscompare(in_fn1 = tmp_gwas_t2d,
             in_fn2 = tmp_pqtl_t2d,
             title1 = 'GWAS',
             title2 = 'pQTL',
             snp = "rs4680458")
dev.off()



## part 5. =====
##===== 5) EA - Bubble plot visualization针对20个SIG protein&PP.H4>0.8 结果绘制 ===========================
p5_sig_aric_pqtl_t2d_ivs_coloc <- sig_aric_pqtl_t2d_ivs_coloc %>% dplyr::select(protein,SNP.PP.H4)
p5_sig_aric_pqtl_t2d_ivs_coloc$outcome <- "T2D_EA"
pdf("pqtl_coloc_bubble_plot(EA)_0427.pdf",
    width = 6,
    height = 3
)
p5_sig_aric_pqtl_t2d_ivs_coloc %>% 
  ggplot(aes(x = protein, y=outcome, size = SNP.PP.H4, na.rm = T, color="#0099B4ff", alpha=0.8)) +
  geom_point() +
  scale_color_manual(values = "#0099B4ff") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=8),
        legend.position = "right",
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #labs(colour = "Beta") +
  #coord_fixed(ratio = 1.2) +
  xlab("Proteins") + ylab(" ") +
  ggtitle("Colocalization analysis") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#################Summary：coloc-localization ARIC(EA 67-proteins) & T2D ==========9个pqtl========

###
# Perform: MR - analysis [ARIC(EA 67-proteins) & risk factors] ==================
###
# risk factors (SMK,SBP,DBP,LDL,HDL,TC,Fasting Glu, Fasting Insu, BMI/Obesity) 
# select sig IVs for the previous 20 proteins (sig MR & sig Coloc)

sig_aric_pqtl_riskfactor <- sig_aric_pqtl_exp_clump_nopleiotropy %>% 
  filter(exposure %in% intersect(sig_coloc_aric_pqtl_t2d$protein,consistent_proteins_robust_and_probable))

# BMI(ukb-b-19953);hba1c(ieu-b-4842);Fasting glucose(ebi-a-GCST90002232);Fasting insulin(ebi-a-GCST90002238);
# HDL(ieu-b-109/ieu-a-299);LDL(ieu-a-300/ieu-b-5089);SBP(ebi-a-GCST90018972);DBP(ebi-a-GCST90018952); 
# obesity(ukb-b-15541); sleeplessness/insomnia(ukb-b-3957)

riskfactors <- c('ukb-b-19953','ieu-b-4842','ebi-a-GCST90002232','ebi-a-GCST90002238',
                 'ieu-b-109','ieu-a-299','ieu-a-300','ieu-b-5089','ebi-a-GCST90018972',
                 'ebi-a-GCST90018952','ukb-b-15541','ukb-b-3957')

riskfactor_H_data <- NULL
riskfactor_mr_results <- NULL


for (riskfactor in riskfactors) {
  tmp_riskfactor_dat <- extract_outcome_data(
    snps = sig_aric_pqtl_riskfactor$SNP,
    outcomes = riskfactor,
    proxies = FALSE,
    access_token = NULL)
  print(paste0('extract outcome information====', riskfactor))
  
  tmp_riskfactor_H_data <- harmonise_data(
    exposure_dat = sig_aric_pqtl_riskfactor, 
    outcome_dat = tmp_riskfactor_dat
  )
  print('harmonise data')
  if(nrow(tmp_riskfactor_H_data) != 0){
    
    tmp_riskfactor_mr_results <- tmp_riskfactor_H_data %>% mr_singlesnp(
    parameters = default_parameters(),
    single_method = "mr_wald_ratio",
    all_method = c('mr_ivw_mre', 'mr_egger_regression'))
    
    riskfactor_H_data <- rbind(riskfactor_H_data, tmp_riskfactor_H_data)
  
    print('perform MR')
    
    riskfactor_mr_results <- rbind(riskfactor_mr_results,tmp_riskfactor_mr_results)
  }
} 

write.csv(riskfactor_H_data, file = 'aric_pqtl_riskfactor_H_data_clump_nopleiotropy(all)_0424.csv', row.names = F)
write.csv(riskfactor_mr_results, file = 'aric_pqtl_riskfactors_nopleiotropy_0424.csv')

## selecting significant risk-factors/

sig_riskfactor_mr_results <- riskfactor_mr_results %>% filter(
  SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
             'All - MR Egger') | p < 0.05/(20*10)
)

sig_riskfactor_mr_results <- sig_riskfactor_mr_results[!is.na(sig_riskfactor_mr_results$b),]


# ARIC(EA) add nsnp to the significant results (965-p1_sig_proteins)
riskfactor_proteins <- unique(sig_riskfactor_mr_results$exposure)

for (riskfactor_protein in riskfactor_proteins) {
  print(paste0('step 0:=========', riskfactor_protein))
  
  for (riskfactor in riskfactors) {
    tmp_sig_pqtl_riskfactor <- sig_riskfactor_mr_results %>% 
      filter(exposure == riskfactor_protein & id.outcome == riskfactor) %>%
      filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
               SNP != "All - MR Egger")
    print(paste0('step 1: filter snp results=====', riskfactor))
    
    if (nrow(tmp_sig_pqtl_riskfactor) != 0) {
      sig_riskfactor_mr_results$nSNP <- dim(tmp_sig_pqtl_riskfactor)[1]
      print(paste0('step 2: add snp numbers======', dim(tmp_sig_pqtl_riskfactor)[1]))
    }
  }
}




# 重新赋值方法，sig_riskfactor_mr_results 包含所有proteins的MR结果
sig_riskfactor_mr_results$mr_method <- ifelse(sig_riskfactor_mr_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                sig_riskfactor_mr_results$SNP != "All - MR Egger" & sig_riskfactor_mr_results$nSNP == 1, 
                                              "Wald_ratio", ifelse(sig_riskfactor_mr_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                     sig_riskfactor_mr_results$SNP != "All - MR Egger" & sig_riskfactor_mr_results$nSNP != 1, "Confused", 
                                                                   ifelse(sig_riskfactor_mr_results$SNP == "All - MR Egger", "MR_egger", "IVW(random_models)")))

sig_riskfactor_mr_results <- sig_riskfactor_mr_results %>% filter(sig_riskfactor_mr_results$mr_method != "Confused")
sig_riskfactor_mr_results <- sig_riskfactor_mr_results %>% generate_odds_ratios()

sig_riskfactor_mr_results <- sig_riskfactor_mr_results %>% filter(SNP != 'All - MR Egger')
sig_riskfactor_mr_results <- sig_riskfactor_mr_results[!is.na(sig_riskfactor_mr_results$b),]

# select significant MR results of 20 proteins and 10 risk factors (0.05/(48*10))
sig_riskfactor_mr_results <- sig_riskfactor_mr_results %>% filter(p < 0.05/(48*10))

# save file (MR-risk factors) 16个protein通过P<0.05/(20*10) threshold
write.csv(sig_riskfactor_mr_results, file = 'Sig_aric_pqtl_riskfactors_0424.csv', row.names = F)

# BMI(ukb-b-19953);hba1c(ieu-b-4842);Fasting glucose(ebi-a-GCST90002232);Fasting insulin(ebi-a-GCST90002238);
# HDL(ieu-b-109/ieu-a-299);LDL(ieu-a-300/ieu-b-5089);SBP(ebi-a-GCST90018972);DBP(ebi-a-GCST90018952); 
# obesity(ukb-b-15541); sleeplessness/insomnia(ukb-b-3957)

riskfactor_lab <- data.frame(
  rbind(c("ukb-b-19953","BMI"),
        c("ieu-b-4842","HAb1c"),
        c("ebi-a-GCST90002232", "Fasting glucose"),
        c("ebi-a-GCST90002238", "Fasting insulin"),
        c("ieu-b-109", "HDL"),
        c("ieu-a-299","HDL-C"),
        c("ieu-a-300",'LDL'),
        c('ieu-b-5089','LDL-C'),
        c('ebi-a-GCST90018972','SBP'),
        c('ebi-a-GCST90018952','DBP'),
        c('ukb-b-15541','Obesity'),
        c("ukb-b-3957",'Sleep disorder'))
)
colnames(riskfactor_lab) <- c('id.outcome','Outcomes')

p3_sig_riskfactor_mr_results <- sig_riskfactor_mr_results %>% left_join(riskfactor_lab, by="id.outcome")

# delete results of obesity of all missing values
p3_sig_riskfactor_mr_results <- p3_sig_riskfactor_mr_results %>% filter(Outcomes != "Obesity")

###
# Perform: coloc of risk factors [16 proteins * 9 risk factors]================================
###
# convert hg38 to hg37
## liftover sumstats_dt为HG37坐标系的ARIC中的SIG(MR-67) PQTL
sumstats_dt <- MungeSumstats::liftover(sumstats_dt = sig_tmp_aric_pqtl_exp/sig_aric_pqtl_exp_clump_nopleiotropy,
                                       ref_genome = "hg38",
                                       convert_ref_genome = "hg19")
knitr::kable(head(sumstats_dt))

# 
#sig_tmp_aric_pqtl_exp_riskfactors <- sumstats_dt %>% left_join(outc, by = c("SNP"="rsID"))

# 
sumstats_dt_exp <- format_data(
  sumstats_dt,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  effect_allele_col = "REF",
  other_allele_col = "A1",
  phenotype_col = "entrezgenesymbol"
)

sumstats_dt_exp <- sumstats_dt_exp %>% filter(exposure %in% riskfactor_proteins)
write.csv(sumstats_dt_exp, file = "sumstats_dt_exp.csv", row.names = F)

riskfactors <- c('ukb-b-19953','ieu-b-4842','ebi-a-GCST90002232','ebi-a-GCST90002238',
                 'ieu-b-109','ieu-b-5089','ebi-a-GCST90018972',
                 'ebi-a-GCST90018952','ukb-b-3957')
riskfactor_proteins <- unique(sig_aric_pqtl_riskfactor$exposure)



riskfactor_H_data <- NULL

for (riskfactor_protein in riskfactor_proteins) {
  
  tmp_sig_aric_pqtl_exp <- sig_aric_pqtl_exp_nopleiotropy %>% filter(exposure == riskfactor_protein)
  
  for (riskfactor in riskfactors) {
    tmp_riskfactor_outcome_dat <- extract_outcome_data(
      snps = tmp_sig_aric_pqtl_exp$SNP,
      outcomes = riskfactor,
      access_token = NULL,
      proxies = F)
    
    print(paste0(riskfactor_protein,"====", riskfactor))
    
    if(nrow(tmp_riskfactor_outcome_dat) !=0){
      tmp_riskfactor_H_data <- harmonise_data(
        exposure_dat = tmp_sig_aric_pqtl_exp, 
        outcome_dat = tmp_riskfactor_outcome_dat) 
    }
    print("harmonise data")
    
    if(nrow(tmp_riskfactor_H_data) != 0){
      riskfactor_H_data <- rbind(riskfactor_H_data, tmp_riskfactor_H_data)
    }
  }
}

#
write.csv(riskfactor_H_data, 'coloc_proteins_riskfactor_H_data_0424.csv',row.names = F)

# add additional info for proteins_riskfactors_H_data to perfrom the coloc analysis
coloc_aric_riskfactor_proteins <- unique(riskfactor_H_data$exposure)
riskfactors <- unique(riskfactor_H_data$id.outcome)

riskfactor_H_data$samplesize.outcome[riskfactor_H_data$id.outcome == "ieu-b-109"] <- 403943
riskfactor_H_data$samplesize.outcome[riskfactor_H_data$id.outcome == "ieu-b-5089"] <- 201678

coloc_aric_pqtl_riskfactor <- NULL

for (i in coloc_aric_riskfactor_proteins){
  # dataset 1: pqtl summary; dataset 2: t2d
  print(paste0('step 1: coloc-analysis=====',i))
  
  for (riskfactor in riskfactors) {
    tmp_coloc_aric_pqtl_riskfactor <- riskfactor_H_data %>% filter(exposure == i & id.outcome == riskfactor)
    print(paste0('exposure===',i,'; outcome====',riskfactor))
    
    
    tmp_coloc_aric_pqtl_riskfactor <- tmp_coloc_aric_pqtl_riskfactor[order(tmp_coloc_aric_pqtl_riskfactor[,'SNP'], tmp_coloc_aric_pqtl_riskfactor[,'pval.exposure']),]
    tmp_coloc_aric_pqtl_riskfactor <- tmp_coloc_aric_pqtl_riskfactor[!duplicated(tmp_coloc_aric_pqtl_riskfactor$SNP), ]
    print('step 2: remove duplicates')
      
    tmp_coloc_aric_pqtl_riskfactor <- tmp_coloc_aric_pqtl_riskfactor %>% filter(!is.na(eaf.outcome))
    tmp_coloc_aric_pqtl_riskfactor <- tmp_coloc_aric_pqtl_riskfactor %>% filter(!is.na(eaf.outcome))
      #tmp_coloc_aric_pqtl_riskfactor <- tmp_coloc_aric_pqtl_riskfactor[tmp_coloc_aric_pqtl_riskfactor$pval.exposure != 0]
    print('step 3: remove missing values')
      
    if(nrow(tmp_coloc_aric_pqtl_riskfactor)!=0){
         tmp <- coloc.abf(dataset2 = list(pvalues=tmp_coloc_aric_pqtl_riskfactor$pval.outcome, 
                                       type = "quant",  
                                       N = tmp_coloc_aric_pqtl_riskfactor$samplesize.outcome, 
                                       snp=tmp_coloc_aric_pqtl_riskfactor$SNP),
                       dataset1 = list(pvalues=tmp_coloc_aric_pqtl_riskfactor$pval.exposure, 
                                       type ="quant", 
                                       N = 7213, 
                                       snp = tmp_coloc_aric_pqtl_riskfactor$SNP), 
                       MAF = tmp_coloc_aric_pqtl_riskfactor$eaf.outcome)
      
          print(paste0('step 4: coloc-analysis outcome=====',riskfactor))
          tmp <- tmp$results
          if (nrow(tmp) != 0) {
            tmp$outcome <- riskfactor
            tmp$protein <- i
            coloc_aric_pqtl_riskfactor <- rbind(coloc_aric_pqtl_riskfactor, tmp)
      }
    }
  }
}

write.csv(coloc_aric_pqtl_riskfactor, file = "coloc_aric_20_coloc_proteins_riskfactor_0424.csv", row.names = F)

sig_coloc_aric_pqtl_riskfactor <- coloc_aric_pqtl_riskfactor %>% dplyr::filter(SNP.PP.H4>=0.8)
sig_coloc_aric_pqtl_riskfactor <- sig_coloc_aric_pqtl_riskfactor %>% dplyr::rename(exposure=protein,id.outcome=outcome)

# sig_aric_pqtl_riskfactor_ivs_coloc 为pass PP.H4 和 MR(proteins - riskfactors) analysis且 共定位SNP在mr analysis中
sig_aric_pqtl_riskfactor_ivs_coloc <- sig_coloc_aric_pqtl_riskfactor %>% 
  filter(protein %in% sig_aric_pqtl_t2d_ivs_coloc$protein) %>% 
  filter(snp %in% sig_aric_pqtl_exp_clump_nopleiotropy$SNP)
sig_aric_pqtl_riskfactor_ivs_coloc <- sig_aric_pqtl_riskfactor_ivs_coloc %>% 
  left_join(riskfactor_lab, by=c('outcome' = 'id.outcome'))

sig_aric_pqtl_riskfactor_ivs_coloc <- sig_aric_pqtl_riskfactor_ivs_coloc %>% filter(protein %in% p3_sig_riskfactor_mr_results.1$exposure)

# p3_sig_riskfactor_mr_results 为MR 结果显著 且 有共定位证据的
p3_sig_riskfactor_mr_results.1 <- p3_sig_riskfactor_mr_results %>% left_join(sig_coloc_aric_pqtl_riskfactor[,c('id.outcome','exposure','SNP.PP.H4')],
                                                                           by = c('exposure','id.outcome'))
# p3_sig_riskfactor_mr_results.2 为MR显著，且有共定位证据显示 pqtl与iv重叠,TEK无显著PPH4且MR显著的结果
p3_sig_riskfactor_mr_results.1$exposure <- ifelse(p3_sig_riskfactor_mr_results.1$exposure %in% unique(sig_aric_pqtl_riskfactor_ivs_coloc$protein)[1:6],
                                                  paste0(p3_sig_riskfactor_mr_results.1$exposure,"*"),
                                                  p3_sig_riskfactor_mr_results.1$exposure)

## Part 3. ============
##===== 3) Bubble plot visualization针对20个SIG protein&PP.H4>0.8中的16个MR sig(risk factors)结果绘制 ===========================
p3_sig_riskfactor_mr_results.2 <- p3_sig_riskfactor_mr_results.1 %>% filter(Outcomes != "LDL-C")
pdf("pqtl_riskfactor_bubble_plot(EA)_0428.pdf",
    width = 6,
    height = 6
)
p3_sig_riskfactor_mr_results.2 %>% 
  ggplot(aes(x = exposure, y = Outcomes, color=b, na.rm = T)) +
  geom_point() +
  geom_point(shape = 8, aes(size = SNP.PP.H4)) +
  scale_color_continuous_diverging(palette = "Blue-Red2") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=8),
        legend.position = "right",
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(colour = "Beta") +
  #coord_fixed(ratio = 1.2) +
  xlab(" ") + ylab(" ") +
  ggtitle(expression("MR: Proteins"  %->% "Risk factors(European Ancestry)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
############## 2024-04-27 summary: CDNF/GCKR/MDGA2/SERPINA1 为进行多COLOC分析的四个蛋白质========================================

###
# perform multi-coloc analysis
###

# protein - T2D [sig_tmp_aric_pqtl_exp_outc]
# protein - risk factors [riskfactor_H_data]

multi_coloc_pqtl_t2d_riskfactor <- riskfactor_H_data %>% left_join(outc, by = c("SNP"="rsID"))


multi_coloc_pqtl_t2d_riskfactor_H_data <- NULL
multi_coloc_proteins <- c('GCKR','SERPINA1')
multi_coloc_riskfactors <- c('ebi-a-GCST90002232','ebi-a-GCST90002238',"ukb-b-19953",
                             'ebi-a-GCST90018972','ebi-a-GCST90018952')

sig_tmp_aric_pqtl_exp_outc_format <- sig_tmp_aric_pqtl_exp_outc %>% 
  filter(entrezgenesymbol %in% multi_coloc_proteins) %>%
  format_data(
    type = "exposure",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P",
    effect_allele_col = "REF",
    other_allele_col = "A1",
    phenotype_col = "entrezgenesymbol"
  )

for (riskfactor_protein in multi_coloc_proteins) {
  
  tmp_sig_aric_pqtl_exp <- sig_tmp_aric_pqtl_exp_outc_format %>% filter(exposure  == riskfactor_protein)
  
  print(riskfactor_protein)
  
  for (riskfactor in multi_coloc_riskfactors) {
    tmp_riskfactor_outcome_dat <- extract_outcome_data(
      snps = tmp_sig_aric_pqtl_exp$SNP,
      outcomes = multi_coloc_riskfactors,
      access_token = NULL,
      proxies = F)
    
    print(paste0(riskfactor_protein,"====", multi_coloc_riskfactors))
    
    if(nrow(tmp_riskfactor_outcome_dat) !=0){
      tmp_riskfactor_H_data <- harmonise_data(
        exposure_dat = tmp_sig_aric_pqtl_exp, 
        outcome_dat = tmp_riskfactor_outcome_dat) 
    }
    print("harmonise data")
    
    if(nrow(tmp_riskfactor_H_data) != 0){
      multi_coloc_pqtl_t2d_riskfactor_H_data <- rbind(multi_coloc_pqtl_t2d_riskfactor_H_data, tmp_riskfactor_H_data)
    }
  }
}

multi_coloc_pqtl_t2d_riskfactor <- multi_coloc_pqtl_t2d_riskfactor_H_data %>% 
  dplyr::select(SNP, beta.exposure,se.exposure,beta.outcome,se.outcome,exposure,id.outcome) %>% 
  left_join(riskfactor_lab, by="id.outcome")

# 2024-04-27 初步测试可能有显著意义的结果为：GCKR/SERPINA1 
# 分别建立 执行MULTI_COLOC分析的数据集：'GCKR','SERPINA1'
multi_coloc_riskfactor_outcome_dat <- multi_coloc_pqtl_t2d_riskfactor %>%
  filter(exposure == "SERPINA1") %>% dplyr::select(SNP,beta.exposure,se.exposure)
multi_coloc_riskfactor_outcome_dat <- unique(multi_coloc_riskfactor_outcome_dat)

#
multi_coloc_riskfactors <- c('ebi-a-GCST90002232','ebi-a-GCST90002238',"ukb-b-19953",
                             'ebi-a-GCST90018972','ebi-a-GCST90018952')
multi_coloc_riskfactors<-multi_coloc_riskfactors[5]
for (riskfactor in multi_coloc_riskfactors) {
    tmp_riskfactor_outcome_dat <- multi_coloc_pqtl_t2d_riskfactor %>% 
      filter(id.outcome == multi_coloc_riskfactors) %>%
      dplyr::select(SNP,beta.outcome,se.outcome)
    print(paste0('outcome====', riskfactor))
    tmp_riskfactor_outcome_dat <- unique(tmp_riskfactor_outcome_dat)
    multi_coloc_riskfactor_outcome_dat <- left_join(multi_coloc_riskfactor_outcome_dat,
                                                    tmp_riskfactor_outcome_dat, by="SNP", 
                                                    suffix =c("", paste0("_",riskfactor)))
}

riskfactor_lab <- data.frame(
  rbind(c("ukb-b-19953","BMI"),
        c("ieu-b-4842","HAb1c"),
        c("ebi-a-GCST90002232", "Fasting glucose"),
        c("ebi-a-GCST90002238", "Fasting insulin"),
        c("ieu-b-109", "HDL"),
        c("ieu-a-299","HDL-C"),
        c("ieu-a-300",'LDL'),
        c('ieu-b-5089','LDL-C'),
        c('ebi-a-GCST90018972','SBP'),
        c('ebi-a-GCST90018952','DBP'),
        c('ukb-b-15541','Obesity'),
        c("ukb-b-3957",'Sleep disorder'))
)
colnames(riskfactor_lab) <- c('id.outcome','Outcomes')

multi_coloc_riskfactor_outcome_dat_GCKR <- multi_coloc_riskfactor_outcome_dat %>% 
  dplyr::rename(beta.fg=beta.outcome,se.fg=se.outcome,
                beta.fi=`beta.outcome_ebi-a-GCST90002238`,
                se.fi= `se.outcome_ebi-a-GCST90002238`,
                beta.bmi=`beta.outcome_ukb-b-19953`,
                se.bmi=`se.outcome_ukb-b-19953`)

multi_coloc_riskfactor_outcome_dat_GCKR <- multi_coloc_riskfactor_outcome_dat_GCKR %>%
  left_join(outc[,c('rsID','Fixed.effects_beta','Fixed.effects_SE')], by = c("SNP"="rsID"))

multi_coloc_riskfactor_outcome_dat_GCKR <- multi_coloc_riskfactor_outcome_dat_GCKR[!duplicated(multi_coloc_riskfactor_outcome_dat_GCKR$SNP),]
multi_coloc_riskfactor_outcome_dat_GCKR <- multi_coloc_riskfactor_outcome_dat_GCKR%>% na.omit()

# define betas
betas <- data.frame(SNPs = multi_coloc_riskfactor_outcome_dat_GCKR$SNP,
                    T1 = multi_coloc_riskfactor_outcome_dat_GCKR$beta.exposure,
                    T2 = multi_coloc_riskfactor_outcome_dat_GCKR$Fixed.effects_beta,
                    T3 = multi_coloc_riskfactor_outcome_dat_GCKR$beta.fg,
                    T4 = multi_coloc_riskfactor_outcome_dat_GCKR$beta.fi,
                    T5 = multi_coloc_riskfactor_outcome_dat_GCKR$beta.bmi)
rownames(betas) <- betas$SNPs
betas=subset(betas,select=-SNPs)
# define ses
ses <- data.frame(SNPs = multi_coloc_riskfactor_outcome_dat_GCKR$SNP,
                  T1 = multi_coloc_riskfactor_outcome_dat_GCKR$se.exposure,
                  T2 = multi_coloc_riskfactor_outcome_dat_GCKR$Fixed.effects_SE,
                  T3 = multi_coloc_riskfactor_outcome_dat_GCKR$se.fg,
                  T4 = multi_coloc_riskfactor_outcome_dat_GCKR$se.fi,
                  T5 = multi_coloc_riskfactor_outcome_dat_GCKR$se.bmi)
rownames(ses)<-ses$SNPs
ses=subset(ses,select=-SNPs)
# define traits
traits <-  paste0("T", 1:dim(betas)[2])
rsid <- rownames(betas)
betas <- data.matrix(betas)
ses <- data.matrix(ses)

res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
res

# perform multi coloc for SERPINA1
multi_coloc_riskfactor_outcome_dat_SERPINA1 <- multi_coloc_riskfactor_outcome_dat %>% 
  dplyr::rename(beta.fg=beta.outcome,se.fg=se.outcome,
                beta.fi=`beta.outcome_ebi-a-GCST90002238`,
                se.fi= `se.outcome_ebi-a-GCST90002238`,
                beta.bmi=`beta.outcome_ukb-b-19953`,
                se.bmi=`se.outcome_ukb-b-19953`,
                beta.sbp=`beta.outcome_ebi-a-GCST90018972`,
                se.sbp=`se.outcome_ebi-a-GCST90018972`,
                beta.dbp=`beta.outcome_ebi-a-GCST90018952`,
                se.dbp=`se.outcome_ebi-a-GCST90018952`)
# combine with ourcome(T2D)
multi_coloc_riskfactor_outcome_dat_SERPINA1 <- multi_coloc_riskfactor_outcome_dat_SERPINA1 %>%
  left_join(outc[,c('rsID','Fixed.effects_beta','Fixed.effects_SE')], by = c("SNP"="rsID"))

multi_coloc_riskfactor_outcome_dat_SERPINA1 <- multi_coloc_riskfactor_outcome_dat_SERPINA1[!duplicated(multi_coloc_riskfactor_outcome_dat_SERPINA1$SNP),]
multi_coloc_riskfactor_outcome_dat_SERPINA1 <- multi_coloc_riskfactor_outcome_dat_SERPINA1%>% na.omit()

# define betas
betas_SERPINA1 <- data.frame(SNPs = multi_coloc_riskfactor_outcome_dat_SERPINA1$SNP,
                    T1 = multi_coloc_riskfactor_outcome_dat_SERPINA1$beta.exposure,
                    T2 = multi_coloc_riskfactor_outcome_dat_SERPINA1$Fixed.effects_beta,
                    T3 = multi_coloc_riskfactor_outcome_dat_SERPINA1$beta.fg,
                    T4 = multi_coloc_riskfactor_outcome_dat_SERPINA1$beta.fi,
                    T5 = multi_coloc_riskfactor_outcome_dat_SERPINA1$beta.bmi)
rownames(betas_SERPINA1) <- betas_SERPINA1$SNPs
betas_SERPINA1=subset(betas_SERPINA1,select=-SNPs)
# define ses
ses_SERPINA1 <- data.frame(SNPs = multi_coloc_riskfactor_outcome_dat_SERPINA1$SNP,
                  T1 = multi_coloc_riskfactor_outcome_dat_SERPINA1$se.exposure,
                  T2 = multi_coloc_riskfactor_outcome_dat_SERPINA1$Fixed.effects_SE,
                  T3 = multi_coloc_riskfactor_outcome_dat_SERPINA1$se.fg,
                  T4 = multi_coloc_riskfactor_outcome_dat_SERPINA1$se.fi,
                  T5 = multi_coloc_riskfactor_outcome_dat_SERPINA1$se.bmi)
rownames(ses_SERPINA1)<-ses_SERPINA1$SNPs
ses_SERPINA1=subset(ses_SERPINA1,select=-SNPs)
# define traits
traits_SERPINA1 <-  paste0("T", 1:dim(betas_SERPINA1)[2])
rsid_SERPINA1 <- rownames(betas_SERPINA1)
betas_SERPINA1 <- data.matrix(betas_SERPINA1)
ses_SERPINA1 <- data.matrix(ses_SERPINA1)

res_SERPINA1 <- hyprcoloc(betas_SERPINA1, ses_SERPINA1, trait.names=traits_SERPINA1, snp.id=rsid_SERPINA1)
res_SERPINA1


###
# perform: MR - mediation analysis (step2---risk factors & T2D) (European Ancestry)==============================================
###
riskfactors <- c('ukb-b-19953','ieu-b-4842','ebi-a-GCST90002232','ebi-a-GCST90002238',
                 'ieu-b-109','ieu-b-5089','ebi-a-GCST90018972',
                 'ebi-a-GCST90018952','ukb-b-3957')

riskfactor_t2d_mr_results <- NULL
riskfactor_exposure_dat <- NULL
t2d_outcome_dat <- NULL
riskfactor_t2d_H_data <- NULL

for (riskfactor in riskfactors) {
  tmp_riskfactor_exposure_dat <- extract_instruments(
    outcomes = riskfactor,
    access_token = NULL)
  print(paste0('extract exposure snps=====',riskfactor))
  
  t2d_outcome_dat <- format_data(
    type = "outcome",
    dat = outc,
    snps = tmp_riskfactor_exposure_dat$SNP,
    header=T,
    snp_col = 'rsID',
    beta_col = 'Fixed.effects_beta',
    se_col = 'Fixed.effects_SE',
    pval_col = 'Fixed.effects_p.value',
    effect_allele_col = 'effect_allele',
    other_allele_col = 'other_allele',
    chr_col = 'chromosome',
    pos_col = 'position',
    eaf_col = 'effect_allele_frequency'
  )
  print('format outcome')
  tmp_riskfactor_t2d_H_data <- harmonise_data(
    exposure_dat = tmp_riskfactor_exposure_dat, 
    outcome_dat = t2d_outcome_dat
  )
  print('harmonise data')
  riskfactor_t2d_H_data <- rbind(riskfactor_t2d_H_data, tmp_riskfactor_t2d_H_data)
  
  tmp_riskfactor_t2d_mr_results <- tmp_riskfactor_t2d_H_data %>% mr_singlesnp(
    parameters = default_parameters(),
    single_method = "mr_wald_ratio",
    all_method = c('mr_ivw_mre', 'mr_egger_regression'))
  
  print('perform MR')
  
  if(!is.null(tmp_riskfactor_t2d_mr_results)){
    riskfactor_t2d_mr_results <- rbind(riskfactor_t2d_mr_results,tmp_riskfactor_t2d_mr_results)
  }
} 

riskfactor_t2d_mr_results$outcome <- "T2D_EA"
riskfactor_t2d_mr_results <- riskfactor_t2d_mr_results %>% generate_odds_ratios() 
# save data
write.csv(riskfactor_t2d_mr_results, file = "F:/SH_desktop/pQTL-analysis/riskfactor_t2d_mr_results.csv")

sig_riskfactor_t2d_mr_results <- riskfactor_t2d_mr_results %>% 
  filter(SNP %in% c("All - Inverse variance weighted (multiplicative random effects)",
                    "All - MR Egger")) %>% filter(p < 0.05/9)


riskfactor_t2d_exposure <- unique(riskfactor_t2d_mr_results$id.exposure)

for (riskfactor in riskfactor_t2d_exposure) {
  
  tmp_riskfactor_t2d <- riskfactor_t2d_mr_results %>% 
    filter(id.exposure == riskfactor) %>%
    filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
             SNP != "All - MR Egger")
  print(paste0('step 1: filter snp results=====', riskfactor))
  
  if (nrow(tmp_riskfactor_t2d) != 0) {
    riskfactor_t2d_mr_results$nSNP[riskfactor_t2d_mr_results$id.exposure == riskfactor] <- dim(tmp_riskfactor_t2d)[1]
    print(paste0('step 2: add snp numbers======', dim(tmp_riskfactor_t2d)[1]))
  }
}
write.csv(riskfactor_t2d_mr_results, file = "riskfactor_t2d_mr_results_0424.csv")


# 重新赋值方法，riskfactor_t2d_mr_results_AA 包含所有proteins的MR结果
riskfactor_t2d_mr_results$mr_method <- ifelse(riskfactor_t2d_mr_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                riskfactor_t2d_mr_results$SNP != "All - MR Egger" & riskfactor_t2d_mr_results$nSNP == 1, 
                                                 "Wald_ratio", ifelse(riskfactor_t2d_mr_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                        riskfactor_t2d_mr_results$SNP != "All - MR Egger" & riskfactor_t2d_mr_results$nSNP != 1, "Confused", 
                                                                      ifelse(riskfactor_t2d_mr_results$SNP == "All - MR Egger", "MR_egger", "IVW(random_models)")))

riskfactor_t2d_mr_results <- riskfactor_t2d_mr_results %>% filter(riskfactor_t2d_mr_results$mr_method != "Confused")
#riskfactor_t2d_mr_results <- riskfactor_t2d_mr_results %>% generate_odds_ratios()

#riskfactor_t2d_mr_results_AA <- riskfactor_t2d_mr_results_AA %>% filter(SNP != 'All - MR Egger')
riskfactor_t2d_mr_results <- riskfactor_t2d_mr_results[!is.na(riskfactor_t2d_mr_results$b),]


##===== 5) EA-MR (riskfactors - T2D) EA forest plot visualization 结果绘制 ===========================
# forest plot (EA)====
p4_riskfactor_t2d_mr_results <- riskfactor_t2d_mr_results %>% filter(mr_method %in% c('IVW(random_models)',"Wald_ratio"))

p4_riskfactor_t2d_mr_results$orlab <- paste0("",
                                                ifelse(sprintf("%.2f",p4_riskfactor_t2d_mr_results$or)<0.0051,
                                                       format(p4_riskfactor_t2d_mr_results$or,scientific = TRUE,digits=3),
                                                       sprintf("%.2f",p4_riskfactor_t2d_mr_results$or)),
                                                "(", 
                                                ifelse(sprintf("%.2f",p4_riskfactor_t2d_mr_results$or_lci95)<0.0051,
                                                       format(p4_riskfactor_t2d_mr_results$or_lci95,scientific = TRUE,digits=3),
                                                       sprintf("%.2f",p4_riskfactor_t2d_mr_results$or_lci95)),
                                                "-",
                                                ifelse(sprintf("%.2f",p4_riskfactor_t2d_mr_results$or_uci95)<0.0051,
                                                       format(p4_riskfactor_t2d_mr_results$or_uci95,scientific = TRUE,digits=3),
                                                       sprintf("%.2f",p4_riskfactor_t2d_mr_results$or_uci95)),
                                                ")")

p4_riskfactor_t2d_mr_results <- p4_riskfactor_t2d_mr_results[order(p4_riskfactor_t2d_mr_results[,'or']),]
p4_riskfactor_t2d_mr_results$exposure <- factor(p4_riskfactor_t2d_mr_results$exposure, levels = p4_riskfactor_t2d_mr_results$exposure)

riskfactor_lab <- data.frame(
  rbind(c("ukb-b-19953","BMI"),
        c("ieu-b-4842","HAb1c"),
        c("ebi-a-GCST90002232", "Fasting glucose"),
        c("ebi-a-GCST90002238", "Fasting insulin"),
        c("ieu-b-109", "HDL"),
        c('ieu-b-5089','LDL-C'),
        c('ebi-a-GCST90018972','SBP'),
        c('ebi-a-GCST90018952','DBP'),
        c("ukb-b-3957",'Sleep disorder'))
)
colnames(riskfactor_lab) <- c('id.exposure','Outcomes')

p4_riskfactor_t2d_mr_results <- p4_riskfactor_t2d_mr_results %>% left_join(riskfactor_lab, by="id.exposure")
p4_riskfactor_t2d_mr_results <- p4_riskfactor_t2d_mr_results %>% filter(!(Outcomes %in% c("LDL-C",'Sleep disorder')))

p4_riskfactor_t2d_mr_results$x <- 0
p4_riskfactor_t2d_mr_results$y <- seq(1,7,1)

# annotation 
riskfactor_annotation <- NULL
riskfactor_annotation0 <- p4_riskfactor_t2d_mr_results[c('x','y','nSNP')] 
riskfactor_annotation1 <- p4_riskfactor_t2d_mr_results[c('x','y','orlab')] %>% dplyr::rename(nSNP=orlab) %>% rbind(riskfactor_annotation0)
riskfactor_annotation2 <- p4_riskfactor_t2d_mr_results[c('x','y','Outcomes')] %>% dplyr::rename(nSNP=Outcomes) 
riskfactor_annotation <- riskfactor_annotation1 %>% rbind(riskfactor_annotation2) 

riskfactor_annotation[1:7,1] <- -0.5
riskfactor_annotation[8:14,1] <- -1.9
riskfactor_annotation[15:21,1] <- -3.2


riskfactor_annotation <- rbind(c(-0.5,7.4,"OR(95%CI)"),riskfactor_annotation)
riskfactor_annotation <- rbind(c(-1.9,7.4,"nSNP"),riskfactor_annotation)
riskfactor_annotation <- rbind(c(-3.2,7.4,"Risk factors"),riskfactor_annotation)

riskfactor_annotation$x <- as.numeric(riskfactor_annotation$x)
riskfactor_annotation$y <- as.numeric(riskfactor_annotation$y)
# AA - riskfactor & T2D forest plot====
pdf("riskfactor_t2d_forest_plot(EA)_0424.pdf",
    width = 8,
    height = 8
)
p4_riskfactor_t2d_mr_results %>% 
  ggplot(aes(or, exposure)) +
  geom_point(aes(col=or, size = nSNP)) +
  #scale_color_gradient(low = "#5E4FA2", high = "#9E0142") +
  scale_color_continuous_sequential(palette = "ag_Sunset") +
  #scale_color_manual(values = c("#8F2618",'#5B3A91')) +
  geom_errorbarh(aes(xmax=or_lci95, xmin=or_uci95, col = or), height=0.2, size=1) +
  scale_x_continuous(limits = c(-4,6), breaks = seq(1,6,1)) +
  geom_vline(aes(xintercept = 1)) + theme_bw()+
  theme(axis.text.y = element_blank(), axis.ticks.y =element_blank()) +
  theme(axis.line.x = element_line(color = "#808080")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), plot.background = element_rect(color="black",size = 0.5)) +
  theme(legend.text = element_text(size = 4), legend.title = element_text(size = 6)) +
  geom_text(data = riskfactor_annotation, aes(x=x, y=y, label=nSNP)) +
  #labs(colour = "P value") +
  #coord_fixed(ratio = 1.2) +
  xlab("Odds Ratio") + ylab(" ") +
  ggtitle(expression("MR: Risk factors"  %->% "T2D(EA)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###
# Perform: PWAS - comparison results between protome-MR and PWAS results from FUSION (European Ancestry)=====================================================
###
# 1) INTERVAL (blood/plasma proteins)
# Define INTERVAL data path
interval_path <- "FUSHION_pwas/IGAP_PWAS/INTERVAL/"
interval_data <- ".top"
intervals <- list.files(path=interval_path,pattern = paste0("*",interval_data))

# 
interval_pwas_top <- NULL
for (i in intervals){
  tmp <- suppressWarnings(fread(paste0(interval_path,i),
                                stringsAsFactors = F,
                                data.table = F))
  if (nrow(tmp) != 0) {
    tmp$info <- gsub(interval_data,"",i)
    interval_pwas_top <- rbind(interval_pwas_top, tmp)
  }
}
# 2024-05-08 adjusted the top pwas results (0.05/total proteins in interval)
interval_pwas_top <- fread('/Volumes/Echo_SSD/SH_desktop/Proteome_MR/t2d.INTERVAL.PWAS.combined.top')

interval_pwas_top$proteins <- sapply(strsplit(interval_pwas_top$ID, "-"), function(x) x[2])

intersect(interval_pwas_top$proteins, sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$exposure)

# 2) ROSMAP (brain proteins)
rosmap_path <- "F:/SH_desktop/pQTL-analysis/FUSHION_pwas/IGAP_PWAS/ROSMAP/"
rosmap_data <- ".top"
rosmaps <- list.files(path=rosmap_path,pattern = paste0("*",rosmap_data))

# 
rosmap_pwas_top <- NULL
for (i in rosmaps){
  tmp <- suppressWarnings(fread(paste0(rosmap_path,i),
                                stringsAsFactors = F,
                                data.table = F))
  if (nrow(tmp) != 0) {
    tmp$info <- gsub(rosmap_data,"",i)
    rosmap_pwas_top <- rbind(rosmap_pwas_top, tmp)
  }
}

rosmap_pwas_top <- rosmap_pwas_top %>% separate(col = ID, into = c("gene_id", "proteins"), sep = "\\.")

intersect(sig_rosmap_brain_pqtl_t2d_singlesnp_results_clump$exposure,rosmap_pwas_top$proteins)

###
# perform literature triangulation analysis=======================================================
###
library("magrittr")
library("dplyr")
library("purrr")
library("glue")
library("ggplot2")
library("igraph")
library("epigraphdb")

# set the trait
STARTING_TRAIT <- "Sleep duration"

get_mr <- function(trait) {
  endpoint <- "/mr"
  params <- list(
    exposure_trait = trait,
    pval_threshold = 1e-10
  )
  mr_df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  mr_df
}

mr_df <- get_mr(STARTING_TRAIT)
mr_df %>% glimpse()

# map outcome traits to disease

trait_to_disease <- function(trait) {
  endpoint <- "/ontology/gwas-efo-disease"
  params <- list(trait = trait)
  disease_df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  if (nrow(disease_df) > 0) {
    res <- disease_df %>% pull(`disease.label`)
  } else {
    res <- c()
  }
  res
}

outcome_disease_df <- mr_df %>%
  select(`outcome.trait`) %>%
  distinct() %>%
  mutate(disease = map(`outcome.trait`, trait_to_disease)) %>%
  filter(map_lgl(`disease`, function(x) !is.null(x)))
outcome_disease_df


mr_df <- mr(outcome_trait = "Type 2 diabetes")
mr_df %>% glimpse()

exposure_disease_df <- mr_df %>%
  dplyr::select(`exposure.trait`) %>%
  distinct() %>%
  mutate(disease = map(`exposure.trait`, trait_to_disease)) %>%
  filter(map_lgl(`disease`, function(x) !is.null(x)))

mr_df[grep("Total protein",mr_df$exposure.trait, value = T),]

# look for literature evidence
get_gwas_pair_literature <- function(gwas_id, assoc_gwas_id) {
  endpoint <- "/literature/gwas/pairwise"
  # NOTE in this example we blacklist to semmentic types
  params <- list(
    gwas_id = gwas_id,
    assoc_gwas_id = assoc_gwas_id,
    by_gwas_id = TRUE,
    pval_threshold = 1e-1,
    semmantic_types = "nusq",
    semmantic_types = "dsyn",
    blacklist = TRUE,
    limit = 1000
  )
  lit_df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  lit_df
}

# GCKR - prot-a-1184; prot-a-2172 -
GWAS_ID_X <- "prot-a-2172"
GWAS_ID_Y <- "ieu-a-25"

lit_df <- get_gwas_pair_literature(GWAS_ID_X, GWAS_ID_Y)
glimpse(lit_df)

lit_df %>%
  count(`s1.predicate`) %>%
  arrange(desc(n))

lit_df %>%
  count(`s2.predicate`) %>%
  arrange(desc(n))

# filter out some predicates that are not informative
pred_filter <- c("COEXISTS_WITH", "ASSOCIATED_WITH")

lit_df_filter <- lit_df %>%
  filter(
    !`s1.predicate` %in% pred_filter,
    !`s2.predicate` %in% pred_filter
  )
# literature results
lit_counts <- lit_df_filter %>%
  count(`st.type`, `st.name`) %>%
  arrange(`st.type`, desc(`n`))

lit_counts %>% print(n = 30)

# plotting the overlapping term frequencies
lit_counts %>%
  filter(n < 300) %>%
  {
    ggplot(.) +
      aes(x = `st.name`, y = `n`) +
      geom_col() +
      geom_text(
        aes(label = `n`),
        position = position_dodge(0.9),
        hjust = 0
      ) +
      coord_flip()
  }

# look in detail at one overlapping term 
focus_term <- "Enzymes"
lit_detail <- lit_df_filter %>% filter(`st.name` == focus_term)
lit_detail %>% head()

lit_detail <- lit_detail %>%
  mutate_at(vars(`gwas.trait`, `assoc_gwas.trait`), stringr::str_to_upper)

nodes <- bind_rows(
  lit_detail %>% dplyr::select(node = `gwas.trait`) %>% distinct() %>% mutate(node_type = 1),
  lit_detail %>% dplyr::select(node = `assoc_gwas.trait`) %>% distinct() %>% mutate(node_type = 1),
  lit_detail %>% dplyr::select(node = `st1.name`) %>% distinct() %>% mutate(node_type = 2),
  lit_detail %>% dplyr::select(node = `st2.name`) %>% distinct() %>% mutate(node_type = 2),
  lit_detail %>% dplyr::select(node = `st.name`) %>% distinct() %>% mutate(node_type = 3),
) %>% distinct()

edges <- bind_rows(
  # exposure -> s1 subject
  lit_detail %>%
    dplyr::select(node = `gwas.trait`, assoc_node = `st1.name`) %>%
    distinct(),
  # s2 object -> outcome
  lit_detail %>%
    dplyr::select(node = `st2.name`, assoc_node = `assoc_gwas.trait`) %>%
    distinct(),
  # s1 subject - s1 predicate -> s1 object
  lit_detail %>%
    dplyr::select(
      node = `st1.name`, assoc_node = `st.name`,
      label = `s1.predicate`
    ) %>%
    distinct(),
  # s2 subject - s2 predicate -> s2 object
  lit_detail %>%
    dplyr::select(
      node = `st.name`, assoc_node = `st2.name`,
      label = `s2.predicate`
    ) %>%
    distinct()
) %>%
  distinct()

plot_network <- function(edges, nodes) {
  graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  graph$layout <- layout_with_kk
  
  # generate colors based on node type
  colors <- c("tomato", "lightblue", "gold")
  V(graph)$color <- colors[V(graph)$node_type]
  
  # Configure canvas
  default_mar <- par("mar")
  new_mar <- c(0, 0, 0, 0)
  par(mar = new_mar)
  
  plot.igraph(
    graph,
    vertex.size = 13,
    vertex.label.color = "black",
    vertex.label.family = "Helvetica",
    vertex.label.cex = 0.8,
    edge.arrow.size = 0.4,
    edge.label.color = "black",
    edge.label.family = "Helvetica",
    edge.label.cex = 0.5
  )
  par(mar = default_mar)
}

plot_network(edges, nodes)

###
# perform eQTL based MR analysis (eQTLGen genes - T2D)==========================
###
ieu_list <- read.csv('D:/SH_desktop/Proteome_MR/IEU_opengwas_list.csv')

ieu_list_eqtl <- ieu_list[grep('^eqtl-a',ieu_list$id),]
druggable_gene_list <- read.csv('D:/SH_desktop/Proteome_MR/Druggable_genes.csv')
druggable_gene_list <- druggable_gene_list$ensembl_gene_id

ieu_list_eqtl_ids <- ieu_list_eqtl$trait

matachs <- ieu_list_eqtl_ids %in% druggable_gene_list
ieu_list_eqtl_ids <- ieu_list_eqtl_ids[matachs]
ieu_list_eqtl_ids_ <- paste0('eqtl-a-',ieu_list_eqtl_ids)


eqtl_t2d_mr_results <- NULL
eqtl_exposure_dat <- NULL
t2d_outcome_dat <- NULL
eqtl_t2d_H_data <- NULL

for (ieu_list_eqtl_id in ieu_list_eqtl_ids_) {
  tmp_eqtl_exposure_dat <- extract_instruments(
    outcomes = ieu_list_eqtl_id,
    access_token = NULL)
  print(paste0('extract exposure snps=====',ieu_list_eqtl_id))
  
  if(!is.null(tmp_eqtl_exposure_dat) && any(tmp_eqtl_exposure_dat$SNP %in% outc$rsID)){
    t2d_outcome_dat <- format_data(
      type = "outcome",
      dat = outc,
      snps = tmp_eqtl_exposure_dat$SNP,
      header=T,
      snp_col = 'rsID',
      beta_col = 'Fixed.effects_beta',
      se_col = 'Fixed.effects_SE',
      pval_col = 'Fixed.effects_p.value',
      effect_allele_col = 'effect_allele',
      other_allele_col = 'other_allele',
      chr_col = 'chromosome',
      pos_col = 'position',
      eaf_col = 'effect_allele_frequency'
    )
    print('format outcome')
    
      tmp_eqtl_t2d_H_data <- harmonise_data(
      exposure_dat = tmp_eqtl_exposure_dat, 
      outcome_dat = t2d_outcome_dat
    )
      if(nrow(tmp_eqtl_t2d_H_data) !=0){
        print('harmonise data')
        eqtl_t2d_H_data <- rbind(eqtl_t2d_H_data, tmp_eqtl_t2d_H_data)
    
        print('perform MR')
        tmp_eqtl_t2d_mr_results <- tmp_eqtl_t2d_H_data %>% mr_singlesnp(
          parameters = default_parameters(),
          single_method = "mr_wald_ratio",
          all_method = c('mr_ivw_mre', 'mr_egger_regression',"mr_weighted_median", "mr_weighted_mode"))
        if(!is.null(tmp_eqtl_t2d_mr_results)){
          eqtl_t2d_mr_results <- rbind(eqtl_t2d_mr_results,tmp_eqtl_t2d_mr_results)
        }
    }
  }
} 
write.csv(eqtl_t2d_H_data, file = "eqtl_t2d_H_data_0502.csv", row.names = F)
write.csv(eqtl_t2d_mr_results, file = "eqtl_t2d_mr_results_0502.csv", row.names = F)

eqtl_t2d_mr_results$outcome <- "T2D_EA"
# trimws() function is then used to remove any leading or trailing whitespace
eqtl_t2d_mr_results <- eqtl_t2d_mr_results %>%
  separate(exposure, into = c('exposure','exposure.id'), sep = "\\|\\|", extra = "merge") %>%
  mutate(exposure = trimws(exposure))

# 2024-08-13 extract the instruments for Ziyao
IEU_opengwas_list <- read.csv("/Volumes/Echo_SSD/SH_desktop/Proteome_MR/IEU_opengwas_list.csv")
eqtl <- IEU_opengwas_list[grep("^eqtl-a", IEU_opengwas_list$id), ]
Druggable_genes <- read.csv("/Volumes/Echo_SSD/SH_desktop/Proteome_MR/Druggable_genes.csv",header = TRUE,sep = ",")
Druggable_genes_filtered <- Druggable_genes %>%
  filter(ensembl_gene_id %in% intersect(eqtl$trait, Druggable_genes$ensembl_gene_id)) %>%
  mutate(id = paste0("eqtl-a", ensembl_gene_id))

ieu_list_eqtl_ids_ <- Druggable_genes_filtered$id

for (ieu_list_eqtl_id in ieu_list_eqtl_ids_) {
  tmp_eqtl_exposure_dat <- extract_instruments(
    outcomes = ieu_list_eqtl_id,
    access_token = NULL)
  print(paste0('extract exposure snps=====',ieu_list_eqtl_id))
  eqtl_exposure_dat <- rbind(eqtl_exposure_dat, tmp_eqtl_exposure_dat)
  }


require('biomaRt')
## gene ID converation===========
# 1)设置需要连接的biomart
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl = useMart(biomart = "ensembl",
                  dataset = "hsapiens_gene_ensembl")
# 设置参数，建立检索
#filters - input
#attributes - output
#values - 应用于FILTER的值，用于检索
goids = getBM(attributes = c('hgnc_symbol','ensembl_gene_id'), 
              filters = "ensembl_gene_id", 
              values = eqtl_t2d_mr_results$exposure, 
              mart = ensembl)
# merge
eqtl_t2d_mr_results <- eqtl_t2d_mr_results %>% left_join(goids, by=c('exposure' ="ensembl_gene_id"))

# add snp number
eqtl_t2d_mr_results <- eqtl_t2d_mr_results %>% dplyr::select(-id.outcome)
eqtl_t2d_mr_results$id.outcome <- "L9LByT"
eqtl_t2d_mr_results$outcome <- "T2D_EA"
eqtl_t2d_mr_results <- unique(eqtl_t2d_mr_results)
#
sig_eqtl_t2d_mr_results <- eqtl_t2d_mr_results %>% filter(
  SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
             'All - MR Egger',
             "All - Weighted median",
             "All - Weighted mode") | p < 0.05/4479)

sig_eqtl_t2d_mr_results <- sig_eqtl_t2d_mr_results[!is.na(sig_eqtl_t2d_mr_results$b),]


# ARIC add nsnp to the significant results (968-sig_aric_proteins)
sig_eqtl_genes <- unique(sig_eqtl_t2d_mr_results$exposure)

for (gene in sig_eqtl_genes) {
  tmp_eqtl_t2d <- eqtl_t2d_mr_results %>% 
    filter(exposure == gene) %>%
    filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
             SNP != "All - MR Egger" &
             SNP != "All - Weighted median" &
             SNP != "All - Weighted mode")
  print(paste0('ARIC step1: filter snp results=====', gene))
  
  sig_eqtl_t2d_mr_results$nSNP[sig_eqtl_t2d_mr_results$exposure == gene]<- dim(tmp_eqtl_t2d)[1]
  print(paste0('ARIC step 2: add snp numbers======', dim(tmp_eqtl_t2d)[1]))
}
# ARIC 根据SNP列的信息，重新赋值是IVW,Wald ratio,MR_egger方法，筛选结果
sig_eqtl_t2d_mr_results$mr_method <- ifelse(sig_eqtl_t2d_mr_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                              sig_eqtl_t2d_mr_results$SNP != "All - MR Egger" & 
                                              sig_eqtl_t2d_mr_results$SNP != "All - Weighted median" & 
                                              sig_eqtl_t2d_mr_results$SNP != "All - Weighted mode" &
                                              sig_eqtl_t2d_mr_results$nSNP == 1, "Wald_ratio", 
                                                                           ifelse(sig_eqtl_t2d_mr_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                                    sig_eqtl_t2d_mr_results$SNP != "All - MR Egger" & 
                                                                                    sig_eqtl_t2d_mr_results$SNP != "All - Weighted median" &
                                                                                    sig_eqtl_t2d_mr_results$SNP != "All - Weighted mode" &
                                                                                    sig_eqtl_t2d_mr_results$nSNP != 1, "Confused", 
                                                                                  ifelse(sig_eqtl_t2d_mr_results$SNP != "All - MR Egger" &
                                                                                           sig_eqtl_t2d_mr_results$SNP != "All - Weighted median" &
                                                                                           sig_eqtl_t2d_mr_results$SNP != "All - Weighted mode", "IVW(random_models)", "Sensitivity analysis")))

sig_eqtl_t2d_mr_results_update <- sig_eqtl_t2d_mr_results %>%
  filter(mr_method != "Confused")


# ARIC 筛选SIG的结果 0.05/ 1563 (968)(sig_eqtl_proteins) - 67-sig_aric_pqtl_t2d_singlesnp_results_clump sig-results (1 duplicated proteins)
sig_eqtl_proteins <- unique(sig_eqtl_t2d_mr_results_update$exposure)
sig_eqtl_t2d_mr_results_withevidence <- NULL

for (protein in sig_eqtl_proteins) {
  tmp_sig_eqtl_t2d_mr_results <- sig_eqtl_t2d_mr_results_update %>% filter(exposure == protein)
  tmp_sig_eqtl_t2d_mr_results <- tmp_sig_eqtl_t2d_mr_results[!duplicated(tmp_sig_eqtl_t2d_mr_results$SNP),]
  tmp_sig_eqtl_t2d_mr_results$evidence1 <- ifelse(tmp_sig_eqtl_t2d_mr_results$mr_method == "IVW(random_models)" & 
                                                    tmp_sig_eqtl_t2d_mr_results$p < 0.05/4479, "robust",
                                                                    ifelse(tmp_sig_eqtl_t2d_mr_results$mr_method == "IVW(random_models)" &
                                                                             tmp_sig_eqtl_t2d_mr_results$p < 0.05, "probable","no_evidence"))
  print(paste0("Step 1:", protein))
  condition1 <- tmp_sig_eqtl_t2d_mr_results$SNP %in% c("All - MR Egger","All - Weighted median","All - Weighted mode") & 
    tmp_sig_eqtl_t2d_mr_results$p < 0.05
  #condition2 <- tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$SNP == "All - Weighted median" & tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05
  #condition3 <- tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$SNP == "All - Weighted mode" & tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$p < 0.05
  tmp_sig_eqtl_t2d_mr_results$evidence2[condition1] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[condition2] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[condition3] <- "robust"
  #tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2[is.na(tmp_sig_aric_pqtl_t2d_singlesnp_results_clump$evidence2)] <- "non_evidence"
  #
  sig_eqtl_t2d_mr_results_withevidence <- rbind(sig_eqtl_t2d_mr_results_withevidence, tmp_sig_eqtl_t2d_mr_results)
}

sig_eqtl_t2d_mr_results_withevidence <- unique(sig_eqtl_t2d_mr_results_withevidence)

sig_eqtl_t2d_mr_results_withevidence_update <- NULL

for (protein in sig_eqtl_proteins) {
  tmp_sig_eqtl_t2d_mr_results_withevidence_update <- sig_eqtl_t2d_mr_results_withevidence %>% filter(exposure == protein)
  if(sum(is.na(tmp_sig_eqtl_t2d_mr_results_withevidence_update[tmp_sig_eqtl_t2d_mr_results_withevidence_update$mr_method != "IVW(random_models)",]$evidence2)) ==3){
    tmp_sig_eqtl_t2d_mr_results_withevidence_update$evidence3 <- "non_evidence"
  }
  if(sum(is.na(tmp_sig_eqtl_t2d_mr_results_withevidence_update[tmp_sig_eqtl_t2d_mr_results_withevidence_update$mr_method != "IVW(random_models)",]$evidence2)) %in% c(1,2)){
    tmp_sig_eqtl_t2d_mr_results_withevidence_update$evidence3 <- "robust"
  }
  if(sum(is.na(tmp_sig_eqtl_t2d_mr_results_withevidence_update[tmp_sig_eqtl_t2d_mr_results_withevidence_update$mr_method != "IVW(random_models)",]$evidence2)) ==0){
    tmp_sig_eqtl_t2d_mr_results_withevidence_update$evidence3 <- "no_sensitivity_test"
  }
  print(paste0("step1:", protein))
  sig_eqtl_t2d_mr_results_withevidence_update <- rbind(sig_eqtl_t2d_mr_results_withevidence_update,tmp_sig_eqtl_t2d_mr_results_withevidence_update)
}

sig_eqtl_t2d_mr_results_withevidence_update$evidence_combined <- ifelse(sig_eqtl_t2d_mr_results_withevidence_update$evidence1 == "robust" & 
                                                                          sig_eqtl_t2d_mr_results_withevidence_update$evidence3 == "robust", "robust",
                                                                                                 ifelse(sig_eqtl_t2d_mr_results_withevidence_update$evidence1 == "robust"&
                                                                                                          sig_eqtl_t2d_mr_results_withevidence_update$evidence3 == "non_evidence","suggestive",
                                                                                                        ifelse(sig_eqtl_t2d_mr_results_withevidence_update$evidence1 == "probable"&
                                                                                                                 sig_eqtl_t2d_mr_results_withevidence_update$evidence3 == "robust","probable",
                                                                                                               ifelse(sig_eqtl_t2d_mr_results_withevidence_update$evidence1 == "probable"&
                                                                                                                        sig_eqtl_t2d_mr_results_withevidence_update$evidence3 == "non_evidence","suggestive","non_evidence"))))

write.csv(sig_eqtl_t2d_mr_results_withevidence_update, file = '/Volumes/Echo_SSD/SH_desktop/Proteome_MR/pQTL-analysis/sig_eqtl_t2d_mr_results_withevidence_update_0507.csv', row.names = F)
sig_eqtl_t2d_mr_results_withevidence_update <- sig_eqtl_t2d_mr_results_withevidence_update %>% filter(evidence_combined %in% c('robust', 'probable'))

intersect(sig_aric_pqtl_t2d_singlesnp_results_clump_nopleiotropy$exposure,sig_eqtl_t2d_mr_results_withevidence_update$hgnc_symbol)

sig_eqtl_t2d_mr_results_withevidence_update_druggable_gene_list <- sig_eqtl_t2d_mr_results_withevidence_update %>% filter(exposure %in% druggable_gene_list)


###
# perform: coloc analysis (exposure/eqtl-gene - outcome/T2D)
###

# ARIC(EA) add nsnp to the significant results (965-p1_sig_proteins)
eQTL_genes <- unique(sig_eqtl_t2d_mr_results_withevidence_update_druggable_gene_list$id.exposure)
eqtl_druggable_gene_list_exp <- NULL
eqtl_druggable_gene_list_H_data <- NULL

for (eQTL_gene in eQTL_genes) {
  print(paste0('step 0:=========', eQTL_gene))
  tmp_eqtl_druggable_gene_list_exp <- extract_instruments(outcomes=eQTL_gene, p1 = 1e-04, clump = F)
  if(nrow(tmp_eqtl_druggable_gene_list_exp!= 0)){
    tmp_eqtl_outcome_dat <- format_data(
      dat = outc,
      type = "outcome", 
      snps = tmp_eqtl_druggable_gene_list_exp$SNP,
      header = TRUE,
      phenotype_col = "phenotype",  
      snp_col = "rsID",
      beta_col = "Fixed.effects_beta",
      se_col = "Fixed.effects_SE",
      effect_allele_col = "effect_allele",                                                                                                                    
      other_allele_col = "other_allele",
      pval_col = "Fixed.effects_p.value",
      eaf_col = 'effect_allele_frequency'
    )
    print("step 1: extract outcome data")
    tmp_eqtl_druggable_gene_list_H_data <- harmonise_data(tmp_eqtl_druggable_gene_list_exp,
                                                          tmp_eqtl_outcome_dat)
    eqtl_druggable_gene_list_H_data <- rbind(eqtl_druggable_gene_list_H_data, tmp_eqtl_druggable_gene_list_H_data)
  }
}

write.csv(eqtl_druggable_gene_list_H_data, file = '/Volumes/Echo_SSD/SH_desktop/Proteome_MR/pQTL-analysis/eqtl_druggable_gene_list_H_data_0507.csv', row.names = F)

# 2024-05-08 共定位准备
# coloc analysis
eqtl_druggable_gene_list_H_data$eaf.exposure <- ifelse(eqtl_druggable_gene_list_H_data$eaf.exposure < 0.5, 
                                                       eqtl_druggable_gene_list_H_data$eaf.exposure,
                                                       1- eqtl_druggable_gene_list_H_data$eaf.exposure)

# 
library(stringr)
eqtl_druggable_gene_list_H_data$genes <- str_extract(eqtl_druggable_gene_list_H_data$id.exposure, "ENSG\\d+")

# run coloc analysis based on eqtl_druggable_gene_list_H_data
coloc_eqtl_genes <- unique(eqtl_druggable_gene_list_H_data$genes)
coloc_eqtl_t2d <- NULL

for (i in coloc_eqtl_genes){
  # dataset 1: pqtl summary; dataset 2: t2d
  tmp_coloc_eqtl_t2d <- eqtl_druggable_gene_list_H_data %>% filter(genes == i)
  
  print(paste0("step 1:", i))
  
  tmp_coloc_eqtl_t2d <- tmp_coloc_eqtl_t2d %>% na.omit()
  tmp_coloc_eqtl_t2d <- tmp_coloc_eqtl_t2d[tmp_coloc_eqtl_t2d$pval.exposure != 0 & !(duplicated(tmp_coloc_eqtl_t2d$SNP)),]
  
  print(paste0("step 2: remove duplicated SNPs and P missing values"))
  
  tmp <- coloc.abf(dataset2 = list(pvalues=tmp_coloc_eqtl_t2d$pval.outcome, 
                                   type = "cc", s = 0.16 , 
                                   N = 429191, snp=tmp_coloc_eqtl_t2d$SNP),
                   dataset1 = list(pvalues=tmp_coloc_eqtl_t2d$pval.exposure, type ="quant", 
                                   N = 31684, snp = tmp_coloc_eqtl_t2d$SNP), 
                   MAF = tmp_coloc_eqtl_t2d$eaf.exposure)
  print(paste0('step 1: coloc-analysis=====',i))
  tmp <- tmp$results
  if (nrow(tmp) != 0) {
    tmp$exposure <- i
    coloc_eqtl_t2d <- rbind(coloc_eqtl_t2d, tmp)
  }
}

write.csv(coloc_eqtl_t2d, file = "coloc_eqtl_t2d_0508.csv", row.names = F)

sig_coloc_eqtl_t2d <- coloc_eqtl_t2d %>% dplyr::filter(SNP.PP.H4>0.8)

# save file [MR info & coloc info] 16 genes pass the PP.H4 test 其中 [15 geness shared the same eqtl(IVs) with MR analysis]
# 10/29 proteins share the same pqtl(IVs) with MR anslysis
# sig_aric_pqtl_t2d_ivs_coloc (pqtl-analysis)/sig_coloc_eqtl_t2d (eqtl-analysis) 为 pass PP.H4 & pass unique pqtl(IVs) check 15 proteins(采用p<1e-04)
# rs7904519 don't overlap any IVs in MR analysis
write.csv(sig_coloc_eqtl_t2d, file = 'sig_coloc_eqtl_t2d_pph40.8_0508.csv', row.names = F)

# delete rs7904519 & merge
sig_coloc_eqtl_t2d <- sig_coloc_eqtl_t2d %>% left_join(goids, by=c('exposure' ="ensembl_gene_id"))


###
# perform TWAS: whole blood
###
YFS_twas_top <- fread('/Volumes/Echo_SSD/SH_desktop/Proteome_MR/t2d.YFS.BLOOD.TWAS.combined.top')
intersect(YFS_twas_top$ID, sig_eqtl_t2d_mr_results_withevidence_update_druggable_gene_list$hgnc_symbol)

###
# 2024-05-26 perform cox regression (SU_drug * breast/gastric cancer risk)
###
library("survival")
library("survminer")
# read data include covariables and outcomes
pheno_data <- read.csv('/Users/qianhe/Downloads/pheno_data_whi_baseline_outcome_0526.csv')
res.cox <- coxph(Surv(time_to_gastric, gastric_cancer) ~ diabetes_baseline, data = pheno_data)
summary(res.cox)
# diabetes baseline status 
# univariable coxph function applied to several covariates
covariates <- c('age','smk_now','logHDL','logTG','logBMI','logGLU','logHSCRP','logSBP','diabetes_baseline')
univ_formulas <- sapply(covariates, 
                       function(x) as.formula(paste('Surv(time_to_breast, breast_cancer)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = pheno_data)})

# extract values and make tables
univ_results <- lapply(univ_models, function(x) {
  x <-summary(x)
  p.value <- signif(x$wald['pvalue'], digits=2)
  wald.test <- signif(x$wald['test'], digits=2)
  beta <- signif(x$coef[1], digits=2) 
  HR <- signif(x$coef[2], digits=2)
  HR.confint.lower <- signif(x$conf.int[,'lower .95'], 2)
  HR.confint.upper <- signif(x$conf.int[,'upper .95'], 2)
  HR <- paste0(HR, "(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res <- c(beta, HR, wald.test, p.value)
  names(res) <- c("beta", "HR(95%CI)", "wald.test","p.value")
  return(res)
})
# summary the univariable cox model
res <- t(as.data.frame(univ_results, check.names=F))
res <- as.data.frame(res)

# multi-variables cox model 
res.cox <- coxph(Surv(time_to_gastric, gastric_cancer) ~ age+logBMI+logHSCRP+diabetes_baseline, 
                 data = pheno_data)
summary(res.cox)

# draw survival plot
diabetes_baseline_df <- with(pheno_data, data.frame(diabetes_baseline =c(0,1),
                                                    age = rep(mean(age, na.rm = T),2),
                                                    logBMI = rep(mean(logBMI, na.rm = T),2),
                                                    logHSCRP = rep(mean(logHSCRP, na.rm = T),2)))

fit <- survfit(res.cox, newdata = diabetes_baseline_df)
ggsurvplot(fit, data = pheno_data, conf.int = T, 
           legend.labs = c('Without Diabetes =0', 'With Diabetes=1'),
           ggtheme = theme_minimal())

# taking SU drugs
# univariable coxph function applied to several covariates
covariates1 <- c('age','smk_now','logHDL','logTG','logBMI','logGLU','logHSCRP','logSBP','SU_drug')
univ_formulas1 <- sapply(covariates1, 
                        function(x) as.formula(paste('Surv(time_to_breast, breast_cancer)~', x)))

univ_models1 <- lapply(univ_formulas1, function(x){coxph(x, 
                                                         data = subset(pheno_data, diabetes_baseline ==1))})

# extract values and make tables
univ_results1 <- lapply(univ_models1, function(x) {
  x <-summary(x)
  p.value <- signif(x$wald['pvalue'], digits=2)
  wald.test <- signif(x$wald['test'], digits=2)
  beta <- signif(x$coef[1], digits=2) 
  HR <- signif(x$coef[2], digits=2)
  HR.confint.lower <- signif(x$conf.int[,'lower .95'], 2)
  HR.confint.upper <- signif(x$conf.int[,'upper .95'], 2)
  HR <- paste0(HR, "(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res1 <- c(beta, HR, wald.test, p.value)
  names(res1) <- c("beta", "HR(95%CI)", "wald.test","p.value")
  return(res1)
})
# summary the univariable cox model
res1 <- t(as.data.frame(univ_results1, check.names=F))
res1 <- as.data.frame(res1)

# print tables for demographic information 
require(tableone)
# we need to convert the category variables as factor before making tables 
pheno_data$diabetes_baseline <- as.factor(pheno_data$diabetes_baseline)
pheno_data$SU_drug <- as.factor(pheno_data$SU_drug)

# breast cancer demographic table

table1 <- CreateTableOne(data = pheno_data, includeNA = T,
                         strata = "breast_cancer", test = T,
                         vars = c('age','smk_now','sbp','bmi',
                                  'ldl','hdl','hsCRP','diabetes_baseline',
                                  'SU_drug'))
print(table1, showAllLevels = F)
table1$CatTable
smmary(tableCatTable)

# gastric cancer demographic table
table2 <- CreateTableOne(data = pheno_data, includeNA = T,
                         strata = "gastric_cancer", test = T,
                         vars = c('age','smk_now','sbp','bmi',
                                  'ldl','hdl','hsCRP','diabetes_baseline',
                                  'SU_drug'))
print(table2, showAllLevels = F)


