graphics.off()
library("openxlsx")
library(pacman)
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
library(devtools)
# BiocManager::install("GenomicRanges")
# BiocManager::install("Rsamtools")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("VariantAnnotation")
# install_github("mrcieu/gwasvcf")
# install_github("mrcieu/gwasglue")
library(gwasglue)
library(gwasvcf)
suppressWarnings(suppressPackageStartupMessages({
  library(gwasvcf)
  library(VariantAnnotation)
  library(dplyr)
  library(magrittr)
}))
library(ieugwasr)
# install_github("phenoscanner/phenoscanner")
library(phenoscanner)
library(coloc)
library(colorspace)
library(hyprcoloc)
library(locuscomparer)

# 2025-03-27 add single cell analysis=========================================================

# Set input data path
sceqtl_path <- "/Volumes/Echo_SSD/QTL_summary_data/onek1k/sc_eQTL/"
sceqtl_data <- ".tsv"
sceqtl_esnps <- list.files(path = sceqtl_path, pattern = paste0("*","_esnp_table", sceqtl_data))

sc_name <- NULL

# For Loop of single-cell eQTL input (contain all the analysis steps)
suppressWarnings(for(i in sceqtl_esnps){
  # 1) read input data
  sceqtl_exp <- NULL
  sceqtl_exp <- suppressWarnings(fread(paste0(sceqtl_path, i), stringsAsFactors = FALSE, data.table = FALSE))
  #sceqtl_exp <- sceqtl_exp %>% head()
  sceqtl_exp$samplesize <- 982
  
  # read and format input data
  sceqtl_exp$se <- sqrt((1 - sceqtl_exp$SPEARMANS_RHO^2) / (sceqtl_exp$samplesize - 1))
  sceqtl_exp <- format_data(
    sceqtl_exp,
    type = "exposure",
    snp_col = "RSID",
    beta_col = "SPEARMANS_RHO",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "A2_FREQ_ONEK1K",
    pval_col = "P_VALUE",
    phenotype_col = "GENE"
  )
  print(paste0('step 0: format single cell data'))
  
  # Add the cancer phenotype information
  sceqtl_exp$phenotype.exposure <- sub("_.*", "", i)
  tmp_sceqtl_name <- sub("_.*", "", i)
  print(paste0("step 1: load single cell data======", tmp_sceqtl_name))
  
  # step 2) clump sceQTL input data
  sig_sceqtl_exp_clump <- NULL
  sc_genes <- unique(sceqtl_exp$exposure)
  
  for (sc_gene in sc_genes) {
    tmp_sceqtl_exp_clump <- sceqtl_exp %>%
      dplyr::filter(exposure == sc_gene)
    print(paste0("clump gene=============", sc_gene))
    #clump exp data
    tmp_sceqtl_exp_clump <- tryCatch({
      ld_clump(dat=tibble(rsid=tmp_sceqtl_exp_clump$SNP,
                          pval=tmp_sceqtl_exp_clump$pval.exposure,
                          id=tmp_sceqtl_exp_clump$id.exposure),
               clump_kb=5000,
               clump_r2 = 0.01,
               clump_p = 1,
               bfile = "/Volumes/Echo_SSD/QTL_summary_data/1kg.v3/EUR",
               plink_bin = "/Volumes/Echo_SSD/QTL_summary_data/1kg.v3/plink.exe")
    }, error = function(e) {
      message("Ignoring error: ", e$message)
      return(NULL)
    })
    sig_sceqtl_exp_clump <- rbind(sig_sceqtl_exp_clump, tmp_sceqtl_exp_clump) 
  }
  # filter clumped SNPs
  sig_sceqtl_exp_clump <- sceqtl_exp %>% dplyr::filter(SNP %in% sig_sceqtl_exp_clump$rsid)
  print("Filter clumpped SNPs")
  
  # save exposure input
  file_name <- paste0('/Volumes/Echo_SSD/QTL_summary_data/',tmp_sceqtl_name,'_exp_clump.csv')
  write.csv(sig_sceqtl_exp_clump, file = file_name, row.names = F)
  print("save exposure input")
  
  # 2) read outcome
  outc <- suppressWarnings(fread("/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/ORAL_OLK_CAVITY/summary_stats_release_finngen_R12_C3_ORALCAVITY_EXALLC.gz", 
                                 stringsAsFactors = FALSE, data.table = FALSE))
  
  oral_cavity_out <- format_data(dat = outc,
                                 type = "outcome",
                                 snps = sig_sceqtl_exp_clump$SNP,
                                 header = T,
                                 phenotype_col = "Oral_cavity_cancer",
                                 snp_col = "rsids",
                                 beta_col = "beta",
                                 se_col = "sebeta",
                                 effect_allele_col = "ref",
                                 other_allele_col = "alt",
                                 pval_col = "pval")
  print("step 2: Format outcome data")
  
  # 3) Perform MR (with clump)
  sceqtl_H_data <- NULL
  sceqtl_cancer_singlesnp_results <- NULL
  sig_sceqtl_gene_name <- NULL
  
  sig_sceqtl_genes <- unique(sig_sceqtl_exp_clump$exposure)
  
  for (sig_sceqtl_gene in sig_sceqtl_genes) {
    tmp_sig_sceqtl_exp_clump <- sig_sceqtl_exp_clump %>%
      dplyr::filter(exposure == sig_sceqtl_gene)
    sig_sceqtl_gene_name <- c(sig_sceqtl_gene_name, sig_sceqtl_gene)
    mm <- length(sig_sceqtl_gene_name)
    print(paste0("step1: filter single-cell genes======", sig_sceqtl_gene_name[mm]))
    
    # harmonise data
    tmp_sceqtl_H_data <- harmonise_data(
      exposure_dat = tmp_sig_sceqtl_exp_clump,
      outcome_dat = oral_cavity_out)
    sceqtl_H_data <- rbind(sceqtl_H_data,tmp_sceqtl_H_data)
    print("step2: Harmonize data")
    
    # perform MR
    if(!(is.data.frame(tmp_sceqtl_H_data) && nrow(tmp_sceqtl_H_data) == 0)){
      tmp_sceqtl_cancer_singlesnp_results <- tmp_sceqtl_H_data %>%
        mr_singlesnp(
          parameters = default_parameters(),
          single_method = "mr_wald_ratio",
          all_method = c('mr_ivw_mre', 'mr_egger_regression',"mr_weighted_median", "mr_weighted_mode"))
      tmp_sceqtl_cancer_singlesnp_results$genes <- sig_sceqtl_gene
    }
    print("step3: perform MR")
    
    sceqtl_cancer_singlesnp_results <- rbind(sceqtl_cancer_singlesnp_results,tmp_sceqtl_cancer_singlesnp_results)
    # save H_data
    file_name <- paste0('/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/',tmp_sceqtl_name, '_H_data.csv')
    write.csv(sceqtl_H_data, file = file_name, row.names = F)
    print("save H_data")
    
    # save MR results
    file_name <- paste0('/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/',tmp_sceqtl_name,'_singlesnp_results.csv')
    write.csv(sceqtl_cancer_singlesnp_results, file = file_name, row.names = F)
    print("save singlesnp_results")
  }
  
  # ============================================ PART 2 sc_eQTL quality control =========================
  
  # select significant results === sig genes number: length(unique(sig_sceqtl_genes))
  sig_sceqtl_cancer_singlesnp_results <- sceqtl_cancer_singlesnp_results %>% dplyr::filter(
    SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
               'All - MR Egger',
               "All - Weighted median",
               "All - Weighted mode") | p < 0.05)
  sig_sceqtl_cancer_singlesnp_results <- sig_sceqtl_cancer_singlesnp_results[!is.na(sig_sceqtl_cancer_singlesnp_results$b),]
  print("Quality control step1==============remove missing values")
  
  # add nsnp to the significant results === bin: 47 sig_sceqtl_genes
  sig_sceqtl_genes <- unique(sig_sceqtl_cancer_singlesnp_results$genes)
  
  for (gene in sig_sceqtl_genes) {
    tmp_sig_sceqtl_cancer <- sceqtl_cancer_singlesnp_results %>% 
      dplyr::filter(genes == gene) %>%
      dplyr::filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
               SNP != "All - MR Egger" &
               SNP != "All - Weighted median" &
               SNP != "All - Weighted mode")
    print(paste0('Quality control step2: filter snp results=====', gene))
    
    sig_sceqtl_cancer_singlesnp_results$nSNP[sig_sceqtl_cancer_singlesnp_results$genes == gene] <- dim(tmp_sig_sceqtl_cancer)[1]
    print(paste0('Quality control step 3: add snp numbers======', dim(tmp_sig_sceqtl_cancer)[1]))
  }
  # Reassign SNP methods based on SNP column info and filter results
  sig_sceqtl_cancer_singlesnp_results$mr_method <- ifelse(sig_sceqtl_cancer_singlesnp_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                            sig_sceqtl_cancer_singlesnp_results$SNP != "All - MR Egger" & 
                                                            sig_sceqtl_cancer_singlesnp_results$SNP != "All - Weighted median" & 
                                                            sig_sceqtl_cancer_singlesnp_results$SNP != "All - Weighted mode" &
                                                            sig_sceqtl_cancer_singlesnp_results$nSNP == 1, "Wald_ratio", 
                                                          ifelse(sig_sceqtl_cancer_singlesnp_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                   sig_sceqtl_cancer_singlesnp_results$SNP != "All - MR Egger" & 
                                                                   sig_sceqtl_cancer_singlesnp_results$SNP != "All - Weighted median" &
                                                                   sig_sceqtl_cancer_singlesnp_results$SNP != "All - Weighted mode" &
                                                                   sig_sceqtl_cancer_singlesnp_results$nSNP != 1, "Confused", 
                                                                 ifelse(sig_sceqtl_cancer_singlesnp_results$SNP != "All - MR Egger" &
                                                                          sig_sceqtl_cancer_singlesnp_results$SNP != "All - Weighted median" &
                                                                          sig_sceqtl_cancer_singlesnp_results$SNP != "All - Weighted mode", "IVW(random_models)", "Sensitivity analysis")))
  
  sig_sceqtl_cancer_singlesnp_results <- sig_sceqtl_cancer_singlesnp_results %>% 
    dplyr::filter(sig_sceqtl_cancer_singlesnp_results$mr_method != "Confused")
  
  # filter significant results === sig genes number: length(unique(sig_sceqtl_genes))
  sig_sceqtl_genes <- unique(sig_sceqtl_cancer_singlesnp_results$genes)
  sig_sceqtl_cancer_singlesnp_results_withevidence <- NULL
  
  for (gene in sig_sceqtl_genes) {
    tmp_sig_sceqtl_cancer_singlesnp_results <- sig_sceqtl_cancer_singlesnp_results %>% dplyr::filter(genes == gene)
    tmp_sig_sceqtl_cancer_singlesnp_results$evidence1 <- ifelse(tmp_sig_sceqtl_cancer_singlesnp_results$mr_method %in% c("IVW(random_models)", "Wald_ratio") & 
                                                                  tmp_sig_sceqtl_cancer_singlesnp_results$p < 0.05 / length(unique(sig_sceqtl_genes)), "robust",
                                                                ifelse(tmp_sig_sceqtl_cancer_singlesnp_results$mr_method %in% c("IVW(random_models)", "Wald_ratio") &
                                                                         tmp_sig_sceqtl_cancer_singlesnp_results$p < 0.05, "probable", "no_evidence"))
    print(paste0("Step 4 define sensitivity analysis:", gene))
    
    condition1 <- tmp_sig_sceqtl_cancer_singlesnp_results$SNP %in% c("All - MR Egger", "All - Weighted median", "All - Weighted mode") & 
      tmp_sig_sceqtl_cancer_singlesnp_results$p < 0.05
    tmp_sig_sceqtl_cancer_singlesnp_results$evidence2[condition1] <- "robust"
    sig_sceqtl_cancer_singlesnp_results_withevidence <- rbind(sig_sceqtl_cancer_singlesnp_results_withevidence, tmp_sig_sceqtl_cancer_singlesnp_results)
  }
  
  sig_sceqtl_cancer_singlesnp_results_withevidence <- unique(sig_sceqtl_cancer_singlesnp_results_withevidence)
  
  sig_sceqtl_cancer_singlesnp_results_withevidence_update <- NULL
  
  for (gene in sig_sceqtl_genes) {
    tmp_sig_sceqtl_cancer_singlesnp_results_withevidence <- sig_sceqtl_cancer_singlesnp_results_withevidence %>% 
      dplyr::filter(genes == gene)
    if(sum(is.na(tmp_sig_sceqtl_cancer_singlesnp_results_withevidence[!(tmp_sig_sceqtl_cancer_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) >= 3){
      tmp_sig_sceqtl_cancer_singlesnp_results_withevidence$evidence3 <- "non_evidence"
    }
    if(sum(is.na(tmp_sig_sceqtl_cancer_singlesnp_results_withevidence[!(tmp_sig_sceqtl_cancer_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) %in% c(1,2)){
      tmp_sig_sceqtl_cancer_singlesnp_results_withevidence$evidence3 <- "robust"
    }
    if(sum(is.na(tmp_sig_sceqtl_cancer_singlesnp_results_withevidence[!(tmp_sig_sceqtl_cancer_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) == 0){
      tmp_sig_sceqtl_cancer_singlesnp_results_withevidence$evidence3 <- "no_sensitivity_test"
    }
    print(paste0("step5 add sensitivity analysis:", gene))
    sig_sceqtl_cancer_singlesnp_results_withevidence_update <- rbind(sig_sceqtl_cancer_singlesnp_results_withevidence_update, tmp_sig_sceqtl_cancer_singlesnp_results_withevidence)
  }
  
  sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence_combined <- ifelse(sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust" & 
                                                                                        sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "robust", "robust",
                                                                                      ifelse(sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust" & 
                                                                                               sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "non_evidence", "suggestive",
                                                                                             ifelse(sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust" & 
                                                                                                      sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "no_sensitivity_test" & 
                                                                                                      sig_sceqtl_cancer_singlesnp_results_withevidence_update$mr_method == "Wald_ratio", "probable",
                                                                                                    ifelse(sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "probable" & 
                                                                                                             sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "robust", "probable",
                                                                                                           ifelse(sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "probable" & 
                                                                                                                    sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "non_evidence", "suggestive",
                                                                                                                  ifelse(sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "non_evidence" & 
                                                                                                                           sig_sceqtl_cancer_singlesnp_results_withevidence_update$evidence3 != "robust", "non_evidence", "suggestive"))))))
  
  print("Save quality control exp")
  file_name <- paste0("/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/sig_sceqtl_", tmp_sceqtl_name, "_mr_results_clump_withevidence.csv")
  write.csv(sig_sceqtl_cancer_singlesnp_results_withevidence_update, file = file_name, row.names = FALSE)
  
  sig_sceqtl_cancer_singlesnp_results <- sig_sceqtl_cancer_singlesnp_results_withevidence_update %>% 
    dplyr::filter(evidence_combined %in% c("robust", "probable", "suggestive"))
  
  print("Save quality control singlesnp_results")
  file_name <- paste0("/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/sig_sceqtl_", tmp_sceqtl_name, "_singlesnp_results_clump.csv")
  write.csv(sig_sceqtl_cancer_singlesnp_results, file = file_name, row.names = FALSE)
  
  
  # =========================== PART 3 Perform: coloc analysis sc-genes - oral cavity cancer ===================
  
  # filter MR results with proable and robust evidence
  tmp_sceqtl_exp <- sig_sceqtl_cancer_singlesnp_results %>% dplyr::select(exposure, genes, evidence_combined)
  # merge with the eQTL summary statics
  print("COLOC==========read single cell eqtl data")
  sceqtl_eqtl_exp <- suppressWarnings(fread(paste0(sceqtl_path, tmp_sceqtl_name, "_eqtl_table.tsv.gz"), 
                                            stringsAsFactors = FALSE, data.table = FALSE))
  tmp_sceqtl_exp <- tmp_sceqtl_exp %>% left_join(sceqtl_eqtl_exp, by = c('genes'='GENE'))
  
  # 合并sceqtl和GWAS的数据 [sig_tmp_sceqtl_exp_outc] 为[sceqtl - cancer]进行COLOC分析的输入数据
  sig_tmp_sceqtl_exp_outc <- tmp_sceqtl_exp %>% left_join(outc, by = c("RSID"="rsids"))
  
  print("COLOC==========save exp file")
  file_name <- paste0('/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/',tmp_sceqtl_name,'_exp_out_coloc.csv')
  write.csv(sig_tmp_sceqtl_exp_outc, file = file_name, row.names = F)
  
  # remove duplicated SNPs
  sig_tmp_sceqtl_exp_outc <- sig_tmp_sceqtl_exp_outc[order(sig_tmp_sceqtl_exp_outc[,'RSID'], sig_tmp_sceqtl_exp_outc[,'P_VALUE']),]
  sig_tmp_sceqtl_exp_outc <- sig_tmp_sceqtl_exp_outc[!duplicated(sig_tmp_sceqtl_exp_outc$RSID), ]
  
  # run coloc analysis based on sig_tmp_sceqtl_exp_outc
  sc_coloc_genes <- unique(sig_tmp_sceqtl_exp_outc$exposure)
  
  coloc_scqtl_cancer <- NULL
  
  for (i in sc_coloc_genes){
    # dataset 1: sceqtl summary || dataset 2: cancer || manually adjust s (case samplesize/total samplesize) and N (samplesize) 
    tmp_coloc_sceqtl_cancer <- sig_tmp_sceqtl_exp_outc %>% dplyr::filter(exposure == i)
    
    tmp_coloc_sceqtl_cancer <- tmp_coloc_sceqtl_cancer %>% na.omit()
    tmp_coloc_sceqtl_cancer <- tmp_coloc_sceqtl_cancer[tmp_coloc_sceqtl_cancer$P_VALUE != 0,]
    
    tmp <- coloc.abf(dataset2 = list(pvalues= tmp_coloc_sceqtl_cancer$pval, 
                                     type = "cc", s = 0.004261397, 
                                     N = 380363, snp= tmp_coloc_sceqtl_cancer$RSID),
                     dataset1 = list(pvalues= tmp_coloc_sceqtl_cancer$P_VALUE, type ="quant", 
                                     N = 982, snp = tmp_coloc_sceqtl_cancer$RSID), 
                     MAF = tmp_coloc_sceqtl_cancer$A2_FREQ_ONEK1K)
    print(paste0('COLOC step 1: coloc-analysis=====',i))
    tmp <- tmp$results
    if (nrow(tmp) != 0) {
      tmp$genes <- i
      coloc_scqtl_cancer <- rbind(coloc_scqtl_cancer, tmp)
    }
  }
  print("COLOC==========save colcoc_resules file")
  file_name <- paste0('/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/',tmp_sceqtl_name,'_results_coloc.csv')
  write.csv(coloc_scqtl_cancer, file = file_name, row.names = F)
  
  sig_coloc_scqtl_cancer <- coloc_scqtl_cancer %>% dplyr::filter(SNP.PP.H4>0.8)
  
  print("COLOC==========save colcoc_resules file (P>0.8)")
  file_name <- paste0('/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/sig_',tmp_sceqtl_name,'_results_coloc.csv')
  write.csv(sig_coloc_scqtl_cancer, file = file_name, row.names = F)
})

# =============================== 20250504 LDmatrix ======================================
# 1.将sig_results_coloc中的SNP挑选出来，选出它所在的gene
# 2.再在exp_clump中选择这个gene的IVs
# 3.用LDmatrix计算R2，结果存为LDinfo
# 4.选出LDinfo中除了1以外最大的R2，添加信息在sig_results_coloc中新的一列LD_R2
# 5.把sig_results_coloc_R2保存起来，把所有的sig_results_coloc_R2 merge在一起

# 先执行ORAL_LEUCOPLACIA；再oral cavity

# LD API: 3b284b35d031
library(dplyr)
library(LDlinkR)
library(httr)

sceqtl_path <- "C:/Single Cell MR/sc_eQTL/"
sceqtl_data <- ".tsv"
sceqtl_esnps <- list.files(path = sceqtl_path, pattern = paste0("*","_esnp_table", sceqtl_data))

suppressWarnings(for(i in sceqtl_esnps){
  tmp_sceqtl_name <- sub("_.*", "", i)
  
  file_name <- paste0("C:/Single Cell MR/SCRNA_MR/sceqtl_ORALCAVITY_EXALLC/sig_", tmp_sceqtl_name, "_results_coloc.csv")
  sig_results_coloc <- read.csv(file_name)
  
  sig_results_coloc$cell_type <- tmp_sceqtl_name
  file_name <- paste0("sig_", tmp_sceqtl_name, "_results_coloc_0405.csv")
  write.csv(sig_results_coloc, file = file_name, row.names = F)
  
  file_name <- paste0("C:/Single Cell MR/SCRNA_MR/sceqtl_ORALCAVITY_EXALLC/", tmp_sceqtl_name, "_exp_clump.csv")
  exp_clump <- read.csv(file_name)
  
  # 添加一列用于存储 R² 最大值
  sig_results_coloc$LD_R2 <- NA
  
  # 关闭 SSL 证书验证
  set_config(config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))
  
  for (i in 1:nrow(sig_results_coloc)) {
    coloc_SNP <- sig_results_coloc$snp[i]
    gene <- sig_results_coloc$genes[i]
    
    # 找出该 gene 的 IV SNPs
    IV_SNPs <- exp_clump %>%
      filter(exposure == gene) %>%
      pull(SNP)
    
    # 合并 SNP 列表
    calculate_SNPs <- unique(c(coloc_SNP, IV_SNPs))
    
    # 判断是否只有一个 SNP
    if (length(calculate_SNPs) == 1) {
      sig_results_coloc$LD_R2[i] <- NA
    } else {
      # 查询 LD 信息
      try({
        LDinfo <- LDmatrix(snps = calculate_SNPs,
                           pop = "EUR",
                           r2d = "r2",
                           token = "3b284b35d031",
                           file = FALSE)
        
        # 转换成矩阵并去掉对角线
        LD_mat <- LDinfo[, -1]
        LD_mat[LD_mat == "NA"] <- NA
        LD_mat_num <- suppressWarnings(apply(LD_mat, 2, as.numeric))
        diag(LD_mat_num) <- NA
        
        # 提取最大 R² 值
        max_R2 <- max(LD_mat_num, na.rm = TRUE)
        sig_results_coloc$LD_R2[i] <- max_R2
      }, silent = TRUE)
    }
    
    # 根据 LD_R2 分类
    sig_results_coloc$LD_confidence[i] <- ifelse(
      is.na(sig_results_coloc$LD_R2[i]), "oneSNP",
      ifelse(sig_results_coloc$LD_R2[i] > 0.8, "in_LD", "not_in_LD")
    )
    
    # 判断 coloc SNP 是否出现在 IV 列表中
    sig_results_coloc$confidence[i] <- ifelse(
      coloc_SNP %in% IV_SNPs, 
      "high_confidence", 
      "no"
    )
  }
  
  file_name <- paste0("sig_", tmp_sceqtl_name, "_results_coloc_confidence.csv")
  write.csv(sig_results_coloc, file = file_name, row.names = F)
  
})




# merge sig_results_coloc_confidence
sig_results_coloc_confidence <- NULL

suppressWarnings(for(i in sceqtl_esnps){
  tmp_sceqtl_name <- sub("_.*", "", i)
  
  file_name <- paste0("C:/Single Cell MR/SCRNA_MR/sceqtl_ORALCAVITY_EXALLC/confidence/sig_", tmp_sceqtl_name, "_results_coloc_confidence.csv")
  tmp_sig_results_coloc_confidence <- read.csv(file_name)
  
  sig_results_coloc_confidence <- rbind(sig_results_coloc_confidence, tmp_sig_results_coloc_confidence)
})

write.csv(sig_results_coloc_confidence, file = "ORALCAVITY_EXALLC_sig_results_coloc_confidence.csv", row.names = F)

