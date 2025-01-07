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

### Perform: primary Transcriptome - MR analysis [eQTL on 41 cancers]

# ============================= PART 1 Primary MR of eQTL (in druggable genes) and 41 cancers ====================================================

# Read the input list of genetic instruments for the analysis
ieu_list <- read.csv('D:/SH_desktop/Proteome_MR/IEU_opengwas_list.csv')

# Filter the list for eQTL data
ieu_list_eqtl <- ieu_list[grep('^eqtl-a',ieu_list$id),]

# load druggable genes
druggable_gene_list <- read.csv('D:/SH_desktop/Proteome_MR/Druggable_genes.csv')
druggable_gene_list <- druggable_gene_list$ensembl_gene_id

# Extract the protein identifiers for eQTL data
ieu_list_eqtl_ids <- ieu_list_eqtl$protein

# Extract the protein identifiers for eQTL data
ieu_list_eqtl_ids <- ieu_list_eqtl$protein

# Match eQTL proteins with druggable genes
matachs <- ieu_list_eqtl_ids %in% druggable_gene_list
ieu_list_eqtl_ids <- ieu_list_eqtl_ids[matachs]
ieu_list_eqtl_ids_ <- paste0('eqtl-a-',ieu_list_eqtl_ids)

# Extract associated genetic instruments
eqtl <- NULL

for (ieu_list_eqtl_id in ieu_list_eqtl_ids_){
  tmp_eqtl <- extract_instruments(
    outcomes = ieu_list_eqtl_id,
    access_token = NULL)
  print(paste0('extract exposure snps=====',ieu_list_eqtl_id))
  eqtl <- rbind(eqtl,tmp_eqtl)
} 
write.csv(eqtl, file = "ieu_extracted_druggable_genes.csv", row.names = F)

eqtl <- read.csv("C:/40cancers_analysis/ieu_extracted_druggable_genes.csv")

# filter druggable gene without hgnc_names
druggable_gene_list$ensembl_gene_id <- paste0('eqtl-a-',druggable_gene_list$ensembl_gene_id)
druggable_gene_list <- druggable_gene_list %>% filter(!hgnc_names=="noresult")

# Merge eQTL data with druggable gene list
eqtl <- eqtl %>% left_join(druggable_gene_list, by = c('id.exposure'='ensembl_gene_id'))

# Format the eQTL data
eqtl_exp <- format_data(
  eqtl,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  eaf_col = "eaf.exposure",
  phenotype_col = "hgnc_names",
)

eqtl_exp_protein <- unique(eqtl_exp$exposure)

# clump eQTL input data
sig_eqtl_exp_clump <- NULL

for (protein in eqtl_exp_protein) {
  tmp_eqtl_exp_clump <- eqtl_exp %>%
    filter(exposure == protein)
  #clump exp data
  tmp_eqtl_exp_clump <- tryCatch({
    ld_clump(dat=tibble(rsid=tmp_eqtl_exp_clump$SNP,
                        pval=tmp_eqtl_exp_clump$pval.exposure,
                        id=tmp_eqtl_exp_clump$id.exposure),
             clump_kb=5000,
             clump_r2 = 0.01,
             clump_p = 1,
             bfile = "C:/40cancers_analysis/1kg.v3/EUR",
             plink_bin = "C:/40cancers_analysis/1kg.v3/plink.exe")
  }, error = function(e) {
    message("Ignoring error: ", e$message)
    return(NULL)
  })
  sig_eqtl_exp_clump <- rbind(sig_eqtl_exp_clump, tmp_eqtl_exp_clump) 
}
sig_eqtl_exp_clump <- eqtl_exp %>% filter(eqtl_exp$SNP %in% sig_eqtl_exp_clump$rsid)
sig_eqtl_protein <- unique(sig_eqtl_exp_clump$exposure)

write.csv(sig_eqtl_exp_clump, file = "sig_eqtl_exp_clump.csv", row.names = F)

# Set outcome cancer GWAS data path
cancer_path <- "C:/40cancers_analysis/GWAS_cancer/"
cancer_data <- ".tsv"
cancers <- list.files(path = cancer_path, pattern = paste0("*",cancer_data))
cancer_name <- NULL

# For Loop of 41 cancers (contain all the analysis steps)
suppressWarnings(for(i in cancers){
  # 1) load outcome cancers data
  outc <- NULL
  outc <- suppressWarnings(fread(paste0(cancer_path,i),
                                        stringsAsFactors = F,
                                        data.table = F))
  
  if("effect_allele_frequency" %in% names(outc)) {
    outc <- outc %>% 
      dplyr::select(variant_id,
                    beta,
                    standard_error,
                    effect_allele,
                    other_allele,
                    p_value,
                    effect_allele_frequency)
  } else {
    outc <- outc %>% 
      dplyr::select(variant_id,
                    beta,
                    standard_error,
                    effect_allele,
                    other_allele,
                    p_value)
  }
  
  # Add the cancer phenotype info
  outc$phenotype <- gsub(cancer_data,"",i)
  tmp_cancer_name <- gsub(cancer_data,"",i)
  cancer_name <- c(cancer_name,tmp_cancer_name)

  # 2) Format outcome data 
  if("effect_allele_frequency" %in% names(outc)) {
    outc_dat <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_eqtl_exp_clump$SNP,
      header = TRUE,
      phenotype_col = "phenotype",
      snp_col = "variant_id",
      beta_col = "beta",
      se_col = "standard_error",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value",
      eaf_col = 'effect_allele_frequency' # Use the effect allele frequency column
    )
  } else {
    outc_dat <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_eqtl_exp_clump$SNP,
      header = TRUE,
      phenotype_col = "phenotype",
      snp_col = "variant_id",
      beta_col = "beta",
      se_col = "standard_error",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value"
    )
  }
  outc_dat$outcome <- tmp_cancer_name
  print("step0: format outcome")
  
  # 3) perform MR
  eqtl_H_data <- NULL
  eqtl_cancer_singlesnp_results <- NULL
  sig_eqtl_protein_name <- NULL

  for (protein in sig_eqtl_protein) {
    tmp_sig_eqtl_exp_clump <- sig_eqtl_exp_clump %>%
      filter(exposure == protein)
    sig_eqtl_protein_name <- c(sig_eqtl_protein_name, protein)
    mm <- length(sig_eqtl_protein_name)
    print(paste0("step1: filter proteins======", sig_eqtl_protein_name[mm]))
    
    # harmonise data
    tmp_eqtl_H_data <- harmonise_data(
      exposure_dat = tmp_sig_eqtl_exp_clump,
      outcome_dat = outc_dat)
    eqtl_H_data <- rbind(eqtl_H_data,tmp_eqtl_H_data)
    print("step2: Harmonize data")
    
    # perform MR
    if(!(is.data.frame(tmp_eqtl_H_data) && nrow(tmp_eqtl_H_data) == 0)){
      tmp_eqtl_cancer_singlesnp_results <- tmp_eqtl_H_data %>%
        mr_singlesnp(
          parameters = default_parameters(),
          single_method = "mr_wald_ratio",
          all_method = c('mr_ivw_mre', 'mr_egger_regression',"mr_weighted_median", "mr_weighted_mode"))
      tmp_eqtl_cancer_singlesnp_results$proteins <- protein
    }
    print("step3: perform MR")
    
    eqtl_cancer_singlesnp_results <- rbind(eqtl_cancer_singlesnp_results,tmp_eqtl_cancer_singlesnp_results)
  }
  
  file_name <- paste0('eqtl_', tmp_cancer_name, '_H_data.csv')
  write.csv(eqtl_H_data, file = file_name, row.names = F)

  file_name <- paste0('eqtl_',tmp_cancer_name,'_singlesnp_results.csv')
  write.csv(eqtl_cancer_singlesnp_results, file = file_name, row.names = F)

  # =========================================== PART 2 eQTL quality control =========================================
  
  eqtl_cancer_singlesnp_results$proceed[eqtl_cancer_singlesnp_results$exposure != eqtl_cancer_singlesnp_results$proteins] <- 0
  eqtl_cancer_singlesnp_results$proceed[eqtl_cancer_singlesnp_results$exposure == eqtl_cancer_singlesnp_results$proteins] <- 1

  # Prepare for secondary eQTL identification
  # input eQTL data
  cis_eqtl <- read.csv("C:/40cancers_analysis/Transcriptome-wide Association Analysis/eqtl_results.csv")
  cis_eqtl$ensembl_gene_id <- gsub("eqtl-a-", "", cis_eqtl$id.exposure)
  
  # add gene SYMBOL
  library(clusterProfiler)
  library(org.Hs.eg.db)
  ensembl_gene_id <- unique(cis_eqtl$ensembl_gene_id)
  ensembl_gene_id <- read.csv("C:/40cancers_analysis/Transcriptome-wide Association Analysis/ensemble_gene_id.csv")
  ensembl_gene_id <- bitr(ensembl_gene_id$ensembl_gene_id,
                        fromType = "ENSEMBL",
                        toType = "SYMBOL",
                        OrgDb = "org.Hs.eg.db")
  
  cis_eqtl <- cis_eqtl %>% left_join(ensembl_gene_id, by = c("ensembl_gene_id"="ENSEMBL"))
  cis_eqtl <- cis_eqtl %>% na.omit()

  # Select significant eQTL results
  nn <- length(unique(eqtl_cancer_singlesnp_results$proteins))

  sig_eqtl_cancer_singlesnp_results <- eqtl_cancer_singlesnp_results %>% 
    filter(SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
                      'All - MR Egger',
                      "All - Weighted median",
                      "All - Weighted mode") | p < 0.05/nn)

  sig_eqtl_cancer_singlesnp_results <- sig_eqtl_cancer_singlesnp_results[!is.na(sig_eqtl_cancer_singlesnp_results$b),]

  # eQTL add nsnp to the significant results
  sig_eqtl_proteins <- unique(sig_eqtl_cancer_singlesnp_results$proteins)

  for (protein in sig_eqtl_proteins) {
    tmp_eqtl_cancer_singlesnp_results <- eqtl_cancer_singlesnp_results %>%
      filter(proteins == protein) %>%
      filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" &
               SNP != "All - MR Egger" &
               SNP != "All - Weighted median" &
               SNP != "All - Weighted mode")
    print(paste0('eQTL step1: filter snp results=====', protein))

    sig_eqtl_cancer_singlesnp_results$nSNP[sig_eqtl_cancer_singlesnp_results$proteins == protein]<- dim(tmp_eqtl_cancer_singlesnp_results)[1]
    print(paste0('eQTL step 2: add snp numbers======', dim(tmp_eqtl_cancer_singlesnp_results)[1]))
  }
  
  # Update the method of MR based on SNP information (e.g., IVW, Wald ratio, MR Egger, etc.)
  sig_eqtl_cancer_singlesnp_results$mr_method <- ifelse(sig_eqtl_cancer_singlesnp_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" &
                                                          sig_eqtl_cancer_singlesnp_results$SNP != "All - MR Egger" &
                                                          sig_eqtl_cancer_singlesnp_results$SNP != "All - Weighted median" &
                                                          sig_eqtl_cancer_singlesnp_results$SNP != "All - Weighted mode" &
                                                          sig_eqtl_cancer_singlesnp_results$nSNP == 1, "Wald_ratio",
                                                        ifelse(sig_eqtl_cancer_singlesnp_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" &
                                                                 sig_eqtl_cancer_singlesnp_results$SNP != "All - MR Egger" &
                                                                 sig_eqtl_cancer_singlesnp_results$SNP != "All - Weighted median" &
                                                                 sig_eqtl_cancer_singlesnp_results$SNP != "All - Weighted mode" &
                                                                 sig_eqtl_cancer_singlesnp_results$nSNP != 1, "Confused",
                                                               ifelse(sig_eqtl_cancer_singlesnp_results$SNP != "All - MR Egger" &
                                                                        sig_eqtl_cancer_singlesnp_results$SNP != "All - Weighted median" &
                                                                        sig_eqtl_cancer_singlesnp_results$SNP != "All - Weighted mode", "IVW(random_models)", "Sensitivity analysis")))

  sig_eqtl_cancer_singlesnp_results <- sig_eqtl_cancer_singlesnp_results %>% filter(sig_eqtl_cancer_singlesnp_results$mr_method != "Confused")

  # Filter significant results
  sig_eqtl_proteins <- unique(sig_eqtl_cancer_singlesnp_results$proteins)
  mm <- length(unique(sig_eqtl_cancer_singlesnp_results$proteins))
  sig_eqtl_cancer_singlesnp_results_withevidence <- NULL

  for (protein in sig_eqtl_proteins) {
    tmp_sig_eqtl_cancer_singlesnp_results <- sig_eqtl_cancer_singlesnp_results %>% filter(proteins == protein)
    tmp_sig_eqtl_cancer_singlesnp_results$evidence1 <- ifelse(tmp_sig_eqtl_cancer_singlesnp_results$mr_method %in% c("IVW(random_models)", "Wald_ratio") &
                                                                tmp_sig_eqtl_cancer_singlesnp_results$p < 0.05/mm, "robust",
                                                              ifelse(tmp_sig_eqtl_cancer_singlesnp_results$mr_method %in% c("IVW(random_models)", "Wald_ratio") &
                                                                       tmp_sig_eqtl_cancer_singlesnp_results$p < 0.05, "probable","no_evidence"))
    print(paste0("Step 1:", protein))
    condition1 <- tmp_sig_eqtl_cancer_singlesnp_results$SNP %in% c("All - MR Egger","All - Weighted median","All - Weighted mode") &
      tmp_sig_eqtl_cancer_singlesnp_results$p < 0.05
    tmp_sig_eqtl_cancer_singlesnp_results$evidence2[condition1] <- "robust"
    sig_eqtl_cancer_singlesnp_results_withevidence <- rbind(sig_eqtl_cancer_singlesnp_results_withevidence, tmp_sig_eqtl_cancer_singlesnp_results)
  }

  # Further update evidence based on sensitivity tests and missing values
  sig_eqtl_cancer_singlesnp_results_withevidence_update <- NULL

  for (protein in sig_eqtl_proteins) {
    tmp_sig_eqtl_cancer_singlesnp_results_withevidence <- sig_eqtl_cancer_singlesnp_results_withevidence %>%
      filter(proteins == protein)

    # Count the number of missing values
    na_count <- sum(is.na(tmp_sig_eqtl_cancer_singlesnp_results_withevidence[!(tmp_sig_eqtl_cancer_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2))

    # Update evidence3 column based on NA counts
    if(na_count > 2){
      tmp_sig_eqtl_cancer_singlesnp_results_withevidence$evidence3 <- "non_evidence"
    } else if(na_count %in% c(1, 2)) {
      tmp_sig_eqtl_cancer_singlesnp_results_withevidence$evidence3 <- "robust"
    } else if(na_count == 0) {
      tmp_sig_eqtl_cancer_singlesnp_results_withevidence$evidence3 <- "no_sensitivity_test"
    }
    print(paste0("step1:", protein))

    # Ensure the update is initialized for the first protein
    if(protein == sig_eqtl_proteins[1]) {
      sig_eqtl_cancer_singlesnp_results_withevidence_update <- tmp_sig_eqtl_cancer_singlesnp_results_withevidence
    } else {
      sig_eqtl_cancer_singlesnp_results_withevidence_update <- rbind(sig_eqtl_cancer_singlesnp_results_withevidence_update, tmp_sig_eqtl_cancer_singlesnp_results_withevidence)
    }
  }

  # Update evidence_combined based on multiple conditions from evidence1 and evidence3
  sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence_combined <- ifelse(sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust" &
                                                                                      sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "robust", "robust",
                                                                                    ifelse(sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust"&
                                                                                             sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "non_evidence","suggestive",
                                                                                           ifelse(sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust"&
                                                                                                    sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "no_sensitivity_test" &
                                                                                                    sig_eqtl_cancer_singlesnp_results_withevidence_update$mr_method == "Wald_ratio","probable",
                                                                                                  ifelse(sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "probable"&
                                                                                                           sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "robust","probable",
                                                                                                         ifelse(sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "probable"&
                                                                                                                  sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "non_evidence","suggestive",
                                                                                                                ifelse(sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "non_evidence"&
                                                                                                                         sig_eqtl_cancer_singlesnp_results_withevidence_update$evidence3 != "robust","non_evidence", "suggestive"))))))


  # Filter significant results based on the combined evidence
  sig_eqtl_cancer_singlesnp_results <- sig_eqtl_cancer_singlesnp_results_withevidence_update %>%
    filter(evidence_combined %in% c("robust","probable"))
  
  file_name <- paste0("sig_eqtl_",tmp_cancer_name,"_singlesnp_results.csv")
  write.csv(sig_eqtl_cancer_singlesnp_results, file = file_name, row.names = F)
  
  #================================================= PART 4 secondary eQTL identification (association between eQTLs) =============================================
  
  # Filter significant IVs
  sig_eqtl_cancer_ivs <- eqtl_cancer_singlesnp_results %>% 
    filter(eqtl_cancer_singlesnp_results$exposure %in% sig_eqtl_cancer_singlesnp_results$exposure) %>%
    filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
             SNP != "All - MR Egger" & SNP != "All - Weighted median" & SNP != "All - Weighted mode")
  
  # Extract secondary eQTL data based on the significant SNPs
  secondary_eqtl <- cis_eqtl %>% filter(SNP %in% sig_eqtl_cancer_ivs$SNP & !(SYMBOL %in% sig_eqtl_cancer_ivs$exposure)) 
  
  # Filter secondary eQTL data based on p-value threshold for horizontal pleiotropy
  secondary_eqtl <- secondary_eqtl %>% filter(pval.exposure < 5e-08)
  secondary_eqtl <- data.frame(secondary_eqtl)
  
  # Group by gene symbol and filter for genes associated with more than one SNP
  secondary_eqtl <- secondary_eqtl %>%
    group_by(SYMBOL) %>%
    filter(n() >= 2)
  
  # Format the secondary eQTL data for exposure
  secondary_eqtl_exp <- format_data(
    secondary_eqtl,
    type = "exposure",
    snp_col = "SNP",
    beta_col = 'beta.exposure',
    se_col = 'se.exposure',
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    eaf_col = 'eaf.exposure',
    phenotype_col = "SYMBOL")
  
  secondary_eqtl_proteins <- unique(secondary_eqtl_exp$exposure)
  
  # Clump the secondary eQTL data
  sig_secondary_eqtl_exp_clump <- NULL
  
  for (protein in secondary_eqtl_proteins) {
    tmp_secondary_eqtl_exp_clump <- secondary_eqtl_exp %>% 
      filter(exposure == protein)
  
    tmp_secondary_eqtl_exp_clump <- tryCatch({
      ld_clump(dat=tibble(rsid= tmp_secondary_eqtl_exp_clump$SNP,
                          pval= tmp_secondary_eqtl_exp_clump$pval.exposure,
                          id= tmp_secondary_eqtl_exp_clump$id.exposure),
               clump_kb=5000,
               clump_r2 = 0.01,
               clump_p = 1,
               bfile = "C:/40cancers_analysis/1kg.v3/EUR",
               plink_bin = "C:/40cancers_analysis/1kg.v3/plink.exe")
    }, error = function(e) {
      message("Ignoring error: ", e$message)
      return(NULL)
    })
    sig_secondary_eqtl_exp_clump <- rbind(sig_secondary_eqtl_exp_clump, tmp_secondary_eqtl_exp_clump) 
  }
  sig_secondary_eqtl_exp_clump <- secondary_eqtl_exp %>% filter(secondary_eqtl_exp$SNP %in% sig_secondary_eqtl_exp_clump$rsid)
  
  # format outcome data
  if("effect_allele_frequency" %in% names(outc)) {
    # Use effect_allele_frequency if present
    secondary_outcome_dat <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_secondary_eqtl_exp_clump$SNP,
      header = TRUE,
      phenotype_col = "phenotype",
      snp_col = "variant_id",
      beta_col = "beta",
      se_col = "standard_error",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value",
      eaf_col = 'effect_allele_frequency'  # 使用effect_allele_frequency列
    )
  } else {
    # If effect_allele_frequency is not available
    secondary_outcome_dat <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_secondary_eqtl_exp_clump$SNP,
      header = TRUE,
      phenotype_col = "phenotype",
      snp_col = "variant_id",
      beta_col = "beta",
      se_col = "standard_error",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value"
    )
  }
  secondary_outcome_dat$outcome <- tmp_cancer_name
  print("step1: format outcome")
  
  # perform MR (with clumped data)
  secondary_H_data_clump <- NULL
  secondary_cancer_singlesnp_results_clump <- NULL
  secondary_protein_name <- NULL
  
  secondary_proteins <- unique(sig_secondary_eqtl_exp_clump$exposure)
  
  for (protein in secondary_proteins) {
    tmp_sig_secondary_eqtl_exp_clump <- sig_secondary_eqtl_exp_clump %>% 
      filter(exposure == protein)
    
    secondary_protein_name <- c(secondary_protein_name, protein)
    mm <- length(secondary_protein_name)
    print(paste0("step1: filter proteins======", secondary_protein_name[mm]))
    
    tmp_secondary_H_data_clump <- harmonise_data(exposure_dat = tmp_sig_secondary_eqtl_exp_clump,
                                                 outcome_dat = secondary_outcome_dat)
    secondary_H_data_clump <- rbind(secondary_H_data_clump, tmp_secondary_H_data_clump)
    print("step 2: Harnomize name")
    
    if (!(is.data.frame(tmp_secondary_H_data_clump) && nrow(tmp_secondary_H_data_clump) == 0)){
      tmp_secondary_cancer_singlesnp_results_clump <- tmp_secondary_H_data_clump %>% 
        mr_singlesnp(
          parameters = default_parameters(),
          single_method = "mr_wald_ratio",
          all_method = c('mr_ivw_mre', 'mr_egger_regression',"mr_weighted_median", "mr_weighted_mode"))
      tmp_secondary_cancer_singlesnp_results_clump$proteins <- protein
    }
    print("step 3: perform MR")
    
    secondary_cancer_singlesnp_results_clump <- rbind(secondary_cancer_singlesnp_results_clump,
                                                      tmp_secondary_cancer_singlesnp_results_clump)
  }
  
  file_name <- paste0('secondary_', tmp_cancer_name, '_H_data_clump.csv')
  write.csv(secondary_H_data_clump, file = file_name)

  file_name <- paste0('secondary_',tmp_cancer_name,'_singlesnp_results.csv')
  write_csv(secondary_cancer_singlesnp_results_clump,file = file_name)
 
  # =================================================== PART 5 Delete horizontal pleiotropy eQTL and run MR again ===================================
  
  # Filter out SNPs associated with the specified methods 
  secondary_eqtl_snp_list <- secondary_cancer_singlesnp_results_clump %>% filter(!(SNP %in% c("All - Inverse variance weighted (multiplicative random effects)",
                                                                                                                  "All - MR Egger","All - Weighted median",
                                                                                                                  "All - Weighted mode")))
  secondary_eqtl_snp_list <- unique(secondary_eqtl_snp_list$SNP)
  
  # Filter out SNPs associated with horizontal pleiotropy from the eQTL data
  sig_eqtl_exp_clump_nopleiotropy <- sig_eqtl_exp_clump %>% filter(!(SNP %in% secondary_eqtl_snp_list))
  sig_eqtl_exp_nopleiotropy <- eqtl_exp %>% filter(!(SNP %in% secondary_eqtl_snp_list))

  # 1) format outcome data
  if("effect_allele_frequency" %in% names(outc)) {
    outcome_dat_nopleiotropy <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_eqtl_exp_clump_nopleiotropy$SNP,
      header = TRUE,
      phenotype_col = "phenotype",
      snp_col = "variant_id",
      beta_col = "beta",
      se_col = "standard_error",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value",
      eaf_col = 'effect_allele_frequency'
      )
  } else {
    outcome_dat_nopleiotropy <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_eqtl_exp_clump_nopleiotropy$SNP,
      header = TRUE,
      phenotype_col = "phenotype",
      snp_col = "variant_id",
      beta_col = "beta",
      se_col = "standard_error",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value"
    )
  }
outcome_dat_nopleiotropy$outcome <- tmp_cancer_name

  # 2) perform MR (with clump)
  eqtl_H_data_clump_nopleiotropy <- NULL
  eqtl_cancer_singlesnp_results_clump_nopleiotropy <- NULL
  eqtl_protein_name_nopleiotropy <- NULL

  eqtl_proteins_nopleiotropy <- unique(sig_eqtl_exp_clump_nopleiotropy$exposure)

  for (protein in eqtl_proteins_nopleiotropy) {
    tmp_sig_eqtl_exp_clump_nopleiotropy <- sig_eqtl_exp_clump_nopleiotropy %>% 
      filter(exposure == protein)
  
    eqtl_protein_name_nopleiotropy <- c(eqtl_protein_name_nopleiotropy, protein)
    mm <- length(eqtl_protein_name_nopleiotropy)
    print(paste0("step1: filter proteins======", eqtl_protein_name_nopleiotropy[mm]))
  
    tmp_eqtl_H_data_clump_nopleiotropy <- harmonise_data(exposure_dat = tmp_sig_eqtl_exp_clump_nopleiotropy,
                                                         outcome_dat = outcome_dat_nopleiotropy)
    eqtl_H_data_clump_nopleiotropy <- rbind(eqtl_H_data_clump_nopleiotropy, tmp_eqtl_H_data_clump_nopleiotropy)
    print("step 2: Harnomize name")
  
    if (!(is.data.frame(tmp_eqtl_H_data_clump_nopleiotropy) && nrow(tmp_eqtl_H_data_clump_nopleiotropy) == 0)){
      tmp_eqtl_cancer_singlesnp_results_clump_nopleiotropy <- tmp_eqtl_H_data_clump_nopleiotropy %>% 
        mr_singlesnp(
          parameters = default_parameters(),
          single_method = "mr_wald_ratio",
          all_method = c('mr_ivw_mre', 'mr_egger_regression',"mr_weighted_median", "mr_weighted_mode"))
      tmp_eqtl_cancer_singlesnp_results_clump_nopleiotropy$proteins <- protein
    }
    print("step 3: perform MR")
  
    eqtl_cancer_singlesnp_results_clump_nopleiotropy <- rbind(eqtl_cancer_singlesnp_results_clump_nopleiotropy,
                                                                             tmp_eqtl_cancer_singlesnp_results_clump_nopleiotropy)
  }
  
  file_name <- paste0('eqtl_',tmp_cancer_name,'_H_data_clump_nopleiotropy(all).csv')
  write.csv(eqtl_H_data_clump_nopleiotropy, file = file_name, row.names = F)


  file_name <- paste0('eqtl_',tmp_cancer_name,'_mr_results_clump_nopleiotropy(all).csv')
  write.csv(eqtl_cancer_singlesnp_results_clump_nopleiotropy, file = file_name)

  # eQTL nopleiotropy quality control
  # eQTL data selecting significant results/最终剩 1030(aric_proteins)个蛋白==== 把上一个step中识别的robust/probable protein过滤到nopleiotropy中
  sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy <- eqtl_cancer_singlesnp_results_clump_nopleiotropy %>% filter(
    SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
               'All - MR Egger',
               "All - Weighted median",
               "All - Weighted mode") | p < 0.05)

  sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy <- sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy[!is.na(sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy$b),]
  
  sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy <- sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy %>% 
    filter(exposure %in% sig_eqtl_cancer_singlesnp_results_clump$exposure)
  
  file_name <- paste0("sig_eqtl_",tmp_cancer_name,"_singlesnp_results_clump_nopleiotropy.csv")
  write.csv(sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy, file = file_name, row.names = F)
  })

# =============================================== PART 6  overlap pQTL & eQTL MR result ====================================
overlap_nopleiotropy <- NULL
sig_pqtl_nopleiotropy <- NULL
sig_eqtl_nopleiotropy <- NULL

suppressWarnings(for (i in cancer_name) {
  print(i)
  
  # Load pQTL results
  path_sig_pqtl_nopleiotropy <- paste0("C:/40cancers_analysis/Transcriptome-wide Association Analysis/Overlap_identification/sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy/sig_aric_pqtl_",i,"_singlesnp_results_clump_nopleiotropy.csv")
  sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy <- read.csv(file = path_sig_pqtl_nopleiotropy)
  unique(sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$exposure)
  
  # Load eQTL results
  path_sig_eqtl_nopleiotropy <- paste0("C:/40cancers_analysis/Transcriptome-wide Association Analysis/2SMR(remove_secondary_eqtl)_1105/sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy/sig_eqtl_",i,"_singlesnp_results_clump_nopleiotropy.csv")
  sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy <- read.csv(file = path_sig_eqtl_nopleiotropy)
  unique(sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy$exposure)
  
  # overlap pQTL & eQTL nopleiotropy
  tmp_overlap_nopleiotropy <- intersect(sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$exposure,sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy$exposure)
  
  # filter for sig_pqtl_nopleiotropy
  tmp_sig_pqtl_nopleiotropy <- sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy %>% 
    filter(exposure %in% tmp_overlap_nopleiotropy)
  sig_pqtl_nopleiotropy <- rbind(sig_pqtl_nopleiotropy,tmp_sig_pqtl_nopleiotropy)
  
  # filter for sig_eqtl_nopleiotropy
  tmp_sig_eqtl_nopleiotropy <- sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy %>% 
    filter(exposure %in% tmp_overlap_nopleiotropy)
  sig_eqtl_nopleiotropy <- rbind(sig_eqtl_nopleiotropy,tmp_sig_eqtl_nopleiotropy)
  }
)

write.csv(sig_pqtl_nopleiotropy, file = "sig_pqtl_nopleiotropy.csv")
write.csv(sig_eqtl_nopleiotropy, file = "sig_eqtl_nopleiotropy.csv")

# ===================================================== PART 7 MR forest plot visualization ===========================

p2_sig_eqtl_singlesnp_results <- sig_eqtl_nopleiotropy %>% filter(
  SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
             'All - MR Egger') | p < 0.05)

p2_sig_eqtl_singlesnp_results <- p2_sig_eqtl_singlesnp_results[!is.na(p2_sig_eqtl_singlesnp_results$b),]

# Add SNP count for significant results
p2_sig_proteins <- unique(p2_sig_eqtl_singlesnp_results$proteins)

for (protein in p2_sig_proteins) {
  tmp_sig_eqtl_nopleiotropy <- sig_eqtl_nopleiotropy %>%
    filter(proteins == protein) %>%
    filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" &
             SNP != "All - MR Egger")
  print(paste0('step1: filter snp results=====', protein))
  
  p2_sig_eqtl_singlesnp_results$nSNP[p2_sig_eqtl_singlesnp_results$proteins == protein]<- dim(tmp_sig_eqtl_nopleiotropy)[1]
  print(paste0('step 2: add snp numbers======', dim(tmp_sig_eqtl_nopleiotropy)[1]))
}

# Assign MR method to each result and filter out confusing results
p2_sig_eqtl_singlesnp_results$mr_method <- ifelse(p2_sig_eqtl_singlesnp_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                    p2_sig_eqtl_singlesnp_results$SNP != "All - MR Egger" & p2_sig_eqtl_singlesnp_results$nSNP == 1, 
                                                  "Wald_ratio", ifelse(p2_sig_eqtl_singlesnp_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                         p2_sig_eqtl_singlesnp_results$SNP != "All - MR Egger" & p2_sig_eqtl_singlesnp_results$nSNP != 1, "Confused", 
                                                                       ifelse(p2_sig_eqtl_singlesnp_results$SNP == "All - MR Egger", "MR_egger", "IVW(random_models)")))

# Filter out "Confused" results, calculate odds ratios, and ensure only results with at least 1 SNP are included
p2_sig_eqtl_singlesnp_results <- p2_sig_eqtl_singlesnp_results %>% 
  filter(p2_sig_eqtl_singlesnp_results$mr_method != "Confused") %>% 
  generate_odds_ratios() %>% 
  filter(SNP != 'All - MR Egger') %>%
  filter(nSNP >=1)

# Format odds ratio labels (OR with 95% confidence intervals) for annotation
p2_sig_eqtl_singlesnp_results$orlab <- paste0("",
                                              ifelse(sprintf("%.2f",p2_sig_eqtl_singlesnp_results$or)<0.0051,
                                                     format(p2_sig_eqtl_singlesnp_results$or,scientific = TRUE,digits=3),
                                                     sprintf("%.2f",p2_sig_eqtl_singlesnp_results$or)),
                                              "(", 
                                              ifelse(sprintf("%.2f",p2_sig_eqtl_singlesnp_results$or_Oral_cavity_canceri95)<0.0051,
                                                     format(p2_sig_eqtl_singlesnp_results$or_Oral_cavity_canceri95,scientific = TRUE,digits=3),
                                                     sprintf("%.2f",p2_sig_eqtl_singlesnp_results$or_Oral_cavity_canceri95)),
                                              "-",
                                              ifelse(sprintf("%.2f",p2_sig_eqtl_singlesnp_results$or_uci95)<0.0051,
                                                     format(p2_sig_eqtl_singlesnp_results$or_uci95,scientific = TRUE,digits=3),
                                                     sprintf("%.2f",p2_sig_eqtl_singlesnp_results$or_uci95)),
                                              ")")
# Sort results by odds ratio
p2_sig_eqtl_singlesnp_results <- p2_sig_eqtl_singlesnp_results[order(p2_sig_eqtl_singlesnp_results[,'or']),]

# Set exposure factor for plotting
p2_sig_eqtl_singlesnp_results$exposure <- factor(p2_sig_eqtl_singlesnp_results$exposure, levels = unique(p2_sig_eqtl_singlesnp_results$exposure))

# Set coordinates for points in plot
p2_sig_eqtl_singlesnp_results$x <- 5.7
p2_sig_eqtl_singlesnp_results$y <- seq(1,28,1)

# Create annotation data frame
annotation<- NULL
annotation <- data.frame(matrix(ncol = 3, nrow = 0))
names(annotation) <- c("x", "y", "nSNP")

# Add relevant columns from p2_sig_eqtl_singlesnp_results to annotation
annotation0 <- p2_sig_eqtl_singlesnp_results[c('x','y','nSNP')]
annotation1 <- p2_sig_eqtl_singlesnp_results[c('x','y','outcome')] %>% dplyr::rename(nSNP=outcome)
annotation2 <- p2_sig_eqtl_singlesnp_results[c('x','y','exposure')] %>% dplyr::rename(nSNP=exposure)
annotation3 <- p2_sig_eqtl_singlesnp_results[c('x','y','orlab')] %>% dplyr::rename(nSNP=orlab)

# Merge the data frames to create the final annotation
annotation <- rbind(annotation0, annotation1, annotation2, annotation3)

# Adjust 'x' values for annotation positions
annotation[1:28, "x"] <- 0
annotation[29:56, "x"] <- -0.3
annotation[57:84, "x"] <- -0.8
annotation[85:112, "x"] <- 0.3

# Add additional annotation lines for labels
annotation <- rbind(annotation, c(0.3, 29, "OR(95%CI)"))
annotation <- rbind(annotation, c(0, 29, "nSNP"))
annotation <- rbind(annotation, c(-0.3, 29, "Cancer"))
annotation <- rbind(annotation, c(-0.8, 29, "Protein"))

# Ensure all columns are characters for consistency
annotation[] <- lapply(annotation, as.character)

# Convert x and y to numeric for proper positioning
annotation$x <- as.numeric(annotation$x)
annotation$y <- as.numeric(annotation$y)

# Create a PDF file to save the forest plot
pdf("eqtl_cancer_forest_plot_1106_2.pdf",
    width =10,
    height = 10
)

# Create forest plot using ggplot
library(colorspace)
p2_sig_eqtl_singlesnp_results %>% 
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
  ggtitle(expression("MR: eQTL"  %->% "Cancer(European Ancestry)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
