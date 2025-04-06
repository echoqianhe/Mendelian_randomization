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

suppressWarnings(for (tmp_cancer_name in cancer_name) {
  print(tmp_cancer_name)
  
  # Load pQTL results
  path_sig_pqtl_nopleiotropy <- paste0("C:/40cancers_analysis/Proteome-wide Association Analysis/sig_aric_pqtl_",tmp_cancer_name,"_singlesnp_results_clump_nopleiotropy_coloc.csv")
  sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy <- read.csv(file = path_sig_pqtl_nopleiotropy)
  unique(sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$exposure)
  
  # Load eQTL results
  path_sig_eqtl_nopleiotropy <- paste0("C:/40cancers_analysis/Transcriptome-wide Association Analysis/eQTL_optimization/sig_eqtl_",tmp_cancer_name,"_singlesnp_results_clump_nopleiotropy.csv")
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
p2_sig_eqtl_singlesnp_results$y <- seq(1,14,1)

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
annotation[1:14, "x"] <- 0
annotation[15:28, "x"] <- -0.3
annotation[29:42, "x"] <- -0.8
annotation[43:56, "x"] <- 0.3

# Add additional annotation lines for labels
annotation <- rbind(annotation, c(0.3, 14.4, "OR(95%CI)"))
annotation <- rbind(annotation, c(0, 14.4, "nSNP"))
annotation <- rbind(annotation, c(-0.3, 14.4, "Cancer"))
annotation <- rbind(annotation, c(-0.8, 14.4, "Protein"))

# Ensure all columns are characters for consistency
annotation[] <- lapply(annotation, as.character)

# Convert x and y to numeric for proper positioning
annotation$x <- as.numeric(annotation$x)
annotation$y <- as.numeric(annotation$y)

# Create a PDF file to save the forest plot
pdf("eqtl_cancer_forest_plot_1122.pdf",
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

# ================================== PART 8 Perform: MR - analysis [eQTL(gene) & risk factors] ==================

###
# risk factors (SMK,SBP,DBP,LDL,HDL,TC,Fasting Glu, Fasting Insu, BMI/Obesity) 
# select sig IVs for the previous genes

sig_eqtl_exp_clump <- read.csv("C:/40cancers_analysis/Transcriptome-wide Association Analysis/eQTL_optimization/sig_eqtl_exp_clump.csv")

suppressWarnings(for(tmp_cancer_name in cancer_name){
  path_secondary_cancer_singlesnp_results_clump <- paste0('C:/40cancers_analysis/Transcriptome-wide Association Analysis/eQTL_optimization/secondary_',tmp_cancer_name,'_singlesnp_results.csv')
  secondary_cancer_singlesnp_results_clump <- read.csv(path_secondary_cancer_singlesnp_results_clump)
  
  # Filter out SNPs associated with the specified methods 
  secondary_eqtl_snp_list <- secondary_cancer_singlesnp_results_clump %>% filter(!(SNP %in% c("All - Inverse variance weighted (multiplicative random effects)",
                                                                                              "All - MR Egger","All - Weighted median",
                                                                                              "All - Weighted mode")))
  secondary_eqtl_snp_list <- unique(secondary_eqtl_snp_list$SNP)
  
  # Filter out SNPs associated with horizontal pleiotropy from the eQTL data
  sig_eqtl_exp_clump_nopleiotropy <- sig_eqtl_exp_clump %>% filter(!(SNP %in% secondary_eqtl_snp_list))
  
  file_name <- paste0("sig_eqtl_",tmp_cancer_name,"_exp_clump_nopleiotropy.csv")
  write.csv(sig_eqtl_exp_clump_nopleiotropy, file = file_name, row.names = F)
  
  path_sig_eqtl_nopleiotropy <- paste0("C:/40cancers_analysis/Transcriptome-wide Association Analysis/eQTL_optimization/sig_eqtl_",tmp_cancer_name,"_singlesnp_results_clump_nopleiotropy.csv")
  sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy <- read.csv(file = path_sig_eqtl_nopleiotropy)
  
  sig_eqtl_riskfactor <- sig_eqtl_exp_clump_nopleiotropy %>% filter(SNP %in% sig_eqtl_cancer_singlesnp_results_clump_nopleiotropy$SNP)
  
  file_name <- paste0("sig_eqtl_",tmp_cancer_name,"_riskfactor.csv")
  write.csv(sig_eqtl_riskfactor, file = file_name)
})

# ======== The above is the preparation work =======

# ======== perform MR analysis [genes - riskfactors] =========

suppressWarnings(for(tmp_cancer_name in cancer_name){
  
  path_sig_eqtl_riskfactor <- paste0("./sig_eqtl_",tmp_cancer_name,"_riskfactor.csv")
  
  sig_eqtl_riskfactor <- read.csv(file = path_sig_eqtl_riskfactor)
  
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
      snps = sig_eqtl_riskfactor$SNP,
      outcomes = riskfactor,
      proxies = FALSE,
      access_token = NULL)
    print(paste0('extract outcome information====', riskfactor))
    
    tmp_riskfactor_H_data <- harmonise_data(
      exposure_dat = sig_eqtl_riskfactor, 
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
  
  file_name <- paste0("eqtl__",tmp_cancer_name,"_riskfactor_H_data_clump_nopleiotropy(all).csv")
  write.csv(riskfactor_H_data, file = file_name, row.names = F)
  
  file_name <- paste0("eqtl_",tmp_cancer_name,"riskfactors_nopleiotropy.csv")
  write.csv(riskfactor_mr_results, file = file_name, row.names = F)
  
})

# ======== riskfactors quality control ============

for(tmp_cancer_name in cancer_name){
  
  # load riskfactor_mr_results
  path_riskfactor_mr_results <- paste0("C:/Users/Kevin/Desktop/Table_csv/Echo2me/riskfactor_eqtl_analysis_results/eqtl_",tmp_cancer_name,"riskfactors_nopleiotropy.csv")
  riskfactor_mr_results <- read.csv(file = path_riskfactor_mr_results)
  
  ## selecting significant risk-factors
  sig_riskfactor_mr_results <- riskfactor_mr_results %>% filter(
    SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
               'All - MR Egger') | p < 0.05)
  
  sig_riskfactor_mr_results <- sig_riskfactor_mr_results[!is.na(sig_riskfactor_mr_results$b),]
  
  # eQTL add nsnp to the significant results
  riskfactor_proteins <- unique(sig_riskfactor_mr_results$exposure)
  
  for (riskfactor_protein in riskfactor_proteins) {
    print(paste0('step 0:=========', riskfactor_protein))
    
    for (riskfactor in riskfactors) {
      tmp_sig_eqtl_riskfactor <- sig_riskfactor_mr_results %>% 
        filter(exposure == riskfactor_protein & id.outcome == riskfactor) %>%
        filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                 SNP != "All - MR Egger")
      print(paste0('step 1: filter snp results=====', riskfactor))
      
      if (nrow(tmp_sig_eqtl_riskfactor) != 0) {
        sig_riskfactor_mr_results$nSNP <- dim(tmp_sig_eqtl_riskfactor)[1]
        print(paste0('step 2: add snp numbers======', dim(tmp_sig_eqtl_riskfactor)[1]))
      }
    }
  }
  
  # reassign method (sig_riskfactor_mr_results contain MR results of all genes)
  sig_riskfactor_mr_results$mr_method <- ifelse(sig_riskfactor_mr_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                  sig_riskfactor_mr_results$SNP != "All - MR Egger" & sig_riskfactor_mr_results$nSNP == 1, 
                                                "Wald_ratio", ifelse(sig_riskfactor_mr_results$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                       sig_riskfactor_mr_results$SNP != "All - MR Egger" & sig_riskfactor_mr_results$nSNP != 1, "Confused", 
                                                                     ifelse(sig_riskfactor_mr_results$SNP == "All - MR Egger", "MR_egger", "IVW(random_models)")))
  sig_riskfactor_mr_results <- sig_riskfactor_mr_results %>% 
    filter(sig_riskfactor_mr_results$mr_method != "Confused") %>% 
    generate_odds_ratios() %>% 
    filter(SNP != 'All - MR Egger')
  sig_riskfactor_mr_results <- sig_riskfactor_mr_results[!is.na(sig_riskfactor_mr_results$b),]
  
  
  # select significant MR results of genes and 12 risk factors
  sig_riskfactor_mr_results <- sig_riskfactor_mr_results %>% filter(p < 0.05/12)
  
  # save file (MR-risk factors) genes pass P< 0.05/12 threshold
  file_name <- paste0("sig_eqtl_",tmp_cancer_name,"_riskfactors.csv")
  write.csv(sig_riskfactor_mr_results, file = file_name, row.names = F)
}

# ===========2025-03-02 prepration of supplementary tables=================
# eqtl IVs -- stable 2.1
eqtl_ivs <- read.csv('/Users/qianhe/Downloads/Supplementary Result_0301/stable2.1_eqtl.csv')

eqtl_ivs$Tissue <- gsub(" \\(.*\\)", "", eqtl_ivs$Tissue)

eqtl_ivs <- eqtl_ivs %>%
  dplyr::rename(Gene = Tissue)
# remove duplicated IVs
eqtl_ivs <- eqtl_ivs[!duplicated(eqtl_ivs),]

write.csv(eqtl_ivs, file = "/Users/qianhe/Downloads/Supplementary Result_0301/stable2.1_eqtl_update.csv", row.names = F)

# pqtl IVs -- stable 2.2
pqtl_ivs <- read.csv('/Users/qianhe/Downloads/Supplementary Result_0301/stable2.2_pqtl.csv')

pqtl_ivs$Tissue <- gsub(" \\(.*\\)", "", pqtl_ivs$Tissue)

pqtl_ivs <- pqtl_ivs %>%
  dplyr::rename(Proteins = Tissue)
# remove duplicated IVs
pqtl_ivs <- pqtl_ivs[!duplicated(pqtl_ivs),]

write.csv(pqtl_ivs, file = "/Users/qianhe/Downloads/Supplementary Result_0301/stable2.2_pqtl_update.csv", row.names = F)


# eqtl & MR -- stable 3.1
eqtl_mr <- read.csv("/Users/qianhe/Downloads/Supplementary Result_0301/stable3.1_eqtl_mr.csv")

eqtl_mr$Outcome <- gsub("_", " ", eqtl_mr$Outcome)

eqtl_mr$OR <- round(eqtl_mr$OR, digits = 2)
eqtl_mr$X95.LCI <- round(eqtl_mr$X95.LCI, digits = 2)
eqtl_mr$X95.UCI <- round(eqtl_mr$X95.UCI, digits = 2)

eqtl_mr <- eqtl_mr %>% filter(P.value < 0.05)

write.csv(eqtl_mr, file = "/Users/qianhe/Downloads/Supplementary Result_0301/stable3.1_eqtl_mr_update.csv", row.names = F)


# pqtl & MR -- stable 3.2
pqtl_mr <- read.csv("/Users/qianhe/Downloads/Supplementary Result_0301/stable3.2_pqtl_mr.csv")

pqtl_mr$Outcome <- gsub("_", " ", pqtl_mr$Outcome)

pqtl_mr$OR <- round(pqtl_mr$OR, digits = 2)
pqtl_mr$X95.LCI <- round(pqtl_mr$X95.LCI, digits = 2)
pqtl_mr$X95.UCI <- round(pqtl_mr$X95.UCI, digits = 2)

pqtl_mr <- pqtl_mr %>% filter(P.value < 0.05)

write.csv(pqtl_mr, file = "/Users/qianhe/Downloads/Supplementary Result_0301/stable3.2_pqtl_mr_update.csv", row.names = F)

# identify overlapped results
overlapped_eqtl_pqtl <- eqtl_mr %>% dplyr::inner_join(pqtl_mr, 
                                                      by=c('Exposure','Outcome'), 
                                                      suffix=c('_eqtl',"_pqtl"))
# 2025-03-27 add single cell analysis=========================================================

# Set input data path
sceqtl_path <- "/Volumes/Echo_SSD/QTL_summary_data/onek1k/sc_eQTL/"
sceqtl_data <- ".tsv.gz"
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
  print(paste0('step 1: format single cell data'))
  # Add the cancer phenotype information
  sceqtl_exp$phenotype.exposure <- sub("_.*", "", i)
  tmp_sceqtl_name <- sub("_.*", "", i)
  #sc_name <- c(sc_name, tmp_sceqtl_name)
  #nn <- length(sc_name)
  print(paste0("step1: load single cell data======",tmp_sceqtl_name))
  
  # step 2) clump sceQTL input data
  sig_sceqtl_exp_clump <- NULL
  sc_genes <- unique(sceqtl_exp$exposure)
  
  for (sc_gene in sc_genes) {
    tmp_sceqtl_exp_clump <- sceqtl_exp %>%
      filter(exposure == sc_gene)
    print(paste0("clump gene=============", sc_gene))
    #clump exp data
    tmp_sceqtl_exp_clump <- tryCatch({
      ld_clump(dat=tibble(rsid=tmp_sceqtl_exp_clump$SNP,
                          pval=tmp_sceqtl_exp_clump$pval.exposure,
                          id=tmp_sceqtl_exp_clump$id.exposure),
               clump_kb=5000,
               clump_r2 = 0.01,
               clump_p = 1,
               bfile = "/Volumes/Echo_SSD/QTL_summary_data/10K_LD_ref/1kg.v3/EUR",
               plink_bin = "/Volumes/Echo_SSD/QTL_summary_data/10K_LD_ref/1kg.v3/plink_mac_20241022/plink")
    }, error = function(e) {
      message("Ignoring error: ", e$message)
      return(NULL)
    })
    sig_sceqtl_exp_clump <- rbind(sig_sceqtl_exp_clump, tmp_sceqtl_exp_clump) 
  }
  # filter clumped SNPs
  sig_sceqtl_exp_clump <- sceqtl_exp %>% filter(SNP %in% sig_sceqtl_exp_clump$rsid)
  print("Filter clumpped SNPs")
  # save exposure input
  file_name <- paste0('/Volumes/Echo_SSD/SH_desktop/SCRNA_MR/',tmp_sceqtl_name,'_exp_clump.csv')
  write.csv(sig_sceqtl_exp_clump, file = file_name, row.names = F)
  print("save exposure input")
  
  # 2) read outcome
  skin_melanoma_out <- format_data(dat = outc,
                                   type = "outcome",
                                   snps = sig_sceqtl_exp_clump$SNP,
                                   header = T,
                                   phenotype_col = "skin_melanona",
                                   snp_col = "variant_id",
                                   beta_col = "beta",
                                   se_col = "standard_error",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",
                                   pval_col = "p_value")
  print("step 2: Format outcome data")
  
  # 3) Perform MR (with clump)
  sceqtl_H_data <- NULL
  sceqtl_cancer_singlesnp_results <- NULL
  sig_sceqtl_gene_name <- NULL
  
  sig_sceqtl_genes <- unique(sig_sceqtl_exp_clump$exposure)
  
  for (sig_sceqtl_gene in sig_sceqtl_genes) {
    tmp_sig_sceqtl_exp_clump <- sig_sceqtl_exp_clump %>%
      filter(exposure == sig_sceqtl_gene)
    sig_sceqtl_gene_name <- c(sig_sceqtl_gene_name, sig_sceqtl_gene)
    mm <- length(sig_sceqtl_gene_name)
    print(paste0("step1: filter single-cell genes======", sig_sceqtl_gene_name[mm]))
    
    # harmonise data
    tmp_sceqtl_H_data <- harmonise_data(
      exposure_dat = tmp_sig_sceqtl_exp_clump,
      outcome_dat = skin_melanoma_out)
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
})




  
# =================================== 2025-03-02 preparation for the forest plot ======================================
# eqtl - MR Format odds ratio labels (OR with 95% confidence intervals) for annotation
p2_eqtl_mr <- eqtl_mr
p2_eqtl_mr$orlab <- paste0("",ifelse(sprintf("%.2f",p2_eqtl_mr$OR)<0.0051,
                                                     format(p2_eqtl_mr$OR,scientific = TRUE,digits=3),
                                                     sprintf("%.2f",p2_eqtl_mr$OR)),
                                              "(", 
                                              ifelse(sprintf("%.2f",p2_eqtl_mr$X95.LCI)<0.0051,
                                                     format(p2_eqtl_mr$X95.LCI,scientific = TRUE,digits=3),
                                                     sprintf("%.2f",p2_eqtl_mr$X95.LCI)),
                                              "-",
                                              ifelse(sprintf("%.2f",p2_eqtl_mr$X95.UCI)<0.0051,
                                                     format(p2_eqtl_mr$X95.UCI,scientific = TRUE,digits=3),
                                                     sprintf("%.2f",p2_eqtl_mr$X95.UCI)),
                                              ")")
# Sort results by odds ratio
p2_eqtl_mr <- p2_eqtl_mr[order(p2_eqtl_mr[,'OR']),]

# Set exposure factor for plotting
p2_eqtl_mr$Exposure <- factor(p2_eqtl_mr$Exposure, levels = unique(p2_eqtl_mr$Exposure))

# Set coordinates for points in plot
p2_eqtl_mr$x <- 2.5
p2_eqtl_mr$y <- seq(1,8,1)

# Create annotation data frame
annotation<- NULL
annotation <- data.frame(matrix(ncol = 3, nrow = 0))
names(annotation) <- c("x", "y", "nSNP")

# Add relevant columns from eqtl_mr to annotation
annotation0 <- p2_eqtl_mr[c('x','y','nSNP')]
annotation1 <- p2_eqtl_mr[c('x','y','Outcome')] %>% dplyr::rename(nSNP=Outcome)
annotation2 <- p2_eqtl_mr[c('x','y','Exposure')] %>% dplyr::rename(nSNP=Exposure)
annotation3 <- p2_eqtl_mr[c('x','y','orlab')] %>% dplyr::rename(nSNP=orlab)

# Merge the data frames to create the final annotation
annotation <- rbind(annotation0, annotation1, annotation2, annotation3)

# Adjust 'x' values for annotation positions
annotation[1:8, "x"] <- -0.5
annotation[9:16, "x"] <- -1.2
annotation[17:24, "x"] <- -2.2
annotation[25:32, "x"] <- 0.1

# Add additional annotation lines for labels
annotation <- rbind(annotation, c(0.1, 8.4, "OR(95%CI)"))
annotation <- rbind(annotation, c(-0.5, 8.4, "nSNP"))
annotation <- rbind(annotation, c(-1.2, 8.4, "Cancer"))
annotation <- rbind(annotation, c(-2.2, 8.4, "Protein"))

# Ensure all columns are characters for consistency
annotation[] <- lapply(annotation, as.character)

# Convert x and y to numeric for proper positioning
annotation$x <- as.numeric(annotation$x)
annotation$y <- as.numeric(annotation$y)

# Create a PDF file to save the forest plot
pdf("/Users/qianhe/Downloads/Supplementary Result_0301/eqtl_cancer_forest_plot_0302.pdf",
    width =10,
    height = 10
)

# Create forest plot using ggplot
library(colorspace)
p2_eqtl_mr %>% 
  ggplot(aes(OR, Exposure)) +
  geom_point(aes(col=OR,size = nSNP)) +
  scale_color_continuous_sequential(palette = "ag_Sunset", l1 = 20, c2 = 70, p1 = 1) +
  #scale_color_continuous_diverging(palette = "Tofino") +
  geom_errorbarh(aes(xmax=X95.UCI, xmin=X95.LCI, col = OR), height=0.2, size=1.0) +
  scale_x_continuous(limits = c(-2.5,1.5), breaks = seq(0,2.0,0.5)) +
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


# Format odds ratio labels (OR with 95% confidence intervals) for annotation
p2_eqtl_mr <- eqtl_mr
p2_eqtl_mr$orlab <- paste0("",ifelse(sprintf("%.2f",p2_eqtl_mr$OR)<0.0051,
                                     format(p2_eqtl_mr$OR,scientific = TRUE,digits=3),
                                     sprintf("%.2f",p2_eqtl_mr$OR)),
                           "(", 
                           ifelse(sprintf("%.2f",p2_eqtl_mr$X95.LCI)<0.0051,
                                  format(p2_eqtl_mr$X95.LCI,scientific = TRUE,digits=3),
                                  sprintf("%.2f",p2_eqtl_mr$X95.LCI)),
                           "-",
                           ifelse(sprintf("%.2f",p2_eqtl_mr$X95.UCI)<0.0051,
                                  format(p2_eqtl_mr$X95.UCI,scientific = TRUE,digits=3),
                                  sprintf("%.2f",p2_eqtl_mr$X95.UCI)),
                           ")")
# Sort results by odds ratio
p2_eqtl_mr <- p2_eqtl_mr[order(p2_eqtl_mr[,'OR']),]

# Set exposure factor for plotting
p2_eqtl_mr$Exposure <- factor(p2_eqtl_mr$Exposure, levels = unique(p2_eqtl_mr$Exposure))

# Set coordinates for points in plot
p2_eqtl_mr$x <- 2.5
p2_eqtl_mr$y <- seq(1,8,1)

# Create annotation data frame
annotation<- NULL
annotation <- data.frame(matrix(ncol = 3, nrow = 0))
names(annotation) <- c("x", "y", "nSNP")

# Add relevant columns from eqtl_mr to annotation
annotation0 <- p2_eqtl_mr[c('x','y','nSNP')]
annotation1 <- p2_eqtl_mr[c('x','y','Outcome')] %>% dplyr::rename(nSNP=Outcome)
annotation2 <- p2_eqtl_mr[c('x','y','Exposure')] %>% dplyr::rename(nSNP=Exposure)
annotation3 <- p2_eqtl_mr[c('x','y','orlab')] %>% dplyr::rename(nSNP=orlab)

# Merge the data frames to create the final annotation
annotation <- rbind(annotation0, annotation1, annotation2, annotation3)

# Adjust 'x' values for annotation positions
annotation[1:8, "x"] <- -0.5
annotation[9:16, "x"] <- -1.2
annotation[17:24, "x"] <- -2.2
annotation[25:32, "x"] <- 0.1

# Add additional annotation lines for labels
annotation <- rbind(annotation, c(0.1, 8.4, "OR(95%CI)"))
annotation <- rbind(annotation, c(-0.5, 8.4, "nSNP"))
annotation <- rbind(annotation, c(-1.2, 8.4, "Cancer"))
annotation <- rbind(annotation, c(-2.2, 8.4, "Protein"))

# Ensure all columns are characters for consistency
annotation[] <- lapply(annotation, as.character)

# Convert x and y to numeric for proper positioning
annotation$x <- as.numeric(annotation$x)
annotation$y <- as.numeric(annotation$y)

# Create a PDF file to save the forest plot
pdf("/Users/qianhe/Downloads/Supplementary Result_0301/eqtl_cancer_forest_plot_0302.pdf",
    width =10,
    height = 10
)
# Create forest plot using ggplot
library(colorspace)
p2_eqtl_mr %>% 
  ggplot(aes(OR, Exposure)) +
  geom_point(aes(col=OR,size = nSNP)) +
  scale_color_continuous_sequential(palette = "ag_Sunset", l1 = 20, c2 = 70, p1 = 1) +
  #scale_color_continuous_diverging(palette = "Tofino") +
  geom_errorbarh(aes(xmax=X95.UCI, xmin=X95.LCI, col = OR), height=0.2, size=1.0) +
  scale_x_continuous(limits = c(-2.5,1.5), breaks = seq(0,2.0,0.5)) +
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


# pqtl - MR Format odds ratio labels (OR with 95% confidence intervals) for annotation==============
p3_pqtl_mr <- pqtl_mr
p3_pqtl_mr$orlab <- paste0("",ifelse(sprintf("%.2f",p3_pqtl_mr$OR)<0.0051,
                                     format(p3_pqtl_mr$OR,scientific = TRUE,digits=3),
                                     sprintf("%.2f",p3_pqtl_mr$OR)),
                           "(", 
                           ifelse(sprintf("%.2f",p3_pqtl_mr$X95.LCI)<0.0051,
                                  format(p3_pqtl_mr$X95.LCI,scientific = TRUE,digits=3),
                                  sprintf("%.2f",p3_pqtl_mr$X95.LCI)),
                           "-",
                           ifelse(sprintf("%.2f",p3_pqtl_mr$X95.UCI)<0.0051,
                                  format(p3_pqtl_mr$X95.UCI,scientific = TRUE,digits=3),
                                  sprintf("%.2f",p3_pqtl_mr$X95.UCI)),
                           ")")
# Sort results by odds ratio
p3_pqtl_mr <- p3_pqtl_mr[order(p3_pqtl_mr[,'OR']),]

# Set exposure factor for plotting
p3_pqtl_mr$Exposure <- factor(p3_pqtl_mr$Exposure, levels = unique(p3_pqtl_mr$Exposure))

# Set coordinates for points in plot
p3_pqtl_mr$x <- 2.5
p3_pqtl_mr$y <- seq(1,10,1)

# Create annotation data frame
annotation<- NULL
annotation <- data.frame(matrix(ncol = 3, nrow = 0))
names(annotation) <- c("x", "y", "nSNP")

# Add relevant columns from eqtl_mr to annotation
annotation0 <- p3_pqtl_mr[c('x','y','nSNP')]
annotation1 <- p3_pqtl_mr[c('x','y','Outcome')] %>% dplyr::rename(nSNP=Outcome)
annotation2 <- p3_pqtl_mr[c('x','y','Exposure')] %>% dplyr::rename(nSNP=Exposure)
annotation3 <- p3_pqtl_mr[c('x','y','orlab')] %>% dplyr::rename(nSNP=orlab)

# Merge the data frames to create the final annotation
annotation <- rbind(annotation0, annotation1, annotation2, annotation3)

# Adjust 'x' values for annotation positions
annotation[1:10, "x"] <- -1.0
annotation[11:20, "x"] <- -2.0
annotation[21:30, "x"] <- -2.8
annotation[31:40, "x"] <- -0.5

# Add additional annotation lines for labels
annotation <- rbind(annotation, c(-0.5, 10.4, "OR(95%CI)"))
annotation <- rbind(annotation, c(-1.0, 10.4, "nSNP"))
annotation <- rbind(annotation, c(-2.0, 10.4, "Cancer"))
annotation <- rbind(annotation, c(-2.8, 10.4, "Protein"))

# Ensure all columns are characters for consistency
annotation[] <- lapply(annotation, as.character)

# Convert x and y to numeric for proper positioning
annotation$x <- as.numeric(annotation$x)
annotation$y <- as.numeric(annotation$y)

# Create a PDF file to save the forest plot
pdf("/Users/qianhe/Downloads/Supplementary Result_0301/pqtl_cancer_forest_plot_0302.pdf",
    width =10,
    height = 10
)
# Create forest plot using ggplot
library(colorspace)
p3_pqtl_mr %>% 
  ggplot(aes(OR, Exposure)) +
  geom_point(aes(col=OR,size = nSNP)) +
  scale_color_continuous_sequential(palette = "ag_Sunset", l1 = 20, c2 = 70, p1 = 1) +
  #scale_color_continuous_diverging(palette = "Tofino") +
  geom_errorbarh(aes(xmax=X95.UCI, xmin=X95.LCI, col = OR), height=0.2, size=1.0) +
  scale_x_continuous(limits = c(-3.0,2.0), breaks = seq(0,2.0,0.5)) +
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
  ggtitle(expression("MR: pQTL"  %->% "Cancer(European Ancestry)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#### 2025-03-06 zoom locus plot preparation ================
CFD_zoomlocus <- read.csv('/Users/qianhe/Downloads/coloc_aric_pqtl_Meta_Melanoma_skin_1e-02_CFD.csv')
gwas_cfd <- CFD_zoomlocus %>% dplyr::select(snp,pvalues.df2) %>%
  dplyr::rename(rsid=snp, pval = pvalues.df2)
pqtl_cfd <- CFD_zoomlocus %>% dplyr::select(snp,pvalues.df1) %>%
  dplyr::rename(rsid=snp, pval = pvalues.df1)

pdf(file = "/Users/qianhe/Downloads/Supplementary Result_0301/CFD_skin_melanoma_pqtl_zoomlocus_plot_EA_0306.pdf",
    width = 6,
    height = 6
)
locuscompare(in_fn1 = gwas_cfd,
             in_fn2 = pqtl_cfd,
             title1 = 'Skin melanoma GWAS',
             title2 = 'pQTL',
             snp = "rs117381356")
dev.off() 

##### 2025-03-18 bubble plot ====================
# perform meta-analysis for the mr results of eqtl and pqtl
library(metafor)
# mpo
data <- data.frame(
  study = c("eqtl", "pqtl"),
  or = c(0.89, 0.88),
  se = c(0.01441638,0.023207086)
)
# 使用IVW方法进行Meta分析
data$log_or <- log(data$or)
data$var <- data$se^2
# Perform IVW meta-analysis
# meta_result <- rma(yi = log_or, sei = se, method = "FE", data = data) # fixed-effect model

meta_result <- rma(yi = log_or, vi = var, data = data, method = "REML")

# Print the meta-analysis results
summary(meta_result)

# Create a forest plot
forest(meta_result)


# read overlapped exposure and outcome
p3_sig_overlapped_mr_results <- read.csv('/Users/qianhe/Downloads/Supplementary Result_0301/overlapped_mr_results_0318.csv')

pdf("/Users/qianhe/Downloads/Supplementary Result_0301/overlapped_proteins_cancers_bubble_plot(EA)_0318.pdf",
    width = 6,
    height = 6
)
p3_sig_overlapped_mr_results %>% 
  ggplot(aes(x = exposure, y = Outcome, color=meta_beta, na.rm = T)) +
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
  ggtitle(expression("MR: Genes(Proteins)"  %->% "Cancer susceptibility(European Ancestry)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
 
##### 2025-03-20 heatmap for TWMR and PWMR raw results========================
Atlas_sigtargets_forestplot_eqtl <- read.csv('/Users/qianhe/Downloads/Supplementary Result_0301/stable31_eqtl_mr.csv')

pdf("/Users/qianhe/Downloads/Supplementary Result_0301/Atlas_sigtargets_expression_heatmap_eqtl.pdf", # 文件名称
    width = 20,           # 宽
    height = 10)   
ggplot(Atlas_sigtargets_forestplot_eqtl, aes(y = Outcome, x = Exposure )) +
  labs(y = " ", x = "", fill = "Beta\n") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  #facet_grid(Atlas_sigtargets_forestplot_2$drugs~., scales = "free", space = "free", switch = "both") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = "right",
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 0, size=4),
        strip.text.y = element_text(angle = 90, hjust = 1, size=10),
        strip.background = element_blank()) +
  geom_tile(aes(fill = beta)) +
  #geom_point(aes(shape = factor(ex_null), color = factor(ex_null)), size = 1.8) +
  scale_shape_discrete(name  = " ",
                       breaks=c("1", "0"),
                       labels=c("P value < 0.00125", "0.00125<P value<0.05")) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       guide = "colourbar",
                       na.value = NA)
dev.off()


Atlas_sigtargets_forestplot_pqtl <- read.csv('/Users/qianhe/Downloads/Supplementary Result_0301/stable32_pqtl_mr.csv')

pdf("/Users/qianhe/Downloads/Supplementary Result_0301/Atlas_sigtargets_expression_heatmap_pqtl.pdf", # 文件名称
    width = 20,           # 宽
    height = 10)   
ggplot(Atlas_sigtargets_forestplot_pqtl, aes(y = Outcome, x = Exposure )) +
  labs(y = " ", x = "", fill = "Beta\n") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  #facet_grid(Atlas_sigtargets_forestplot_2$drugs~., scales = "free", space = "free", switch = "both") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = "right",
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 0, size=8),
        strip.text.y = element_text(angle = 90, hjust = 1, size=10),
        strip.background = element_blank()) +
  geom_tile(aes(fill = beta)) +
  #geom_point(aes(shape = factor(ex_null), color = factor(ex_null)), size = 1.8) +
  scale_shape_discrete(name  = " ",
                       breaks=c("1", "0"),
                       labels=c("P value < 0.00125", "0.00125<P value<0.05")) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       guide = "colourbar",
                       na.value = NA)
dev.off()

