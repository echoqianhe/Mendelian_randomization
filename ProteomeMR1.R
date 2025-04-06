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

### Perform: primary proteome - MR analysis [ARIC (EA -hg38) pQTL(proteins) on 41 cancers]

# ========================================= PART 1 Primary MR of pQTL and 41 cancers ====================================================

#  pQTL - input data
aric_path <- "C:/40cancers_analysis/EA/EA/"
aric_data <- ".PHENO1.glm.linear"
protein_names <- list.files(path=aric_path, pattern = paste0("*", aric_data))

tmp_aric_pqtl_exp <- NULL

# Process each protein file
for (i in protein_names){
  tmp <- suppressWarnings(fread(paste0(aric_path, i),
                                stringsAsFactors = FALSE,
                                data.table = FALSE))
  tmp$P <- as.numeric(tmp$P)
  print(paste0('step 1:', i))
  tmp <- tmp[tmp$P < 5e-08,] # Filter significant p-values
  
  if (nrow(tmp) != 0) {
    tmp$tissue <- gsub(aric_data, "", i)
    tmp_aric_pqtl_exp <- rbind(tmp_aric_pqtl_exp, tmp)
    print(paste0('step 2: combine the sig pQTL', gsub(aric_data, "", i)))
  }
}

# Load sequence data and rename columns
seq <- read.table("C:/40cancers_analysis/EA/seqid.txt", sep = "\t", header = TRUE)
seq <- seq %>% dplyr::rename(chr_hg38 = chromosome_name, transcription_start_site_hg38 = transcription_start_site)

# Merge pQTL data with sequence data
aric_pqtl_exp <- tmp_aric_pqtl_exp %>% left_join(seq, by = c('tissue' = 'seqid_in_sample'))

# Format input data for MR analysis
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
write.csv(aric_pqtl_exp, file = "aric_pqtl_exp.csv", row.names = FALSE)

aric_proteins <- unique(aric_pqtl_exp$exposure)

# Clump input data
sig_aric_pqtl_exp_clump <- NULL

for (protein in aric_proteins) {
  tmp_aric_pqtl_exp_clump <- aric_pqtl_exp %>%
    filter(exposure == protein)
  # Clump exp data
  tmp_aric_pqtl_exp_clump <- tryCatch({
    ld_clump(dat = tibble(rsid = tmp_aric_pqtl_exp_clump$SNP,
                          pval = tmp_aric_pqtl_exp_clump$pval.exposure,
                          id = tmp_aric_pqtl_exp_clump$id.exposure),
             clump_kb = 5000,
             clump_r2 = 0.01,
             clump_p = 1,
             bfile = "C:/40cancers_analysis/1kg.v3/EUR",
             plink_bin = "C:/40cancers_analysis/1kg.v3/plink.exe")
  }, error = function(e) {
    message("Error ignored: ", e$message)
    return(NULL)
  })
  sig_aric_pqtl_exp_clump <- rbind(sig_aric_pqtl_exp_clump, tmp_aric_pqtl_exp_clump)
}
sig_aric_pqtl_exp_clump <- aric_pqtl_exp %>% filter(aric_pqtl_exp$SNP %in% sig_aric_pqtl_exp_clump$rsid)

write.csv(sig_aric_pqtl_exp_clump, file = "sig_aric_pqtl_exp_clump.csv", row.names = FALSE)

# Set outcome data path
cancer_path <- "C:/40cancers_analysis/GWAS_cancer/"
cancer_data <- ".tsv"
cancers <- list.files(path = cancer_path, pattern = paste0("*", cancer_data))

# Supplementary table: s (case samplesize/total samplesize) and N (samplesize) of 41 cancers
supplementary_table <- read.csv(file = "C:/40cancers_analysis/Proteome-wide Association Analysis/supplementary_table.csv")

cancer_name <- NULL

# For Loop of 41 cancers (contain all the analysis steps)
suppressWarnings(for(i in cancers){
  # 1) load outcome cancers data
  outc <- NULL
  outc <- suppressWarnings(fread(paste0(cancer_path, i), stringsAsFactors = FALSE, data.table = FALSE))
  
  if("effect_allele_frequency" %in% names(outc)) {
    # If effect_allele_frequency column exists, select it along with other specified columns
    outc <- outc %>% 
      dplyr::select(variant_id,
                    beta,
                    standard_error,
                    effect_allele,
                    other_allele,
                    p_value,
                    effect_allele_frequency)
  } else {
    # If effect_allele_frequency column does not exist, select other specified columns
    outc <- outc %>% 
      dplyr::select(variant_id,
                    beta,
                    standard_error,
                    effect_allele,
                    other_allele,
                    p_value)
  }
  
  # Add the cancer phenotype information
  outc$phenotype <- gsub(cancer_data, "", i)
  tmp_cancer_name <- gsub(cancer_data, "", i)
  cancer_name <- c(cancer_name, tmp_cancer_name)
  
  # 2) Format outcome data
  if("effect_allele_frequency" %in% names(outc)) {
    # If effect_allele_frequency column exists, use it in format_data function
    aric_outcome_dat <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_aric_pqtl_exp_clump$SNP,
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
    # If effect_allele_frequency column does not exist, do not use it in format_data function
    aric_outcome_dat <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_aric_pqtl_exp_clump$SNP,
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
  aric_outcome_dat$outcome <- tmp_cancer_name
  
  # 3) Perform MR (with clump)
  aric_H_data_clump <- NULL
  aric_pqtl_cancer_singlesnp_results_clump <- NULL
  aric_protein_name <- NULL
  
  aric_proteins <- unique(sig_aric_pqtl_exp_clump$exposure)
  
  # Note: Not all proteins have results in outcome data, if not, cannot add PROTEINS
  for (protein in aric_proteins) {
    tmp_sig_aric_pQTL_exp_clump <- sig_aric_pqtl_exp_clump %>% 
      filter(exposure == protein)
    
    aric_protein_name <- c(aric_protein_name, protein)
    mm <- length(aric_protein_name)
    
    print(paste0("step1: filter proteins======", aric_protein_name[mm]))
    
    tmp_aric_H_data_clump <- harmonise_data(exposure_dat = tmp_sig_aric_pQTL_exp_clump,
                                            outcome_dat = aric_outcome_dat)
    aric_H_data_clump <- rbind(aric_H_data_clump, tmp_aric_H_data_clump)
    print("step 2: Harmonize name")
    
    if (!(is.data.frame(tmp_aric_H_data_clump) && nrow(tmp_aric_H_data_clump) == 0)){
      tmp_aric_pqtl_cancer_singlesnp_results_clump <- tmp_aric_H_data_clump %>% 
        mr_singlesnp(
          parameters = default_parameters(),
          single_method = "mr_wald_ratio",
          all_method = c('mr_ivw_mre', 'mr_egger_regression', "mr_weighted_median", "mr_weighted_mode"))
      tmp_aric_pqtl_cancer_singlesnp_results_clump$proteins <- protein
    }
    print("step 3: perform MR")
    
    aric_pqtl_cancer_singlesnp_results_clump <- rbind(aric_pqtl_cancer_singlesnp_results_clump,tmp_aric_pqtl_cancer_singlesnp_results_clump)
  }
  
  file_name <- paste0("aric_pqtl_", tmp_cancer_name, "_H_data_clump(all).csv")
  write.csv(aric_H_data_clump, file = file_name, row.names = FALSE)
  
  # ============================================ PART 2 ARIC(EA) quality control =========================
  
  aric_pqtl_cancer_singlesnp_results_clump$proceed[aric_pqtl_cancer_singlesnp_results_clump$exposure != aric_pqtl_cancer_singlesnp_results_clump$proteins] <- 0
  aric_pqtl_cancer_singlesnp_results_clump$proceed[aric_pqtl_cancer_singlesnp_results_clump$exposure == aric_pqtl_cancer_singlesnp_results_clump$proteins] <- 1
  
  file_name <- paste0("aric_pqtl_", tmp_cancer_name, "_mr_results_clump(all).csv")
  write.csv(aric_pqtl_cancer_singlesnp_results_clump, file = file_name, row.names = FALSE)
  
  # ARIC data selecting significant results/final remaining 1563 (aric_proteins) proteins ===
  sig_aric_pqtl_cancer_singlesnp_results_clump <- aric_pqtl_cancer_singlesnp_results_clump %>% filter(
    SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
               'All - MR Egger',
               "All - Weighted median",
               "All - Weighted mode") | p < 0.05 / 1563)
  sig_aric_pqtl_cancer_singlesnp_results_clump <- sig_aric_pqtl_cancer_singlesnp_results_clump[!is.na(sig_aric_pqtl_cancer_singlesnp_results_clump$b),]
  
  # ARIC add nsnp to the significant results ("1054"-sig_aric_proteins)
  sig_aric_proteins <- unique(sig_aric_pqtl_cancer_singlesnp_results_clump$proteins)
  
  for (protein in sig_aric_proteins) {
    tmp_sig_aric_pqtl_cancer <- aric_pqtl_cancer_singlesnp_results_clump %>% 
      filter(proteins == protein) %>%
      filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
               SNP != "All - MR Egger" &
               SNP != "All - Weighted median" &
               SNP != "All - Weighted mode")
    print(paste0('ARIC step1: filter snp results=====', protein))
    
    sig_aric_pqtl_cancer_singlesnp_results_clump$nSNP[sig_aric_pqtl_cancer_singlesnp_results_clump$proteins == protein] <- dim(tmp_sig_aric_pqtl_cancer)[1]
    print(paste0('ARIC step 2: add snp numbers======', dim(tmp_sig_aric_pqtl_cancer)[1]))
  }
  # ARIC Reassign SNP methods based on SNP column info and filter results
  sig_aric_pqtl_cancer_singlesnp_results_clump$mr_method <- ifelse(sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                     sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - MR Egger" & 
                                                                     sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - Weighted median" & 
                                                                     sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - Weighted mode" &
                                                                     sig_aric_pqtl_cancer_singlesnp_results_clump$nSNP == 1, "Wald_ratio", 
                                                                   ifelse(sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - Inverse variance weighted (multiplicative random effects)" & 
                                                                            sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - MR Egger" & 
                                                                            sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - Weighted median" &
                                                                            sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - Weighted mode" &
                                                                            sig_aric_pqtl_cancer_singlesnp_results_clump$nSNP != 1, "Confused", 
                                                                          ifelse(sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - MR Egger" &
                                                                                   sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - Weighted median" &
                                                                                   sig_aric_pqtl_cancer_singlesnp_results_clump$SNP != "All - Weighted mode", "IVW(random_models)", "Sensitivity analysis")))
  
  sig_aric_pqtl_cancer_singlesnp_results_clump <- sig_aric_pqtl_cancer_singlesnp_results_clump %>% filter(sig_aric_pqtl_cancer_singlesnp_results_clump$mr_method != "Confused")
  
  # ARIC filter significant results 
  sig_aric_proteins <- unique(sig_aric_pqtl_cancer_singlesnp_results_clump$proteins)
  sig_aric_pqtl_cancer_singlesnp_results_withevidence <- NULL
  
  for (protein in sig_aric_proteins) {
    tmp_sig_aric_pqtl_cancer_singlesnp_results_clump <- sig_aric_pqtl_cancer_singlesnp_results_clump %>% filter(proteins == protein)
    tmp_sig_aric_pqtl_cancer_singlesnp_results_clump$evidence1 <- ifelse(tmp_sig_aric_pqtl_cancer_singlesnp_results_clump$mr_method %in% c("IVW(random_models)", "Wald_ratio") & 
                                                                           tmp_sig_aric_pqtl_cancer_singlesnp_results_clump$p < 0.05 / 1054, "robust",
                                                                         ifelse(tmp_sig_aric_pqtl_cancer_singlesnp_results_clump$mr_method %in% c("IVW(random_models)", "Wald_ratio") &
                                                                                  tmp_sig_aric_pqtl_cancer_singlesnp_results_clump$p < 0.05, "probable", "no_evidence"))
    print(paste0("Step 1:", protein))
    
    condition1 <- tmp_sig_aric_pqtl_cancer_singlesnp_results_clump$SNP %in% c("All - MR Egger", "All - Weighted median", "All - Weighted mode") & 
      tmp_sig_aric_pqtl_cancer_singlesnp_results_clump$p < 0.05
    tmp_sig_aric_pqtl_cancer_singlesnp_results_clump$evidence2[condition1] <- "robust"
    sig_aric_pqtl_cancer_singlesnp_results_withevidence <- rbind(sig_aric_pqtl_cancer_singlesnp_results_withevidence, tmp_sig_aric_pqtl_cancer_singlesnp_results_clump)
  }
  
  sig_aric_pqtl_cancer_singlesnp_results_withevidence <- unique(sig_aric_pqtl_cancer_singlesnp_results_withevidence)
  
  sig_aric_pqtl_cancer_singlesnp_results_withevidence_update <- NULL
  
  for (protein in sig_aric_proteins) {
    tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence <- sig_aric_pqtl_cancer_singlesnp_results_withevidence %>% filter(proteins == protein)
    if(sum(is.na(tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence[!(tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) >= 3){
      tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence$evidence3 <- "non_evidence"
    }
    if(sum(is.na(tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence[!(tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) %in% c(1,2)){
      tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence$evidence3 <- "robust"
    }
    if(sum(is.na(tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence[!(tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence$mr_method %in% c("IVW(random_models)", "Wald_ratio")),]$evidence2)) == 0){
      tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence$evidence3 <- "no_sensitivity_test"
    }
    print(paste0("step1:", protein))
    sig_aric_pqtl_cancer_singlesnp_results_withevidence_update <- rbind(sig_aric_pqtl_cancer_singlesnp_results_withevidence_update, tmp_sig_aric_pqtl_cancer_singlesnp_results_withevidence)
  }
  
  sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence_combined <- ifelse(sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust" & 
                                                     sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "robust", "robust",
                                                   ifelse(sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust" & 
                                                            sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "non_evidence", "suggestive",
                                                          ifelse(sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "robust" & 
                                                                   sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "no_sensitivity_test" & 
                                                                   sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$mr_method == "Wald_ratio", "probable",
                                                                 ifelse(sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "probable" & 
                                                                          sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "robust", "probable",
                                                                        ifelse(sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "probable" & 
                                                                                 sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence3 == "non_evidence", "suggestive",
                                                                               ifelse(sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence1 == "non_evidence" & 
                                                                                        sig_aric_pqtl_cancer_singlesnp_results_withevidence_update$evidence3 != "robust", "non_evidence", "suggestive"))))))
  
  file_name <- paste0("sig_aric_pqtl_", tmp_cancer_name, "_mr_results_clump_withevidence.csv")
  write.csv(sig_aric_pqtl_cancer_singlesnp_results_withevidence_update, file = file_name, row.names = FALSE)
  
  sig_aric_pqtl_cancer_singlesnp_results_clump <- sig_aric_pqtl_cancer_singlesnp_results_withevidence_update %>% 
    filter(evidence_combined %in% c("robust", "probable"))
  
  file_name <- paste0("sig_aric_pqtl_", tmp_cancer_name, "_singlesnp_results_clump.csv")
  write.csv(sig_aric_pqtl_cancer_singlesnp_results_clump, file = file_name, row.names = FALSE)
  
  
  # =========================== PART 3 Perform: coloc-localization [ARIC(EA-hg38) significant MR proteins & cancer]  ===================
  
  # Only select significant protein GWAS files and filter robust and probable proteins
  tmp_aric_pqtl_exp <- sig_aric_pqtl_cancer_singlesnp_results_clump %>% dplyr::select(exposure, proteins, evidence_combined)
  
  # Delete duplicated probes' proteins based on QC principles
  duplicated_probes <- seq[duplicated(seq$entrezgenesymbol) == TRUE, ]
  duplicated_probes_protein <- intersect(duplicated_probes$entrezgenesymbol, sig_aric_pqtl_cancer_singlesnp_results_clump$exposure)
  tmp_aric_pqtl_exp <- tmp_aric_pqtl_exp %>% filter(!(exposure %in% duplicated_probes_protein))
  
  tmp_aric_pqtl_exp <- tmp_aric_pqtl_exp %>% left_join(seq, by = c('exposure' = 'entrezgenesymbol'))
  
  # Only read proteins GWAS in the following step/ tmp_aric_pqtl_exp_proteins is the vector include proteins names
  tmp_aric_pqtl_exp_proteins <- tmp_aric_pqtl_exp$seqid_in_sample
  tmp_aric_pqtl_exp_proteins <- paste(tmp_aric_pqtl_exp_proteins, ".PHENO1.glm.linear", sep = "")
  
  ## COLOC: this part of the code is for colocalization analysis
  tmp_aric_pqtl_exp <- NULL
  
  for (i in tmp_aric_pqtl_exp_proteins){
    tmp <- suppressWarnings(fread(paste0(aric_path, i),
                                  stringsAsFactors = FALSE,
                                  data.table = FALSE))
    tmp$P <- as.numeric(tmp$P)
    print(paste0('step 1:', i))
    tmp <- tmp[tmp$P < 1e-04,]
    
    if (nrow(tmp) != 0) {
      tmp$tissue <- gsub(aric_data, "", i)
      tmp_aric_pqtl_exp <- rbind(tmp_aric_pqtl_exp, tmp)
      print(paste0('step 2: combine the sig pQTL', gsub(aric_data, "", i)))
    }
  }
  
  # sig_tmp_aric_pqtl_exp is the PQTL (P<5e-04) for significant proteins used for COLOC analysis
  sig_tmp_aric_pqtl_exp <- tmp_aric_pqtl_exp %>% left_join(seq, by = c('tissue' = 'seqid_in_sample'))
  
  # Rename columns
  sig_tmp_aric_pqtl_exp <- sig_tmp_aric_pqtl_exp %>% dplyr::rename(CHR = "#CHROM", SNP = ID, BP = POS)
  sig_tmp_aric_pqtl_exp <- sig_tmp_aric_pqtl_exp %>% dplyr::select(SNP, CHR, BP, REF, ALT, A1, A1_FREQ, BETA, SE, P, tissue, uniprot_id, entrezgenesymbol)
  
  # Convert hg38 to hg37
  # Liftover sumstats_dt to hg37 coordinates for ARIC significant (MR-67) PQTL
  # sumstats_dt is hg37 coordinates for sig_tmp_aric_pqtl_exp
  library(MungeSumstats)
  sumstats_dt <- MungeSumstats::liftover(sumstats_dt = sig_tmp_aric_pqtl_exp,
                                         ref_genome = "hg38",
                                         convert_ref_genome = "hg19")
  knitr::kable(head(sumstats_dt))
  
  # Merge pqtl and GWAS data for COLOC analysis input
  sig_tmp_aric_pqtl_exp_outc <- sig_tmp_aric_pqtl_exp %>% left_join(outc, by = c("SNP" = "variant_id"))
  
  # Remove duplicated SNPs
  sig_tmp_aric_pqtl_exp_outc <- sig_tmp_aric_pqtl_exp_outc[order(sig_tmp_aric_pqtl_exp_outc[, 'SNP'], sig_tmp_aric_pqtl_exp_outc[, 'P']), ]
  sig_tmp_aric_pqtl_exp_outc <- sig_tmp_aric_pqtl_exp_outc[!duplicated(sig_tmp_aric_pqtl_exp_outc$SNP), ]
  
  file_name <- paste0("sig_tmp_aric_pqtl_exp_outc_", tmp_cancer_name, ".csv")
  write.csv(sig_tmp_aric_pqtl_exp_outc, file = file_name, row.names = FALSE)
  
  # Run coloc analysis based on sig_tmp_aric_pqtl_exp_outc
  coloc_aric_proteins <- unique(sig_tmp_aric_pqtl_exp_outc$entrezgenesymbol)
  
  coloc_aric_pqtl_cancer <- NULL
  
  # Get N and s values from the supplementary table
  N_value <- supplementary_table %>%
    filter(cancer_name == tmp_cancer_name) %>%
    dplyr::select(N) %>%
    pull()
  
  s_value <- supplementary_table %>%
    filter(cancer_name == tmp_cancer_name) %>%
    dplyr::select(s) %>%
    pull()
  
  # Check if the effect_allele_frequency column exists in sig_tmp_aric_pqtl_exp_outc
  if ("effect_allele_frequency" %in% names(sig_tmp_aric_pqtl_exp_outc)) {
    # Check if the effect_allele_frequency column is all NA
    all_na <- sig_tmp_aric_pqtl_exp_outc %>% 
      summarise(all_na = all(is.na(effect_allele_frequency))) %>% 
      pull(all_na)
    
    # If all NA, delete the column
    if (all_na) {
      sig_tmp_aric_pqtl_exp_outc <- sig_tmp_aric_pqtl_exp_outc %>% 
        dplyr::select(-effect_allele_frequency)
    }
  }
  
  # Perform coloc analysis for each protein
  for (i in coloc_aric_proteins){
    # Dataset 1: pqtl summary: N = 7213; dataset 2: cancer: s (case samplesize/total samplesize) and N (samplesize)
    tmp_coloc_aric_pqtl_cancer <- sig_tmp_aric_pqtl_exp_outc %>% filter(entrezgenesymbol == i)
    
    tmp_coloc_aric_pqtl_cancer <- tmp_coloc_aric_pqtl_cancer %>% na.omit()
    
    tmp_coloc_aric_pqtl_cancer <- tmp_coloc_aric_pqtl_cancer[tmp_coloc_aric_pqtl_cancer$P != 0,]
    
    tmp <- coloc.abf(dataset2 = list(pvalues = tmp_coloc_aric_pqtl_cancer$p_value, type = "cc", s = s_value, 
                                     N = N_value, snp = tmp_coloc_aric_pqtl_cancer$SNP),
                     dataset1 = list(pvalues = tmp_coloc_aric_pqtl_cancer$P, type = "quant", 
                                     N = 7213, snp = tmp_coloc_aric_pqtl_cancer$SNP), 
                     MAF = tmp_coloc_aric_pqtl_cancer$A1_FREQ)
    
    print(paste0('step 1: coloc-analysis=====', i))
    
    tmp <- tmp$results
    
    if (nrow(tmp) != 0) {
      tmp$protein <- i
      coloc_aric_pqtl_cancer <- rbind(coloc_aric_pqtl_cancer, tmp)
    }
  }
  
  file_name <- paste0("coloc_aric_pqtl_", tmp_cancer_name, ".csv")
  write.csv(coloc_aric_pqtl_cancer, file = file_name, row.names = FALSE)
  
  # Filter significant coloc results with PP.H4 > 0.8
  sig_coloc_aric_pqtl_cancer <- coloc_aric_pqtl_cancer %>% dplyr::filter(SNP.PP.H4 > 0.8)
  
  file_name <- paste0('sig_coloc_aric_pqtl_', tmp_cancer_name, '_pph40.8.csv')
  write.csv(sig_coloc_aric_pqtl_cancer, file = file_name, row.names = FALSE)
  
  
  # ========================================== PART 4 Perform: check association between pQTL & eQTL ======================================
  
  # 1) Select IVs(pQTL) info from aric_pqtl_cancer_singlesnp_results_clump data
  sig_aric_pqtl_cancer_ivs <- aric_pqtl_cancer_singlesnp_results_clump %>% 
    filter(aric_pqtl_cancer_singlesnp_results_clump$exposure %in% sig_aric_pqtl_cancer_singlesnp_results_clump$exposure) %>%
    filter(SNP != "All - Inverse variance weighted (multiplicative random effects)" & SNP != "All - MR Egger")
  
  # Horizontal pleiotropy test: check pQTL pleiotropy in (robust & probable MR results)
  # check if cis-pqtl associated with other proteins or genes
  sig_aric_pqtl_cancer_ivs_qc <- sig_aric_pqtl_cancer_ivs %>% 
    filter(!(exposure %in% duplicated_probes_protein)) 
  
  # 2) filter secondary genes (in the same pathway as original genes) (STRING database)
  secondary_pqtl <- sig_aric_pqtl_exp_clump %>%
    filter(SNP %in% sig_aric_pqtl_cancer_ivs_qc$SNP & !(exposure %in% sig_aric_pqtl_cancer_ivs_qc$exposure)) %>%
    dplyr::select(SNP, exposure) %>%
    left_join(sig_aric_pqtl_cancer_ivs_qc[, c("SNP", "proteins")], by = "SNP") %>%
    dplyr::rename(secondary_protein = exposure, original_protein = proteins)
  
  secondary_pqtl %>% distinct(secondary_protein, original_protein)
  
  unique(secondary_pqtl$secondary_protein)
  
  # load horizontal pleiotropy pQTL file (previously mapped secondary genes)
  secondary_pqtl_pleiotropy <- read.csv("C:/40cancers_analysis/Proteome-wide Association Analysis/secondary_pqtl_pleiotropy.csv")
  
  secondary_pqtl_pleiotropy <- secondary_pqtl_pleiotropy %>%
    filter(cancer_type == tmp_cancer_name) %>% 
    pull(gene_list)
  
  secondary_pqtl_cancer_singlesnp_results_clump_pqtl <- aric_pqtl_cancer_singlesnp_results_clump %>%
    filter(exposure %in% secondary_pqtl_pleiotropy)
  
  # Check if cis-pqtl associated with other genes
  # Read eQTL dataset 
  eqtl <- fread("C:/40cancers_analysis/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
                sep = "\t", header = TRUE)
  eqtl <- eqtl %>% filter(SNP %in% sig_aric_pqtl_cancer_ivs_qc$SNP & !(GeneSymbol %in% sig_aric_pqtl_cancer_ivs_qc$exposure))
  eqtl <- eqtl %>% filter(Pvalue < 5e-08 & FDR < 0.05)

  # Check how many proteins corresponding to these pQTL genes have protein abundance measured
  eqtl_seq <- eqtl %>% left_join(seq, by = c("GeneSymbol" = "entrezgenesymbol"))
  eqtl_seq <- eqtl_seq[!is.na(eqtl_seq$seqid_in_sample)]
  eqtl_seq <- eqtl_seq[order(eqtl_seq$SNP, eqtl_seq$Pvalue),]
  eqtl_seq <- eqtl_seq[!duplicated(eqtl_seq$SNP),]
  eqtl_seq <- eqtl_seq[order(eqtl_seq$seqid_in_sample, eqtl_seq$Pvalue),]
  
  # Identify the potential pleiotropy pQTL from ARIC
  tmp_eqtl_seq_proteins <- eqtl_seq$seqid_in_sample
  tmp_eqtl_seq_proteins <- paste(tmp_eqtl_seq_proteins, ".PHENO1.glm.linear", sep = "")
  
  # 3) identify the secondary proteins/genes
  tmp_eqtl_seq <- NULL
  
  for (i in tmp_eqtl_seq_proteins){
    tmp <- suppressWarnings(fread(paste0(aric_path, i),
                                  stringsAsFactors = FALSE,
                                  data.table = FALSE))
    print(paste0('step 1:', i))
    tmp <- tmp[tmp$ID %in% eqtl_seq$SNP,]
    
    if (nrow(tmp) != 0) {
      tmp$tissue <- gsub(aric_data, "", i)
      tmp_eqtl_seq <- rbind(tmp_eqtl_seq, tmp)
      print(paste0('step 2: combine the sig pQTL', gsub(aric_data, "", i)))
    }
  }
  tmp_eqtl_seq <- tmp_eqtl_seq %>% left_join(seq, by = c('tissue' = 'seqid_in_sample')) 
  
  # 4) perform two-sample MR analysis for secondary proteins and cancer
  # 4.1) Format input data
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
  file_name <- paste0("secondary_pqtl_eqtl_seq_exp_", tmp_cancer_name, ".csv")
  write.csv(secondary_pqtl_eqtl_seq_exp, file = file_name, row.names = FALSE)
  
  secondary_pqtl_eqtl_proteins <- unique(secondary_pqtl_eqtl_seq_exp$exposure)
  
  # 4.2) Clump input data 
  sig_secondary_pqtl_eqtl_exp_clump <- NULL
  
  for (protein in secondary_pqtl_eqtl_proteins) {
    tmp_secondary_pqtl_eqtl_seq_exp_clump <- secondary_pqtl_eqtl_seq_exp %>% 
      filter(exposure == protein)
    # Clump exp data
    tmp_secondary_pqtl_eqtl_seq_exp_clump <- tryCatch({
      ld_clump(dat = tibble(rsid = tmp_secondary_pqtl_eqtl_seq_exp_clump$SNP,
                            pval = tmp_secondary_pqtl_eqtl_seq_exp_clump$pval.exposure,
                            id = tmp_secondary_pqtl_eqtl_seq_exp_clump$id.exposure),
               clump_kb = 5000,
               clump_r2 = 0.01,
               clump_p = 1,
               bfile = "C:/40cancers_analysis/1kg.v3/EUR",
               plink_bin = "C:/40cancers_analysis/1kg.v3/plink.exe")
    }, error = function(e) {
      message("Error ignored: ", e$message)
      return(NULL)
    })
    sig_secondary_pqtl_eqtl_exp_clump <- rbind(sig_secondary_pqtl_eqtl_exp_clump, tmp_secondary_pqtl_eqtl_seq_exp_clump) 
  }
  sig_secondary_pqtl_eqtl_exp_clump <- secondary_pqtl_eqtl_seq_exp %>% filter(secondary_pqtl_eqtl_seq_exp$SNP %in% sig_secondary_pqtl_eqtl_exp_clump$rsid)
  
  file_name <- paste0('sig_secondary_pqtl_eqtl_exp_clump_', tmp_cancer_name, '.csv')
  write.csv(sig_secondary_pqtl_eqtl_exp_clump, file = file_name, row.names = FALSE)
  
  # 4.3) Format outcome data
  if ("effect_allele_frequency" %in% names(outc)) {
    # If effect_allele_frequency column exists, use it in format_data function
    secondary_outcome_dat <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_secondary_pqtl_eqtl_exp_clump$SNP,
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
    # If effect_allele_frequency column does not exist, do not use it in format_data function
    secondary_outcome_dat <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_secondary_pqtl_eqtl_exp_clump$SNP,
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
  
  # 4.4) Perform MR (with clump)
  secondary_H_data_clump <- NULL
  secondary_pqtl_cancer_singlesnp_results_clump <- NULL
  secondary_protein_name <- NULL
  
  secondary_proteins <- unique(sig_secondary_pqtl_eqtl_exp_clump$exposure)
  
  # Note: Not all proteins have results in outcome data, if not, cannot add PROTEINS
  for (protein in secondary_proteins) {
    tmp_sig_secondary_pqtl_eqtl_exp_clump <- sig_secondary_pqtl_eqtl_exp_clump %>% 
      filter(exposure == protein)
    
    secondary_protein_name <- c(secondary_protein_name, protein)
    mm <- length(secondary_protein_name)
    
    print(paste0("step1: filter proteins======", secondary_protein_name[mm]))
    
    tmp_secondary_H_data_clump <- harmonise_data(exposure_dat = tmp_sig_secondary_pqtl_eqtl_exp_clump,
                                                 outcome_dat = secondary_outcome_dat)
    secondary_H_data_clump <- rbind(secondary_H_data_clump, tmp_secondary_H_data_clump)
    print("step 2: Harmonize name")
    
    if (!(is.data.frame(tmp_secondary_H_data_clump) && nrow(tmp_secondary_H_data_clump) == 0)){
      tmp_secondary_pqtl_cancer_singlesnp_results_clump <- tmp_secondary_H_data_clump %>% 
        mr_singlesnp(
          parameters = default_parameters(),
          single_method = "mr_wald_ratio",
          all_method = c('mr_ivw_mre', 'mr_egger_regression', "mr_weighted_median", "mr_weighted_mode"))
      tmp_secondary_pqtl_cancer_singlesnp_results_clump$proteins <- protein
    }
    print("step 3: perform MR")
    
    secondary_pqtl_cancer_singlesnp_results_clump <- rbind(secondary_pqtl_cancer_singlesnp_results_clump,
                                                           tmp_secondary_pqtl_cancer_singlesnp_results_clump)
  }
  
  file_name <- paste0('secondary_pqtl_eqtl_', tmp_cancer_name, '_H_data_clump(all).csv')
  write.csv(secondary_H_data_clump, file = file_name, row.names = FALSE)
  
  # Merge secondary pQTL & eQTL MR results (horizontal pleiotropy proteins)
  secondary_pqtl_cancer_singlesnp_results_clump <- rbind(secondary_pqtl_cancer_singlesnp_results_clump, secondary_pqtl_cancer_singlesnp_results_clump_pqtl)
  
  file_name <- paste0('secondary_pqtl_eqtl_', tmp_cancer_name, '_singlesnp_results_clump(all).csv')
  write.csv(secondary_pqtl_cancer_singlesnp_results_clump, file = file_name, row.names = FALSE)
  
  # ============================================== PART 5 Delete the horizontal pQTL and run the MR again ================================================
 
   secondary_pqtl_eqtl_snp_list <- secondary_pqtl_cancer_singlesnp_results_clump %>% filter(!(SNP %in% c("All - Inverse variance weighted (multiplicative random effects)",
                                                                                                        "All - MR Egger", "All - Weighted median",
                                                                                                        "All - Weighted mode")))
  secondary_pqtl_eqtl_snp_list <- unique(secondary_pqtl_eqtl_snp_list$SNP)
  sig_aric_pqtl_exp_clump_nopleiotropy <- sig_aric_pqtl_exp_clump %>% filter(!(SNP %in% secondary_pqtl_eqtl_snp_list))
  
  file_name <- paste0("sig_aric_pqtl_exp_clump_nopleiotropy_", tmp_cancer_name, ".csv")
  write.csv(sig_aric_pqtl_exp_clump_nopleiotropy, file = file_name, row.names = FALSE)
  
  sig_aric_pqtl_exp_nopleiotropy <- aric_pqtl_exp %>% filter(!(SNP %in% secondary_pqtl_eqtl_snp_list))
  
  file_name <- paste0("sig_aric_pqtl_exp_nopleiotropy_", tmp_cancer_name, ".csv")
  write.csv(sig_aric_pqtl_exp_nopleiotropy, file = file_name, row.names = FALSE)
  
  # 1) Format outcome data
  if ("effect_allele_frequency" %in% names(outc)) {
    # If effect_allele_frequency column exists, use it in format_data function
    aric_outcome_dat_nopleiotropy <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_aric_pqtl_exp_clump_nopleiotropy$SNP,
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
    # If effect_allele_frequency column does not exist, do not use it in format_data function
    aric_outcome_dat_nopleiotropy <- format_data(
      dat = outc,
      type = "outcome",
      snps = sig_aric_pqtl_exp_clump_nopleiotropy$SNP,
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
  aric_outcome_dat_nopleiotropy$outcome <- tmp_cancer_name
  
  # 2) Perform MR (with clump)
  aric_H_data_clump_nopleiotropy <- NULL
  aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy <- NULL
  aric_protein_name_nopleiotropy <- NULL
  
  aric_proteins_nopleiotropy <- unique(sig_aric_pqtl_exp_clump_nopleiotropy$exposure)
  
  # Note: Not all proteins have results in outcome data, if not, cannot add PROTEINS
  for (protein in aric_proteins_nopleiotropy) {
    tmp_sig_aric_pQTL_exp_clump_nopleiotropy <- sig_aric_pqtl_exp_clump_nopleiotropy %>% 
      filter(exposure == protein)
    
    aric_protein_name_nopleiotropy <- c(aric_protein_name_nopleiotropy, protein)
    mm <- length(aric_protein_name_nopleiotropy)
    
    print(paste0("step1: filter proteins======", aric_protein_name_nopleiotropy[mm]))
    
    tmp_aric_H_data_clump_nopleiotropy <- harmonise_data(exposure_dat = tmp_sig_aric_pQTL_exp_clump_nopleiotropy,
                                                         outcome_dat = aric_outcome_dat_nopleiotropy)
    aric_H_data_clump_nopleiotropy <- rbind(aric_H_data_clump_nopleiotropy, tmp_aric_H_data_clump_nopleiotropy)
    print("step 2: Harmonize name")
    
    if (!(is.data.frame(tmp_aric_H_data_clump_nopleiotropy) && nrow(tmp_aric_H_data_clump_nopleiotropy) == 0)){
      tmp_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy <- tmp_aric_H_data_clump_nopleiotropy %>% 
        mr_singlesnp(
          parameters = default_parameters(),
          single_method = "mr_wald_ratio",
          all_method = c('mr_ivw_mre', 'mr_egger_regression', "mr_weighted_median", "mr_weighted_mode"))
      tmp_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$proteins <- protein
    }
    print("step 3: perform MR")
    
    aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy <- rbind(aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy,
                                                                   tmp_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy)
  }
  
  file_name <- paste0('aric_pqtl_', tmp_cancer_name, '_H_data_clump_nopleiotropy(all).csv')
  write.csv(aric_H_data_clump_nopleiotropy, file = file_name, row.names = FALSE)
  
  file_name <- paste0('aric_pqtl_', tmp_cancer_name, '_mr_results_clump_nopleiotropy(all).csv')
  write.csv(aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy, file = file_name, row.names = FALSE)
  
  # 3) ARIC(EA) quality control
  # Selecting significant results - proteins remaining after filtering based on previous step's robust/probable protein identification
  sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy <- aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy %>% filter(
    SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
               'All - MR Egger',
               "All - Weighted median",
               "All - Weighted mode") | p < 0.05)
  
  sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy <- sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy[!is.na(sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$b),]
  
  sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy <- sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy %>% 
    filter(exposure %in% sig_aric_pqtl_cancer_singlesnp_results_clump$exposure)
  
  sig_aric_proteins <- unique(sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$proteins)
  
  # 4) Check how many SNPs passing the PPH4 test are in the MR analysis IVS
  sig_aric_pqtl_cancer_ivs_coloc <- sig_coloc_aric_pqtl_cancer %>% 
    filter(protein %in% intersect(sig_coloc_aric_pqtl_cancer$protein, sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$exposure)) %>%
    filter(snp %in% aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$SNP)
  
  file_name <- paste0("sig_aric_pqtl_", tmp_cancer_name, "_ivs_coloc.csv")
  write.csv(sig_aric_pqtl_cancer_ivs_coloc, file = file_name, row.names = FALSE)
  
  sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy_coloc <- sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy %>% 
    filter(proteins %in% intersect(sig_coloc_aric_pqtl_cancer$protein, sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$exposure)) %>% 
    filter(SNP %in% aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy$SNP)
  
  file_name <- paste0("sig_aric_pqtl_", tmp_cancer_name, "_singlesnp_results_clump_nopleiotropy_coloc.csv")
  write.csv(sig_aric_pqtl_cancer_singlesnp_results_clump_nopleiotropy_coloc, file = file_name, row.names = FALSE)
  
})

# ============================================== PART 6 Perform: MR - analysis [ARIC(EA proteins) & risk factors] ==================
###
# risk factors (SMK,SBP,DBP,LDL,HDL,TC,Fasting Glu, Fasting Insu, BMI/Obesity) 
# select sig IVs for the previous proteins (sig MR & sig Coloc)

for(tmp_cancer_name in cancer_name){
  
  # load sig_aric_pqtl_exp_clump_nopleiotropy
  path_sig_aric_pqtl_exp_clump_nopleiotropy <- paste0("./riskfactor_pqtl_analysis/sig_aric_pqtl_exp_clump_nopleiotropy_", tmp_cancer_name, ".csv")
  sig_aric_pqtl_exp_clump_nopleiotropy <- read.csv(file = path_sig_aric_pqtl_exp_clump_nopleiotropy)
  
  #load sig_aric_pqtl_cancer_ivs_coloc
  path_sig_aric_pqtl_cancer_ivs_coloc <- paste0("./riskfactor_pqtl_analysis/sig_aric_pqtl_", tmp_cancer_name, "_ivs_coloc.csv")
  sig_aric_pqtl_cancer_ivs_coloc <- read.csv(file = path_sig_aric_pqtl_cancer_ivs_coloc)
  
  sig_aric_pqtl_riskfactor <- sig_aric_pqtl_exp_clump_nopleiotropy %>% 
    filter(exposure %in% sig_aric_pqtl_cancer_ivs_coloc$protein)
  
  print(sig_aric_pqtl_riskfactor)
  
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
  
  file_name <- paste0("aric_pqtl__",tmp_cancer_name,"_riskfactor_H_data_clump_nopleiotropy(all).csv")
  write.csv(riskfactor_H_data, file = file_name, row.names = F)
  
  file_name <- paste0("aric_pqtl_",tmp_cancer_name,"_riskfactors_nopleiotropy.csv")
  write.csv(riskfactor_mr_results, file = file_name, row.names = F)
}

# ========== riskfactor quality control =========

for(tmp_cancer_name in cancer_name){
  # load riskfactor_H_data
  path_riskfactor_H_data <- paste0("C:/Users/Kevin/Desktop/Table_csv/riskfactor_outcomes/aric_pqtl_",tmp_cancer_name,"_riskfactor_H_data_clump_nopleiotropy(all).csv")
  riskfactor_H_data <- read.csv(file = path_riskfactor_H_data)
  
  # load riskfactor_mr_results
  path_riskfactor_mr_results <- paste0("C:/Users/Kevin/Desktop/Table_csv/riskfactor_outcomes/aric_pqtl_",tmp_cancer_name,"_riskfactors_nopleiotropy.csv")
  riskfactor_mr_results <- read.csv(file = path_riskfactor_mr_results)
  
  ## selecting significant risk-factors
  sig_riskfactor_mr_results <- riskfactor_mr_results %>% filter(
    SNP %in% c('All - Inverse variance weighted (multiplicative random effects)',
               'All - MR Egger') | p < 0.05)
  
  sig_riskfactor_mr_results <- sig_riskfactor_mr_results[!is.na(sig_riskfactor_mr_results$b),]
  
  # ARIC(EA) add nsnp to the significant results
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
  
  # reassign method (sig_riskfactor_mr_results contain MR results of all proteins)
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
  
  
  # select significant MR results of 9 proteins and 12 risk factors
  sig_riskfactor_mr_results <- sig_riskfactor_mr_results %>% filter(p < 0.05/12)
  
  # save file (MR-risk factors) proteins pass P< 0.05/12 threshold
  file_name <- paste0("sig_aric_pqtl_",tmp_cancer_name,"_riskfactors.csv")
  write.csv(sig_riskfactor_mr_results, file = file_name, row.names = F)
}

# ================================ PART 7 Coloc-localization Analysis (Supplementary calculation) ================================

# ======= step 1: Coloc-localization Analysis of ARIC(EA proteins) - cancers =====================

sig_coloc_aric_pqtl_cancer <- read.csv(file = "C:/Users/Kevin/Desktop/Table_csv/pp.h4.csv")

df <- sig_coloc_aric_pqtl_cancer %>% filter(protein == "CFD", outcome == "Oral_cavity_cancer")

# build dataset1(pQTL) and dataset2(cancers)
dataset1 <- list(
  pvalues = df$pvalues.df1,
  N = df$N.df1,
  MAF = df$MAF.df1,
  beta = df$z.df1,
  varbeta = df$V.df1,
  snp = df$snp,
  type = "quant" # "quant" for continuous traits, "cc" for case-control studies
)

dataset2 <- list(
  pvalues = df$pvalues.df2,
  N = df$N.df2,
  MAF = df$MAF.df2,
  beta = df$z.df2,
  varbeta = df$V.df2,
  snp = df$snp,
  type = "cc" # Modify to an appropriate type
)

# calculate sdY = sqrt(V * N)
dataset1$sdY <- sqrt(df$V.df1 * df$N.df1)
dataset2$sdY <- sqrt(df$V.df2 * df$N.df2)

# perform colocalisation analysis
coloc_result <- coloc.abf(dataset1, dataset2, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)

coloc_result$summary


# ========= step 2: Coloc-localization ARIC(EA proteins) & riskfactor ===========

# perform coloc analysis

for(tmp_cancer_name in cancer_name){
  path_riskfactor_H_data <- paste0("C:/Users/Kevin/Desktop/Table_csv/Echo2me/riskfactor_outcomes/aric_pqtl_",tmp_cancer_name,"_riskfactor_H_data_clump_nopleiotropy(all).csv")
  riskfactor_H_data <- read.csv(path_riskfactor_H_data)
  
  path_sig_riskfactor_mr_results <- paste0("C:/Users/Kevin/Desktop/Table_csv/table5_riskfactor/table5_pqtl_riskfactor/sig_aric_pqtl_",tmp_cancer_name,"_riskfactors.csv")
  sig_riskfactor_mr_results <- read.csv(path_sig_riskfactor_mr_results)
  
  # filter sig_riskfactor_H_data for coloc analysis
  sig_riskfactor_H_data <- riskfactor_H_data %>% filter(exposure %in% sig_riskfactor_mr_results$exposure)
  
  # add additional info for proteins_riskfactors_H_data to perfrom the coloc analysis
  coloc_aric_riskfactor_proteins <- unique(sig_riskfactor_H_data$exposure)
  riskfactors <- unique(sig_riskfactor_H_data$id.outcome)
  
  sig_riskfactor_H_data$samplesize.outcome[sig_riskfactor_H_data$id.outcome == "ieu-b-109"] <- 403943
  sig_riskfactor_H_data$samplesize.outcome[sig_riskfactor_H_data$id.outcome == "ieu-b-5089"] <- 201678
  
  coloc_aric_pqtl_riskfactor <- NULL
  
  for (i in coloc_aric_riskfactor_proteins){
    # dataset 1: pqtl summary; dataset 2: riskfactors
    print(paste0('step 1: coloc-analysis=====',i))
    
    for (riskfactor in riskfactors) {
      tmp_coloc_aric_pqtl_riskfactor <- sig_riskfactor_H_data %>% filter(exposure == i & id.outcome == riskfactor)
      print(paste0('exposure===',i,'; outcome====',riskfactor))
      
      
      tmp_coloc_aric_pqtl_riskfactor <- tmp_coloc_aric_pqtl_riskfactor[order(tmp_coloc_aric_pqtl_riskfactor[,'SNP'], tmp_coloc_aric_pqtl_riskfactor[,'pval.exposure']),]
      tmp_coloc_aric_pqtl_riskfactor <- tmp_coloc_aric_pqtl_riskfactor[!duplicated(tmp_coloc_aric_pqtl_riskfactor$SNP), ]
      print('step 2: remove duplicates')
      
      tmp_coloc_aric_pqtl_riskfactor <- tmp_coloc_aric_pqtl_riskfactor %>% filter(!is.na(eaf.outcome))
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
  
  file_name <- paste0("coloc_aric_",tmp_cancer_name,"_coloc_proteins_riskfactor.csv")
  write.csv(coloc_aric_pqtl_riskfactor, file = file_name, row.names = F)
  
  if(!is.null(coloc_aric_pqtl_riskfactor)){
    sig_coloc_aric_pqtl_riskfactor <- coloc_aric_pqtl_riskfactor %>% dplyr::filter(SNP.PP.H4>=0.8)
    sig_coloc_aric_pqtl_riskfactor <- sig_coloc_aric_pqtl_riskfactor %>% dplyr::rename(exposure=protein,id.outcome=outcome)
    
    file_name <- paste0("sig_coloc_aric_",tmp_cancer_name,"_coloc_proteins_riskfactor.csv")
    write.csv(sig_coloc_aric_pqtl_riskfactor, file = file_name, row.names = F)
  }
  
}

# ======= step 3: coloc analysis for (H0, H1, H2, H3, H4) ARIC(EA proteins) - cancer ==========

sig_coloc_aric_pqtl_cancer <- read.csv(file = "C:/Users/Kevin/Desktop/Table_csv/riskfactor_pp.h4.csv")

df <- sig_coloc_aric_pqtl_cancer %>% filter(exposure == "CFD", id.outcome == "ebi-a-GCST90002238")

# build dataset1(pQTL) and dataset2(riskfactor)  
dataset1 <- list(
  pvalues = df$pvalues.df1,
  N = df$N.df1,
  MAF = df$MAF.df1,
  beta = df$z.df1,
  varbeta = df$V.df1,
  snp = df$snp,
  type = "quant" # "quant" for continuous traits, "cc" for case-control studies
)

dataset2 <- list(
  pvalues = df$pvalues.df2,
  N = df$N.df2,
  MAF = df$MAF.df2,
  beta = df$z.df2,
  varbeta = df$V.df2,
  snp = df$snp,
  type = "quant" # Modify the appropriate type
)

# calculate sdY = sqrt(V * N)
dataset1$sdY <- sqrt(df$V.df1 * df$N.df1)
dataset2$sdY <- sqrt(df$V.df2 * df$N.df2)

# perform coloc-analysis
coloc_result <- coloc.abf(dataset1, dataset2, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)

coloc_result$summary

df

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

riskfactor_pp.h4 <- riskfactor_pp.h4 %>% left_join(riskfactor_lab, by="id.outcome")
write.csv(riskfactor_pp.h4, file = "riskfactor_pp.h4.csv", row.names = F)