###GSMR
##按照gsmr默认参数调用gsmr
gsmr_own<-function(dat_gsmr,mat,gwas_thres)
{
dat_gsmr = dat_gsmr[dat_gsmr$mr_keep==T,]
dat_gsmr = dat_gsmr[order(dat_gsmr$pval.exposure),]
dat_gsmr = dat_gsmr[!duplicated(dat_gsmr$SNP),]

gsmr_data = dat_gsmr%>%rename('a1'='effect_allele.exposure', 'a2'='other_allele.exposure',
                         'bzx_pval'='pval.exposure', 'bzy_pval'='pval.outcome',
                         'bzx'='beta.exposure', 'bzy'='beta.outcome',
                         'bzx_se'='se.exposure', 'bzy_se'='se.outcome',
                         'bzx_n'='samplesize.exposure', 'bzy_n'='samplesize.outcome',
                         'a1_freq'='eaf.exposure')

snpfreq = gsmr_data$a1_freq             # allele frequencies of the SNPs
bzx = gsmr_data$bzx     # effects of the instruments on risk factor
bzx_se = gsmr_data$bzx_se       # standard errors of bzx
bzx_n = gsmr_data$bzx_n          # GWAS sample size for the risk factor
std_zx = std_effect(snpfreq, bzx, as.numeric(bzx_se), as.numeric(bzx_n))    # perform standardisation
gsmr_data$std_bzx = std_zx$b    # standardized bzx
gsmr_data$std_bzx_se = std_zx$se    # standardized bzx_se

ldrho = mat[rownames(mat) %in% gsmr_data$SNP, colnames(mat) %in% gsmr_data$SNP]
snp_coeff_id = rownames(ldrho)
gsmr_data = gsmr_data[match(unique(gsmr_data$SNP), gsmr_data$SNP),] # sometime has duplicate SNP

n_ref <- 489    # Sample size of the 1000g eur (nrow of fam)
gwas_thresh <- gwas_thres    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
single_snp_heidi_thresh <- 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis | default is 0.01
multi_snp_heidi_thresh <- 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis | default is 0.01
nsnps_thresh <- 10   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag <- F    # flag for HEIDI-outlier analysis
ld_r2_thresh <- 0.05    # LD r2 threshold to remove SNPs in high LD 
ld_fdr_thresh <- 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr2_beta <- 1

gsmr_results = gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval,
                    gsmr_data$bzy, gsmr_data$bzy_se, gsmr_data$bzy_pval,
                    ldrho, snp_coeff_id, n_ref,
                    heidi_outlier_flag, gwas_thresh,
                    single_snp_heidi_thresh, multi_snp_heidi_thresh,
                    nsnps_thresh, ld_r2_thresh, ld_fdr_thresh,
                    gsmr2_beta)
filtered_index=gsmr_results$used_index
bzy = gsmr_data$bzy    # SNP effects on the disease
bzy_se = gsmr_data$bzy_se    # standard errors of bzy
bzy_pval = gsmr_data$bzy_pval
re<-c(gsmr_results$bxy,gsmr_results$bxy_se,gsmr_results$bxy_pval,length(gsmr_results$used_index))
return(re)
}
