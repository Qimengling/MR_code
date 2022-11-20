MR_own<-function(exp_file,out_file,MR_threshold,mat,exp_ori,out_ori)
{
    ###融合数据
df1 = exp_file
df2 = out_file
df1 = df1 %>% rename('pval.exposure' = 'P', 'effect_allele.exposure' = 'A1', 'other_allele.exposure' = 'A2',
                     'samplesize.exposure' = 'N',
                     'beta.exposure' = 'BETA', 'se.exposure' = 'SE', 'eaf.exposure' = 'FRQ') %>% mutate(exposure=rep("exposure",dim(df1)[1]),id.exposure=rep("drug",dim(df1)[1]))
df2 = df2 %>% rename('pval.outcome' = 'P', 'effect_allele.outcome' = 'A1', 'other_allele.outcome' = 'A2',
                     'samplesize.outcome' = 'N',
                     'beta.outcome' = 'BETA', 'se.outcome' = 'SE', 'eaf.outcome' = 'FRQ') %>% mutate(outcome=rep("outcome",dim(df2)[1]),id.outcome=rep("ecg",dim(df2)[1]))
dat = harmonise_data(df1, df2)
###直接去掉outcome中显著的SNP，不做异质性检验
dat<-dat[dat$mr_keep==TRUE,]
dat<-dat[which(dat$pval.outcome>1e-5),]

temp<-mr_pleiotropy_test(dat)
try({
if (temp$pval<0.05)
{
    temp<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  
          SignifThreshold = 0.05)
    dat<-dat[-temp$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,]
}
})
mr_re<-mr(dat)

MR_egger_re<-mr_egger_regression(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
single<-mr_leaveoneout(dat)
res_single <- mr_singlesnp(dat)

try({
gsmr_re<-gsmr_own(dat,mat,MR_threshold)
})
mr_re<-mr_re[-4,-c(1:2)]

re<-mr_re
try({
re<-rbind(re,data.frame(exposure="exposure",outcome="outcome",method="GSMR",nsnp=gsmr_re[4],b=gsmr_re[1],se=gsmr_re[2],pval=gsmr_re[3]))
})
###MRlap
source("/home/qiml/soft/MRlap-master/R/get_correction.R")
source("/home/qiml/soft/MRlap-master/R/run_LDSC.R")
source("/home/qiml/soft/MRlap-master/R/tidy_inputGWAS.R")

verbose<-TRUE
ld<-"/home/qiml/work/ukb_ecg/ref_ld/eur_w_ld_chr/"
hm3<-"/home/qiml/work/ukb_ecg/ldsc/w_hm3.snplist"
save_logfiles<-FALSE

exposure_data = tidy_inputGWAS(exp_ori, verbose)
outcome_data = tidy_inputGWAS(out_ori, verbose)
LDSC_results = run_LDSC(exposure_data, "exposure",outcome_data, "outcome", ld, hm3, save_logfiles, verbose)

###dat是MR-PRESSO的filter后的结果
###标准化beta,se,跑MR
dat$Z<-dat$beta.exposure/dat$se.exposure
dat$Z.exposure<-dat$beta.exposure/dat$se.exposure
dat$Z.outcome<-dat$beta.outcome/dat$se.outcome
dat$std_beta.exp<-dat$Z.exposure/sqrt(dat$samplesize.exposure)
dat$std_beta.out<-dat$Z.outcome/sqrt(dat$samplesize.outcome)
dat$std_SE.exp<-1/sqrt(dat$samplesize.e)
dat$std_SE.exp<-1/sqrt(dat$samplesize.exposure)
dat$std_SE.out<-1/sqrt(dat$samplesize.outcome)
dat$p.exp<-2*stats::pnorm(-abs(dat$Z.exposure))
dat$p.out<-2*stats::pnorm(-abs(dat$Z.outcome))
dat1<-dat[which(dat$mr_keep=="TRUE"),]
res_MR_TwoSampleMR<-TwoSampleMR::mr_ivw(dat1$std_beta.exp, dat1$std_beta.out,
                      dat1$std_SE.exp, dat1$std_SE.out)

MR_results<-list("alpha_obs" = res_MR_TwoSampleMR$b,"alpha_obs_se" = res_MR_TwoSampleMR$se,"n_exp" = round(mean(dat1$samplesize.exposure)),"n_out" = round(mean(dat1$samplesize.outcome)),"IVs" = dat1[,c(32,34)],"IVs_rs" = dat1$SNP)
try({
correction_results = with(c(MR_results, LDSC_results),
    get_correction(IVs, lambda, lambda_se, h2_LDSC, h2_LDSC_se,
                                      alpha_obs, alpha_obs_se,
                                      n_exp, n_out, MR_threshold, verbose))
re<-rbind(re,data.frame(exposure="exposure",outcome="outcome",method="MRlap",nsnp=nrow(MR_results$IVs),b=correction_results$alpha_corrected,se=correction_results$alpha_corrected_se,pval=2*stats::pnorm(-abs(correction_results$alpha_corrected/correction_results$alpha_corrected_se))))
})
return(re)
}
