library(MRPRESSO)
library(TwoSampleMR)
library(dplyr)
library(gsmr)
library(MRlap)
library(data.table)

source("/home/yanglab_data/user/qiml/ukb/mr/code/GSMR_own.R")
source("/home/yanglab_data/user/qiml/ukb/mr/code/MR_own_all.R")

mr_list<-read.table("../list_mr",header=F)
file_exp_mat<-"/home/yanglab_data/user/qiml/ukb/heart_dis_gwas/finngen/5e-8"
file_exp<-"/home/yanglab_data/user/qiml/ukb/heart_dis_gwas/finngen/clean"
file_out<-"/home/yanglab_data/user/qiml/ukb/metabolites/ieu/clean_data"
exp_index<-2
out_index<-1

p_threshold<-5e-8

MR_re<-NULL
for (i in 1:nrow(mr_list))
{
if (file.exists(sprintf("./%s_%s",mr_list[i,exp_index],mr_list[i,out_index])))
{
dat_re<-read.table(sprintf("./%s_%s",mr_list[i,exp_index],mr_list[i,out_index]),header=T,sep="\t")
MR_re<-rbind(MR_re,dat_re)
}
if (!file.exists(sprintf("./%s_%s",mr_list[i,exp_index],mr_list[i,out_index])))
{
try({
exp_mat<-read.table(sprintf("%s/finngen_R6_%s.R2Matrix.ld",file_exp_mat,mr_list[i,exp_index]))
snplist<-read.table(sprintf("%s/finngen_R6_%s.R2Matrix.snplist",file_exp_mat,mr_list[i,exp_index]))[,1]
colnames(exp_mat)<-snplist
rownames(exp_mat)<-snplist
exp_ori<-fread(sprintf("%s/finngen_R6_%s_clean",file_exp,mr_list[i,exp_index]))
out_ori<-fread(sprintf("%s/%s_clean",file_out,mr_list[i,out_index]))

exp_file<-exp_ori[exp_ori$SNP %in% snplist,]
out_file<-out_ori[out_ori$SNP %in% snplist,]
exp_ori = exp_ori%>%rename(rsid=SNP,a1=A1,a2=A2,beta=BETA,se=SE,p=P)
out_ori = out_ori%>%rename(rsid=SNP,a1=A1,a2=A2,beta=BETA,se=SE,p=P)
results<-MR_own(exp_file,out_file,p_threshold,exp_mat,exp_ori,out_ori)
results$exposure<-mr_list[i,exp_index]
results$outcome<-mr_list[i,out_index]
write.table(results,sprintf("./%s_%s",mr_list[i,exp_index],mr_list[i,out_index]),row.names=F,sep="\t",quote=F)
MR_re<-rbind(MR_re,results)
})
}
}
write.table(MR_re,"../MR_results_all",row.names=F,sep="\t",quote=F)
