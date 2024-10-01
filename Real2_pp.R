##tcga_lung_cox

rm(list = ls())
#source("~/Documents/OneDrive - Texas A&M University/MLFDR/semifinal_fun.R")
source("~/Downloads/dact.R")
source("~/Documents/OneDrive - Texas A&M University/MLFDR/EM_general_fun.R")
library(readr)
library(tidyverse)
library(dplyr)
library(janitor)
library(UCSCXenaTools)
library(HDMT)
library(qvalue)
library(forecast)
library(compositions)
library(dplyr)
#library(DACT)


showTCGA(project = "LUNG")
XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "LUNG") -> xenalung

XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "LUNG") %>%
  XenaFilter(filterDatasets = "TCGA.LUNG.sampleMap/HumanMethylation450") -> meth450

##Use datasets(xenalung)/cohort(xenalung), whichever applicable to find datasets within a Xenahub object.

XenaQuery(meth450) %>%
  XenaDownload() -> xe_download 

df_meth = XenaPrepare(xe_download)
#samples = samples(meth450)


##exposures and confounders
XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "LUNG") %>%
  XenaFilter(filterDatasets = "clinical") -> cli
XenaQuery(cli) %>%
  XenaDownload() -> xe_cli 

clinical = XenaPrepare(xe_cli)


#XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
 # XenaFilter(filterDatasets = "LUNG") %>%
#  XenaFilter(filterDatasets = "survival") -> surv

#XenaQuery(surv) %>%
#  XenaDownload() -> xe_surv 


#survival = XenaPrepare(xe_surv)


survival = data.frame(clinical$sampleID, clinical$days_to_initial_pathologic_diagnosis, clinical$days_to_death, clinical$days_to_last_followup, clinical$number_pack_years_smoked, clinical$age_at_initial_pathologic_diagnosis)
#attach(survival)
#detach(survival)
##filter out NA rows
survival [!(is.na(survival[,3]) * is.na(survival[,4])),] -> surv2
#filterout NA confounders
surv2[!is.na(surv2[,6]),] -> surv2

##Calculate survival time and censoring information
surv_time = c()
censoring = c()
dtd = surv2[,3]
dtc = surv2[,4]
index = 0
for(i in 1:nrow(surv2))
{
  if(is.na(dtd[i]) && !is.na(dtc[i])) {
    surv_time[i] = dtc[i]
    censoring[i] = 1
  }else if (!is.na(dtd[i]) && is.na(dtc[i])){
    surv_time[i] = dtd[i]
    censoring[i] = 0
  }else if(!is.na(dtd[i]) && !is.na(dtc[i])){
    surv_time[i] = min(dtc[i], dtd[i])
    censoring[i] = ifelse(dtc[i] >= dtd[i], 0, 1)
  }
}
surv2 = cbind(surv2, surv_time, censoring)

##remove rows with zero survival times
flags = which(surv_time <= 0)
surv2 %>% filter(!row_number() %in% flags ) -> surv2

##remove duplicate samples
surv2[,1] = substr(surv2[,1], 1, 12)
surv2 <- surv2[!duplicated(surv2),]

flags2 = which(is.na(surv2[,5]))
##Survival data is ready at this step.
surv2 %>% filter(!row_number() %in% flags2 ) -> surv2

colnames(surv2) = c("sample", "diag.init", "daystodeath", "lastfollowup", "smoking","age_at_diag", "surv_time", "censoring")
##Survival+exposure data is now ready


#Processing the cpg data
cpg_sites = df_meth$sample
sample = substr(names(df_meth)[-1], 1, 12)
names(df_meth) = c("sample", sample)
common_samples = intersect(sample, surv2[,1])
df_meth2 <- df_meth[,common_samples]
df_meth2 <- cbind(cpg_sites, df_meth2)
df_meth3 <- na.omit(df_meth2)
cpg_sites <- df_meth3$cpg_sites
##df_meth is ready


df_meth4 <- as.data.frame(t(df_meth3[,-1]))
names(df_meth4) <- cpg_sites
df_meth4$sample = rownames(df_meth4)

data = inner_join(surv2, df_meth4, by = "sample")
#write_csv(data, "tcga_cox_data.csv")
#rm(list = ls())
#.rs.restartR()
#data = read_csv("tcga_cox_data.csv")
data_sub = data[,1:1000]
write_csv(data_sub, "tcga_cox_data_subset2.csv")
##data is ready
data_sub = read_csv("tcga_cox_data_subset.csv")
n = nrow(data_sub)

meth = data_sub[,8:690]
surv = data_sub$surv_time
cens = data_sub$censoring
smoking = as.numeric(data_sub$smoking)
rm(data_sub)
gc()
#meth = scale(meth)
smoking = scale(smoking)
p = ncol(meth)
alpha <- beta <- var_alpha <- var_beta <- c()
p1 <- p2 <- c()
for(i in 1:p)
{
  obj1 = lm(scale(meth[,i])~ -1 + smoking)
  alpha[i] = obj1$coefficients
  var_alpha[i] = coef(summary(obj1))[1,2]^2
  p1[i] = coef(summary(obj1))[1,4]
  #print(i)
  #on.exit(gc())
}
obj2 = cox_inference(meth, surv, cens, kk = 5)
beta_hat = obj$bhat
var_beta = obj$var
sum((beta - beta_hat)^2)