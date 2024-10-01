rm(list = ls())
#source("~/Documents/OneDrive - Texas A&M University/MLFDR/semifinal_fun.R")
source("~/Downloads/dact.R")
source("~/Documents/OneDrive - Texas A&M University/MLFDR/EM_general_fun.R")

library(tidyverse)
library(dplyr)
library(janitor)
library(UCSCXenaTools)
library(HDMT)
library(qvalue)
library(forecast)
library(compositions)
#library(DACT)
library(stringi)

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
samples = samples(meth450)

XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "LUNG") %>%
  XenaFilter(filterDatasets = "TCGA.LUNG.sampleMap/HiSeqV2") -> genex

XenaQuery(genex) %>%
  XenaDownload() ->xe_genex

##exposures and confounders
XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "LUNG") %>%
  XenaFilter(filterDatasets = "clinical") -> cli
XenaQuery(cli) %>%
  XenaDownload() -> xe_cli 

clinical = XenaPrepare(xe_cli)


df_genex = XenaPrepare(xe_genex)
df_genexv2 = df_genex[[2]]

probemap = read_tsv("genecode.tsv")
probemap2 = probemap[,1:2]
names(probemap2) = c("id", "gene")
x = unlist(df_meth[,1])
y =unlist(probemap2[,1])
common_id = intersect(x,y)
u = as.factor(common_id)


names(df_meth)[1] = "id"



probemap2 %>%
  filter(gene != ".") -> probemap3
sthing = sub(",.*", " ", probemap3$gene)
probemap3$gene = sthing

names(df_genexv2)[1] = "id"

df1 = df_meth %>% inner_join(probemap3, by = "id")
df1$gene = trimws(df1$gene)
df1$gene = as.factor(df1$gene)
length(levels(df1$gene))

df1[,-1] %>%
  group_by(gene) %>%
  summarise_all(mean, na.rm = TRUE) -> testdf

#summarise(across(everything(), mean(.,na.rm = TRUE))) ->testdf
names(df_genexv2)[1] = "gene"
common_samples = intersect(names(testdf), names(df_genexv2))
common_samples = intersect(common_samples, c("gene", clinical$sampleID))

X = clinical$tobacco_smoking_history
Z1 = clinical$age_at_initial_pathologic_diagnosis
Z2 = clinical$gender
expo = data.frame(X, Z1, Z2, row.names = clinical$sampleID)

expo_new = na.omit(expo[common_samples,])
common_samples = c("gene", rownames(expo_new))

df_genexv2[,common_samples] -> df2  ##geneexpression data
testdf[,common_samples] -> df3


#common_genes = intersect(testdf$gene, df_genexv2$gene)

df4 = inner_join(df2, df3, by = "gene") ##contains genexpression|DNAmethylation data

common_genes = df4$gene
#write.csv(df4, file = "alldata.csv")
write.csv(expo_new, file = "exposure.csv")
##Analysis of the data

m = length(common_genes) #16544
n = nrow(expo_new)  #798

##extracting Y-linked genes
missing = vector()



##Computing the coefficients

alpha = vector()
beta = vector()
p1 = vector()
p2 = vector()
sd1 = vector()
sd2 = vector()

M = df4[,(n+2):(2*n + 1)]
M_new = replace(M, is.na(M), 0)
Y = df4[,2:(n+1)]
M = M_new
#M = na.omit(M)
#X = factor(expo_new$X)
##permute the variables, break the dependence
X = expo_new$X[sample(1:n)]
M = M[sample(1:m),]
Y = Y[sample(1:m),]
for(i in 1:m)
{
  
  obj1 = lm(unlist(M[i,])~ expo_new$X + expo_new$Z1 + expo_new$Z2)
  p1[i] = summary(obj1)$coefficients[2,4]
  obj2 = lm(unlist(Y[i,]) ~ unlist(M[i,]) + expo_new$X + expo_new$Z1 + expo_new$Z2)
  p2[i] = summary(obj2)$coefficients[2,4]
  sd1[i] = summary(obj1)$coefficients[2,2]
  sd2[i] = summary(obj2)$coefficients[2,2]
  
  alpha[i] = obj1$coefficients[2]
  beta[i] = obj2$coefficients[2]
  rm(obj1)
  rm(obj2)
  print(i)
}

#holandric = common_genes[missing >= 324]

coeff_pvals = na.omit(data.frame(common_genes,alpha,beta, p1, p2, sd1, sd2))

#write.csv(coeff_pvals, file = "coeff_pvals.csv")

counter = !is.na(p2)
Y_new = Y[counter,]
M_new = M[counter,]
common_genes = common_genes[counter]
alpha = alpha[counter]
beta = beta[counter]
p1 = p1[counter]
p2 = p2[counter]




pi10 = qvalue(p1)$pi0
pi01 = qvalue(p2)$pi0
input_pvalues = na.omit(cbind(p1, p2))

##DACT
pmax = apply(input_pvalues, 1, max)
p_dact = DACT(p1, p2)
##null estimation
rej_dact = qvalue(p_dact)$qvalues<= size
sum(rej_dact)

##HDMT
nullprop = null_estimation(input_pvalues)
fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                         nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
threshhold = max(pmax[fdr_hdmt<= size])
rej_HDMT = pmax <= threshhold
sum(rej_HDMT)


##MLFDR
Y = matrix(unlist(Y_new), ncol = 798)
M = matrix(unlist(M_new), ncol = 798)
index_mat = rowSums(Y == 0) < 80
M = M[index_mat,]
Y = Y[index_mat,]
Y = as.matrix(clr(Y + 0.5))
M = as.matrix(clr(M + 0.5))
m = nrow(Y)
pi_start = c(nullprop$alpha00, nullprop$alpha10, nullprop$alpha01, 1-(nullprop$alpha00+ nullprop$alpha10 + nullprop$alpha01))


alpha = c()
beta = c()
sd1 = c()
sd2 = c()
lambda1 = vector()
lambda2 = vector()
for(i in 1:m)
{
  #lambda1[i] = BoxCox.lambda(M[i,]+5, lower = -2, upper = 5)
  #lambda2[i] = BoxCox.lambda(Y[i,] + 5, lower = -2, upper = 5)
  
  obj1 = lm(M[i,]~ expo_new$X + expo_new$Z1)
  p1[i] = summary(obj1)$coefficients[2,4]
  obj2 = lm(Y[i,] ~ M[i,] + expo_new$X + expo_new$Z1)
  p2[i] = summary(obj2)$coefficients[2,4]
  sd1[i] = summary(obj1)$coefficients[2,2]
  sd2[i] = summary(obj2)$coefficients[2,2]
  
  
  alpha[i] = obj1$coefficients[2]
  beta[i] = obj2$coefficients[2]
  rm(obj1)
  rm(obj2)
  print(i)
}

#write.csv(data.frame(alpha, beta), "coeff_mlfdr.csv")
#write.csv(data.frame(expo_new$X, expo_new$Z1), "expo_mlfdr.csv")
#write.csv(Y, "Y_mlfdr.csv")
#write.csv(M, "M_mlfdr.csv")
#obj = maximization(alpha, beta, expo_new$X, Y, M, pi_start, maxiter = 1000)


vector1 = alpha/sd1 
vector2 = beta/sd2

