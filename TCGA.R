##TCGA data analysis
library(HDMT)
library(DACT)
library(qvalue)
library(ggvenn)
library(xtable)
#Run from line 22
load("TCGA_mediation_result.RData")

alpha = TCGA_mediation_results$meqtl_beta
beta = TCGA_mediation_results$methy_beta

t_alpha = TCGA_mediation_results$meqtl_tstat
t_beta = TCGA_mediation_results$methy_tstat

var_alpha = (alpha/t_alpha)^2
var_beta = (beta/t_beta)^2
p_alpha = TCGA_mediation_results$meqtl_pvalue
p_beta = TCGA_mediation_results$methy_pvalue


obj1 = EM_fun(cbind(alpha, beta), k = 4, var_alpha, var_beta, lambda.init = c(0.4, 0.1, 0.4, 0.1), kappa.init = 0.001, psi.init = 10)
obj1
beepr::beep(4)

saveRDS(obj1, file = "TCGA_EM_fit.RDS")

###Run from here in future
fit = readRDS("TCGA_EM_fit.RDS")
pi = fit$lambda
mu = fit$mu
k = length(mu)
sigma = fit$sigma
lfdr = vector()
m = length(alpha)
k = 4
t = matrix(nrow = m, ncol = k)

load("TCGA_mediation_result.RData")

alpha = TCGA_mediation_results$meqtl_beta
beta = TCGA_mediation_results$methy_beta

x = cbind(alpha, beta)
for(i in 1:m)
{
  for(j in 1:k)
  {
    t[i,j] = pi[j]*emdbook::dmvnorm(x[i,], mu[[j]], sigma[j,i,,])
  }
  
  lfdr[i] = (t[i,1] + t[i,2] + t[i,3])/(t[i,1] + t[i,2] + t[i,3] + t[i,4])
  
}
st.lfdr<-sort(lfdr)
k=1
size = 0.01
while(k<m && ((1/k)*sum(st.lfdr[1:k])) <= size){
  k=k+1
}

k<-k-1
lfdrk<-st.lfdr[k]
snp_data = cbind(TCGA_mediation_results[,c(1:4,19)], lfdr)
reject<- snp_data[lfdr<=lfdrk,]
library(dplyr)
rej_MLFDR = reject %>% arrange(pmax)
rej_MLFDR

##HDMT
input_pvalues = cbind(p_alpha, p_beta)
pmax = apply(input_pvalues, 1, max)
nullprop = null_estimation(input_pvalues)
fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                         nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)
threshhold = max(pmax[fdr_hdmt<= size])
rej_HDMT = snp_data[pmax <= threshhold,] %>% arrange(pmax)


##DACT
p_dact = DACT(p_alpha, p_beta, correction = "JC")
p_dact_corr = qvalue(p_dact, pi0 = 1)$qvalues
rej_DACT = snp_data[ p_dact_corr <= size,]


saveRDS(rej_MLFDR, file = "rej_MLFDR.RDS")
saveRDS(rej_HDMT, file = "rej_HDMT.RDS")
saveRDS(rej_DACT, file = "rej_HDMT.RDS")

temp = layer_data(pp,1)
temp$colour

# Assuming df1 and df2 are your DataFrames

# Rows in df1 but not in df2
unique_to_df1 <- anti_join(rej_MLFDR, rej_HDMT, by = c("SNP", "cpg", "cpg_annogene", "genexp", "pmax", "lfdr"))

unique_to_MLFDR <- anti_join(unique_to_df1, rej_DACT, by = c("SNP", "cpg", "cpg_annogene", "genexp", "pmax", "lfdr"))
print(xtable(unique_to_MLFDR, digits = c(0,0,0,0,0,4,4)), type = "latex")
# Rows in df2 but not in df1
unique_to_df2 <- anti_join(rej_HDMT, rej_MLFDR, by = c("SNP", "cpg", "cpg_annogene", "genexp", "pmax", "lfdr"))

# View results
print("Rows in MLFDR but not in HDMT:")
print(unique_to_df1)

print("Rows in HDMT but not in MLFDR:")
print(unique_to_df2)

#Get rejection indices
ind_MLFDR = which(lfdr<= lfdrk)
ind_HDMT = which(pmax <= threshhold)
ind_DACT = which(p_dact_corr <= size)

venn_list = list(DACT = ind_DACT, HDMT = ind_HDMT, MLFDR = ind_MLFDR)
tcga_venn = ggvenn::ggvenn(venn_list,
               fill_color = c("#F8766D","#00BA38","#619CFF"),
               stroke_size = 0.5)+
  labs(title = "SNP-CPG-Gene Expression triplets identified by all methods")+
  theme(plot.title = element_text(hjust = 0.5))  
ggsave("TCGA_venn.png", plot = tcga_venn,bg = "white", width = 7, height = 7, dpi = 1200)  # Adjust width and height as needed




##Data #2, Hutch Mediation Data

load("Hutch_Mediation_output.RData")

apply(outdat, 2, function(x)sum(x>0 & x<=1))

apply(outdat, 2, function(x) sum(x>0))
p1 = outdat[,1]
p2 = outdat[,3]
input_pvalues = cbind(p1,p2)
pmax = apply(input_pvalues, 1, max)
nullprop = null_estimation(input_pvalues)
fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                         nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)
threshhold = max(pmax[fdr_hdmt<= size])
rej_HDMT = snp_data[pmax <= threshhold,] %>% arrange(pmax)