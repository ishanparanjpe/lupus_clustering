source("/load_data_meqtl.R")
additive.p<-function(snp, y, args){
  model<-glm(y~as.numeric(snp) + ., data = args,family=binomial)
  coef(summary(model))[2,4] # i.e. the p-value for association
}
setwd("/Users/iparanjpe/Documents/clues/integrative/meqtl/complete_clues_imputed/")

load(file="analysis/all_acr_clusters_linear_cis_meqtl_fdr_0.05.RData" )
### load clinical data 
load("analysis/all_clues_acr_data_with_kmeans_no_hem_no_imm_clusters.RData")
clinical_data<-output

clinical_data<- clinical_data[which(rownames(clinical_data) %in% metadata$subject_id),]
rownames(clinical_data)<-metadata$chip_id[match(rownames(clinical_data), metadata$subject_id)]
rownames(clinical_data)==rownames(bed)

### find SNPs directly associated with cluster
me_cis_sig$snps <- trimws(me_cis_sig$snps)
snps <- unique(me_cis_sig$snps)
cpgs <- unique(me_cis_sig$gene)

table(snps %in% colnames(bed))   ### check that all snps are in the snp data matrix

pvals <- sapply(snps, function(x){
  additive.p(bed[,x], clinical_data$clusters , covariates_meqtl)
  
})
sum(adj.pvals<0.05)
adj.pvals <- p.adjust(pvals, "fdr")
sig_snps <- names(pvals)[which(pvals<.05)]
sig_meqtl <- (me_cis_sig[which(me_cis_sig$snps %in% sig_snps),])

library(cit)
cit_results <- list()
for (i in 1:nrow(sig_meqtl)){
  

  L<- bed[,as.character(sig_meqtl$snps)[i]]
  G <- t(meth_data[as.character(sig_meqtl$gene)[i],])
  T <- clinical_data$clusters
  
  
  cit_results[[i]] <- cit.cp(L, G, T, n.perm = 20)
  
}
cit_results_fdr <- fdr.cit(cit_results)
sum(cit_results_fdr$q.cit<0.05)
cit_meqtl <-sig_meqtl[which(cit_results_fdr$q.cit < .05),]
cit_meqtl$cit_pval <- cit_results_fdr$p.cit[which(cit_results_fdr$q.cit < .05)]
cit_meqtl$cit_fdr <- cit_results_fdr$q.cit[which(cit_results_fdr$q.cit < .05)]
length(unique(cit_meqtl$snps))
length(unique(cit_meqtl$gene))
cit_meqtl[cit_meqtl$snp_gene=="GAB2",]



### Adjusted linear models
i <- 1
y <- clinical_data$clusters
snp_name <- as.character(cit_meqtl$snps)[i]
cpg_name <- as.character(cit_meqtl$gene)[i]
snp <- bed[,snp_name]

meth <- t(meth_data[cpg_name,])
meth_beta <-  t(data_beta[cpg_name,])
model_unadj <- glm(y~as.numeric(snp), family=binomial)

model_adj <- glm(y~as.numeric(snp) + meth, family=binomial)
coef(summary(model_unadj))[2,4] # i.e. the p-value for association
coef(summary(model_adj))[2,4] # i.e. the p-value for association
unadj_beta <- coef(summary(model_unadj))[2,1]
unadj_beta_std <- coef(summary(model_unadj))[2,2]
adj_beta <- coef(summary(model_adj))[2,1]
adj_beta_std <- coef(summary(model_adj))[2,2]


### Make plots
theme <-   theme(axis.title = element_text(size = 30),
                 plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
                 axis.text.x  = element_text(size = 30, colour = "black"),
                 axis.text.y =element_text(size = 30, colour ="black"),
                 legend.text = element_text(size = 35),
                 legend.title=element_text(face = "bold", size = 25),
                 title =element_text(size=20, face='bold'),
                 plot.title = element_text(size = 30,hjust = 0.5))


ref_allele <- dat$bim$V6[match(gsub("[.]",":", snp_name), dat$bim$V2)]
alt_allele <- dat$bim$V5[match(gsub("[.]",":", snp_name), dat$bim$V2)]

snp[snp==0] <- paste0(alt_allele, alt_allele)
snp[snp==1] <- paste0(ref_allele, alt_allele)
snp[snp==2] <- paste0(ref_allele, ref_allele)
### Remove all na samples
na_samples <- which(is.na(snp))
snp<- snp[-na_samples]
meth<- meth[-na_samples]
meth_beta<- meth_beta[-na_samples]
clusters <- clinical_data$clusters[-na_samples]

library(ggplot2)
ggplot(data= NULL, aes(x= clusters, y= meth_beta)) +geom_point(position = "jitter",size= 5)+geom_boxplot(width = .1)+
  ylab("% Methylation")+
  xlab("Cluster")+
  theme

ggplot(data= NULL, aes(x= snp, y= meth_beta))+geom_point(position = "jitter", size= 5)+geom_boxplot(width = .1)+
  ylab("% Methylation")+
  xlab("SNP")+
  theme


snp_df <- data.frame(snp=snp, cluster = clusters)
library(reshape)
melted_data <-plyr::count(snp_df)
melted_data <- melted_data[!is.na(melted_data$snp),]
require(dplyr)    
freq_sum <- melted_data %>% group_by(cluster) %>%summarise(sum = sum(freq))

melted_data$percent <- signif(100*melted_data$freq/freq_sum$sum, 2)

ggplot(melted_data, aes(fill=factor(snp), y= percent, x=cluster)) + 
  geom_bar( stat="identity")+
  ylab("% Patients in Cluster")+
  xlab("Cluster")+
  guides(fill= guide_legend("SNP"))+
  theme




coef_df <- data.frame(Model = c("Unadjusted", "Adjusted for M"), Beta = c(unadj_beta, adj_beta), std = c(unadj_beta_std, adj_beta_std))
ggplot(data= coef_df , aes(x= Model, y = Beta))+geom_point(size = 14, fill = "white", shape = 21)+
  geom_errorbar(aes(ymin = Beta - std, ymax = Beta+std, width = .3, size = 2))+
  geom_hline(yintercept = 0, linetype="dashed", size = 3)+
  guides(size = F)+
  theme


