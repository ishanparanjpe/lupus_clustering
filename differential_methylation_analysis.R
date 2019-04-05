############ Load data
#### Clinical data with K-Means clusters, Covariates, EPIC annotation, 
####  chip metadata, methylation data after Combat
############ 
load("data/chip_metadata.RData")
load("data/clues_metadata.RData")

### Methylation EPIC chip annotation file
load("data/epic_anno.RData") 

### Methylation M values with CpGs in rows and subject IDs in columns
load("data/data_m_no_sex.RData")

### Top principal components (PCs) from drug matrix
### rows are sample IDs and columns are PCs
drug_data <- read.csv("data/common_meds_pca.csv", row.names  = 1, stringsAsFactors = F)

covariates<-read.csv("data/data_cell.samples.csv")
load("data/genetic_pc.RData")

meth_data<- data_m
rownames(covariates)<-covariates$sampleid
metadata<- cbind(clues_metadata,chip_metadata )
rownames(metadata)<-metadata$chip_id
drug_data <- drug_data[which(rownames(drug_data) %in%metadata$subject_id ),]

rownames(output_pca)<-metadata$chip_id[match(rownames(output_pca), metadata$subject_id)]
rownames(drug_data)<-metadata$chip_id[match(rownames(drug_data), metadata$subject_id)]


### load clusters from clustering ACR data
load("analysis/all_clues_acr_data_with_kmeans_no_hem_no_imm_clusters.RData")

clinical_data<-output

clinical_data<- clinical_data[which(rownames(clinical_data) %in% metadata$subject_id),]
rownames(clinical_data)<-metadata$chip_id[match(rownames(clinical_data), metadata$subject_id)]

### find samples which are common between all the data files
common_subjects<-  intersect(rownames(drug_data), intersect(rownames(clinical_data),intersect(rownames(output_pca), 
                                                                                              intersect(covariates$sampleid, intersect(metadata$chip_id, colnames(meth_data))))))

#### Only take the subset of samples that have methylation, genetic, and clinical data

metadata<- metadata[common_subjects,]
meth_data<- meth_data[,common_subjects]
rownames(covariates)<-covariates$sampleid
covariates<- covariates[common_subjects,]
output_pca <- output_pca[common_subjects, ]
clinical_data<-clinical_data[common_subjects,]

drug_data <- drug_data[common_subjects,]

all(rownames(clinical_data)== rownames(metadata))


############ 
####   Differential Methylation Analysis with LIMMA
############

library(limma)

metadata$smokeever[is.na(metadata$smokeever)]<-0
steroids <- metadata$steroidsoral
steroids[is.na(steroids)]<-0
design <- model.matrix(~0+factor(clinical_data$clusters)+factor(covariates$female)+factor(metadata$smokeever)+
                         covariates$age+covariates$rc1+covariates$rc2+
                         covariates$rc3+output_pca[,1]+output_pca[,2]+output_pca[,3]+factor(metadata$alcohol) +drug_data[,1]+drug_data[,2]+drug_data[,3]
)


colnames(design)<-c("clusterM", "clusterS1","clusterS2","female","smoke","age","rc1","rc2","rc3","geno1","geno2","geno3","alcohol", paste0("drug_data",seq(1,3)))
head(design)

fit <-lmFit(meth_data, design)
contrasts<-makeContrasts(clusterS1-clusterM, clusterS2 - clusterS1, clusterS2 - clusterM,levels = design)


fit2<-contrasts.fit(fit, contrasts)
fit2<- eBayes(fit2)

### Determine cluster-associated CpGs (FDR<0.1) using F test
ttable <- topTable(fit2, coef =1:3, n= Inf, p.value = 0.1)

#### Use nestedF to determine which clusters contribute to significant F statistic
results<- decideTests(fit2,p.value = 0.1, lfc = 0, method="nestedF")
cluster_S1_M <- names(which(results[,1]!=0))
cluster_S2_S1 <- names(which(results[,2]!=0))
cluster_S2_M<- names(which(results[,3]!=0))


#######
## Annotate results
#######

### get cpg annotations
cpg <- as.character(rownames(ttable))

### cpg associated gene name
cpg_gene <- epic_anno[match(cpg, epic_anno$IlmnID),"UCSC_RefGene_Name"]
cpg_gene <- as.character(cpg_gene)
cpg_gene <- sapply(cpg_gene, FUN=function(x) strsplit(x,split=";")[[1]][1] )
cpg_gene <- unname(cpg_gene)

### position of cpg relative to gene
cpg_position<-epic_anno[match(cpg, epic_anno$IlmnID),"UCSC_RefGene_Group"]
cpg_position<-as.character(cpg_position)
cpg_position <- sapply(cpg_position, FUN=function(x) strsplit(x,split=";")[[1]][1] )
cpg_position <- unname(cpg_position)
### decription of cpg gene
description <- gene_anno$description[match(cpg_gene, gene_anno$external_gene_name)]
ttable$cpg_gene<-cpg_gene
ttable$cpg_description <- description
ttable$cpg_position <- cpg_position

### look at mean beta values for each cpg in each cluster
beta_values <- data.frame(matrix(nrow= dim(ttable)[1], ncol = 3))
rownames(beta_values)<-rownames(ttable)
colnames(beta_values)<-c("clusterM","clusterS1","clusterS2")
beta_values[,1]<- apply((data_beta[rownames(beta_values), which(clinical_data$clusters=="M"),]),
                        1, function(x) mean(x, na.rm = T))
beta_values[,2]<- apply((data_beta[rownames(beta_values), which(clinical_data$clusters=="S1"),]),
                        1, function(x) mean(x, na.rm = T))
beta_values[,3]<- apply((data_beta[rownames(beta_values), which(clinical_data$clusters=="S2"),]),
                        1, function(x) mean(x, na.rm = T))
beta_values$variance <- apply(data_beta[rownames(beta_values),], 1, function(x) var(x, na.rm = T))
beta_values$clusterS2_hypo<- (beta_values[,3]<beta_values[,2] & beta_values[,3]<beta_values[,1])  ### cluster 3 hypomethylated CpGs
beta_values$clusterS2_hyper<- (beta_values[,3]>beta_values[,2] & beta_values[,3]>beta_values[,1])

beta_values$clusterS1_hypo<- (beta_values[,2]<beta_values[,1] & beta_values[,2]<beta_values[,3])  ### cluster 2 hypomethylated CpGs
beta_values$clusterS1_hyper<- (beta_values[,2]>beta_values[,1] & beta_values[,2]>beta_values[,3])

beta_values$clusterM_hypo<- (beta_values[,1]<beta_values[,2] & beta_values[,1]<beta_values[,3])  ### cluster 2 hypomethylated CpGs
beta_values$clusterM_hyper<- (beta_values[,1]>beta_values[,2] & beta_values[,1]>beta_values[,3])

### check overlap with IFN-alpha and IFN-gamma genes
ifn_alpha_genes <- read.table("reference_gene_sets/HALLMARK_INTERFERON_ALPHA_RESPONSE.txt", skip = 2, header = F)
beta_values$ifn_alpha<- ttable$cpg_gene %in% ifn_alpha_genes$V1

ifn_gamma_genes <- read.table("reference_gene_sets/HALLMARK_INTERFERON_GAMMA_RESPONSE.txt", skip = 2, header = F)
beta_values$ifn_gamma<- ttable$cpg_gene %in% ifn_gamma_genes$V1

ttable_beta<-cbind(ttable, beta_values)

### add results of nestedF test to ttable_beta
ttable_beta$clusterS2_S1_sig <- ifelse(rownames(ttable_beta) %in% cluster_S2_S1, 1, 0)
ttable_beta$clusterS2_M_sig <- ifelse(rownames(ttable_beta) %in% cluster_S2_M, 1, 0)
ttable_beta$clusterS1_M_sig <- ifelse(rownames(ttable_beta) %in% cluster_S1_M, 1, 0)


########
###  Make volcano plot of pairwise comparisons
########
source("helper_diff_meth.R")
delta_beta_cutoff <- 0.05
filename = "acr_clinical_clusters_S2_M_volcano.png"
volcano_contrast(ttable_beta,"clusterS2","clusterM",cpgs = cluster_S2_M, delta_beta_cutoff  = 0.1, padj_cutoff = 0.1, filename=filename)


sample_ <- sample(nrow(ttable_beta),size = 100000 )
delta_beta <- ttable_beta$clusterS2[sample_ ]-ttable_beta$clusterM[sample_ ] 
var_beta <-ttable_beta$variance[sample_ ] 
p <- -log10(ttable_beta$P.Value[sample_ ])
plot((delta_beta), p )
plot(density(var_beta))

### Pathway Analysis:
library(clusterProfiler)
library(ReactomePA)
entrez_genes = bitr(ttable$cpg_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
enriched_pathways<-enrichPathway(gene=entrez_genes$ENTREZID,pvalueCutoff=0.05, readable=T)
x<-barplot(enriched_pathways, showCategory = 10, font.size = 20)
descriptions <- as.character(x$data$Description)
wrapped_text <-list()
for (i in 1:length(descriptions)){
  str1<- descriptions[i]
  # wrapped_text[[i]]<-paste(sapply(seq(1, nchar(str1), 16), function(i) paste0(substring(str1, i, min(i + 15, nchar(str1))), '\n')), collapse='')
  wrapped_text[[i]]<-paste(strwrap(str1,width=40), collapse = "\n")
  
  
}
x$data$Description<- unlist(wrapped_text)
x$data$Description<- factor(x$data$Description, levels = x$data$Description[order(-x$data$p.adjust)])
x
library(RColorBrewer)
library(ggplot2)
ggplot(data=x$data, aes(x = Description, y = -log10(p.adjust), fill =-log10(p.adjust) ))+ geom_bar(stat='identity' ) +
  coord_flip()+
  scale_x_discrete(name = "Pathway")+
  scale_y_continuous(name = "-log10(q)")+
  scale_fill_gradient(low = "purple", high = "red", guide = "none")+
  theme(axis.title = element_text(size = 45),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        axis.text.x  = element_text(size = 25, colour = "blue"),
        axis.text.y =element_text(size = 35, colour ="blue"),
        legend.text = element_text(size = 20),
        legend.title=element_text(face = "bold", size = 25),
        title =element_text(size=20, face='bold'),
        plot.title = element_text(size = 30,hjust = 0.5))



### Heatmap of significant CpGs 
library(pheatmap)
library(plyr)

sig_meth_data <- data_beta[rownames(ttable)[which(ttable$adj.P.Val<0.1)],]
race<-revalue(as.character(metadata$raceeth), c('1'='White','2'='Hispanic','3'='African American','4'='Asian', '5'="Other"))


anno <- data.frame(Cluster= factor(clinical_data$clusters ), Race= factor(race=="White", labels = c("Non-White", "White")))
rownames(anno)<- colnames(sig_meth_data)
anno_colors <- list(Cluster = c('S1'="darkgoldenrod1",'M'="seagreen2",'S2'="royalblue1"),
                    Race = c("White"="red2","Non-White"="gold"))

anno<- anno[order(anno$Cluster, decreasing = F),]
anno_rows <- data.frame("IFN-A Responsive" = factor(as.numeric(ttable_beta$ifn_alpha), labels  = c("No","Yes")), check.names = F
)
rownames(anno_rows)<-rownames(ttable)

legend_titles <- c("Race","Cluster","IFN-alpha Responsive")

library(RColorBrewer)
hmcol<-brewer.pal(11,"RdBu")
cols <-  colorRampPalette(rev(brewer.pal(n = 11, name =
                                           "RdYlBu")))(100)
scaled_data <- t(scale(t(sig_meth_data[, rownames(anno)])))
scaled_data[scaled_data>3]<-3
scaled_data[scaled_data< -3]<- -3

pheatmap(scaled_data, annotation = anno,color = cols, cluster_cols = F, cluster_rows = T,annotation_colors = anno_colors,
         scale = "none", show_rownames = F, show_colnames = F, annotation_row = anno_rows, annotation_names_row = T,border_color = NA, legend = F,annotation_legend = T )


#########
#####   QQ plot
########
library(qqman)
qqplot <- qq(as.numeric(ttable$P.Value)[1:10000],pch = 18, col = "blue4", cex = 1, las = 1, cex.axis = 2 )
qqplot
pvals <- as.numeric(ttable$P.Value)

observed = -log10(sort(pvals,decreasing=FALSE))
expected = -log10( ppoints(length(pvals) ))

library(ggplot2)
ggplot(data = NULL , aes(x= expected, y = observed ))+geom_point(size = 3, color = "blue")+
  xlab(bquote('Expected'~ -log[10] ~ '(p-value)'))+
  ylab(bquote('Observed'~ -log[10] ~ '(p-value)'))+
  # geom_abline(mapping = NULL, data = NULL, 1,
  #             na.rm = FALSE, show.legend = NA)+
  
  geom_abline(slope = 1, linetype = "dashed", size = 2, color = "red")+
  theme(axis.title = element_text(size = 45),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        axis.text.x  = element_text(size = 35, colour = "blue"),
        axis.text.y =element_text(size = 35, colour ="blue"),
        legend.text = element_text(size = 20),
        legend.title=element_text(face = "bold", size = 25),
        title =element_text(size=20, face='bold'),
        plot.title = element_text(size = 30,hjust = 0.5))

plot((observed), (expected))

#######
### Enrichment of race associated CpGs
#######
cluster_cpg <- rownames(ttable)[which(ttable$adj.P.Val<0.1)]


sample_size <- dim(meth_data)[2]
num_iter <- 1000

#### for each iteration randomly shuffle self-reported race and determine number of race-associated CpGs
results <- sapply(seq(1, num_iter), function(iter){
  race_pvals <-double()
  # sample <- sample(rownames(meth_data), sample_size)
  sample <- sample(colnames(meth_data), sample_size)
  race_pvals <- apply(meth_data[cluster_cpg, ],1, function(x){
    metadata_shuffle <- metadata[sample,]
    covariates_shuffle <- covariates[sample,]
    output_pca_shuffle <- output_pca[sample,]
    drug_data_shuffle <- drug_data[sample,]
    model <- aov(as.numeric(x) ~factor(metadata_shuffle$raceeth)+factor(covariates_shuffle$female)+factor(metadata_shuffle$smokeever)+covariates_shuffle$age+
                   covariates_shuffle$rc1+covariates_shuffle$rc2+covariates_shuffle$rc3+ output_pca_shuffle[,1]+output_pca_shuffle[,2]+output_pca_shuffle[,3]+ 
                   factor(metadata_shuffle$alcohol)+drug_data_shuffle[,1]+drug_data_shuffle[,2]+drug_data_shuffle[,3] )
    summary(model)[[1]][1,5]
    
  } )
  
  
})

results_summary <- apply(results, MARGIN = 2,FUN = function(x){
  sum(x<0.05)
})
plot(density(results_summary))

#### Check how many cluster associated CpGs are also associated with self-reported race
race_pvals_clues <- apply(meth_data[cluster_cpg, ],1, function(x){
  model <- aov(as.numeric(x) ~factor(metadata$raceeth)+factor(covariates$female)+factor(metadata$smokeever)+covariates$age+
                 covariates$rc1+covariates$rc2+covariates$rc3+ output_pca[,1]+output_pca[,2]+output_pca[,3]+ 
                 factor(metadata$alcohol)+drug_data[,1]+drug_data[,2]+drug_data[,3] )
  summary(model)[[1]][1,5]
  
} )

enrichment_pval <- sum(results_summary >=num_sig_race)/num_iter