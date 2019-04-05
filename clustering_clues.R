library(clustMixType)
library(Rtsne)
library(magrittr)
library(cluster)
library(dplyr)
library(Rtsne)
library(compareGroups)
library(ggfortify)
library(ggplot2)
library(plyr)
library(mice)
library(klaR)
library(missForest)
library(vegan)
library(FactoMineR)

### Read clinical data
data<-read.csv("./cluestime1_updated.csv", sep="\t")
all_data<-data

### Convert variable names
variable_conversion<-read.csv("../clinical analysis/variable_conversion.csv")
clinical_vars <- read.csv("ACR_clustering_vars.csv", stringsAsFactors = F)


colnames(data)<-variable_conversion$new[match(colnames(data),variable_conversion$original)]
data<- data[,clinical_vars$x]
data <- data.frame(lapply(data,as.factor), check.names = F)


which(is.na(data), arr.ind = T) # check which values are NA 
data<-data[-313,]
all_data<-all_data[-313,]
rownames(data)<-all_data$subjectid
rownames(all_data)<-all_data$subjectid

### remove ACR from variable names 
colnames(data)<- unlist(unname(sapply(colnames(data), function(x) strsplit(x, "ACR ")[[1]][2])))

data <- data.frame(data, check.names = F)
all_data <- data.frame(all_data, check.names = F)


###########################
##### Clustering with PCAMix
###########################
library(Rtsne)
# data[is.na(data)]<-0
library(PCAmixdata)
split <- splitmix(data)
X1 <- split$X.quanti 
X2 <- split$X.quali 

res.pcamix <- PCAmix(X.quanti=X1, X.quali=X2,rename.level=TRUE,
                     graph=FALSE, ndim = 200)

loadings <- res.pcamix$ind$coord

### check number of MCA components to retain using missMDA with K-fold CV
library(missMDA)
num_comp <- estim_ncpMCA(data, method.cv = "Kfold", ncp.max = 10, ncp.min = 1)
num_comp


### Performe clustering and evaluate cluster significance with clusterboot
library(fpc)
clust_boot <- clusterboot(loadings[,1:2],clustermethod = kmeansCBI,krange=3, seed = 1)


# rename clusters to match Mild (M), Severe1 (S1), Severe2 (S2)
clusters[clusters==1]<- "S1"
clusters[clusters==2]<-"M"
clusters[clusters==3]<- "S2"
stability_scores <- data.frame(t(clust_boot$bootresult))
colnames(stability_scores)<- c("Cluster S1", "Cluster M","Cluster S2")
stability_scores<- melt(stability_scores)

###########################
###### Compare each variable across clusters 
###### Chi squared tests for variables that are factors (categorical)
###### non parametric kruskal wallist test for numeric variables
###########################


output<-cbind(data, clusters = clusters, subject_id= all_data$subjectid)
save(output, file = "/Users/iparanjpe/Documents/clues/clinical analysis/clinical_analysis_output/all_clues_acr_data_with_kmeans_no_hem_no_imm_clusters.RData")
data_race <- cbind(data, race=factor(all_data$raceeth), lsi = all_data$lupusseverityindex, SLEDAI = all_data$sledaiscore, Female = factor(all_data$female))
temp <- cbind(data_race, clusters)
save(temp,file = "/Users/iparanjpe/Documents/clues/clinical analysis/clinical_analysis_output/all_clues_acr_data_with_kmeans_no_hem_no_imm_clusters_with_race.RData")
pvals<-c()
for (i in 1:ncol(data_race)){ 
  print(i)
  if (class(data_race[,i])=="factor"){
    print("factor")
    #pvals[i] <-cv.test(data_imputed[,i], data_imputed$raceeth)
    pvals[i]<-chisq.test(data_race[,i], clusters)$p.value
    
    
  }
  else{
    print("numeric")
    aov_test <- anova(aov(data_race[,i]~clusters, data = NULL))
    # kw_test <- kruskal.test(data_race[,i]~clusters, data = NULL)
    pvals[i] <- aov_test$`Pr(>F)`[1]
  }
}
pvals <- pvals
adj.p <- p.adjust(pvals, method ="fdr")

significant_vars <- data.frame(colnames(data_race)[adj.p<1], pvals[pvals<=1],adj.p[adj.p<=1])
colnames(significant_vars)<-c("Variable","pval", "Adj-p-val")
significant_vars<-significant_vars[order(significant_vars$`Adj-p-val`),]
significant_vars$Variable<-as.character(significant_vars$Variable)


####### 
### Train random forest classifier 
#######

library(randomForest)

### Remove NA values from data
data_no_na <- data
data_no_na[is.na(data_no_na)]<-0
oob_errors <- c()
set.seed(1)

### Find optimal number of variables at each split(mtry) by minizming out of bag (oob) error
for(mtry in seq(1,dim(data_no_na)[2] )){
  model<- randomForest(x = data_no_na, y = factor(clusters),classwt = 1/table(clusters), mtry = mtry, ntree =5000,strata = factor(clusters))
  oob_errors[[mtry]]<- model$err.rate[500,1]
}
best_mtry <- which.min(oob_errors)

### Plot mtry vs OOB
ggplot(data=NULL, aes(x= seq(1,dim(data_no_na)[2]),y = oob_errors))+geom_point(size = 8)+
  guides(fill = FALSE)+
  scale_x_continuous(name = "Num. Variables at each node\n (mtry)")+
  scale_y_continuous(name = "OOB Misclassification Error")+
  theme(axis.title = element_text(size = 25),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        axis.text = element_text(size = 20, colour = "red"))



model<- randomForest(x = data_no_na, y = factor(clusters),classwt = 1/table(clusters), mtry = best_mtry, ntree = 5000,strata = factor(clusters))
model
save(model, file="all_clues_random_forest.RData")


