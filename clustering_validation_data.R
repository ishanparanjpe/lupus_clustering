### Load validation data
data <- read.csv("data/validation_clinical_data.csv", stringsAsFactors = F)

### remove samples that overlap with CLUES
overlap <- read.csv("data/clues_validation_overlap.csv")

### Convert variable names for consistency
variable_conversion <- read.csv("data/variable_conversion_validation.csv", stringsAsFactors = F,
                                sep= "\t",header = 1)
data_acr<- data[, variable_conversion$Vars]

colnames(data_acr)<-variable_conversion$NewVars[match(colnames(data_acr), variable_conversion$Vars)]
data_acr$`ACR Photosensitivity`<- data$PhotosensitivityMD |data$PhotosensitivityS
data_acr$`ACR APLA` <- data$CLIP_MD | data$RVVT_MD | data$FalseSyphilisMD

clinical_vars <- read.csv("data/ACR_clustering_vars.csv", stringsAsFactors = F)
### there was a misspelling in the CLUES data so to keep variable names consistent, change the clinical_vars name
clinical_vars$x[clinical_vars$x=='ACR Lymphopenia'] <- 'ACR Lymphoenia'
data_acr <- data_acr[,clinical_vars$x]
for (i in 1:ncol(data_acr)){
  data_acr[,i]<-(as.numeric(data_acr[,i]) ) 
}

head(data_acr)

### assign clusters based on random forest model trained on CLUES data
library(randomForest)
load("analysis/all_clues_random_forest.RData")
clusters <- predict(model, newdata = data_acr)

clusters[clusters==2]<- "S1"
clusters[clusters==3]<-"M"
clusters[clusters==1]<- "S2"
output<-cbind(data_acr, clusters = clusters, subject_id= data$subjectid)

###########################
###### Compare each variable across clusters 
###### Chi squared tests for variables that are factors (categorical)
###### non parametric kruskal wallist test for numeric variables
###########################
pvals<-c()
for (i in 1:ncol(data_acr)){
  data_acr[,i]<-as.factor(as.numeric(data_acr[,i]) ) 
}


for (i in 1:ncol(data_acr)){ 
  print(i)
  if (class(data_acr[,i])=="factor"){
    print("factor")
    #pvals[i] <-cv.test(data_imputed[,i], data_imputed$raceeth)
    pvals[i]<-chisq.test(data_acr[,i], clusters)$p.value
    
    
  }
  else{
    print("numeric")
    #aov_test <- anova(aov(data_imputed[,i]~raceeth, data = data_imputed))
    kw_test <- kruskal.test(data_acr[,i]~clusters, data = clusters)
    pvals[i] <- kw_test$p.value
  }
}
pvals <- pvals
adj.p <- p.adjust(pvals, method ="fdr")

significant_vars <- data.frame(colnames(data_acr)[adj.p<1], pvals[adj.p<1], adj.p[adj.p<1])
colnames(significant_vars)<-c("Variable","pval", "Adj-p-val")
significant_vars<-significant_vars[order(significant_vars$`Adj-p-val`),]
significant_vars$Variable<-as.character(significant_vars$Variable)



