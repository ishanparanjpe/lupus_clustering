######
#### preload data
######
load("data/chip_metadata.RData")
load("data/clues_metadata.RData")
load("dataepic_anno.RData")
load("data/data_m_no_sex.RData")
load("data/data_beta_no_sex.RData")

covariates<-read.csv("data/data_cell.samples.csv")
load("data/genetic_pc.RData")
load("data/biomart_human_gene_anno.RData")
load("data/affy_6.RData")

drug_data <- read.csv("data/common_meds_pca.csv", row.names  = 1, stringsAsFactors = F)



meth_data<- data_m
rownames(covariates)<-covariates$sampleid
metadata<- cbind(clues_metadata,chip_metadata )
rownames(metadata)<-metadata$chip_id

drug_data <- drug_data[which(rownames(drug_data) %in%metadata$subject_id ),]

rownames(output_pca)<-metadata$chip_id[match(rownames(output_pca), metadata$subject_id)]
rownames(drug_data)<-metadata$chip_id[match(rownames(drug_data), metadata$subject_id)]


### find samples which are common between all the data files
common_subjects<-   intersect(rownames(drug_data), intersect(rownames(output_pca), 
                                                             intersect(covariates$sampleid, intersect(metadata$chip_id, colnames(meth_data)))))


#### Only take the subset of samples that have methylation, genetic, and clinical data
metadata<- metadata[common_subjects,]
meth_data<- meth_data[,common_subjects]
rownames(covariates)<-covariates$sampleid
covariates<- covariates[common_subjects,]
output_pca <- output_pca[common_subjects, ]
drug_data <- drug_data[common_subjects,]

### read in imputed genotyping data
library(plink2R)
dat<-read_plink("data/final_merged_ld_pruned")

bed <-data.frame(dat$bed)
##### convert genotyping sample names
sample_names <- unname(unlist(sapply(rownames(dat$bed), function(x) strsplit(x, "_")[[1]][1])))
gene_id_conversion <- read.csv("data/clues-IDmap.txt",sep = "\t")

gene_id_conversion$SLE_GeneticsStudyID<-as.character(gene_id_conversion$SLE_GeneticsStudyID)

genID <- match(sample_names[181:282], gene_id_conversion$SLE_GeneticsStudyID)
sample_names[181:282]<-gene_id_conversion[genID,"SubjectID"]
rownames(bed)<-sample_names
bed<-bed[which(rownames(bed) %in% metadata$subjectid),]
rownames(bed)<-metadata$chip_id[match(rownames(bed), metadata$subject_id)]


#filter genotyping data to match the subset which has methylation data
snp_meth_common<- intersect(rownames(metadata), rownames(bed))

metadata<- metadata[snp_meth_common,]
meth_data<- meth_data[,snp_meth_common]
covariates<- covariates[snp_meth_common,]
output_pca <- output_pca[snp_meth_common, ]
bed<-bed[snp_meth_common,]
drug_data <- drug_data[snp_meth_common,]


### remove methylation probes with cpgs
probes_with_snps <- read.table("data/snp.in.probe.txt", header = F,stringsAsFactors = F)
meth_data<- meth_data[which(!rownames(meth_data) %in% probes_with_snps[,1]),]

### get SNP names (colnames of bed)
snp_names <- colnames(bed)
snp_names <- sapply(snp_names, function(x) strsplit(x,"X")[[1]][2])
snp_names <- unname(unlist(snp_names))

chr <- unname(unlist(sapply(snp_names, function(x) strsplit(x, "[.]")[[1]][1])))
pos <- unname(unlist(sapply(snp_names, function(x) strsplit(x, "[.]")[[1]][2])))

colnames(bed) <- as.character(snp_names)
## create dataframes of snp location and meth location
snp_loc <- data.frame(snp = as.character(snp_names), chr = paste0('chr', as.character(chr)), pos = as.integer(pos))
# save(snp_loc, file="/Users/iparanjpe/Documents/clues/integrative/meqtl/snp_loc.RData")
meth_loc<- data.frame(cpg= as.character(epic_anno$Name), chr= paste0('chr', as.character(epic_anno$CHR)), s1 = epic_anno$MAPINFO, s2 = epic_anno$MAPINFO+1) 
meth_loc <- meth_loc[which(!is.na(meth_loc[,3])),]
#save(meth_loc, file = "/Users/iparanjpe/Documents/clues/integrative/meqtl/meth_loc.RData")




### remove methylation probes on sex chromosomes
sex_probes <- meth_loc$cpg[which(meth_loc$chr %in% c("chrX","chrY"))]
meth_data<-meth_data[which(!rownames(meth_data) %in% sex_probes),]

metadata$smokeever[is.na(metadata$smokeever)]<-0
steroids <- metadata$steroidsoral
steroids[is.na(steroids)]<-0

### create covariates matrix for meqtl

covariates_meqtl<-(data.frame(metadata$female, metadata$age, covariates$rc1, covariates$rc2,
                              covariates$rc3, output_pca[,1:3],
                              (metadata$smokeever), (metadata$alcohol), drug_data[,1:3]))


rownames(covariates_meqtl)<-colnames(meth_data)


data_beta <- data_beta[rownames(meth_data), colnames(meth_data)]


