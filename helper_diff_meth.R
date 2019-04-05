library(VennDiagram)
library(ggplot2)
library(calibrate)
makeVenn<- function(itemList, nameList,filename, title){
  # list1<- list1[!is.na(list1)]
  # list2<-list2[!is.na(list2)]
  # list3<-list3[!is.na(list3)]
  plotList<-itemList
  print(plotList)
  (venn.diagram(x=plotList,category.names = nameList,main = title,
                
                lty = 'blank',
                fill= c("green","red"),
                cex = 2.4,
                main.cex=2,
                fontface = "bold",
                fontfamily = "sans",
                cat.cex = 1.3,
                cat.fontface = "bold",
                alpha = c(0.5,0.5),
                cat.fontfamily = "sans",
                width = 5000,
                height = 5000,
                
                filename= filename
                
                
  ))
}

volcano <- function(limma_table, FC_cutoff, padj_cutoff, filename){
  # limma_table<- limma_table[cpgs,]
  png(filename, width = 10, height = 10)
  # Make a basic volcano plot
  with(limma_table, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2), ylim = c(4,9)))
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(limma_table, adj.P.Val<padj_cutoff), points(logFC, -log10(P.Value), pch=20, col="red"))
  with(subset(limma_table, abs(logFC)>1 & adj.P.Val<padj_cutoff), points(logFC, -log10(P.Value), pch=20, col="blue"))
  #with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

  ### add gene names for significant CpGs
  with(subset(limma_table, abs(logFC)>FC_cutoff & adj.P.Val<padj_cutoff), textxy(logFC, -log10(P.Value), labs=cpg_gene, cex=.8))
  dev.off()
}
volcano_contrast <- function(limma_table,group1, group2, delta_beta_cutoff, cpgs, padj_cutoff, filename){
  # limma_table<- limma_table[cpgs,]
  limma_table$delta_beta <- limma_table[,group1]-limma_table[,group2]
  
  
  png(filename, width = 15, height = 15, units = "in", res = 200)
  par(ps = 12, cex = 1.8, cex.axis = 2, cex.lab = 4)
  # Make a basic volcano plot
  with(limma_table, plot(delta_beta, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-0.4,0.4), xlab= "", ylab = "", ylim = c(0,9)))
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(limma_table[cpgs,], adj.P.Val<padj_cutoff), points(delta_beta, -log10(P.Value), pch=20, col="red"))
  with(subset(limma_table[cpgs,], abs(delta_beta)>delta_beta_cutoff & adj.P.Val<padj_cutoff), points(delta_beta, -log10(P.Value), pch=20, col="blue"))
  #with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  
  ### add gene names for significant CpGs
  # with(subset(limma_table[cpgs,], abs(delta_beta)>delta_beta_cutoff & adj.P.Val<padj_cutoff), textxy(delta_beta, -log10(P.Value), labs=cpg_gene, cex=.8))
  num_sig <- dim(subset(limma_table[cpgs,], abs(delta_beta)>delta_beta_cutoff & adj.P.Val<padj_cutoff))[1]
  num_sig <- ifelse(num_sig<10, num_sig, 10)
  # with(subset(limma_table[cpgs,], abs(delta_beta)>delta_beta_cutoff & adj.P.Val<padj_cutoff)[1:num_sig,], textxy(delta_beta, -log10(P.Value), labs=cpg_gene, cex=.8))
  with(subset(limma_table[cpgs,], abs(delta_beta)>delta_beta_cutoff & adj.P.Val<padj_cutoff), textxy(delta_beta, -log10(P.Value), labs=cpg_gene, cex=.8))
  
    dev.off()
}
