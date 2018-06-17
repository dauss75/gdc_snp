maf_to_matrix<-function(maf_dir_path, clinical_data, data_matrix){
  # data_preprocessing(raw_dir_path, clinical_data, data)
  #data_matrix<-file.path(data, "gdc_tcga")
  #data_matrix<-data
  dir.create(data_matrix,showWarnings = FALSE)
  maf_files<-list.files(maf_dir_path, full.names = TRUE)
  
  clinical<-fread(clinical_data)
  colnames(clinical)[2]<-"Tumor_Sample_Barcode"
  
  # create an aggregated data matrix
  y<-NULL; y_sample<-NULL
  
  for (i in 1:length(maf_files)){
    laml<-read.maf(maf=maf_files[i],isTCGA = FALSE, clinicalData=clinical,removeDuplicatedVariants=F)
    snp_cnt<-mutCountMatrix(maf=laml, includeSyn=FALSE, countOnly = "Missense_Mutation")
    laml@data[["Hugo_Symbol"]]<-paste(laml@data[["Hugo_Symbol"]],":",laml@data[["HGVSp_Short"]],sep="")
    snp_aa<-mutCountMatrix(laml, includeSyn = FALSE, countOnly = "Missense_Mutation")
    y<-c(y,rep(i-1,ncol(snp_aa))); y_sample<-c(y_sample,colnames(snp_aa))
    
    sample_name<-paste(unlist(strsplit(basename(maf_files[i]),'[.]'))[1:2],collapse="-")
    fwrite(snp_cnt,paste(data_matrix,"/",sample_name,".snp_cnt.txt",sep=""),row.names = TRUE, col.names=TRUE)
    fwrite(snp_aa,paste(data_matrix,"/",sample_name,".snp_gene_aa.txt",sep=""),row.names = TRUE, col.names=TRUE)
    rm(laml,snp_cnt,snp_aa)
  }
  y<-as.data.frame(y); y_sample<-as.data.frame(y_sample)
  fwrite(y,paste(data_matrix,"/y_label.txt",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
  fwrite(y_sample,paste(data_matrix,"/y_samplename.txt",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
}

gene_devide_exon_size<-function(exprs){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
  colnames(exonic.gene.sizes)[1]<-"size"
  exonic.gene.sizes$genesymbol<-getSYMBOL(na.omit(rownames(exonic.gene.sizes)), data='org.Hs.eg')
  exonic.gene.sizes<-na.omit(exonic.gene.sizes)
  rownames(exonic.gene.sizes)<-exonic.gene.sizes$genesymbol
  
  # y<-as.data.frame(exprs[,1]); 
  x<-as.data.frame(exprs); rownames(x)<-rownames(exprs)
  
  com_gene<-intersect(rownames(x),rownames(exonic.gene.sizes))
  x2<-x[com_gene,]; x2<-na.omit(x2)
  exonic.gene.sizes2<-exonic.gene.sizes[com_gene,] # sort by the common gene
  exprs_devided_by_exon_size<-x2/exonic.gene.sizes2$size
  exprs_devided_by_exon_size
}

combine_matrix<-function(data_path, gene=NULL){
  library(reshape2)
  snp_cnt_files<-list.files(data_path, pattern = "_cnt.txt$", full.names = TRUE)
  y_label<-fread(paste(data_path, "/y_label.txt",sep=""))
  if (is.null(gene)){
    snp_cnt_data <- Reduce(function(x, y) merge(x, y, all=TRUE), lapply(snp_cnt_files,fread)); 
    row.names(snp_cnt_data)<-snp_cnt_data$V1;snp_cnt_data$V1<-NULL
    snp_cnt_data[is.na(snp_cnt_data)] <- 0
    fwrite(snp_cnt_data,paste(data,"/snp_cnt_merged.txt",sep=""), row.names=T, col.names = T);
  } else {
    gene_list<-fread(gene); 
    snp_cnt_data<-gene_list
    for (i in 1:length(snp_cnt_files)){
      t_data<-as.data.frame(fread(snp_cnt_files[i])); rownames(t_data)<-t_data$V1; t_data$V1<-NULL
      t_data2<-t_data[gene_list$HGNC_SYMBOL,]; t_data2[is.na(t_data2)] <- 0
      snp_cnt_data<-cbind(snp_cnt_data,t_data2)
    }
    rm(t_data,t_data2)
    rownames(snp_cnt_data)<-snp_cnt_data$HGNC_SYMBOL;snp_cnt_data$HGNC_SYMBOL<-NULL
    
    snp_cnt_data2<-snp_cnt_data[rowSums(snp_cnt_data)!=0,];
    snp_cnt_data2<-as.data.frame(snp_cnt_data2)
    rownames(snp_cnt_data2)<-rownames(snp_cnt_data)[which(rowSums(snp_cnt_data)!=0)];
    snp_cnt_data3<-snp_cnt_data2[,colSums(snp_cnt_data2) != 0];
    y_label2<-y_label[colSums(snp_cnt_data2)!=0,]
    rownames(y_label2)<-colnames(snp_cnt_data3)
    fwrite(snp_cnt_data3,paste(data,"/",unlist(strsplit(basename(gene),"[.]"))[1],".snp_cnt_merged.txt",sep=""), row.names=T, col.names = T);
    fwrite(y_label2,paste(data,"/",unlist(strsplit(basename(gene),"[.]"))[1],".snp_cnt_merged_label.txt",sep=""), row.names=T, col.names = T);
  }
  rm(snp_cnt_data,snp_cnt_data2,snp_cnt_data3,y_label,y_label2)
  
  snp_aa_files<-list.files(data_path, pattern = "_aa.txt$", full.names = TRUE)
  y_label<-fread(paste(data_path, "/y_label.txt",sep=""))
  if (is.null(gene)){
    snp_aa_data <- Reduce(function(x, y) merge(x, y, all=TRUE), lapply(snp_aa_files,fread)); 
    row.names(snp_aa_data)<-snp_aa_data$V1;snp_aa_data$V1<-NULL
    snp_aa_data[is.na(snp_aa_data)] <- 0
    fwrite(snp_cnt_data,paste(data,"/snp_aa_merged.txt",sep=""), row.names=T, col.names = T);
  } else {
    gene_list<-fread(gene); 
    snp_aa_data<-gene_list;
    x<-NULL
    for (i in 1:length(snp_aa_files)){
    # for (i in 1:3){
      t_data<-as.data.frame(fread(snp_aa_files[i])); rownames(t_data)<-t_data$V1; 
      t_data$V1<-colsplit(string=t_data$V1, pattern=":", names=c("Gene", "AA"))[,1]
      tmp<-t_data[which(t_data$V1 %in% gene_list$HGNC_SYMBOL),];tmp$V1<-NULL
      x<-merge(x,tmp,by="row.names",all=TRUE)
      rownames(x)<-x$Row.names; x$Row.names<-NULL
      x[is.na(x)] <- 0
      print(i)
    }
    rm(t_data,tmp,snp_aa_data)
    snp_aa_data<-x;
    snp_aa_data2<-snp_aa_data[rowSums(snp_aa_data)!=0,]
    snp_aa_data3<-snp_aa_data2[,colSums(snp_aa_data2) != 0];
    y_label2<-y_label[which(colSums(snp_aa_data2)!=0),]
    rownames(y_label2)<-colnames(snp_aa_data3)
    fwrite(snp_aa_data3,paste(data,"/",unlist(strsplit(basename(gene),"[.]"))[1],".snp_aa_merged.txt",sep=""), row.names=T, col.names = T);
    fwrite(y_label2,paste(data,"/",unlist(strsplit(basename(gene),"[.]"))[1],".snp_aa_merged_label.txt",sep=""), row.names=T, col.names = T);
  }
  rm(snp_aa_data,snp_aa_data2,snp_aa_data3,y_label,y_label2,x)
}

#------------------------------
#   build a pathway matrix (c2)
#------------------------------
build_pathway_matrix <- function(adj_exprs, p_file, y_label){
  library(qusage)
  library (plyr)
  pathway<-read.gmt(p_file)
  # pathways<-do.call(c, list(c2_pathway, c5_pathway))

  df2 <- ldply (pathway, data.frame)
  colnames(df2)<-c("pathway","gene")
  reshape_pathway <- dcast(df2, pathway ~ gene, function(x) 1, fill = 0)
  rownames(reshape_pathway)<-reshape_pathway$pathway; reshape_pathway$pathway<-NULL
  
  com_pathway_gene<-intersect(rownames(adj_exprs),colnames(reshape_pathway))
  genes_pathway<-as.matrix(t(adj_exprs[com_pathway_gene,]))
  reshape_pathway2<-as.matrix(t(reshape_pathway[,com_pathway_gene]))
  
  pathway_mutation<-genes_pathway %*% reshape_pathway2; dim(pathway_mutation)
  
  pathway_mutation2<-pathway_mutation[rowSums(pathway_mutation)!=0,]; dim(pathway_mutation2)
  y_adjusted<-y_adjusted[rowSums(pathway_mutation)!=0,]
  pathway_mutation2<-as.data.frame(pathway_mutation2[,colSums(pathway_mutation2)!=0]); dim(pathway_mutation2)
  
  
  fwrite(pathway_mutation2, paste(data,"/pathway_mutation_data.txt",sep=""), row.names = TRUE)
  fwrite(y_adjusted, paste(data,"/pathway_mutation_data_label.txt",sep=""), row.names = TRUE)
}
