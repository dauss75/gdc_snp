maf_to_matrix<-function(maf_dir_path, clinical_data, data){
  # data_preprocessing(raw_dir_path, clinical_data, data)
  dir.create(data, showWarnings = FALSE)
  maf_files<-list.files(maf_dir_path, full.names = TRUE)
  
  clinical<-fread(clinical_data)
  colnames(clinical)[2]<-"Tumor_Sample_Barcode"
  
  # create an aggregated data matrix
  y_cnt<-NULL; y_aa_mis<-NULL; y_aa_all<-NULL; 
  x_cnt<-NULL; x_aa_mis<-NULL; x_aa_all<-NULL;
  y_map<-NULL; y_summary<-NULL; tcga_type<-NULL; 
  
  for (i in 1:length(maf_files)){
    # sink(file=paste(data,"/",basename(maf_files[i]),"_summary.txt",sep=""))
    laml<-read.maf(maf=maf_files[i],isTCGA = FALSE, clinicalData=clinical,removeDuplicatedVariants=F)
    # sink()
    
    mut_cnt<-mutCountMatrix(maf=laml, includeSyn=F)
    # x_cnt<-merge(x_cnt,mut_cnt, by="row.names",all=TRUE)
    # rownames(x_cnt)<-x_cnt$Row.names; x_cnt$Row.names<-NULL
    y_cnt<-c(y_cnt,rep(i-1,ncol(mut_cnt)))
    # y_map<-rbind(y_map,paste(i,unlist(strsplit(basename(maf_files[i]),"[.]"))[2],sep="\t"))
    # y_summary<-rbind(y_summary,paste(nrow(mut_cnt),ncol(mut_cnt),sep="\t"))
    
    # png(file = paste("/Users/sjung/Documents/GitHub/gdc_snp/figs/",unlist(strsplit(basename(maf_files[i]),"[.]"))[2],".png",sep=""), width = 1024, height = 768)
    # plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, titleSize = c(14,12), fs=14, statFontSize=5)
    # dev.off()
    
    # laml@data[["Hugo_Symbol"]]<-paste(laml@data[["Hugo_Symbol"]],":",laml@data[["Protein_Change"]],sep="")
    laml@data[["Hugo_Symbol"]]<-paste(laml@data[["Hugo_Symbol"]],":",laml@data[["HGVSp_Short"]],sep="")
    #snp_aa_all<-mutCountMatrix(laml, includeSyn = FALSE)
    #x_aa_all<-merge(x_aa_all,snp_aa_all, by="row.names",all=TRUE)
    #rownames(x_aa_all)<-x_aa_all$Row.names; x_aa_all$Row.names<-NULL
    #y_aa_all<-c(y_aa_all,rep(i-1,ncol(snp_aa_all)))
    
    snp_aa_mis<-mutCountMatrix(laml, includeSyn = FALSE, countOnly = "Missense_Mutation")
    # x_aa_mis<-merge(x_aa_mis,snp_aa_mis, by="row.names",all=TRUE)
    # rownames(x_aa_mis)<-x_aa_mis$Row.names; x_aa_mis$Row.names<-NULL
    y_aa_mis<-c(y_aa_mis,rep(i-1,ncol(snp_aa_mis)))
    
    write.table(mut_cnt,paste(data,"/",basename(maf_files[i]),"_mut_cnt",sep=""),row.names = TRUE, col.names=TRUE,quote = FALSE)
    fwrite(mut_cnt,paste(data,"/",basename(maf_files[i]),"mut_cnt",sep=""),row.names = TRUE, col.names=TRUE)
    write.table(snp_aa_mis,paste(data,"/",basename(maf_files[i]),".snp_aa",sep=""),row.names = TRUE, col.names=TRUE,quote = FALSE)
    fwrite(snp_aa_mis,paste(data,"/",basename(maf_files[i]),"_snp_aa",sep=""),row.names = TRUE, col.names=TRUE)
    
    # rm(laml)
    
    # sprintf("%f is done", i)
  }
  # rm(snp_aa_mis, mut_cnt,laml,clinical)
  # write.table(rownames(x),paste(data,"/mut_gene_list.txt",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
  # write.table(colnames(x),paste(data,"/mut_sample_list.txt",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
  # save(x,y_label,y_map,y_summary,file=paste(data,"/combined_mut_data.RData",sep=""))
  # x_cnt[is.na(x_cnt)] <- 0; x_aa_mis[is.na(x_aa_mis)] <- 0; #x_aa_all[is.na(x_aa_all)] <- 0; 
  
  # x_cnt<-as.data.frame(t(x_cnt)); x_aa_mis<-as.data.frame(t(x_aa_mis)); #x_aa_all<-as.data.frame(t(x_aa_all))
  y_cnt<-as.data.frame(y_cnt); y_aa_mis<-as.data.frame(y_aa_mis); #y_aa_all<-as.data.frame(y_aa_all)
  # write.table(x_cnt,paste(data,"/x_cnt",sep=""),row.names = TRUE, col.names=TRUE,quote = FALSE)
  # fwrite(x_cnt,paste(data,"/x_cnt_2",sep=""),sep=" ", row.names=T, col.names = T)
  # write.table(x_aa_mis,paste(data,"/x_aa_mis",sep=""),row.names = TRUE, col.names=TRUE,quote = FALSE)
  # fwrite(x_aa_mis,paste(data,"/x_aa_mis_2",sep=""),sep=" ", row.names=T, col.names = T)
  #write.table(x_aa_all,paste(data,"/x_aa_all",sep=""),row.names = TRUE, col.names=TRUE,quote = FALSE)
  
  write.table(y_cnt,paste(data,"/y_cnt_mis",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
  write.table(y_aa_mis,paste(data,"/y_aa_mis",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
  #write.table(y_aa_all,paste(data,"/y_aa_all",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
  # write.table(y_map,paste(data,"/y_mut_count.map",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
  # write.table(y_summary,paste(data,"/y_mut_count.summary",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
} 


gene_devide_exon_size<-function(exprs){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
  colnames(exonic.gene.sizes)[1]<-"size"
  exonic.gene.sizes$genesymbol<-getSYMBOL(na.omit(rownames(exonic.gene.sizes)), data='org.Hs.eg')
  exonic.gene.sizes<-na.omit(exonic.gene.sizes)
  rownames(exonic.gene.sizes)<-exonic.gene.sizes$genesymbol
  
  y<-as.data.frame(exprs[,1]); x<-as.data.frame(t(exprs[,-1]))
  
  # retrieve pathway
  com_gene<-intersect(rownames(x),rownames(exonic.gene.sizes))
  x2<-x[com_gene,]; x2<-na.omit(x2)
  exonic.gene.sizes2<-exonic.gene.sizes[com_gene,]
  exprs_devided_by_exon_size<-x2/exonic.gene.sizes2$size
  exprs_devided_by_exon_size
}

combine_matrix<-function(data_path){
#  cnt_files<-list.files(data_path, pattern = "_cnt$", full.names = TRUE)
#  cnt_data <- Reduce(function(x, y) merge(x, y, all=TRUE), lapply(cnt_files,fread)); row.names(cnt_data)<-cnt_data$V1;cnt_data$V1<-NULL 
#  cnt_data[is.na(cnt_data)] <- 0
#  fwrite(cnt_data,paste(data,"/cnt_merged_data",sep=""),sep=" ", row.names=T, col.names = T);
#  rm(cnt_data)
  aa_files<-list.files(data_path, pattern = "_aa$", full.names = TRUE)
  aa_data2 <- Reduce(function(x, y) merge(x, y, all=TRUE), lapply(aa_files,fread)); 
  row.names(aa_data)<-aa_data$V1;
  aa_data$V1<-NULL 
  aa_data[is.na(aa_data)] <- 0; 
  fwrite(aa_data,paste(data,"/aa_merged_data",sep=""),sep=" ", row.names=T, col.names = T)
}