library(data.table)
library(maftools)
library(dplyr)

# read maf and clinical files
maf_dir_path="/Users/sjung/Documents/GitHub/gdc-snp/mutect/masked_somatic_mutation/maf_files"
data="/Users/sjung/Documents/GitHub/gdc-snp/mutect/data"

dir.create(data, showWarnings = FALSE)

maf_files<-list.files(maf_dir_path, full.names = TRUE)
clinical_data<-"/Users/sjung/Documents/GitHub/gdc-snp/mutect/masked_somatic_mutation/clinical.cart.2018-02-26/clinical.tsv"
clinical<-fread(clinical_data)
colnames(clinical)[2]<-"Tumor_Sample_Barcode"

# create an aggregated data matrix

y_label<-NULL; y_map<-NULL; y_summary<-NULL; tcga_type<-NULL; 
x<-NULL;

for (i in 1:length(maf_files)){
  laml<-read.maf(maf=maf_files[i],isTCGA=TRUE,clinicalData=clinical)
  tmp<-mutCountMatrix(maf=laml, removeNonMutated=TRUE)
  x<-merge(x,tmp, by="row.names",all=TRUE)
  rownames(x)<-x$Row.names; x$Row.names<-NULL
  y_label<-c(y_label,rep(i,ncol(tmp)))
  y_map<-rbind(y_map,paste(i,unlist(strsplit(basename(maf_files[i]),"[.]"))[2],sep="\t"))
  y_summary<-rbind(y_summary,paste(nrow(tmp),ncol(tmp),sep="\t"))
  sprintf("%f is done", i)
}
save(x,y_label,y_map,y_summary,file=paste(data,"/Combined.RData",sep=""))
x[is.na(x)] <- 0


# y_label<-as.data.frame(y_label)

X<-as.data.frame(t(x))
y_label<-as.data.frame(y_label)
write.table(X,paste(data,"/X",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
write.table(y_label,paste(data,"/y",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
write.table(y_map,paste(data,"/y.map",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)
write.table(y_summary,paste(data,"/y.summary",sep=""),row.names = FALSE, col.names=FALSE,quote = FALSE)


#----------------------------------
# data split
#----------------------------------

library(caret)

X<-fread("/Users/sjung/Documents/GitHub/gdc-snp/mutect/data/X")
y<-fread("/Users/sjung/Documents/GitHub/gdc-snp/mutect/data/y")
Z<-cbind(y,X)

set.seed(3456)
trainIndex <- createDataPartition(Z$V1, p = .8, 
                                  list = FALSE, 
                                  times = 1)

snpTrain <- Z[ trainIndex,]
snpTest  <- Z[-trainIndex,]
fwrite(snpTrain,"/Users/sjung/Documents/GitHub/Benchmarks/Data/Pilot1/snp_train.csv",col.names = F)
fwrite(snpTest,"/Users/sjung/Documents/GitHub/Benchmarks/Data/Pilot1/snp_test.csv",col.names = F)

unique(snpTrain[,1])
