## Rscript script/main_run.R -m input/maf_files -c input/clinical.tsv -d data/gdc_tcga -g 

suppressPackageStartupMessages(library(optparse))

source("/Users/sjung/Documents/GitHub/gdc_snp/script/library.R")

# file and directory path
option_list = list(
  make_option(c("-m", "--maf_dir"),  type="character", help="maf file directory path"),
  make_option(c("-c", "--clinical"),  type="character", help="clinical data"),
  make_option(c("-d", "--data"),  type="character", help="data directory path"),
  make_option(c("-g", "--gene"),  type="character", help="gene list")
  make_option(c("-p", "--pathway"),  type="character", help="pathway file")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(OptionParser(option_list=option_list))

if ( is.null(opt$m) & is.null(opt$c) & is.null(opt$d)) {
    print_help(opt_parser)
    stop("Missing input files path!\n")
} else {
    maf_dir_path<-opt$m;
    clinical_data<-opt$c;
    data<-opt$d;
    print(data)
}

# maf_dir_path="/Users/sjung/Documents/GitHub/gdc_snp/input/maf_files"
# clinical_data<-"/Users/sjung/Documents/GitHub/gdc_snp/input/clinical.tsv"
# data="/Users/sjung/Documents/GitHub/gdc_snp/data/gdc_tcga"
# gene<-"/Users/sjung/Documents/GitHub/gdc_snp/data/cosmic_ostlund_cancer_gene.csv"

## data preprocessing
source("/Users/sjung/Documents/GitHub/gdc_snp/script/maf_data_process.R")

## main
maf_to_matrix(maf_dir_path, clinical_data, data)
if ( is.null(opt$g)){
  combine_matrix(data)
} else {
  combine_matrix(data,opt$g) #combine only the provided genes
}

# normalize snp count by gene length

input_data_matrix<-paste(data,"/",unlist(strsplit(basename(opt$g),"[.]"))[1],".snp_cnt_merged.txt",sep="")
input_data_label<-paste(data,"/",unlist(strsplit(basename(opt$g),"[.]"))[1],".snp_cnt_merged_label.txt",sep="")
# input_data_matrix<-fread("/Users/sjung/Documents/GitHub/gdc_snp/data/gdc_tcga/cosmic_ostlund_cancer_gene.snp_cnt_merged.txt")
# input_data_label<-fread("/Users/sjung/Documents/GitHub/gdc_snp/data/gdc_tcga/cosmic_ostlund_cancer_gene.snp_cnt_merged_label.txt")
rownames(input_data_label)<-input_data_label$V1; input_data_label$V1<-NULL
rownames(input_data_matrix)<-input_data_matrix$V1;input_data_matrix$V1<-NULL
norm_gene_count<-gene_devide_exon_size(input_data_matrix)

norm_gene_count_label<-input_data_label[colSums(norm_gene_count)!=0,]
norm_gene_count<-norm_gene_count[,colSums(norm_gene_count)!=0]
rownames(norm_gene_count_label)<-colnames(norm_gene_count)


fwrite(norm_gene_count,paste(data,"/",unlist(strsplit(basename(opt$g),"[.]"))[1],".snp_cnt_merged_normalized.txt",sep=""), row.names = TRUE)
fwrite(norm_gene_count_label,paste(data,"/",unlist(strsplit(basename(opt$g),"[.]"))[1],".snp_cnt_merged_normalized_label.txt",sep=""), row.names = TRUE)
# fwrite(norm_gene_count,"/Users/sjung/Documents/GitHub/gdc_snp/data/gdc_tcga/cosmic_ostlund_cancer_gene.snp_cnt_merged_normalized.txt", row.names = TRUE)
# fwrite(norm_gene_count_label,"/Users/sjung/Documents/GitHub/gdc_snp/data/gdc_tcga/cosmic_ostlund_cancer_gene.snp_cnt_merged_normalized_label.txt", row.names = TRUE)

# calcuate pathway
build_pathway_matrix(norm_gene_count, opt$p, norm_gene_count_label)
