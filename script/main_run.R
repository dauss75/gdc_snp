source("/Users/sjung/Documents/GitHub/gdc_snp/script/library.R")

# file and directory path
raw_dir_path="/Users/sjung/Documents/GitHub/gdc_snp/input/maf_files"
clinical_data<-"/Users/sjung/Documents/GitHub/gdc_snp/input/clinical.tsv"
data="/Users/sjung/Documents/GitHub/gdc_snp/data"
figs="/Users/sjung/Documents/GitHub/gdc_snp/figs"

## data preprocessing
source("/Users/sjung/Documents/GitHub/gdc_snp/script/maf_data_process.R")

## main
maf_to_matrix(raw_dir_path, clinical_data, data)
# combine_matrix(data)
