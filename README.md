# GDC SNP data analysis using deep learning

## Table of Contents

- [Data Preparation](#data-preparation)
   - [Feature generation](#feature-generation)
   - [Feature evaluation](#feature-evaluation)
      - [PCA and Covariates](#pca-and-covariates)
- [Simulation](#simulation)
     - [Trial 1 with Mutect Data](#trial-1-with-mutect-data)
- [Supplementary Information](#supplementary-information)
     - [List of mutations](#list-of-mutations)
     - [Somatic Variant Calling Workflow](#somatic-variant-calling-workflow)

## Data Preparation

Masked somatic mutation data for 33 cancer types (WGS) are available which are processed by 4 somatic calling tools: [MuSE (MS)](http://www.biorxiv.org/content/early/2016/05/25/055467.abstract), [MuTect2 (MK)](https://www.nature.com/articles/nbt.2514), [VarScan2 (VS)](https://genome.cshlp.org/content/22/3/568.short), and [SomaticSniper (SS)](http://bioinformatics.oxfordjournals.org/content/28/3/311.short). The data is from whole genome sequencing.

### Feature generation

Generate a count matrix of mutation that is based on the number of mutation per gene per sample.

- First, we include all mutations - [sysnonymous](https://en.wikipedia.org/wiki/Synonymous_substitution) and [non-sysnonymous](https://en.wikipedia.org/wiki/Nonsynonymous_substitution). Note that sysnonymous mutation does not change encoded amino acid where as non-sysnonymous mutation does.
- Second, we include only non-sysnonymous mutation.


TBD

### Feature evaluation

TBD

#### PCA and Covariates

TBD


## Simulation

### Test 1 with Muteck Data

Dataset: Somatic mutation data from MuTeck2 (MT) in total of 8132 samples

| Label | Cancer Type  | # of features (MT) | # of Samples (MT)|
| ------|------------- | -------------------|------------------|
|   1 	|    ACC 	     |          4752	    |        92        |
|   2	  |    BLCA 	   |          16551	    |       412        |
|   3		|    BRCA      |          16131	    |       985        |
|   4		|    CESC 	   |          15342	    |       289        |
|   5		|    CHOL      |          2333 	    |        50        |
|   6		|    COAD 	   |          18037	    |       399        |
|   7		|    DLBC 	   |          2976 	    |        37        |
|   8		|    ESCA 	   |          10719	    |       184        |
|   9		|    GBM 	     |          14778	    |       392        |
|  10		|    HNSC 	   |          15436	    |       508        |
|  11		|    KICH 	   |          1634 	    |        66        |
|  12		|    KIRC 	   |          9444 	    |       336        |
|  13		|    KIRP 	   |          8689 	    |       281        |
|  14		|    LAML 	   |          4810 	    |       140        |
|  15		|    LGG 	     |          10728 	  |       507        |
|  16		|    LIHC 	   |          12704 	  |       363        |
|  17		|    LUAD 	   |          16980 	  |       565        |
|  18		|    LUSC 	   |          16784 	  |       491        |
|  19		|    MESO 	   |          2027 	    |        81        |
|  20		|    OV 	     |          14226 	  |       436        |
|  21		|    PAAD 	   |          9967  	  |       175        |
|  22		|    PCPG 	   |          1416  	  |       178        |
|  23		|    PRAD 	   |          9829  	  |       494        |
|  24		|    READ 	   |          13333 	  |       137        |
|  25		|    SARC 	   |          8119  	  |       237        |
|  26		|    SKCM 	   |          17509 	  |       467        |
|  27		|    STAD 	   |          17431	    |       437        |
|  28		|    TGCT 	   |          1749 	    |       144        |
|  29		|    THCA 	   |          4643 	    |       490        |
|  30		|    THYM 	   |          2084 	    |       123        |
|  31		|    UCEC 	   |          18998	    |       530        |
|  32		|    UCS 	     |          5110  	  |        57        |
|  33		|    UVM   	   |          1082  	  |        80        |

**Summary of DNN**

- Input data: training (8132 samples x 19384 features) and test (2031 samples x 19384 features)
- class labels: 33
- primarily results

![Screenshot](figs/trial1_dropout_0.1.png)  
*Figure 1. accuracy (A) and loss (B) scores are shown over epochs with dropout=0.1. The training accuracy continues to improve whereas the test accuracy
stabilizes around 20 epochs that is explanined by the loss score where the validation loss does not decrease at around 20 epochs.*

![Screenshot](figs/trial1_dropout_0.2.png)  
*Figure 2. accuracy (A) and loss (B) scores are shown over epochs with dropout=0.2. With the change of dropout to 0.2, slight improvements have been observed, but not the overfitting problem.*

**Evaluation with cancer data with >100 samples**

| Label | Cancer Type  | # of features (MT) | # of Samples (MT)|
| ------|------------- | -------------------|------------------|
|   2	  |    BLCA 	   |          16551	    |       412        |
|   3		|    BRCA      |          16131	    |       985        |
|   4		|    CESC 	   |          15342	    |       289        |
|   6		|    COAD 	   |          18037	    |       399        |
|   8		|    ESCA 	   |          10719	    |       184        |
|   9		|    GBM 	     |          14778	    |       392        |
|  10		|    HNSC 	   |          15436	    |       508        |
|  12		|    KIRC 	   |          9444 	    |       336        |
|  13		|    KIRP 	   |          8689 	    |       281        |
|  14		|    LAML 	   |          4810 	    |       140        |
|  15		|    LGG 	     |          10728 	  |       507        |
|  16		|    LIHC 	   |          12704 	  |       363        |
|  17		|    LUAD 	   |          16980 	  |       565        |
|  18		|    LUSC 	   |          16784 	  |       491        |
|  20		|    OV 	     |          14226 	  |       436        |
|  21		|    PAAD 	   |          9967  	  |       175        |
|  22		|    PCPG 	   |          1416  	  |       178        |
|  23		|    PRAD 	   |          9829  	  |       494        |
|  24		|    READ 	   |          13333 	  |       137        |
|  25		|    SARC 	   |          8119  	  |       237        |
|  26		|    SKCM 	   |          17509 	  |       467        |
|  27		|    STAD 	   |          17431	    |       437        |
|  28		|    TGCT 	   |          1749 	    |       144        |
|  29		|    THCA 	   |          4643 	    |       490        |
|  30		|    THYM 	   |          2084 	    |       123        |
|  31		|    UCEC 	   |          18998	    |       530        |


----------

## Supplementary Information

### List of mutations

The variant classifications from TCGA MAF.

- Intron -- variant lies between exons within the bounds of the chosen transcript.
- 5'UTR -- variant is on the 5'UTR for the chosen transcript
- 3'UTR -- variant is on the 3'UTR for the chosen transcript
- IGR -- intergenic region. Does not overlap any transcript.
- 5'Flank -- the variant is upstream of the chosen transcript (within 3kb)
- [Missense_Mutation](https://en.wikipedia.org/wiki/Missense_mutation) -- the point mutation alters the protein structure by one amino acid.
- [Nonsense_Mutation](https://en.wikipedia.org/wiki/Nonsense_mutation) -- a premature stop codon is created by the variant.
- Nonstop_Mutation -- variant removes stop codon.
- Silent -- variant is in coding region of the chosen transcript, but protein structure is identical. I.e. a synonymous mutation
- Splice_Site -- the variant is within two bases of a splice site. See the secondary classification to determine if it lies on the exon or intron side.
- In_Frame_Del -- deletion that keeps the sequence in frame.
- In_Frame_Ins -- insertion that keeps the sequence in frame.
- [Frame_Shift_Ins](https://en.wikipedia.org/wiki/Frameshift_mutation) -- insertion that moves the coding sequence out of frame.
- [Frame_Shift_Del](https://en.wikipedia.org/wiki/Frameshift_mutation) -- deletion that moves the coding sequence out of frame.
- Start_Codon_SNP -- point mutation that overlaps the start codon.
- Start_Codon_Ins -- insertion that overlaps the start codon.
- Start_Codon_Del -- deletion that overlaps the start codon.
- De_novo_Start_InFrame -- New start codon is created by the given variant using the chosen transcript. However, it is in frame relative to the coded protein.
- De_novo_Start_OutOfFrame -- New start codon is created by the given variant using the chosen transcript. However, it is out of frame relative to the coded protein.
- RNA -- variant lies on one of the RNA transcripts.
- lincRNA -- variant lies on one of the lincRNAs.

### Somatic Variant Calling Workflow

Aligned tumor-normal BAM pairs are processed through the Somatic Mutation Calling Workflow.

![Screenshot](figs/gdc-alignment.png)

*DNA-seq alignment pipeline-Courtesy of GDC*

![Screenshot](figs/gdc-Broad_MuTect.png)

*Mutect somatic variant calling pipeline-Courtesy of GDC*
