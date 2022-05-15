#Authors: David J. Bradshaw
#Emails: dbradshaw2015@fau.edu, dbradshaw3366@gmail.com

#Importing information into R was done using the following tximport vignette
#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

#Analysis using DESeq2 was done using the following vignette
#http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

#Setting the working directory
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/R")


#Install packages#
BiocManager::install("tximportData")
BiocManager::install("BUSpaRse")
BiocManager::install("tximport")
BiocManager::install("edgeR")
BiocManager::install("apeglm")
install.packages("tidyverse")
install.packages("readxl")
BiocManager::install("topGO")
BiocManager::install("org.Dr.eg.db")
BiocManager::install("org.Mn.eg.db")
BiocManager::install("org.Hs.eg.db")
install.packages("dendextend")

BiocManager::install("RnaSeqGeneEdgeRQL")


#Load libraries#
library(RnaSeqGeneEdgeRQL)
library(topGO)
library(phyloseq)
library(DESeq2)
library(tximportData)
library(tximport)
library(edgeR)
library(apeglm)
library(tidyverse)
library(readxl)
library(phyloseq)
library(tximeta)
library(gtools)
library(pheatmap)
library(dendextend)
library(FSA)
library(biomaRt)
library(GO.db)
library(clusterProfiler)
library(rcompanion)
library(ggpubr)
library(pathview)
library(ggvenn)

###PREPARE FOR IMPORTING OF CLC GENOMICS DATA AND PHYLOSEQ###

##Prepare your working environment

#In your working directory you should have:
#CLC Genomics Excel file - each sheet should be a different sample
#sample.txt - tab deliminated file that has your metadata for your samples including a column named "sample" that lists the names that you want your samples to have in the ultimate abundance table
#tx2gene.txt - tab deliminated file with the first column being transcript names (TXNAME) and the second column being gene names (GENEID), these are based upon the CLC Genomics columns of "Name" and "Gene name" or if you want to get technical "Transcript ID" and "GeneID" 

#Set your working directory
getwd()
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis")
getwd()

#Upload your metadata file
samples <- read.csv("salinity_transcriptomics_metadata.csv", row.names = 1)

#Create a new folder  in that work directory
dir.create("CLC_Genomics_Files")

#Set your working directory to that new folder for this next section
getwd()
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/CLC_Genomics_Files")
getwd()

##Prep abundance files##

#Load function that will write every sheet in an Excel file (designated path here) to a series of tsv files#
read_then_tsv <- function(sheet, path) {
  pathbase <- path %>%
    basename() %>%
    tools::file_path_sans_ext()
  path %>%
    read_excel(sheet = sheet) %>%
    write_tsv(paste0(pathbase, "-", sheet, ".tsv"))
}

#Set path to equal the Excel with all your sample sheets
path <- "C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/Salinity_TE_all.xlsx"

#Make a tsv file for every sheet in the Excel file
path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_then_tsv, path = path)

##Rename all the files using your metadata table

#list.files() - just lists files in your present working directory (pwd) similar to ls in Linux
myFiles <- list.files()

#See what this list looks like
myFiles

#mixed sort allows numbers to be sorted in correct order instead of treating them like a factor, rev reverses the order
myFiles <- rev(mixedsort(myFiles))

#Check that it worked
myFiles

#use my files to change the name of all files in the pwd
file.rename(myFiles, paste0(samples$sample,".tsv"))

#Check to see that it worked
list.files()

#Set your working directory back to the other folder
getwd()
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis")
getwd()

##Make a vector pointing to quantification files

#Make a value that points to folder with your new set of tsv files for each sample 
dir <- "C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/CLC_Genomics_Files"
list.files(dir)

#Create a dataframe of the file names to use later
names <- as.data.frame(list.files(dir))
names$names <- names$`list.files(dir)`
names

#Create a list of files and their associated paths
files <- file.path(dir, names$names)
files
all(file.exists(files))

#Use the metadata table to add labels to each of these files
names(files) <- paste0(samples$sample)
names(files)

#Load in your transcripts to genes file
tx2gene <- read.csv("tx2gene.csv")
head(tx2gene)

#Use tximport to import all the files into a list of data frames that have transcript information

#You are using the type = "none" because tximport does not have a built in way to deal with CLC Genomics outputs. That means you must define the following information by telling the tool the names of the columns in the tables that have this information: column with transcript name that must match whats in your tx2gene table (txIdCol), column with TPM or FPKM abundances (abundanceCol), column with estimate counts (countsCol), and column with length information (lengthCol)

#You must also define the importer used to get your tables into R, use read_tsv here since it is opposite of the write_tsv used earlier

#countsFromAbundance - allows you to get counts scaled to library size (scaledTPM), or scaled using average transcript length over samples and then by library size(lengthScaledTPM), or scaled using the median transcript length among isoforms fo the gene and then by library size (dtuScaledTPM) that is generally used with txOut=TRUE and if you want to do differential transcript usage analysis

#txOut = TRUE - means that you are just generating transcript level outputs

#Make a set of tables with transcript information
txi.CLC.tx <- tximport(files, type="none", txIdCol = "Name", tx2gene=tx2gene, abundanceCol = "TPM", countsCol = "Total transcript reads", lengthCol = "Transcript length", importer = read_tsv, txOut=TRUE)

names(txi.CLC.tx)
#"abundance" "counts" "length" "countsFromAbundance"

#Extract tables of each of the elements of this tximport object and save them
head(txi.CLC.tx$counts)
tx_counts_table <- as.data.frame(txi.CLC.tx$counts)
write.csv(tx_counts_table, "salinity_tx_counts.csv")

head(txi.CLC.tx$abundance)
tx_abundance_table <- as.data.frame(txi.CLC.tx$abundance)
write.csv(tx_abundance_table, "salinity_tx_abundance.csv")

head(txi.CLC.tx$length)
tx_length_table <- as.data.frame(txi.CLC.tx$length)
write.csv(tx_length_table, "salinity_tx_length.csv")

#Make a gene-level object
txi.CLC.gn <- summarizeToGene(txi.CLC.tx, tx2gene)

names(txi.CLC.gn)
#"abundance" "counts" "length" "countsFromAbundance"

#Make a set of tables with gene information
head(txi.CLC.gn$counts)
gn_counts_table <- as.data.frame(txi.CLC.gn$counts)
write.csv(gn_counts_table, "salinity_gn_counts.csv")

head(txi.CLC.gn$abundance)
gn_abundance_table <- as.data.frame(txi.CLC.gn$abundance)
write.csv(gn_abundance_table, "salinity_gn_abundance.csv")

#Take just one column of the length table to be used later and relate it length
head(txi.CLC.gn$length)
gn_length_table <- as.data.frame(txi.CLC.gn$length)
gn_length_table <- gn_length_table[,1, drop=FALSE]
colnames(gn_length_table)=c("length")
head(gn_length_table)
write.csv(gn_length_table, "salinity_gn_length.csv")

##Get sequence information for transcripts and genes
tx_Totals <- as.data.frame(colSums(tx_counts_table))
colnames(tx_Totals) <- "Sum"
max(tx_Totals$Sum)/min(tx_Totals$Sum) #3.57964
sum(tx_Totals$Sum) #553,366,960
dim(tx_counts_table) #34893    57

gene_Totals <- as.data.frame(colSums(gn_counts_table))
colnames(gene_Totals) <- "Sum"
max(gene_Totals$Sum)/min(gene_Totals$Sum) #3.579674
sum(gene_Totals$Sum) #553366960
dim(gn_counts_table) #23800    57

write.csv(tx_Totals, "transcrip_totals.csv")

##Prepare transcript-level data for phyloseq
#Convert counts and metadata into phyloseq friendly versions#
META <- import_qiime_sample_data("salinity_transcriptomics_metadata.txt")

#Make a table of the transcript counts into phyloseq format
TXCTS = otu_table(tx_counts_table, taxa_are_rows = TRUE)

#Create a phyloseq object
TX_dataA <- merge_phyloseq(TXCTS, META)

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
nofactor_tx_dds <- DESeqDataSetFromTximport(txi.CLC.tx, samples, ~1)

#Conduct DESEQ2 test#
nofactor_tx_dds = DESeq(nofactor_dds, fitType="local")

#Make a copy of Fdata so you can have a designated DESeq transformed copy
TX_dataA
VST_TX_dataA = TX_dataA
VST_TX_dataA

#Switch the asv table with the DESeq2 transformed data
otu_table(VST_TX_dataA) <- otu_table(getVarianceStabilizedData(nofactor_tx_dds), taxa_are_rows = TRUE)

#Extract out VST transcript table
VST_TX_Table <- as.data.frame(otu_table(VST_TX_dataA))

#Print the VST transcript table
write.csv(otu_table(VST_TX_dataA), "VST_TX_Table.csv")

#Make a PCoA
TX_PCoA_Ord <- ordinate(VST_TX_dataA,"PCoA")
TX_PCoA_Plot = plot_ordination(VST_TX_dataA, TX_PCoA_Ord, color="Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"), labels=c("3", "6", "12", "15", "18", "20", "24"))+
  #ggtitle("Tissue PCoA by DPH and Salinity") +
  scale_shape_discrete(name="Salinity (ppt)", labels=c("10", "20", "30"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13)) 
TX_PCoA_Plot


##Prepare gene-level data for phyloseq
GNCTS = otu_table(gn_counts_table, taxa_are_rows = TRUE)
GN_dataA <- merge_phyloseq(GNCTS, META)

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
nofactor_gn_dds <- DESeqDataSetFromTximport(txi.CLC.gn, samples, ~1)

#Conduct DESEQ2 test#
nofactor_gn_dds = DESeq(nofactor_gn_dds, fitType="local")

#Make a copy of Fdata so you can have a designated DESeq transformed copy
GN_dataA
VST_GN_dataA = GN_dataA
VST_GN_dataA

#Switch the asv table with the DESeq2 transformed data
otu_table(VST_GN_dataA) <- otu_table(getVarianceStabilizedData(nofactor_gn_dds), taxa_are_rows = TRUE)

#Extract out VST gene table
VST_GN_Table <- as.data.frame(otu_table(VST_GN_dataA))

#Print the VST gene table for PRIMER7/PERMANOVA+ analysis (TABLE T4)
write.csv(VST_GN_Table, "VST_GN_Table.csv")

#Extract out the VST abundance data from the phyloseq object
data <- as.data.frame(otu_table(VST_GN_dataA))

#Check how different the sample sums are between samples
max(rowSums(t(data)))/min(rowSums(t(data))) #1.046653

#Visualize the samples sums
hist(rowSums(t(data)))

#Make a PCoA

##SUPP FIGURE 1

GN_PCoA_Ord <- ordinate(VST_GN_dataA,"PCoA")
GN_PCoA_Plot = plot_ordination(VST_GN_dataA, GN_PCoA_Ord, color="Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"), labels=c("3", "6", "12", "15", "18", "20", "24"))+
  #ggtitle("Tissue PCoA by DPH and Salinity") +
  scale_shape_discrete(name="Salinity (ppt)", labels=c("10", "20", "30"))+
  theme_bw() +
  geom_point(size=4)+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13)) 
GN_PCoA_Plot

#Publish your gene PCoA
tiff('Supp Fig 1 Gene-level PCoA.tiff', units="in", width=8, height=8, res=300)
GN_PCoA_Plot
dev.off()

###SUMMARIZING THE POOLED DEGS

##TABLE 2

#Upload all the DEG lists
All_S10_S30 <- read.csv("All DPH 10 vs. 30 ppt.csv")
All_S10_S20 <- read.csv("All DPH 10 vs. 20 ppt.csv")
All_S20_S30 <- read.csv("All DPH 20 vs. 30 ppt.csv")

##Subset DEGs and determine their number

DEGS_All_S10_S20 <- filter(All_S10_S20, FDR.p.value < 0.05 & FDR.p.value != "#N/A")#5
DEGS_All_S10_S30 <- filter(All_S10_S30, FDR.p.value < 0.05 & FDR.p.value != "#N/A")#12
DEGS_All_S20_S30 <- filter(All_S20_S30, FDR.p.value < 0.05 & FDR.p.value != "#N/A")#0

#Copy results to Excel
clipr::write_clip(DEGS_All_S10_S20)
clipr::write_clip(DEGS_All_S10_S30)
clipr::write_clip(DEGS_All_S20_S30)

##Pooled S10 vs S20
#DEGs associated with S10
as.numeric(length( which( All_S10_S20$FDR.p.value < 0.05 & All_S10_S20$FDR.p.value!= "#N/A" & All_S10_S20$Fold.change >0) ))#3

#DEGs associated with S20
as.numeric(length( which( All_S10_S20$FDR.p.value < 0.05 & All_S10_S20$FDR.p.value!= "#N/A" & All_S10_S20$Fold.change <0) )) #2

##Pooled S10 vs S30
#DEGs associated with S10
as.numeric(length( which( All_S10_S30$FDR.p.value < 0.05 & All_S10_S30$FDR.p.value!= "#N/A" & All_S10_S30$Fold.change >0) ))#2

#DEGs associated with S30
as.numeric(length( which( All_S10_S30$FDR.p.value < 0.05 & All_S10_S30$FDR.p.value!= "#N/A" & All_S10_S30$Fold.change <0) )) #10

##Pooled S20 vs S30
#DEGs associated with S20
as.numeric(length( which( All_S20_S30$FDR.p.value < 0.05 & All_S20_S30$FDR.p.value!= "#N/A" & All_S20_S30$Fold.change >0) ))#0

#DEGs associated with S30
as.numeric(length( which( All_S20_S30$FDR.p.value < 0.05 & All_S20_S30$FDR.p.value!= "#N/A" & All_S20_S30$Fold.change <0) )) #0


###SUMMARIZING THE DPH DEGS

##TABLE 3

#Set your working directory to the location of your DEG lists that have been saved as csvs from the original Excel sheet
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/Salinity DEG Lists")

#Upload all the DEG lists
D3_S10_S30 <- read.csv("3 DPH 10 vs. 30 ppt.csv")
D3_S10_S20 <- read.csv("3 DPH 10 vs. 20 ppt.csv")
D3_S20_S30 <- read.csv("3 DPH 20 vs. 30 ppt.csv")
D6_S10_S30 <- read.csv("6 DPH 10 vs. 30 ppt.csv")
D6_S10_S20 <- read.csv("6 DPH 10 vs. 20 ppt.csv")
D6_S20_S30 <- read.csv("6 DPH 20 vs. 30 ppt.csv")
D12_S10_S30 <- read.csv("12 DPH 10 vs. 30 ppt.csv")
D12_S10_S20 <- read.csv("12 DPH 10 vs. 20 ppt.csv")
D12_S20_S30 <- read.csv("12 DPH 20 vs. 30 ppt.csv")
D15_S10_S30 <- read.csv("15 DPH 10 vs. 30 ppt.csv")
D15_S10_S20 <- read.csv("15 DPH 10 vs. 20 ppt.csv")
D15_S20_S30 <- read.csv("15 DPH 20 vs. 30 ppt.csv")
D18_S10_S30 <- read.csv("18 DPH 10 vs. 30 ppt.csv")
D18_S10_S20 <- read.csv("18 DPH 10 vs. 20 ppt.csv")
D18_S20_S30 <- read.csv("18 DPH 20 vs. 30 ppt.csv")
D20_S10_S30 <- read.csv("20 DPH 10 vs. 30 ppt.csv")
D20_S10_S20 <- read.csv("20 DPH 10 vs. 20 ppt.csv")
D20_S20_S30 <- read.csv("20 DPH 20 vs. 30 ppt.csv")
D24_S10_S30 <- read.csv("24 DPH 10 vs. 30 ppt.csv")
D24_S10_S20 <- read.csv("24 DPH 10 vs. 20 ppt.csv")
D24_S20_S30 <- read.csv("24 DPH 20 vs. 30 ppt.csv")

#convert working directory back to the base one
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis")

#Make subsets of DEGS that are up- (positive) and down- (negative) regulated
Up_D3_S10_S30 <- subset(D3_S10_S30, Log..fold.change >0)
Down_D3_S10_S30 <- subset(D3_S10_S30, Log..fold.change <0)
Up_D3_S10_S20 <- subset(D3_S10_S20, Log..fold.change >0)
Down_D3_S10_S20 <- subset(D3_S10_S20, Log..fold.change <0)
Up_D3_S20_S30 <- subset(D3_S20_S30, Log..fold.change >0)
Down_D3_S20_S30 <- subset(D3_S20_S30, Log..fold.change <0)

Up_D6_S10_S30 <- subset(D6_S10_S30, Log..fold.change >0)
Down_D6_S10_S30 <- subset(D6_S10_S30, Log..fold.change <0)
Up_D6_S10_S20 <- subset(D6_S10_S20, Log..fold.change >0)
Down_D6_S10_S20 <- subset(D6_S10_S20, Log..fold.change <0)
Up_D6_S20_S30 <- subset(D6_S20_S30, Log..fold.change >0)
Down_D6_S20_S30 <- subset(D6_S20_S30, Log..fold.change <0)

Up_D12_S10_S30 <- subset(D12_S10_S30, Log..fold.change >0)
Down_D12_S10_S30 <- subset(D12_S10_S30, Log..fold.change <0)
Up_D12_S10_S20 <- subset(D12_S10_S20, Log..fold.change >0)
Down_D12_S10_S20 <- subset(D12_S10_S20, Log..fold.change <0)
Up_D12_S20_S30 <- subset(D12_S20_S30, Log..fold.change >0)
Down_D12_S20_S30 <- subset(D12_S20_S30, Log..fold.change <0)

Up_D15_S10_S30 <- subset(D15_S10_S30, Log..fold.change >0)
Down_D15_S10_S30 <- subset(D15_S10_S30, Log..fold.change <0)
Up_D15_S10_S20 <- subset(D15_S10_S20, Log..fold.change >0)
Down_D15_S10_S20 <- subset(D15_S10_S20, Log..fold.change <0)
Up_D15_S20_S30 <- subset(D15_S20_S30, Log..fold.change >0)
Down_D15_S20_S30 <- subset(D15_S20_S30, Log..fold.change <0)

Up_D18_S10_S30 <- subset(D18_S10_S30, Log..fold.change >0)
Down_D18_S10_S30 <- subset(D18_S10_S30, Log..fold.change <0)
Up_D18_S10_S20 <- subset(D18_S10_S20, Log..fold.change >0)
Down_D18_S10_S20 <- subset(D18_S10_S20, Log..fold.change <0)
Up_D18_S20_S30 <- subset(D18_S20_S30, Log..fold.change >0)
Down_D18_S20_S30 <- subset(D18_S20_S30, Log..fold.change <0)

Up_D20_S10_S30 <- subset(D20_S10_S30, Log..fold.change >0)
Down_D20_S10_S30 <- subset(D20_S10_S30, Log..fold.change <0)
Up_D20_S10_S20 <- subset(D20_S10_S20, Log..fold.change >0)
Down_D20_S10_S20 <- subset(D20_S10_S20, Log..fold.change <0)
Up_D20_S20_S30 <- subset(D20_S20_S30, Log..fold.change >0)
Down_D20_S20_S30 <- subset(D20_S20_S30, Log..fold.change <0)

Up_D24_S10_S30 <- subset(D24_S10_S30, Log..fold.change >0)
Down_D24_S10_S30 <- subset(D24_S10_S30, Log..fold.change <0)
Up_D24_S10_S20 <- subset(D24_S10_S20, Log..fold.change >0)
Down_D24_S10_S20 <- subset(D24_S10_S20, Log..fold.change <0)
Up_D24_S20_S30 <- subset(D24_S20_S30, Log..fold.change >0)
Down_D24_S20_S30 <- subset(D24_S20_S30, Log..fold.change <0)


#Create a list of all these dataframes
x <- list(Up_D3_S10_S20, Down_D3_S10_S20, Up_D3_S10_S30, Down_D3_S10_S30, Up_D3_S20_S30, Down_D3_S20_S30,Up_D6_S10_S20, Down_D6_S10_S20, Up_D6_S10_S30, Down_D6_S10_S30, Up_D6_S20_S30, Down_D6_S20_S30,Up_D12_S10_S20, Down_D12_S10_S20, Up_D12_S10_S30, Down_D12_S10_S30, Up_D12_S20_S30, Down_D12_S20_S30,Up_D15_S10_S20, Down_D15_S10_S20, Up_D15_S10_S30, Down_D15_S10_S30, Up_D15_S20_S30, Down_D15_S20_S30,Up_D18_S10_S20, Down_D18_S10_S20, Up_D18_S10_S30, Down_D18_S10_S30, Up_D18_S20_S30, Down_D18_S20_S30, Up_D20_S10_S20, Down_D20_S10_S20, Up_D20_S10_S30, Down_D20_S10_S30, Up_D20_S20_S30, Down_D20_S20_S30,Up_D24_S10_S20, Down_D24_S10_S20, Up_D24_S10_S30, Down_D24_S10_S30, Up_D24_S20_S30, Down_D24_S20_S30)

#Add names to each one
names(x) <- list("Up_D3_S10_S20", "Down_D3_S10_S20", "Up_D3_S10_S30", "Down_D3_S10_S30", "Up_D3_S20_S30", "Down_D3_S20_S30", "Up_D6_S10_S20", "Down_D6_S10_S20", "Up_D6_S10_S30", "Down_D6_S10_S30", "Up_D6_S20_S30", "Down_D6_S20_S30",    "Up_D12_S10_S20", "Down_D12_S10_S20", "Up_D12_S10_S30", "Down_D12_S10_S30", "Up_D12_S20_S30", "Down_D12_S20_S30","Up_D15_S10_S20", "Down_D15_S10_S20", "Up_D15_S10_S30", "Down_D15_S10_S30", "Up_D15_S20_S30", "Down_D15_S20_S30",
                 "Up_D18_S10_S20", "Down_D18_S10_S20", "Up_D18_S10_S30", "Down_D18_S10_S30", "Up_D18_S20_S30", "Down_D18_S20_S30", "Up_D20_S10_S20", "Down_D20_S10_S20", "Up_D20_S10_S30", "Down_D20_S10_S30", "Up_D20_S20_S30", "Down_D20_S20_S30",
                 "Up_D24_S10_S20", "Down_D24_S10_S20", "Up_D24_S10_S30", "Down_D24_S10_S30", "Up_D24_S20_S30", "Down_D24_S20_S30")

#Create a blank list
DEGS_Table<-c()

#Create a for loop that will go across all uploaded DEG files and extract out the number of DEGS
for (i in names(x)) {
  #Use the dataframe name to recall the dataframe and store it
  df <- x[[i]]
  #Extract the number of DEGS for which the FDR p value is less than 0.05 and not    N/A
  FDR_DEGS <- as.numeric(length( which( df$FDR.p.value < 0.05 & df$FDR.p.value!= "#N/A") ))
  #Extract the number of DEGS for which the Bonerroni p value is less than 0.05      and not N/A
  Bon_DEGS <- as.numeric(length( which( df$Bonferroni < 0.05 & df$Bonferroni!= "#N/A") ))
  #Add a row to the DEGS_Table that will have three columns: dataframe name and      the two values above
  DEGS_Table <- rbind(DEGS_Table,c(i,FDR_DEGS,Bon_DEGS))
}
#Define the column names for the DEGS_Table, make it a dataframe, and both sets of values numeric
colnames(DEGS_Table) <- c("df", "FDR_p_value", "Bonferonni_p_value")
DEGS_Table <- as.data.frame(DEGS_Table)
DEGS_Table$FDR_p_value <- as.numeric(DEGS_Table$FDR_p_value)
DEGS_Table$Bonferonni_p_value <- as.numeric(DEGS_Table$Bonferonni_p_value)

###PREPARE FOR GO AND KEGG ENRICHMENT ANALYSIS###

##Note: Although we focused on the Gene Set Enrichment Analysis results for this study, the tables made to prepare for GO and KEGG Analysis (Found in Exploratory_Transcriptome_Salinity_R_Scripts_Not_In_Manuscript.R) were still useful as references points, especially the Full_GO table which allowed us to search for relevant genes to our questions

##Get Gene Ontology (GO) Term information using BioMart and Ensembl database for target organism

#List the available biomarts
listEnsembl()

#Save the genes one in an object
ensembl <- useEnsembl(biomart = "genes")

#Search for your dataset using first letter of genus name followed by species name
searchDatasets(mart = ensembl, pattern = "sdumerili")

#Add the dataset name to your ensembl object
ensembl <- useDataset(dataset = "sdumerili_gene_ensembl", mart = ensembl)

#Previous steps in one step if know ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "sdumerili_gene_ensembl")

#use the mirrors object to go to a different one if its taking too long, options include: "useast", "uswest", "asia", and "www"

#Make a table of the attributes that can be extracted from this object
attributes = listAttributes(ensembl)
attributes[1:5,]

##Make some tables for summarizing multiple attributes for use later
#Make a table of the gene symbols associated with the entrez gene ids associated with your gene table
#Note you can only populate 3 attributes at a time
GENEID_SYM <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol', 'zfin_id_symbol' ), filters = 'entrezgene_id', mart = ensembl, values = rownames(gn_length_table))

#Make a table of the GO IDs and Names/Terms associated with the entrez gene ids associated with your gene table
GENEID_GO <- getBM(attributes = c('entrezgene_id','go_id', 'name_1006', 'definition_1006' ), filters = 'entrezgene_id', mart = ensembl, values = rownames(gn_length_table))

#Upload a table based upon the CLC Genomics transcripts file that compares the entrez gene name to the gene id
name2id <- read.csv("Genename_GeneID.csv")

#Extract out the unique rows
name2id <- unique(name2id)

#Make a copy of the GENEID_GO table so you can manipulate the column names to match the name2id table
TEST_GO <- GENEID_GO
colnames(TEST_GO) <- c("GeneID", "go_id", "name_1006", "definition_1006")

#Merge the two tables to add the gene ids based upon the CLC Genomics output to the GO table which will provide more information than if you just merged the GENEID_SYM table instead (which is based upon a different set of gene ids/symbols)
Full_GO <-  merge(TEST_GO, name2id, all = TRUE, by = "GeneID")

##Make a series of tables to use for GO and KEGG enrichment analysis

#Extract out the Entrez gene ids and their associated GO Names associated with the target organism's three letter designation for KEGG Enrichment Analysis

#The three letter designation can be found by going to the following website and searching for the scientific or common name 
#https://www.genome.jp/kegg/catalog/org_list.html
gene_pathway <- getGeneKEGGLinks("sdu")

#Get another table of the names associated with each of the pathways for KEGG Enrichment Analysis
pathway_names <- getKEGGPathwayNames("sdu", remove=TRUE)

#Create a full table of all GO IDs and their associated terms for GO Enrichment Analysis
GO.Name <- toTable(GOTERM)
GO.Name <- GO.Name[, c("go_id","Term")]

#Get a table of just the gene ids to the go ids for GO Enrichment Analysis
GENEID_GO_Terms <- getBM(attributes = c('entrezgene_id', 'go_id'), filters = 'entrezgene_id', mart = ensembl, values = rownames(gn_length_table))

#Get a table of just the gene ids to the go names/terms for your reference
GENEID_GO_Names <- getBM(attributes = c('entrezgene_id', 'name_1006' ), filters = 'entrezgene_id', mart = ensembl, values = rownames(gn_length_table))

###Gene Set Enrichment Analysis
#Because there was such a variation in the number of DEGs at the dph level and there were little DEGs at the pooled level, this lead to similar varied results in the GO and KEGG analysis (See Other Transcriptome Exploration Scripts Not In Manuscript.R). Thus it was decided to test all genes at once instead of just the DEGs using GSEA and jus to focus on the pooled trends.

##S10 vs S30
#Upload the full list of genes post CLC-Genomics DEGs analysis, get rid of lines with #N/A as the Fold Change, and make the fold change numeric
full_S10_S30 <- read.csv("Total All DPH 10 vs. 30 ppt.csv", stringsAsFactors=FALSE)
full_S10_S30 <- full_S10_S30[which(full_S10_S30$Fold.change!="#N/A"),]
full_S10_S30$Fold.change <- as.numeric(full_S10_S30$Fold.change)

#Prep for GSEA
# feature 1: create a numeric vector using the column with Fold Change values
geneList = full_S10_S30[,6]
# feature 2: name the values in this vector with the Gene Id column
names(geneList) = as.character(full_S10_S30[,10])
# feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

#Run the GSEA with normal parameters make a dataframe from the results
GSEA_10vs30 <- gseKEGG(geneList     = geneList,
                       organism     = 'sdu',
                       minGSSize    = 120,
                       pvalueCutoff = 0.05)
head(GSEA_10vs30)

df_GSEA_10vs30 <- as.data.frame(GSEA_10vs30)

##S10 vs S20
#Upload the full list of genes post CLC-Genomics DEGs analysis, get rid of lines with #N/A as the Fold Change, and make the fold change numeric
full_S10_S20 <- read.csv("Total All DPH 10 vs. 20 ppt.csv", stringsAsFactors=FALSE)
full_S10_S20 <- full_S10_S20[which(full_S10_S20$Fold.change!="#N/A"),]
full_S10_S20$Fold.change <- as.numeric(full_S10_S20$Fold.change)

#Prep for GSEA
# feature 1: create a numeric vector using the column with Fold Change values
geneList = full_S10_S20[,6]
# feature 2: name the values in this vector with the Gene Id column
names(geneList) = as.character(full_S10_S20[,10])
# feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

#Run the GSEA with normal parameters make a dataframe from the results
GSEA_10vs20 <- gseKEGG(geneList     = geneList,
                       organism     = 'sdu',
                       minGSSize    = 120,
                       pvalueCutoff = 0.05)
head(GSEA_10vs20)

df_GSEA_10vs20 <- as.data.frame(GSEA_10vs20)

##S20 vs S30
#Upload the full list of genes post CLC-Genomics DEGs analysis, get rid of lines with #N/A as the Fold Change, and make the fold change numeric

full_S20_S30 <- read.csv("Total All DPH 20 vs. 30 ppt.csv", stringsAsFactors=FALSE)
full_S20_S30 <- full_S20_S30[which(full_S20_S30$Fold.change!="#N/A"),]
full_S20_S30$Fold.change <- as.numeric(full_S20_S30$Fold.change)

#Prep for GSEA
# feature 1: create a numeric vector using the column with Fold Change values
geneList = full_S20_S30[,6]
# feature 2: name the values in this vector with the Gene Id column
names(geneList) = as.character(full_S20_S30[,10])
# feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

#Run the GSEA with normal parameters make a dataframe from the results
GSEA_20vs30 <- gseKEGG(geneList     = geneList,
                       organism     = 'sdu',
                       minGSSize    = 120,
                       pvalueCutoff = 0.05)
head(GSEA_20vs30)

df_GSEA_20vs30 <- as.data.frame(GSEA_20vs30)

##Plot the results

#S10 vs S30

#Create a new column denoting activated/suppressed status
df_GSEA_10vs30$Direction <- "Activated"

df_GSEA_10vs30$Direction[which(df_GSEA_10vs30$NES < 0)] <- "Suppressed"

##FIGURE 1

#Plot and publish the results
S10vsS30_GSEA_point_plot <- ggplot(df_GSEA_10vs30)+
  geom_point(mapping = aes(y= reorder(Description, NES), x=(rank), color = p.adjust, size=abs(NES), shape=Direction, fill = p.adjust),
             data = head(df_GSEA_10vs30[which(df_GSEA_10vs30$p.adjust < 0.05),], n = 37))+
  theme_bw()+
  scale_shape_manual(values=c(24,25))+
  scale_fill_continuous(low="red", high="blue") +
  scale_color_continuous(low="red", high="blue", guide="none") +
  labs(y = "Enriched Gene Sets", x = "rank", color = "p.adjust", size = "NES", shape = "Expression") +
  theme(axis.text=element_text(size=8)) +
  #ggtitle('24 DPH 10ppt vs 30 ppt') +
  theme(legend.key.size = unit(1.2, "line"))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
S10vsS30_GSEA_point_plot

tiff('Fig 1 Gene Sets Enriched between 10 ppt and 30 ppt.tiff', units="in", width=6, height=6, res=300)
S10vsS30_GSEA_point_plot
dev.off()

#S20 vs S30

#Create a new column denoting activated/suppressed status
df_GSEA_20vs30$Direction <- "Activated"

df_GSEA_20vs30$Direction[which(df_GSEA_20vs30$NES < 0)] <- "Suppressed"

##SUPP FIGURE 2

#Plot and publish the results
S20vsS30_GSEA_point_plot <- ggplot(df_GSEA_20vs30)+
  geom_point(mapping = aes(y= reorder(Description, NES), x=(rank), color = p.adjust, size=abs(NES), shape=Direction, fill = p.adjust),
             data = head(df_GSEA_20vs30[which(df_GSEA_20vs30$p.adjust < 0.05),], n = 37))+
  theme_bw() +
  scale_shape_manual(values=c(24,25))+
  scale_fill_continuous(low="red", high="blue") +
  scale_color_continuous(low="red", high="blue", guide="none") +
  labs(y = "Enriched Gene Sets", x = "rank", color = "p.adjust", size = "NES", shape = "Expression") +
  theme(axis.text=element_text(size=8)) +
  #ggtitle('24 DPH 10ppt vs 30 ppt') +
  theme(legend.key.size = unit(1.2, "line"))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
S20vsS30_GSEA_point_plot

tiff('Supp Fig 2 Gene Sets Enriched between 20 ppt and 30 ppt.tiff', units="in", width=6, height=6, res=300)
S20vsS30_GSEA_point_plot
dev.off()

#S10 vs S20

#Create a new column denoting activated/suppressed status
df_GSEA_10vs20$Direction <- "Activated"

df_GSEA_10vs20$Direction[which(df_GSEA_10vs20$NES < 0)] <- "Suppressed"

##SUPP FIGURE 3

#Plot and publish the results
S10vsS20_GSEA_point_plot <- ggplot(df_GSEA_10vs20)+
  geom_point(mapping = aes(y= reorder(Description, NES), x=(rank), color = p.adjust, size=abs(NES), shape=Direction, fill = p.adjust),
             data = head(df_GSEA_10vs20[which(df_GSEA_10vs20$p.adjust < 0.05),], n = 37))+
  theme_bw()+
  scale_shape_manual(values=c(24,25))+
  scale_fill_continuous(low="red", high="blue") +
  scale_color_continuous(low="red", high="blue", guide="none") +
  labs(y = "Enriched Gene Sets", x = "rank", color = "p.adjust", size = "NES", shape = "Expression") +
  theme(axis.text=element_text(size=8)) +
  #ggtitle('24 DPH 10ppt vs 30 ppt') +
  theme(legend.key.size = unit(1.2, "line"))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
S10vsS20_GSEA_point_plot

tiff('Supp Fig 3 Gene Sets Enriched between 10 ppt and 20 ppt.tiff', units="in", width=6, height=6, res=300)
S10vsS20_GSEA_point_plot
dev.off()






###GENE TREND STATISTICS AND FIGURES

#Copy data to manipulate it, change rownames to from gene ids to gene names, Transpose table, and add in metadata
VST_cts_boxplot <- data
rownames(VST_cts_boxplot) <- name2id$Gene.name[match(rownames(VST_cts_boxplot), name2id$GeneID)]
VST_cts_boxplot <- as.data.frame(t(VST_cts_boxplot))
VST_cts_boxplot <- cbind(VST_cts_boxplot, samples)

##Perform statistics on each target

##TABLE 5

##Osmoregulation

#atp1b4
#Salinity
kruskal.test(atp1b4 ~ Salinity, data=VST_cts_boxplot)
dunnTest(atp1b4 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(atp1b4 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(atp1b4 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab ab         ab
# 2   D15     ab ab         ab
# 3   D18      a b         a 
# 4    D2      a b         a 
# 5   D24     ab ab         ab
# 6    D3     ab ab         ab
# 7    D6      b a          b 

#"ab", "a", "ab", "ab", "b", "b", "ab"

##Figure 2A

atp1b4 <- ggplot(VST_cts_boxplot, aes(x=Salinity, y=atp1b4)) +   
  geom_boxplot(fill="gray90") + 
  #scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  theme_bw() +
  ggtitle(expression(atp1b4~-~Na^{"+"}:K^{"+"}~exchanger))+
  xlab("Salinity (ppt)") +
  scale_x_discrete(labels=c("10", "20", "30"))+
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(10.5, 11.5, 10.75), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
atp1b4

tiff('atp1b4 VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
atp1b4
dev.off()


#cftr
#Salinity
kruskal.test(cftr ~ Salinity, data=VST_cts_boxplot)
dunnTest(cftr ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(cftr ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(cftr ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab         ab
# 2   D15     ab         ab
# 3   D18      a         a 
# 4    D2      a         a 
# 5   D24      a         a 
# 6    D3      b          b
# 7    D6     ab         ab

#"ab", "a", "ab", "ab", "b", "b", "ab"

##Figure 2B

cftr <- ggplot(VST_cts_boxplot, aes(x=Salinity, y=cftr)) +   
  geom_boxplot(fill="gray90") + 
  #scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  theme_bw()+
  ggtitle(expression(cftr~-~Cl^{"-"}~transporter))+
  xlab("Salinity (ppt)") +
  scale_x_discrete(labels=c("10", "20", "30"))+
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(7, 7.75, 7.75), label = c("a", "ab", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
cftr

tiff('cftr VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
cftr
dev.off()

#slc26a3
#Salinity
kruskal.test(slc26a3 ~ Salinity, data=VST_cts_boxplot)
dunnTest(slc26a3 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(slc26a3 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(slc26a3 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab         ab
# 2   D15     ab         ab
# 3   D18      a         a 
# 4    D2      a         a 
# 5   D24      a         a 
# 6    D3      b          b
# 7    D6      b          b

#"ab", "a", "ab", "ab", "b", "b", "ab"

##Figure 2C

slc26a3 <- ggplot(VST_cts_boxplot, aes(x=Salinity, y=slc26a3)) +   
  geom_boxplot(fill="gray90") + 
  #scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  theme_bw()+
  ggtitle("slc26a3 - anion:anion transporter") +
  xlab("Salinity (ppt)") +
  scale_x_discrete(labels=c("10", "20", "30"))+
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(6, 7.3, 7.15), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
slc26a3

tiff('slc26a3 VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
slc26a3
dev.off()

#slc9a7
#Salinity
kruskal.test(slc9a7 ~ Salinity, data=VST_cts_boxplot)
dunnTest(slc9a7 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(slc9a7 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(slc9a7 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab        ab 
# 2   D15    abc        abc
# 3   D18     ab        ab 
# 4    D2      a        a  
# 5   D24      a        a  
# 6    D3      c          c
# 7    D6     bc         bc

#"ab", "a", "ab", "ab", "b", "b", "ab"

##Figure 2D

slc9a7 <- ggplot(VST_cts_boxplot, aes(x=Salinity, y=slc9a7)) +   
  geom_boxplot(fill="gray90") + 
  #scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  theme_bw()+
  ggtitle(expression(slc9a7~-~Na^{"+"}:H^{"+"}~antiporter))+
  xlab("Salinity (ppt)") +
  scale_x_discrete(labels=c("10", "20", "30"))+
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(6.6, 6.65, 6.7), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
slc9a7

tiff('slc9a7 VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
slc9a7
dev.off()

#slc26a6
#Salinity
kruskal.test(slc26a6 ~ Salinity, data=VST_cts_boxplot)
dunnTest(slc26a6 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(slc26a6 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(slc26a6 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12      a        a  
# 2   D15    abc        abc
# 3   D18     bc         bc
# 4    D2      b         b 
# 5   D24    abc        abc
# 6    D3    abc        abc
# 7    D6     ac        a c

#"ab", "a", "ab", "ab", "b", "b", "ab"

##Figure 2E

slc26a6 <- ggplot(VST_cts_boxplot, aes(x=Salinity, y=slc26a6)) +   
  geom_boxplot(fill="gray90") + 
  #scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  theme_bw()+
  ggtitle("slc26a6 - sulfate transporter") +
  xlab("Salinity (ppt)") +
  scale_x_discrete(labels=c("10", "20", "30"))+
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(6.5, 8.2, 7.1), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
slc26a6

tiff('slc26a6 VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
slc26a6
dev.off()

#atp6v0a1
#Salinity
kruskal.test(atp6v0a1 ~ Salinity, data=VST_cts_boxplot)
dunnTest(atp6v0a1 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(atp6v0a1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(atp6v0a1 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

##Figure 2F

atp6v0a1 <- ggplot(VST_cts_boxplot, aes(x=Salinity, y=atp6v0a1)) +   
  geom_boxplot(fill="gray90") + 
  #scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  theme_bw()+
  ggtitle(expression(atp6v0a1~-~V-ATPase~H^{"+"}~transporter))+
  xlab("Salinity (ppt)") +
  scale_x_discrete(labels=c("10", "20", "30"))+
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(8.35, 8.15, 8.2), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
atp6v0a1

tiff('atp6v0a1 VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
atp6v0a1
dev.off()

##Figure 2

#Make a multi-part figure, removing the axis labels from each separate figure
target_osmoregulation_gene_figure <- ggarrange(atp1b4 + rremove("ylab") + rremove("xlab"), cftr + rremove("ylab") + rremove("xlab"), slc26a3 + rremove("ylab") + rremove("xlab"), slc9a7 + rremove("ylab") + rremove("xlab"), slc26a6 + rremove("ylab") + rremove("xlab"), atp6v0a1 + rremove("ylab") + rremove("xlab"),
                                               labels = c("A", "B", "C", "D", "E", "F"),
                                               ncol = 2, nrow = 3,
                                               common.legend = TRUE)
target_osmoregulation_gene_figure

#Add annotations to the above figure
target_osmoregulation_gene_figure <- annotate_figure(target_osmoregulation_gene_figure, 
                                                     bottom = text_grob("Salinity (ppt)", size = 16),
                                                     left = text_grob("VST Gene Count", rot = 90, size =16))
target_osmoregulation_gene_figure

tiff('Figure 2 Target Osmoregulation Genes Figure.tiff', units="in", width=10, height=10, res=300)
target_osmoregulation_gene_figure
dev.off()

##Fatty Acid targets

#elovl5
#Salinity
kruskal.test(elovl5 ~ Salinity, data=VST_cts_boxplot)
dunnTest(elovl5 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(elovl5 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(elovl5 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab        ab 
# 2   D15    abc        abc
# 3   D18     ac        a c
# 4    D2      c          c
# 5   D24    abc        abc
# 6    D3      b         b 
# 7    D6      b         b 

##FIGURE T4A

elovl5 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=elovl5, fill=Feeding)) +   
  geom_boxplot() +
  # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("elovl5 - elongase") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(7.3, 9, 9, 10.1, 10.3, 10.1, 10), label = c("a", "ab", "bc", "bc", "c", "c", "bc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))  +
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
elovl5

tiff('elovl5 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
elovl5
dev.off()

#acadm
#Salinity
kruskal.test(acadm ~ Salinity, data=VST_cts_boxplot)
dunnTest(acadm ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(acadm ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(acadm ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12      a         a 
# 2   D15     ab         ab
# 3   D18      b          b
# 4    D2      b          b
# 5   D24     ab         ab
# 6    D3     ab         ab
# 7    D6      a         a

#"ab", "a", "a", "ab", "b", "b", "ab"

##FIGURE T4B

acadm <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=fads2, fill=Feeding)) +   geom_boxplot() +
  # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("acadm - fatty acid beta-oxidation") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(8.7,8.35, 10.5, 11.4, 11.7, 12, 11.25), label = c("ab", "a", "a", "ab", "b", "b", "ab"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
acadm

tiff('acadm VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
acadm
dev.off()
#fabp1
#Salinity
kruskal.test(fabp1 ~ Salinity, data=VST_cts_boxplot)
dunnTest(fabp1 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(fabp1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(fabp1 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab         ab
# 2   D15     ab         ab
# 3   D18     ab         ab
# 4    D2     ab         ab
# 5   D24      a         a 
# 6    D3      b          b
# 7    D6      b          b

##FIGURE T4C

fabp1 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=fabp1, fill=Feeding)) +   geom_boxplot() +
  # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("fabp1 - fatty acid binding") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(13, 12.65, 12.15, 13.35, 12.1, 11.6, 11.6), label = c("a", "a", "ab", "ab", "ab", "ab", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
fabp1

tiff('fabp1 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
fabp1
dev.off()

#fasn
#Salinity
kruskal.test(fasn ~ Salinity, data=VST_cts_boxplot)
dunnTest(fasn ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(fasn ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(fasn ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab        ab 
# 2   D15     ab        ab 
# 3   D18      a        a  
# 4    D2      a        a  
# 5   D24     ab        ab 
# 6    D3      c          c
# 7    D6     bc         bc

##FIGURE T4D

fasn <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=fasn, fill=Feeding)) +   geom_boxplot() +
  # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("fasn - fatty acid synthesis") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(8.1, 10.6, 11.55, 11.75, 12.3, 12.9, 12.3), label = c("a", "ab", "bc", "bc", "c", "c", "cb"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
fasn

tiff('fasn VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
fasn
dev.off()

#fads2
#Salinity
kruskal.test(fads2 ~ Salinity, data=VST_cts_boxplot)
dunnTest(fads2 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(fads2 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(fads2 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12    abc        abc
# 2   D15      a        a  
# 3   D18      a        a  
# 4    D2      a        a  
# 5   D24     ab        ab 
# 6    D3     bc         bc
# 7    D6      c          c

#"ab", "a", "abc", "c", "c", "c", "bc"

##FIGURE T4E

fads2 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=fads2, fill=Feeding)) +   geom_boxplot() +
  # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("fads2 - fatty acid desaturation") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(8.7, 8.4, 10.5, 11.4, 11.7, 12, 11.25), label = c("ab", "a", "abc", "c", "c", "c", "bc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
fads2

tiff('fads2 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
fads2
dev.off()

#elovl6
#Salinity
kruskal.test(elovl6 ~ Salinity, data=VST_cts_boxplot)
dunnTest(elovl6 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(elovl6 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(elovl6 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12      a         a 
# 2   D15     ab         ab
# 3   D18     ab         ab
# 4    D2      a         a 
# 5   D24      b          b
# 6    D3      b          b
# 7    D6      b          b

##FIGURE T4F

elovl6 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=elovl6, fill=Feeding)) +   geom_boxplot() +
  # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("elovl6 - elongase") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(8.1, 8.1, 8.8, 8.85, 9.25, 9.55, 8.2), label = c("a", "a", "b", "ab", "ab", "b", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
elovl6

tiff('elovl6 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
elovl6
dev.off()

##FIGURE T4

#Make a multi-part figure, removing the axis labels from each separate figure
target_fatty_acid_genes_figure <- ggarrange(elovl5 + rremove("ylab") + rremove("xlab"), acadm + rremove("ylab") + rremove("xlab"), fabp1 + rremove("ylab") + rremove("xlab"), fasn + rremove("ylab") + rremove("xlab"), fads2 + rremove("ylab") + rremove("xlab"), elovl6 + rremove("ylab") + rremove("xlab"),
                                            labels = c("A", "B", "C", "D", "E", "F"),
                                            ncol = 2, nrow = 3,
                                            common.legend = TRUE)
target_fatty_acid_genes_figure

#Add annotations to the above figure
target_fatty_acid_genes_figure <- annotate_figure(target_fatty_acid_genes_figure, 
                                                  bottom = text_grob("Days Post Hatch", size = 16),
                                                  left = text_grob("VST Gene Count", rot = 90, size =16))
target_fatty_acid_genes_figure

tiff('Figure 3 Target Fatty Acid Genes Figure.tiff', units="in", width=10, height=10, res=300)
target_fatty_acid_genes_figure
dev.off()

#Adjust Kruskal Wallis p-values for TABLE M5
gene_trends = c(
  "elovl5_Sal", "elovl5_DPH", 
  "elovl6_Sal", "elovl6_DPH", 
  "fapb1_Sal", "fapb1_DPH", 
  "fasn_Sal", "fasn_DPH", 
  "acadm_Sal", "acadm_DPH", 
  "fads2_Sal", "fads2_DPH", 
  "atp1b4_Sal", 
  "cftr_Sal", 
  "slc26a6_Sal",
  "slc26a3_Sal",
  "slc9a7_Sal", 
  "atp6v0a1_Sal",
  "sod1_Sal", "sod1_DPH", 
  "tlr3_Sal", "tlr3_DPH", 
  "nod1_Sal", "nod1_DPH", 
  "lyz_Sal", "lyz_DPH",
  "tnfaip3_Sal", "tnfaip3_DPH",
  "cat_Sal", "cat_DPH")
gene_trends=data.frame(gene_trends)
gene_trends$p.value <- c(0.7289, 2.298e-06, 0.8474, 2.205e-05, 0.858, 0.004525, 0.7896, 1.093e-07, 0.1417, 3.801e-05, 0.7958, 4.848e-08, 0.4114, 0.04045, 0.5025, 0.7175, 0.617, 0.854, 0.3055, 0.01173, 0.4264, 0.0004088, 0.01419, 1.889e-05, 0.9674, 0.000146, 0.4354, 7.815e-05, 0.4054, 5.033e-05)

gene_trends$padj <- p.adjust(gene_trends$p.value, method = "BH")
gene_trends

###VENN DIAGRAMS
##Transcript Level Salinity Venn Diagrams

#Subset each salinity and filter out transcripts without any sequences and make a list of transcripts
S10_TX_dataA <- subset_samples(TX_dataA, Salinity=="S10")
S10_TX_dataA # 34893 x 18
S10_TX_dataA <- filter_taxa(S10_TX_dataA, function(x) sum(x) >0, TRUE)
S10_TX_dataA # 31338 x 18
S10_txs <- as.list(rownames(otu_table(S10_TX_dataA)))

S20_TX_dataA <- subset_samples(TX_dataA, Salinity=="S20")
S20_TX_dataA # 34893 x 18
S20_TX_dataA <- filter_taxa(S20_TX_dataA, function(x) sum(x) >0, TRUE)
S20_TX_dataA # 31302 x 18
S20_txs <- as.list(rownames(otu_table(S20_TX_dataA)))

S30_TX_dataA <- subset_samples(TX_dataA, Salinity=="S30")
S30_TX_dataA # 34893 x 21
S30_TX_dataA <- filter_taxa(S30_TX_dataA, function(x) sum(x) >0, TRUE)
S30_TX_dataA # 31395 x 18
S30_txs <- as.list(rownames(otu_table(S30_TX_dataA)))

#Make a list of lists
TX_Salinities_Lists <- list('S10 (31338)' = S10_txs,
                            'S20 (31302)' = S20_txs,
                            'S30 (31395)' = S30_txs)

#Create the venn diagram
ggvenn(TX_Salinities_Lists, fill_color = c("deeppink", "midnightblue", "blue"), show_percentage = FALSE)+
  ggtitle("txs by Salinity") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Transcript level Salinity venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

##Gene Level Salinity Venn Diagrams

#Subset each salinity and filter out genes without any sequences and make a list of genes
S10_GN_dataA <- subset_samples(GN_dataA, Salinity=="S10")
S10_GN_dataA # 23800 x 18
S10_GN_dataA <- filter_taxa(S10_GN_dataA, function(x) sum(x) >0, TRUE)
S10_GN_dataA # 22157 x 18
S10_genes <- as.list(rownames(otu_table(S10_GN_dataA)))

S20_GN_dataA <- subset_samples(GN_dataA, Salinity=="S20")
S20_GN_dataA # 23800 x 18
S20_GN_dataA <- filter_taxa(S20_GN_dataA, function(x) sum(x) >0, TRUE)
S20_GN_dataA # 22143 x 18
S20_genes <- as.list(rownames(otu_table(S20_GN_dataA)))

S30_GN_dataA <- subset_samples(GN_dataA, Salinity=="S30")
S30_GN_dataA # 23800 x 21
S30_GN_dataA <- filter_taxa(S30_GN_dataA, function(x) sum(x) >0, TRUE)
S30_GN_dataA # 22172 x 18
S30_genes <- as.list(rownames(otu_table(S30_GN_dataA)))

Gene_Salinities_Lists <- list('S10 (22157)' = S10_genes,
                              'S20 (22143)' = S20_genes,
                              'S30 (22172)' = S30_genes)

ggvenn(Gene_Salinities_Lists, fill_color = c("deeppink", "midnightblue", "blue"), show_percentage = FALSE)+
  ggtitle("Genes by Salinity") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Gene level Salinity venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

###EXPLORING VST COUNT AND FATTY ACID PERCENTAGES CORRELATIONS

##Create a VST dataframe focused on the target fatty acid synthesis genes

#Extract out rows of target activity into new dataframe
fatty_acid_targets <- filter(Full_GO, Gene.name %in% c("acadm", "elovl6", "elovl5", "fasn", "fabp1", "fads2"))

#Extract out variance stabilized transformed counts from phyloseq object
data <- as.data.frame(otu_table(VST_GN_dataA))

#Keep only the rows that are found in the target dataframe and change rownames to Gene names based upon name2id table
fatty_acid_VST_cts <- data[rownames(data) %in% fatty_acid_targets$GeneID,]
rownames(fatty_acid_VST_cts) <- name2id$Gene.name[match(rownames(fatty_acid_VST_cts), name2id$GeneID)]
fatty_acid_VST_cts <- as.data.frame(t(fatty_acid_VST_cts))

write.csv(fatty_acid_VST_cts, "fatty_acid_VST_cts.csv")

#Upload the fatty acid data
NFA <- read.csv("NFA.csv", row.names = 1)

#The row names for these dataframes do not match
rownames(NFA)
rownames(fatty_acid_VST_cts)

#Remove the two samples from the VST counts dataframe without fatty acid data
fatty_acid_VST_cts <- subset(fatty_acid_VST_cts, rownames(fatty_acid_VST_cts) != "LST15DP3a")
fatty_acid_VST_cts <- subset(fatty_acid_VST_cts, rownames(fatty_acid_VST_cts) != "LST24DP1a")


##Take preliminary look at EPA
#Combine into one dataframe
fatty_acid_VST_cts_NFA <- cbind(fatty_acid_VST_cts, NFA)

#Plot trends with correlations
ggscatter(data = fatty_acid_VST_cts_NFA, x = "elovl5", y = "C20.5..EPA.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "elovl5 expression", ylab = "EPA")

ggplot(fatty_acid_VST_cts_NFA, aes(x=elovl5, y=C20.5..EPA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)

ggscatter(data = fatty_acid_VST_cts_NFA, x = "fads2", y = "C20.5..EPA.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "elovl5 expression", ylab = "EPA")

ggplot(fatty_acid_VST_cts_NFA, aes(x=fads2, y=C20.5..EPA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)

##Take preliminary look at DHA

#Plot trends with correlations
ggscatter(data = fatty_acid_VST_cts_NFA, x = "elovl5", y = "C22.6..DHA.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "elovl5 expression", ylab = "DHA")

ggplot(fatty_acid_VST_cts_NFA, aes(x=elovl5, y=C22.6..DHA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)

ggscatter(data = fatty_acid_VST_cts_NFA, x = "fads2", y = "C22.6..DHA.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "elovl5 expression", ylab = "DHA")

ggplot(fatty_acid_VST_cts_NFA, aes(x=fads2, y=C22.6..DHA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)

##Take preliminary look at ARA

#Plot trends with correlations
ggscatter(data = fatty_acid_VST_cts_NFA, x = "elovl5", y = "C20.4N6..ARA.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "elovl5 expression", ylab = "EPA")

ggplot(fatty_acid_VST_cts_NFA, aes(x=elovl5, y=C20.4N6..ARA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)

ggscatter(data = fatty_acid_VST_cts_NFA, x = "fads2", y = "C20.4N6..ARA.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "elovl5 expression", ylab = "EPA")

ggplot(fatty_acid_VST_cts_NFA, aes(x=fads2, y=C20.4N6..ARA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)


#Note: We noted a trend that the amounts of fatty acids associated with the 24dph samples were typically major outliers, this is likely because the larvae were likely just feed pellets at that time which are designed to have a large amount of fatty acids thus skewing results. The decision was made to continue analysis without those samples. 

#Remove the 24 dph samples due to them being an outlier
fatty_acid_VST_cts_no24DPH <- subset(fatty_acid_VST_cts, rownames(fatty_acid_VST_cts) != "LST24DP4a")
fatty_acid_VST_cts_no24DPH <- subset(fatty_acid_VST_cts_no24DPH, rownames(fatty_acid_VST_cts_no24DPH) != "LST24DP5a")
fatty_acid_VST_cts_no24DPH <- subset(fatty_acid_VST_cts_no24DPH, rownames(fatty_acid_VST_cts_no24DPH) != "LST24DP6a")
fatty_acid_VST_cts_no24DPH <- subset(fatty_acid_VST_cts_no24DPH, rownames(fatty_acid_VST_cts_no24DPH) != "LST24DP7a")
fatty_acid_VST_cts_no24DPH <- subset(fatty_acid_VST_cts_no24DPH, rownames(fatty_acid_VST_cts_no24DPH) != "LST24DP8a")
fatty_acid_VST_cts_no24DPH <- subset(fatty_acid_VST_cts_no24DPH, rownames(fatty_acid_VST_cts_no24DPH) != "LST24DP9a")

#Remove the 24 dph samples due to them being an outlier
NFA_no24DPH <- subset(NFA, rownames(NFA) != "LST24DP4a")
NFA_no24DPH <- subset(NFA_no24DPH, rownames(NFA_no24DPH) != "LST24DP5a")
NFA_no24DPH <- subset(NFA_no24DPH, rownames(NFA_no24DPH) != "LST24DP6a")
NFA_no24DPH <- subset(NFA_no24DPH, rownames(NFA_no24DPH) != "LST24DP7a")
NFA_no24DPH <- subset(NFA_no24DPH, rownames(NFA_no24DPH) != "LST24DP8a")
NFA_no24DPH <- subset(NFA_no24DPH, rownames(NFA_no24DPH) != "LST24DP9a")

#Combine into one dataframe
fatty_acid_VST_cts_no24DPH_NFA <- cbind(fatty_acid_VST_cts_no24DPH, NFA_no24DPH)


###Correlations###

#Upload table of relative abundances of target fatty acids. Here we were interested in total fatty acid (TFA), neutral fatty acid (NFA), and polar fatty acid (PFA) fractions of linoleic acid (LA), linolenic acid (LNA), arachidonic acid (ARA), (eicosapentaenoic acid EPA), and docosahexaenoic acid (DHA).
FA <- read.csv("Target_Fatty_Acids.csv", row.names = 1)

#Preliminary testing revealed that 24 dph was an outlier, thus it will not be tested
FA_no24DPH <- subset(FA, rownames(FA) != "LST24DP4a")
FA_no24DPH <- subset(FA_no24DPH, rownames(FA_no24DPH) != "LST24DP5a")
FA_no24DPH <- subset(FA_no24DPH, rownames(FA_no24DPH) != "LST24DP6a")
FA_no24DPH <- subset(FA_no24DPH, rownames(FA_no24DPH) != "LST24DP7a")
FA_no24DPH <- subset(FA_no24DPH, rownames(FA_no24DPH) != "LST24DP8a")
FA_no24DPH <- subset(FA_no24DPH, rownames(FA_no24DPH) != "LST24DP9a")

##Conduct correlation tests between two fatty acid genes and the five target fatty acids from three different fractions

#Record rho and p-value in an Excel sheet to be used for p-value adjustment later

#TABLE M6

##elovl5
##NFAs vs elovl5 
#NFA-EPA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Neutral.C20.5..EPA., method ="spearman") #R 0.6866096  p=5.132e-08

#NFA-DHA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Neutral.C22.6..DHA., method ="spearman") #R -0.4063265 p=0.004025

#NFA-ARA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Neutral.C20.4N6..ARA., method ="spearman") #R 0.5747959  p=2.235e-05

#NFA-LA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Neutral.C18.2..LA., method ="spearman") #R -0.4473469 p=0.001413

#NFA-LNA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Neutral.C18.3..LNA., method ="spearman") #R -0.5788776  p=1.906e-05

##PFAs vs elovl5 
#PFA-EPA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Polar.C20.5..EPA., method ="spearman") #R 0.6316017   p=1.136e-06

#PFA-DHA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Polar.C22.6..DHA., method ="spearman") #R -0.2456373 p=0.0889

#PFA-ARA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Polar.C20.4N6..ARA., method ="spearman") #R 0.3336969   p=0.01912

#PFA-LA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Polar.C18.2..LA., method ="spearman") #R 0.3674299  p=0.0094

#PFA-LNA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Polar.C18.3..LNA., method ="spearman") #R 0.1701009   p=0.2426

##TFAs vs elovl5 
#TFA-EPA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Total.C20.5..EPA., method ="spearman") #R 0.709064   p=1.185e-08

#TFA-DHA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Total.C22.6..DHA., method ="spearman") #R -0.4380612 p=0.00181

#TFA-ARA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Total.C20.4N6..ARA., method ="spearman") #R 0.4226531   0.002692

#TFA-LA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Total.C18.2..LA., method ="spearman") #R -0.4302041 p=0.002221

#TFA-LNA vs elovl5
cor.test(fatty_acid_VST_cts_no24DPH$elovl5, FA_no24DPH$Total.C18.3..LNA., method ="spearman") #R -0.6378571  p=1.581e-06


##fads2
##NFAs vs fads2 
#NFA-EPA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Neutral.C20.5..EPA., method ="spearman") #R 0.7761181   p=5.702e-11

#NFA-DHA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Neutral.C22.6..DHA., method ="spearman") #R -0.4308163 p=0.002186

#NFA-ARA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Neutral.C20.4N6..ARA., method ="spearman") #R 0.6702041   p=3.643e-07

#NFA-LA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Neutral.C18.2..LA., method ="spearman") #R -0.5280612 p=0.000122

#NFA-LNA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Neutral.C18.3..LNA., method ="spearman") #R -0.7338776  p=1.217e-08

##PFAs vs fads2 
#PFA-EPA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Polar.C20.5..EPA., method ="spearman") #R 0.756844   p=3.146e-10

#PFA-DHA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Polar.C22.6..DHA., method ="spearman") #R -0.1647107  p=0.2581

#PFA-ARA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Polar.C20.4N6..ARA., method ="spearman") #R 0.5335722   p=7.886e-05

#PFA-LA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Polar.C18.2..LA., method ="spearman") #R 0.5752157  p=1.537e-05

#PFA-LNA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Polar.C18.3..LNA., method ="spearman") #R 0.3056301   p=0.03271

##TFAs vs fads2 
#TFA-EPA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Total.C20.5..EPA., method ="spearman") #R 0.7825042   p=3.119e-11

#TFA-DHA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Total.C22.6..DHA., method ="spearman") #R -0.472551  p=0.0006977

#TFA-ARA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Total.C20.4N6..ARA., method ="spearman") #R 0.4882653   p=0.0004378

#TFA-LA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Total.C18.2..LA., method ="spearman") #R -0.5511224 p=0.001413

#TFA-LNA vs fads2
cor.test(fatty_acid_VST_cts_no24DPH$fads2, FA_no24DPH$Total.C18.3..LNA., method ="spearman") #R -0.7968367  p=2.2e-16

##Make Correlation Figures

#Combine target VST counts and target fatty acid percentages dataframes into one dataframe
fatty_acid_VST_cts_FA_RPs_no24DPH <- cbind(fatty_acid_VST_cts_no24DPH, FA_no24DPH)

##elovl5
##NFAs vs elovl5

##SUPP FIGURE T4A

#NFA-LA vs elovl5

NFA_LA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Neutral.C18.2..LA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("NFA LA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("NFA - C18:2 (LA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_LA_elovl5_JA

tiff('NFA LA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
NFA_LA_elovl5_JA
dev.off()

##SUPP FIGURE T4B

#NFA-LNA vs elovl5

NFA_LNA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Neutral.C18.3..LNA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("LNA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("NFA - C20:5 (LNA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_LNA_elovl5_JA

tiff('NFA LNA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
NFA_LNA_elovl5_JA
dev.off()

##SUPP FIGURE T4C

#NFA-ARA vs elovl5

NFA_ARA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Neutral.C20.4N6..ARA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("NFA ARA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("NFA - C20:4 (ARA)") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_ARA_elovl5_JA

tiff('NFA ARA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
NFA_ARA_elovl5_JA
dev.off()

##SUPP FIGURE T4D

#NFA-EPA vs elovl5
NFA_EPA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Neutral.C20.5..EPA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("NFA EPA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("NFA - C20:5 (EPA)") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_EPA_elovl5_JA

tiff('NFA EPA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
NFA_EPA_elovl5_JA
dev.off()

##SUPP FIGURE T4E

#NFA-DHA vs elovl5

NFA_DHA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Neutral.C22.6..DHA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("NFA DHA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("NFA - C22:6 (DHA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_DHA_elovl5_JA

tiff('NFA DHA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
NFA_DHA_elovl5_JA
dev.off()

#SUPP FIGURE T4

elovl5_vs_NFAs_figure_JA <- ggarrange(NFA_LA_elovl5_JA , NFA_LNA_elovl5_JA , NFA_ARA_elovl5_JA , NFA_EPA_elovl5_JA , NFA_DHA_elovl5_JA,
                                      labels = c("A", "B", "C", "D", "E"),
                                      ncol = 2, nrow = 3,
                                      common.legend = TRUE)
elovl5_vs_NFAs_figure_JA

tiff('Supp Fig 4 elovl5 vs NFAs JA.tiff', units="in", width=10, height=10, res=300)
elovl5_vs_NFAs_figure_JA
dev.off()

##PFAs vs elovl5

#SUPP FIGURE T5A

#PFA-LA vs elovl5

PFA_LA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Polar.C18.2..LA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("PFA LA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("PFA - C18:2 (LA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_LA_elovl5_JA

tiff('PFA LA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
PFA_LA_elovl5_JA
dev.off()

#SUPP FIGURE T5B

#PFA-LNA vs elovl5

PFA_LNA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Polar.C18.3..LNA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("LNA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("PFA - C20:5 (LNA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_LNA_elovl5_JA

tiff('PFA LNA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
PFA_LNA_elovl5_JA
dev.off()

#SUPP FIGURE T5C

#PFA-ARA vs elovl5

PFA_ARA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Polar.C20.4N6..ARA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("PFA ARA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("PFA - C20:4 (ARA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_ARA_elovl5_JA

tiff('PFA ARA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
PFA_ARA_elovl5_JA
dev.off()

#SUPP FIGURE T5D

#PFA-EPA vs elovl5
PFA_EPA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Polar.C20.5..EPA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("PFA EPA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("PFA - C20:5 (EPA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_EPA_elovl5_JA

tiff('PFA EPA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
PFA_EPA_elovl5_JA
dev.off()

#SUPP FIGURE T5E

#PFA-DHA vs elovl5

PFA_DHA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Polar.C22.6..DHA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("PFA DHA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("PFA - C22:6 (DHA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_DHA_elovl5_JA

tiff('PFA DHA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
PFA_DHA_elovl5_JA
dev.off()

#SUPP FIGURE T5

elovl5_vs_PFAs_figure_JA <- ggarrange(PFA_LA_elovl5_JA , PFA_LNA_elovl5_JA , PFA_ARA_elovl5_JA , PFA_EPA_elovl5_JA , PFA_DHA_elovl5_JA,
                                      labels = c("A", "B", "C", "D", "E"),
                                      ncol = 2, nrow = 3,
                                      common.legend = TRUE)
elovl5_vs_PFAs_figure_JA

tiff('Supp Fig 5 elovl5 vs PFAs JA.tiff', units="in", width=10, height=10, res=300)
elovl5_vs_PFAs_figure_JA
dev.off()

##Total Fatty Acids vs elovl5

#SUPP FIGURE T6A

#TFA-LA vs elovl5

TFA_LA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Total.C18.2..LA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("TFA LA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("TFA - C18:2 (LA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_LA_elovl5_JA

tiff('TFA LA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
TFA_LA_elovl5_JA
dev.off()

#SUPP FIGURE T6B

#TFA-LNA vs elovl5

TFA_LNA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Total.C18.3..LNA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("LNA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("TFA - C20:5 (LNA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_LNA_elovl5_JA

tiff('TFA LNA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
TFA_LNA_elovl5_JA
dev.off()

#SUPP FIGURE T6C

#TFA-ARA vs elovl5

TFA_ARA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Total.C20.4N6..ARA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("TFA ARA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("TFA - C20:4 (ARA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_ARA_elovl5_JA

tiff('TFA ARA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
TFA_ARA_elovl5_JA
dev.off()

#SUPP FIGURE T6D

#TFA-EPA vs elovl5
TFA_EPA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Total.C20.5..EPA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("TFA EPA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("TFA - C20:5 (EPA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_EPA_elovl5_JA

tiff('TFA EPA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
TFA_EPA_elovl5_JA
dev.off()

#SUPP FIGURE T6E

#TFA-DHA vs elovl5

TFA_DHA_elovl5_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=elovl5, y=Total.C22.6..DHA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("TFA DHA vs elovl5") +
  xlab("elovl5 - elongase") +
  ylab("TFA - C22:6 (DHA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_DHA_elovl5_JA

tiff('TFA DHA vs elovl5 JA.tiff', units="in", width=6, height=6, res=300)
TFA_DHA_elovl5_JA
dev.off()

#SUPP FIGURE T6

elovl5_vs_TFAs_figure_JA <- ggarrange(TFA_LA_elovl5_JA , TFA_LNA_elovl5_JA , TFA_ARA_elovl5_JA , TFA_EPA_elovl5_JA , TFA_DHA_elovl5_JA,
                                      labels = c("A", "B", "C", "D", "E"),
                                      ncol = 2, nrow = 3,
                                      common.legend = TRUE)
elovl5_vs_TFAs_figure_JA

tiff('Supp Fig 6 elovl5 vs TFAs JA.tiff', units="in", width=10, height=10, res=300)
elovl5_vs_TFAs_figure_JA
dev.off()

##fads2
##NFAs vs fads2

#SUPP FIGURE T7D

#NFA-EPA vs fads2
NFA_EPA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Neutral.C20.5..EPA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("NFA EPA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("NFA - C20:5 (EPA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_EPA_fads2_JA

tiff('NFA EPA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
NFA_EPA_fads2_JA
dev.off()

#SUPP FIGURE T7E

#NFA-DHA vs fads2

NFA_DHA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Neutral.C22.6..DHA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("NFA DHA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("NFA - C22:6 (DHA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_DHA_fads2_JA

tiff('NFA DHA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
NFA_DHA_fads2_JA
dev.off()

#SUPP FIGURE T7C

#NFA-ARA vs fads2

NFA_ARA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Neutral.C20.4N6..ARA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("NFA ARA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("NFA - C20:4 (ARA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_ARA_fads2_JA

tiff('NFA ARA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
NFA_ARA_fads2_JA
dev.off()

#SUPP FIGURE T7A

#NFA-LA vs fads2

NFA_LA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Neutral.C18.2..LA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("NFA LA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("NFA - C18:2 (LA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_LA_fads2_JA

tiff('NFA LA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
NFA_LA_fads2_JA
dev.off()

#SUPP FIGURE T7B

#NFA-LNA vs fads2

NFA_LNA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Neutral.C18.3..LNA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("LNA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("NFA - C20:5 (LNA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
NFA_LNA_fads2_JA

tiff('NFA LNA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
NFA_LNA_fads2_JA
dev.off()

#SUPP FIGURE T7

fads2_vs_NFAs_figure_JA <- ggarrange(NFA_LA_fads2_JA , NFA_LNA_fads2_JA , NFA_ARA_fads2_JA , NFA_EPA_fads2_JA , NFA_DHA_fads2_JA,
                                     labels = c("A", "B", "C", "D", "E"),
                                     ncol = 2, nrow = 3,
                                     common.legend = TRUE)
fads2_vs_NFAs_figure_JA

tiff('Supp Fig 7 fads2 vs NFAs JA.tiff', units="in", width=10, height=10, res=300)
fads2_vs_NFAs_figure_JA
dev.off()

##PFAs vs fads2

#SUPP FIGURE T8D

#PFA-EPA vs fads2
PFA_EPA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Polar.C20.5..EPA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("PFA EPA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("PFA - C20:5 (EPA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_EPA_fads2_JA

tiff('PFA EPA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
PFA_EPA_fads2_JA
dev.off()

#SUPP FIGURE T8E

#PFA-DHA vs fads2

PFA_DHA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Polar.C22.6..DHA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("PFA DHA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("PFA - C22:6 (DHA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_DHA_fads2_JA

tiff('PFA DHA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
PFA_DHA_fads2_JA
dev.off()

#SUPP FIGURE T8C

#PFA-ARA vs fads2

PFA_ARA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Polar.C20.4N6..ARA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("PFA ARA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("PFA - C20:4 (ARA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_ARA_fads2_JA

tiff('PFA ARA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
PFA_ARA_fads2_JA
dev.off()

#SUPP FIGURE T8A

#PFA-LA vs fads2

PFA_LA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Polar.C18.2..LA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("PFA LA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("PFA - C18:2 (LA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_LA_fads2_JA

tiff('PFA LA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
PFA_LA_fads2_JA
dev.off()

#SUPP FIGURE T8B

#PFA-LNA vs fads2

PFA_LNA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Polar.C18.3..LNA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("LNA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("PFA - C20:5 (LNA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
PFA_LNA_fads2_JA

tiff('PFA LNA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
PFA_LNA_fads2_JA
dev.off()

#SUPP FIGURE T8

fads2_vs_PFAs_figure_JA <- ggarrange(PFA_LA_fads2_JA , PFA_LNA_fads2_JA , PFA_ARA_fads2_JA , PFA_EPA_fads2_JA , PFA_DHA_fads2_JA,
                                     labels = c("A", "B", "C", "D", "E"),
                                     ncol = 2, nrow = 3,
                                     common.legend = TRUE)
fads2_vs_PFAs_figure_JA

tiff('Supp Fig 8 fads2 vs PFAs JA.tiff', units="in", width=10, height=10, res=300)
fads2_vs_PFAs_figure_JA
dev.off()

##Total Fatty Acids vs fads2

#SUPP FIGURE T9D

#TFA-EPA vs fads2
TFA_EPA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Total.C20.5..EPA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("TFA EPA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("TFA - C20:5 (EPA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_EPA_fads2_JA

tiff('TFA EPA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
TFA_EPA_fads2_JA
dev.off()

#SUPP FIGURE T9E

#TFA-DHA vs fads2

TFA_DHA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Total.C22.6..DHA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("TFA DHA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("TFA - C22:6 (DHA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_DHA_fads2_JA

tiff('TFA DHA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
TFA_DHA_fads2_JA
dev.off()

#SUPP FIGURE T9C

#TFA-ARA vs fads2

TFA_ARA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Total.C20.4N6..ARA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("TFA ARA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("TFA - C20:4 (ARA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_ARA_fads2_JA

tiff('TFA ARA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
TFA_ARA_fads2_JA
dev.off()

#SUPP FIGURE T9A

#TFA-LA vs fads2

TFA_LA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Total.C18.2..LA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("TFA LA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("TFA - C18:2 (LA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_LA_fads2_JA

tiff('TFA LA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
TFA_LA_fads2_JA
dev.off()

#SUPP FIGURE T9B

#TFA-LNA vs fads2

TFA_LNA_fads2_JA <- ggplot(fatty_acid_VST_cts_FA_RPs_no24DPH, aes(x=fads2, y=Total.C18.3..LNA.)) + geom_point(aes(color=Days_Post_Hatch))+
  geom_smooth(method=lm)+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("3", "6", "12", "15", "18", "20", "24"))+
  ggtitle("LNA vs fads2") +
  xlab("fads2 - desaturase") +
  ylab("TFA - C20:5 (LNA)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank())
TFA_LNA_fads2_JA

tiff('TFA LNA vs fads2 JA.tiff', units="in", width=6, height=6, res=300)
TFA_LNA_fads2_JA
dev.off()

#SUPP FIGURE T9

fads2_vs_TFAs_figure_JA <- ggarrange(TFA_LA_fads2_JA , TFA_LNA_fads2_JA , TFA_ARA_fads2_JA , TFA_EPA_fads2_JA , TFA_DHA_fads2_JA,
                                     labels = c("A", "B", "C", "D", "E"),
                                     ncol = 2, nrow = 3,
                                     common.legend = TRUE)
fads2_vs_TFAs_figure_JA

tiff('Supp Fig 9 fads2 vs TFAs JA.tiff', units="in", width=10, height=10, res=300)
fads2_vs_TFAs_figure_JA
dev.off()

#Upload the Excel sheet you recorded all the data in and adjust the pvalues
correlations <- read.csv("Correlations.csv")
head(correlations)
correlations$padj <- p.adjust(correlations$p.value, method = "BH")
head(correlations)

write.csv(correlations, "correlations_padj.csv")
