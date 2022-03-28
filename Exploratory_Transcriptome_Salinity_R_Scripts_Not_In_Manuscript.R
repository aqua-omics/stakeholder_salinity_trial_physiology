#Extract out all information associated with a particular GO name 
elongases <- filter(Full_GO, name_1006=="fatty acid elongase activity")

#Extract out target genes and change rownames to gene names
fatty_acid_VST_cts_boxplot <- data[rownames(data) %in% elongases$GeneID,]
rownames(fatty_acid_VST_cts_boxplot) <- elongases$Gene.name


#elovl7
#Salinity
kruskal.test(elovl7 ~ Salinity, data=VST_cts_boxplot)
dunnTest(elovl7 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(elovl7 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(elovl7 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

#None

elovl7 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=elovl7, color=Feeding)) +   geom_boxplot() +
  # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  ggtitle("elovl7") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
elovl7

tiff('elovl7 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
elovl7
dev.off()


###MAKING A HEATMAP OF TARGET GENES###

##Make a heatmap for a particular GO name

#Extract out Entrez Gene Ids for a particular GO name
lipid_GENEIDs <- as.data.frame(GENEID_GO$entrezgene_id[GENEID_GO$name_1006=="fatty acid elongase activity"])
colnames(lipid_GENEIDs) <- "GENEID"

#Extract out the VST abundance data from the phyloseq object
data <- as.data.frame(otu_table(VST_GN_dataA))

#Use the previous table to extract out rows in data for particular GO name
data_go_lipids <- data[rownames(data) %in% lipid_GENEIDs$GENEID,]

#Create a clustering object based upon a distance matrix for target genes and plot horizontally
my_hclust_gene <- hclust(dist(data_go_lipids), method = "complete")

#Visualize the clustering as a dendrogram
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)
#Create a list of integers corresponding to the cluster number
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)

#Convert this data to a dataframe and use ifelse to relabel the information by add a second column denoting what each result is equal to
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))

#Save the metadata table as a dataframe and keep only target columns
META2 <- as.data.frame(META)
META_sub <- select(samples, Salinity, Days_Post_Hatch)

#Create a heatmap of the target genes with rows annotated with the clustering information and columns annotated with the metadata information
pheatmap(data_go_lipids, annotation_row = my_gene_col, annotation_col = META_sub)

#Do not cluster the heatmap, instead it will be ordered by sample name
pheatmap(data_go_lipids, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = META_sub)


##Make a heatmap for a certain list of genes
#Extract out rows of target activity into new dataframe
fatty_acid_targets <- filter(Full_GO, Gene.name %in% c("acadm", "elovl6", "elovl5", "fasn", "fabp1", "fads2"))

#Extract out variance stabilized transformed counts from phyloseq object
data <- as.data.frame(otu_table(VST_GN_dataA))

#Keep only the rows that are found in the target dataframe and change rownames to Gene names based upon name2id table
fatty_acid_VST_cts_heatmap <- data[rownames(data) %in% fatty_acid_targets$GeneID,]
rownames(fatty_acid_VST_cts_heatmap) <- name2id$Gene.name[match(rownames(fatty_acid_VST_cts_heatmap), name2id$GeneID)]

#Select columns from the metadata table to annotate the heatmap and change column names
META_sub <- select(samples, Salinity, Days_Post_Hatch)
colnames(META_sub) <- c("Salinity", "Days Post Hatch")

#Make a list of lists to tell it what colors to use
my_colour = list(
  "Salinity" = c(S10 = "pink", S20 = "lightblue", S30 = "green"),
  "Days Post Hatch" = c(D3 = "rosybrown", D6 = "red", D12 = "forestgreen", D15 = "brown", D18 = "purple", D20 = "dark gray", D24 = "orange")
)

#Make your heat map
fatty_acid_heatmap <- pheatmap(fatty_acid_VST_cts_heatmap, annotation_col = META_sub, cluster_cols = FALSE, cluster_rows = FALSE, annotation_colors = my_colour)
fatty_acid_heatmap

#Publish your heat map
tiff('Fatty Acid Gene Heatmap.tiff', units="in", width=12, height=8, res=300)
fatty_acid_heatmap
dev.off()

###GO TERM TRENDS###

##Summarize the total VST Counts for a particular GO Term (Ex: fatty acid elongase activity)

#Extract out Entrez Gene Ids for a particular GO name
lipid_GENEIDs <- as.data.frame(GENEID_GO$entrezgene_id[GENEID_GO$name_1006=="fatty acid elongase activity"])
colnames(lipid_GENEIDs) <- "GENEID"

#Extract out all information associated with a particular GO name 
elongases <- filter(Full_GO, name_1006=="fatty acid elongase activity")

#Extract out abundances of genes using their Gene Id and transpose it
lipid_abundances <- data[rownames(data) %in% elongases$GeneID,]
lipid_abundances_V2 <- as.data.frame(t(lipid_abundances))

#Add column with totals of the genes
lipid_abundances_V2 <- cbind(lipid_abundances_V2, Total = rowSums(lipid_abundances_V2))

#Add in metadata
lipid_abundances_V2 <- cbind(lipid_abundances_V2, samples)

#Make basic bar plot of Total associated with particular GO Name
ggplot(data=lipid_abundances_V2, aes(x=rownames(lipid_abundances_V2), y=Total)) +
  geom_bar(stat="identity") +
  theme(axis.text.x  = element_text(angle=45, hjust=1))

#Make boxplot exploring the Total by Metadata
ggplot(lipid_abundances_V2, aes(x=Salinity, y=Total, color=Salinity)) +   geom_boxplot() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels = c("10", "20", "30"))+
  scale_x_discrete(labels=c("10", "20", "30"))+
  #ggtitle("Tissue PPG% by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Total Variance Stabilized Transformed Elongase Counts") +
  theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))

###KEGG ENRICHMENT ANALYSIS WITH LIMMA LIBRARY###

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

#Conduct KEGG enrichment analysis using each set of these objects
keg_Up_D3_S10_S30 <- kegga(Up_D3_S10_S30$GeneID, species.KEGG ="sdu")
keg_Down_D3_S10_S30 <- kegga(Down_D3_S10_S30$GeneID, species.KEGG ="sdu")
keg_Up_D3_S10_S20 <- kegga(Up_D3_S10_S20$GeneID, species.KEGG ="sdu")
keg_Down_D3_S10_S20 <- kegga(Down_D3_S10_S20$GeneID, species.KEGG ="sdu")
keg_Up_D3_S20_S30 <- kegga(Up_D3_S20_S30$GeneID, species.KEGG ="sdu")
keg_Down_D3_S20_S30 <- kegga(Down_D3_S20_S30$GeneID, species.KEGG ="sdu")

keg_Up_D6_S10_S30 <- kegga(Up_D6_S10_S30$GeneID, species.KEGG ="sdu")
keg_Down_D6_S10_S30 <- kegga(Down_D6_S10_S30$GeneID, species.KEGG ="sdu")
keg_Up_D6_S10_S20 <- kegga(Up_D6_S10_S20$GeneID, species.KEGG ="sdu")
keg_Down_D6_S10_S20 <- kegga(Down_D6_S10_S20$GeneID, species.KEGG ="sdu")
keg_Up_D6_S20_S30 <- kegga(Up_D6_S20_S30$GeneID, species.KEGG ="sdu")
keg_Down_D6_S20_S30 <- kegga(Down_D6_S20_S30$GeneID, species.KEGG ="sdu")

keg_Up_D12_S10_S30 <- kegga(Up_D12_S10_S30$GeneID, species.KEGG ="sdu")
keg_Down_D12_S10_S30 <- kegga(Down_D12_S10_S30$GeneID, species.KEGG ="sdu")
keg_Up_D12_S10_S20 <- kegga(Up_D12_S10_S20$GeneID, species.KEGG ="sdu")
keg_Down_D12_S10_S20 <- kegga(Down_D12_S10_S20$GeneID, species.KEGG ="sdu")
keg_Up_D12_S20_S30 <- kegga(Up_D12_S20_S30$GeneID, species.KEGG ="sdu")
keg_Down_D12_S20_S30 <- kegga(Down_D12_S20_S30$GeneID, species.KEGG ="sdu")

keg_Up_D15_S10_S30 <- kegga(Up_D15_S10_S30$GeneID, species.KEGG ="sdu")
keg_Down_D15_S10_S30 <- kegga(Down_D15_S10_S30$GeneID, species.KEGG ="sdu")
keg_Up_D15_S10_S20 <- kegga(Up_D15_S10_S20$GeneID, species.KEGG ="sdu")
keg_Down_D15_S10_S20 <- kegga(Down_D15_S10_S20$GeneID, species.KEGG ="sdu")
keg_Up_D15_S20_S30 <- kegga(Up_D15_S20_S30$GeneID, species.KEGG ="sdu")
keg_Down_D15_S20_S30 <- kegga(Down_D15_S20_S30$GeneID, species.KEGG ="sdu")

keg_Up_D18_S10_S30 <- kegga(Up_D18_S10_S30$GeneID, species.KEGG ="sdu")
keg_Down_D18_S10_S30 <- kegga(Down_D18_S10_S30$GeneID, species.KEGG ="sdu")
keg_Up_D18_S10_S20 <- kegga(Up_D18_S10_S20$GeneID, species.KEGG ="sdu")
keg_Down_D18_S10_S20 <- kegga(Down_D18_S10_S20$GeneID, species.KEGG ="sdu")
keg_Up_D18_S20_S30 <- kegga(Up_D18_S20_S30$GeneID, species.KEGG ="sdu")
keg_Down_D18_S20_S30 <- kegga(Down_D18_S20_S30$GeneID, species.KEGG ="sdu")

keg_Up_D20_S10_S30 <- kegga(Up_D20_S10_S30$GeneID, species.KEGG ="sdu")
keg_Down_D20_S10_S30 <- kegga(Down_D20_S10_S30$GeneID, species.KEGG ="sdu")
keg_Up_D20_S10_S20 <- kegga(Up_D20_S10_S20$GeneID, species.KEGG ="sdu")
keg_Down_D20_S10_S20 <- kegga(Down_D20_S10_S20$GeneID, species.KEGG ="sdu")
keg_Up_D20_S20_S30 <- kegga(Up_D20_S20_S30$GeneID, species.KEGG ="sdu")
keg_Down_D20_S20_S30 <- kegga(Down_D20_S20_S30$GeneID, species.KEGG ="sdu")

keg_Up_D24_S10_S30 <- kegga(Up_D24_S10_S30$GeneID, species.KEGG ="sdu")
keg_Down_D24_S10_S30 <- kegga(Down_D24_S10_S30$GeneID, species.KEGG ="sdu")
keg_Up_D24_S10_S20 <- kegga(Up_D24_S10_S20$GeneID, species.KEGG ="sdu")
keg_Down_D24_S10_S20 <- kegga(Down_D24_S10_S20$GeneID, species.KEGG ="sdu")
keg_Up_D24_S20_S30 <- kegga(Up_D24_S20_S30$GeneID, species.KEGG ="sdu")
keg_Down_D24_S20_S30 <- kegga(Down_D24_S20_S30$GeneID, species.KEGG ="sdu")

#Adjust the p values for each of these analysis for multiple testing
keg_Up_D3_S10_S30$padj <- p.adjust(keg_Up_D3_S10_S30$P.DE, method = "BH")
keg_Down_D3_S10_S30$padj <- p.adjust(keg_Down_D3_S10_S30$P.DE, method = "BH")
keg_Up_D3_S10_S20$padj <- p.adjust(keg_Up_D3_S10_S20$P.DE, method = "BH")
keg_Down_D3_S10_S20$padj <- p.adjust(keg_Down_D3_S10_S20$P.DE, method = "BH")
keg_Up_D3_S20_S30$padj <- p.adjust(keg_Up_D3_S20_S30$P.DE, method = "BH")
keg_Down_D3_S20_S30$padj <- p.adjust(keg_Down_D3_S20_S30$P.DE, method = "BH")

keg_Up_D6_S10_S30$padj <- p.adjust(keg_Up_D6_S10_S30$P.DE, method = "BH")
keg_Down_D6_S10_S30$padj <- p.adjust(keg_Down_D6_S10_S30$P.DE, method = "BH")
keg_Up_D6_S10_S20$padj <- p.adjust(keg_Up_D6_S10_S20$P.DE, method = "BH")
keg_Down_D6_S10_S20$padj <- p.adjust(keg_Down_D6_S10_S20$P.DE, method = "BH")
keg_Up_D6_S20_S30$padj <- p.adjust(keg_Up_D6_S20_S30$P.DE, method = "BH")
keg_Down_D6_S20_S30$padj <- p.adjust(keg_Down_D6_S20_S30$P.DE, method = "BH")

keg_Up_D12_S10_S30$padj <- p.adjust(keg_Up_D12_S10_S30$P.DE, method = "BH")
keg_Down_D12_S10_S30$padj <- p.adjust(keg_Down_D12_S10_S30$P.DE, method = "BH")
keg_Up_D12_S10_S20$padj <- p.adjust(keg_Up_D12_S10_S20$P.DE, method = "BH")
keg_Down_D12_S10_S20$padj <- p.adjust(keg_Down_D12_S10_S20$P.DE, method = "BH")
keg_Up_D12_S20_S30$padj <- p.adjust(keg_Up_D12_S20_S30$P.DE, method = "BH")
keg_Down_D12_S20_S30$padj <- p.adjust(keg_Down_D12_S20_S30$P.DE, method = "BH")

keg_Up_D15_S10_S30$padj <- p.adjust(keg_Up_D15_S10_S30$P.DE, method = "BH")
keg_Down_D15_S10_S30$padj <- p.adjust(keg_Down_D15_S10_S30$P.DE, method = "BH")
keg_Up_D15_S10_S20$padj <- p.adjust(keg_Up_D15_S10_S20$P.DE, method = "BH")
keg_Down_D15_S10_S20$padj <- p.adjust(keg_Down_D15_S10_S20$P.DE, method = "BH")
keg_Up_D15_S20_S30$padj <- p.adjust(keg_Up_D15_S20_S30$P.DE, method = "BH")
keg_Down_D15_S20_S30$padj <- p.adjust(keg_Down_D15_S20_S30$P.DE, method = "BH")

keg_Up_D18_S10_S30$padj <- p.adjust(keg_Up_D18_S10_S30$P.DE, method = "BH")
keg_Down_D18_S10_S30$padj <- p.adjust(keg_Down_D18_S10_S30$P.DE, method = "BH")
keg_Up_D18_S10_S20$padj <- p.adjust(keg_Up_D18_S10_S20$P.DE, method = "BH")
keg_Down_D18_S10_S20$padj <- p.adjust(keg_Down_D18_S10_S20$P.DE, method = "BH")
keg_Up_D18_S20_S30$padj <- p.adjust(keg_Up_D18_S20_S30$P.DE, method = "BH")
keg_Down_D18_S20_S30$padj <- p.adjust(keg_Down_D18_S20_S30$P.DE, method = "BH")

keg_Up_D20_S10_S30$padj <- p.adjust(keg_Up_D20_S10_S30$P.DE, method = "BH")
keg_Down_D20_S10_S30$padj <- p.adjust(keg_Down_D20_S10_S30$P.DE, method = "BH")
keg_Up_D20_S10_S20$padj <- p.adjust(keg_Up_D20_S10_S20$P.DE, method = "BH")
keg_Down_D20_S10_S20$padj <- p.adjust(keg_Down_D20_S10_S20$P.DE, method = "BH")
keg_Up_D20_S20_S30$padj <- p.adjust(keg_Up_D20_S20_S30$P.DE, method = "BH")
keg_Down_D20_S20_S30$padj <- p.adjust(keg_Down_D20_S20_S30$P.DE, method = "BH")

keg_Up_D24_S10_S30$padj <- p.adjust(keg_Up_D24_S10_S30$P.DE, method = "BH")
keg_Down_D24_S10_S30$padj <- p.adjust(keg_Down_D24_S10_S30$P.DE, method = "BH")
keg_Up_D24_S10_S20$padj <- p.adjust(keg_Up_D24_S10_S20$P.DE, method = "BH")
keg_Down_D24_S10_S20$padj <- p.adjust(keg_Down_D24_S10_S20$P.DE, method = "BH")
keg_Up_D24_S20_S30$padj <- p.adjust(keg_Up_D24_S20_S30$P.DE, method = "BH")
keg_Down_D24_S20_S30$padj <- p.adjust(keg_Down_D24_S20_S30$P.DE, method = "BH")

#Make dot plots to summarize each of these datasets
ggplot(keg_Up_D3_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/sum(DE)), color=padj, size=DE),
             data = head(keg_Up_D3_S10_S30[which(keg_Up_D3_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Down_D3_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(keg_Down_D3_S10_S30[which(keg_Down_D3_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))


ggplot(keg_Up_D6_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(keg_Up_D6_S10_S30[which(keg_Up_D6_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Down_D6_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(keg_Down_D6_S10_S30[which(keg_Down_D6_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Up_D12_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(keg_Up_D12_S10_S30[which(keg_Up_D12_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Down_D12_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(keg_Down_D12_S10_S30[which(keg_Down_D12_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Up_D15_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(keg_Up_D15_S10_S30[which(keg_Up_D15_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Down_D15_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(keg_Down_D15_S10_S30[which(keg_Down_D15_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Up_D18_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(keg_Up_D18_S10_S30[which(keg_Up_D18_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Down_D18_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(keg_Down_D18_S10_S30[which(keg_Down_D18_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Up_D20_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(keg_Up_D20_S10_S30[which(keg_Up_D20_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Down_D20_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(keg_Down_D20_S10_S30[which(keg_Down_D20_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Up_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(keg_Up_D24_S10_S30[which(keg_Up_D24_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(keg_Down_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(keg_Down_D24_S10_S30[which(keg_Down_D24_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "KEGG pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

###GO ENRICHMENT ANALYSIS WITH LIMMA LIBRARY###

#Conduct GO enrichment analysis using each set of these objects, kegga can do GO enrichment if the appropriate dataframes are provided. 
go_Up_D3_S10_S30 <- kegga(Up_D3_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D3_S10_S30 <- kegga(Down_D3_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D3_S10_S20 <- kegga(Up_D3_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D3_S10_S20 <- kegga(Down_D3_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D3_S20_S30 <- kegga(Up_D3_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D3_S20_S30 <- kegga(Down_D3_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)

go_Up_D6_S10_S30 <- kegga(Up_D6_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D6_S10_S30 <- kegga(Down_D6_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D6_S10_S20 <- kegga(Up_D6_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D6_S10_S20 <- kegga(Down_D6_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D6_S20_S30 <- kegga(Up_D6_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D6_S20_S30 <- kegga(Down_D6_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)

go_Up_D12_S10_S30 <- kegga(Up_D12_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D12_S10_S30 <- kegga(Down_D12_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D12_S10_S20 <- kegga(Up_D12_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D12_S10_S20 <- kegga(Down_D12_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D12_S20_S30 <- kegga(Up_D12_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D12_S20_S30 <- kegga(Down_D12_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)

go_Up_D15_S10_S30 <- kegga(Up_D15_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D15_S10_S30 <- kegga(Down_D15_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D15_S10_S20 <- kegga(Up_D15_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D15_S10_S20 <- kegga(Down_D15_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D15_S20_S30 <- kegga(Up_D15_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D15_S20_S30 <- kegga(Down_D15_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)

go_Up_D18_S10_S30 <- kegga(Up_D18_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D18_S10_S30 <- kegga(Down_D18_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D18_S10_S20 <- kegga(Up_D18_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D18_S10_S20 <- kegga(Down_D18_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D18_S20_S30 <- kegga(Up_D18_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D18_S20_S30 <- kegga(Down_D18_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)

go_Up_D20_S10_S30 <- kegga(Up_D20_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D20_S10_S30 <- kegga(Down_D20_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D20_S10_S20 <- kegga(Up_D20_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D20_S10_S20 <- kegga(Down_D20_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D20_S20_S30 <- kegga(Up_D20_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D20_S20_S30 <- kegga(Down_D20_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)

go_Up_D24_S10_S30 <- kegga(Up_D24_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D24_S10_S30 <- kegga(Down_D24_S10_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D24_S10_S20 <- kegga(Up_D24_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D24_S10_S20 <- kegga(Down_D24_S10_S20$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Up_D24_S20_S30 <- kegga(Up_D24_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)
go_Down_D24_S20_S30 <- kegga(Down_D24_S20_S30$GeneID, gene.pathway=GENEID_GO_Terms, pathway.names=GO.Name)

#Adjust the p values for each of these analysis for multiple testing
go_Up_D3_S10_S30$padj <- p.adjust(go_Up_D3_S10_S30$P.DE, method = "BH")
go_Down_D3_S10_S30$padj <- p.adjust(go_Down_D3_S10_S30$P.DE, method = "BH")
go_Up_D3_S10_S20$padj <- p.adjust(go_Up_D3_S10_S20$P.DE, method = "BH")
go_Down_D3_S10_S20$padj <- p.adjust(go_Down_D3_S10_S20$P.DE, method = "BH")
go_Up_D3_S20_S30$padj <- p.adjust(go_Up_D3_S20_S30$P.DE, method = "BH")
go_Down_D3_S20_S30$padj <- p.adjust(go_Down_D3_S20_S30$P.DE, method = "BH")

go_Up_D6_S10_S30$padj <- p.adjust(go_Up_D6_S10_S30$P.DE, method = "BH")
go_Down_D6_S10_S30$padj <- p.adjust(go_Down_D6_S10_S30$P.DE, method = "BH")
go_Up_D6_S10_S20$padj <- p.adjust(go_Up_D6_S10_S20$P.DE, method = "BH")
go_Down_D6_S10_S20$padj <- p.adjust(go_Down_D6_S10_S20$P.DE, method = "BH")
go_Up_D6_S20_S30$padj <- p.adjust(go_Up_D6_S20_S30$P.DE, method = "BH")
go_Down_D6_S20_S30$padj <- p.adjust(go_Down_D6_S20_S30$P.DE, method = "BH")

go_Up_D12_S10_S30$padj <- p.adjust(go_Up_D12_S10_S30$P.DE, method = "BH")
go_Down_D12_S10_S30$padj <- p.adjust(go_Down_D12_S10_S30$P.DE, method = "BH")
go_Up_D12_S10_S20$padj <- p.adjust(go_Up_D12_S10_S20$P.DE, method = "BH")
go_Down_D12_S10_S20$padj <- p.adjust(go_Down_D12_S10_S20$P.DE, method = "BH")
go_Up_D12_S20_S30$padj <- p.adjust(go_Up_D12_S20_S30$P.DE, method = "BH")
go_Down_D12_S20_S30$padj <- p.adjust(go_Down_D12_S20_S30$P.DE, method = "BH")

go_Up_D15_S10_S30$padj <- p.adjust(go_Up_D15_S10_S30$P.DE, method = "BH")
go_Down_D15_S10_S30$padj <- p.adjust(go_Down_D15_S10_S30$P.DE, method = "BH")
go_Up_D15_S10_S20$padj <- p.adjust(go_Up_D15_S10_S20$P.DE, method = "BH")
go_Down_D15_S10_S20$padj <- p.adjust(go_Down_D15_S10_S20$P.DE, method = "BH")
go_Up_D15_S20_S30$padj <- p.adjust(go_Up_D15_S20_S30$P.DE, method = "BH")
go_Down_D15_S20_S30$padj <- p.adjust(go_Down_D15_S20_S30$P.DE, method = "BH")

go_Up_D18_S10_S30$padj <- p.adjust(go_Up_D18_S10_S30$P.DE, method = "BH")
go_Down_D18_S10_S30$padj <- p.adjust(go_Down_D18_S10_S30$P.DE, method = "BH")
go_Up_D18_S10_S20$padj <- p.adjust(go_Up_D18_S10_S20$P.DE, method = "BH")
go_Down_D18_S10_S20$padj <- p.adjust(go_Down_D18_S10_S20$P.DE, method = "BH")
go_Up_D18_S20_S30$padj <- p.adjust(go_Up_D18_S20_S30$P.DE, method = "BH")
go_Down_D18_S20_S30$padj <- p.adjust(go_Down_D18_S20_S30$P.DE, method = "BH")

go_Up_D20_S10_S30$padj <- p.adjust(go_Up_D20_S10_S30$P.DE, method = "BH")
go_Down_D20_S10_S30$padj <- p.adjust(go_Down_D20_S10_S30$P.DE, method = "BH")
go_Up_D20_S10_S20$padj <- p.adjust(go_Up_D20_S10_S20$P.DE, method = "BH")
go_Down_D20_S10_S20$padj <- p.adjust(go_Down_D20_S10_S20$P.DE, method = "BH")
go_Up_D20_S20_S30$padj <- p.adjust(go_Up_D20_S20_S30$P.DE, method = "BH")
go_Down_D20_S20_S30$padj <- p.adjust(go_Down_D20_S20_S30$P.DE, method = "BH")

go_Up_D24_S10_S30$padj <- p.adjust(go_Up_D24_S10_S30$P.DE, method = "BH")
go_Down_D24_S10_S30$padj <- p.adjust(go_Down_D24_S10_S30$P.DE, method = "BH")
go_Up_D24_S10_S20$padj <- p.adjust(go_Up_D24_S10_S20$P.DE, method = "BH")
go_Down_D24_S10_S20$padj <- p.adjust(go_Down_D24_S10_S20$P.DE, method = "BH")
go_Up_D24_S20_S30$padj <- p.adjust(go_Up_D24_S20_S30$P.DE, method = "BH")
go_Down_D24_S20_S30$padj <- p.adjust(go_Down_D24_S20_S30$P.DE, method = "BH")

#Make dot plots to summarize each of these datasets
ggplot(go_Up_D3_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(go_Up_D3_S10_S30[which(go_Up_D3_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Down_D3_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(go_Down_D3_S10_S30[which(go_Down_D3_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))


ggplot(go_Up_D6_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(go_Up_D6_S10_S30[which(go_Up_D6_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Down_D6_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(go_Down_D6_S10_S30[which(go_Down_D6_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Up_D12_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(go_Up_D12_S10_S30[which(go_Up_D12_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Down_D12_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(go_Down_D12_S10_S30[which(go_Down_D12_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Up_D15_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(go_Up_D15_S10_S30[which(go_Up_D15_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Down_D15_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(go_Down_D15_S10_S30[which(go_Down_D15_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Up_D18_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(go_Up_D18_S10_S30[which(go_Up_D18_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Down_D18_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(go_Down_D18_S10_S30[which(go_Down_D18_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Up_D20_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(go_Up_D20_S10_S30[which(go_Up_D20_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Down_D20_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(go_Down_D20_S10_S30[which(go_Down_D20_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Up_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=padj, size=DE),
             data = head(go_Up_D24_S10_S30[which(go_Up_D24_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

ggplot(go_Down_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Pathway, -padj), x=(DE/N), color=-padj, size=DE),
             data = head(go_Down_D24_S10_S30[which(go_Down_D24_S10_S30$padj < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Names", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(legend.key.size = unit(1.2, "line"))

###GO AND KEEG ENRICHMENT ANALYSIS WITH CLUSTERPROFILER

#Switch around the Gene ID and GO Terms
GO_TERM_2_GENE_ID <- GENEID_GO_Terms[c(2,1)]

#Create a list of all these dataframes
x <- list(Up_D3_S10_S20, Down_D3_S10_S20, Up_D3_S10_S30, Down_D3_S10_S30, Up_D3_S20_S30, Down_D3_S20_S30,Up_D6_S10_S20, Down_D6_S10_S20, Up_D6_S10_S30, Down_D6_S10_S30, Up_D6_S20_S30, Down_D6_S20_S30,Up_D12_S10_S20, Down_D12_S10_S20, Up_D12_S10_S30, Down_D12_S10_S30, Up_D12_S20_S30, Down_D12_S20_S30,Up_D15_S10_S20, Down_D15_S10_S20, Up_D15_S10_S30, Down_D15_S10_S30, Up_D15_S20_S30, Down_D15_S20_S30,Up_D18_S10_S20, Down_D18_S10_S20, Up_D18_S10_S30, Down_D18_S10_S30, Up_D18_S20_S30, Down_D18_S20_S30, Up_D20_S10_S20, Down_D20_S10_S20, Up_D20_S10_S30, Down_D20_S10_S30, Up_D20_S20_S30, Down_D20_S20_S30,Up_D24_S10_S20, Down_D24_S10_S20, Up_D24_S10_S30, Down_D24_S10_S30, Up_D24_S20_S30, Down_D24_S20_S30)

#Add names to each one
names(x) <- list("Up_D3_S10_S20", "Down_D3_S10_S20", "Up_D3_S10_S30", "Down_D3_S10_S30", "Up_D3_S20_S30", "Down_D3_S20_S30", "Up_D6_S10_S20", "Down_D6_S10_S20", "Up_D6_S10_S30", "Down_D6_S10_S30", "Up_D6_S20_S30", "Down_D6_S20_S30",    "Up_D12_S10_S20", "Down_D12_S10_S20", "Up_D12_S10_S30", "Down_D12_S10_S30", "Up_D12_S20_S30", "Down_D12_S20_S30","Up_D15_S10_S20", "Down_D15_S10_S20", "Up_D15_S10_S30", "Down_D15_S10_S30", "Up_D15_S20_S30", "Down_D15_S20_S30",
                 "Up_D18_S10_S20", "Down_D18_S10_S20", "Up_D18_S10_S30", "Down_D18_S10_S30", "Up_D18_S20_S30", "Down_D18_S20_S30", "Up_D20_S10_S20", "Down_D20_S10_S20", "Up_D20_S10_S30", "Down_D20_S10_S30", "Up_D20_S20_S30", "Down_D20_S20_S30",
                 "Up_D24_S10_S20", "Down_D24_S10_S20", "Up_D24_S10_S30", "Down_D24_S10_S30", "Up_D24_S20_S30", "Down_D24_S20_S30")

#Determine if there are any overlapping genes across all dph at each up or down regulated set of DEGS for each salinity

#S10 vs S20 Common Genes
Up_S10_S20 <- semi_join(Up_D3_S10_S20, Up_D6_S10_S20, by = "Name")
Up_S10_S20 <- semi_join(Up_S10_S20, Up_D12_S10_S20, by = "Name")
Up_S10_S20 <- semi_join(Up_S10_S20, Up_D15_S10_S20, by = "Name")
Up_S10_S20 <- semi_join(Up_S10_S20, Up_D18_S10_S20, by = "Name")
Up_S10_S20 <- semi_join(Up_S10_S20, Up_D20_S10_S20, by = "Name")
Up_S10_S20 <- semi_join(Up_S10_S20, Up_D24_S10_S20, by = "Name")
nrow(Up_S10_S20) #0


Down_S10_S20 <- semi_join(Down_D3_S10_S20, Down_D6_S10_S20, by = "Name")
Down_S10_S20 <- semi_join(Down_S10_S20, Down_D12_S10_S20, by = "Name")
Down_S10_S20 <- semi_join(Down_S10_S20, Down_D15_S10_S20, by = "Name")
Down_S10_S20 <- semi_join(Down_S10_S20, Down_D18_S10_S20, by = "Name")
Down_S10_S20 <- semi_join(Down_S10_S20, Down_D20_S10_S20, by = "Name")
Down_S10_S20 <- semi_join(Down_S10_S20, Down_D24_S10_S20, by = "Name")
nrow(Down_S10_S20) #0


#S10 vs S30 Common Genes
Up_S10_S30 <- semi_join(Up_D3_S10_S30, Up_D6_S10_S30, by = "Name")
Up_S10_S30 <- semi_join(Up_S10_S30, Up_D12_S10_S30, by = "Name")
Up_S10_S30 <- semi_join(Up_S10_S30, Up_D15_S10_S30, by = "Name")
Up_S10_S30 <- semi_join(Up_S10_S30, Up_D18_S10_S30, by = "Name")
Up_S10_S30 <- semi_join(Up_S10_S30, Up_D20_S10_S30, by = "Name")
Up_S10_S30 <- semi_join(Up_S10_S30, Up_D24_S10_S30, by = "Name")
nrow(Up_S10_S30) #1
#LOC111231300 111231300 uncharacterized

Down_S10_S30 <- semi_join(Down_D3_S10_S30, Down_D6_S10_S30, by = "Name")
Down_S10_S30 <- semi_join(Down_S10_S30, Down_D12_S10_S30, by = "Name")
Down_S10_S30 <- semi_join(Down_S10_S30, Down_D15_S10_S30, by = "Name")
Down_S10_S30 <- semi_join(Down_S10_S30, Down_D18_S10_S30, by = "Name")
Down_S10_S30 <- semi_join(Down_S10_S30, Down_D20_S10_S30, by = "Name")
Down_S10_S30 <- semi_join(Down_S10_S30, Down_D24_S10_S30, by = "Name")
nrow(Down_S10_S30) #2
#LOC111226922 111226922 uncharacterized
#LOC111239098 111239098 transmembrane protease serine 2-like

#S20 vs S30 Common Genes
Up_S20_S30 <- semi_join(Up_D3_S20_S30, Up_D6_S20_S30, by = "Name")
Up_S20_S30 <- semi_join(Up_S20_S30, Up_D12_S20_S30, by = "Name")
Up_S20_S30 <- semi_join(Up_S20_S30, Up_D15_S20_S30, by = "Name")
Up_S20_S30 <- semi_join(Up_S20_S30, Up_D18_S20_S30, by = "Name")
Up_S20_S30 <- semi_join(Up_S20_S30, Up_D20_S20_S30, by = "Name")
Up_S20_S30 <- semi_join(Up_S20_S30, Up_D24_S20_S30, by = "Name")
nrow(Up_S20_S30) #0


Down_S20_S30 <- semi_join(Down_D3_S20_S30, Down_D6_S20_S30, by = "Name")
Down_S20_S30 <- semi_join(Down_S20_S30, Down_D12_S20_S30, by = "Name")
Down_S20_S30 <- semi_join(Down_S20_S30, Down_D15_S20_S30, by = "Name")
Down_S20_S30 <- semi_join(Down_S20_S30, Down_D18_S20_S30, by = "Name")
Down_S20_S30 <- semi_join(Down_S20_S30, Down_D20_S20_S30, by = "Name")
Down_S20_S30 <- semi_join(Down_S20_S30, Down_D24_S20_S30, by = "Name")
nrow(Down_S20_S30) #0

#Create a function to prepare for a run clusterProfiler's enricher tool for GO Enrichment Analysis

#Three elements go into the function, the dataframe of CLC Genomics differential expression results with Log Fold Change and Gene Id information (df), a dataframe pairing GO Terms to each Gene Id (y, here a default dataframe is given), and a dataframe associating the GO TERM with a GO Name

enricher_prep <- function(df, y = GO_TERM_2_GENE_ID, z = GO.Name) {
  ## feature 1: numeric vector, column with Fold Change values
  geneList = df[,6]
  ## feature 2: add names to the vector based upon the Gene Id column
  names(geneList) = as.character(df[,10])
  #Keep values for which there is actually a fold change, to allow conversion to      numeric
  df <- df[which(df$Fold.change!="#N/A"),]
  #Make Fold change numeric
  df$Fold.change <- as.numeric(df$Fold.change)
  ## feature 3: Make the list in decreasing order
  geneList = sort(geneList, decreasing = TRUE)
  #Keep DEGS with FC > 2
  gene <- names(geneList)#[abs(geneList) > 2]
  #run enrichr and save the results as a data frame
  ewp <- enricher(gene, TERM2GENE = y, TERM2NAME = z)
  df_ewp <- as.data.frame(ewp)
  #Create another column using the already existing GeneRatio column to create a one value ratio instead of a fraction which is easier to plot
  df_ewp$GeneRatio2 <- sapply(strsplit(df_ewp$GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  #Return this dataframe
  return(df_ewp)
}

#Create a for loop that will save the results of the above function as a distinctly named dataframe
for (i in names(x)) {
  #Use the dataframe name to recall the dataframe and store it
  df <- x[[i]]
  #run the above function and save as a different dataframe
  df2 <- enricher_prep(df)
  #assign that dataframe a different name based upon the original dataframe
  assign(paste("enrichr_go", i, sep="_"), df2)
}

enrichKEGG_prep <- function(df) {
  ## feature 1: numeric vector, column with Fold Change values
  geneList = df[,6]
  ## feature 2: add names to the vector based upon the Gene Id column
  names(geneList) = as.character(df[,10])
  #Keep values for which there is actually a fold change, to allow conversion to      numeric
  df <- df[which(df$Fold.change!="#N/A"),]
  #Make Fold change numeric
  df$Fold.change <- as.numeric(df$Fold.change)
  ## feature 3: Make the list in decreasing order
  geneList = sort(geneList, decreasing = TRUE)
  #Keep DEGS with FC > 2
  gene <- names(geneList)#[abs(geneList) > 2]
  #run enrichKEGG and save the results as a data frame
  kk <- enrichKEGG(gene = gene, organism = 'sdu', pvalueCutoff = 0.05)
  df_kk <- as.data.frame(kk)
  #Create another column using the already existing GeneRatio column to create a one value ratio instead of a fraction which is easier to plot
  df_kk$GeneRatio2 <- sapply(strsplit(df_kk$GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  return(df_kk)
}

#Create a for loop that will save the results of the above function as a distinctly named dataframe
for (i in names(x)) {
  #Use the dataframe name to recall the dataframe and store it
  df <- x[[i]]
  #run the above function and save as a different dataframe
  df2 <- enrichKEGG_prep(df)
  #assign that dataframe a different name based upon the original dataframe
  assign(paste("enrichr_kegg", i, sep="_"), df2)
}

##Upload and do same thing for pooled samples
All_S10_S30 <- read.csv("All DPH 10 vs. 30 ppt.csv")
All_S10_S20 <- read.csv("All DPH 10 vs. 20 ppt.csv")
All_S20_S30 <- read.csv("All DPH 20 vs. 30 ppt.csv")

Up_All_S10_S30 <- subset(All_S10_S30, Log..fold.change >0)
Down_All_S10_S30 <- subset(All_S10_S30, Log..fold.change <0)
Up_All_S10_S20 <- subset(All_S10_S20, Log..fold.change >0)
Down_All_S10_S20 <- subset(All_S10_S20, Log..fold.change <0)
Up_All_S20_S30 <- subset(All_S20_S30, Log..fold.change >0)
Down_All_S20_S30 <- subset(All_S20_S30, Log..fold.change <0)


enrichr_go_Up_All_S10_S30 <- enricher_prep(Up_All_S10_S30)
enrichr_go_Down_All_S10_S30 <- enricher_prep(Down_All_S10_S30)
enrichr_go_Up_All_S10_S20 <- enricher_prep(Up_All_S10_S20)
enrichr_go_Down_All_S10_S20 <- enricher_prep(Down_All_S10_S20)
enrichr_go_Up_All_S20_S30 <- enricher_prep(Up_All_S20_S30)
enrichr_go_Down_All_S20_S30 <- enricher_prep(Down_All_S20_S30)

enrichr_kegg_Up_All_S10_S30 <- enrichKEGG_prep(Up_All_S10_S30)
enrichr_kegg_Down_All_S10_S30 <- enrichKEGG_prep(Down_All_S10_S30)
enrichr_kegg_Up_All_S10_S20 <- enrichKEGG_prep(Up_All_S10_S20)
enrichr_kegg_Down_All_S10_S20 <- enrichKEGG_prep(Down_All_S10_S20)
enrichr_kegg_Up_All_S20_S30 <- enrichKEGG_prep(Up_All_S20_S30)
enrichr_kegg_Down_All_S20_S30 <- enrichKEGG_prep(Down_All_S20_S30)

#Determine if there are overlapping up or downregulated pathways across all dph at each salinity

#S10 vs S20 Common Pathways
enrichr_kegg_Up_S10_S20 <- semi_join(enrichr_kegg_Up_D3_S10_S20, enrichr_kegg_Up_D6_S10_S20, by = "ID")
enrichr_kegg_Up_S10_S20 <- semi_join(enrichr_kegg_Up_S10_S20, enrichr_kegg_Up_D12_S10_S20, by = "ID")
enrichr_kegg_Up_S10_S20 <- semi_join(enrichr_kegg_Up_S10_S20, enrichr_kegg_Up_D15_S10_S20, by = "ID")
enrichr_kegg_Up_S10_S20 <- semi_join(enrichr_kegg_Up_S10_S20, enrichr_kegg_Up_D18_S10_S20, by = "ID")
enrichr_kegg_Up_S10_S20 <- semi_join(enrichr_kegg_Up_S10_S20, enrichr_kegg_Up_D20_S10_S20, by = "ID")
enrichr_kegg_Up_S10_S20 <- semi_join(enrichr_kegg_Up_S10_S20, enrichr_kegg_Up_D24_S10_S20, by = "ID")
nrow(enrichr_kegg_Up_S10_S20) #0


enrichr_kegg_Down_S10_S20 <- semi_join(enrichr_kegg_Down_D3_S10_S20, enrichr_kegg_Down_D6_S10_S20, by = "ID")
enrichr_kegg_Down_S10_S20 <- semi_join(enrichr_kegg_Down_S10_S20, enrichr_kegg_Down_D12_S10_S20, by = "ID")
enrichr_kegg_Down_S10_S20 <- semi_join(enrichr_kegg_Down_S10_S20, enrichr_kegg_Down_D15_S10_S20, by = "ID")
enrichr_kegg_Down_S10_S20 <- semi_join(enrichr_kegg_Down_S10_S20, enrichr_kegg_Down_D18_S10_S20, by = "ID")
enrichr_kegg_Down_S10_S20 <- semi_join(enrichr_kegg_Down_S10_S20, enrichr_kegg_Down_D20_S10_S20, by = "ID")
enrichr_kegg_Down_S10_S20 <- semi_join(enrichr_kegg_Down_S10_S20, enrichr_kegg_Down_D24_S10_S20, by = "ID")
nrow(enrichr_kegg_Down_S10_S20) #0


#S10 vs S30 Common Pathways
enrichr_kegg_Up_S10_S30 <- semi_join(enrichr_kegg_Up_D3_S10_S30, enrichr_kegg_Up_D6_S10_S30, by = "ID")
enrichr_kegg_Up_S10_S30 <- semi_join(enrichr_kegg_Up_S10_S30, enrichr_kegg_Up_D12_S10_S30, by = "ID")
enrichr_kegg_Up_S10_S30 <- semi_join(enrichr_kegg_Up_S10_S30, enrichr_kegg_Up_D15_S10_S30, by = "ID")
enrichr_kegg_Up_S10_S30 <- semi_join(enrichr_kegg_Up_S10_S30, enrichr_kegg_Up_D18_S10_S30, by = "ID")
enrichr_kegg_Up_S10_S30 <- semi_join(enrichr_kegg_Up_S10_S30, enrichr_kegg_Up_D20_S10_S30, by = "ID")
enrichr_kegg_Up_S10_S30 <- semi_join(enrichr_kegg_Up_S10_S30, enrichr_kegg_Up_D24_S10_S30, by = "ID")
nrow(enrichr_kegg_Up_S10_S30) #1
#Cytokine-cytokine receptor interaction

enrichr_kegg_Down_S10_S30 <- semi_join(enrichr_kegg_Down_D3_S10_S30, enrichr_kegg_Down_D6_S10_S30, by = "ID")
enrichr_kegg_Down_S10_S30 <- semi_join(enrichr_kegg_Down_S10_S30, enrichr_kegg_Down_D12_S10_S30, by = "ID")
enrichr_kegg_Down_S10_S30 <- semi_join(enrichr_kegg_Down_S10_S30, enrichr_kegg_Down_D15_S10_S30, by = "ID")
enrichr_kegg_Down_S10_S30 <- semi_join(enrichr_kegg_Down_S10_S30, enrichr_kegg_Down_D18_S10_S30, by = "ID")
enrichr_kegg_Down_S10_S30 <- semi_join(enrichr_kegg_Down_S10_S30, enrichr_kegg_Down_D20_S10_S30, by = "ID")
enrichr_kegg_Down_S10_S30 <- semi_join(enrichr_kegg_Down_S10_S30, enrichr_kegg_Down_D24_S10_S30, by = "ID")
nrow(enrichr_kegg_Down_S10_S30) #0

#S20 vs S30 Common Pathways
enrichr_kegg_Up_S20_S30 <- semi_join(enrichr_kegg_Up_D3_S20_S30, enrichr_kegg_Up_D6_S20_S30, by = "ID")
enrichr_kegg_Up_S20_S30 <- semi_join(enrichr_kegg_Up_S20_S30, enrichr_kegg_Up_D12_S20_S30, by = "ID")
enrichr_kegg_Up_S20_S30 <- semi_join(enrichr_kegg_Up_S20_S30, enrichr_kegg_Up_D15_S20_S30, by = "ID")
enrichr_kegg_Up_S20_S30 <- semi_join(enrichr_kegg_Up_S20_S30, enrichr_kegg_Up_D18_S20_S30, by = "ID")
enrichr_kegg_Up_S20_S30 <- semi_join(enrichr_kegg_Up_S20_S30, enrichr_kegg_Up_D20_S20_S30, by = "ID")
enrichr_kegg_Up_S20_S30 <- semi_join(enrichr_kegg_Up_S20_S30, enrichr_kegg_Up_D24_S20_S30, by = "ID")
nrow(enrichr_kegg_Up_S20_S30) #0


enrichr_kegg_Down_S20_S30 <- semi_join(enrichr_kegg_Down_D3_S20_S30, enrichr_kegg_Down_D6_S20_S30, by = "ID")
enrichr_kegg_Down_S20_S30 <- semi_join(enrichr_kegg_Down_S20_S30, enrichr_kegg_Down_D12_S20_S30, by = "ID")
enrichr_kegg_Down_S20_S30 <- semi_join(enrichr_kegg_Down_S20_S30, enrichr_kegg_Down_D15_S20_S30, by = "ID")
enrichr_kegg_Down_S20_S30 <- semi_join(enrichr_kegg_Down_S20_S30, enrichr_kegg_Down_D18_S20_S30, by = "ID")
enrichr_kegg_Down_S20_S30 <- semi_join(enrichr_kegg_Down_S20_S30, enrichr_kegg_Down_D20_S20_S30, by = "ID")
enrichr_kegg_Down_S20_S30 <- semi_join(enrichr_kegg_Down_S20_S30, enrichr_kegg_Down_D24_S20_S30, by = "ID")
nrow(enrichr_kegg_Down_S20_S30) #0

#Test out two ways of visualizing the results

##Version 1

#Plot the S10 upregulated GO Terms
D24_Up_S10vsS30_GO_point_plot <- ggplot(enrichr_go_Up_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Description, -p.adjust), x=(GeneRatio2), color=p.adjust, size=Count),
             data = head(enrichr_go_Up_D24_S10_S30[which(enrichr_go_Up_D24_S10_S30$p.adjust < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Terms", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('Upregulated') +
  theme(legend.key.size = unit(1.2, "line"))
D24_Up_S10vsS30_GO_point_plot

#Plot the downregulated GO Terms
D24_Down_S10vsS30_GO_point_plot <- ggplot(enrichr_go_Down_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Description, -p.adjust), x=(GeneRatio2), color=-p.adjust, size=Count),
             data = head(enrichr_go_Down_D24_S10_S30[which(enrichr_go_Down_D24_S10_S30$p.adjust < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "GO Terms", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('Downregulated') +
  theme(legend.key.size = unit(1.2, "line"))
D24_Down_S10vsS30_GO_point_plot

#Combine the up- and downregulated GO Term graphs
GO_figure <- ggarrange(D24_Up_S10vsS30_GO_point_plot, D24_Down_S10vsS30_GO_point_plot,
                       labels = c("A", "B"),
                       ncol = 1, nrow = 2)
GO_figure

tiff('24 DPH Enrich GO Terms between 10 ppt and 20 ppt.tiff', units="in", width=6, height=8, res=300)
GO_figure
dev.off()

#Do the same thing for KEGG pathways
D24_Up_S10vsS30_KEGG_point_plot <- ggplot(enrichr_kegg_Up_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Description, -p.adjust), x=(GeneRatio2), color=p.adjust, size=Count),
             data = head(enrichr_kegg_Up_D24_S10_S30[which(enrichr_kegg_Up_D24_S10_S30$p.adjust < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "Enriched KEGG Pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('Upregulated') +
  theme(legend.key.size = unit(1.2, "line"))
D24_Up_S10vsS30_KEGG_point_plot

D24_Down_S10vsS30_KEGG_point_plot <- ggplot(enrichr_kegg_Down_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Description, -p.adjust), x=(GeneRatio2), color=-p.adjust, size=Count),
             data = head(enrichr_kegg_Down_D24_S10_S30[which(enrichr_kegg_Down_D24_S10_S30$p.adjust < 0.05),], n = 30))+
  scale_color_continuous(low="red", high="blue") +
  labs(y = "Enriched KEGG Pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('Downregulated') +
  theme(legend.key.size = unit(1.2, "line"))
D24_Down_S10vsS30_KEGG_point_plot

KEGG_figure <- ggarrange(D24_Up_S10vsS30_KEGG_point_plot, D24_Down_S10vsS30_KEGG_point_plot,
                         labels = c("A", "B"),
                         ncol = 1, nrow = 2)
KEGG_figure

tiff('24 DPH Enrich KEGG Terms between 10 ppt and 20 ppt.tiff', units="in", width=6, height=8, res=300)
KEGG_figure
dev.off()

ggarrange(
  lp,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(bxp, dp, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A"       # Label of the line plot
) 



##Version 2

#Add a new column to each of the up and downregualted dataframes denoting direction and combine them into one dataframe
enrichr_kegg_Down_D24_S10_S30$Direction <- "Down"

enrichr_kegg_Up_D24_S10_S30$Direction <- "Up"

enrichr_kegg_D24_S10_S30 <- rbind(enrichr_kegg_Down_D24_S10_S30, enrichr_kegg_Up_D24_S10_S30)

#Plot the dataframe with direction dictated by shape
D24_S10vsS30_KEGG_point_plot <- ggplot(enrichr_kegg_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Description, -p.adjust), x=(GeneRatio2), color = p.adjust, size=Count, shape=Direction, fill = p.adjust),
             data = head(enrichr_kegg_D24_S10_S30[which(enrichr_kegg_D24_S10_S30$p.adjust < 0.05),], n = 37))+
  scale_shape_manual(values=c(25,24))+
  scale_fill_continuous(low="red", high="blue") +
  scale_color_continuous(low="red", high="blue", guide="none") +
  labs(y = "Enriched KEGG Pathways", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('24 DPH 10ppt vs 30 ppt') +
  theme(legend.key.size = unit(1.2, "line"))
D24_S10vsS30_KEGG_point_plot

tiff('24 DPH Enriched KEGG Terms between 10 ppt and 30 ppt V2.tiff', units="in", width=6, height=6, res=300)
D24_S10vsS30_KEGG_point_plot
dev.off()

#Do same thing for GO Terms
enrichr_go_Down_D24_S10_S30$Direction <- "Down"

enrichr_go_Up_D24_S10_S30$Direction <- "Up"

enrichr_go_D24_S10_S30 <- rbind(enrichr_go_Down_D24_S10_S30, enrichr_go_Up_D24_S10_S30)

enrichr_go_D24_S10_S30["Description"][enrichr_go_D24_S10_S30["Description"] == "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"] <- "oxidoreductase acitvity version 2"


enrichr_go_D24_S10_S30 <- enrichr_go_D24_S10_S30[order(enrichr_go_D24_S10_S30$p.adjust),]

D24_S10vsS30_GO_point_plot <- ggplot(enrichr_go_D24_S10_S30)+
  geom_point(mapping = aes(y= reorder(Description, -p.adjust), x=(GeneRatio2), color = p.adjust, size=Count, shape=Direction, fill = p.adjust),
             data = head(enrichr_go_D24_S10_S30[which(enrichr_go_D24_S10_S30$p.adjust < 0.05),], n = 37))+
  scale_shape_manual(values=c(25,24))+
  scale_fill_continuous(low="red", high="blue") +
  scale_color_continuous(low="red", high="blue", guide="none") +
  labs(y = "Enriched GO Terms", x = "Gene Ratio", color = "p value", size = "Counts") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('24 DPH 10ppt vs 30 ppt') +
  theme(legend.key.size = unit(1.2, "line"))
D24_S10vsS30_GO_point_plot

tiff('24 DPH Enrich GO Terms between 10 ppt and 30 ppt V2.tiff', units="in", width=6, height=6, res=300)
D24_S10vsS30_GO_point_plot
dev.off()


###GSEA ANALYSIS AT THE DPH LEVEL

#Upload the total genes for each dataset separated by dph
All_D3_S10_S30 <- read.csv("All 3 DPH 10 vs. 30 ppt.csv")
All_D3_S10_S20 <- read.csv("All 3 DPH 10 vs. 20 ppt.csv")
All_D3_S20_S30 <- read.csv("All 3 DPH 20 vs. 30 ppt.csv")
All_D6_S10_S30 <- read.csv("All 6 DPH 10 vs. 30 ppt.csv")
All_D6_S10_S20 <- read.csv("All 6 DPH 10 vs. 20 ppt.csv")
All_D6_S20_S30 <- read.csv("All 6 DPH 20 vs. 30 ppt.csv")
All_D12_S10_S30 <- read.csv("All 12 DPH 10 vs. 30 ppt.csv")
All_D12_S10_S20 <- read.csv("All 12 DPH 10 vs. 20 ppt.csv")
All_D12_S20_S30 <- read.csv("All 12 DPH 20 vs. 30 ppt.csv")
All_D15_S10_S30 <- read.csv("All 15 DPH 10 vs. 30 ppt.csv")
All_D15_S10_S20 <- read.csv("All 15 DPH 10 vs. 20 ppt.csv")
All_D15_S20_S30 <- read.csv("All 15 DPH 20 vs. 30 ppt.csv")
All_D18_S10_S30 <- read.csv("All 18 DPH 10 vs. 30 ppt.csv")
All_D18_S10_S20 <- read.csv("All 18 DPH 10 vs. 20 ppt.csv")
All_D18_S20_S30 <- read.csv("All 18 DPH 20 vs. 30 ppt.csv")
All_D20_S10_S30 <- read.csv("All 20 DPH 10 vs. 30 ppt.csv")
All_D20_S10_S20 <- read.csv("All 20 DPH 10 vs. 20 ppt.csv")
All_D20_S20_S30 <- read.csv("All 20 DPH 20 vs. 30 ppt.csv")
All_D24_S10_S30 <- read.csv("All 24 DPH 10 vs. 30 ppt.csv")
All_D24_S10_S20 <- read.csv("All 24 DPH 10 vs. 20 ppt.csv")
All_D24_S20_S30 <- read.csv("All 24 DPH 20 vs. 30 ppt.csv")

#Create a list of the above dataframes and add names to them
y <- list(All_D3_S10_S20, All_D3_S10_S30, All_D3_S20_S30, All_D6_S10_S20, All_D6_S10_S30, All_D6_S20_S30, All_D12_S10_S20, All_D12_S10_S30, All_D12_S20_S30, All_D15_S10_S20, All_D15_S10_S30, All_D15_S20_S30, All_D18_S10_S20,  All_D18_S10_S30, All_D18_S20_S30, All_D20_S10_S20,  All_D20_S10_S30, All_D20_S20_S30, All_D24_S10_S20,  All_D24_S10_S30, All_D24_S20_S30)

names(y) <- list("D3_S10_S20", "D3_S10_S30", "D3_S20_S30", "D6_S10_S20", "D6_S10_S30", "D6_S20_S30",    "D12_S10_S20", "D12_S10_S30", "D12_S20_S30", "D15_S10_S20", "D15_S10_S30",  "D15_S20_S30", "D18_S10_S20", "D18_S10_S30", "D18_S20_S30",  "D20_S10_S20", "D20_S10_S30", "D20_S20_S30", "D24_S10_S20",  "D24_S10_S30", "D24_S20_S30")


GSEA_KEGG_prep <- function(df) {
  #Get rid of rows with N/A in Fold Change column and make numeric
  df <- df[which(df$Fold.change!="#N/A"),]
  df$Fold.change <- as.numeric(df$Fold.change)
  ## feature 1: numeric vector, column with Fold Change values
  geneList = df[,6]
  ## feature 2: name the above vector using column with gene names
  names(geneList) = as.character(df[,10])
  ## feature 3: decreasing order
  geneList = sort(geneList, decreasing = TRUE)
  #run GSEA
  gk <- gseKEGG(geneList     = geneList, organism     = 'sdu',minGSSize    = 120,pvalueCutoff = 0.05)
  #Save results as a dataframe
  df_gk <- as.data.frame(gk)
  return(df_gk)
}

#Run a for loop on all the samples above and save as a uniquely named dataframe
for (i in names(y)) {
  df <- x[[i]]
  df2 <- GSEA_KEGG_prep(df)
  assign(paste("gsea_kegg", i, sep="_"), df2)
}
