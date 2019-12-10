library(installr) # for is.empty function
#BiocManager::install("pcaMethods")
library(pcaMethods)
library(ggplot2)
library(plotly)

#browseVignettes("pcaMethods")

#### Upload Data ####

palette(c("#5ABCB2", "#D6F8D6","#D6F8D6" ,"mediumorchid2","darkgoldenrod1", 
          "hotpink3", "red2", "sienna2","slategray4","mediumturquoise", 
          "deepskyblue", "orangered", "midnightblue"))

# Read in MaxQuant output proteingroups, I filtered out the contaminants and the fasta header overspill

# data.frame of Eclipse 45 min runs using HP column. Searched with new aliqouts to confirm the strange intensity trends 
proteinGroups_Zymo_Eclipse45min <- read.csv("P:/EAT_20190926_Zymomona/txt_Eclipse_45min/proteinGroups.csv", 
                                              header = TRUE, sep = ",", stringsAsFactors = FALSE)
peptides_zymo <- read.csv("P:/EAT_20190926_Zymomona/txt_Eclipse_45min/peptides.csv", 
                          header = TRUE, sep = ",",stringsAsFactors = FALSE)

rownames(DRB_imputed_data) <- DRB_imputed_data$Majority.protein.ID

#### ISPG peptides #####
#subset the protein groups that are associated to the ISPG gene
ISPG_data <- proteinGroups_Zymo_Eclipse45min[grep("GN=ispG", proteinGroups_Zymo_Eclipse45min$Fasta.headers),]

ISPG_lfq <- ISPG_data[3,c(1,grep("LFQ",colnames(ISPG_data)))]

Sample <- sub("LFQ.intensity.", "", colnames(ISPG_lfq[,-1]))
Intensity <- unlist(ISPG_lfq[2:ncol(ISPG_lfq)])
ISPG_lfq <- data.frame("Sample" = Sample, "Intensity" = Intensity)

ISPG_lfq <- ISPG_lfq[-grep("new",ISPG_lfq$Sample),]
ISPG_lfq$LOG2 <- log2(ISPG_lfq$Intensity)
ISPG_lfq$Type <- c("ISPG","ISPG","ISPG","ISPH","ISPH","ISPH","WT","WT","WT")


type.colors = as.numeric(factor(ISPG_lfq$Type))
q <- ggplot(ISPG_lfq, aes(x=Sample, y=LOG2, fill=Type)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Log2(ISPG protein)", 
                                                                    subtitle = "Across sample types in triplicate",
                                                                    y = "Log2(Intensity)")

#subset zymo peptides for protein Q5NR50
ISPG_peptides <- peptides_zymo[grep("Q5NR50", peptides_zymo$Proteins),]

ISPG_peptides_averages <- ISPG_peptides[,grep("IspG", colnames(ISPG_peptides))]
sub_ISPG_peptides_avg <- ISPG_peptides_averages[,10:12]
rownames(sub_ISPG_peptides_avg) <- ISPG_peptides$Sequence

sub_ISPG_peptides_avg$Avg <- rowMeans(sub_ISPG_peptides_avg)
sub_ISPG_peptides_avg$Log2Avg <- log2(sub_ISPG_peptides_avg$Avg)

#subset peptides with abundance greater than log2
ISPG_peptides_abundance <- ISPG_peptides[which(log2(ISPG_peptides$Intensity.IspG_2b__20191004) > 34),]



# Plot all peptide abundances for protein "Q5NR50" aka ISPG
ggplot(ISPG_peptides, aes(x=Sequence, y= log2(Intensity.IspG_2b__20191004))) + 
  geom_point() +
  ylim(20,36) +
  theme_light()+
  labs(title = "Peptide Abundances" ,subtitle = "IspG_2b__20191004 Log2 Transformed", y = "Log2(Intensity)")

# Label peptides with abundance greater than 34
ggplot(ISPG_peptides, aes(x=Sequence, y= log2(Intensity.IspG_2b__20191004))) + 
  geom_point() +
  ylim(20,36) +
  scale_y_continuous(breaks = seq(20,36,2)) +
  geom_text(data=subset(ISPG_peptides, log2(Intensity.IspG_2b__20191004) > 34),
            aes(x=Sequence, y= log2(Intensity.IspG_2b__20191004),label= Sequence)) +
  theme_light()+
  labs(title = "Peptide Abundances greater than 34" ,subtitle = "IspG_2b__20191004 Log2 Transformed", y = "Log2(Intensity)")

# averaged intensities across the three replicates for the peptides
ggplot(sub_ISPG_peptides_avg, aes(x=rownames(sub_ISPG_peptides_avg), y= Log2Avg)) + 
  geom_point() +
  ylim(20,36) +
  theme_light()+
  labs(title = "Average Peptide Abundances" ,subtitle = "IspG Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides")


# Averaged intensities with intensities greater than 33
sub_ISPG_peptides_avg$sequence <- rownames(sub_ISPG_peptides_avg)

ggplot(sub_ISPG_peptides_avg, aes(x=rownames(sub_ISPG_peptides_avg), y= Log2Avg)) + 
  geom_point() +
  #ylim(20,36) +
  geom_text(data = subset(sub_ISPG_peptides_avg, Log2Avg > 33),
            aes(x=sequence, y= Log2Avg, label= sequence)) +
  theme_light()+
  labs(title = "Average Peptide Abundances" ,subtitle = "IspG Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides")


#### ISPH peptides #####
#subset the protein groups that are associated to the ISPH gene
ISPH_data <- proteinGroups_Zymo_Eclipse45min[grep("GN=ispH", proteinGroups_Zymo_Eclipse45min$Fasta.headers),]

ISPH_lfq <- ISPH_data[,c(1,grep("LFQ",colnames(ISPH_data)))]

Sample <- sub("LFQ.intensity.", "", colnames(ISPH_lfq[,-1]))
Intensity <- unlist(ISPH_lfq[2:ncol(ISPH_lfq)])
ISPH_lfq <- data.frame("Sample" = Sample, "Intensity" = Intensity)

ISPH_lfq <- ISPH_lfq[-grep("new",ISPH_lfq$Sample),]
ISPH_lfq$LOG2 <- log2(ISPH_lfq$Intensity)
ISPH_lfq$Type <- c("ISPG","ISPG","ISPG","ISPH","ISPH","ISPH","WT","WT","WT")


type.colors = as.numeric(factor(ISPH_lfq$Type))
q <- ggplot(ISPH_lfq, aes(x=Sample, y=LOG2, fill=Type)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Log2(ISPH protein)", 
                                                                    subtitle = "Across sample types in triplicate",
                                                                    y = "Log2(Intensity)")
#subset zymo peptides for protein Q5NR50
ISPH_peptides <- peptides_zymo[grep("Q5NP61", peptides_zymo$Proteins),]

ISPH_peptides_averages <- ISPH_peptides[,grep("IspH", colnames(ISPH_peptides))]
sub_ISPH_peptides_avg <- ISPH_peptides_averages[,8:10]
rownames(sub_ISPH_peptides_avg) <- ISPH_peptides$Sequence

sub_ISPH_peptides_avg$Avg <- rowMeans(sub_ISPH_peptides_avg)
sub_ISPH_peptides_avg$Log2Avg <- log2(sub_ISPH_peptides_avg$Avg)

sub_ISPH_peptides_avg <- sub_ISPH_peptides_avg[-which(sub_ISPH_peptides_avg$Avg == 0),]

#subset peptides with abundance greater than log2
ISPH_peptides_abundance <- sub_ISPH_peptides_avg[which(sub_ISPH_peptides_avg$Log2Avg > 32),]

# ISPH averaged intensities across the three replicates for the peptides
ggplot(sub_ISPH_peptides_avg, aes(x=rownames(sub_ISPH_peptides_avg), y= Log2Avg)) + 
  geom_point() +
  ylim(20,34) +
  theme_light()+
  labs(title = "Average Peptide Abundances" ,subtitle = "IspH Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides")


# ISPH averaged intensities with intensities greater than 33
sub_ISPH_peptides_avg$sequence <- rownames(sub_ISPH_peptides_avg)

ggplot(sub_ISPH_peptides_avg, aes(x=rownames(sub_ISPH_peptides_avg), y= Log2Avg)) + 
  geom_point() +
  #ylim(20,36) +
  geom_text(data = subset(sub_ISPH_peptides_avg, Log2Avg > 32),
            aes(x=sequence, y= Log2Avg, label= sequence)) +
  theme_light()+
  labs(title = "Average Peptide Abundances" ,subtitle = "IspH Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides")





