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

DRB_imputed_data <- read.csv("T:/File_exchange/Dain/From Anji/DRB_imputed_data.csv")
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




#### Format Data ####

# subset specific columns and bind function
subsetLFQ <- function(x){
  y<- x[,which(names(x) %in% c("Protein.IDs", 
                               "Majority.protein.IDs",
                               "Fasta.headers",
                               "Peptides",
                               "Razor...unique.peptides",
                               "Unique.peptides",
                               "Sequence.coverage....",
                               "Unique.sequence.coverage....",
                               "Mol..weight..kDa.",
                               "Q.value",
                               "Score",
                               "Intensity",
                               "MS.MS.count",
                               "id"))]
  z <- x[,grep("LFQ.intensity.",names(x))]
  z[z == 0] <- NA
  x <- cbind(y,z)
  return(x)
}

# Perform subset function
subset_df <- subsetLFQ(proteinGroups_Zymo_Eclipse45min)

# remove the new aliqout runs - these were purely for confirming the crazy inverted trend
subset_df <- subset_df[,-grep("new",colnames(subset_df))]
# convert Protein.IDs to characters 
subset_df$Protein.IDs <- as.character(subset_df$Protein.IDs)
# convert Majority.protein.IDs to characters
subset_df$Majority.protein.IDs <- as.character(subset_df$Majority.protein.IDs)

# use id as identifier in row names
rownames(subset_df) <- subset_df$id

# Isolate only the first protein in protein group for Majority.protein.IDs and Protein.IDs
majorityProtein <- lapply(strsplit(subset_df$Majority.protein.IDs, ";"), '[', 1)
protein.ID <- lapply(strsplit(subset_df$Protein.IDs, ";"), '[', 1)

# Replace with reduced array
subset_df$Protein.IDs <- protein.ID
subset_df$Majority.protein.IDs <- majorityProtein


# remove protein groups that are missing over 50% of the data across samples
remove.features.50percentcuttoff <- function (x) {
  
  # Calculate how many missing values per feature
  features.missing = rowMeans(is.na(x)) 
  print(paste0("Number of protein groups that have over 50% missing measurements: ",sum(features.missing > 0.50))) 
  features.missing.50more = rownames(x)[features.missing > 0.50] 
  
  # create a vector of the features that meet criteria
  keep.features = which(features.missing <= 0.50) 
  print(paste0("Protein groups that pass the 50% filteration: ", length(keep.features)))
  names(keep.features) = keep.features 
  
  # create a vector of the features that will be removed from data.frame
  remove.features = which(features.missing > 0.50)
  print(paste0("Number of protein groups removed from dataset: ", length(remove.features)))
  names(remove.features) = remove.features
  
  # If there isn't any features to remove keep all features
  if(is.empty(remove.features) == TRUE ){
    filtered = x[which(rownames(x) %in% keep.features),]
  # otherwise filter out the features that contain over 50% of data missing
  } else{
      filtered = x[-which(rownames(x) %in% remove.features),]
    }
  return(filtered)
}

#### 50 %  Cut-off ####

# Perform feature removal on only the LFQ data only
clean_df <- remove.features.50percentcuttoff(subset_df[15:ncol(subset_df)])
# Remove Intensity from column names
colnames(clean_df) <- sub(".*intensity.", "", colnames(clean_df))

# Create a meta object with only the meta data from MaxQuant
meta <- subset_df[,1:14]

# Filtered meta data frame with the correct protein groups after 50% data cut off
filtered = meta[which(rownames(meta) %in% rownames(clean_df)),]

# DRB imputed data will be reduced to the protein groups left after 50% data cut off
drb_merged <- DRB_imputed_data[which(DRB_imputed_data$Majority.protein.ID %in% filtered$Majority.protein.IDs),]

# Bind together meta data and data
ready_to_impute <- cbind(filtered,clean_df)
merged_drb_impute <- cbind(filtered, drb_merged)

#### Statistics ####

#### IspG Stats ####
IspG_average <- drb_merged[,grep("IspG", colnames(drb_merged))]
IspG_average$Average.IspG <- rowMeans(IspG_average)
IspG_average$Std.IspG <- apply(IspG_average[1:3],1,sd)


#### IspH Stats ####
IspH_average <- drb_merged[,grep("IspH", colnames(drb_merged))]
IspH_average$Average.IspH <- rowMeans(IspH_average)
IspH_average$Std.IspH <- apply(IspH_average[1:3],1,sd)

#### WT Stats ####
WT_average <- drb_merged[,grep("WT", colnames(drb_merged))]
WT_average$Average.WT <- rowMeans(WT_average)
WT_average$Std.WT <- apply(WT_average[1:3],1,sd)

#### Fold Change ####
IspG_average$FC_IspG_WT <- IspG_average$Average.IspG-WT_average$Average.WT

IspH_average$FC_IspH_WT <- IspH_average$Average.IspH-WT_average$Average.WT


#### Append columns to imputed data set ####
merged_drb_impute$Average.IspG <- IspG_average$Average.IspG
merged_drb_impute$Average.IspH <- IspH_average$Average.IspH
merged_drb_impute$Average.WT <- WT_average$Average.WT


merged_drb_impute$FC_IspG_WT <- IspG_average$FC_IspG_WT
merged_drb_impute$FC_IspH_WT <- IspH_average$FC_IspH_WT

merged_drb_impute$Std.IspG <- IspG_average$Std.IspG
merged_drb_impute$Std.IspH <- IspH_average$Std.IspH
merged_drb_impute$Std.WT <- WT_average$Std.WT

#### t.test #####
with(IspG_average[1:3], rep)
merged_drb_impute$t.test.Ispg <- apply(c(IspG_average[1:3], WT_average[1:3]),1,function(x) {t.test(x[1:3], x[4:6])$p.value})

#### Transform lists into vectors for a complete data frame
merged_drb_impute$Protein.IDs <- vapply(merged_drb_impute$Protein.IDs, paste, collapse = ",", character(1L))
merged_drb_impute$Majority.protein.IDs <- vapply(merged_drb_impute$Majority.protein.IDs, paste, collapse = ",", character(1L))

write.csv(merged_drb_impute, "P:/EAT_20190926_Zymomona/txt_Eclipse_45min/Zymo20191008_Mean_STD_FC_proteinGroups.csv")
write.csv(ready_to_impute, "P:/EAT_20190926_Zymomona/txt_Eclipse_45min/EAT_50percentfiltered_proteinGroups.csv")


p <- ggplot(merged_drb_impute, aes(y=Average.IspH, x = Average.IspG)) +
  geom_point()

ggplotly(p)

merged_drb_impute[grep(32.16088, merged_drb_impute$Average.IspG),]

p <- ggplot(merged_drb_impute, aes(y=Average.IspH, x = Average.WT)) +
  geom_point()

ggplotly(p)


p <- ggplot(merged_drb_impute, aes(y=Average.IspG, x = Average.WT)) +
  geom_point()

ggplotly(p)

merged_drb_impute[grep(37.704, merged_drb_impute$Average.IspG),]

ggplot(merged_drb_impute, aes(x=pvalues)) + 
        geom_histogram(binwidth = .05, color="black", fill="white")

# t.test on all the data, not by row
tt <- with(df, 
           t.test(merged_drb_impute$Average.IspG, merged_drb_impute$Average.WT,
           var.equal = TRUE))




# myFun <- function(data) {
#   temp1 <- sapply(data, is.list)
#   temp2 <- do.call(
#     cbind, lapply(data[temp1], function(x) 
#       data.frame(do.call(rbind, x), check.names=FALSE)))
#   cbind(data[!temp1], temp2)
# }

# temp1 <- sapply(proteinGroups_Zymo_Eclipse45min, is.list)
# temp2 <- do.call(
#   cbind, lapply(drb_merged[temp1], function(x) 
#     data.frame(do.call(rbind, x), check.names=FALSE)))


df <- data.frame(matrix(unlist(drb_merged[,-1]), nrow=9, byrow=T),stringsAsFactors=FALSE)
df <- t(df)
colnames(df) <- colnames(drb_merged[,-1])
rownames(df) <- drb_merged[,1]

pvals=apply(df,1,function(x) {t.test(x[1:3],x[7:9])$p.value})

merged_drb_impute$pvalues <- pvals

var(df[1:3])
var(df[7:9])

hist(merged_drb_impute[16:24], breaks = 50)
boxplot(log2(clean_df[,1:ncol(clean_df)]))

matrix_clean_df <- matrix(unlist(clean_df),ncol = 9, byrow = FALSE)
colnames(matrix_clean_df) <- colnames(clean_df)
rownames(matrix_clean_df) <- rownames(clean_df)

# no missing values allowed #resPCA <- pca(matrix_clean_df, method="svd", center=FALSE, nPcs=5)
resPPCA <- pca(matrix_clean_df, method="ppca", center=FALSE, nPcs=5)
resBPCA <- pca(matrix_clean_df, method="bpca", center=FALSE, nPcs=5)
resSVDI <- pca(matrix_clean_df, method="svdImpute", center=FALSE, nPcs=5)
resNipals <- pca(matrix_clean_df, method="nipals", center=FALSE, nPcs=5)
resNLPCA <- pca(matrix_clean_df, method="nlpca", center=FALSE, nPcs=5, maxSteps=300)

