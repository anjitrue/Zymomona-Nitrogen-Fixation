library(pcaMethods)
library(ggplot2)
library(plotly)
library(rapportools)
library(pheatmap)
library()

#browseVignettes("pcaMethods")

#### Upload Data ####

palette(c("#5ABCB2", "#D6F8D6","#D6F8D6" ,"mediumorchid2","darkgoldenrod1", 
          "hotpink3", "red2", "sienna2","slategray4","mediumturquoise", 
          "deepskyblue", "orangered", "midnightblue"))

# Read in MaxQuant output proteingroups, I filtered out the contaminants and the fasta header overspill

# data.frame of Eclipse 45 min runs using HP column. Searched with new aliqouts to confirm the strange intensity trends 
proteinGroups_Zymo_NitrogenFixation_Zm4Tag <- read.csv("H:/Projects/Proteomics/Zymomona/Martein,Julia_NitrogenFixation/combined_ZM4Tag/txt/proteinGroups_NitrogenFixation_ZM4Tag.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

meta_Zymo_NitrogenFixation <- read.csv("H:/Projects/Proteomics/Zymomona/Martein,Julia_NitrogenFixation/Meta_Zymo_NitrogenFixation.csv", 
                                       header = TRUE, sep = ",", stringsAsFactors = FALSE)
meta_Zymo_NitrogenFixation <- meta_Zymo_NitrogenFixation[-grep("A10",meta_Zymo_NitrogenFixation$Index),]

imputed_Zymo_Upshifted_NitrogenFixation_ZM4tag <- read.csv("H:/Projects/Proteomics/Zymomona/Martein,Julia_NitrogenFixation/R_outputs/NitrogenFixation_UpshiftImputed_ZM4tag_proteinGroups.csv", 
                                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)
imputed_Zymo_Downshifted_NitrogenFixation_ZM4tag <- read.csv("H:/Projects/Proteomics/Zymomona/Martein,Julia_NitrogenFixation/R_outputs/NitrogenFixation_DownshiftImputed_ZM4tag_proteinGroups.csv", 
                                                           header = TRUE, sep = ",", stringsAsFactors = FALSE)
#### Data Wrangling #####
# Perform subset function
subset_NitrogenFixation <- subsetLFQ(proteinGroups_Zymo_NitrogenFixation_Zm4Tag)

# use id as identifier in row names
rownames(subset_NitrogenFixation) <- subset_NitrogenFixation$id

# remove the blank runs 
#subset_NitrogenFixation <- subset_NitrogenFixation[,-grep("lank",colnames(subset_NitrogenFixation))]

# Isolate only the first protein in protein group for Majority.protein.IDs and Protein.IDs
majorityProtein <- lapply(strsplit(subset_NitrogenFixation$Majority.protein.IDs, ";"), '[', 1)
protein.ID <- lapply(strsplit(subset_NitrogenFixation$Protein.IDs, ";"), '[', 1)

# Replace with reduced array
subset_NitrogenFixation$Protein.IDs <- protein.ID
subset_NitrogenFixation$Majority.protein.IDs <- majorityProtein


# remove protein groups that are missing over 50% of the data across samples
remove.features.50percentcuttoff <- function (x) {
  
  # Calculate how many missing values per feature
  features.missing = rowMeans(is.na(x)) 
  print(paste0("Number of protein groups that have over 50% missing measurements: ",sum(features.missing > 0.50))) 
  features.missing.50more = rownames(x)[features.missing > 0.50] 
  
  # create a vector of the features that meet criteria
  keep.features = rownames(x)[features.missing <= 0.50]
  print(paste0("Protein groups that pass the 50% filteration: ", length(keep.features)))
  #names(keep.features) = keep.features 
  
  # create a vector of the features that will be removed from data.frame
  remove.features = rownames(x)[features.missing > 0.50]
  print(paste0("Number of protein groups removed from dataset: ", length(remove.features)))
  #names(remove.features) = remove.features
  
  # If there isn't any features to remove keep all features
  if(sum(features.missing) == 0 ){
    filtered = x[which(rownames(x) %in% keep.features),]
    # otherwise filter out the features that contain over 50% of data missing
  } else{
    filtered = x[-which(rownames(x) %in% remove.features),]
  }
  return(filtered)
}


# Perform feature removal on only the LFQ data only
clean_df_NitrogenFixation <- remove.features.50percentcuttoff(subset_NitrogenFixation[,15:ncol(subset_NitrogenFixation)])
# Remove Intensity from column names
colnames(clean_df_NitrogenFixation) <- sub(".*_EAT_", "", colnames(clean_df_NitrogenFixation))

# Create a meta object with only the meta data from MaxQuant
meta_proteins_NitrogenFixation <- subset_NitrogenFixation[,1:14]

# Filtered meta data frame with the correct protein groups after 50% data cut off
filtered_proteinsMeta_NitrogenFixation = meta_proteins_NitrogenFixation[which(rownames(meta_proteins_NitrogenFixation) %in% rownames(clean_df_NitrogenFixation)),]

# Transpose cleaned up dataframe with intensity values
t_clean_df_NitrogenFixation <- t(clean_df_NitrogenFixation)

sample_number <- vector()
for(i in 1:nrow(t_clean_df_NitrogenFixation)){
  
  row_number <- meta_Zymo_NitrogenFixation[grep(rownames(t_clean_df_NitrogenFixation)[i], meta_Zymo_NitrogenFixation$BarCode),2]
  
  sample_number <- append(sample_number, row_number)
}

rownames(t_clean_df_NitrogenFixation) <- sample_number
clean_df_NitrogenFixation <- t(t_clean_df_NitrogenFixation)
clean_df_NitrogenFixation <- clean_df_NitrogenFixation[,order(as.numeric(as.character(colnames(clean_df_NitrogenFixation))))]
clean_df_NitrogenFixation <- cbind(filtered_proteinsMeta_NitrogenFixation, clean_df_NitrogenFixation)
clean_df_NitrogenFixation <- as.data.frame(clean_df_NitrogenFixation)

clean_df_NitrogenFixation$Protein.IDs <- vapply(clean_df_NitrogenFixation$Protein.IDs, paste, collapse = ", ", character(1L))
clean_df_NitrogenFixation$Majority.protein.IDs <- vapply(clean_df_NitrogenFixation$Majority.protein.IDs, paste, collapse = ", ", character(1L))

colnames(clean_df_NitrogenFixation) <- c(colnames(clean_df_NitrogenFixation)[1:14], meta_Zymo_NitrogenFixation$sample.name)

#write.csv(clean_df_NitrogenFixation, "H:/Projects/Proteomics/Zymomona/Martein,Julia_NitrogenFixation/R_outputs/EAT_50PercentFilter_NitrogenFixation_UpAndDownShiftData_ZM4tag.csv")

# DRB algorithm imputed data matched with cleaned up 50% filtered data (I don't think this is an important step but will keep it)
NitrogenFixation_upshift <- imputed_Zymo_Upshifted_NitrogenFixation_ZM4tag[which(imputed_Zymo_Upshifted_NitrogenFixation_ZM4tag$Majority.protein.ID %in% clean_df_NitrogenFixation$Majority.protein.IDs),]
NitrogenFixation_downshift <- imputed_Zymo_Downshifted_NitrogenFixation_ZM4tag[which(imputed_Zymo_Downshifted_NitrogenFixation_ZM4tag$Majority.protein.ID %in% clean_df_NitrogenFixation$Majority.protein.IDs),]

# Bind together meta data and data
NitrogenFixation_upshift_merged <- cbind(filtered_proteinsMeta_NitrogenFixation,NitrogenFixation_upshift[,-1])

NitrogenFixation_downshift_merged <- cbind(filtered_proteinsMeta_NitrogenFixation, NitrogenFixation_downshift[,-1])


#### Data Analysis ####
#Upshift Data
NitrogenFixation_upshift_merged_log2_matrix <- log2(as.matrix(NitrogenFixation_upshift_merged[,-c(1:14)]))
rownames(NitrogenFixation_upshift_merged_log2_matrix) <- NitrogenFixation_upshift_merged$Protein.IDs
pheatmap::pheatmap(NitrogenFixation_upshift_merged_log2_matrix,
                   scale = "row",
                   cluster_cols = FALSE)
#DownShift Data
NitrogenFixation_downshift_merged_log2_matrix <- log2(as.matrix(NitrogenFixation_downshift_merged[,-c(1:14)]))
rownames(NitrogenFixation_downshift_merged_log2_matrix) <- NitrogenFixation_downshift_merged$Protein.IDs
pheatmap::pheatmap(NitrogenFixation_downshift_merged_log2_matrix,
                   scale = "row",
                   cluster_cols = FALSE)


#### Statistics ####
# transforms data frames such that samples are rows and features are columns
t_NF_upshift_noMeta <- t(NitrogenFixation_upshift_merged[,-c(1:14)])
rownames(t_NF_upshift_noMeta) <- meta_Zymo_NitrogenFixation$sample.name[1:23]

# in order to aggregate by sample name, use numeric matrix as R object, by Sample Name from a different object
# and the function of choice
average_time.point_NF_upshift <- aggregate(t_NF_upshift_noMeta,
                                              list(rownames(t_NF_upshift_noMeta)), mean)

# rename the colnames and rown names for future re-ordering
colnames(average_time.point_NF_upshift) <- c("Sample.Name", NitrogenFixation_upshift_merged$Majority.protein.IDs)
rownames(average_time.point_NF_upshift) <- average_time.point_NF_upshift$Sample.Name

# order samples by substring "time"
reorder_samples_time <- function(x){
  m = x[,1] # Extract sample names
  
  r= strsplit(m,"-") # Split the column names based on the -
  
  # Create a data.frame with four columns for each of the string split characters
  # X1 = Upshift(U)/Downshift(D), X2 = N2/NH4 , X3 = Replicate (A,B,C), X4 = time (min) 
  q = data.frame(matrix(unlist(r),nrow=length(m),byrow=T)) 
  q$X4 <- as.numeric(as.character(q$X4))
  
  q = q[order(q$X2,as.numeric(as.character(q$X4))),] # order by X2 and X4
  
  # Reformat names
  m = paste0(q$X1,"-",q$X2,"-",q$X3, "-", q$X4)
  
  # Rearrange columns
  x = x[m,]
  
  return(x)
}


