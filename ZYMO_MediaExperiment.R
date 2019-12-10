#### Read data ####
# 45 min gradients for 5 different conditions, each with 3 replicates. Ran on Boudicca Lumos Fusion Instrument
# MaxQuant Output
proteinGroups_Zymo_MediaExperiment <- read.csv("P:/EAT_20190926_Zymomona/txt_Zymmo_MediaExperiment/proteinGroups_MediaExperiment.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Ran Imputation Algorithm developed my Dain. C# code output was filtered to 50% 
EAT_imputed_data <- read.csv("P:/EAT_20190926_Zymomona/txt_Zymmo_MediaExperiment/MediaExperiment_imputed_proteinGroups.csv")

#### Format data ####

# Rownames are the majority protein id, column names are cleaned up by removing everything before .*intensity
rownames(EAT_imputed_data) <- EAT_imputed_data$Majority.protein.ID 
colnames(EAT_imputed_data) <- c("Majority.protein.ID",sub(".*intensity.","",colnames(EAT_imputed_data[,-1])))

# Run subsetLFQ function from Data_Formating_Functions.R file
subset_MediaExp <- subsetLFQ(proteinGroups_Zymo_MediaExperiment)

# use id as identifier in row names
rownames(subset_MediaExp) <- subset_MediaExp$id

#### 50 %  Cut-off ####

# Perform feature removal on only the LFQ data only
clean_df_MediaExp <- remove.features.50percentcuttoff(subset_MediaExp[15:ncol(subset_MediaExp)])
# Remove Intensity from column names
colnames(clean_df_MediaExp) <- sub(".*intensity.", "", colnames(clean_df_MediaExp))

# All NA enteries will be replaced with an empty
clean_df_MediaExp[is.na(clean_df_MediaExp)] <- ""

# Create a meta object with only the meta data from MaxQuant
# [1] "Protein.IDs"                  "Majority.protein.IDs"         "Fasta.headers"                "Peptides"                    
# [5] "Razor...unique.peptides"      "Unique.peptides"              "Sequence.coverage...."        "Unique.sequence.coverage...."
# [9] "Mol..weight..kDa."            "Q.value"                      "Score"                        "Intensity"                   
# [13] "MS.MS.count"                  "id"   
meta_MediaExp <- subset_MediaExp[,1:14]

# Filtered meta data frame with the correct protein groups after 50% data cut off
filtered_MediaExp = meta_MediaExp[which(rownames(meta_MediaExp) %in% rownames(clean_df_MediaExp)),]

# Bind together meta data and data
ready_to_impute <- cbind(filtered_MediaExp,clean_df_MediaExp)

# object ready for imputation, output into the P folder 
write.csv(ready_to_impute, "P:/EAT_20190926_Zymomona/txt_Zymmo_MediaExperiment/EAT_50Percentfiltered_proteinGroups.csv")

# This object contains imputed data that is reduced to the 50% data cut off
MediaExp_merged <- EAT_imputed_data[which(EAT_imputed_data$Majority.protein.ID %in% filtered_MediaExp$id),]

# Bind together meta data and data
merged_MediaExp_impute <- cbind(filtered_MediaExp, MediaExp_merged[,-1])


#### Statistics ####

# reduce samples to the basic condition either 100 | 25 | 50 | Ctr | MM
replicates <- sapply(strsplit(colnames(clean_df_MediaExp), "_"), '[', 1) 
conditions <- unique(replicates) # reduce further to only unique conditions


# For loop that makes list of dataframes for each of the conditions
# First create an empty list then iterate through each unique condition
# Note I am using the most current imputed dataset with all the correct filters
# Save to each condition list the replicates
AllConditions_list <- list()
for(i in 1:length(conditions)){
  index <- which(grepl(conditions[i],colnames(merged_MediaExp_impute)))
  cond <- paste0("df_",conditions[i],collapse = "") 
  tmp <- merged_MediaExp_impute[,index]
  
  AllConditions_list[[cond]] <- tmp

}

# Create a function to calculate averages across replicates
average_fun <-  function(dataframe){
  rowMeans(dataframe,na.rm=TRUE) 
}
# apply average function to each dataframe in the list and add the output as a new column
average <- lapply(AllConditions_list, function(x) transform(x,Average=average_fun(x)))
# transform list of dataframe values to log 2 scale
average_log2 <- lapply(average, function(x) log2(x))
# standard deviation function
sd_fun <-  function(dataframe){
  apply(dataframe[,1:3],1,sd) 
}

# apply standard deviation to all data frames and add the output as a new column
sd_log2 <- lapply(average_log2, function(x) transform(x,SD = sd_fun(x)))

# fold change function
fc_fun <- function(dataframe){
  dataframe[,4]-sd_log2$df_Ctrl["Average"]
}

# applyl fold change funtion to all data frames and add the output as a new column
fc_log2 <- lapply(sd_log2, function(x) transform(x, FoldChange=fc_fun(x)))

# cbind all dataframes into a mega dataframe
allConditions_df <- cbind(fc_log2$df_MM, fc_log2$df_25, fc_log2$df_50, fc_log2$df_100, fc_log2$df_Ctrl)
samples <- sub("Average.1", "FoldChange",colnames(allConditions_df))
samples <- sub("X","", samples)
colnames(allConditions_df) <- samples

#### t.test #####
MM_ttest <-apply(allConditions_df[,c(1:3,25:27)],1,function(x) {t.test(x[1:3],x[4:6])$p.value})
MM_fdr <- p.adjust(MM_ttest, method = "fdr")
X25_ttest <- apply(allConditions_df[,c(7:9,25:27)],1,function(x) {t.test(x[1:3],x[4:6])$p.value})
X25_fdr <- p.adjust(X25_ttest,method = "fdr")
X50_ttest <- apply(allConditions_df[,c(13:15,25:27)],1,function(x) {t.test(x[1:3],x[4:6])$p.value})
X50_fdr <- p.adjust(X50_ttest, method = "fdr")
X100_ttest <- apply(allConditions_df[,c(19:21,25:27)],1,function(x) {t.test(x[1:3],x[4:6])$p.value})
X100_fdr <- p.adjust(X100_ttest, method = "fdr")

allConditions_df$MM_ttest <- MM_ttest
allConditions_df$MM_fdr <- MM_fdr
allConditions_df$'25_ttest' <- X25_ttest
allConditions_df$'25_fdr' <- X25_fdr
allConditions_df$'50_ttest' <- X50_ttest
allConditions_df$'50_fdr' <- X50_fdr
allConditions_df$'100_ttest' <- X100_ttest
allConditions_df$'100_fdr' <- X100_fdr

Collaborator_pdf <- cbind(filtered_MediaExp,allConditions_df)

# object ready for imputation, output into the P folder 
write.csv(Collaborator_pdf, "P:/EAT_20190926_Zymomona/txt_Zymmo_MediaExperiment/EAT_MediaExperiment_Average_SD_FC_TTEST_20191205proteinGroups.csv")
