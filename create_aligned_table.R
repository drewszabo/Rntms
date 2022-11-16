create_aligened_features <- function(fGroups = fGroupsXCMS, csv_path = "aligned_feature_table.csv"){

# From NeatMS Documentation
# Feature table export from XCMS (aligned peaks with gapfilling)


# The first part of the code is the same as for unaligned peaks, you can jump to the feature information addition 

# This code assumes that the xdata variable corresponds 
# to the XCMSnExp object that contains the detected peaks 

# Load dplyr (required for left_join())
library(tidyverse)
# Load tibble (required for rownames_to_column())
library(tibble)
# Load XCMS
library(xcms)

#--------------
# Extract XCMS object from fGroups (Drew)
#--------------
xdata_filled <- fGroups@xdata

# Create dataframe containing adjusted retention times 
df_chrom_peaks_fill <- as.data.frame(chromPeaks(xdata_filled))
# Create dataframe with raw retention times 
feature_dataframe <- as.data.frame(chromPeaks(dropAdjustedRtime(xdata_filled)))

# Get the peaks that have been recovered by the gapfilling step
df_filled_peaks <- df_chrom_peaks_fill[!row.names(df_chrom_peaks_fill) %in% c(rownames(feature_dataframe)),]

# Add them to the raw retention time dataframe (filled peaks do not have raw retention times so we use the adjusted ones)
feature_dataframe <- bind_rows(feature_dataframe, df_filled_peaks)

# Rename the retention time columns of the adjusted rt dataframe
df_chrom_peaks_fill <- df_chrom_peaks_fill %>%
  dplyr::rename(rt_adjusted = rt, rtmin_adjusted = rtmin, rtmax_adjusted = rtmax)

# Add the adjusted rt columns to the dataframe containing raw rt
# To have a dataframe containing both raw and adjusted rt information
feature_dataframe <- left_join(rownames_to_column(feature_dataframe), rownames_to_column(df_chrom_peaks_fill[,c("rt_adjusted","rtmin_adjusted","rtmax_adjusted")]), by="rowname")

# Remove the rownames as we won't need them
feature_dataframe$rowname <- NULL

# Retrieve the sample names and store them as a dataframe
sample_names_df <- as.data.frame(sampleNames(xdata_filled))

#------------------
# Add sample names manually from patRoon data (Drew)
#------------------
sample_names_df <- as.data.frame(fGroups@analysisInfo$analysis)
sample_names_df$`fGroups@analysisInfo$analysis` <- paste0(sample_names_df$`fGroups@analysisInfo$analysis`, ".mzML")

# Rename the unique column "sample_name"
colnames(sample_names_df) <- c("sample_name")

# Generate the correct sample ids for matching purposes
# XCMS sampleNames() function returns sample names ordered by their ids
sample_names_df$sample <- seq.int(nrow(sample_names_df))

# Attach the sample names to the main dataframe by matching ids (sample column)
feature_dataframe <- left_join(feature_dataframe,sample_names_df, by="sample")

### Feature information addition ###

# Here we will bring the feature alignment information stored in the XCMSnExp object to the dataframe that we have already created

featuresDef <- featureDefinitions(xdata_filled)
featuresDef_df = data.frame(featuresDef)

#--------------------
# Added as.numeric() to fix the column_index assignment (Drew)
#--------------------

# Adjust variable
# Only keep the information we need (column named 'peakidx')
# Get the index of the peakidx column
column_index <- as.numeric(which(colnames(featuresDef_df)=="peakidx"))
# Extract the peakidx column
features_df <- as.list(featuresDef_df$peakidx)

for(x in 1:length(features_df)){
  length(features_df[[x]]) <- length(fGroups@groups[[x]])
}

features_df = data.frame(features_df)

features_df = data.frame(t(features_df))

# Rename the column
peak_colummn_name <- colnames(features_df)
features_df = rename(features_df, "peak_id"=all_of(peak_colummn_name))

#------------------
# Add feature names from patRoon to features_df (Drew)
#------------------

feature_names <- colnames(fGroups@groups)

features_df <- cbind(feature_id = feature_names, features_df)

features_df <- features_df %>%
  pivot_longer(!feature_id, values_to = "peak_id") %>%
  select(feature_id, peak_id) %>%
  drop_na(peak_id)


# We'll use data.table for the next step
require(data.table)

# Get all the peak_id for each feature_id
features_df <- data.table(features_df)
features_df = features_df[, list(peak_id = unlist(peak_id)), by=feature_id]

# Bring the feature_id to the original peak dataframe
feature_dataframe = cbind(peak_id= row.names(feature_dataframe),feature_dataframe)
feature_dataframe$peak_id = as.character(feature_dataframe$peak_id)
features_df$peak_id = as.character(features_df$peak_id)
feature_dataframe = left_join(feature_dataframe, features_df, by="peak_id")

# Note: The dataframe contains an extra column called peak_id, but this won't affect NeatMS and will simply be ignored (as would any other column not present in the list above).

# Export the data as csv. 
# Note: Set row.names to FALSE as NeatMS does not need them
write.csv(feature_dataframe, csv_path, row.names = FALSE)

return(feature_dataframe)

}
