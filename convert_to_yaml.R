convert_to_yaml <- function(ntms_results = "ntms_export.csv", yaml_path = "model_session.yml"){

# Import results
neatms_export <- read.csv(ntms_results)

# Re-allocate feature_id from feature_dataframe based on "into" variable
neatms_export <- neatms_export %>%
  dplyr::rename(mz = m.z,
                maxo = height,
                into = area,
                rt = retention.time,
                analysis = sample) %>%
  dplyr::mutate(rt = rt*60)

neatms_export <- neatms_export %>%
  dplyr::left_join(anaInfo) %>%
  dplyr::select(-path, -blank)

neatms_export$feature_id <- 
  feature_dataframe$feature_id[match(round(neatms_export$into,2), round(feature_dataframe$into,2))]

# Categorize features by quality
neatms_export <- neatms_export %>%
  dplyr::group_by(feature.ID, group) %>%
  dplyr::summarise(quality = ifelse("High_quality" %in% label, TRUE, FALSE),
            feature = feature_id[1])

# Filter features with <1 "high-quality"
removeFully <- neatms_export %>%
  dplyr::group_by(feature) %>%
  dplyr::filter(sum(quality)/length(quality) == 0) %>%
  dplyr::distinct(feature) %>%
  dplyr::select(feature) %>%
  dplyr::rename(removeFully = feature)


# Create YAML file for removal of noisy features
library(yaml)

list <- as.list(removeFully)

# Set the version field
version <- 1.0 # Should be written without quotation marks

# Set the type field
type <- "featureGroups"

# Set the removePartially field
removePartially <- list()
featureGroups <- list()

# Add the version, type, and removePartially fields to the list
list$removePartially <- removePartially
list$version <- version
list$type <- type

featureGroups <- patRoon::as.data.table(fGroups) %>%
  dplyr::select(group, ret, mz) %>%
  dplyr::filter(group %in% removeFully$removeFully)
featureGroups <- list(featureGroups = split(replace(featureGroups, "group", NULL), featureGroups$group))
list <- append(list, featureGroups)

# Convert the list to a YAML object
yaml_object <- as.yaml(list)

# Write the YAML object to a file
write(yaml_object, file = yaml_path)

}
