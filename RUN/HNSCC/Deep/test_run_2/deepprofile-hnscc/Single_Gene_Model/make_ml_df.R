


# Generalized function to process gene list and count data
process_gene_data <- function(gene_list, annotated_file, count_file) {
  
  # Read the CSV file with annotated results
  find_id <- read.csv(annotated_file)
  
  # Check if the file was loaded correctly
  print("Annotated File Preview:")
  print(head(find_id))
  
  # Find matching rows where `external_gene_name` matches any of the `gene_list`
  matched_genes <- find_id[find_id$external_gene_name %in% gene_list, ]
  
  # Extract corresponding `Given_id` for those matched genes
  given_ids <- matched_genes$Given_id
  
  # Read the count data
  data <- read.csv(count_file)
  
  # Find matching rows where `gene_id` in count data matches `Given_id`
  data_ids <- data[data$gene_id %in% given_ids, ]
  
  # Check if the data was filtered correctly
  print("Filtered Count Data Preview:")
  print(head(data_ids))
  
  # Replace `gene_id` with the corresponding gene name from `matched_genes`
  data_ids$gene_id <- matched_genes$external_gene_name[match(data_ids$gene_id, matched_genes$Given_id)]
  
  # Display the processed data
  print("Processed Data with Gene Names:")
  print(head(data_ids))
  
  # Return the processed data
  return(data_ids)
}

# Example usage:
gene_list <- c("CCND1", "HSD17B10", "LPAR3", "SNRPC", "HSD17B4", "TBC1D4", "CTSD", "ZNF451", 
               "GATA6", "PLAGL1", "ICAM5", "RRAD", "ZNF212", "DHX38", "GDF9", "BCL10", "MGRN1", 
               "SLC9A8", "PHTF1", "NUP58", "PIK3CD", "NFAT5", "RFPL3S", "CARS1", "UBXN7", 
               "DIO2", "PKNOX1", "MSI1", "ZNF266", "SZT2", "DMPK", "NEDD4L", "FOSL1", "SIN3B")

annotated_file <- "annotated_7_LFC_Monocyte_Mice_FROM_36_Samples.csv"
count_file <- "2_count_Data_Normalized_124_Imputed_MICE.csv"

# Call the function
processed_data <- process_gene_data(gene_list, annotated_file, count_file)

# Create a new filename using the original count file name, adding "processed_" as a prefix
output_file <- paste0("processed_", count_file)

# Save the processed data to a CSV file
write.csv(processed_data, file = output_file, row.names = FALSE)











# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


# take gene list
gene_ids <- c("CCND1", "HSD17B10", "LPAR3", "SNRPC", "HSD17B4", "TBC1D4", "CTSD", "ZNF451", "GATA6",
              "PLAGL1", "ICAM5", "RRAD", "ZNF212", "DHX38", "GDF9", "BCL10", "MGRN1", "SLC9A8",
              "PHTF1", "NUP58", "PIK3CD", "NFAT5", "RFPL3S", "CARS1", "UBXN7", "DIO2", "PKNOX1",
              "MSI1", "ZNF266", "SZT2", "DMPK", "NEDD4L", "FOSL1", "SIN3B")

# Create a data frame
gene_df <- data.frame(Gene = gene_ids)
gene_df <- read.csv("genes_list.csv",)

# write.csv(gene_df, "genes_list.csv", row.names = FALSE)
find_id <- read.csv ("annotated_resLFC_Treatment_MPXV.clade.I.infected_vs_mock.csv")
head(find_id)
#find id <- find_id$external_gene_name[ gene_ids,]
# i have gene_ids in the column of find_id named external_gene_name, i want to find out the corresponding Given_ids

given_id <- find_id[find_id$external_gene_name %in% find_id$Given_id ]




# Read your gene list (if it's coming from a CSV)
gene_name <- c("CCND1", "HSD17B10", "LPAR3", "SNRPC", "HSD17B4", "TBC1D4", "CTSD", "ZNF451",
               "GATA6", "PLAGL1", "ICAM5", "RRAD", "ZNF212", "DHX38", "GDF9", "BCL10", "MGRN1",
               "SLC9A8", "PHTF1", "NUP58", "PIK3CD", "NFAT5", "RFPL3S", "CARS1", "UBXN7",
               "DIO2", "PKNOX1", "MSI1", "ZNF266", "SZT2", "DMPK", "NEDD4L", "FOSL1", "SIN3B")

# Read the CSV file with annotated results
find_id <- read.csv("annotated_7_LFC_Monocyte_Mice_FROM_36_Samples.csv")

# Check the first few rows to ensure the data is loaded correctly
head(find_id)

# Find matching rows where `external_gene_name` matches any of the gene_ids
matched_genes <- find_id[find_id$external_gene_name %in% gene_name, ]

# Extract corresponding `Given_id` for those matched genes
given_ids <- matched_genes$Given_id
# Remove the second value
given_ids <- given_ids[-2]

data <- read.csv("2_count_Data_Normalized_124_Imputed_MICE.csv")

data_ids <- data[data$X %in% given_ids,]

head(data_ids)

# Replace `gene_id` with the corresponding gene name from `matched_genes`
data_ids$X <- matched_genes$external_gene_name[match(data_ids$X, matched_genes$Given_id)]
head(data_ids)


# Save the processed data to a CSV file
write.csv(data_ids, file = "Processed_2_count_Data_Normalized_124_Imputed_MICE.csv", row.names = FALSE)
