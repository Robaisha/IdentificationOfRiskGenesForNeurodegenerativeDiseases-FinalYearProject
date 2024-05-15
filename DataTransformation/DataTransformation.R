if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("tidyverse")
BiocManager::install("preprocessCore")
BiocManager::install("openxlsx")
BiocManager::install("dplyr")

library(affy)
library(tidyverse)
library(openxlsx)
library(dplyr)

setwd("D:/FYP")
gzipped_txt_file_path <- "C:/Users/Masoo/Desktop/FYP_Material/GSE118553_non-normalized_data.txt.gz"
con <- gzfile(gzipped_txt_file_path, "rt")
raw_data <- read.table(con, header = TRUE, sep = "\t", comment.char = "", quote = "")
numeric_columns <- sapply(raw_data, is.numeric)
if (any(!numeric_columns)) {
  print("Warning: Non-numeric columns were detected and not included in normalization.")
}
numeric_data <- raw_data[, numeric_columns]
normalized_numeric_data <- normalize.quantiles(as.matrix(numeric_data))
normalized_data <- cbind(ID_REF = raw_data[, "ID_REF"], as.data.frame(normalized_numeric_data))
colnames(normalized_data) <- colnames(raw_data)
col_names <- colnames(normalized_data)


# Transpose the data frame
df_transposed <- t(normalized_data)

# Convert the transposed matrix back to a data frame
df_transposed_normalized_data <- as.data.frame(df_transposed)

file_path_excel <- "C:/Users/Masoo/Desktop/FYP_Material/BrainPart_GENE_IDS.xlsx"
excel_data_row <- readxl::read_excel(file_path_excel)
df_brainPartGeneIDs_excel <- data.frame(excel_data_row)

# Select columns 1, 2, and 3 from df_brainPartGeneIDs_excel Dataframe 
selected_columns_brainPartGeneIDs_excel <- df_brainPartGeneIDs_excel[, c(1, 2, 3)]

# Specify new column names Brain Parts Dataframe 
selected_column_names_brainPartGeneIDs_excel <- c("Brain_Part", "Type", "GENE_ID")

# Create a new data frame Brain Parts Dataframe 
df_SelectedColumns_brainPartGeneID_excel <- data.frame(selected_columns_brainPartGeneIDs_excel)

# Set column names for the selected columns Brain Parts Dataframe 
colnames(df_SelectedColumns_brainPartGeneID_excel) <- selected_column_names_brainPartGeneIDs_excel

#Select excel file of transposed data without ILMN Gene
file_path_excel <- "C:/Users/Masoo/Desktop/FYP_Material/df_normalized_data.xlsx"
excel_data_row <- readxl::read_excel(file_path_excel)
df_Excel_Transposed_NormalizedData<- data.frame(excel_data_row)

#Checking col names
colnames(df_Excel_Transposed_NormalizedData)
colnames(df_SelectedColumns_brainPartGeneID_excel)

#Checking column types
str(df_Excel_Transposed_NormalizedData)
str(df_SelectedColumns_brainPartGeneID_excel)

#Perform Join on Brain Parts
df_join_brainPartGeneID_transposed_normalized_data <- merge(x = df_SelectedColumns_brainPartGeneID_excel, y =df_Excel_Transposed_NormalizedData, by = "GENE_ID", all.x = TRUE)
write.csv(df_join_brainPartGeneID_transposed_normalized_data, "C:/Users/Masoo/Desktop/FYP_Material/df_join_brainPartGeneID_transposed_normalized_data.txt", row.names=FALSE)

#s_no<- 1:47323
#normalized_data <- cbind(s_no, normalized_data )
#print (normalized_data)
#df_normalized_data <- normalized_data[,-1 ]
#df_normalized_data <-t(df_normalized_data)
#selected_columns_normalized <- df_normalized_data[,c(0) ]
#df_rownames_transpose_data <- data.frame(rownames(df_normalized_data))
#df_rownames_transpose_data <- df_rownames_transpose_data[-c(1,2), ]


#Reading from files
excel_file_path <- "C:/Users/Masoo/Desktop/FYP_Material/Book1.xlsx"
# Read data from Excel sheet
excel_data <- readxl::read_excel(excel_file_path)
df_Book1_excel <- data.frame(excel_data)

df_normalized <- data.frame(normalized_data)

#removing detection_p value
normalized_data_without_col <- normalized_data
cols_to_remove <- grep("Detection_Pval", colnames(normalized_data_without_col))

# Remove the columns using column indexing
normalized_data_without_col <- normalized_data_without_col[, -cols_to_remove]

colnames(df_normalized)[1] ="ID_REF_Probe_ID"
# Select columns 4, 5, and 14 from df_Book1_excel Dataframe 1
selected_columns_Book1_excel <- df_Book1_excel[, c(4, 5, 14)]

# Specify new column names Dataframe 1
selected_column_names_Book1_excel <- c("ID_REF", "ILMN_GENE", "ID_REF_Probe_ID")

# Create a new data frame Dataframe 1
df_SelectedColumns_Book1_excel <- data.frame(selected_columns_Book1_excel)

# Set column names for the selected columns Dataframe 1
colnames(df_SelectedColumns_Book1_excel) <- selected_column_names_Book1_excel

#Removing row 1 containing names of col from excel
df_SelectedColumns_Book1_excel <- df_SelectedColumns_Book1_excel[-1, ]

# Select columns 1 from df_Book1_excel Dataframe 2
selected_columns_normalized <- df_normalized[, c(1)]

# Specify new column names Dataframe 2
selected_column_names_normalized <- c("ID_REF_Probe_ID")

# Create a new data frame Dataframe 2
df_SelectedColumns_normalized <- data.frame(selected_columns_normalized)

# Set column names for the selected columns Dataframe 2
colnames(df_SelectedColumns_normalized) <- selected_column_names_normalized

#Perform Join
df_join <- merge(x = df_SelectedColumns_Book1_excel, y = df_SelectedColumns_normalized, by = "ID_REF_Probe_ID", all.x = TRUE)

#Checking col names
colnames(df_SelectedColumns_Book1_excel)
colnames(df_SelectedColumns_normalized)

#Checking column types
str(df_SelectedColumns_Book1_excel)
str(df_SelectedColumns_normalized)


#Join on normalized data
df_Join_NormalizedData_df_join <- merge(x = df_join , y = df_normalized, by = "ID_REF_Probe_ID", all.x = TRUE)

#removing detection_p value
df_Join_NormalizedData_df_join_copy <- df_Join_NormalizedData_df_join
cols_to_remove <- grep("Detection_Pval", colnames(df_Join_NormalizedData_df_join_copy))

# Remove the columns using column indexing
df_Join_NormalizedData_df_join_copy <- df_Join_NormalizedData_df_join_copy[, -cols_to_remove]
df_Join_NormalizedData_df_join_copy<-df_Join_NormalizedData_df_join_copy[,-2]

#convert this dataframe into excel sheet
write.xlsx(df_Join_NormalizedData_df_join, "C:/Users/Masoo/Desktop/FYP_Material/df_Join_NormalizedData_df_join.xlsx")


# Select column 3 ILMN GENE For Trrust database 1
selected_columns_Trust_df1 <- df_Join_NormalizedData_df_join[, c(3)]

# Specify column name ILMN GENE For Trrust database 1
selected_columns_name_Trust_df1 <- c("ILMN_GENE")

# Create a new data frame ILMN_GENE For Trust database 1
df_Normalized_Trust_df1 <- data.frame(selected_columns_Trust_df1)

# Set column names for df_ILMN_GENE_Trrust_df1
colnames(df_Normalized_Trust_df1) <- selected_columns_name_Trust_df1

#Specify file path for Trrust database
file_path <- "C:/Users/Masoo/Desktop/FYP_Material/trrust_rawdata.human.tsv"

#Read file
Trrust_database<- read.delim(file_path)

#Create dataframe
df_Trrust_database<- data.frame(Trrust_database)

# Select columns 1,2 For Trrust database 2
selected_columns_Trust_df2 <- df_Trrust_database[, c(1,2)]

# Specify column name ILMN GENE,ALIAS For Trrust database 2
selected_columns_name_Trust_df2 <- c("ILMN_GENE","TF")

# Create a new data frame df_Trrust_df2 For Trust database 2
df_Trust_df2 <- data.frame(selected_columns_Trust_df2)

# Set column names for df_Trrust_df2
colnames(df_Trust_df2) <- selected_columns_name_Trust_df2

#Join on Trrust database
df_Join_Trrust_df1_df2<- merge(x = df_Trust_df2 , y = df_Normalized_Trust_df1, by = "ILMN_GENE", all.x = TRUE)

df_Join_NormalizedData_Trrust_resultant <- merge(x =df_Join_Trrust_df1_df2 , y = df_Join_NormalizedData_df_join , by = "ILMN_GENE", all.x = TRUE)
df_Join_NormalizedData_Trrust_resultant_copy<-df_Join_NormalizedData_Trrust_resultant
#removing detection_p value
cols_to_remove <- grep("Detection_Pval", colnames(df_Join_NormalizedData_Trrust_resultant_copy))

# Remove the columns using column indexing
df_Join_NormalizedData_Trrust_resultant_copy <- df_Join_NormalizedData_Trrust_resultant_copy[, -cols_to_remove]
df_Join_NormalizedData_Trrust_resultant_copy<-df_Join_NormalizedData_df_join_copy[,-4]
write.csv(df_Join_NormalizedData_Trrust_resultant, "C:/Users/Masoo/Desktop/FYP_Material/df_Join_Trrust_df1_df2.txt", row.names=FALSE)


