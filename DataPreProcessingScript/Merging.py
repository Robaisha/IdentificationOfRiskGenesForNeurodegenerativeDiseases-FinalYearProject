import sys
import numpy as np
import pandas as pd
# Load gene expression data
gene_expression_df = pd.read_csv('output_matrix.csv')

# Load trait data, skip the first row as it contains column headers
trait_data_df = pd.read_csv('output_trait.csv')

# Merge datasets on common identifier (e.g., sample ID)
merged_df = pd.concat([trait_data_df,gene_expression_df], axis=0)

merged_df.to_csv('merged_data.csv', index=False)
