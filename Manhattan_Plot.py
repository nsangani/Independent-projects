#!/usr/bin/env python
# coding: utf-8

# In[1]:
# Import necessary libraries for data manipulation and visualization
from bioinfokit import analys, visuz  # bioinfokit is used for analysis and plotting
import pandas as pd  # pandas for data handling

# In[2]:
# Read the data from the CSV file. Ensure the file 'A_H_genes_gtf.csv' exists in your working directory
df = pd.read_csv('A_H_genes_gtf.csv')

# In[3]:
# Data preprocessing: Remove the 'chr' prefix from chromosome names if present
# This step ensures that 'chr' is not part of the chromosome number when plotting
df['chr'] = df['chr'].str.replace('chr', "")  # Replacing 'chr' prefix
df.head()  # Display the first few rows of the dataframe to confirm changes

# In[4]:
# Visualizing the Manhattan plot using bioinfokit's visuz.marker.mhat method
# Parameters:
# - df: DataFrame containing the input data
# - chr: Column containing chromosome information (numeric or string)
# - pv: Column containing p-values (used to calculate -log10(p-value))
# - show: Whether to display the plot (True = show plot)
# - gwas_sign_line: Display a significance line based on GWAS p-value threshold (typically 5e-8)
# - markeridcol: Column containing gene or marker names for labeling
# - gstyle: A style parameter for marker style (set to 2 here)
# - dim: Dimensions of the plot (in this case, width=10 and height=2)

visuz.marker.mhat(
    df=df,                    # DataFrame to be used
    chr='chr',                # Column name with chromosome information
    pv='pvalue',              # Column with p-values (will be used for -log10(p-value))
    show=True,                # Display the plot
    gwas_sign_line=True,      # Add a GWAS significance line at p-value threshold
    markeridcol='gene',       # Column with gene identifiers for the plot markers
    gstyle=2,                 # Style for the plot markers (visual style)
    dim=(10, 2)               # Set the dimensions of the plot (width, height)
)
