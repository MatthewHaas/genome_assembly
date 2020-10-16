# Python code for merging Zizania palustris-specific duplications with their genomic positions
# This was intially run on my (Matthew's) MacBook using Jupyter Notebook.

# Load required packages
import pandas as pd

# Read in data
zpalustris_bed = pd.read_table('rice.gene_structures_post_PASA_updates.21917.bed', delimiter = '\t', header = None)

# Keep only the first four columns (these have the info we need)
zpalustris_bed = zpalustris_bed.iloc[:,0:4]

# Split the third column by semicolon so that it is split according to gene names
gene_names_expanded = zpalustris_bed[3].str.split(";", expand = True)

# Keep only the column with standard format gene names
gene_names = gene_names_expanded[[2]]

# Rename column names so that they can be merged
gene_names = gene_names.rename(columns={2: "gene_name"})

# Merge position information with gene names
gene_positions = pd.concat([zpalustris_bed.iloc[:,0:3], gene_names], axis=1)

# Assign names to first three columns
gene_positions = gene_positions.rename(columns={0: "scaffold"})
gene_positions = gene_positions.rename(columns={1: "start"})
gene_positions = gene_positions.rename(columns={2: "end"})

# Define the major scaffolds of interest and filter the table so that we only keep these.
scaffolds_of_interest = ["Scaffold_1", "Scaffold_3", "Scaffold_7", "Scaffold_9", "Scaffold_13", "Scaffold_18", "Scaffold_48", "Scaffold_51", "Scaffold_70", "Scaffold_93", "Scaffold_415", "Scaffold_693", "Scaffold_1062", "Scaffold_1063", "Scaffold_1064", "Scaffold_1065"]

gene_positions_filtered = gene_positions[gene_positions['scaffold'].isin(scaffolds_of_interest)]

# Now sort the dataframe based on the start position of each gene (as well as by scaffold, but it is already sorted that way)
gene_pos_filt_sort = gene_positions_filtered.sort_values(by=['scaffold', 'start'])

# Remove rows with duplicate values
gene_pos_filt_sort = gene_pos_filt_sort.drop_duplicates()

# Read in duplicated genes from the file Duplications.tsv which was parsed with 'parse_dulications_file'
duplications = pd.read_csv('zizania_duplications_simplified.csv', delimiter=',')

# Remove the string 'Zizania_palustris_' (positions 0-18)from each cell in columns 'gene_1' and 'gene_2'
# Note: only run this once because it replaces by position values, not by the string itself
duplications['gene_1'] = duplications['gene_1'].str.slice_replace(0, 18, '')
duplications['gene_2'] = duplications['gene_2'].str.slice_replace(0, 18, '')

# Rename the 'gene_1' column in 'duplications' to allow for merging with 'gene_pos_filt_sort'
duplications = duplications.rename(columns={'gene_1': 'gene_name'})

# Merge
duplications = gene_pos_filt_sort.merge(duplications, on = 'gene_name')

# Rename columns in preparation for merging again based on 'gene_2'
duplications = duplications.rename(columns={'scaffold':'scaffold_1', 'start':'start_1', 'end':'end_1', 'gene_name':'gene_1', 'gene_2':'gene_name'})

# Merge again, based on the duplicated gene column (formerly 'gene_2')
duplications = gene_pos_filt_sort.merge(duplications, on = 'gene_name')

# Rename columns to make it clear which are the positions based on 'gene_2'
duplications = duplications.rename(columns={'scaffold':'scaffold_2', 'start':'start_2', 'end':'end_2', 'gene_name':'gene_2'})

# Reorder dataframe (also, remove unnnamed column; those values correspond to row numbers in original Duplications.tsv file)
duplications = duplications[['gene_1', 'scaffold_1', 'start_1', 'end_1', 'gene_2', 'scaffold_2', 'start_2', 'end_2']]

# Write to csv file
duplications.to_csv('duplications_with_positions.csv')