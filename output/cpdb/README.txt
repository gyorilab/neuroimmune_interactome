This folder contains files representing the interactome as follows:
- ligand_receptor_unique_interactions.csv: ligand-receptor interactions (using UniProt IDs)
- ligand_receptor_unique_interactions_hgnc.csv: ligand-receptor interactions (using gene names)
- ligand_ion_channel_unique_interactions.csv: ligand-ion channel interactions (using UniProt IDs)
- ligand_ion_channel_unique_interactions_hgnc.csv: ligand-ion channel interactions (using gene names)
- enzyme_product_unique_interactions.csv: enzyme-receptor/ion channel interactions through enzyme products (using UniProt IDs)
- enzyme_product_unique_interactions_hgnc.csv: enzyme-receptor/ion channel interactions through enzyme products (using gene names)
- all_interactions.csv: all interactions merged (using UniProt IDs)
- all_interactions_hgnc.csv: all interactions merged (using gene names)

All files using UniProt IDs are in a format compatible with CellPhoneDB.
The input file used for CellPhoneDB database construction and analysis was
all_interactions.csv.
