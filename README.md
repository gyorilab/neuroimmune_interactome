# Neuro-immune interactome

This repository accompanies the manuscript "Nociceptor neuroimmune interactomes
reveal common and injury-specific inflammatory pain pathways" and contains
code to generate the interactome.

The Python script to generate the interactome is at `interactome/core.py`.
Output files, including gene lists and generated CellPhoneDB input interaction tables
are in the `output` folder.

In the `output/cpdb` folder, interactome files
using UniProt IDs are provided:
- `ligand_receptor_unique_interactions.csv`
- `ligand_ion_channel_interactions.csv`
- `enzyme_product_unique_interactions.csv`

which are used as input to CellPhoneDB. In addition, versions of these
interactome files that use gene names 
- `enzyme_product_unique_interactions_hgnc.csv`
- `ligand_ion_channel_unique_interactions_hgnc.csv`
- `ligand_receptor_unique_interactions_hgnc.csv`

are also provided for convenience. Finally, the merged interactome is available
both using UniProt IDs and using gene names:
- `all_interactions.csv`
- `all_interactions_hgnc.csv`
