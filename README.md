# Neuro-immune interactome

This repository accompanies the manuscript "Nociceptor neuroimmune interactomes
reveal common and injury-specific inflammatory pain pathways" and contains
code to generate the interactome.

The Python script to generate the interactome is at `interactome/core.py`.
Output files, including gene lists and generated CellPhoneDB input interaction tables
are in the `output` folder. In the `output/cpdb` folder, interactome files
using UniProt IDs are provided (`enzyme_product_interactions.csv`,
`ligand_ion_channel_interactions.csv`, `ligand_receptor_interactions.csv`)
which are used as input to CellPhoneDB. In addition, versions of these
interactome files are using HGNC gene names (`enzyme_product_interactions_hgnc.csv`,
`ligand_ion_channel_interactions_hgnc.csv`, `ligand_receptor_interactions_hgnc.csv`)
are also provided for convenience.
