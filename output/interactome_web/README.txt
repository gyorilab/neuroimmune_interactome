This file contains three HTML pages that can be opened locally
in a web browser. These pages provide Statements and their supporting evidence
as assembled by INDRA, constituting the interactome.
- ligand_receptor_statements.html.zip: Statements representing ligand/receptor
  interactions. Note that this file is very large and can take a while
  to load in a browser, it is available as a zipped file due to size.
- ligand_ion_channel_statements.html: Statements representing ligand/ion channel
  interactions.
- enzyme_product_statements.html: Statements representing enzyme product
  interactions with receptors or ion channels.

Note that enzyme_product_statements.html only contains enzyme product, not
enzymes themselves. The relationships between enzymes and the products they
control the production of is provided separately in the
tab-separated spreadsheet enzyme_product_relations.tsv. This is derived
directly from https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.Detailed.hgnc.txt.gz
and contains the following columns:
- enzyme: the official gene symbol of the enzyme
- product_id: the CHEBI ID of the product
- product_name: the standard CHEBI name of the product
- pubmed_ids: the list of PubMed IDs (separated by semicolons) supporting the
              relationship
