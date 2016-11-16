This is a collection of scripts suitable for estimating the evolutionary age of genes. (and recreating my honors thesis)

Download the data, then run the appropriate scripts to perform the desired analyses

Downloading Data
================
`bash download_3rdparty_data.sh`

Age Genes
=========

Ensembl Compara
---------------
`perl lage_pantaxa.pl Compara.phyloxml_aa_trees.22 pantaxa_ages.txt`

HOGENOM
-------
`perl time_species_code.pl species_code`

`perl age_hogenom.pl`

`perl lage_hogenom.pl hogenom_aged.txt acc_id5_id6`

Homologene
----------
`perl time_taxid_taxname.pl taxid_taxname`

`perl age_homologene.pl`

`perl lage_homologene.pl homologene_aged.txt homologene.data.txt`
