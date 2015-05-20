#!usr/bin/bash

##ENSEMBL COMPARA
#Use perl lage_pantaxa.pl

wget ftp://ftp.ensemblgenomes.org/pub/pan_ensembl/release-22/emf/ensembl-compara/homologies/Compara.phyloxml_aa_trees.22.tar.gz
wget ftp://ftp.ensemblgenomes.org/pub/pan_ensembl/release-22/mysql/ensembl_compara_pan_homology_22_75/member.txt.gz
wget ftp://ftp.ensemblgenomes.org/pub/pan_ensembl/release-22/mysql/ensembl_compara_pan_homology_22_75/gene_tree_node.txt.gz
wget ftp://ftp.ensemblgenomes.org/pub/pan_ensembl/release-22/mysql/ensembl_compara_pan_homology_22_75/gene_tree_root.txt.gz
wget ftp://ftp.ensemblgenomes.org/pub/pan_ensembl/release-22/mysql/ensembl_compara_pan_homology_22_75/gene_tree_node_attr.txt.gz

##HOGENOM
#Use perl time_species_code.pl species_code
#then perl age_hogenom.pl
#then perl lage_hogenom.pl hogenom_aged.txt acc_id5_id6

wget ftp://pbil.univ-lyon1.fr/pub/hogenom/release_06/acc_id5_id6
wget ftp://pbil.univ-lyon1.fr/pub/hogenom/release_06/FAM_HOGENOMID.gz
wget ftp://pbil.univ-lyon1.fr/pub/hogenom/release_06/species_code
wget ftp://pbil.univ-lyon1.fr/pub/hogenom/release_06/hogenom6.phyml

##HOMOLOGENE
#Use perl time_taxid_taxname.pl then
#perl age_homologene.pl then
#perl lage_homologene.pl

wget ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data
wget ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/build_inputs/taxid_taxname
