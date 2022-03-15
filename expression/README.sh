#!/bin/bash

# get expression data from : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75748

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75748/suppl/GSE75748_bulk_cell_type_ec.csv.gz
gunzip GSE75748_bulk_cell_type_ec.csv.gz

