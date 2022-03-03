#!/bin/bash

## STEP 1: Get raw data

### https://data.4dnucleome.org/experiment-set-replicates/4DNESEVQ1V15/
### Micro-C on H1 embryonic stem cells from Maehr lab with 1% Formaldehyde and 3mM DSG crosslinking

wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/24d35a88-5b24-478a-ac0d-d06a1ed52c40/4DNFIJLK5WML.mcool


### https://data.4dnucleome.org/experiment-set-replicates/4DNESP4MARXG/
### Micro-C on definitive endoderm derived from H1 embryonic stem cells (Maehr lab protocol) with 1% Formaldehyde and 3mM DSG crosslinking

wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/0a1b3d24-d714-4687-9a82-1b232bda9736/4DNFI9E222YJ.mcool


## STEP 2: Call loops with [peakachu](https://github.com/tariks/peakachu)

#### time = 2.5 minutes  ||  depth = 173,850,883
peakachu depth -p 4DNFIJLK5WML.mcool::/resolutions/5000
#### time = 2.5 minutes  ||  depth = 154,174,393
peakachu depth -p 4DNFI9E222YJ.mcool::/resolutions/5000
#### Using the 150 million, 5kb loop set (10% 5kb)
wget https://dl.dropboxusercontent.com/s/i2axj1ij4vbhcha/high-confidence.200million.5kb.pkl?dl=0

./run_peakachu.sh = call loops

ln -s 4DNFI9E222YJ_scores_1kb.loops.final endoderm.loops
ln -s 4DNFIJLK5WML_scores_1kb.loops.final h1.loops


## STEP 3: Get chromHMM calls from [epimap](http://compbio.mit.edu/epimap/)

# endoderm
wget https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/BSS00285_18_CALLS_segments.bed.gz
# H1-hESC
wget https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/BSS00478_18_CALLS_segments.bed.gz


## STEP 4: Get scores for loops

`get_scores_mcool.py`




## STEP 5: Get contact labels associated with loops

`get_contact_labels`
`./summarize.sh`



## STEP 5: Create matrix

`create_matrix.sh`

## STEP 6: run NMF

`RcppML_nmf/run_nmf.R`



