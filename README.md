# Tripartition

## Dependencies
- Python 2.7.15
  -DendroPy
  -Scipy
  -Numpy
- RAxML 8.2.11


## Overview

Given a gene family multiple sequence alignment, partition the alignment in three approximately equal sized alignments, create maximum likelihood trees, create bootstrap trees, and conduct a statistical analysis for presence/absence of sub-gene/partial gene transfer event.

1. 01_splitFasta.py - Tripartition the multiple sequence alignment:  

   command: `python 01_splitFasta.py`  
 
2. 02_raxml.sh - Generate maximum likelihood trees with RAxML:  
   
   command: `bash 02_raxml.sh /path/to/tripartitioned/alignment/directory/ model`

3. 03_raxml_bs.sh - Generate bootstrap tree (not rapid bootstraps) with RAxML:

   command: `bash 03_raxml_bs.sh /path/to/tripartitioned/alignment/directory/ model `

4. 04_hist_intersection.py - Conducts the statistical analysis for determining whether gene family is contains sub-gene transfer events:
   
   command: `python 04_hist_intersection.py --thresh 0.5 --bs_sample /path/to/tripartitioned/bootstrap/trees/ --ml_sample /path/to/tripartitioned/maximum/likelihood/trees/`

## Data Set Availability

Simulated data sets and the scripts used to generate such data are available in the sim_data directory.