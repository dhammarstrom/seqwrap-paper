# seqwrap-paper


## Introduction

This repository contains details of the relevant files, data and code
needed to replicate the results and figures presented in the
[seqwrap](https://github.com/trainome/seqwrap) paper.

## Relevant folders

- `/R` This directory contains the scripts used to produce the results
  present in the paper. Below are short descriptions of what each file
  does

  1.  `data-prep.R` . If the libraries in this file are installed,
      running this script will download the raw gene counts from the
      [Pillion study](https://pubmed.ncbi.nlm.nih.gov/36070371/) , use
      seqwrap to build and save a preliminary model .

  2.  `simulation-functions.R` Contains R functions used at different
      points of the analysis such as the function to simulate data sets
      and summary functions. The other files (except data-prep.R) loads
      and executes this file to make the functions and variables in the
      script available in the R environment

  iii\. `m1-m5-pillion-data.R` This file builds ands save the models
  based on real world data from the [Pillion
  study](https://pubmed.ncbi.nlm.nih.gov/36070371/)

  iv\. `simulations.R`

  v\. `simulations2.R`

## Recommended steps to reproduce analysis
