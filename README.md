# CONIPHER

## CONIPHER phylogenetic tree building 
### R package

This is the official github repository for the R package to perform tumour phylogenetic tree building using CONIPHER. For full details of the method, please refer to our manuscript and protocol XXX. For details on how to run mutation clustering and phylogenetic tree builing consecutively with one wrapper script, please refer to the github repository at https://github.com/FrancisCrickInstitute/TRACERx-clustering-treeBuilding. 

#### Software
The current implementation of CONIPHER tree building is written in `r=3.6.1` and is distributed as an R package.

---
### Quick start

To get start quickly, you can install CONIPHER tree building and perform phylogenetic tree reconstruction on the example data provided using the following instructions.

#### CONIPHER installation and run instructions
To install and be able to run CONIPHER tree building, you must have R package `devtools` installed. 

**Step 1.** Start R and install the 'CONIPHER' R package from this github repository using the following command:
```
devtools::install_github("McGranahanLab/CONIPHER")
library(CONIPHER)
```

**Step 2.** Specify your input data. To run CONIPHER treebuilding successfully, we require an input table, a sample and tumour case prefix and a desired output directory. For example:
```
data(input_table)
prefix <- 'CRUK'
out_dir <- 'data/results'
```

**Step 3.** Run data preprocessing:
```
input_list <- treebuilding_preprocess(input_table, prefix, out_dir)
```

**Step 4.** Run the main tree building function:
```
tree_output_list <- treebuilding_run(input_list)
```
Full details of all input parameters in `treebuilding_run()` can be found in our protocol [here](CONIPHER_arxiv)


**Step 5.** Plot tree output, if desired:
```
treebuilding_plot(tree_output_list)
```
--- 

### Anticipated results

The tree building will return 3 output files:
- <CASE_ID>.tree.RDS: an R list object containing tree building output information
- pytree_and_bar.pdf: a plot of the default reconstructed tree and barplot
![Pytree and bar](https://github.com/McGranahanLab/CONIPHER/files/10381097/pytree_and_bar.pdf)
- pytree_multipletrees.pdf: a plot showing all possible alternative phylogenetic trees found by CONIPHER
![Multiple trees](https://github.com/McGranahanLab/CONIPHER/files/10381106/pytree_multipletrees.pdf)



