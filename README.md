# CONIPHER

## CONIPHER phylogenetic tree building 
#### R package

This is the official github repository for the R package to perform tumour phylogenetic tree building using CONIPHER. For full details of the method, please refer to our manuscript and protocol XXX. For details on how to run mutation clustering and phylogenetic tree builing consecutively with one wrapper script, please refer to the github repository at https://github.com/FrancisCrickInstitute/TRACERx-clustering-treeBuilding. 

---
#### CONIPHER installation and run instructions
**Step 1.** To install and be able to run CONIPHER tree building, start R and install the 'CONIPHER' R package from this github repository using the following command:

```
devtools::install_github("McGranahanLab/CONIPHER")
library(CONIPHER)
```
(NOTE: Currently github repo is private so use `devtools::install_github("McGranahanLab/CONIPHER", auth_token='YOUR_AUTH_TOKEN_HERE'))`)


**Step 2.** Specify your input data. To run CONIPHER treebuilding successfully, we require an input table, a sample and tumour case prefix and a desired output directory. For example:
```
input_table <- read.delim2('path_to_your_input_table_here.tsv', stringsAsFactors = FALSE)
prefix <- 'CRUK'
out_dir <- 'path_to_your_output_dir_here'
```

**Step 3.** Run data preprocessing:
```
input_list <- treebuilding_preprocess(input_table, prefix, out_dir)
```

**Step 4.** Run the main tree building function:
```
tree_output_list <- treebuilding_run(input_list)
```
Full details of all input parameters in `treebuilding_run()` can be found in our protocol [here](https://www.example.com)


**Step 5.** Plot tree output, if desired:
```
treebuilding_plot(tree_output_list)
```


---
(might need to load in conda env - check this)
```
conda create -n pyclone -c bioconda -c conda-forge pyclone r-base=3.6.1 r-essentials r-tidyverse r-cowplot r-ggpubr r-fst r-biocmanager r-devtools r-seqminer
```
