# CONIPHER

## CONIPHER mutation clustering and phylogenetic tree building R package

This is the official github repository for the R package to perform mutation clustering and phylogenetic tree building using CONIPHER. For details on how to run mutation clustering and phylogenetic tree builing consecutively with one wrapper script from the command line, please refer to the github repository [CONIPHER-wrapper](https://github.com/McGranahanLab/CONIPHER-wrapper). For full details of all the inputs into CONIPHER clustering and tree building, refer to our protocol (https://doi.org/10.21203/rs.3.pex-2158/v1).

### Software
The current implementation of CONIPHER is written in `R>=3.6.1` and is distributed as an R package.

---
## CONIPHER installation 


CONIPHER can be installed and run in multiple modes.

1) To run mutation clustering + phylogenetic tree building, install the `conipher` conda environment from bioconda using the command below. This environment contains R version 4.1.3, the CONIPHER R package, PyClone v0.13.1 (which is required for mutation clustering), and all other dependencies.

```
conda create -n conipher -c conda-forge -c bioconda conipher
```

2) To run phylogenetic tree building only, the conda environment can also be used but is not strictly necessary. Alternatively, just the CONIPHER R package can be installed, provided the R package `devtools >= 2.4.1` is installed.

To install the CONIPHER R package from an R console, run the following command:

```
library(devtools)
devtools::install_github("McGranahanLab/CONIPHER")
```

---
## Quick start

### Running clustering + tree building end-to-end 
To get start quickly, you can install CONIPHER and perform mutation clustering and phylogenetic tree reconstruction on the example data provided using the following instructions.


**Step 1.** 
Install the `conipher` conda environment using the instructions above.

**Step 2.**
Start R and load the 'CONIPHER' and 'tidyverse' R packages using the following command:
```
library(CONIPHER)
library(tidyverse)
```

**Step 3.**
Specify a parent output directory where the clustering and tree building results will be saved in individual subfolders, for example using the following command:
```
out_dir <- "conipher_results/"
```

**Step 4.**
Specify the location of the input table .tsv file. For example, the file path of the toy input table provided in this package is specified using the following command:
```
input_tsv_loc <- system.file("extdata", "input_table.tsv", package = "CONIPHER", mustWork = TRUE)
```

**Step 5.**
Run clustering + tree building end-to-end (interactively) using the following command:
```
conipher_run(case_id = "CRUKTOY001",
             prefix = "CRUK",
             out_dir = out_dir,
             input_tsv_loc = input_tsv_loc)
```

### Running clustering only
Run steps 1 - 4 as described in "Running clustering + tree building end-to-end" above. 

**Step 5a.**
Run clustering (interactively) using the following command:

```
conipher_clustering(case_id = "CRUKTOY001", 
                    out_dir = out_dir, 
                    input_tsv_loc = input_tsv_loc)
```

### Running tree building only
Run steps 1 - 4 as described in "Running clustering + tree building end-to-end" above. 

**Step 1b.** 
Alternatively, instead of installing the `conipher` conda environment, install the CONIPHER R package only using the instructions described above.


**Step 5b.**
Run tree building (interactively) using the following command:

```
conipher_treebuilding(prefix = "CRUK",
                      out_dir = out_dir,
                      input_tsv_loc = tree_input_tsv_loc)
```


---
### Anticipated results
The clustering output will include the following output files (examples are in "conipher_results/Clustering"):
- <CASE_ID>.SCoutput.CLEAN.tsv: A mutation table in the same format as input.tsv, including columns for the
following: mutation cluster assignments (CLUSTER), mutation PhyloCCF estimates (CCF_PHYLO); observed CCF (CCF_OBS); mutation copy number estimates for each mutation in each tumor sample (MUT_COPY); mutation identifier <CASE_ID>:<CHR>:<POS>:<REF>:<ALT> (mutation_id)
- <CASE_ID>.SCoutput.FULL.tsv: A mutation table in same format as above except with additional column CLEAN, a logical flag indicating whether this mutation was removed (CLEAN==FALSE)
- <CASE_ID>.removed.muts.txt: A mutation table comprising mutations removed during the clustering stage. This file will not be generated if no mutations are removed.
- <CASE_ID>_pyclone_cluster_assignment_copynumber_clean.pdf: Genome plots of each mutation plotted at its genomic position (x-axis) vs it's estimated mutation copy number (y-axis) for each tumor sample

The clustering output will include the following output plots:
- <CASE_ID>.removedCPN.muts.pdf: Similar genome plots of mutation copy numbers plotted for mutations that were removed during mutation clustering post-processing
- <CASE_ID>.pyclone_cluster_assignment_phylo_clean.pdf: Scatter plot of estimated PhyloCCF of each (non-removed) mutation between each pair of samples, coloured by assigned mutation cluster

The tree building output will include the following output files (examples are in "conipher_results/Trees"):
- allTrees.txt: A text file containing all inferred phylogenetic trees
- alternativeTreeMetrics.txt: A tab-delimted file containing summary metrics of all alternative phylogenetic trees
- clusterInfo.txt: A tab-delimted text file containing information about each inferred mutation cluster
- cloneProportionsMinErrorTrees.txt: A tab-delimted file containing inferred clone proportion tables for the phylogenetic trees that had the lowest error
- subclonalExpansionScoreMinErrorTrees.txt: A tab-delimted file containing subclonal expansion scores for the phylogenetic trees that had the lowest error
- consensusBranches.txt: A text file containing all ancestor cluster - descendant cluster pairs that were shared amongst all inferred alternative phylogenetic trees
- consensusRelationships.txt: A text file containing all parent cluster - child cluster edge pairs that were shared amongst all inferred alternative phylogenetic trees
- treeTable.tsv: A tab-separated mutation table in the format of input.tsv, including columns originalCLUSTER (cluster ID that the mutation was assigned to in CONIPHER clustering stage) and treeCLUSTER (the mutation's final cluster assignment after CONIPHER tree building stage)
- <CASE_ID>.tree.RDS: an R list object containing all tree building output information

The tree building output will include the following output plots:
- pytree_and_bar.pdf: a plot of the default reconstructed tree and barplot
- pytree_multipletrees.pdf: a plot showing all possible alternative phylogenetic trees found by CONIPHER



