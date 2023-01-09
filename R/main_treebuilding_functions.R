################################################################################
##################               MAIN FUNCTIONS               ##################
################################################################################

#' Input data preprocessing function
#'
#' This function takes the input tsv and formats the data to be compatible with
#' the main CONIPHER tree building function. NOTE: it is assumed that
#' clustering has been carried out prior to running tree building.
#' @param input_table An dataframe of the input mutation table in the correct
#' format. For more information on the input table format, please see our
#' tree building protocol.
#' @param prefix A tumour case and sample prefix, e.g. 'CRUK'.
#' @param out_dir A file path to the desired output directory
#' @export treebuilding_preprocess

treebuilding_preprocess <- function(input_table, prefix, out_dir)
{
  # check if the correct columns are included
  required_cols <- c("CASE_ID", "SAMPLE", "CHR", "POS", "REF", "ALT", "CLUSTER", "CCF_PHYLO", "CCF_OBS", "MUT_COPY", "COPY_NUMBER_A", "COPY_NUMBER_B")
  if (FALSE%in% (required_cols %in% colnames(input_table)))
  {
    print('\nThe following columns are required in input_tsv:\n')
    cat(required_cols)
    stop()
  }

  # add mutation id column
  input_table$mutation_id <- paste(input_table$CASE_ID,
                                   input_table$CHR,
                                   input_table$POS,
                                   input_table$REF,
                                   input_table$ALT,
                                   sep=":")

  nr_unique_muts           <- length(unique(input_table$mutation_id))
  nr_regions               <- length(unique(input_table$SAMPLE))
  regions                  <- unique(input_table$SAMPLE)

  # Raise an error if prefix is not specified, or incorrectly specified
  if (is.null(prefix)){
    stop("No prefix specified. Please indicate a prefix for the current tumour case.")
  } else if (!grepl(prefix, input_table$CASE_ID[1])){
    stop("Incorrect prefix specified. Please input the correct prefix for the current tumour case.")
  }

  if(!'CLUSTER'%in%colnames(input_table$CASE_ID[1]))
  {
    if('originalCLUSTER'%in%colnames(input_table))
    {
      colnames(input_table) <- gsub("originalCLUSTER","CLUSTER",colnames(input_table))
    }
  }

  # make sure all columns are the correct class
  input_table$POS <- as.numeric(input_table$POS)
  input_table$CCF_PHYLO <- as.numeric(input_table$CCF_PHYLO)
  input_table$CCF_OBS <- as.numeric(input_table$CCF_OBS)
  input_table$MUT_COPY <- as.numeric(input_table$MUT_COPY)
  input_table$COPY_NUMBER_A <- as.numeric(input_table$COPY_NUMBER_A)
  input_table$COPY_NUMBER_B <- as.numeric(input_table$COPY_NUMBER_B)
  input_table$REF_COUNT <- as.numeric(input_table$REF_COUNT)
  input_table$VAR_COUNT <- as.numeric(input_table$VAR_COUNT)
  input_table$DEPTH <- as.numeric(input_table$DEPTH)
  input_table$ACF <- as.numeric(input_table$ACF)
  input_table$PLOIDY <- as.numeric(input_table$PLOIDY)


  # make sure all mutations have a cluster assigned
  tmp <- input_table[!is.na(input_table$CLUSTER),]
  removed_mutations <- c()
  for (mutation_id in unique(input_table$mutation_id))
  {
    if(length(unique(tmp[tmp$mutation_id%in%mutation_id,]$CLUSTER))==0)
    {
      # warning
      removed_mutations <- c(removed_mutations,mutation_id)
      cat('\nwarning:')
      cat('', paste(mutation_id),'does not have a CLUSTER assigned, will remove')
      next;
    }

    input_table[input_table$mutation_id%in%mutation_id,]$CLUSTER <- unique(tmp[tmp$mutation_id%in%mutation_id,]$CLUSTER)
  }

  if(length(removed_mutations)>=1)
  {
    cat(paste('\nwarning: ',length(removed_mutations), ' mutations removed due to lack of cluster assignment',sep=""))
  }

  input_table <- input_table[!is.na(input_table$CLUSTER),,drop=FALSE]

  # check again:
  nr_unique_muts           <- length(unique(input_table$mutation_id))
  nr_regions               <- length(unique(input_table$SAMPLE))
  regions                  <- unique(input_table$SAMPLE)


  # Next convert the input_table into a sample_input_list
  input_list <- list()
  names_input_list <- c("pyclone",
                        "pyclone_absolute",
                        "sampleID",
                        "dirName",
                        "pycloneType",
                        "prefix",
                        "releaseDir",
                        "saveDir",
                        "generalSave",
                        "correct_subclonality",
                        "mutTable",
                        "merged_clusters")

  #create the pyclone table
  input_format           <- data.frame(matrix(data = NA,
                                              nrow = nr_unique_muts,
                                              ncol = 11*nr_regions+2),
                                       stringsAsFactors = FALSE)
  colnames(input_format)  <- c(paste(regions, "_cov", sep = "")
                               ,paste(regions, "_var_count", sep = "")
                               ,paste(regions, "_VAF", sep = "")
                               ,paste(regions,"_PhyloCCF",sep="")
                               ,paste(regions,"_PycloneCCF",sep="")
                               ,paste(regions,"_Pyclone_0.05",sep="")
                               ,paste(regions,"_Pyclone_0.95",sep="")
                               ,paste(regions,"_cpn.copies",sep="")
                               ,paste(regions,"_mut.cpn.num",sep="")
                               ,paste(regions,"_nAraw",sep="")
                               ,paste(regions,"_nBraw",sep="")
                               ,"PycloneCluster"
                               ,"CleanCluster")
  rownames(input_format) <- unique(input_table$mutation_id)
  # next populate the table
  for (mutation_id in rownames(input_format))
  {
    spec_mut_table <- input_table[input_table$mutation_id%in%mutation_id,,drop=FALSE]
    for (i in 1:nrow(spec_mut_table))
    {
      region_spec_mut <- spec_mut_table[i,,drop=FALSE]
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_cov',sep="")] <- region_spec_mut$DEPTH
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_var_count',sep="")] <- region_spec_mut$VAR_COUNT
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_VAF',sep="")] <- region_spec_mut$VAR_COUNT / region_spec_mut$DEPTH

      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_PhyloCCF',sep="")] <- region_spec_mut$CCF_PHYLO
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_PycloneCCF',sep="")] <- region_spec_mut$CCF_PHYLO
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_Pyclone_0.05',sep="")] <- region_spec_mut$CCF_PHYLO
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_Pyclone_0.95',sep="")] <- region_spec_mut$CCF_PHYLO
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_cpn.copies',sep="")] <- 1

      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_mut.cpn.num',sep="")] <- region_spec_mut$MUT_COPY
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_nAraw',sep="")] <- region_spec_mut$COPY_NUMBER_A
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_nBraw',sep="")] <- region_spec_mut$COPY_NUMBER_B
      input_format[mutation_id,"PycloneCluster"] <- region_spec_mut$CLUSTER
      input_format[mutation_id,"CleanCluster"]  <- 1
    }
  }

  # Do the same thing for pyclone_absolute (non-subclonal-copy-number-corrected version)
  #create the pyclone table
  input_format_absolute           <- data.frame(matrix(data = NA,
                                                       nrow = nr_unique_muts,
                                                       ncol = 11*nr_regions+2),
                                                stringsAsFactors = FALSE)
  colnames(input_format_absolute)  <- c(paste(regions, "_cov", sep = "")
                                        ,paste(regions, "_var_count", sep = "")
                                        ,paste(regions, "_VAF", sep = "")
                                        ,paste(regions,"_PhyloCCF",sep="")
                                        ,paste(regions,"_PycloneCCF",sep="")
                                        ,paste(regions,"_Pyclone_0.05",sep="")
                                        ,paste(regions,"_Pyclone_0.95",sep="")
                                        ,paste(regions,"_cpn.copies",sep="")
                                        ,paste(regions,"_mut.cpn.num",sep="")
                                        ,paste(regions,"_nAraw",sep="")
                                        ,paste(regions,"_nBraw",sep="")
                                        ,"PycloneCluster"
                                        ,"CleanCluster")
  rownames(input_format_absolute) <- unique(input_table$mutation_id)
  # next populate the table
  for (mutation_id in rownames(input_format_absolute))
  {
    spec_mut_table <- input_table[input_table$mutation_id%in%mutation_id,,drop=FALSE]

    for (i in 1:nrow(spec_mut_table))
    {
      region_spec_mut <- spec_mut_table[i,,drop=FALSE]
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_cov',sep="")] <- region_spec_mut$DEPTH
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_var_count',sep="")] <- region_spec_mut$VAR_COUNT
      input_format[mutation_id,paste(region_spec_mut$SAMPLE,'_VAF',sep="")] <- region_spec_mut$VAR_COUNT / region_spec_mut$DEPTH

      input_format_absolute[mutation_id,paste(region_spec_mut$SAMPLE,'_PhyloCCF',sep="")] <- region_spec_mut$CCF_OBS
      input_format_absolute[mutation_id,paste(region_spec_mut$SAMPLE,'_PycloneCCF',sep="")] <- region_spec_mut$CCF_OBS
      input_format_absolute[mutation_id,paste(region_spec_mut$SAMPLE,'_Pyclone_0.05',sep="")] <- region_spec_mut$CCF_OBS
      input_format_absolute[mutation_id,paste(region_spec_mut$SAMPLE,'_Pyclone_0.95',sep="")] <- region_spec_mut$CCF_OBS
      input_format_absolute[mutation_id,paste(region_spec_mut$SAMPLE,'_cpn.copies',sep="")] <- 1

      input_format_absolute[mutation_id,paste(region_spec_mut$SAMPLE,'_mut.cpn.num',sep="")] <- region_spec_mut$MUT_COPY
      input_format_absolute[mutation_id,paste(region_spec_mut$SAMPLE,'_nAraw',sep="")] <- region_spec_mut$COPY_NUMBER_A
      input_format_absolute[mutation_id,paste(region_spec_mut$SAMPLE,'_nBraw',sep="")] <- region_spec_mut$COPY_NUMBER_B
      input_format_absolute[mutation_id,"PycloneCluster"] <- region_spec_mut$CLUSTER
      input_format_absolute[mutation_id,"CleanCluster"]  <- 1

    }

  }

  # Now create list object for input to treebuilding
  input_list$pyclone          <- input_format
  input_list$pyclone_absolute <- input_format_absolute
  input_list$sampleID         <- input_table$CASE_ID[1]
  input_list$dirName          <- NA
  input_list$pycloneType      <- NA
  input_list$prefix           <- prefix
  input_list$releaseDir       <- NA
  input_list$saveDir          <- out_dir
  input_list$generalSave      <- out_dir
  input_list$correct_subclonality <- TRUE
  input_list$mutTable           <- NA
  input_list$merged_clusters   <- NA

  if(!file.exists(input_list$generalSave))
  {
    dir.create(input_list$generalSave,showWarnings = TRUE, recursive = TRUE, mode = "0775")
  }

  return(input_list)
}


#' TRACERx tree building function
#'
#' This function is the main CONIPHER wrapper function to run phylogenetic
#' tree building from mutation clustering output. NOTE: it is assumed that
#' clustering has been carried out prior to running tree building.
#' @param sample_input_list An R list object. This object contains information
#' about each mutation in each tumour region sampled, including which cluster
#' each mutation was assigned to in mutation clustering.
#' This object can be generated by running the data preprocessing function
#' @param ccf_buffer PhyloCCF buffer allowance for testing tree level issue
#' (default=10)
#' @param pval_cutoff P-value cut off for testing cluster nesting (default=0.01)
#' @param use_boot Should bootstrapping be used to compute confidence interval?
#' (default=TRUE)
#' @param merge_clusters Should similar clusters be merged if possible?
#' (default=TRUE)
#' @param correct_cpn_clusters Should clusters driven by copy number errors be
#' removed? (default=TRUE)
#' @param adjust_noisy_clusters Should noisy clusters be adjusted? (default=TRUE)
#' @param adjust_noisy_clusters_prop What is the minimum proportion of mutations
#' required to be present in a region to avoid cluster adjustment? (default=0.05)
#' @param min_ccf What is the minimum CCF threshold to consider a mutation as
#' present? (default=0.01)
#' @param min_cluster_size What is the minimum number of mutations required in a
#' cluster to be included in analysis? (default=5)
#' @param plotting Should inferred phylogenetic trees be plotted? (default=TRUE)
#' @param run.multi.trees Should alternative tumour phylogenies be explored?
#' (default=TRUE)
#' @returns sample_pyclone_tree, an R list object containing output information
#' from CONIPHER tree building
#' @export treebuilding_run

treebuilding_run <- function(sample_input_list
                                  , ccf_buffer = 10
                                  , pval_cutoff = 0.01
                                  , use_boot = TRUE
                                  , merge_clusters = TRUE
                                  , correct_cpn_clusters = TRUE
                                  , adjust_noisy_clusters = FALSE
                                  , adjust_noisy_clusters_prop = 0.05
                                  , min_ccf = 0.01
                                  , min_cluster_size = 5
                                  , plotting = TRUE
                                  , run.multi.trees = TRUE
)
{
  suppressPackageStartupMessages(require(igraph))
  suppressPackageStartupMessages(require(mapplots))

  cat('\n\nStarting Tree Building')
  #first keep track of parameters used for this
  input_parameter_list <- list()
  input_parameter_list$sampleID                      <- sampleID               <- sample_input_list$sampleID
  input_parameter_list$dirName                       <- dirName                <- sample_input_list$dirName
  input_parameter_list$pycloneType                   <- pycloneType            <- sample_input_list$pycloneType
  input_parameter_list$prefix                        <- prefix                 <- sample_input_list$prefix
  input_parameter_list$releaseDir                    <- releaseDir             <- sample_input_list$releaseDir
  input_parameter_list$saveDir                       <- saveDir                <- sample_input_list$saveDir
  input_parameter_list$generalSave                   <- generalSave            <- sample_input_list$generalSave
  input_parameter_list$correct_subclonality          <- correct_subclonality   <- sample_input_list$correct_subclonality
  input_parameter_list$mutTable                      <- mutTable               <- sample_input_list$mutTable
  input_parameter_list$ccf_buffer                    <- ccf_buffer
  input_parameter_list$pval_cutoff                   <- pval_cutoff
  input_parameter_list$use_boot                      <- use_boot
  input_parameter_list$merge_clusters                <- merge_clusters
  input_parameter_list$correct_cpn_clusters          <- correct_cpn_clusters
  input_parameter_list$adjust_noisy_clusters         <- adjust_noisy_clusters
  input_parameter_list$adjust_noisy_clusters_prop    <- adjust_noisy_clusters_prop
  input_parameter_list$min_ccf                       <- min_ccf
  input_parameter_list$min_cluster_size              <- min_cluster_size

  #
  cat('\nFollowing parameters used for tree building:\n')
  print(do.call(rbind,input_parameter_list))

  # prepare the output
  output_list  <- list()
  output_list$ccf_table_pyclone                      <- sample_input_list$pyclone
  output_list$ccf_table_absolute                     <- sample_input_list$pyclone_absolute
  output_list$ccf_table_pyclone_clean                <- sample_input_list$pyclone
  output_list$ccf_table_absolute_clean               <- sample_input_list$pyclone_absolute

  output_list$merged_clusters                        <- sample_input_list$merged_clusters

  output_list$noisy_clusters_adjusted                <- NA
  output_list$cpn_removed_clusters                   <- NA
  output_list$tree_removed_clusters                  <- NA

  test_pyclone          <- sample_input_list$pyclone
  test_pyclone_absolute <- sample_input_list$pyclone_absolute

  #check how many regions you have...

  #state parameters for use below.
  if(adjust_noisy_clusters)
  {
    cat('\nAdjusting noisy clusters\n')
    pyclone_adj  <- clean.noisy.clusters(pyclone = test_pyclone,max.absent.prop = adjust_noisy_clusters_prop)
    test_pyclone <- pyclone_adj$corrected_pyclone

    output_list$ccf_table_pyclone_clean <- test_pyclone

    if(!is.na(pyclone_adj$corrected_cluster)[1])
    {
      tmp <- pyclone_adj$corrected_cluster
      colnames(tmp) <- c('region','cluster')
      output_list$noisy_clusters_adjusted <- tmp
      cat('\nThe following clusters were adjusted:\n')
      print(tmp)
    }

    pyclone_adj_absolute  <- clean.noisy.clusters(test_pyclone_absolute,max.absent.prop = adjust_noisy_clusters_prop)
    test_pyclone_absolute <- pyclone_adj_absolute$corrected_pyclone

    output_list$ccf_table_absolute_clean <- test_pyclone_absolute
  }

  # make sure you only use clean clusters
  test_pyclone          <- test_pyclone[test_pyclone[, "CleanCluster"] %in% 1, -ncol(test_pyclone)]
  test_pyclone_absolute <- test_pyclone_absolute[test_pyclone_absolute[, "CleanCluster"] %in% 1, -ncol(test_pyclone_absolute)]

  if (nrow(test_pyclone) < min_cluster_size) { stop('too few mutations to run tree building') }
  if (sort(table(test_pyclone[, 'PycloneCluster']), decreasing = T)[1] < min_cluster_size) { stop('too few mutations to run tree building') }

  clusters_with_min_cluster_sizeations <- table(test_pyclone[,'PycloneCluster'])[table(test_pyclone[,'PycloneCluster'])>=min_cluster_size]


  cat('\n\n\nDetermining nesting of clusters\n')
  nested_pyclone <- determine.cluster.nesting(pyclone = test_pyclone
                                              , prefix = prefix
                                              , min_cluster_size = max(c(2, min_cluster_size))
                                              , pval_cutoff = pval_cutoff
                                              , use_boot =use_boot
                                              , min_ccf =  min_ccf
  )

  # NM additional step, taking clonality into account when nesting [05/04/2022]
  directedGraph_input_full <- matrix(0, 0, 2)
  colsums <- colSums(nested_pyclone$nestedclust)
  rowsums <- rowSums(nested_pyclone$nestedclust)
  trunk_cluster <- names(colsums[which(colsums == max(colsums))])
  if (length(trunk_cluster) > 1) {
    trunk_cluster <- names(sort(rowMeans(nested_pyclone$ccf_cluster_table[trunk_cluster,, drop = F]), decreasing = T))[1]
  }

  clonality_table <- clonality.function(pyclone = test_pyclone
                                        ,trunk =trunk_cluster
                                        ,prefix = prefix
                                        , min_cluster_size = max(c(2, min_cluster_size))
                                        ,pval_cutoff = pval_cutoff
                                        ,use_boot =use_boot )


  nested_pyclone  <- correct.clonality.nesting(nestedlist = nested_pyclone
                                               , pyclone = test_pyclone
                                               , clonality_table = clonality_table
                                               , pval_cutoff = pval_cutoff
                                               , min_cluster_size = min_cluster_size
                                               , min_ccf = min_ccf
                                               , prefix = prefix
  )
  # finish additional step NM [05/04/2022]
  cat('\nThe following nesting identified:\n')
  print(nested_pyclone$nestedclust[,order(colSums(nested_pyclone$nestedclust),decreasing=T)])
  cat('\n')

  nested_pyclone_absolute <- determine.cluster.nesting(pyclone = test_pyclone_absolute
                                                       , prefix = prefix
                                                       , min_cluster_size = max(c(2, min_cluster_size))
                                                       , pval_cutoff = pval_cutoff
                                                       , use_boot =use_boot
                                                       , min_ccf =  min_ccf
  )


  # we only merge clusters if we're not using absolute.
  if(merge_clusters %in% TRUE &
     nrow(nested_pyclone$nestedclust) > 1 &
     nrow(nested_pyclone_absolute$nestedclust) > 1)
  {
    cat('\nChecking for cluster merging\n')
    out <- merge.clusters.full(test_pyclone = test_pyclone
                               ,test_pyclone_absolute = test_pyclone_absolute
                               ,nested_pyclone = nested_pyclone
                               ,nested_pyclone_absolute = nested_pyclone_absolute
                               ,min_ccf = min_ccf
                               ,p_value_cut = pval_cutoff
                               ,min_cluster_size = min_cluster_size
                               ,use_boot = use_boot
    )

    nested_pyclone   <- out$nested_pyclone
    test_pyclone     <- out$test_pyclone

    if(!is.na(out$mergedclusters[1]) & !is.na(output_list$merged_clusters))
    {
      output_list$merged_clusters <- rbind(out$mergedclusters,output_list$merged_clusters)
      cat('\nThe following clusters were merged:\n')
      print(out$mergedclusters)
      cat('\n')
    }
    else if (!is.na(out$mergedclusters[1])) {
      output_list$merged_clusters <- out$mergedclusters
      cat('\nThe following clusters were merged:\n')
      print(out$mergedclusters)
      cat('\n')
    }
    output_list$ccf_table_pyclone_clean <- test_pyclone
  }


  # remove the genomically clustered clones which may be driven by undetected subclonal copy number
  if(correct_cpn_clusters)
  {
    cat('\nChecking for chromosome clustered clusters')
    # determine clonal/trunk cluster using same method as for tree building
    directedGraph_input_full <- matrix(0, 0, 2)
    colsums <- colSums(nested_pyclone$nestedclust)
    rowsums <- rowSums(nested_pyclone$nestedclust)
    trunk_cluster <- names(colsums[which(colsums == max(colsums))])
    if (length(trunk_cluster) > 1) {
      trunk_cluster <- names(sort(rowMeans(nested_pyclone$ccf_cluster_table[trunk_cluster,, drop = F]), decreasing = T))[1]
    }
    #remove the genomically clustered clones which may be driven by undetected subclonal copy number
    new_test_pyclone <- remove_clustered_clones(test_pyclone,
                                                clonal_cluster = trunk_cluster,
                                                p_value_cut = 0.01,
                                                clustering_estimate_cut = 2 )

    if(identical(new_test_pyclone,test_pyclone))
    {
      cat('\nNo clusters removed\n')
    }
    if(!identical(new_test_pyclone,test_pyclone))
    {
      #
      cat('\nThe following clusters removed due to genomic clustering:\n')
      # sort out if copy number cluster removed.
      cpn_removed_clusters <- names(nested_pyclone$cluster_qc[,'ClusterName'])[!nested_pyclone$cluster_qc[,'ClusterName']%in%unique(new_test_pyclone[,'PycloneCluster'])]
      nested_pyclone$cluster_qc[nested_pyclone$cluster_qc[,'ClusterName']%in%cpn_removed_clusters,'CopyNumRemove'] <- 1
      output_list$cpn_removed_clusters <- cpn_removed_clusters



      cat(cpn_removed_clusters)
      cat('\n')
    }

    test_pyclone    <- new_test_pyclone
    output_list$ccf_table_pyclone_clean <- test_pyclone
  }

  cat('\nBuilding trees...')
  # check whether this means all the clusters are removed.
  graph_pyclone  <- grow.trees( nestedlist = nested_pyclone
                                , pyclone = test_pyclone
                                , min_cluster_size = min_cluster_size
                                , force_trunk = TRUE
                                , ccf_buffer = ccf_buffer
  )

  output_list$tree_removed_clusters <-  graph_pyclone$Clusters_with_issues

  cat('\n------------------\n')
  cat('\nTree identified\n')

  clonality_table <- clonality.function(pyclone = test_pyclone
                                        ,trunk =graph_pyclone$trunk
                                        ,prefix = prefix
                                        , min_cluster_size = max(c(2, min_cluster_size))
                                        ,pval_cutoff = pval_cutoff
                                        ,use_boot =use_boot )
  clonality_out   <- correct.clonality.table(clonality_table = clonality_table,graph_pyclone = graph_pyclone,trunk_cluster = graph_pyclone$trunk) #TODO still may need correcting for one region cases #EC 20210509

  ### AH edit set CCF in ccf cluster table as well as upper and lower CIs to 0 if cluster is defined as absent
  for (region in colnames(clonality_table)) {
    tmp.absentClust <- rownames(clonality_table)[clonality_table[,region] == "absent"]
    if (any(nested_pyclone$ccf_cluster_table[tmp.absentClust, region] != 0)) {
      print("Absent clusters with meanCCF > 0. Resetting to 0 in ccf cluster table")
      nested_pyclone$ccf_cluster_table[tmp.absentClust, region] <- 0
      nested_pyclone$ccf_ci_lower[tmp.absentClust, region] <- 0
      nested_pyclone$ccf_ci_upper[tmp.absentClust, region] <- 0
    }
  }
  ### AH edit done

  # EC 20210509 check why printing NULL here
  # clonality_list[[i]] <- clonality_out #removed for one region cases

  if (run.multi.trees) {
    cat('\nExploring presence of multiple alternate trees')
    multi.trees     <- grow.multi.trees(nestedlist = nested_pyclone
                                        ,pyclone = test_pyclone
                                        ,graph_pyclone = graph_pyclone
                                        ,ccf_buffer = ccf_buffer
    )

  } else {
    multi.trees <- NULL
  }

  # pyclone_tree   <- prep_pyclone_tree.fn(graph_pyclone,nested_pyclone)
  graph_pyclone$alt_trees          <- multi.trees$good.trees


  if(length(multi.trees)==0)
  {
    graph_pyclone$consensus_branches <- paste(graph_pyclone$Corrected_tree[,1],graph_pyclone$Corrected_tree[,2],sep=":")
    graph_pyclone$nested_clust       <- nested_pyclone[[1]]

    # list all clone - clone relationships which are common to all alternative trees
    # This captures some tree info for clones where the exact tree position is uncertain
    graph_pyclone$consensus_relationships <- extract_consensus_relationships( list(graph_pyclone$Corrected_tree ) )
    graph_pyclone$alt_trees <- list(graph_pyclone$Corrected_tree)
  }

  if(length(multi.trees)!=0)
  {
    # check whether any repeats in alt_trees [this can happen due to level issue]
    tree_vector <- c()
    for (i in 1:length(graph_pyclone$alt_trees))
    {
      tree_vector <- c(tree_vector,PasteVector(sort(paste(graph_pyclone$alt_trees[[i]][,1],graph_pyclone$alt_trees[[i]][,2],sep=":")),sep=","))

    }

    alt_trees    <- list()
    trees_to_use <-  c(1:length(graph_pyclone$alt_trees))[!duplicated(tree_vector)]
    for (i in 1:length(trees_to_use))
    {
      alt_trees[[i]] <- graph_pyclone$alt_trees[[trees_to_use[i]]]
    }

    graph_pyclone$alt_trees <- alt_trees

    graph_pyclone$consensus_branches <- multi.trees$consensus.branches
    graph_pyclone$nested_clust       <- multi.trees$consensus.nestedclust

    # list all clone - clone relationships which are common to all alternative trees
    # This captures some tree info for clones where the exact tree position is uncertain
    graph_pyclone$consensus_relationships <- extract_consensus_relationships( alt_trees )

  }


  ### Compute alternative tree metrics:
  # 1) Compute sum condition error SCE for each alternative tree + find trees with lowest SCE:

  cat('\n\nComputing sum condition error for each alternative tree')
  graph_pyclone$alt_trees_sum_condition_error <- compute_sum_condition_error(tree_list = graph_pyclone$alt_trees, ccf_cluster_table = nested_pyclone$ccf_cluster_table, trunk = trunk_cluster)
  graph_pyclone$min_sce_trees <- names(which(graph_pyclone$alt_trees_sum_condition_error == min(graph_pyclone$alt_trees_sum_condition_error)))
  cat('\nTrees with minimum sum condition error: ', graph_pyclone$min_sce_trees, '\n')

  # 2) Compute edge probability for each alternative tree + find trees with highest edge probability:

  cat('\n\nComputing edge probability score for each alternative tree\n')
  graph_pyclone$alt_trees_edge_probability <- compute_tree_edge_probability(tree_list = graph_pyclone$alt_trees, edgelength = graph_pyclone$edgelength, trunk = trunk_cluster)
  graph_pyclone$max_edge_probability_trees <- names(which(graph_pyclone$alt_trees_edge_probability == max(graph_pyclone$alt_trees_edge_probability)))
  cat('\nTrees with maximum edge probability: ', graph_pyclone$max_edge_probability_trees, '\n')


  ### Compute clone proportions output:
  # 1) Compute subclone proportions from default tree:
  cat('\n\nComputing clone proportions from default tree\n')
  clone_proportion_table <- compute_subclone_proportions(tree_list = graph_pyclone$alt_trees,
                                                         ccf_cluster_table = nested_pyclone$ccf_cluster_table,
                                                         clonality_table = clonality_out$clonality_table_corrected,
                                                         trunk = trunk_cluster,
                                                         force_clonal_100 = TRUE,
                                                         tree_id = 1)

  # 2) Compute subclone proportions from lowest error tree:
  cat('\n\nComputing clone proportions from tree with lowest sum condition error\n')
  clone_proportions_min_sce_trees <- lapply(graph_pyclone$min_sce_trees, function(i){
    compute_subclone_proportions(tree_list = graph_pyclone$alt_trees,
                                 ccf_cluster_table = nested_pyclone$ccf_cluster_table,
                                 clonality_table = clonality_out$clonality_table_corrected,
                                 trunk = trunk_cluster,
                                 force_clonal_100 = TRUE,
                                 tree_id = as.numeric(i))
  })
  names(clone_proportions_min_sce_trees) <- graph_pyclone$min_sce_trees
  clone_proportion_out <- list(clone_proportion_table = clone_proportion_table, clone_proportions_min_sce_trees = clone_proportions_min_sce_trees)


  ### Finally, save all tree output
  # Save sample ID
  graph_pyclone$sampleID  <- sampleID
  graph_pyclone$long_sampleID  <- trx_rename.fn(sampleID, trialID = prefix)

  # Saving all output to list
  output_list$graph_pyclone                          <- graph_pyclone
  output_list$parameters                             <- input_parameter_list
  output_list$nested_pyclone                         <- nested_pyclone
  output_list$clonality_table                        <- clonality_table
  output_list$clonality_out                          <- clonality_out
  output_list$clone_proportion_out                   <- clone_proportion_out

  #let's save the output_list
  saveRDS(output_list,file=paste(generalSave,sampleID,".tree.RDS",sep = ""))

  return(output_list)
}


#' TRACERx tree plotting function
#'
#' This function is the CONIPHER function to plot the inferred phylogenetic tree.
#' @param sample_pyclone_tree A list containing all information about the
#' tree inferred using function tracerx.tree.building()
#' @importFrom grDevices "colorRampPalette" "dev.off" "pdf"
#' @importFrom graphics "abline" "axis" "barplot" "layout" "legend"
#' "par" "plot.new" "segments" "text" "title"
#' @importFrom igraph "get.edgelist"
#' @export treebuilding_plot

treebuilding_plot <- function(sample_pyclone_tree)
{
  require(mapplots)
  sampleID  <-  sample_pyclone_tree$parameters$sampleID
  dirName <- sample_pyclone_tree$parameters$dirName
  pycloneType <- sample_pyclone_tree$parameters$pycloneType
  prefix <- sample_pyclone_tree$parameters$prefix
  releaseDir <- sample_pyclone_tree$parameters$releaseDir
  saveDir <- sample_pyclone_tree$parameters$saveDir
  generalSave <- sample_pyclone_tree$parameters$generalSave
  ccf_buffer <- sample_pyclone_tree$parameters$ccf_buffer
  pval_cutoff <- sample_pyclone_tree$parameters$pval_cutoff
  use_boot <- sample_pyclone_tree$parameters$use_boot
  merge_clusters <- sample_pyclone_tree$parameters$merge_clusters
  correct_cpn_clusters <- sample_pyclone_tree$parameters$correct_cpn_clusters
  adjust_noisy_clusters <- sample_pyclone_tree$parameters$adjust_noisy_clusters
  adjust_noisy_clusters_prop <- sample_pyclone_tree$parameters$adjust_noisy_clusters_prop
  min_ccf <- sample_pyclone_tree$parameters$min_ccf
  min_cluster_size           <- sample_pyclone_tree$parameters$min_cluster_size

  nested_pyclone <- sample_pyclone_tree$nested_pyclone
  pyclone_tree   <- sample_pyclone_tree$graph_pyclone
  clonality_table <- sample_pyclone_tree$clonality_out$clonality_table_corrected
  clonality_out <- sample_pyclone_tree$clonality_out
  test_pyclone  <- sample_pyclone_tree$ccf_table_pyclone_clean
  cpn_removed_clusters <- sample_pyclone_tree$cpn_removed_clusters
  if(length(cpn_removed_clusters)==0)
  {
    cpn_removed_clusters <- NA
  }

  merged_clusters <- sample_pyclone_tree$merged_clusters

  setwd(generalSave)

  ### Plot trees -- AUTOMATIC
  date <- gsub('-', '', substr(Sys.time(), 1, 10))

  pdfname <- 'pytree_and_bar.pdf'

  height.mult.factor <- ceiling(nrow(nested_pyclone$ccf_cluster_table)/25)
  width.mult.factor  <- ceiling(nrow(nested_pyclone$ccf_cluster_table)/25)


  pdf(pdfname, width=22*width.mult.factor, height=12*height.mult.factor)
  {

    #par(mar=c(0.1,5,0.1,2),lend=1)
    par(mar=c(0,0,0,0))
    layout(cbind(1:(nrow(nested_pyclone$ccf_cluster_table)+2),rep(nrow(nested_pyclone$ccf_cluster_table)+3,nrow(nested_pyclone$ccf_cluster_table)+2),rep(nrow(nested_pyclone$ccf_cluster_table)+3,nrow(nested_pyclone$ccf_cluster_table)+2)))
    require(beeswarm)

    tmp <- nested_pyclone$ccf_cluster_table
    main <- paste(substr(colnames(tmp)[1], 1, 8), '\ Phylo CCF values', sep = '')
    colnames(tmp) <- gsub(paste0(substr(colnames(tmp)[1], 1, 8), "_"), "", colnames(tmp))
    suppressPackageStartupMessages(require(gplots))
    #textplot(tmp)
    plot.new()
    par(mar=c(2,2,2,2))
    title(main, cex = 2)


    #lend(1)
    colours.to.use <- color.tree(1:nrow(nested_pyclone$ccf_cluster_table))

    par(mar=c(0.1,5,0.1,2),lend=1)

    for (j in 1:nrow(nested_pyclone$ccf_cluster_table))
    {

      if(j==1)
      {
        border.col <- ifelse(clonality_table[j,]=='clonal','black','grey')
        bp <- barplot(nested_pyclone$ccf_cluster_table[j,],las=1,col=colours.to.use[j],border=border.col,names="",ylab=paste("Cl",rownames(nested_pyclone$ccf_cluster_table)[j],sep=" "),ylim=c(0,115),yaxt='n',cex.axis=1.25)

      }
      if(j!=1)
      {
        border.col <- ifelse(clonality_table[j,]=='clonal','black','grey')
        #border.col <- ifelse(nested_pyclone$ccf_ci_upper[j,]>=nested_pyclone$ccf_ci_lower[trunk_cluster,]-ccf_buffer,'black','grey')
        bp <- barplot(nested_pyclone$ccf_cluster_table[j,],las=1,col=colours.to.use[j],border=border.col,names="",ylab=paste("Cl",rownames(nested_pyclone$ccf_cluster_table)[j],sep=" "),ylim=c(0,115),yaxt='n',cex.axis=1.25)

      }
      axis(side = 2,at = c(0,50,100),labels=c(c(0,50,100)),las=1)
      if(j ==nrow(nested_pyclone$ccf_cluster_table))
      {
        axis(side=1,at=bp,labels=gsub(paste0(substr(colnames(nested_pyclone$ccf_cluster_table)[1], 1, 8), "_"), "",colnames(nested_pyclone$ccf_cluster_table))
             ,tick=FALSE
             #,size=2
             ,cex.axis=1.25)



      }
      abline(h=0)
      abline(h=100,lty='dashed')
      abline(h=50,lty='dashed')
      for (bar in 1:length(bp))
      {

        beeswarm(test_pyclone[test_pyclone[,'PycloneCluster']%in%rownames(nested_pyclone$ccf_cluster_table)[j],grep('PhyloCCF',colnames(test_pyclone))[bar]]*100
                 ,at=bp[bar]
                 ,add=TRUE
                 ,corralWidth = 0.5
                 ,method='swarm'
                 ,corral='wrap'
                 ,pch=21
                 ,col=colours.to.use[j]
                 ,bg='grey')
        segments(x0 = bp[bar],x1 = bp[bar],y0 = nested_pyclone$ccf_ci_lower[j,bar],y1 = nested_pyclone$ccf_ci_upper[j,bar],lwd=5)
        text(x=bp[bar],y=25,labels=nested_pyclone$ccf_cluster_table[j,bar],cex =1.5)



      }




    }

    plot.new()
    par(mar=c(2.1, 2.1, 4.1, 38), xpd=TRUE)

    g <- graph.data.frame(pyclone_tree$Corrected_tree,directed = FALSE)
    indx <- V(g)$name
    vcol <- setNames(color.tree(pyclone_tree$edgelength), names(pyclone_tree$edgelength))[indx]

    l <- layout_as_tree(g, root = pyclone_tree$trunk)

    pie.size <- ncol(sample_pyclone_tree$nested_pyclone$ccf_cluster_table)
    node.shape <- setNames(rep('pie', length(vcol)), names(vcol))
    pie.slices <- lapply(1:length(vcol), function(x) rep(1, pie.size))
    empty.col = '#bdbdbd'#'white'

    node_size_factor <- log2(max(pyclone_tree$edgelength)) / 30
    node.size <- log2(pyclone_tree$edgelength) / node_size_factor
    node.size <- node.size[names(node.shape)]

    pie.colors <- sample_pyclone_tree$nested_pyclone$ccf_cluster_table[match(names(vcol), rownames(sample_pyclone_tree$nested_pyclone$ccf_cluster_table)),, drop = F]
    pie.colors <- ifelse(pie.colors>=90,99,pie.colors)
    pie.colors <- ifelse(pie.colors<10&pie.colors>=1,10,pie.colors)
    pie.colors <- lapply(1:nrow(pie.colors), function(x) {
      if(!all(is.na(pie.colors[x,]))){
        tmp     <- pie.colors[x,]
        tmp2    <- tmp
        colfunc <- colorRampPalette(c("white", vcol[rownames(pie.colors)[x]]))
        speccolours <- colfunc(100)
        tmp[tmp>0]  <- speccolours[tmp]
        tmp[tmp2 == 0] <- empty.col
        tmp
      }
    })


    g_dir <- graph.data.frame(pyclone_tree$Corrected_tree,directed = TRUE)
    edges <- get.edgelist(g_dir)
    ecol <- setNames(rep('#bdbdbd', nrow(edges)),edges[,2])# baseline, set edge color to black
    ewidth <- rep(1,length(ecol))

    #label consensus edges in other colour
    ecol[paste(edges[,1],edges[,2],sep=":")%in%pyclone_tree$consensus_relationships] <- '#000000'
    ewidth[paste(edges[,1],edges[,2],sep=":")%in%pyclone_tree$consensus_relationships] <-150

    plot(g
         , layout=l
         , main = sampleID
         , vertex.color = vcol[indx]
         , vertex.frame.color=vcol[indx]
         , vertex.shape = node.shape
         , vertex.lwd=5
         , vertex.pie.lwd=3
         , vertex.pie = pie.slices
         , vertex.pie.color = lapply(pie.colors,rev)
         , vertex.size = node.size
         , edge.color=ecol
         , edge.size=ewidth
         , vertex.label.cex=2
         , vertex.label.pos=2
         , vertex.label.dist=0
         , vertex.label.family='Helvetica'
         , vertex.label.font=2
         , vertex.label.color = 'black')#vcol[indx])


    legend.pie(1,1,labels=gsub(paste0(substr(colnames(tmp)[1], 1, 8), "_"), "", colnames(tmp)), radius=0.2, bty="n", col='#bdbdbd',
               cex=1.25, label.dist=0.8
               ,border='white')

    snv_clusters <- sort(pyclone_tree$edgelength[indx], decreasing = T)

    snv_clusters_removed <- pyclone_tree$edgelength
    snv_clusters_removed <- sort(snv_clusters_removed[!names(snv_clusters_removed) %in% indx], decreasing = T)
    # snv_clusters_removed_cpn <- cpn_removed_clusters
    if(!is.na(cpn_removed_clusters[1]))
    {
      snv_clusters_removed <- c(snv_clusters_removed,table(sample_pyclone_tree$ccf_table_pyclone[,'PycloneCluster'])[cpn_removed_clusters])
    }

    tmp <- legend('topright', inset = c(-0.3, 0), legend = paste(names(snv_clusters), ' (', snv_clusters,' SNVs)', sep = ''), col = vcol[names(snv_clusters)], pch = 19, title = 'Clusters included:', bty = 'n')  ## inset option controls how far from x and y margins
    if (length(snv_clusters_removed) > 0) {
      if(!is.na(cpn_removed_clusters[1]))
      {
        to_plot <- table(sample_pyclone_tree$ccf_table_pyclone[,'PycloneCluster'])[cpn_removed_clusters]
        legend(x=tmp$rect$left,y = 0,inset = c(-0.3, 0),legend = paste(names(to_plot), ' (', to_plot, ' SNVs)', sep = ''), col = vcol[names(to_plot)], pch = 19, title = 'Copy# clusters removed:', bty = 'n')

      }
      legend('bottomright', inset = c(-0.3, 0), legend = paste(names(snv_clusters_removed), ' (', snv_clusters_removed, ' SNVs)', sep = ''), col = vcol[names(snv_clusters_removed)], pch = 19, title = 'Clusters removed:', bty = 'n')
    }




  }
  dev.off()

  #next, plot all the possible trees
  trees.to.plot <- pyclone_tree$alt_trees
  if(length(trees.to.plot)==0)
  {
    #nothing to plot here.
  }
  if(length(trees.to.plot)!=0)
  {

    date <- gsub('-', '', substr(Sys.time(), 1, 10))

    pdfname <- 'pytree_multipletrees.pdf'

    mult.factor <- ceiling(length(trees.to.plot)/50)

    pdf(pdfname, width=12*mult.factor, height=12*mult.factor)
    {
      nr.trees              <- length(trees.to.plot)
      columnnum             <- 1
      rownum                <- nr.trees/columnnum

      if(nr.trees<=50)
      {
        nr.to.use <- nr.trees
      }

      if(nr.trees>50)
      {
        nr.to.use <- signif(nr.trees+5,2)
      }

      for(i in 1: nr.to.use) {
        if((nr.to.use %% i) == 0) {
          if((i+(nr.to.use/i))<(columnnum+rownum))
          {
            columnnum <- i
            rownum    <- nr.to.use/columnnum
          }
        }
      }

      if(columnnum==1)
      {
        columnnum <- ceiling(columnnum*2)
        rownum    <- ceiling(rownum/2)
      }

      par(mfrow=c(rownum,columnnum),xpd=TRUE,mar=c(1, 1,1, 1))
      for (i in 1:nr.trees)
      {
        # pyclone_tree <- all_pyclone_trees[[1]]
        auto_tree    <- trees.to.plot[[i]]
        g <- graph.data.frame(auto_tree,directed = FALSE)
        indx <- V(g)$name
        vcol <- setNames(color.tree(pyclone_tree$edgelength), names(pyclone_tree$edgelength))[indx]

        l <- layout_as_tree(g, root = pyclone_tree$trunk)
        #get.edgelist(g)


        pie.size <- ncol(nested_pyclone$ccf_cluster_table)
        node.shape <- setNames(rep('pie', length(vcol)), names(vcol))
        pie.slices <- lapply(1:length(vcol), function(x) rep(1, pie.size))
        empty.col = 'gray85'#'white'

        node_size_factor <- log2(max(pyclone_tree$edgelength)) / 30
        node.size <- log2(pyclone_tree$edgelength) / node_size_factor
        node.size <- node.size[names(node.shape)]

        pie.colors <- nested_pyclone$ccf_cluster_table[match(names(vcol), rownames(nested_pyclone$ccf_cluster_table)),, drop = F]
        pie.colors <- lapply(1:nrow(pie.colors), function(x) {
          if(!all(is.na(pie.colors[x,]))){
            tmp <- pie.colors[x,]
            tmp[tmp > 0] <- vcol[rownames(pie.colors)[x]]
            tmp[tmp == 0] <- empty.col
            tmp
          }
        })

        g_dir <- graph.data.frame(auto_tree,directed = TRUE)
        edges <- get.edgelist(g_dir)
        ecol <- setNames(rep('#bdbdbd', nrow(edges)),edges[,2])# baseline, set edge color to black
        ewidth <- rep(1,length(ecol))

        #label consensus edges in other colour
        ecol[paste(edges[,1],edges[,2],sep=":")%in%pyclone_tree$consensus_relationships] <- '#000000'
        ewidth[paste(edges[,1],edges[,2],sep=":")%in%pyclone_tree$consensus_relationships] <- 2


        plot(g, main = sampleID
             , layout = l
             , vertex.color = vcol[indx]
             , vertex.shape = node.shape
             , vertex.pie = pie.slices
             , vertex.pie.color = pie.colors
             , vertex.pie.lty = 0
             , vertex.size = node.size
             , edge.width =ewidth
             , edge.color=ecol
             , arrow.size =0
             ,arrow.width=0
             ,arrow.mode=0
        )

      }




    }
    dev.off()
  }
}




