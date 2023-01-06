################################################################################
######################## Tree building helper functions ########################
################################################################################

#' Function to determine cluster nesting structure.
#'
#' This function takes as input an R list containing information about PhyloCCF
#' of each mutation (pyclone) and computes confidence intervals for each mutation
#' cluster, of the PhyloCCF distributions of the mutations in that cluster. If
#' use_boot==TRUE, then confidence intervals are computed using bootstrapping.
#' The function then performs a statistical test (Wilcoxon) for every pair of
#' clusters to determine whether one cluster can be nested within another.
#' @param pyclone An R list object containing information about the PhyloCCF
#' of each mutation in each tumour region.
#' @param prefix A character string indicating the sample and tumour case prefix
#' @param min_cluster_size Threshold for minimum number of mutations required in
#' a mutation cluster
#' @param pval_cutoff A p-value significance threshold for testing whether
#' clusters can be nested. (i.e. a p-value < pval_cutoff is significant)
#' @param min_ccf Minimum threshold for cluster PhyloCCF allowed to be classified
#' as present
#' @param use_boot Whether to use bootstrapping to determine confidence intervals
#' for each mutation cluster
#' @returns 'nestedlist', an R list containing information about the nesting
#' structure of mutation clusters in each region. Elements of the list include:
#' 'nestedclust',  'ccf_ci_lower', 'ccf_ci_upper', 'ccf_cluster_table', 'cluster_qc'.
determine.cluster.nesting <- function( pyclone
                                       , prefix = "LTX"
                                       , min_cluster_size = 5
                                       , pval_cutoff = 0.01
                                       , min_ccf= 0.01
                                       , use_boot = FALSE

) {


  suppressWarnings(require(boot))
  # To implement: Check for CCF per row does not go beyond the number of cells

  # cluster names
  pyclust    <- which(colnames(pyclone) == "PycloneCluster")

  #calc number of clusters
  nclusters  <- table(pyclone[, pyclust])

  # set minimum number of mutations of cluster to be included
  pyclone    <- pyclone[pyclone[, pyclust] %in% names(nclusters[nclusters >= min_cluster_size]), ]
  nclusters  <- nclusters[nclusters >= min_cluster_size]

  # calculate the confidence interval for each cluster
  print('Calculating cluster confidence intervals using bootstrapping')
  pyclone_ci <- calc.pyclone.ci(pyclone = pyclone, pyclust = pyclust, nclusters = nclusters, prefix = prefix)

  # if boot-strapping used, set the confidence intervals using this approach.
  if(use_boot)
  {
    pyclone_ci$ccf_ci_upper   <- pyclone_ci$ccf_ci_boot_upper
    pyclone_ci$ccf_ci_lower   <- pyclone_ci$ccf_ci_boot_lower
    pyclone_ci$mean_phylo_ccf <- pyclone_ci$mean_phylo_ccf
  }


  for (cl in grep("PhyloCCF",colnames(pyclone)))
  {
    pyclone[,cl] <- ifelse(pyclone[,cl]<min_ccf,0,pyclone[,cl])
  }

  nested_dist          <- test.distributions(pyclone = pyclone, nclusters = nclusters, pval_cutoff = pval_cutoff)
  # this simply looks at whether each cluster is significantly smaller than other clusters

  #Correct potentially wrong clusters, where very low frequency will throw off the significance testing
  clusters_with_issues <- identify.issue.clusters(nested_dist, pyclone_ci, pval_cutoff = pval_cutoff)

  #make qc matrix
  cluster_qc           <- matrix(0, nrow = nrow(nested_dist$nestedclust), ncol = 9)
  rownames(cluster_qc) <- rownames(nested_dist$nestedclust)
  colnames(cluster_qc) <- c('ClusterName', 'numSNVs',  'NotDescendFromTrunk', 'IncomingEdges', 'OutgoingEdges', 'CauseCycle', 'CopyNumRemove','SizeRemove','IssueRemove')
  cluster_qc[names(nclusters), 'ClusterName'] <- names(nclusters)
  cluster_qc[names(nclusters), 'numSNVs'] <- nclusters

  nestedlist <- list(nested_dist$nestedclust, pyclone_ci$ccf_ci_lower, pyclone_ci$ccf_ci_upper, pyclone_ci$mean_phylo_ccf, cluster_qc)
  names(nestedlist) <- c('nestedclust', 'ccf_ci_lower', 'ccf_ci_upper', 'ccf_cluster_table', 'cluster_qc')
  return(nestedlist)
}

#' Function to correct nesting based on cluster clonality
correct.clonality.nesting <- function(nestedlist
                                      ,pyclone
                                      ,clonality_table
                                      ,pval_cutoff =0.01)
{
  suppressWarnings(require(boot))
  # To implement: Check for CCF per row does not go beyond the number of cells

  # cluster names
  pyclust    <- which(colnames(pyclone) == "PycloneCluster")

  #calc number of clusters
  nclusters  <- table(pyclone[, pyclust])

  # set minimum number of mutations of cluster to be included
  pyclone    <- pyclone[pyclone[, pyclust] %in% names(nclusters[nclusters >= min_cluster_size]), ]
  nclusters  <- nclusters[nclusters >= min_cluster_size]

  # calculate the confidence interval for each
  pyclone_ci <- calc.pyclone.ci(pyclone = pyclone, pyclust = pyclust, nclusters = nclusters, prefix = prefix)

  # if boot-strapping used, set the confidence intervals using this approach.
  if(use_boot)
  {
    pyclone_ci$ccf_ci_upper   <- pyclone_ci$ccf_ci_boot_upper
    pyclone_ci$ccf_ci_lower   <- pyclone_ci$ccf_ci_boot_lower
    pyclone_ci$mean_phylo_ccf <- pyclone_ci$mean_phylo_ccf
  }


  for (cl in grep("PhyloCCF",colnames(pyclone)))
  {
    pyclone[,cl] <- ifelse(pyclone[,cl]<min_ccf,0,pyclone[,cl])
  }

  nested_dist          <- test.distributions.with.clonality.table(pyclone = pyclone, nclusters = nclusters, clonality_table = clonality_table,pval_cutoff = pval_cutoff)

  clusters_with_issues <- identify.issue.clusters(nested_dist, pyclone_ci, pval_cutoff = pval_cutoff)

  #make qc matrix
  cluster_qc           <- matrix(0, nrow = nrow(nested_dist$nestedclust), ncol = 9)
  rownames(cluster_qc) <- rownames(nested_dist$nestedclust)
  colnames(cluster_qc) <- c('ClusterName', 'numSNVs',  'NotDescendFromTrunk', 'IncomingEdges', 'OutgoingEdges', 'CauseCycle', 'CopyNumRemove','SizeRemove','IssueRemove')
  cluster_qc[names(nclusters), 'ClusterName'] <- names(nclusters)
  cluster_qc[names(nclusters), 'numSNVs'] <- nclusters

  nestedlist <- list(nested_dist$nestedclust, pyclone_ci$ccf_ci_lower, pyclone_ci$ccf_ci_upper, pyclone_ci$mean_phylo_ccf, cluster_qc)
  names(nestedlist) <- c('nestedclust', 'ccf_ci_lower', 'ccf_ci_upper', 'ccf_cluster_table', 'cluster_qc')
  return(nestedlist)

}


#' Function to determine default tree structure
#'
#' This function takes as input an R list containing information about
#' PhyloCCF of each mutation (pyclone) and an R list containing a cluster ccf table
#' and nesting matrix (nestedlist), and returns an R list containing the default tree
#' structure.
#' @param pyclone An R list object containing information about the PhyloCCF
#' of each mutation in each tumour region.
#' @param prefix A character string indicating the sample and tumour case prefix
#' @param min_cluster_size Threshold for minimum number of mutations required in
#' a mutation cluster
#' @param ccf_buffer PhyloCCF buffer permitted when checking tree level issue
grow.trees <- function(nestedlist
                       , pyclone
                       , min_cluster_size = 5
                       , ccf_buffer =10
                       , force_trunk = TRUE
                       , skip_size  =20
)
{
  # nestedlist is the output from nested.clusters.fn
  # normally expressionptions is set to 5000
  expressionsoptions <- options()[["expressions"]]
  options(expressions=1e5)
  suppressPackageStartupMessages(require(gtools))
  suppressPackageStartupMessages(require(igraph))

  nestedclust       <- nestedlist[['nestedclust']]
  ccf_ci_lower      <- nestedlist[['ccf_ci_lower']]
  ccf_ci_upper      <- nestedlist[['ccf_ci_upper']]
  ccf_cluster_table <- nestedlist[['ccf_cluster_table']]
  cluster_qc        <- nestedlist[['cluster_qc']]

  # we can only build trees with clusters that are in the ccf_cluster_table, so make sure they're all there
  clusters_to_explore <- unique(pyclone[,'PycloneCluster'])
  nestedclust        <- nestedclust[rownames(nestedclust)%in%clusters_to_explore,colnames(nestedclust)%in%clusters_to_explore]
  ccf_ci_lower        <- ccf_ci_lower[rownames(ccf_ci_lower)%in%clusters_to_explore,]
  ccf_ci_upper        <- ccf_ci_upper[rownames(ccf_ci_upper)%in%clusters_to_explore,]
  ccf_cluster_table   <- ccf_cluster_table[rownames(ccf_cluster_table)%in%clusters_to_explore,]


  # you have more than 1 cluster, proceed.
  #create directed graph
  directedGraph_input_full <- matrix(0, 0, 2)
  if (class(nestedclust) != "numeric") {
    colsums <- colSums(nestedclust)
    rowsums <- rowSums(nestedclust)
    trunk_cluster <- names(colsums[which(colsums == max(colsums))])
    nestedclust[trunk_cluster, ] <- 0
  } else {trunk_cluster <- "1"
  nestedclust[trunk_cluster] <- 0
  }
  if (length(trunk_cluster) > 1) {
    trunk_cluster <- names(sort(rowMeans(ccf_cluster_table[trunk_cluster,, drop = F]), decreasing = T))[1]
  }

  if (class(nestedclust)!="numeric"){
    if (!all(nestedclust[!rownames(nestedclust) %in% trunk_cluster, trunk_cluster] == 1)) {
      tmp <- rownames(nestedclust[!rownames(nestedclust) %in% trunk_cluster, ])[!nestedclust[!rownames(nestedclust) %in% trunk_cluster, trunk_cluster] == 1]
      #If any here, check if they are included in the trunk
      if (length(tmp) > 0) {
        for (i in 1:length(tmp)) {
          if (all(ccf_ci_upper[trunk_cluster, ] >= ccf_ci_lower[tmp[i], ])) {
            nestedclust[tmp[i], trunk_cluster] <- 1
          }
        }
      }
      tmp <- rownames(nestedclust[!rownames(nestedclust) %in% trunk_cluster, ])[!nestedclust[!rownames(nestedclust) %in% trunk_cluster, trunk_cluster] == 1]
      #If any remains
      if (length(tmp) > 0) {
        warning(paste('Cluster(s) ', tmp, ' is not descended from the trunk!!'))
        cluster_qc[tmp, 'NotDescendFromTrunk'] <- 1
        if(force_trunk) {
          nestedclust[!rownames(nestedclust) %in% trunk_cluster, trunk_cluster] <- 1
        }
      }
    }
  }

  #create new nesting based on trunk.
  pyclust       <- ncol(pyclone)
  nclusters     <- table(pyclone[, pyclust])


  if (class(ccf_ci_upper)!="numeric"){
    max_per_level <- max(max(ccf_ci_upper[trunk_cluster,]) + ccf_buffer, 100+ccf_buffer)
  } else {
    max_per_level <- max(max(ccf_ci_upper[trunk_cluster]) + ccf_buffer, 100+ccf_buffer)
  }

  #Get size of clusters
  pyclust <- ncol(pyclone)
  nclusters <- table(pyclone[, pyclust])

  if (class(nestedclust) !="numeric"){
    for (i in colnames(nestedclust)) {
      tmp <- cbind(i, rownames(nestedclust)[nestedclust[, i] == 1])
      if (ncol(tmp) == 2) {
        directedGraph_input_full <- rbind(directedGraph_input_full, tmp)
      }
    }
  } else {
    directedGraph_input_full <- rbind(directedGraph_input_full,c("1","1")) #TOFIX: EC add #need to adapt where cluster number !=1
  }
  cluster_qc[names(table(directedGraph_input_full[, 2])), 'IncomingEdges'] <- table(directedGraph_input_full[, 2])
  cluster_qc[names(table(directedGraph_input_full[, 1])), 'OutgoingEdges'] <- table(directedGraph_input_full[, 1])

  #Remove multiple incoming edges to same node
  directedGraph_input        <- directedGraph_input_full
  if(class(nestedclust)!="numeric"){
    directedGraph_input_pruned <- prune.tree_test(edgelist_full = directedGraph_input,nestedclust = nestedclust,trunk_cluster = trunk_cluster)
  } else { directedGraph_input_pruned <- directedGraph_input}
  plot(graph.data.frame(directedGraph_input_pruned),layout=layout.reingold.tilford(graph.data.frame(directedGraph_input_pruned), root=trunk_cluster) )
  directedGraph_input_pruned_unfolded <- directedGraph_input_pruned

  # identify if there are issues with the tree
  # are there circles #
  directedGraph <- graph.data.frame(directedGraph_input_pruned_unfolded)
  issue_circles <- girth(directedGraph)$girth != 0
  issue_circles <- ifelse(issue_circles==TRUE,TRUE,ifelse(max(table(directedGraph_input_pruned_unfolded[,2]))>1,TRUE,FALSE))
  cluster_qc_new <- cluster_qc
  cluster_qc_new[names(table(directedGraph_input_pruned_unfolded[, 2])), 'IncomingEdges'] <- table(directedGraph_input_pruned_unfolded[, 2])
  cluster_qc_new[names(table(directedGraph_input_pruned_unfolded[, 1])), 'OutgoingEdges'] <- table(directedGraph_input_pruned_unfolded[, 1])


  # is there a level issue (CCF>100%+buffer) in a region
  directedGraph_input_corrected <- directedGraph_input_pruned#unfold_tree(directedGraph_input_corrected, lower, trunk_cluster)[[1]]
  if (class(ccf_ci_lower) != "numeric") {
    lower                         <- ccf_ci_lower[rownames(ccf_ci_lower) %in% unique(c(directedGraph_input_corrected)),, drop = F]
  } else if (length(unique(directedGraph_input[1,]))==1) {
    ccf_ci_lower <- t(as.matrix(ccf_ci_lower))
    lower <- ccf_ci_lower
  } else {
    lower <- ccf_ci_lower[names(ccf_ci_lower) %in% unique(c(directedGraph_input_corrected)), drop = F]
  }
  tlevels                       <- check.levels.ccf(directedGraph_input_corrected, lower, trunk_cluster)
  original_tree                 <- directedGraph_input_corrected
  issue_level          <- TRUE%in%c(tlevels[[1]]>max_per_level)

  # are clusters internally consistent (this is calculated, but not used fully)
  if (length(unique(directedGraph_input[1,]))!=1) {
    internal_consistency <- check.internally.consistent(tree_small = directedGraph_input_corrected,tree_full = directedGraph_input,trunk_cluster)
    issue_consistency    <- FALSE%in%internal_consistency[,2]

    # let's try to fix any issues
    newlevels            <- tlevels
    nestedclust_to_use   <- nestedclust
    skip_next_steps      <- FALSE
    clusters_to_use      <- names(nclusters)
    clusters_to_remove   <- c()

    if(issue_level|issue_circles|issue_consistency)
    {


      #first, assess whether this can be simply explained by less strict clustering
      test_cycle  <- prune.cycle(egdgelist = directedGraph_input_corrected
                                 ,nestedclust = nestedclust
                                 ,ccf_ci_lower = ccf_ci_lower
                                 ,max_per_level = max_per_level
                                 ,trunk_cluster = trunk_cluster
      )
      #issue_consistency <- check.internally.consistent(test_cycle$directedGraph_input_pruned_unfolded,tree_full = directedGraph_input,trunk_cluster)
      if(test_cycle$levelissue %in% FALSE & test_cycle$circle %in% FALSE & test_cycle$consistencyissue %in% FALSE)
      {
        #report this now, with the latest case....
        nestedclust_to_use <- test_cycle$nestedcluster
        skip_next_steps    <- TRUE
      }

      #let's try removing clusters and see if that can solve this issue (would have been good to integrate with above, but we'll do that for another iteration)
      clusters_decreasing_order               <- names(sort(nclusters[nclusters>=min_cluster_size][names(nclusters[nclusters>=min_cluster_size])%in%clusters_to_use],decreasing=FALSE))
      clusters_decreasing_order_without_trunk <- clusters_decreasing_order[!clusters_decreasing_order%in%trunk_cluster]


      if(length(clusters_decreasing_order_without_trunk)>skip_size&skip_next_steps%in%FALSE)
      {
        clusters_with_multiple_daughters <- colnames(nestedclust_to_use)[colSums(nestedclust_to_use)>0]
        clusters_with_no_daughters       <- colnames(nestedclust_to_use)[colSums(nestedclust_to_use)==0]
        endterminals                     <- directedGraph_input_pruned_unfolded[,2][!directedGraph_input_pruned_unfolded[,2]%in%directedGraph_input_pruned_unfolded[,1]]
        clusters_with_multiple_daughters <- clusters_with_multiple_daughters[!clusters_with_multiple_daughters%in%endterminals]
        clusters_to_remove_for_now       <- unique(c(directedGraph_input_pruned_unfolded))[!unique(c(directedGraph_input_pruned_unfolded))%in%unique(c(clusters_with_multiple_daughters))]

        still_keep                       <- names(which(nclusters[!names(nclusters)%in%trunk_cluster]/sum(nclusters[!names(nclusters)%in%trunk_cluster])>1/length(nclusters)))
        dont_keep                        <- names(which(nclusters[!names(nclusters)%in%trunk_cluster]/sum(nclusters[!names(nclusters)%in%trunk_cluster])<0.5/length(nclusters)))

        clusters_to_remove_for_now       <- clusters_to_remove_for_now[!clusters_to_remove_for_now%in%still_keep]
        clusters_to_remove_for_now       <- unique(c(clusters_to_remove_for_now,dont_keep))


        #check how many clusters we're going to be assessing for permutation
        all_clusters_to_test             <- clusters_to_use[!clusters_to_use%in%clusters_to_remove_for_now]
        test_out                         <- is.there.ccf.issue(nestedclust = nestedclust_to_use,directed_input_graph = directedGraph_input_corrected,ccf_ci_lower = ccf_ci_lower,trunk_cluster = trunk_cluster,clusters_to_remove = clusters_to_remove_for_now,clusters_to_use = clusters_to_use,max_per_level = max_per_level)

        if(length(all_clusters_to_test)>25)
        {
          clusters_to_remove_for_now       <- c(all_clusters_to_test[-c(1:25)],clusters_to_remove_for_now)
          all_clusters_to_test             <- clusters_to_use[!clusters_to_use%in%clusters_to_remove_for_now]
          test_out                         <- is.there.ccf.issue(nestedclust = nestedclust_to_use,directed_input_graph = directedGraph_input_corrected,ccf_ci_lower = ccf_ci_lower,trunk_cluster = trunk_cluster,clusters_to_remove = clusters_to_remove_for_now,clusters_to_use = clusters_to_use,max_per_level = max_per_level)

        }

        if (nrow(test_out$directed_input_corrected) != 1) {
          # try removing clusters
          new_test_out                     <- permute.clusters.to.remove(test_out = test_out,nestedclust = nestedclust_to_use,max_per_level = max_per_level,tlevels = tlevels,trunk_cluster = trunk_cluster,cluster_qc = cluster_qc,ccf_ci_lower = ccf_ci_lower,nclusters=nclusters)
          new_clusters_to_test             <- unique(c(clusters_to_remove_for_now,new_test_out$clusters.to.remove))
        } else {
          new_test_out <- test_out
          new_clusters_to_test <- clusters_to_remove_for_now
        }
        new_clusters_to_test             <- as.character(sort(as.numeric(new_clusters_to_test),decreasing=FALSE))

        # try removed clusters back
        new_test_out                     <- permute.clusters.to.add.back(new_test_out = new_test_out,nestedclust = nestedclust_to_use,max_per_level = max_per_level,tlevels = tlevels,trunk_cluster = trunk_cluster,cluster_qc = cluster_qc,ccf_ci_lower = ccf_ci_lower,ccf_ci_upper = ccf_ci_upper,new_clusters_to_test = new_clusters_to_test)
        clusters_to_remove               <- clusters_to_use[!clusters_to_use%in%as.character(sort(as.numeric(unique(c(new_test_out$directed_input_corrected)))))]
        nestedclust_to_use               <- new_test_out$nestedclust_new

        # try adding rest of clusters back
        new_test_out                     <- permute.clusters.to.add.back(new_test_out = new_test_out,nestedclust = nestedclust_to_use,max_per_level = max_per_level,tlevels = tlevels,trunk_cluster = trunk_cluster,cluster_qc = cluster_qc,ccf_ci_lower = ccf_ci_lower,ccf_ci_upper = ccf_ci_upper,new_clusters_to_test = clusters_to_remove)
        nestedclust_to_use               <- new_test_out$nestedclust_new

        #final list of clusters to remove and use
        clusters_to_remove               <- clusters_to_use[!clusters_to_use%in%as.character(sort(as.numeric(unique(c(new_test_out$directed_input_corrected)))))]
        clusters_to_use                  <- as.character(sort(as.numeric(unique(c(new_test_out$directed_input_corrected)))))


      }

      if(length(clusters_decreasing_order_without_trunk)<=skip_size&skip_next_steps%in%FALSE)
      {
        test_out                         <- is.there.ccf.issue(nestedclust = nestedclust_to_use,directed_input_graph = directedGraph_input_corrected,ccf_ci_lower = ccf_ci_lower,trunk_cluster = trunk_cluster,clusters_to_remove = c(),clusters_to_use = clusters_to_use,max_per_level = max_per_level)

        # try removing clusters
        new_test_out                     <- permute.clusters.to.remove(test_out = test_out
                                                                       ,nestedclust = nestedclust_to_use
                                                                       ,max_per_level = max_per_level
                                                                       ,tlevels = tlevels
                                                                       ,trunk_cluster = trunk_cluster
                                                                       ,cluster_qc = cluster_qc
                                                                       ,ccf_ci_lower = ccf_ci_lower
                                                                       ,nclusters=nclusters
        )


        new_clusters_to_test             <- unique(c(new_test_out$clusters_to_remove))
        new_clusters_to_test             <- c(clusters_decreasing_order_without_trunk,as.character(sort(as.numeric(new_clusters_to_test),decreasing=FALSE)))
        new_clusters_to_test             <- new_clusters_to_test[!new_clusters_to_test%in%unique(c(new_test_out$directed_input_corrected))]
        clusters_to_remove               <- as.character(sort(as.numeric(unique(new_clusters_to_test))))

        # try adding them back in
        new_test_out                     <- permute.clusters.to.add.back(new_test_out = new_test_out,nestedclust = nestedclust_to_use,max_per_level = max_per_level,tlevels = tlevels,trunk_cluster = trunk_cluster,cluster_qc = cluster_qc,ccf_ci_lower = ccf_ci_lower,ccf_ci_upper = ccf_ci_upper,new_clusters_to_test = clusters_to_remove)
        nestedclust_to_use               <- new_test_out$nestedclust_new

        # final list of clusters to remove and include
        clusters_to_remove               <- clusters_to_use[!clusters_to_use%in%as.character(sort(as.numeric(unique(c(new_test_out$directed_input_corrected)))))]
        clusters_to_use                  <- as.character(sort(as.numeric(unique(c(new_test_out$directed_input_corrected)))))


      }

    }


    #now create the final tree output with these all added.
    clusters_to_use <- c(clusters_to_use)
    nestedclust_test     <- nestedclust_to_use[rownames(nestedclust_to_use)%in%clusters_to_use,colnames(nestedclust_to_use)%in%clusters_to_use]
    directedGraph_input_full_test <- matrix(0, 0, 2)
    colsums <- colSums(nestedclust_test)
    rowsums <- rowSums(nestedclust_test)

    nestedclust_test[trunk_cluster, ] <- 0
    if (!all(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1)) {
      tmp <- rownames(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, ])[!nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1]
      #If any here, check if they are included in the trunk
      if (length(tmp) > 0) {
        for (i in 1:length(tmp)) {
          if (all(ccf_ci_upper[trunk_cluster, ] >= ccf_ci_lower[tmp[i], ])) {
            nestedclust_test[tmp[i], trunk_cluster] <- 1
          }
        }
      }
      tmp <- rownames(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, ])[!nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1]
      #If any remains
      if (length(tmp) > 0) {
        warning(paste('Cluster(s) ', tmp, ' is not descended from the trunk!!'))
        cluster_qc[tmp, 'NotDescendFromTrunk'] <- 1
        if(force_trunk) {
          nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] <- 1
        }
      }
    }
    for (i in colnames(nestedclust_test)) {
      tmp <- cbind(i, rownames(nestedclust_test)[nestedclust_test[, i] == 1])
      if (ncol(tmp) == 2) {
        directedGraph_input_full_test <- rbind(directedGraph_input_full_test, tmp)
      }
    }

    directedGraph_input_pruned <- prune.tree_test(edgelist_full = directedGraph_input_full_test,nestedclust = nestedclust_test,trunk_cluster = trunk_cluster)
  } else {clusters_to_remove <- NA}

  directedGraph_input_pruned_unfolded <- directedGraph_input_pruned
  tlevels <- check.levels.ccf(directedGraph_input_pruned_unfolded, ccf_ci_lower, trunk_cluster)

  outlist <- list(directedGraph_input_pruned_unfolded, original_tree, c(clusters_to_remove), tlevels, nclusters, trunk_cluster)
  names(outlist) <- c('Corrected_tree', 'Original_tree', 'Clusters_with_issues', 'ccf_per_level', 'edgelength', 'trunk')

  if (class(ccf_cluster_table) != "numeric") {
    if (nrow(ccf_cluster_table) == 1) {
      outlist <- list(Corrected_tree = matrix(as.character(c(1, 1)), ncol = 2, byrow = T), 'Original_tree' = matrix(as.character(c(1, 1)), ncol = 2, byrow = T), 'Clusters_with_issues' = NA, 'ccf_per_level' = NA, 'edgelength' = setNames(cluster_qc[1, 1], rownames(ccf_cluster_table)[1]), 'trunk' = rownames(ccf_cluster_table)[1])
    }}

  options(expressions=expressionsoptions)

  return(outlist)
}


# this function is unfinished and doesn't work - but highlights need to check which clusters are removed and why.
check.why.cluster.excluded <- function(nestedlist,max_per_level,trunk_cluster,tree_removed_clusters)
{

  expressionsoptions <- options()[["expressions"]]
  options(expressions=1e5)
  suppressPackageStartupMessages(require(gtools))
  suppressPackageStartupMessages(require(igraph))

  nestedclust      <- nestedlist[[1]]
  ccf_ci_lower     <- nestedlist[[2]]
  ccf_ci_upper     <- nestedlist[[3]]
  ccf_cluster_table <- nestedlist[[4]]
  cluster_qc       <- nestedlist[[5]]

  # you have more than 1 cluster, proceed.
  #create directed graph
  directedGraph_input_full <- matrix(0, 0, 2)
  colsums <- colSums(nestedclust)
  rowsums <- rowSums(nestedclust)
  trunk_cluster <- names(colsums[which(colsums == max(colsums))])
  if (length(trunk_cluster) > 1) {
    trunk_cluster <- names(sort(rowMeans(ccf_cluster_table[trunk_cluster,, drop = F]), decreasing = T))[1]
  }
  #Check that all clusters are decended from the trunk, and ensure that the trunk does not descend from any
  nestedclust[trunk_cluster, ] <- 0
  if (!all(nestedclust[!rownames(nestedclust) %in% trunk_cluster, trunk_cluster] == 1)) {
    tmp <- rownames(nestedclust[!rownames(nestedclust) %in% trunk_cluster, ])[!nestedclust[!rownames(nestedclust) %in% trunk_cluster, trunk_cluster] == 1]
    #If any here, check if they are included in the trunk
    if (length(tmp) > 0) {
      for (i in 1:length(tmp)) {
        if (all(ccf_ci_upper[trunk_cluster, ] >= ccf_ci_lower[tmp[i], ])) {
          nestedclust[tmp[i], trunk_cluster] <- 1
        }
      }
    }
    tmp <- rownames(nestedclust[!rownames(nestedclust) %in% trunk_cluster, ])[!nestedclust[!rownames(nestedclust) %in% trunk_cluster, trunk_cluster] == 1]
    #If any remains
    if (length(tmp) > 0) {
      warning(paste('Cluster(s) ', tmp, ' is not descended from the trunk!!'))
      cluster_qc[tmp, 'NotDescendFromTrunk'] <- 1
      if(force_trunk) {
        nestedclust[!rownames(nestedclust) %in% trunk_cluster, trunk_cluster] <- 1
      }
    }
  }

  #create new nesting based on trunk.
  pyclust       <- ncol(pyclone)
  nclusters     <- table(pyclone[, pyclust])
  max_per_level <- max(ccf_ci_upper[trunk_cluster,])+ccf_buffer

  for (i in colnames(nestedclust)) {
    tmp <- cbind(i, rownames(nestedclust)[nestedclust[, i] == 1])
    if (ncol(tmp) == 2) {
      directedGraph_input_full <- rbind(directedGraph_input_full, tmp)
    }
  }

  #Remove multiple incoming edges to same node
  directedGraph_input        <- directedGraph_input_full
  directedGraph_input_pruned <- prune.tree_test(edgelist_full = directedGraph_input,
                                                nestedclust = nestedclust,
                                                trunk_cluster = trunk_cluster)


  ccf.test                   <- is.there.ccf.issue(nestedclust = nestedclust,
                                                   directed_input_graph = directedGraph_input_pruned,
                                                   ccf_ci_lower = ccf_ci_lower,
                                                   trunk_cluster = trunk_cluster,
                                                   clusters_to_remove = c(),
                                                   clusters_to_use = clusters.in.use.test,
                                                   max_per_level = max_per_level)

  nestedclust_to_use <- update.nested.cluster(nestedclust,ccf_ci_lower = ccf_ci_lower,max_per_level = max(ccf_ci_upper[trunk_cluster,]),trunk_cluster)

  plot(graph.data.frame(directedGraph_input_pruned),layout=layout.reingold.tilford(graph.data.frame(directedGraph_input_pruned), root=trunk_cluster) )
  directedGraph_input_pruned_unfolded <- directedGraph_input_pruned


  for(tree_removed_cluster in tree_removed_clusters)
  {
    descendents_of_cluster <- names(nestedclust_to_use[,tree_removed_cluster][nestedclust_to_use[,tree_removed_cluster]==1])
    parents_of_cluster     <- names(nestedclust_to_use[tree_removed_cluster,][nestedclust_to_use[tree_removed_cluster,]==1])

    #check which descendents have other parents
    parent_clash <- c()
    for (descendent_of_cluster in descendents_of_cluster)
    {

    }
  }
}


#' Function to determine cluster clonality in each tumour region
#'
#' This function takes as input an R list containing information about PhyloCCF
#' of each mutation (pyclone) and computes confidence intervals for each mutation
#' cluster, of the PhyloCCF distributions of the mutations in that cluster. If
#' use_boot==TRUE, then confidence intervals are computed using bootstrapping.
#' The function then performs a statistical test (Wilcoxon) for every pair of
#' clusters to determine whether one cluster can be nested within another.
#' @param pyclone An R list object containing information about the PhyloCCF
#' of each mutation in each tumour region.
#' @param trunk truncal cluster name
#' @param ccf_buffer PhyloCCF buffer permitted when checking tree level issue
#' @param min_cluster_size Threshold for minimum number of mutations required in
#' a mutation cluster
#' @param pval_cutoff A p-value significance threshold for testing whether
#' clusters can be nested. (i.e. a p-value < pval_cutoff is significant)
#' @param min_ccf Minimum threshold for cluster PhyloCCF allowed to be classified
#' as present
#' @param use_boot Whether to use bootstrapping to determine confidence intervals
#' for each mutation cluster
#' @returns a matrix of dimensions (n_clusters x n_regions) classifying each
#' cluster as 'clonal', 'subclonal', or 'absent' in each tumour region.
clonality.function <- function(pyclone
                               , trunk
                               , ccf_buffer=10
                               , prefix = "LTX"
                               , min_cluster_size = 5
                               , pval_cutoff = 0.01
                               , min_ccf= 0.05
                               , use_boot = FALSE)
{
  suppressWarnings(require(boot))
  # To implement: Check for CCF per row does not go beyond the number of cells

  pyclust    <- which(colnames(pyclone) == "PycloneCluster")
  nclusters  <- table(pyclone[, pyclust])
  pyclone    <- pyclone[pyclone[, pyclust] %in% names(nclusters[nclusters >= min_cluster_size]), ]
  nclusters  <- nclusters[nclusters >= min_cluster_size]
  pyclone_ci <- calc.pyclone.ci(pyclone = pyclone, pyclust = pyclust, nclusters = nclusters, prefix = prefix)


  if(use_boot)
  {
    pyclone_ci$ccf_ci_upper   <- pyclone_ci$ccf_ci_boot_upper
    pyclone_ci$ccf_ci_lower   <- pyclone_ci$ccf_ci_boot_lower
    pyclone_ci$mean_phylo_ccf <- pyclone_ci$mean_phylo_ccf
  }

  nested_dist          <- test.distributions(pyclone, nclusters, pval_cutoff = pval_cutoff)

  #create a table where everything is clonal
  clonal.matrix <- matrix(1,ncol=length(colnames(pyclone_ci$ccf_cluster_table))
                          ,nrow=length(rownames(pyclone_ci$ccf_cluster_table)))
  colnames(clonal.matrix) <- colnames(pyclone_ci$ccf_cluster_table)
  rownames(clonal.matrix) <- rownames(pyclone_ci$ccf_cluster_table)
  clonal.matrix[trunk,] <- 1
  lower.ccf.limit         <- ifelse(pyclone_ci$ccf_ci_lower[trunk,]-ccf_buffer<90,pyclone_ci$ccf_ci_lower[trunk,]-ccf_buffer,90)
  for (region in colnames(clonal.matrix))
  { if (nrow(clonal.matrix)!=1) {
    if (ncol(clonal.matrix)!=1) {
      clonal.matrix[,region] <- ifelse(nested_dist[[region]][trunk,]<pval_cutoff&pyclone_ci$ccf_ci_upper[,region]<=lower.ccf.limit[region],0,1) #TODO check this for one region case(?)
    } else {
      clonal.matrix[,region] <- ifelse(nested_dist[[region]][trunk,]<pval_cutoff&pyclone_ci$ccf_ci_upper[,region]<=lower.ccf.limit,0,1)
    }} else {
      clonal.matrix[,region] <- ifelse(nested_dist[[region]][trunk,]<pval_cutoff&pyclone_ci$ccf_ci_upper[,region]<=lower.ccf.limit[region],0,1)
    }}

  for (region in colnames(clonal.matrix))
  {
    clonal.matrix[,region] <- ifelse(pyclone_ci$mean_phylo_ccf[,region]<=0,-1,clonal.matrix[,region])
  }

  clonal.matrix <- data.frame(ifelse(clonal.matrix==1,'clonal',ifelse(clonal.matrix==-1,'absent','subclonal')))

  return(clonal.matrix)


}

#' Function to create a list of all tree paths
#' @param pyclone An R matrix describing the tree structure with two columns
#' specifying 'parent' (column 1) and child (column 2)
#' @returns An R list of all tree paths from trunk to leaves
createAllPathsList <- function(tree.structure, trunk)
{
  tree.structure <- as.data.frame(tree.structure, stringsAsFactors = F)
  colnames(tree.structure) <- c("V1", "V2")
  leaves <- tree.structure$V2[!tree.structure$V2%in%tree.structure$V1]
  all.paths <- NULL
  if(length(leaves)==0){
    all.paths <- paste(as.character(tree.structure[1,]), collapse = ",")
  }
  while(length(leaves>0)){
    l <- leaves[1]
    preceding <- tree.structure$V1[tree.structure$V2%in%l]
    path <- paste(c(preceding, l), collapse = ",")
    while(any(preceding!=trunk)){
      preceding <- tree.structure$V1[tree.structure$V2%in%preceding]
      path <- paste(c(preceding, path), collapse = ",")
    }
    all.paths <- c(all.paths, list(path))
    leaves <- leaves[!leaves%in%l]
  }
  return(all.paths)
}

tree.structure.adjust.clonality <- function(tree.structure, trunk, clonality.called.cluster.region.df)
{
  allpaths <- createAllPathsList(tree.structure, trunk)
  allpaths <- unlist(allpaths)
  allpaths <- strsplit(allpaths, split = ",")
  clusters.to.keep <- unique(as.character(tree.structure))
  clonality.called.cluster.region.df <- clonality.called.cluster.region.df[rownames(clonality.called.cluster.region.df)%in%clusters.to.keep,]
  updated.clonality.df <- clonality.called.cluster.region.df
  tree.structure <- as.data.frame(tree.structure, stringsAsFactors = F)
  colnames(tree.structure) <- c("V1", "V2")
  error <- NULL
  if (class(clonality.called.cluster.region.df) != "factor") {
    regions <- colnames(clonality.called.cluster.region.df)
    for(i in 1:length(regions)){
      clusters.in.region <- rownames(clonality.called.cluster.region.df)[clonality.called.cluster.region.df[,i]%in%c("clonal", "subclonal")]
      clusters.in.region <- clusters.in.region[!clusters.in.region%in%trunk]
      if(length(clusters.in.region)==0){
        next()
      }
      paths.to.keep <- NULL
      for(j in 1:length(clusters.in.region)){
        c <- clusters.in.region[j]
        for(k in 1:length(allpaths)){
          p <- allpaths[[k]]
          if(c%in%p){
            paths.to.keep <- c(paths.to.keep, list(p))
          }
        }
      }
      paths.to.keep <- unique(paths.to.keep)

      ## cannot have a region move from subclonal to clonal; or absent to subclonal; or absent to clonal; currently only flags as opposed to corrects
      for(l in 1:length(paths.to.keep)){
        p <- paths.to.keep[[l]]
        if(length(p)==1){
          next()
        }else{
          for(m in 2:length(p)){
            cluster <- p[m]
            preceding.cluster <- p[m-1]
            transition <- paste0(clonality.called.cluster.region.df[rownames(clonality.called.cluster.region.df)%in%preceding.cluster,regions[i]],
                                 ">",
                                 clonality.called.cluster.region.df[rownames(clonality.called.cluster.region.df)%in%cluster,regions[i]])
            if(transition%in%c("absent>subclonal", "absent>clonal", "subclonal>clonal")){
              error <- c(error, paste(regions[i], preceding.cluster, cluster,"badtransition", transition, sep=":"))
            }
          }
        }
        error <- unique(error)
      }

      ## now we have to check that a region that is in multiple branches of the tree cannot have multiple branches that have clonal clusters post MRCA
      if(length(paths.to.keep)==1){
        next()
      }else{
        clonal.clusters.in.region <- rownames(clonality.called.cluster.region.df)[clonality.called.cluster.region.df[ ,regions[i]]%in%"clonal"]
        clonal.clusters.in.region <- clonal.clusters.in.region[!clonal.clusters.in.region%in%trunk]
        all.clusters.in.regions <- rownames(clonality.called.cluster.region.df)[clonality.called.cluster.region.df[ ,regions[i]]%in%c("subclonal", "clonal")]
        all.consensus.paths <- NULL
        ##get the consensus paths for this region
        for(a in 1:length(paths.to.keep)){
          path <- paths.to.keep[[a]]
          tmp.consensus.path <- NULL
          for(b in 1:length(path)){
            if(path[b]%in%all.clusters.in.regions){
              tmp.consensus.path <- c(tmp.consensus.path, path[b])
            }
          }
          all.consensus.paths <- c(all.consensus.paths, list(tmp.consensus.path))
        }
        all.consensus.paths <- unique(all.consensus.paths)

        if(length(all.consensus.paths)>1){
          ##identify the "most recent common ancestor" of the consensus paths and the branching clones
          max.length.consensus.path <- max(unlist(lapply(all.consensus.paths, length)))
          tree.df <- plyr::ldply(all.consensus.paths, rbind)
          tree.df <- data.frame(tree.df, stringsAsFactors = F)
          for(c in 1:ncol(tree.df)){
            if(length(unique(tree.df[,c]))==1){
              MRCA <- as.character(unique(tree.df[,c]))
            }else{
              branching.clones <- as.character(unique(tree.df[,c]))
              break()
            }
          }
          branching.clones <- branching.clones[!is.na(branching.clones)]
          if(length(branching.clones)>1){
            if(any(branching.clones%in%clonal.clusters.in.region)){
              error <- c(error, paste(regions[i],  "multipleclonalbranches", paste(branching.clones, collapse = ","), sep=":"))
            }
          }

        }

      }

    }}
  print(error)
  return(error)
}

correct.clonality.table <- function(clonality_table,graph_pyclone,trunk_cluster)
{
  tmp <- data.frame(lapply(clonality_table, as.character), stringsAsFactors = FALSE)
  rownames(tmp) <- rownames(clonality_table)
  clonality_table <- tmp
  rm(tmp)
  original_clonality_table <- clonality_table

  if (nrow(clonality_table) != 1 & ncol(clonality_table) != 1) {

    tmp <- tree.structure.adjust.clonality(tree.structure = graph_pyclone$Corrected_tree,trunk = trunk_cluster,clonality.called.cluster.region.df = clonality_table)
    if(length(tmp)==0)
    {
      out <- list()
      out$clonality_table_corrected <- clonality_table
      out$clonality_table_original <- original_clonality_table
      return(out)
    }

    tmp <- grep("multipleclonalbranches",tmp,value=TRUE)
    corrected.clusters <- c()

    if(length(tmp)>0)
    {  for (ln in 1:length(tmp))
    {
      info <- unlist(strsplit(tmp[ln],split=":"))
      region <- info[1]
      clusters <- unlist(strsplit(info[3],split="\\,"))
      corrected.clusters <- c(corrected.clusters,clusters)
      clonality_table[clusters,region] <- 'subclonal'

    }

    }

    corrected.clusters <- unique(corrected.clusters)
    test.tree <- tree.structure.adjust.clonality(tree.structure = graph_pyclone$Corrected_tree,trunk = trunk_cluster,clonality.called.cluster.region.df = clonality_table)
    badtransition.present          <- TRUE%in%ifelse(grep('badtransition',test.tree)>0,TRUE,FALSE)
    multipleclonalbranches.present <- TRUE%in%ifelse(grep('multipleclonalbranches',test.tree)>0,TRUE,FALSE)
    k <- 1
    while((badtransition.present%in%TRUE|multipleclonalbranches.present%in%TRUE)&k<50)
    {
      tmp1 <- tree.structure.adjust.clonality(tree.structure = graph_pyclone$Corrected_tree,trunk = trunk_cluster,clonality.called.cluster.region.df = clonality_table)
      tmp2 <- grep("badtransition",tmp1,value=TRUE)
      corrected.clusters2 <- corrected.clusters

      if(length(tmp2)>0)
      {
        for (ln in 1:length(tmp2))
        {
          info       <- unlist(strsplit(tmp2[ln],split=":"))
          region     <- info[1]
          cluster1   <- info[2]
          cluster2   <- info[3]
          #info[5]
          #
          if (cluster1%in%corrected.clusters)
          {
            # assume cluster 1 is correct
            new.labelling <- unlist(strsplit(info[5],split=">"))[1]
            clonality_table[cluster2,region] <- new.labelling


          }
          if (cluster2%in%corrected.clusters)
          {
            new.labelling <- unlist(strsplit(info[5],split=">"))[2]
            clonality_table[cluster1,region] <- new.labelling


          }

          if (!TRUE%in%c(cluster2%in%corrected.clusters,cluster1%in%corrected.clusters))
          {
            #make decision based on cluster size
            if(as.numeric(cluster1)>=as.numeric(cluster2))
            {
              # assume cluster 1 correct
              new.labelling <- unlist(strsplit(info[5],split=">"))[1]
              clonality_table[cluster2,region] <- new.labelling

            }
            if(as.numeric(cluster1)<as.numeric(cluster2))
            {
              # assume cluster 2 correct
              new.labelling <- unlist(strsplit(info[5],split=">"))[2]
              clonality_table[cluster1,region] <- new.labelling

            }


          }

          corrected.clusters2 <- c(corrected.clusters2,cluster1)
          corrected.clusters2 <- c(corrected.clusters2,cluster2)




        }
        corrected.clusters2 <- unique(corrected.clusters2)
      }


      tmp3 <- grep("multipleclonalbranches",tmp1,value=TRUE)
      # corrected.clusters <- c()
      if(length(tmp3)>0)
      {
        for (ln in 1:length(tmp3))
        {
          info <- unlist(strsplit(tmp3[ln],split=":"))
          region <- info[1]
          clusters <- unlist(strsplit(info[3],split="\\,"))
          corrected.clusters2 <- c(corrected.clusters2,clusters)
          clonality_table[clusters,region] <- 'subclonal'

        }
      }


      corrected.clusters <- corrected.clusters2
      test.tree <- tree.structure.adjust.clonality(tree.structure = graph_pyclone$Corrected_tree,trunk = trunk_cluster,clonality.called.cluster.region.df = clonality_table)
      badtransition.present          <- TRUE%in%ifelse(grep('badtransition',test.tree)>0,TRUE,FALSE)
      multipleclonalbranches.present <- TRUE%in%ifelse(grep('multipleclonalbranches',test.tree)>0,TRUE,FALSE)
      k <- k+1

    }
  }

  out <- list()
  out$clonality_table_corrected <- clonality_table
  out$clonality_table_original <- original_clonality_table
  return(out)
}

#' Function to determine all possible alternative phylogenies
#'
#' @param nestedlist An R list containing information about the nesting
#' structure of mutation clusters in each region.
#' @param graph_pyclone An R list containing information about the tree structure
#' @param pyclone An R list containing information about mutation PhyloCCF
#' @param ccf_buffer PhyloCCF buffer permitted when checking tree level issue
#' @returns An R list containing all possible alternative tree structures and
#' information about which branches are consensus across multiple trees
grow.multi.trees <- function(nestedlist,graph_pyclone,pyclone,ccf_buffer=10)
{
  suppressWarnings(require(gtools))
  nestedclust       <- nestedlist[[1]]
  ccf_ci_lower      <- nestedlist[[2]]
  ccf_ci_upper      <- nestedlist[[3]]
  ccf_cluster_table <- nestedlist[[4]]
  cluster_qc        <- nestedlist[[5]]

  nestedclust <- nestedclust[,colnames(nestedclust)%in%unique(c(graph_pyclone$Corrected_tree)),drop=FALSE]
  nestedclust <- nestedclust[rownames(nestedclust)%in%unique(c(graph_pyclone$Corrected_tree)),,drop=FALSE]

  trunk_cluster       <- graph_pyclone$trunk
  max_per_level       <- max(max(ccf_ci_upper[trunk_cluster,])+ccf_buffer,100+ccf_buffer)

  non_trunk_clusters  <- c(graph_pyclone$Corrected_tree)[!c(graph_pyclone$Corrected_tree)%in%trunk_cluster]
  non_trunk_clusters  <- unique(non_trunk_clusters[order(as.numeric(non_trunk_clusters),decreasing = FALSE)])

  # for altering the tree (I don't include the CCF allowance)
  strict_nestedclust  <- update.nested.cluster(nestedclust,ccf_ci_lower = ccf_ci_lower,max_per_level = max(ccf_ci_upper[trunk_cluster,]),trunk_cluster)

  #check which clusters could be changed potentially
  clusters_with_potential_to_move <- c()
  for (cluster in non_trunk_clusters)
  {
    if(identical(nestedclust[cluster,],strict_nestedclust[cluster,]))
    {
      next;
    }
    clusters_with_potential_to_move <- c(clusters_with_potential_to_move,cluster)

  }

  if(length(clusters_with_potential_to_move)==0)
  {
    return(NULL)
    # need to return
  }

  #now, let's assume we can create new trees with all of these cluster
  # first let's determine the possible ways we could do that
  potential.trees <- list()
  k <- 1
  cat('\nCalculating combinations of clusters to remove\n')
  # let's limit to 6 possible combinations of clusters being moved (otherwise we get a segfault)
  if(length(clusters_with_potential_to_move)>6)
  {
    warning('Restricting number of clusters to simultaneously move to 6')
  }
  max_clusters_to_move <- min(length(clusters_with_potential_to_move),6)
  for (r in 1:max_clusters_to_move)
  {
    cat('\n%complete:',round(((r-1)/max_clusters_to_move)*100),'')
    comb.out <- combinations(v = clusters_with_potential_to_move,r=r,n=length(clusters_with_potential_to_move))
    for (i in 1:nrow(comb.out))
    {
      potential.trees[[k]] <- c(comb.out[i,])
      k <- k+1
    }
  }

  cat('\npotential combinations of trees to remove calculated')
  # now we know some of potential number of trees, let's see how many of these don't break the rules
  trees               <-  list()
  trees[["original"]]$tree    <- graph_pyclone$Corrected_tree
  trees[["original"]]$maxccf  <- max(graph_pyclone$ccf_per_level[[1]])
  trees[["original"]]$nesting <- nestedclust
  #next add all the double trees
  trees[["1_cluster_change"]] <- list()
  cat('\n\ntry moving all possible clusters')

  for (cluster in clusters_with_potential_to_move)
  {
    nesting.to.remove <- names(which(c(nestedclust[cluster,]-strict_nestedclust[cluster,])==1))
    expanded_tree     <- expand.tree(graph_pyclone$Corrected_tree,trunk_cluster)

    nesting.to.remove <- nesting.to.remove[paste(nesting.to.remove,cluster,sep=":")%in%paste(expanded_tree[,1],expanded_tree[,2],sep=":")]

    # make sure all is in the tree
    nestedclust_full  <- nestedclust
    nestedclust_full[cluster,] <- 0
    nestedclust_full[cluster,expanded_tree[expanded_tree[,2]==cluster,1]] <- 1
    keycluster        <- cluster
    old.tree          <- graph_pyclone$Corrected_tree

    k <- 1
    for (cluster1 in nesting.to.remove)
    {
      cat('\nmoving cluster:' ,cluster1)
      previous.tree    <- old.tree
      old.tree.expanded <- expand.tree(previous.tree,trunk_cluster)

      new.tree         <- old.tree
      nestedclust_test <- nestedclust
      nestedclust_test[keycluster,cluster1] <- 0
      nesting.removed  <- paste(cluster1,keycluster,sep=":")
      nesting.to.keep  <- cbind(colnames(nestedclust_test)[nestedclust_test[keycluster,]==1],keycluster)
      nesting.to.keep  <- nesting.to.keep[paste(nesting.to.keep[,1],nesting.to.keep[,2],sep=":")%in%paste(old.tree.expanded[,1],old.tree.expanded[,2],sep=":"),]


      new.tree <- new.tree[!paste(new.tree[,1],new.tree[,2],sep=":")%in%nesting.removed,,drop=FALSE]

      new.tree           <- rbind(new.tree,nesting.to.keep)


      rownames(new.tree) <- 1:nrow(new.tree)
      new.tree            <- prune.tree_test(new.tree,nestedclust_test,trunk_cluster)
      new.tree           <- new.tree[!duplicated(paste(new.tree[,1],new.tree[,2],sep=":")),,drop=FALSE]
      new.tree.graph     <- graph.data.frame(new.tree)
      plot(new.tree.graph,layout=layout.reingold.tilford(new.tree.graph, root=trunk_cluster) )

      new.tree <- new.tree[!duplicated(new.tree[,2]),,drop=FALSE]


      #now, let's explore the tree associated with this.
      tlevels          <- check.levels.ccf(new.tree, ccf_ci_lower, trunk_cluster)
      level.issue      <- max(tlevels[[1]])>=max_per_level
      circles_exist    <- girth(new.tree.graph)$girth != 0
      circles_exist    <- ifelse(circles_exist==TRUE,TRUE,ifelse(max(table(new.tree[,2]))>1,TRUE,FALSE))
      new.tree.change  <- !identical(new.tree,previous.tree)

      if(level.issue%in%FALSE&circles_exist%in%FALSE&new.tree.change)
      {
        trees[["1_cluster_change"]][[keycluster]][[k]] <- list()
        trees[["1_cluster_change"]][[keycluster]][[k]][["tree"]]      <- new.tree
        trees[["1_cluster_change"]][[keycluster]][[k]][["maxccf"]]  <- max(tlevels[[1]])
        trees[["1_cluster_change"]][[keycluster]][[k]][["nesting"]]  <- nestedclust_test
        k <- k+1
        previous.tree <- new.tree

      }

      other.clusters <- nesting.to.remove[!nesting.to.remove%in%cluster1]
      for (cluster2 in other.clusters)
      {
        if(on.same.branch(cluster1,cluster2,graph_pyclone$Corrected_tree,trunk_cluster))
        {
          next;
        }
        if(!on.same.branch(cluster1,cluster2,graph_pyclone$Corrected_tree,trunk_cluster))
        {
          nestedclust_test[keycluster,cluster2] <- 0
          nesting.removed <- c(nesting.removed,paste(cluster2,keycluster,sep=":"))
        }
      }

      new.tree         <- old.tree
      nesting.to.keep  <- cbind(colnames(nestedclust_test)[nestedclust_test[keycluster,]==1],keycluster)
      new.tree <- new.tree[!paste(new.tree[,1],new.tree[,2],sep=":")%in%nesting.removed,,drop=FALSE]


      new.tree           <- rbind(new.tree,nesting.to.keep)
      rownames(new.tree) <- 1:nrow(new.tree)
      new.tree           <- prune.tree_test(new.tree,nestedclust_test,trunk_cluster)
      new.tree.graph     <- graph.data.frame(new.tree)
      plot(new.tree.graph,layout=layout.reingold.tilford(new.tree.graph, root=trunk_cluster) )


      #now, let's explore the tree associated with this.
      tlevels          <- check.levels.ccf(new.tree, ccf_ci_lower, trunk_cluster)
      level.issue      <- max(tlevels[[1]])>=max_per_level
      circles_exist    <- girth(new.tree.graph)$girth != 0
      circles_exist    <- ifelse(circles_exist==TRUE,TRUE,ifelse(max(table(new.tree[,2]))>1,TRUE,FALSE))
      new.tree.change  <- !identical(new.tree,previous.tree)


      if(level.issue%in%FALSE&circles_exist%in%FALSE&new.tree.change)
      {
        trees[["1_cluster_change"]][[keycluster]][[k]] <- list()
        trees[["1_cluster_change"]][[keycluster]][[k]][["tree"]]      <- new.tree
        trees[["1_cluster_change"]][[keycluster]][[k]][["maxccf"]]  <- max(tlevels[[1]])
        trees[["1_cluster_change"]][[keycluster]][[k]][["nesting"]]  <- nestedclust_test
        k <- k+1
      }

    }

  }

  if(length(potential.trees[[length(potential.trees)]])>=2)
  {
    cat('\n moving different clusters')
    for (change.type in 2:length(potential.trees[[length(potential.trees)]]))
    {
      cat('\nchecking', change.type, 'clusters')
      nameToUse <- paste(change.type,"_cluster_change",sep="")
      trees[[nameToUse]] <- list()
      clusters.to.test <- potential.trees[lapply(potential.trees,length)==change.type]
      #take first part of the cluster
      seq100 <- round(seq(from=1,length(clusters.to.test),length.out =10))
      names(seq100) <- seq(10,100,by=10)
      for (i in 1:length(clusters.to.test))
      {
        #cat(i)
        if(i%in%c(seq100))
        {
          cat('\n%complete:',names(seq100[seq100%in%i]))
        }
        theclusters <- clusters.to.test[[i]]
        cluster     <- theclusters[length(theclusters)]
        clustersdone <- PasteVector(theclusters[-length(theclusters)],sep=":")
        num_cluster_change  <- paste(length(theclusters[-length(theclusters)]),"_cluster_change",sep="")
        trees_to_add_to <- trees[[num_cluster_change]][[clustersdone]]
        if(length(trees_to_add_to)==0)
        {
          next;
        }

        for (i in 1:length(trees_to_add_to))
        {
          #now we can treat it like we did above for the single cases
          nestedclust_orig <- trees_to_add_to[[i]]$nesting
          nesting.to.remove <- names(which(c(nestedclust_orig[cluster,]-strict_nestedclust[cluster,])==1))
          nestedclust_test  <- nestedclust_orig
          keycluster        <- cluster
          old.tree          <- trees_to_add_to[[i]]$tree
          k <- 1
          for (cluster1 in nesting.to.remove)
          {
            new.tree         <- old.tree
            nestedclust_test <- nestedclust_orig
            nestedclust_test[keycluster,cluster1] <- 0
            nesting.removed  <- paste(cluster1,keycluster,sep=":")
            #which other clusters are internal
            other.clusters <- nesting.to.remove[!nesting.to.remove%in%cluster1]
            for (cluster2 in other.clusters)
            {
              if(on.same.branch(cluster1,cluster2,graph_pyclone$Corrected_tree,trunk_cluster))
              {
                next;
              }
              if(!on.same.branch(cluster1,cluster2,graph_pyclone$Corrected_tree,trunk_cluster))
              {
                nestedclust_test[keycluster,cluster2] <- 0
                nesting.removed <- c(nesting.removed,paste(cluster2,keycluster,sep=":"))
              }
            }

            nesting.to.keep <- cbind(colnames(nestedclust_test)[nestedclust_test[keycluster,]==1],keycluster)
            new.tree <- new.tree[!paste(new.tree[,1],new.tree[,2],sep=":")%in%nesting.removed,,drop=FALSE]

            new.tree           <- rbind(new.tree,nesting.to.keep)
            rownames(new.tree) <- 1:nrow(new.tree)
            new.tree           <- prune.tree_test(new.tree,nestedclust_test,trunk_cluster)
            new.tree           <- new.tree[!duplicated(paste(new.tree[,1],new.tree[,2],sep=":")),,drop=FALSE]
            new.tree.graph     <- graph.data.frame(new.tree)

            #now, let's explore the tree associated with this.
            tlevels          <- check.levels.ccf(new.tree, ccf_ci_lower, trunk_cluster)
            level.issue      <- max(tlevels[[1]])>=max_per_level
            circles_exist    <- girth(new.tree.graph)$girth != 0
            circles_exist    <- ifelse(circles_exist==TRUE,TRUE,ifelse(max(table(new.tree[,2]))>1,TRUE,FALSE))

            if(level.issue%in%FALSE&circles_exist%in%FALSE) # &nesting.issue%in%FALSE)
            {
              trees[[nameToUse]][[PasteVector(c(clustersdone,keycluster),sep=":")]][[k]] <- list()
              trees[[nameToUse]][[PasteVector(c(clustersdone,keycluster),sep=":")]][[k]][["tree"]]      <- new.tree
              trees[[nameToUse]][[PasteVector(c(clustersdone,keycluster),sep=":")]][[k]][["maxccf"]]  <- max(tlevels[[1]])
              trees[[nameToUse]][[PasteVector(c(clustersdone,keycluster),sep=":")]][[k]][["nesting"]]  <- nestedclust_test
              k <- k+1
            }
          }
        }
      }
    }
  }


  reference_ccf <- max(trees$original$maxccf,100)
  trees.out.good <- list()
  trees.out       <- list()

  trees.out[[1]] <- trees$original$tree
  trees.out.good[[1]] <- trees$original$tree
  treecode            <- PasteVector(sort(paste(trees$original$tree[,1],trees$original$tree[,2],sep=":")),sep=":")

  k <- 2
  p <- 2
  for (i in names(trees)[2:length(trees)])
  {
    for (j in names(trees[[i]]))
    {
      for (l in 1:length(trees[[i]][[j]]))
      {


        if(as.numeric(trees[[i]][[j]][[l]]$maxccf)<=reference_ccf)
        {
          treeToAdd <- trees[[i]][[j]][[l]]$tree
          treeToAdd <- PasteVector(sort(paste(treeToAdd[,1],treeToAdd[,2],sep=":")),sep=":")
          if(treeToAdd%in%treecode)
          {
            next;
          }
          if(!treeToAdd%in%treecode)
          {
            trees.out.good[[k]] <-trees[[i]][[j]][[l]]$tree
            trees.out[[p]] <- trees[[i]][[j]][[l]]$tree
            p <- p+1
            k <- k+1
            c <- c(treecode,treeToAdd)

          }
        }
        if(as.numeric(trees[[i]][[j]][[l]]$maxccf)>reference_ccf)
        {
          treeToAdd <- trees[[i]][[j]][[l]]$tree
          treeToAdd <- PasteVector(sort(paste(treeToAdd[,1],treeToAdd[,2],sep=":")),sep=":")
          if(treeToAdd%in%treecode)
          {
            next;
          }
          if(!treeToAdd%in%treecode)
          {
            trees.out[[p]] <- trees[[i]][[j]][[l]]$tree
            p <- p+1
            c <- c(treecode,treeToAdd)
          }

        }


        if(p%in%seq(from = 1,to = 100000,by=1000))
        {
          cat('\npossible trees exceeds',p)
        }

      }
    }
  }

  #get concensus branches
  nr.trees <- length(trees.out.good)
  consensus.branches <- names(which(table(paste(do.call(rbind,trees.out.good)[,1],do.call(rbind,trees.out.good)[,2],sep=":"))/nr.trees==1))
  #add this so the stric_nestclust
  for (branch in consensus.branches)
  {
    tmp <- unlist(strsplit(branch,split=":"))
    strict_nestedclust[tmp[2],tmp[1]] <- 1

  }

  if(length(trees.out.good)<=1)
  {
    return(NULL)
  }

  #strict_nestedclust

  #let's get rid of trees that are worse than the original tree
  # okay, so we have some output
  output <- list()
  output$all.trees  <- trees.out
  output$good.trees <- trees.out.good
  output$tree.info  <- trees
  output$consensus.nestedclust<- strict_nestedclust
  output$consensus.branches<- consensus.branches
  return(output)



}


# cleaning functions ####
merge.clusters.full   <- function(test_pyclone
                                  , test_pyclone_absolute
                                  , nested_pyclone
                                  , nested_pyclone_absolute
                                  , min_ccf=0.05
                                  , p_value_cut=0.01
                                  , min_cluster_size=5
                                  , use_boot)
{

  merge_clusters_phylo    <- merge.clusters(nested_pyclone = nested_pyclone,pyclone = test_pyclone,min_ccf = min_ccf)
  merge_clusters_absolute <- merge.clusters(nested_pyclone = nested_pyclone_absolute,pyclone = test_pyclone_absolute,min_ccf = min_ccf)
  update_pyclone_table    <- FALSE
  merged_clusters         <- c()

  if(!is.na(merge_clusters_absolute[1]))
  {
    merge_clusters_absolute <- merge_clusters_absolute[merge_clusters_absolute[,3]%in%TRUE,,drop=FALSE]

    # this suggests that potentially too much copy number correction took place in j_cluster or i_cluster
    # so let's remove copy number correction for mutations in these clusters (and merge the clusters together)
    if(nrow(merge_clusters_absolute)>=1)
    {
      new.pyclone   <- test_pyclone
      for (i in 1:nrow(merge_clusters_absolute))
      {
        sig_p_value <- FALSE
        regions_phylo <- grep('PhyloCCF',colnames(test_pyclone),value=TRUE)
        for (region in regions_phylo)
        {
          cluster_A_pyclone <- test_pyclone[test_pyclone[,'PycloneCluster']%in%merge_clusters_absolute[i,1],,drop=FALSE]
          cluster_B_pyclone <- test_pyclone[test_pyclone[,'PycloneCluster']%in%merge_clusters_absolute[i,2],,drop=FALSE]
          cluster_A_region_phyloCCF  <- cluster_A_pyclone[,region]
          cluster_B_region_phyloCCF  <- cluster_B_pyclone[,region]

          cluster_A_region_mutCPN  <- cluster_A_pyclone[,gsub("PhyloCCF","mut.cpn.num",region)]
          cluster_B_region_mutCPN  <- cluster_B_pyclone[,gsub("PhyloCCF","mut.cpn.num",region)]

          cluster_A_region_phyloNEW  <- ifelse(cluster_A_region_phyloCCF>cluster_A_region_mutCPN,cluster_A_region_mutCPN,cluster_A_region_phyloCCF)
          cluster_B_region_phyloNEW  <- ifelse(cluster_B_region_phyloCCF>cluster_B_region_mutCPN,cluster_B_region_mutCPN,cluster_B_region_phyloCCF)

          wt <- suppressWarnings(wilcox.test(ifelse(cluster_A_region_phyloNEW<min_ccf,0,cluster_A_region_phyloNEW),ifelse(cluster_B_region_phyloNEW<min_ccf,0,cluster_B_region_phyloNEW)))
          if(is.na(wt$p.value))
          {
            wt$p.value <- 1
          }

          if(wt$p.value<p_value_cut)
          {
            sig_p_value <- TRUE
            next;
          }
          update_pyclone_table <- TRUE
          new.pyclone[rownames(cluster_A_pyclone),region] <- cluster_A_region_phyloNEW
          new.pyclone[rownames(cluster_B_pyclone),region] <- cluster_B_region_phyloNEW
        }

        if(!sig_p_value)
        {
          #lets merge the clusters
          cat('merging clusters:',merge_clusters_absolute[i,2], '-->',merge_clusters_absolute[i,1], '\n')
          merged_clusters <- rbind(merged_clusters,c(merge_clusters_absolute[i,2],'-->',merge_clusters_absolute[i,1]))
          new.pyclone[new.pyclone[,'PycloneCluster']%in%merge_clusters_absolute[i,2],'PycloneCluster'] <- as.numeric(merge_clusters_absolute[i,1])

        }

      }

      remove_rows <- which(duplicated(merged_clusters[,1]))
      if(length(remove_rows)!=0) {
        merged_clusters <- merged_clusters[-remove_rows,]
      }

    }

  }

  if(update_pyclone_table)
  {
    cat('\nupdating CCF values')
    test_pyclone <- new.pyclone
    nested_pyclone <- determine.cluster.nesting(pyclone = test_pyclone
                                                , prefix = prefix
                                                # , sameChrFrac = 0.8
                                                , min_cluster_size = max(c(2, min_cluster_size))
                                                , pval_cutoff = p_value_cut
                                                , use_boot =use_boot
    )

  }

  if(length(merged_clusters)==0)
  {
    merged_clusters <- NA
  }

  outList <- list()
  outList[["nested_pyclone"]] <- nested_pyclone
  outList[["test_pyclone"]]   <- test_pyclone
  outList[["mergedclusters"]] <- merged_clusters
  return(outList)


}

prep_pyclone_tree.fn <- function(graph_pyclone
                                 , nested_pyclone
                                 , correct_tree = NULL
                                 , merge_cluster = NULL
                                 , sampleid = substr(colnames(nested_pyclone$ccf_cluster_table)[1], 1, 8))
{
  auto_edgelength <- graph_pyclone$edgelength
  trunk <- graph_pyclone$trunk
  auto_mean_pyclone_ccf <- nested_pyclone$ccf_cluster_table
  auto_tree <- graph_pyclone$Corrected_tree
  tree <- edgelength <- ccf_cluster_table <- merged_clusters <- NA
  if (!is.null(correct_tree)) {
    if (sampleid %in% names(correct_tree)) {
      tree <- correct_tree[[sampleid]]
      edgelength <- auto_edgelength
      ccf_cluster_table <- auto_mean_pyclone_ccf
    }
    if (sampleid %in% names(merge_cluster)) {
      merged_clusters <- merge_cluster[[sampleid]]
      for (i in 1:length(merge_cluster[[sampleid]])) {
        edgelength <- edgelength
        names(edgelength)[names(edgelength) %in% as.character(merge_cluster[[sampleid]])[i]] <- names(merge_cluster[[sampleid]])[i]
        edgelength <- c(edgelength[!names(edgelength) %in% names(merge_cluster[[sampleid]])[i]], setNames(sum(edgelength[names(edgelength) %in% names(merge_cluster[[sampleid]])[i]]), names(merge_cluster[[sampleid]])[i]))
        edgelength <- edgelength[order(as.numeric(names(edgelength)))]
      }
      #When merging clusters, keep the CCF of the main cluster, to avoid low level CCF in regions otherwise at 0
      ccf_cluster_table <- auto_mean_pyclone_ccf[!rownames(auto_mean_pyclone_ccf) %in% merge_cluster[[sampleid]],, drop = F]
    }
  }
  output <- list('manual_tree' = tree, 'manual_edgelength' = edgelength, 'manual_mean_pyclone_ccf' = ccf_cluster_table, 'trunk' = trunk, 'auto_tree' = auto_tree, 'auto_edgelength' = auto_edgelength, 'auto_mean_pyclone_ccf' = auto_mean_pyclone_ccf, 'merged_clusters' = merged_clusters, 'multiple_trees' = NULL)
  return(output)
}

clean.noisy.clusters  <- function(pyclone,max.absent.prop=0.05,min_ccf=0.05)
{

  pyclone_clusters_clean <-  unique(pyclone[pyclone[,'CleanCluster']%in%1,'PycloneCluster'])
  pyclone_clusters_clean <-  sort(pyclone_clusters_clean)
  ccf_columns <- colnames(pyclone)[grepl('PhyloCCF', colnames(pyclone))]
  corrected_pyclone <- pyclone
  corrected_cluster <- c()

  for (py_cluster in pyclone_clusters_clean)
  {
    cluster_pyclone <- pyclone[pyclone[,'PycloneCluster']%in%py_cluster,,drop=FALSE]

    for (region in ccf_columns)
    {
      if((length(which(cluster_pyclone[,region]>min_ccf))/nrow(cluster_pyclone))>=max.absent.prop)
      {
        next;
      }
      if(sum(cluster_pyclone[,region])==0)
      {
        next;
      }
      print(py_cluster)
      print(region)
      corrected_pyclone[rownames(cluster_pyclone),region] <- 0
      corrected_cluster  <- rbind(corrected_cluster,c(region,py_cluster))

    }

  }

  if(length(corrected_cluster)==0)
  {
    corrected_cluster <- NA
  }
  outlist <- list()
  outlist$corrected_pyclone <- corrected_pyclone
  outlist$corrected_cluster <- corrected_cluster
  return(outlist)

}

#' Function to remove mutation clusters with genomically clustered mutations
remove_clustered_clones <- function(test.pyclone
                                    , clonal_cluster
                                    , p_value_cut = 0.01
                                    , clustering_estimate_cut = 2 )
{

  # determine the clones
  clones <- unique( test.pyclone[, "PycloneCluster" ] )
  clones <- sort(clones)

  # get clustering metrics on each clone
  clustering_data <- do.call( rbind, lapply( clones, function(clone) clusterDistributionAcrossGenome(clone, clonal_cluster, test.pyclone ) ))

  num_muts <- sapply( clones, function(clone) sum( test.pyclone[, "PycloneCluster"] %in% clone, na.rm = T ))

  clustering_data <- suppressWarnings( data.frame( PycloneCluster = clones,
                                                   is_clonal = clones == clonal_cluster,
                                                   num_muts = num_muts,
                                                   clustered.p = clustering_data[, 1 ],
                                                   clustered.or = clustering_data[, 2 ],
                                                   cluster.chr.num = clustering_data[, 3 ] ) )

  bad_clones <- clustering_data[ clustering_data$clustered.p  < p_value_cut &
                                   !is.na(clustering_data$clustered.p) & # NA for the clonal cluster and when only 1 mutation
                                   clustering_data$clustered.or > clustering_estimate_cut, "PycloneCluster" ]


  test.pyclone <- test.pyclone[ !test.pyclone[, "PycloneCluster" ] %in% bad_clones, ]

  return(test.pyclone)

}

#' Function to determine the distribution of location of mutations within a cluster
clusterDistributionAcrossGenome <- function(cluster
                                            , clonal_cluster
                                            , test.pyclone
                                            , iterations =10000
)
{

  clusterMuts    <- rownames(test.pyclone)[ test.pyclone[,"PycloneCluster"] == cluster ]

  if( is.na(cluster) | is.na(clonal_cluster) |
      cluster == clonal_cluster| length(clusterMuts) == 1 ) return( c(p_value = NA,
                                                                      estimate = NA,
                                                                      NumChromosomes_obs = NA) )

  # assess how many chromosomes cluster is on
  NumChromosomes_obs <- length(unique(do.call(rbind,strsplit( clusterMuts ,split=":"))[,2]))
  NumMutations_obs   <- length(clusterMuts)

  # use for clonal clusterfor background
  clonalMuts    <- rownames( test.pyclone[ test.pyclone[,"PycloneCluster"] == clonal_cluster ,] )
  AllChromosomes <- do.call(rbind,strsplit(clonalMuts ,split=":"))[,2]

  #can have subclones with more mutations than the clonal - in which case increase clonal dataset size
  while( length(AllChromosomes) < NumMutations_obs ) AllChromosomes <- c(AllChromosomes, AllChromosomes)


  #now sample the same number of mutations repeated
  k <- 1
  numChromosomes_sim <- c()
  while(k<iterations)
  {
    numChromosomes_sim <- c(numChromosomes_sim,length(unique(sample(AllChromosomes,replace = F,size=NumMutations_obs))))
    k <- k+1
  }

  numSimsFewer <- numChromosomes_sim[numChromosomes_sim<NumChromosomes_obs]
  p_value <- length(numSimsFewer)/iterations
  estimate <- mean(numChromosomes_sim) / NumChromosomes_obs

  if(length(p_value)==0 | p_value == 0)
  {
    p_value <- 1/iterations
  }

  return( c(p_value, estimate, NumChromosomes_obs) )
}


# accessory functions ####
# clustergraph ####
update.nested.cluster <- function(nestedclust,ccf_ci_lower,max_per_level,trunk_cluster)
{
  # this function changes to nesting is only enforced if pigeonhole principle is otherwise violated.
  # e.g. 20 < 30, thus cluster could be nested, but need not be.

  new_nestedclust <- nestedclust
  clusters.without.trunk <- colnames(new_nestedclust)[!colnames(new_nestedclust)%in%trunk_cluster]
  for (cluster in clusters.without.trunk)
  {
    which.clusters.nest <- which(new_nestedclust[,cluster]==1)
    for (nested.clust in which.clusters.nest)
    {
      if(!TRUE%in%c(ccf_ci_lower[nested.clust,] + ccf_ci_lower[cluster,]>max_per_level))
      {
        #can remove nesting
        new_nestedclust[nested.clust,cluster] <- 0
      }
    }
  }
  return(new_nestedclust)
}

#' Function to prune the ancestral graph
prune.tree <- function(edgelist
                       , nestedclust)
{

  # attempts to remove cycles, if allowed.
  edgelist_new <- edgelist
  rownames(edgelist_new) <- rownames(edgelist_new)<- 1:nrow(edgelist_new)
  trunk_node <- unique(edgelist_new[,1][!edgelist_new[,1] %in% edgelist_new[,2]])
  receiving_nodes <- table(edgelist_new[,2])
  #Nodes with multiple inputs
  receiving_nodes <- receiving_nodes[receiving_nodes > 1]
  finalOrder <- c()
  for (i in nrow(edgelist):1)
  {
    finalOrder <- c(finalOrder,edgelist[i,2],edgelist[i,1])
  }
  finalOrder <- finalOrder[!duplicated(finalOrder)]
  finalOrder <- finalOrder[finalOrder%in%names(receiving_nodes)]

  receiving_nodes <- receiving_nodes[order(as.numeric(names(receiving_nodes)),decreasing=T)]
  receiving_nodes <- receiving_nodes[finalOrder]

  remove_edges <- c()
  if(length(receiving_nodes) > 0){

    tmpNames <- names(receiving_nodes[order(colSums(nestedclust[,names(receiving_nodes),drop=FALSE]),decreasing=F)])
    receiving_nodes <- receiving_nodes[tmpNames]
    for(i in 1:length(receiving_nodes)){
      test_nets <- edgelist_new[edgelist_new[,2] %in% names(receiving_nodes)[i],,drop=F]
      #This trick ensures that if two nodes are both single transmitters (causing a circle) it is not resolved here, as that would be corrected without testing if it is correct
      tmp <- sort(table(edgelist_new[edgelist_new[,2] %in% test_nets[,1],2]),decreasing=T)
      transmitting_node <- names(tmp[tmp %in% tmp[1]])
      probably_remove_edges <- c(remove_edges,rownames(test_nets)[!test_nets[,1] %in% transmitting_node])
      probably_edgelist_new <- edgelist_new[!rownames(edgelist_new) %in% probably_remove_edges,,drop=F]


      #Test if all nodes are still descended from trunk, and trunk must still be present in the nodes
      if(trunk_node %in% c(probably_edgelist_new)){
        directedGraph <- graph.data.frame(probably_edgelist_new)
        if(all(shortest.paths(directedGraph,mode='in')[-1,trunk_node] != 'Inf')){
          #				if(length(get.shortest.paths(directedGraph,  trunk,receiving_nodes[i])[[1]][[1]]) >1){
          remove_edges<-probably_remove_edges
          edgelist_new<-probably_edgelist_new
        }
      }
    }
  }

  edgelistOut <-edgelist_new[!rownames(edgelist_new) %in% remove_edges,,drop=F]
  edgelist_new <- edgelistOut
  receiving_nodes <- table(edgelist_new[,2])
  #Nodes with multiple inputs
  receiving_nodes <- receiving_nodes[receiving_nodes > 1]
  receiving_nodes <- receiving_nodes[order(as.numeric(names(receiving_nodes)),decreasing=T)]

  if(length(receiving_nodes)>0)
  { tmpNames <- names(receiving_nodes[order(colSums(nestedclust[,names(receiving_nodes),drop=FALSE]),decreasing=F)])
  receiving_nodes <- receiving_nodes[tmpNames]

  for (i in 1:length(receiving_nodes))
  {
    tmpTest <- edgelist_new[edgelist_new[,2]%in%receiving_nodes[i],,drop=FALSE]
    test_nets <- edgelist_new[edgelist_new[,2] %in% names(receiving_nodes)[i],,drop=F]
    tmp <- sort(table(edgelist_new[edgelist_new[,2] %in% test_nets[,1],2]),decreasing=T)

    #are the events all descended?
    transmitting_node <- names(tmp[tmp %in% tmp[1]])
    if(nrow(edgelist_new[edgelist_new[,2]%in%transmitting_node,,drop=FALSE])>=length(transmitting_node))
    {#these must be on same branch, therefore remove the earlier branches
      if(length(transmitting_node)==1)
      {
        probably_remove_edges <- c(remove_edges,rownames(test_nets)[!test_nets[,1] %in% transmitting_node])
        probably_edgelist_new <- edgelist_new[!rownames(edgelist_new) %in% probably_remove_edges,,drop=F]
        #Test if all nodes are still descended from trunk, and trunk must still be present in the nodes
        if(trunk_node %in% c(probably_edgelist_new)){
          directedGraph <- graph.data.frame(probably_edgelist_new)
          if(all(shortest.paths(directedGraph,mode='in')[-1,trunk_node] != 'Inf')){
            #				if(length(get.shortest.paths(directedGraph,  trunk,receiving_nodes[i])[[1]][[1]]) >1){
            remove_edges<-probably_remove_edges
            edgelist_new<-probably_edgelist_new
          }
        }
      }

      if(length(transmitting_node)==2)
      {
        if(TRUE%in%c(edgelist_new[,1]%in%transmitting_node[1]&edgelist_new[,2]%in%transmitting_node[2]))
        {

          probably_remove_edges <- c(remove_edges,rownames(test_nets)[test_nets[,1]%in%transmitting_node[1]])
          probably_edgelist_new <- edgelist_new[!rownames(edgelist_new) %in% probably_remove_edges,,drop=F]
          if(trunk_node %in% c(probably_edgelist_new)){
            directedGraph <- graph.data.frame(probably_edgelist_new)
            if(all(shortest.paths(directedGraph,mode='in')[-1,trunk_node] != 'Inf')){
              #				if(length(get.shortest.paths(directedGraph,  trunk,receiving_nodes[i])[[1]][[1]]) >1){
              remove_edges<-probably_remove_edges
              edgelist_new<-probably_edgelist_new
            }
          }
        }

        if(TRUE%in%c(edgelist_new[,2]%in%transmitting_node[1]&edgelist_new[,1]%in%transmitting_node[2]))
        {

          probably_remove_edges <- c(remove_edges,rownames(test_nets)[test_nets[,1]%in%transmitting_node[2]])
          probably_edgelist_new <- edgelist_new[!rownames(edgelist_new) %in% probably_remove_edges,,drop=F]
          if(trunk_node %in% c(probably_edgelist_new)){
            directedGraph <- graph.data.frame(probably_edgelist_new)
            if(all(shortest.paths(directedGraph,mode='in')[-1,trunk_node] != 'Inf')){
              #				if(length(get.shortest.paths(directedGraph,  trunk,receiving_nodes[i])[[1]][[1]]) >1){
              remove_edges<-probably_remove_edges
              edgelist_new<-probably_edgelist_new
            }
          }
        }



      }



    }


  }


  }
  edgelist_new[order(as.numeric(edgelist_new[,1]),as.numeric(edgelist_new[,2]),decreasing=FALSE),]
  edgelist_new<-edgelist_new[!rownames(edgelist_new) %in% remove_edges,,drop=F]

  return(edgelist_new)
}

prune.tree_test <- function(edgelist_full,nestedclust,trunk_cluster)
{

  edgelist_new   <- c()
  all_nodes      <- unique(c(edgelist_full))
  all_nodes      <- all_nodes[all_nodes%in%colnames(nestedclust)]
  all_nodes_without_trunk <- all_nodes[!all_nodes%in%trunk_cluster]

  for (terminal_node in all_nodes_without_trunk)
  {

    nestedclust[,terminal_node]
    parental_nodes <- edgelist_full[edgelist_full[,2]%in%terminal_node,,drop=F]
    if(nrow(parental_nodes)==1)
    {
      edgelist_new <-rbind(edgelist_new,parental_nodes)
      next;
    }
    # sort clusters by number of mutations
    order_nodes    <- names(sort(colSums(nestedclust[parental_nodes[,1],parental_nodes[,1]]),decreasing=TRUE))
    #if single lineage then should be straightforward
    if_single_lineage <- matrix(0,nrow=length(order_nodes),ncol=length(order_nodes))
    colnames(if_single_lineage) <- order_nodes
    rownames(if_single_lineage) <- order_nodes
    observed_lineage            <- nestedclust[order_nodes,order_nodes]



    parentals <- observed_lineage



    for (i in 1:(ncol(if_single_lineage)-1))
    {
      if_single_lineage[(i+1):length(order_nodes),i] <- 1
    }

    if(identical(observed_lineage,if_single_lineage))
    {
      #can remove all but last descendent
      last_ancestor <- colnames(if_single_lineage)[ncol(if_single_lineage)]
      edgelist_new  <- rbind(edgelist_new,cbind(last_ancestor,terminal_node))
    }

    if(!identical(observed_lineage,if_single_lineage))
    {

      #work out which combination is needed
      essential_parents <- colnames(observed_lineage)[colSums(observed_lineage)==0]
      essential_parents <- unique(c(essential_parents,names(which(observed_lineage[essential_parents,]==0))))

      edgelist_new  <- rbind(edgelist_new,cbind(essential_parents,terminal_node))

    }
  }

  #next, let's check the next set of nodes
  return(edgelist_new)

}

expand.tree <- function(tree.table,trunk_cluster)
{

  non_trunk_clusters <- unique(c(tree.table))[!unique(c(tree.table))%in%trunk_cluster]
  full_tree <- c()
  for (cluster in non_trunk_clusters)
  {
    last_ancestor <- cluster
    while(last_ancestor!=trunk_cluster)
    {
      branch_to_add  <- tree.table[tree.table[,2]==last_ancestor,]
      full_tree      <- rbind(full_tree,branch_to_add)
      #what what last ancestor was
      last_ancestor <- branch_to_add[1]
      if(last_ancestor!=cluster)
      {
        branch_to_add <- c(last_ancestor,cluster)
        full_tree      <- rbind(full_tree,branch_to_add)

      }

    }

  }

  full_tree  <- full_tree[!duplicated(paste(full_tree[,1],full_tree[,2],sep=":")),,drop=FALSE]
  rownames(full_tree) <- 1:nrow(full_tree)
  return(full_tree)
}


get.direct.descendent <- function(cluster,edgeList)
{
  return(edgeList[edgeList[,1]==cluster,2])
}

unfold_tree <- function(edgelist, lower, trunk)
{
  tree_levels <- check.levels.ccf(edgelist, lower, trunk)[[2]]
  tree_values <- unlist(tree_levels)
  for (i in 1:length(tree_levels)) {
    names(tree_values)[tree_values %in% tree_levels[[i]]] <- i
  }
  tree_values <- tree_values[!duplicated(tree_values)]
  edgelist_new <- edgelist
  rownames(edgelist_new) <- rownames(edgelist) <- 1:nrow(edgelist_new)
  trunk_node <- unique(edgelist[, 1][!edgelist[, 1] %in% edgelist[, 2]])
  receiving_nodes <- table(edgelist[, 2])
  #Nodes with multiple inputs
  receiving_nodes <- receiving_nodes[receiving_nodes > 1]
  receiving_nodes <- receiving_nodes[order(as.numeric(names(receiving_nodes)), decreasing = T)]
  remove_edges <- c()
  cause_cicles <- c()
  if (length(receiving_nodes) > 0) {
    for (i in 1:length(receiving_nodes)) {
      test_nets <- edgelist_new[edgelist_new[, 2] %in% names(receiving_nodes)[i],, drop = F]
      #This trick ensures that if two nodes are both single transmitters (causing a circle) it is not resolved here, as that would be corrected without testing if it is correct
      receiving_level <- tree_values[tree_values %in% names(receiving_nodes)[i]]
      transmitting_level <- tree_values[tree_values %in% test_nets[, 1]]
      tmp <- transmitting_level[order(as.numeric(names(transmitting_level)), decreasing = T)]
      transmitting_node <- tmp[names(tmp) %in% names(tmp)[1]]
      if (any(transmitting_node %in% cause_cicles)) next
      if (length(transmitting_node) > 1) {
        cause_cicles <- c(cause_cicles, transmitting_node)
      }
      #Test if all nodes are still descended from trunk
      probably_remove_edges <- c(remove_edges, rownames(test_nets)[!test_nets[, 1] %in% transmitting_node])
      if (trunk %in% c(edgelist[!rownames(edgelist) %in% probably_remove_edges,, drop = F])) {
        directedGraph <- graph.data.frame(edgelist[!rownames(edgelist) %in% probably_remove_edges,, drop = F])
        if (all(shortest.paths(directedGraph, mode = 'in')[-1, trunk_node] != 'Inf')) {
          remove_edges <- c(remove_edges, rownames(test_nets)[!test_nets[, 1] %in% transmitting_node])
        }
        if (!all(shortest.paths(directedGraph)[-1, trunk_node] != 'Inf')) {
          cause_cicles <- c(cause_cicles, transmitting_node)
        }
      }
    }
  }
  edgelist <- edgelist[!rownames(edgelist) %in% remove_edges,, drop = F]

  return(list(edgelist, cause_cicles))
}

check.levels.ccf <- function(tree, lower, trunk)
{

  #Identify tree levels

  tlevels <- list(trunk)
  clusters <- unique(c(tree))
  clusters <- clusters[!clusters %in% trunk]
  count <- 1
  while (length(clusters) > 0) {
    count <- count + 1
    previous_level <- tlevels[[count - 1]]
    tlevels[[count]] <- unique(tree[tree[, 1] %in% previous_level, 2])
    clusters <- clusters[!clusters %in% unlist(tlevels)]
  }

  lower[lower > 100] <- 100
  if (class(lower) != "numeric") {
    level_ccf <- matrix(NA, nrow = length(tlevels), ncol = ncol(lower))
  } else {level_ccf <- matrix(NA, nrow = length(tlevels), ncol = 1)}
  for (i in 1:nrow(level_ccf)) {
    if (class(lower) != "numeric") {
      level_ccf[i, ] <- colSums(lower[rownames(lower) %in% tlevels[[i]],, drop = F])
    } else {level_ccf[i,1] <- lower[names(lower) %in% tlevels[[i]], drop = F]
    }}
  return(list(level_ccf, tlevels))
}

prune.cycle <- function(egdgelist,nestedclust,ccf_ci_lower,max_per_level,trunk_cluster)
{
  unfolded_tree             <- egdgelist
  clusters_with_multi_input <- sort(unique(unfolded_tree[duplicated(unfolded_tree[,2]),2]),decreasing=FALSE)
  strictclust               <- update.nested.cluster(nestedclust,ccf_ci_lower = ccf_ci_lower,max_per_level = max_per_level,trunk_cluster)
  nestedclust_new           <- nestedclust
  clusters_to_use           <- as.character(unique(c(unfolded_tree))[order(as.numeric(unique(c(unfolded_tree))),decreasing=FALSE)])

  for (clust in clusters_with_multi_input)
  {
    #perhaps should include full tree here
    observed_nesting <- unfolded_tree[unfolded_tree[,2]%in%clust,,drop=FALSE]

    strict_nesting   <- cbind(names(which(strictclust[clust,]=='1')),clust)
    unstrict_nesting <- cbind(names(which(nestedclust[clust,]=='1')),clust)
    unstrict_nesting <- unstrict_nesting[!paste(unstrict_nesting[,1],unstrict_nesting[,2],sep=":")%in%paste(strict_nesting[,1],strict_nesting[,2],sep=":"),,drop=FALSE]

    #now check whether the nesting observed in tree is strict or unstrict
    # is any of the nesting strictly necessary?
    if(TRUE%in%paste(observed_nesting[,1],observed_nesting[,2],sep=":")%in%c(paste(strict_nesting[,1],strict_nesting[,2],sep=":")))
    {
      consider_removing <- observed_nesting[!paste(observed_nesting[,1],observed_nesting[,2],sep=":")%in%c(paste(strict_nesting[,1],strict_nesting[,2],sep=":")),,drop=FALSE]
      for (i in 1:nrow(consider_removing))
      {
        ancestral_to_remove <- consider_removing[i,1]
        nestedclust_new[clust,ancestral_to_remove] <- 0
      }

    }

    if(!TRUE%in%paste(observed_nesting[,1],observed_nesting[,2],sep=":")%in%c(paste(strict_nesting[,1],strict_nesting[,2],sep=":")))
    {
      #can choose which nesting to use (choose larger cluster in terms of mutations)
      keepthis    <- observed_nesting[1,,drop=FALSE]
      removethis  <- observed_nesting[2:nrow(observed_nesting),,drop=FALSE]

      nestedclust_new[clust,] <- nestedclust_new[keepthis[1,1],]
      nestedclust_new[clust,keepthis[1,1]] <- 1
      for (ancestor in names(which(strictclust[clust,]==1)))
      {
        nestedclust_new[clust,ancestor] <- 1
      }

      for (i in 1:nrow(removethis))
      {
        ancestral_to_remove <- removethis[i,1]
        nestedclust_new[clust,ancestral_to_remove] <- 0
      }

    }


  }

  #create a new tree with this data
  nestedclust_test     <- nestedclust_new[rownames(nestedclust_new)%in%clusters_to_use,colnames(nestedclust_new)%in%clusters_to_use]
  directedGraph_input_full_test <- matrix(0, 0, 2)
  colsums <- colSums(nestedclust_test)
  rowsums <- rowSums(nestedclust_test)

  nestedclust_test[trunk_cluster, ] <- 0
  if (!all(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1)) {
    tmp <- rownames(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, ])[!nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1]
    #If any here, check if they are included in the trunk
    if (length(tmp) > 0) {
      for (i in 1:length(tmp)) {
        if (all(ccf_ci_upper[trunk_cluster, ] >= ccf_ci_lower[tmp[i], ])) {
          nestedclust_test[tmp[i], trunk_cluster] <- 1
        }
      }
    }
    tmp <- rownames(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, ])[!nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1]
    #If any remains
    if (length(tmp) > 0) {
      warning(paste('Cluster(s) ', tmp, ' is not descended from the trunk!!'))
      cluster_qc[tmp, 'NotDescendFromTrunk'] <- 1
      if(force_trunk) {
        nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] <- 1
      }
    }
  }
  for (i in colnames(nestedclust_test)) {
    tmp <- cbind(i, rownames(nestedclust_test)[nestedclust_test[, i] == 1])
    if (ncol(tmp) == 2) {
      directedGraph_input_full_test <- rbind(directedGraph_input_full_test, tmp)
    }
  }
  directedGraph_input_pruned <- prune.tree_test(edgelist = directedGraph_input_full_test,nestedclust = nestedclust_test,trunk_cluster = trunk_cluster)#,nestedclust = nestedclust)

  directedGraph_input_pruned_unfolded <- directedGraph_input_pruned

  plot(graph.data.frame(directedGraph_input_pruned_unfolded),layout=layout.reingold.tilford(graph.data.frame(directedGraph_input_pruned_unfolded), root=trunk_cluster) )
  tlevels <- check.levels.ccf(directedGraph_input_pruned_unfolded, ccf_ci_lower, trunk_cluster)
  test_ccf <- is.there.ccf.issue(nestedclust = nestedclust_new,directed_input_graph = directedGraph_input_pruned_unfolded,ccf_ci_lower = ccf_ci_lower,trunk_cluster = trunk_cluster,clusters_to_remove = c(),clusters_to_use = clusters_to_use,max_per_level = max_per_level)

  output <- list()
  output$directedGraph_input_pruned_unfolded <- directedGraph_input_pruned_unfolded
  output$nestedcluster                       <- nestedclust_new
  output$levelissue                          <- test_ccf$levelissue
  output$circle                              <- test_ccf$circle
  output$consistencyissue                    <- test_ccf$consistencyissue
  return(output)
}

#' Function to check whether there is a CCF level issue
is.there.ccf.issue <- function(nestedclust
                               ,directed_input_graph
                               ,ccf_ci_lower
                               ,trunk_cluster
                               ,clusters_to_remove
                               ,clusters_to_use
                               ,max_per_level=115)
{
  revised_nestedclust <- nestedclust
  revised_nestedclust <- revised_nestedclust[rownames(revised_nestedclust)%in%unique(c(directed_input_graph))[!unique(c(directed_input_graph))%in%clusters_to_remove],,drop=FALSE]
  revised_nestedclust <- revised_nestedclust[,colnames(revised_nestedclust)%in%unique(c(directed_input_graph))[!unique(c(directed_input_graph))%in%clusters_to_remove],drop=FALSE]
  revised_nestedclust <- revised_nestedclust[rownames(revised_nestedclust)%in%clusters_to_use,colnames(revised_nestedclust)%in%clusters_to_use]


  directedGraph_input_revised <- matrix(0, 0, 2)

  for (i in colnames(revised_nestedclust)) {
    tmp <- cbind(i, rownames(revised_nestedclust)[revised_nestedclust[, i] == 1])
    if (ncol(tmp) == 2) {
      directedGraph_input_revised <- rbind(directedGraph_input_revised, tmp)
    }
  }
  new_lower                           <- ccf_ci_lower[rownames(ccf_ci_lower) %in% unique(c(directedGraph_input_revised)),, drop = F]


  directedGraph_input_revised_pruned <- prune.tree_test(edgelist = directedGraph_input_revised,nestedclust = revised_nestedclust,trunk_cluster = trunk_cluster)#,nestedclust = nestedclust)
  directedGraph_input_revised_pruned_unfolded   <- directedGraph_input_revised_pruned
  plot(graph.data.frame(directedGraph_input_revised_pruned_unfolded),layout=layout.reingold.tilford(graph.data.frame(directedGraph_input_revised_pruned_unfolded), root=trunk_cluster) )
  newlevels                                     <- check.levels.ccf(directedGraph_input_revised_pruned_unfolded, new_lower, trunk_cluster)
  new.levels.with.issues                        <- which(rowSums(newlevels[[1]]>max_per_level)>0)
  directedGraph_revised                         <- graph.data.frame(directedGraph_input_revised_pruned_unfolded)
  circles_exist                                 <- girth(directedGraph_revised)$girth != 0
  circles_exist    <- ifelse(circles_exist==TRUE,TRUE,ifelse(max(table(directedGraph_input_revised_pruned_unfolded[,2]))>1,TRUE,FALSE))


  internal.consistency                          <- check.internally.consistent(tree_small = directedGraph_input_revised_pruned_unfolded,tree_full = directedGraph_input_revised,trunk_cluster)
  consistency.issue                             <- FALSE%in%c(internal.consistency[,2])

  if(length(new.levels.with.issues)>0)
  {
    outList <- list()
    outList$directed_input_corrected <- directedGraph_input_revised_pruned_unfolded
    outList$tlevels                  <- newlevels
    outList$highestCCF               <- max(newlevels[[1]])
    outList$levelissue               <- TRUE
    outList$circle                   <- circles_exist
    outList$consistencyissue         <- consistency.issue
    return(outList)
  }
  if(length(new.levels.with.issues)==0)
  {
    outList <- list()
    outList$directed_input_corrected <- directedGraph_input_revised_pruned
    outList$tlevels                  <- newlevels
    outList$highestCCF               <- max(newlevels[[1]])
    outList$levelissue               <- FALSE
    outList$circle                   <- circles_exist
    outList$consistencyissue         <- consistency.issue
    return(outList)
  }

}


#' Function to test which clusters are best to remove from tree
#' This function permutes through clusters and checks which are the best to remove
#' from the tree - takes into account number of mutations
permute.clusters.to.remove <- function(test_out
                                       ,nestedclust
                                       ,max_per_level
                                       ,tlevels
                                       ,trunk_cluster
                                       ,cluster_qc
                                       ,ccf_ci_lower
                                       ,nclusters
                                       ,speed_cluster = 15
)
{

  test_out.input <- test_out
  clusters_to_use                           <- as.character(sort(as.numeric(unique(c(test_out$directed_input_corrected)))))

  # first let's remove bad cycles
  cluster_qc_new <- cluster_qc
  cluster_qc_new[names(table(test_out$directed_input_corrected[, 2])), 'IncomingEdges'] <- table(test_out$directed_input_corrected[, 2])
  cluster_qc_new[names(table(test_out$directed_input_corrected[, 1])), 'OutgoingEdges'] <- table(test_out$directed_input_corrected[, 1])
  cycles.to.test <- sort(as.numeric(rownames(cluster_qc_new)[cluster_qc_new[,'CauseCycle']==1]),decreasing=TRUE)
  clusters.to.remove <- c()
  directedGraph_input_corrected <- test_out$directed_input_corrected
  clusters_to_test   <- c()
  circle <- test_out$circle


  if(length(clusters_to_use)>speed_cluster)
  {
    k <- 1
    while (k<=length(cycles.to.test)&circle==TRUE)
    {
      clusters_to_test <- c(clusters_to_test,cycles.to.test[k])
      test.now         <- is.there.ccf.issue(nestedclust = nestedclust,directed_input_graph = directedGraph_input_corrected,ccf_ci_lower = ccf_ci_lower,trunk_cluster = trunk_cluster,clusters_to_remove = clusters_to_test,clusters_to_use = clusters_to_use,max_per_level = max_per_level)
      if(test.now$circle%in%FALSE)
      {
        clusters.to.remove <- c(clusters.to.remove,clusters_to_test)
        circle <- FALSE
        test_out <- test.now
      }
      k <- k+1
    }
  }


  clusters_to_use                           <- as.character(sort(as.numeric(unique(c(test_out$directed_input_corrected)))))
  clusters.decreasing.order.without.trunk   <- clusters_to_use[!clusters_to_use%in%trunk_cluster]


  ncolSize <- 10000

  if(length(clusters.decreasing.order.without.trunk)>40)
  {
    ncolSize <- 1000
  }

  order.matrx     <- matrix(NA,nrow =length(clusters.decreasing.order.without.trunk)-1
                            ,ncol=ncolSize)




  combnList <- list()
  for (i in 1:min((length(clusters.decreasing.order.without.trunk)-1),25))
  {
    # cat(i)
    tmp <- combn(rev(clusters.decreasing.order.without.trunk),i)
    #combinations(length(clusters.decreasing.order.without.trunk),i,clusters.decreasing.order.without.trunk,repeats.allowed = FALSE)

    combnList[[i]] <- tmp
    for (k in 1:min(ncol(tmp),ncolSize))
    {
      order.matrx[i,k] <- sum(nclusters[tmp[,k]])
    }
    rm(tmp)

  }

  muts_to_remove       <- unique(sort(c(order.matrx),decreasing = FALSE))
  level_issue          <- TRUE%in%c(tlevels[[1]]>max_per_level)
  comb                 <- 1
  clusters.to.remove   <- c()
  consistency.issue    <- test_out.input$consistencyissue

  # testing to remove mutations
  cat('\n\nTesting mutation removal:')
  while ((level_issue|circle) & comb<min(10000,length(muts_to_remove)) )
  {
    cat('\n')
    cat(paste(muts_to_remove[comb],"/", sum(nclusters)," mutations (",signif(muts_to_remove[comb]/sum(nclusters)*100,3),"%)",sep=""))


    toTest <- which(order.matrx == muts_to_remove[comb], arr.ind = TRUE)

    if(nrow(toTest)>1)
    {
      for (row in 1:nrow(toTest))
      {
        clusters_to_test <- combnList[[toTest[row,1]]][,toTest[row,2]]
        test_out         <- is.there.ccf.issue(nestedclust = nestedclust
                                               ,directed_input_graph = directedGraph_input_corrected
                                               ,ccf_ci_lower = ccf_ci_lower
                                               ,trunk_cluster = trunk_cluster
                                               ,clusters_to_remove = clusters_to_test
                                               ,clusters_to_use = clusters_to_use
                                               ,max_per_level = max_per_level
        )
        if(test_out$circle%in%FALSE&test_out$levelissue%in%FALSE&test_out$consistencyissue%in%FALSE)
        {
          level_issue <- FALSE
          consistency.issue <- FALSE
          circle     <- FALSE
          clusters.to.remove <- clusters_to_test

          cat('\n----------------')
          cat('\nRemoving the following clusters:\n')
          cat(clusters_to_test)
          cat('\nTotal mutations removed: ')
          cat(sum(nclusters[clusters_to_test]))
          cat(' ')
          cat(' ')
          break;
        }
      }

    }

    if(nrow(toTest)==1)
    {

      clusters_to_test <- combnList[[toTest[1]]][,toTest[2]]
      test_out         <- is.there.ccf.issue(nestedclust = nestedclust,directed_input_graph = directedGraph_input_corrected,ccf_ci_lower = ccf_ci_lower,trunk_cluster = trunk_cluster,clusters_to_remove = clusters_to_test,clusters_to_use = clusters_to_use,max_per_level = max_per_level)
      if(test_out$circle%in%FALSE&test_out$levelissue%in%FALSE&test_out$consistencyissue%in%FALSE)
      {
        level_issue <- FALSE
        consistency.issue <- FALSE
        circle <- FALSE
        clusters.to.remove <- clusters_to_test
        #break;
      }


    }



    comb <- comb+1
  }

  test_out$clusters.to.remove <- clusters.to.remove


  return(test_out)
}

permute.clusters.to.add.back <- function(new_test_out
                                         ,nestedclust
                                         ,max_per_level
                                         ,tlevels
                                         ,trunk_cluster
                                         ,cluster_qc
                                         ,ccf_ci_lower
                                         ,ccf_ci_upper
                                         ,new_clusters_to_test
)
{

  clusters.in.use.start  <- unique(c(new_test_out$directed_input_corrected))
  clusters.in.use.test   <- clusters.in.use.start
  clusters.in.use.finish <- clusters.in.use.start
  nestedclust_new        <- nestedclust
  # now go cluster by cluster and try to add
  for (cluster in new_clusters_to_test)
  {
    clusters.in.use.test <- c(clusters.in.use.finish,cluster)
    nestedclust_test     <- nestedclust_new[rownames(nestedclust_new)%in%clusters.in.use.test,colnames(nestedclust_new)%in%clusters.in.use.test]
    directedGraph_input_full_test <- matrix(0, 0, 2)
    colsums <- colSums(nestedclust_test)
    rowsums <- rowSums(nestedclust_test)

    nestedclust_test[trunk_cluster, ] <- 0
    if (!all(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1)) {
      tmp <- rownames(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, ])[!nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1]
      #If any here, check if they are included in the trunk
      if (length(tmp) > 0) {
        for (i in 1:length(tmp)) {
          if (all(ccf_ci_upper[trunk_cluster, ] >= ccf_ci_lower[tmp[i], ])) {
            nestedclust_test[tmp[i], trunk_cluster] <- 1
          }
        }
      }
      tmp <- rownames(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, ])[!nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1]
      #If any remains
      if (length(tmp) > 0) {
        warning(paste('Cluster(s) ', tmp, ' is not descended from the trunk!!'))
        cluster_qc[tmp, 'NotDescendFromTrunk'] <- 1
        if(force_trunk) {
          nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] <- 1
        }
      }
    }
    for (i in colnames(nestedclust_test)) {
      tmp <- cbind(i, rownames(nestedclust_test)[nestedclust_test[, i] == 1])
      if (ncol(tmp) == 2) {
        directedGraph_input_full_test <- rbind(directedGraph_input_full_test, tmp)
      }
    }
    directedGraph_input_pruned <- prune.tree_test(edgelist = directedGraph_input_full_test,nestedclust = nestedclust_test,trunk_cluster = trunk_cluster)

    directedGraph_input_pruned_unfolded <-directedGraph_input_pruned

    plot(graph.data.frame(directedGraph_input_pruned_unfolded),layout=layout.reingold.tilford(graph.data.frame(directedGraph_input_pruned_unfolded), root=trunk_cluster) )
    ccf.test                         <- is.there.ccf.issue(nestedclust = nestedclust_test,directed_input_graph = directedGraph_input_pruned_unfolded,ccf_ci_lower = ccf_ci_lower,trunk_cluster = trunk_cluster,clusters_to_remove = c(),clusters_to_use = clusters.in.use.test,max_per_level = max_per_level)
    if(ccf.test$levelissue%in%FALSE&ccf.test$circle%in%FALSE)#&ccf.test$consistencyissue%in%FALSE)
    {
      print(cluster)
      clusters.in.use.finish <- clusters.in.use.test
      next;
    }
    if(ccf.test$levelissue%in%TRUE|ccf.test$circle%in%TRUE)
    {
      prune.out <- prune.cycle(egdgelist = directedGraph_input_pruned_unfolded,nestedclust = nestedclust_new,ccf_ci_lower = ccf_ci_lower,max_per_level = max_per_level,trunk_cluster = trunk_cluster)
      if(prune.out$levelissue%in%FALSE&prune.out$circle%in%FALSE&prune.out$consistencyissue%in%FALSE)
      {
        nestedclust_new <- prune.out$nestedcluster
        clusters.in.use.finish <- clusters.in.use.test
        print(cluster)
      }
    }





  }

  #now create the final tree output with these all added.
  clusters.in.use.test <- c(clusters.in.use.finish)
  nestedclust_test     <- nestedclust_new[rownames(nestedclust_new)%in%clusters.in.use.test,colnames(nestedclust_new)%in%clusters.in.use.test]
  directedGraph_input_full_test <- matrix(0, 0, 2)
  colsums <- colSums(nestedclust_test)
  rowsums <- rowSums(nestedclust_test)

  nestedclust_test[trunk_cluster, ] <- 0
  if (!all(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1)) {
    tmp <- rownames(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, ])[!nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1]
    #If any here, check if they are included in the trunk
    if (length(tmp) > 0) {
      for (i in 1:length(tmp)) {
        if (all(ccf_ci_upper[trunk_cluster, ] >= ccf_ci_lower[tmp[i], ])) {
          nestedclust_test[tmp[i], trunk_cluster] <- 1
        }
      }
    }
    tmp <- rownames(nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, ])[!nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] == 1]
    #If any remains
    if (length(tmp) > 0) {
      warning(paste('Cluster(s) ', tmp, ' is not descended from the trunk!!'))
      cluster_qc[tmp, 'NotDescendFromTrunk'] <- 1
      if(force_trunk) {
        nestedclust_test[!rownames(nestedclust_test) %in% trunk_cluster, trunk_cluster] <- 1
      }
    }
  }
  for (i in colnames(nestedclust_test)) {
    tmp <- cbind(i, rownames(nestedclust_test)[nestedclust_test[, i] == 1])
    if (ncol(tmp) == 2) {
      directedGraph_input_full_test <- rbind(directedGraph_input_full_test, tmp)
    }
  }
  directedGraph_input_pruned <- prune.tree_test(edgelist = directedGraph_input_full_test,nestedclust = nestedclust_test,trunk_cluster = trunk_cluster)

  directedGraph_input_pruned_unfolded <- directedGraph_input_pruned#test[[1]]

  plot(graph.data.frame(directedGraph_input_pruned_unfolded),layout=layout.reingold.tilford(graph.data.frame(directedGraph_input_pruned_unfolded), root=trunk_cluster) )
  output                         <- is.there.ccf.issue(nestedclust = nestedclust,directed_input_graph = directedGraph_input_pruned_unfolded,ccf_ci_lower = ccf_ci_lower,trunk_cluster = trunk_cluster,clusters_to_remove = c(),clusters_to_use = clusters.in.use.test,max_per_level = max_per_level)
  clusters.to.remove             <- new_clusters_to_test[!new_clusters_to_test%in%c(output$directed_input_corrected)]
  output$clusters.to.remove      <- clusters.to.remove
  output$nestedclust_new         <- nestedclust_test
  return(output)


}

check.internally.consistent <- function(tree_small,tree_full,trunk)
{
  #Identify tree levels

  tlevels <- list(trunk)
  clusters <- unique(c(tree_small))
  clusters <- clusters[!clusters %in% trunk]
  count <- 1
  while (length(clusters) > 0) {
    count <- count + 1
    previous_level <- tlevels[[count - 1]]
    tlevels[[count]] <- unique(tree_small[tree_small[, 1] %in% previous_level, 2])
    clusters <- clusters[!clusters %in% unlist(tlevels)]
  }

  clusters <- unique(c(tree_small))
  clusters <- clusters[!clusters %in% trunk]
  cluster_out <- c()
  for (cluster in clusters)
  {
    tree_level_cluster <- nrow(tree_full[tree_full[,2]%in%cluster,,drop=FALSE])
    if((tree_level_cluster+1)>length(tlevels))
    {
      out <- cbind(cluster,FALSE)
      cluster_out <- rbind(cluster_out,out)
      next;
    }

    out <- cbind(cluster,cluster%in%tlevels[[tree_level_cluster+1]]|cluster%in%tlevels[[tree_level_cluster]])
    cluster_out <- rbind(cluster_out,out)
  }

  colnames(cluster_out) <- c('cluster','internally_consistent')
  return(cluster_out)
}

prune.tree.old  <- function(edgelist)
{
  edgelist_new <- edgelist
  rownames(edgelist_new) <- rownames(edgelist) <- 1:nrow(edgelist_new)
  trunk_node <- unique(edgelist[, 1][!edgelist[, 1] %in% edgelist[, 2]])
  receiving_nodes <- table(edgelist[, 2])
  #Nodes with multiple inputs
  receiving_nodes <- receiving_nodes[receiving_nodes > 1]
  receiving_nodes <- receiving_nodes[order(as.numeric(names(receiving_nodes)), decreasing = T)]
  remove_edges <- c()
  if (length(receiving_nodes) > 0) {
    for (i in 1:length(receiving_nodes)) {
      test_nets <- edgelist_new[edgelist_new[, 2] %in% names(receiving_nodes)[i],, drop = F]
      #This trick ensures that if two nodes are both single transmitters (causing a circle) it is not resolved here, as that would be corrected without testing if it is correct
      tmp <- sort(table(edgelist_new[edgelist_new[, 2] %in% test_nets[, 1], 2]), decreasing = T)
      transmitting_node <- names(tmp[tmp %in% tmp[1]])
      probably_remove_edges <- c(remove_edges, rownames(test_nets)[!test_nets[, 1] %in% transmitting_node])
      probably_edgelist_new <- edgelist_new[!rownames(edgelist_new) %in% probably_remove_edges,, drop = F]
      #Test if all nodes are still descended from trunk, and trunk must still be present in the nodes
      if (trunk_node %in% c(probably_edgelist_new)) {
        directedGraph <- graph.data.frame(probably_edgelist_new)
        if (all(shortest.paths(directedGraph, mode = 'in')[-1, trunk_node] != 'Inf')) {
          remove_edges <- probably_remove_edges
          edgelist_new <- probably_edgelist_new
        }
      }
    }
  }
  edgelist <- edgelist[!rownames(edgelist) %in% remove_edges,, drop = F]
  return(edgelist)
}

on.same.branch   <- function(cluster1,cluster2,tree.to.test,trunk_cluster)
{
  ancestors         <- c()
  previous.ancestor <- tree.to.test[tree.to.test[,2]==cluster1,1]
  ancestors         <- previous.ancestor
  while(previous.ancestor!=trunk_cluster)
  {
    previous.ancestor <- tree.to.test[tree.to.test[,2]==previous.ancestor,1]
    ancestors         <- c(ancestors,previous.ancestor)
  }
  if(cluster2%in%ancestors)
  {
    return(TRUE)
  }
  if(!cluster2%in%ancestors)
  {
    return(FALSE)
  }

}

check.nesting    <- function (tree.to.test,strict_nestedclust,trunk_cluster)
{
  directedGraph_input_full <- matrix(0, 0, 2)
  colsums <- colSums(strict_nestedclust)
  rowsums <- rowSums(strict_nestedclust)
  for (i in colnames(strict_nestedclust)) {
    tmp <- cbind(i, rownames(strict_nestedclust)[strict_nestedclust[, i] == 1])
    if (ncol(tmp) == 2) {
      directedGraph_input_full <- rbind(directedGraph_input_full, tmp)
    }
  }

  tree.test        <- rbind(directedGraph_input_full,tree.to.test)
  full.tree        <- tree.test[!duplicated(paste(tree.test[,1],tree.test[,2],sep=":")),,drop=FALSE]
  new.tree.out     <- prune.tree_test(full.tree,strict_nestedclust,trunk_cluster)
  potential.issues <- new.tree.out[!paste(new.tree.out[,1],new.tree.out[,2],sep=":")%in%paste(tree.to.test[,1],tree.to.test[,2],sep=":"),,drop=FALSE]

  if(nrow(potential.issues)==0)
  {
    return(FALSE)
  }

  row.without.issue <- c()
  for (i in 1:nrow(potential.issues))
  {
    if(on.same.branch(potential.issues[i,2],potential.issues[i,1],tree.to.test,trunk_cluster))
    {
      row.without.issue <- c(row.without.issue,i)
    }
  }
  if(length(row.without.issue)>=1)
  {
    potential.issues <- potential.issues[-row.without.issue,,drop=FALSE]

  }
  if(nrow(potential.issues)>0)
  {
    return(TRUE)
  }
  if(nrow(potential.issues)==0)
  {
    return(FALSE)
  }
}


#' Function to compute confidence intervals of each cluster in each region
#'
#' This function takes in an R list containing mutation PhyloCCF in each region,
#' and mutation assignments to a cluster and computes bootstrapped
#' confidence intervals.
#' @param pyclone An R list object containing information about the PhyloCCF
#' of each mutation in each tumour region.
#' @returns An R list containing elements: 'ccf_cluster_table', 'mean_phylo_ccf',
#' 'median_pyclone_ccf', 'median_phylo_ccf', 'ccf_ci_upper', 'ccf_ci_lower',
#' 'ccf_ci_boot_upper', 'ccf_ci_boot_lower'

calc.pyclone.ci <- function(pyclone
                            , pyclust
                            , nclusters
                            , prefix = "LTX"
                            , lower_min = 0)
{
  tmp <- colnames(pyclone)[grepl(prefix, colnames(pyclone))]
  phyloCCF <- tmp[grepl('PhyloCCF$', tmp)]
  pycloneCCF <- tmp[grepl('PycloneCCF$', tmp)]

  nestedclust <- matrix(0, ncol = length(nclusters), nrow = length(nclusters))
  rownames(nestedclust) <- colnames(nestedclust) <- names(nclusters)
  #Calculate mean pycloneCCF
  ccf_cluster_table <- round(t(do.call(rbind, lapply(pycloneCCF, function(x) {
    tapply(pyclone[, x], factor(pyclone[, pyclust]), mean)}))) * 100)
  #Calculate mean phyloCCF
  mean_phylo_ccf <- round(t(do.call(rbind, lapply(phyloCCF, function(x) {
    tapply(pyclone[, x], factor(pyclone[, pyclust]), mean)}))) * 100)

  #Calculate median pycloneCCF
  median_pyclone_ccf <- round(t(do.call(rbind, lapply(pycloneCCF, function(x) {
    tapply(pyclone[, x], factor(pyclone[, pyclust]), median)}))) * 100)
  #Calculate median phyloCCF
  median_phylo_ccf <- round(t(do.call(rbind, lapply(phyloCCF, function(x) {
    tapply(pyclone[, x], factor(pyclone[, pyclust]), median)}))) * 100)

  #95% confidence interval
  ccf_error  <- t(do.call(rbind, lapply(phyloCCF, function(x) {
    a <- tapply(pyclone[, x] * 100, factor(pyclone[, pyclust]), mean)
    b <- tapply(pyclone[, x] * 100, factor(pyclone[, pyclust]), sd)
    n <- tapply(pyclone[, x] * 100, factor(pyclone[, pyclust]), length)
    error <- qnorm(0.975) * b / sqrt(n)
    error
  })))

  #95% confidence interval for the median value (using bootstrapping)
  ccf_ci_boot_lower <- ccf_ci_boot_upper <- matrix(NA,nrow = nrow(ccf_error),ncol=ncol(ccf_error))
  colnames(ccf_ci_boot_lower) <- colnames(ccf_ci_boot_upper) <- colnames(ccf_error)
  rownames(ccf_ci_boot_lower) <- rownames(ccf_ci_boot_upper) <- rownames(ccf_error)

  # Set the number of bootstrapping replicates > total_tmb, but minimally 5000.
  total_tmb <- sum(nclusters)
  if (total_tmb < 5000){
    BS_replicates <- 5000
  } else {
    rounding_factor <- ceiling(total_tmb/1000)
    BS_replicates <- rounding_factor * 1000
  }
  print(paste0('Number of bootstrapping replicates: ', BS_replicates))

  for(col in 1:ncol(ccf_ci_boot_lower))
  {
    for(cl in rownames(ccf_ci_boot_lower))
    {
      b <- suppressWarnings(boot(data=pyclone[pyclone[,pyclust]==cl, phyloCCF[col]],statistic = function(x,i) mean(x[i]), R=BS_replicates))

      if(b$t0==1&length(unique(b$data))==1|b$t0==1&length(unique(b$t))==1|b$t0==0)
      {
        if(b$t0==0)
        {
          ccf_ci_boot_lower[cl,col]  <- 0
          ccf_ci_boot_upper[cl,col]  <- 0
        }
        if(b$t0==1&length(unique(b$data))==1|b$t0==1&length(unique(b$t))==1)
        {
          ccf_ci_boot_lower[cl,col]  <- 100
          ccf_ci_boot_upper[cl,col]  <- 100
        }

      }
      else if(sd(b$t)==0&b$t0!=0&b$t0!=1)
      {
        ccf_ci_boot_lower[cl,col]  <- b$t0*100
        ccf_ci_boot_upper[cl,col]  <- b$t0*100
      }
      else{
        b.ci <- suppressWarnings(boot.ci(b,conf = 0.99))
        ccf_ci_boot_lower[cl,col] <- b.ci$normal[2]*100
        ccf_ci_boot_upper[cl,col] <- b.ci$normal[3]*100
      }
    }
  }

  ccf_ci_upper <- mean_phylo_ccf + ccf_error
  ccf_ci_lower <- mean_phylo_ccf - ccf_error

  ccf_ci_lower[is.na(ccf_ci_lower)] <- mean_phylo_ccf[is.na(ccf_ci_lower)]
  ccf_ci_upper[is.na(ccf_ci_upper)] <- mean_phylo_ccf[is.na(ccf_ci_upper)]
  ccf_ci_lower[mean_phylo_ccf     < lower_min]     <- 0 # reset anything under lower_min to 0 (noise level)
  ccf_ci_upper[mean_phylo_ccf     < lower_min]     <- 0 # reset anything under lower_min to 0 (noise level)
  mean_phylo_ccf[mean_phylo_ccf < lower_min]     <- 0 # reset anything under lower_min to 0 (noise level)

  ccf_ci_boot_lower[is.na(ccf_ci_boot_lower)] <- mean_phylo_ccf[is.na(ccf_ci_boot_lower)]
  ccf_ci_boot_upper[is.na(ccf_ci_boot_upper)] <- mean_phylo_ccf[is.na(ccf_ci_boot_upper)]
  ccf_ci_boot_lower[mean_phylo_ccf     < lower_min]     <- 0 # reset anything under lower_min to 0 (noise level)
  ccf_ci_boot_upper[mean_phylo_ccf     < lower_min]     <- 0 # reset anything under lower_min to 0 (noise level)

  ccf_ci_boot_lower[ccf_ci_boot_lower<0] <- 0
  ccf_ci_lower[ccf_ci_lower<0] <- 0


  colnames(ccf_cluster_table) <- colnames(mean_phylo_ccf) <- colnames(median_phylo_ccf) <- colnames(median_pyclone_ccf) <- colnames(ccf_ci_upper) <- colnames(ccf_ci_lower) <- colnames(ccf_ci_boot_lower) <- colnames(ccf_ci_boot_upper) <- sub('_PhyloCCF', '', phyloCCF)

  outlist <- list(ccf_cluster_table, mean_phylo_ccf, median_pyclone_ccf, median_phylo_ccf, ccf_ci_upper,ccf_ci_lower,ccf_ci_boot_upper,ccf_ci_boot_lower)
  names(outlist) <- c('ccf_cluster_table', 'mean_phylo_ccf', 'median_pyclone_ccf', 'median_phylo_ccf', 'ccf_ci_upper', 'ccf_ci_lower','ccf_ci_boot_upper','ccf_ci_boot_lower')

  return(outlist)
}




#' Test difference in PhyloCCF distributions of each pair of mutation clusters
#'
#' This function compares the distributions of the PhyloCCF of each pair of
#' mutation clusters in the dataset and outputs a nesting matrix.
#' @param pyclone An R list object containing information about the PhyloCCF
#' of each mutation in each tumour region.
#' @param nclusters Number of clusters
#' @param pval_cutoff A p-value significance threshold for testing whether
#' clusters can be nested. (i.e. a p-value < pval_cutoff is significant)
#' @returns This function returns list of nesting matrices. Each element of the
#' list is a nesting matrix for one tumour region, that describes whether a
#' cluster A (row) can be nested within a cluster B (column).

test.distributions <- function(pyclone
                               , nclusters
                               , pval_cutoff = 0.05)
{
  # Use a statistical test to test for nesting.
  # Done in each region separately, then looking across all regions to create the final nesting matrix
  # Given that if we use phyloCCF can be greater than 1, let's treat as not significantly greater if it's already greater than 0.99
  # also have changed to a wilcox test, so a bit more robust to outliers. could also perform mad winsorization.

  ccf_columns <- colnames(pyclone)[grepl("PhyloCCF", colnames(pyclone))]

  outlist <- list()
  for (i in 1:length(ccf_columns)) {
    outmat <- matrix(1, ncol = length(nclusters), nrow = length(nclusters))
    rownames(outmat) <- colnames(outmat) <- names(nclusters)

    for (j in colnames(outmat)) {
      for (k in rownames(outmat)) {
        if (as.numeric(j) < as.numeric(k)) {
          indx <- pyclone[, 'PycloneCluster'] %in% c(j,k)
          outmat[k,j] <- suppressWarnings(wilcox.test(pyclone[indx, ccf_columns[i]] ~ pyclone[indx, 'PycloneCluster'], alternative = 'less'))$p.value
          # if significant make sure this isn't simply because both are greater than clonal
          if(is.na(outmat[k,j]))
          {
            outmat[k,j] <- 1
          }
          if(outmat[k,j]<pval_cutoff)
          {
            outmat[k,j] <- suppressWarnings(wilcox.test(ifelse(pyclone[indx, ccf_columns[i]]>1,1,pyclone[indx, ccf_columns[i]])~ pyclone[indx, 'PycloneCluster'], alternative = 'less'))$p.value
          }

        }
        if (as.numeric(j) > as.numeric(k)) {
          indx <- pyclone[, 'PycloneCluster'] %in% c(j,k)
          outmat[k,j] <- suppressWarnings(wilcox.test(pyclone[indx, ccf_columns[i]] ~ pyclone[indx, 'PycloneCluster'], alternative = 'greater'))$p.value
          if(is.na(outmat[k,j]))
          {
            outmat[k,j] <- 1
          }
          if(outmat[k,j]<pval_cutoff)
          {
            outmat[k,j] <- suppressWarnings(wilcox.test(ifelse(pyclone[indx, ccf_columns[i]]>1,1,pyclone[indx, ccf_columns[i]]) ~ pyclone[indx, 'PycloneCluster'], alternative = 'greater'))$p.value

          }
        }
      }
    }
    outmat[is.na(outmat)] <- 1
    outlist[[i]] <- outmat
  }

  # this matrix above explores whether each cluster is significantly greater (rows vs. columns) or significantly less (columns vs. rows)

  outmat <- matrix(0, ncol = length(nclusters), nrow = length(nclusters))
  rownames(outmat) <- colnames(outmat) <- names(nclusters)

  for (j in 1:ncol(outmat)) {
    tmp <- do.call(cbind, lapply(outlist, function(x) x[, j])) # which clusters is this cluster less than
    outmat[, j] <- c(0, 1)[match(apply(tmp, 1, function(x) all(x > pval_cutoff)), c('FALSE', 'TRUE'))]
    outmat[j, j] <- 0
  }

  # the structure of outmat is such that if you go column by column you can see which clusters (by row) could be its ancestor (above in tree, indicated by 0)
  # and which could be its descendents (below on tree, labelled as 1)

  outmat_new <- matrix(0, ncol = length(nclusters), nrow = length(nclusters))
  rownames(outmat_new) <- colnames(outmat_new) <- names(nclusters)

  for (j in 1:ncol(outmat)) {
    tmp         <- do.call(cbind, lapply(outlist, function(x) x[, j])) # which clusters is this cluster less than
    outmat[, j] <- c(0, 1)[match(apply(tmp, 1, function(x) all(x > pval_cutoff)), c('FALSE', 'TRUE'))]
    outmat[j, j] <- 0
  }

  #make sure there is no mirroring (this is biologically impossible)
  for (i in 1:nrow(outmat)) {
    for (j in 1:ncol(outmat)) {
      if(outmat[i,j]==1&outmat[j,i]==1)
      { if (j>i) {
        outmat[i,j] <- 0
      }}
    }
  }


  outlist <- c(list(outmat), outlist)
  names(outlist) <- c('nestedclust', sub(paste('_',"PhyloCCF",sep=""), '', ccf_columns))
  return(outlist)
}

test.distributions.with.clonality.table <- function(pyclone
                                                    , nclusters
                                                    , clonality_table
                                                    , pval_cutoff = 0.05)
{
  # Use a statistical test to test for nesting.
  # Done in each region separately, then looking across all regions to create the final nesting matrix
  # Given that if we use phyloCCF can be greater than 1, let's treat as not significantly greater if it's already greater than 0.99
  # also have changed to a wilcox test, so a bit more robust to outliers. could also perform mad winsorization.

  ccf_columns <- colnames(pyclone)[grepl("PhyloCCF", colnames(pyclone))]

  outlist <- list()
  for (i in 1:length(ccf_columns)) {
    outmat <- matrix(1, ncol = length(nclusters), nrow = length(nclusters))
    rownames(outmat) <- colnames(outmat) <- names(nclusters)

    for (j in colnames(outmat)) {
      for (k in rownames(outmat)) {
        if (as.numeric(j) < as.numeric(k)) {
          indx <- pyclone[, 'PycloneCluster'] %in% c(j,k)
          outmat[k,j] <- suppressWarnings(wilcox.test(pyclone[indx, ccf_columns[i]] ~ pyclone[indx, 'PycloneCluster'], alternative = 'less'))$p.value
          # if significant make sure this isn't simply because both are greater than clonal
          if(is.na(outmat[k,j]))
          {
            outmat[k,j] <- 1
          }
          if(outmat[k,j]<pval_cutoff)
          {
            outmat[k,j] <- suppressWarnings(wilcox.test(ifelse(pyclone[indx, ccf_columns[i]]>1,1,pyclone[indx, ccf_columns[i]])~ pyclone[indx, 'PycloneCluster'], alternative = 'less'))$p.value
            if(clonality_table[k,gsub(paste("_","PhyloCCF",sep=""),"",ccf_columns)[i]]%in%'clonal'&clonality_table[j,gsub(paste("_","PhyloCCF",sep=""),"",ccf_columns)[i]]%in%'clonal')

            {
              outmat[k,j] <- 1

            }

          }

        }
        if (as.numeric(j) > as.numeric(k)) {
          indx <- pyclone[, 'PycloneCluster'] %in% c(j,k)
          outmat[k,j] <- suppressWarnings(wilcox.test(pyclone[indx, ccf_columns[i]] ~ pyclone[indx, 'PycloneCluster'], alternative = 'greater'))$p.value
          if(is.na(outmat[k,j]))
          {
            outmat[k,j] <- 1
          }
          if(outmat[k,j]<pval_cutoff)
          {
            outmat[k,j] <- suppressWarnings(wilcox.test(ifelse(pyclone[indx, ccf_columns[i]]>1,1,pyclone[indx, ccf_columns[i]]) ~ pyclone[indx, 'PycloneCluster'], alternative = 'greater'))$p.value
            if(clonality_table[k,gsub(paste("_","PhyloCCF",sep=""),"",ccf_columns)[i]]%in%'clonal'&clonality_table[j,gsub(paste("_","PhyloCCF",sep=""),"",ccf_columns)[i]]%in%'clonal')
            {
              outmat[k,j] <- 1

            }

          }
        }
      }
    }
    outmat[is.na(outmat)] <- 1
    outlist[[i]] <- outmat
  }

  # this matrix above explores whether each cluster is significantly greater (rows vs. columns) or significantly less (columns vs. rows)

  outmat <- matrix(0, ncol = length(nclusters), nrow = length(nclusters))
  rownames(outmat) <- colnames(outmat) <- names(nclusters)

  for (j in 1:ncol(outmat)) {
    tmp <- do.call(cbind, lapply(outlist, function(x) x[, j])) # which clusters is this cluster less than
    outmat[, j] <- c(0, 1)[match(apply(tmp, 1, function(x) all(x > pval_cutoff)), c('FALSE', 'TRUE'))]
    outmat[j, j] <- 0
  }

  # the structure of outmat is such that if you go column by column you can see which clusters (by row) could be its ancestor (above in tree, indicated by 0)
  # and which could be its descendents (below on tree, labelled as 1)

  outmat_new <- matrix(0, ncol = length(nclusters), nrow = length(nclusters))
  rownames(outmat_new) <- colnames(outmat_new) <- names(nclusters)

  for (j in 1:ncol(outmat)) {
    tmp         <- do.call(cbind, lapply(outlist, function(x) x[, j])) # which clusters is this cluster less than
    outmat[, j] <- c(0, 1)[match(apply(tmp, 1, function(x) all(x > pval_cutoff)), c('FALSE', 'TRUE'))]
    outmat[j, j] <- 0
  }

  #make sure there is no mirroring (this is biologically impossible)
  for (i in 1:nrow(outmat)) {
    for (j in 1:ncol(outmat)) {
      if(outmat[i,j]==1&outmat[j,i]==1)
      { if (j>i) {
        outmat[i,j] <- 0
      }}
    }
  }


  outlist <- c(list(outmat), outlist)
  names(outlist) <- c('nestedclust', sub(paste('_',"PhyloCCF",sep=""), '', ccf_columns))
  return(outlist)
}


verify.nesting <- function(nested_clust
                               , pyclone_ci)
{
  #This function will test that the nesting derived from significance testing makes sense when very low frequency clusters are ignored
  nested_clust_original <- nested_clust
  for (i in colnames(nested_clust)) {
    for (j in rownames(nested_clust)) {
      if (nested_clust[j, i] == 1) {
        #Test if descending cluster also exists in the same region(s)
        if (!any(pyclone_ci[['ccf_ci_lower']][j, pyclone_ci[['ccf_ci_lower']][i, ] > 0] > 0)) {
          nested_clust[j, i] <- 0
        }
        #test if the descending cluster also exists in other regions
        if (any(pyclone_ci[['ccf_ci_lower']][j, pyclone_ci[['ccf_ci_upper']][i, ] <= 0] > 0)) {
          nested_clust[j, i] <- 0
        }
      }
    }
  }

  # just retaining cluster 1 in cases where no nesting remains possible
  if (all(nested_clust) == 0) {
    trunk_cluster <- which(rowSums(nested_clust_original)==0)
    nested_clust <- nested_clust_original[trunk_cluster,trunk_cluster,drop=FALSE] #check how this deals with more than one trunk in original tree
  }

  return(nested_clust)
}

identify.issue.clusters <- function(distributions
                                    , pyclone_ci
                                    , pval_cutoff = 0.025)
{
  #Identify clusters that are incompatible (one larger, then smaller than the other).

  outmat <- matrix(0, 0, 2)
  test_list <- distributions[-c(1)]

  #check if any pair of clusters cannot fit within each other but are also nested (ie, one bigger than the other in one region, reversed in another, and both clusters cannot coexist)
  for (i in 1:ncol(distributions$nestedclust)) {
    a <- do.call(cbind, lapply(test_list, function(x) x[, i]))
    b <- do.call(cbind, lapply(test_list, function(x) x[i, ]))
    for (j in 1:nrow(a)) {
      if (i == j) next
      if (all(a[j, ] > pval_cutoff) & all(b[j, ] > pval_cutoff)) {
        coexists <- colSums(pyclone_ci$ccf_ci_lower[c(i, j),, drop = F])
        if(any(coexists > max(pyclone_ci$ccf_ci_lower))) {
          outmat <- rbind(outmat, c(i, j))
        }
      }
    }
  }
  return(outmat)
}

merge.clusters <- function(nested_pyclone
                           , pyclone
                           , min_ccf = 0.05
                           , p_value_cut = 0.01)
{
  potential.merge <- c()
  for (i in 1:(ncol(nested_pyclone$nestedclust)-1))
  {
    for (j in (i+1):ncol(nested_pyclone$nestedclust))
    {
      i_cluster <- colnames(nested_pyclone$nestedclust)[i]
      j_cluster <- colnames(nested_pyclone$nestedclust)[j]

      i_cluster_nesting_column <- nested_pyclone$nestedclust[,i]
      i_cluster_nesting_column <- i_cluster_nesting_column[!names(i_cluster_nesting_column)%in%c(i,j)]
      j_cluster_nesting_column <- nested_pyclone$nestedclust[,j]
      j_cluster_nesting_column <- j_cluster_nesting_column[!names(j_cluster_nesting_column)%in%c(i,j)]



      i_cluster_nesting_row <- nested_pyclone$nestedclust[i,]
      i_cluster_nesting_row <- i_cluster_nesting_row[!names(i_cluster_nesting_row)%in%c(i,j)]
      j_cluster_nesting_row <- nested_pyclone$nestedclust[j,]
      j_cluster_nesting_row <- j_cluster_nesting_row[!names(j_cluster_nesting_row)%in%c(i,j)]


      if(identical(i_cluster_nesting_column,j_cluster_nesting_column)&identical(i_cluster_nesting_row,j_cluster_nesting_row))
      {
        out <- cbind(i_cluster,j_cluster)
        potential.merge <- rbind(potential.merge,out)
      }

    }
  }


  if(length(potential.merge)==0)
  {
    return(NA)
  }
  potential.merge <- cbind(potential.merge,TRUE)
  colnames(potential.merge)[3] <- 'potential_merge'
  for (i in 1:nrow(potential.merge))
  {

    pyclone_cluster_1 <- pyclone[pyclone[,'PycloneCluster']%in%potential.merge[i,1],]
    pyclone_cluster_2 <- pyclone[pyclone[,'PycloneCluster']%in%potential.merge[i,2],]

    pyclone_regions   <- grep("PhyloCCF",colnames(pyclone),value=TRUE)
    for (region in pyclone_regions)
    {
      wt <- suppressWarnings(wilcox.test(ifelse(pyclone_cluster_1[,region]<min_ccf,0,pyclone_cluster_1[,region]),ifelse(pyclone_cluster_2[,region]<min_ccf,0,pyclone_cluster_2[,region])))
      if(is.na(wt$p.value))
      {
        next;
      }

      if(wt$p.value<p_value_cut)
      {
        potential.merge[i,3] <- FALSE
      }
    }

  }

  return(potential.merge)
}

# additional functions ####
PasteVector <- function(v,sep="")
{

  vt <- v[1];
  if(length(v) > 1){
    for(g in 2:length(v)){
      vt <- paste(vt,v[g],sep=sep)

    }
  }
  vt <- paste(vt," EnD",sep="");
  out.v <- sub(" EnD","",vt);
  out.v <- sub("NA , ","",out.v);
  out.v <- sub(" , NA","",out.v);
  out.v <- sub(" , NA , "," , ",out.v);
  return(out.v);
}

trx_rename.fn <- function(trxid, trialID = 'LTX')
{
  new_trxid <- as.character(trxid)
  #Includes hospital id?
  hospital <- grepl('\\w_', substr(trxid[1], 1, 2))
  if (hospital) {
    new_trxid <- substring(new_trxid, 3)
  }
  new_trxid <- gsub(trialID, paste0(trialID, "0"), new_trxid)
  return(new_trxid)
}



# plotting functions ####

#' Plotting function to colour the nodes on the phylogenetic tree
#'
#' @param edgelength A named vector containing number of mutations of each cluster
color.tree <- function(edgelength, opacity = 255)
{
  suppressPackageStartupMessages(require(RColorBrewer))
  ncols <- length(edgelength)
  max.cols <- 12
  cols     <- paste(brewer.pal(min(max.cols,ncols), name = "Paired"), sep = "")
  cols     <- rep(cols, ceiling(ncols / max.cols))[1:ncols]
  #Make colors with more than 12 clusters more transparent
  nopaq <- round(opacity / ceiling(ncols / 12))
  opacity_val <- c(sapply(1:ceiling(ncols / 12), function(x) rep(opacity - nopaq * (x - 1), 12)))[1:ncols]
  cols.opac <- paste(cols, as.hexmode(opacity_val), sep = "")
  return(cols.opac)
}


#' Function to extract all daughter clones from a parent using a phylogenetic tree
#'
#' @param tree A phylogenetic tree matrix with two columns specifying
#' 'parent' (column 1) and child (column 2)
#' @param parent.clones The name of the parent clone(s) for which you wish to find all daughters
extract_daughters <- function(tree
                              , parent.clones)
{

  daughters <- tree[ tree[,1] %in% parent.clones, 2]
  repeat{
    all.daughters <- unique( c( daughters, tree[ tree[,1] %in% daughters, 2] ) )
    if( all(all.daughters %in% daughters ) ) break
    daughters <- all.daughters
  }

  return( daughters )

}

#' Function to take a list of trees and identify the consensus relationships
extract_consensus_relationships <- function(tree_list
                                            , output_as_table = FALSE )
{

  parent_clones <- lapply( tree_list, function( tree ) unique( tree[, 1] ) )

  # determine all nesting relationships - ie all clones which are subclonal
  # to every parent clone in the tree
  # eg this tree: 1 -> 2, 2 -> 3, 2 -> 4 (where parent -> daughter)
  # has these relationships: 1:2, 1:3, 1:4, 2:3, 2:4
  relationships <- lapply( 1:length(tree_list), function(tree_i){

    tree <- tree_list[[ tree_i ]]
    parents <- parent_clones[[ tree_i ]]

    rels <- unlist( lapply( parents, function(parent){

      daughters <- extract_daughters( tree, parent )

      return( paste( parent, daughters, sep = ':' ) )

    } ) )

    return(rels)

  } )

  # detemine relationhsips that are common across all trees
  all_rels <- unique(unlist( relationships ))
  rel_freq <- table(unlist( relationships ))
  tree_num <- length(tree_list)
  consensus_rels <- names(rel_freq)[ rel_freq == tree_num ]

  #option to output as a table in the future - this could be clearer
  if(output_as_table == TRUE){
    consensus_rels <- data.frame( parent = sapply(strsplit( consensus_rels, split = ':' ), function(x) x[1] ),
                                  subclone = sapply(strsplit( consensus_rels, split = ':' ), function(x) x[2] ) )
  }

  return( consensus_rels )

}


#' Function to compute the Sum Condition Error for a list of trees
#'
#' This function takes a list of phylogenetic tree structures and the PhyloCCF cluster table and
#' computes the sum condition error for each tree.
#' @param tree_list A list of tree matrices
#' @param ccf_cluster_table A matrix of mean PhyloCCF of each cluster in
#' each tumour region
#' @param trunk The name of the truncal cluster
#' @returns sce_vec, A named vector of the sum condition error (SCE) for each
#' tree structure in the input tree list
compute_sum_condition_error <- function(tree_list
                                        , ccf_cluster_table
                                        , trunk)
{

  sce_list <- lapply(1:length(tree_list), function(tree_i){

    tree <- tree_list[[tree_i]]

    clusters_in_tree <- unique(c(tree))
    ccf_cluster_table <- ccf_cluster_table[rownames(ccf_cluster_table) %in% clusters_in_tree, , drop = F]

    # extract all tree paths
    all_paths <- createAllPathsList(tree, trunk)
    all_paths <- lapply(all_paths, function(x) strsplit(x, ',')[[1]])
    leaves <- unlist(lapply(all_paths, function(p) p[length(p)]))
    parents <- clusters_in_tree[!clusters_in_tree %in% leaves]

    # for every parent node in the tree, compute the total ccf error from the sum of its daughters
    total_error <- 0
    for (node in parents){
      daughters <- tree[tree[, 1] == node, 2]
      ccf_node <- ccf_cluster_table[rownames(ccf_cluster_table) == node, , drop = F]
      ccf_daughters <- ccf_cluster_table[rownames(ccf_cluster_table) %in% daughters, , drop = F]
      sum_ccf_daughters <- colSums(ccf_daughters)
      ccf_difference <- ccf_node - sum_ccf_daughters
      ccf_error <- ccf_difference
      ccf_error[ccf_error > 0] <- 0
      ccf_error <- mean(ccf_error) # mean over tumour regions
      total_error = total_error + ccf_error
    }

    mean_error_tree <- total_error/length(parents) # mean over parent nodes
    full_error_tree <- mean_error_tree * -1 # (make it positive)

    return(full_error_tree)
  })
  sce_vec <- unlist(sce_list)
  names(sce_vec) <- 1:length(tree_list)

  return(sce_vec)
}


#' Function to compute the edge probability score for a list of trees
#'
#' This function takes a list of phylogenetic tree structures and the PhyloCCF cluster table and
#' computes the sum condition error for each tree.
#' @param tree_list A list of tree matrices
#' @param edgelength A named vector containing number of mutations of each cluster
#' @param trunk The name of the truncal cluster
#' @returns sce_vec, A named vector of the sum condition error (SCE) for each
#' tree structure in the input tree list
# Function to take a list of trees and edge lengths and edge probabilities
compute_tree_edge_probability <- function(tree_list
                                          , edgelength
                                          , trunk)
{
  print(packageVersion('data.table'))
  n_alt_trees <- length(tree_list)
  clusters_in_tree <- unique(as.numeric(tree_list[[1]]))

  # Create a data table of all edges from all alternative trees
  alt_trees_edges = rbindlist(lapply(seq(tree_list), function(i){
    tree = as.data.table(tree_list[[i]])
    colnames(tree) = c('parent_node', 'node')
    tree$tree_name = i
    return(tree)
  }))

  # get edge lengths
  edgelength = data.frame(edgelength)
  colnames(edgelength) = c('node', 'edgelength')
  edgelength = edgelength[edgelength$node %in% clusters_in_tree, , drop = F]

  # merge edge lengths with all alternative tree edges
  alt_trees_edges = merge.data.table(alt_trees_edges,
                                     edgelength,
                                     by = 'node',
                                     all.x = T,
                                     all.y = F, # note: this makes sure we remove the row for the trunk node
                                     allow.cartesian = T)
  alt_trees_edges <- as.data.table(alt_trees_edges)
  # extract edge level metrics
  alt_trees_edges[, edge_count := .N, by = c('parent_node', 'node')] # number of alternative trees edge is in
  alt_trees_edges[, edge_prevalence := edge_count/n_alt_trees] # fraction of alternative trees edge is in
  alt_trees_edges[, log_edge_prevalence_mut_power := edgelength * log(edge_prevalence)] # log probability of the edge prevalence

  # compute tree edge probability score
  alt_trees_edges[, edge_probability := sum(log_edge_prevalence_mut_power), by = 'tree_name']
  tree_score_df <- unique(alt_trees_edges[, c("tree_name", 'edge_probability')])

  # save named vector of final edge probability score
  edge_probability_score <- tree_score_df$edge_probability
  names(edge_probability_score) <- tree_score_df$tree_name

  return(edge_probability_score)
}


#' Function to compute the tree level of a cluster
#' @param tree_graph A matrix of a tree structure
#' @param cluster Name of a cluster for which you want to get the tree level
get_tree_level <- function(tree_graph
                           , cluster)
{
  # trunk has level 1
  if (length(unique(as.numeric(tree_graph))) == 1){
    return(1)
  } else{
    trunk <- unique(tree_graph[, 1][!tree_graph[, 1] %in% tree_graph[, 2]])
    if (cluster == trunk) return(1)
    else {
      clusters_in_tree <- unique(as.numeric(as.matrix(tree_graph)))
      colnames(tree_graph) <- c('Parent', 'Child')
      tree_graph <- apply(tree_graph, c(1, 2), as.numeric)

      tree_df <- as.data.frame(tree_graph)

      level <- 1
      if (cluster %in% clusters_in_tree){
        current_cluster <- as.numeric(cluster)
        while (current_cluster != trunk) {
          parent <- tree_df[tree_df$Child == current_cluster, 'Parent']
          level <- level + 1
          current_cluster <- parent
        }
        return(level)
      } else return(NA)
    }
  }
}


#' Function to compute clone proportions on a selected alternative tree structure
#' @param tree_list A list of tree matrices
#' @param ccf_cluster_table A matrix of mean PhyloCCF of each cluster in
#' each tumour region
#' @param clonality_table A matrix of clonality calls for each cluster in
#' each tumour region
#' @param trunk The name of the truncal cluster
#' @param force_clonal_100 A logical indicating whether to for clusters that are
#' 'clonal' in a region to have CCF==100
#' @param tree_id The tree index of the selected alternative tree for which you
#' want to compute the clone proportions
compute_subclone_proportions <- function(tree_list
                                         , ccf_cluster_table
                                         , clonality_table
                                         , trunk
                                         , force_clonal_100 = TRUE
                                         , tree_id = 1)
{
  options(stringsAsFactors = F)

  tree <- tree_list[[tree_id]]

  clusters_in_tree <- unique(c(as.matrix(tree)))

  # if tree is just clonal cluster
  if (length(clusters_in_tree) == 1){
    clone_proportion_table <- (ccf_cluster_table > 0)
  } else {

    colnames(tree) <- c('Parent', 'Child')
    tree <- as.data.frame(tree)
    region_IDs <- colnames(ccf_cluster_table)

    # only use clusters in tree
    ccf_cluster_table <- ccf_cluster_table[rownames(ccf_cluster_table) %in% clusters_in_tree, , drop = F]
    clonality_table <- clonality_table[rownames(clonality_table) %in% clusters_in_tree, , drop = F]

    # force all clonal ccfs == 100 if flag is TRUE
    if (force_clonal_100 == T){
      clonal.factor <- (clonality_table == 'clonal')
      ccf_cluster_table <- ccf_cluster_table*(1 - clonal.factor) + 100*clonal.factor
    }

    # set a cap on CCFs above 100
    ccf_cluster_table[ccf_cluster_table > 100] <- 100

    # assign clonal cluster to be == 100 ccf
    ccf_cluster_table[rownames(ccf_cluster_table) == trunk, ] <- 100

    ccf_cluster_df <- as.data.frame(ccf_cluster_table)
    ccf_cluster_df$cluster <- rownames(ccf_cluster_df)


    ##### MAIN #####
    # Create a clone proportions dataframe to populate:
    proportions_df <- as.data.frame(ccf_cluster_table)
    proportions_df$cluster <- rownames(proportions_df)
    proportions_df[, region_IDs] <- 0

    # Compute clone proportions for each region independently
    for (r in region_IDs){
      cols <- c(r, 'cluster')
      region_ccf <- ccf_cluster_df[, cols]
      colnames(region_ccf)[1] <- c('ccf')

      clusters_present <- region_ccf[region_ccf$ccf != 0, 'cluster']
      parents_present <- unique(tree$Parent[tree$Parent %in% clusters_present])

      # order the parent subclones by tree level
      parent_df <- data.frame(parent_node = parents_present)
      parent_df$level <- sapply(parent_df$parent_node, function(p) return(get_tree_level(as.matrix(tree), p)))
      setorder(parent_df, level)

      # For each parent node in tree: compute the difference in CCF between parent node and sum of its children
      for (p in parent_df$parent_node){

        # TOP DOWN APPRAOCH:
        # Fix parental node CCFs, scale children nodes accordingly:
        # i.e. if children CCFs sum to > parental CCF, rescale child CCFs to be proportional values of parent ccf (so that their CCFs sum to parent exactly)

        children_nodes <- tree[tree$Parent == p, 'Child']
        parent_ccf <- as.numeric(region_ccf[region_ccf$cluster == p, 'ccf'])
        sum_children_ccf <- sum(region_ccf[region_ccf$cluster %in% children_nodes, 'ccf'])

        if (sum_children_ccf > parent_ccf){
          parent_proportion <- 0
          region_ccf[region_ccf$cluster %in% children_nodes, 'ccf'] <- parent_ccf * region_ccf[region_ccf$cluster %in% children_nodes, 'ccf'] / sum(region_ccf[region_ccf$cluster %in% children_nodes, 'ccf'])
        } else {
          parent_proportion <- parent_ccf - sum_children_ccf
        }
        proportions_df[proportions_df$cluster == p, r] <- parent_proportion

        # if the clones are terminal, add them to proportions_df as well
        for (d in children_nodes){
          if (d %in% tree[, 'Child'] & !(d %in% tree[, 'Parent'])) {
            proportions_df[proportions_df$cluster == d,r] <- region_ccf[region_ccf$cluster == d, 'ccf']
          }
        }
      }
    }

    # fix final output matrix
    rownames(proportions_df) <- proportions_df$cluster
    clone_proportion_table <- as.matrix(proportions_df[, region_IDs])
  }
  return(clone_proportion_table)
}
