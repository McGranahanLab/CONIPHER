findSimpleClusters <- function(input_list, runType = "WES") {
    ######################################################################################################################################################
    ### summary of what's required for function
    ### phylo.region.list - output as part of the initial steps of pyclone script (this contains the CCF of each mutation in each region prior to cluster)
    ### mut.table - the mutation table for the patient
    ######################################################################################################################################################
    new.dir           <- input_list$new.dir
    phylo.region.list <- input_list$phylo.region.list
    mut.table         <- input_list$mut.table
    seg.mat.phylo     <- input_list$seg.mat.phylo
    patient           <- input_list$patient

    ### define the output
    SmallClusters    <- list()

    regionsToUse     <- names(phylo.region.list)
    CCFtable         <- c()

    sharedMutations  <- unlist(phylo.region.list[[1]]$mutation_id)
    for (region in regionsToUse) {
        sharedMutations <- intersect(sharedMutations, unlist(phylo.region.list[[region]]$mutation_id))
    }
    new.phylo.region.list <- lapply(phylo.region.list, function(x) {
        x$cov <- (unlist(x$var_counts) + unlist(x$ref_counts))
        x$var_count <- unlist(x$var_counts)
        x$VAF <- x$var_count / x$cov * 100
        x$presence <- ifelse(x$var_count >= 1 & x$VAF > 1, TRUE, FALSE)
        return(x)
    })
    for (region in regionsToUse) {
        CCFtable <- cbind(CCFtable, unlist(new.phylo.region.list[[region]][sharedMutations, "presence"]))
    }
    rownames(CCFtable) <- sharedMutations
    colnames(CCFtable) <- regionsToUse
    BinaryTable <- ifelse(CCFtable, 1, 0)
    MutClusters <- apply(BinaryTable, 1, PasteVector, sep = ":")

    ### filter out mutations that are in small clusters
    UniqCluster <- unique(MutClusters)  

    ### set Names for the UniqClusters
    names(UniqCluster) <- as.character(1:length(UniqCluster))
    ### reverse the names as well
    UniqClusterNumber  <- names(UniqCluster)
    names(UniqClusterNumber) <- UniqCluster
    MutClusterNum <- UniqClusterNumber[MutClusters]
    names(MutClusterNum) <- names(MutClusters)  

    ### check whether evidence for multi-modality
    ClusterEvidenceForMultiModal <- rep(0, length(UniqCluster))
    names(ClusterEvidenceForMultiModal) <- names(UniqCluster)

    for (cluster in names(UniqCluster)) {
        SmallClusters[[cluster]]$clusterID            <- cluster
        SmallClusters[[cluster]]$clusterBinary        <- UniqCluster[cluster]
        SmallClusters[[cluster]]$MutationsWithCluster <- names(MutClusterNum[MutClusterNum %in% cluster])

        ### next, make pyclone tables so we can also run pyclone on these samples. 
        ClusterMutationsIDs <- names(MutClusterNum[MutClusterNum %in% cluster])
        specClusterTable    <- BinaryTable[ClusterMutationsIDs, , drop = FALSE]
        RegionsInCluster    <- colnames(specClusterTable)[which(specClusterTable[1, ] == 1)]
        ClusterTable        <- mut.table[mut.table$mutation_id %in% ClusterMutationsIDs, , drop = FALSE] 
        SmallClusters[[cluster]]$RegionsInCluster     <- RegionsInCluster

        ### check which regions are in the cluster 
        for (region in RegionsInCluster) {

            region.mut.table <- ClusterTable
            region.seg.copy  <- seg.mat.phylo[seg.mat.phylo$SampleID %in% region, , drop = FALSE]
            pyclone.table    <- data.frame(t(sapply(1:nrow(region.mut.table), identify.subclonal.mut.copy.number.ascat, region.mut.table, region.seg.copy, region, patient)), stringsAsFactors = FALSE)

            na.mutations     <- pyclone.table[is.na(pyclone.table$minor_cn), , drop = FALSE]
            loss.mutations   <- pyclone.table[as.numeric(pyclone.table$major_cn) == 0 | (as.numeric(pyclone.table$var_counts) + as.numeric(pyclone.table$ref_counts) == 0), ]
            error.muts       <- rbind(na.mutations, loss.mutations)
            error.muts       <- unlist(na.mutations$mutation_id, loss.mutations$mutation_id)
            error.muts.table <- paste0(new.dir, "/", region, ".error.muts.tsv")

            if (runType == "sim") {
                #CMR for sim data, modify NA CN 
                print("Running CN updates for simulation data.")
                pyclone.table$major_cn[is.na(pyclone.table$minor_cn)] <- 1
                pyclone.table$major_raw[is.na(pyclone.table$minor_cn)] <- 1
                pyclone.table$minor_raw[is.na(pyclone.table$minor_cn)] <- 1
                pyclone.table$fracA[is.na(pyclone.table$minor_cn)] <- 1
                pyclone.table$nMaj_A[is.na(pyclone.table$minor_cn)] <- 1
                pyclone.table$nMin_A[is.na(pyclone.table$minor_cn)] <- 1
                pyclone.table$fracB[is.na(pyclone.table$minor_cn)] <- 0
                pyclone.table$nMaj_B[is.na(pyclone.table$minor_cn)] <- 1
                pyclone.table$nMin_B[is.na(pyclone.table$minor_cn)] <- 1
                pyclone.table$minor_cn[is.na(pyclone.table$minor_cn)] <- 1
            }

            ### a few sanity checks
            pyclone.table    <- pyclone.table[!is.na(pyclone.table$minor_cn), ]
            pyclone.table    <- pyclone.table[!is.na(pyclone.table$ref_counts), ]
            pyclone.table    <- pyclone.table[!duplicated(pyclone.table$mutation_id), ]
            pyclone.table    <- pyclone.table[as.numeric(pyclone.table$major_cn) >= 1, ]
            pyclone.table    <- pyclone.table[!is.na(pyclone.table$minor_cn), ]

            ### now, let's check what the cancer cell fraction estimates are for this region
            region.ccf               <- phylo.region.list[[region]]
            region.ccf               <- data.frame(region.ccf, stringsAsFactors = FALSE)
            rownames(region.ccf)     <- region.ccf$mutation_id
            tmp                      <- intersect(unlist(pyclone.table$mutation_id), unlist(region.ccf$mutation_id))
            rownames(pyclone.table)  <- pyclone.table$mutation_id
            pyclone.table            <- pyclone.table[tmp, , drop = FALSE]
            region.ccf               <- region.ccf[tmp, , drop = FALSE]

            if (nrow(pyclone.table) > 0) {
                tmp                      <- round(((unlist(pyclone.table$var_counts)) / (unlist(region.ccf$phyloCCF) / 2)) - unlist(pyclone.table$var_counts))
                tmp[is.na(tmp)]          <- unlist(pyclone.table$ref_counts[(is.na(tmp))])
                pyclone.table$ref_counts <- tmp
                pyclone.table$minor_cn   <- 0
                pyclone.table$major_cn   <- 2
                pyclone.table$ref_counts <- apply(cbind(pyclone.table$ref_counts, 2), 1, max)
            }
            
            SmallClusters[[cluster]]$PyCloneTables[[region]]$pyclone.table <- pyclone.table
        }
    }
    return(SmallClusters)
}

RunPyCloneWithSimpleClusters <- function(clusterName, patientID, seg.mat, SmallClusters, patientDirToUse = new.dir, yamlConfigLoc = template.config.yaml, pyclone.burnin = 1000, pyclone.seed = 1024, run.pyclone = TRUE, pyclone.module = "PyClone/0.12.3-foss-2016b-Python-2.7.12-tkinter") {
    PyClone <- "pyclone-vi"
    ### give a name to the sample
    PyCloneRunName <- paste0(patientID, "_cluster", clusterName)
    ### create a specific subDirectory for this analysis
    ClusterDir <- paste0(patientDirToUse, "/", PyCloneRunName, "/")
    if (!dir.exists(ClusterDir)) {
        dir.create(ClusterDir, recursive = TRUE)
    }

    ### make the mutation files for each region
    RegionsInClustering <- SmallClusters[[clusterName]]$RegionsInCluster
    for (region in RegionsInClustering) {
        pyclone.tsv   <- paste0(ClusterDir, "/", region, ".tsv")
        pyclone.table <- SmallClusters[[clusterName]]$PyCloneTables[[region]]$pyclone.table

        # KT: needs reformatting for pyclone-vi
        pyclone.table <- pyclone.table[, c("mutation_id", "region", "ref_counts", "var_counts", "normal_cn", "major_raw", "minor_raw")]
        colnames(pyclone.table) <- c("mutation_id", "sample_id", "ref_counts", "alt_counts", "normal_cn", "major_cn", "minor_cn")
        pyclone.table$sample_id <- as.character(pyclone.table$sample_id)

        # KT: let's also add rhe tumour_content column
        purities <- unique(seg.mat[, c("SampleID", "ACF")])
        colnames(purities) <- c("sample_id", "tumour_content")

        pyclone.table <- left_join(pyclone.table, purities)

        if (nrow(pyclone.table) == 1) {
            pyclone.out <- matrix(apply(pyclone.table, 2, as.character), nrow = 1)
            colnames(pyclone.out) <- colnames(pyclone.table)
        } else {
            pyclone.out   <- apply(pyclone.table, 2, as.character)
        }
        write.table(pyclone.out, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE, file = pyclone.tsv)

        ### Run PyClone build_mutations_file TSV_FILE where TSV_FILE is the input file you have created.
        

        #run PyClone-vi git command, lets run its with default parametes for now
        pyclone.yaml <- paste0(ClusterDir, "/", region, ".yaml")
        sample.results <- paste0(ClusterDir, "/", patientID, '.results.tsv')
        cmd <- paste0(PyClone
                , " fit"
                , " --in-file ", pyclone.tsv
                , " --out-file ", pyclone.yaml
                , " --num-clusters ", min(max(1, floor(length(SmallClusters[[clusterName]]$MutationsWithCluster) / 5)), 10)
                , " --density beta-binomial")
        cat('\n')

        cmd2 <- paste0(PyClone
                , "  write-results-file"
                , " --in-file ", pyclone.yaml
                , " --out-file ", sample.results)
        cat('\n')

        if (run.pyclone) {
            # run model fitting step
            cat(cmd)
            system(cmd)

            # run step to select best results
            cat(cmd2)
            system(cmd2)
            
            
            ### KT: don't think I need any of this for PyClone-vi
            ### AH edit change so it works with states separated over multiple lines as well
            #exclude state including g_v=AB from yaml file
            # yaml <- readLines(pyclone.yaml)
            # rm.indx <- grep("AB", yaml)
            # if (length(grep("prior_weight", grep("AB", yaml, value = TRUE))) > 0) {
            #     yaml <- yaml[-rm.indx]    
            # } else {
            #     yaml <- yaml[-c(rm.indx-2, rm.indx-1, rm.indx, rm.indx+1)]
            # } 
            # write.table(yaml, file = pyclone.yaml, col.names = FALSE, row.names = FALSE, quote = FALSE)


        }
    }

    # KT: also don't think I need any of this for pyclone-vi
    # pyclone.config.yaml <- paste0(ClusterDir, "/", PyCloneRunName, ".config.yaml")
    # pyclone.config      <- readLines(yamlConfigLoc)
    # start.samples       <- (grep("samples", pyclone.config) + 1)
    # end.samples         <- length(pyclone.config)

    # sample.lines        <- pyclone.config[start.samples:end.samples]
    # pyclone.config      <- pyclone.config[-c(start.samples:end.samples)]
    # pyclone.config      <- c(pyclone.config, "init_method: connected", "", "samples:")
    # pyclone.config      <- gsub("working.directory.location", ClusterDir, pyclone.config)

    # write.table(pyclone.config, file = pyclone.config.yaml, col.names = FALSE, row.names = FALSE, quote = FALSE)

    # RegionsInClustering <- SmallClusters[[clusterName]]$RegionsInCluster
    # for (region in RegionsInClustering) {
    #     sample.config <- gsub("TCGA.barcode", region, sample.lines)
    #     # pyclone.yaml  <- paste0(new.dir, "/", region, ".yaml")
    #     pyclone.yaml  <- paste0(region, ".yaml")
    #     sample.config <- gsub("mutations.yaml", pyclone.yaml, sample.config)
    #     region.purity <- 0.5

    #     sample.config <- gsub("value: 1.0", paste0("value: ", signif(region.purity, 3)), sample.config)
    #     sample.config <- sample.config[1:8]

    #     if (run.pyclone) {
    #         write.table(sample.config, file = pyclone.config.yaml, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
    #     }
    # }

    # ### next, run pyclone
    # cmd <- paste0(PyClone
    #             , " run_analysis --config_file "
    #             , pyclone.config.yaml
    #             , " --seed "
    #             , pyclone.seed)
    # cat('\n')

    # if (run.pyclone) {
    #     cat(cmd)
    #     system(cmd)
    # }
    # cat('\n')

    # sample.results <- paste0(ClusterDir, "/", patientID, '.results.tsv')
    # cmd <- paste0(PyClone
    #            , " build_table --config_file "
    #            , pyclone.config.yaml
    #            , " --table_type old_style --out_file "
    #            , sample.results
    #            , " --max_clusters ", min(max(1, floor(length(SmallClusters[[clusterName]]$MutationsWithCluster) / 5)), 10)
    #            , " --burnin "
    #            , pyclone.burnin)
    # cat('\n')

    # if(run.pyclone) {
    #     cat(cmd)
    #     system(cmd)
    # }
    # cat('\n')
}

CreateOutputNoPyCloneRun <- function(clusterName, patientID, SmallClusters, patientDirToUse = new.dir) {
  
    ### give a name to the sample
    PyCloneRunName <- paste0(patientID, "_cluster", clusterName)
    ### create a specific subDirectory for this analysis
    ClusterDir <- paste0(patientDirToUse, "/", PyCloneRunName, "/")
    if (!dir.exists(ClusterDir)) {
        dir.create(ClusterDir, recursive = TRUE)
    }

    sample.results <- paste0(ClusterDir, "/", patientID, '.results.tsv')
    tmp <- matrix(0, nrow = length(SmallClusters[[clusterName]]$MutationsWithCluster), ncol = length(SmallClusters[[clusterName]]$RegionsInCluster))
    rownames(tmp) <- SmallClusters[[clusterName]]$MutationsWithCluster
    colnames(tmp) <- SmallClusters[[clusterName]]$RegionsInCluster
    tmp <- data.frame(tmp, cluster_id = clusterName)
    write.table(tmp, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, file = sample.results)
}

####################################################################################################
######################################### Helper functions #########################################
####################################################################################################

### get modified VAFs (extracted from findSimpleClustersWithSciClone)
getModifiedVAF <- function(rowNumber, regionClusterTable) {
    cat(rowNumber)
    varCountObs <- unlist(regionClusterTable[rowNumber, ]$var_counts)
    if (varCountObs == 0) {
        Depth <- as.numeric(regionClusterTable[rowNumber, ]$ref_counts)
    }
    if (varCountObs != 0) {
        Depth <- round(varCountObs / unlist(regionClusterTable[rowNumber, ]$phyloCCF / 2))
    }
    RefCount    <- Depth - varCountObs
    mutID       <- unlist(regionClusterTable[rowNumber, ]$mutation_id)
    outPutRow   <- cbind(unlist(strsplit(mutID, split = ":"))[2], as.numeric(unlist(strsplit(mutID, split = ":"))[3]), as.numeric(RefCount), as.numeric(varCountObs), as.numeric(varCountObs / Depth))
    if (as.numeric(varCountObs / Depth) > 0.75) {
        outPutRow <- cbind(unlist(strsplit(mutID, split = ":"))[2], as.numeric(unlist(strsplit(mutID, split = ":"))[3]), as.numeric(varCountObs), as.numeric(varCountObs), 0.5)
    }
    colnames(outPutRow) <- c('chr', 'pos', 'ref_count', 'var_count', 'vaf')
    return(outPutRow)
}



