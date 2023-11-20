create.subclonal.copy.number <- function(seg.mat.copy, min.subclonal = 0.1) {
    seg.out           <- seg.mat.copy
    seg.out$nMaj1     <- floor(as.numeric(seg.out$COPY_NUMBER_A))
    seg.out$nMaj2     <- ceiling(as.numeric(seg.out$COPY_NUMBER_A))
    seg.out$nMin1     <- floor(as.numeric(seg.out$COPY_NUMBER_B))
    seg.out$nMin2     <- ceiling(as.numeric(seg.out$COPY_NUMBER_B))
  
    seg.out$fracMaj1  <- as.numeric(seg.out$nMaj2) - as.numeric(seg.out$COPY_NUMBER_A)
    seg.out$fracMaj2  <- 1 - as.numeric(seg.out$fracMaj1)
    seg.out$fracMin1  <- as.numeric(seg.out$nMin2) - as.numeric(seg.out$COPY_NUMBER_B)  
    seg.out$fracMin2  <- 1 - as.numeric(seg.out$fracMin1)
  
  
    # next, let's deal with the minimum subclonal 
    seg.out$fracMaj1 <- ifelse(seg.out$fracMaj1 < min.subclonal,
                               0,
                               ifelse(seg.out$fracMaj1 > (1 - min.subclonal),
                                      1,
                                      seg.out$fracMaj1))
    seg.out$fracMaj2 <- ifelse(seg.out$fracMaj2 < min.subclonal,
                               0,
                               ifelse(seg.out$fracMaj2 > (1 - min.subclonal),
                                      1,
                                      seg.out$fracMaj2))
    seg.out$fracMin1 <- ifelse(seg.out$fracMin1 < min.subclonal,
                               0,
                               ifelse(seg.out$fracMin1 > (1 - min.subclonal),
                                      1,
                                      seg.out$fracMin1))
    seg.out$fracMin2 <- ifelse(seg.out$fracMin2 < min.subclonal,
                               0,
                               ifelse(seg.out$fracMin2 > (1 - min.subclonal),
                                      1,
                                      seg.out$fracMin2))
  
    # how much of the genome for each tumour region is subject to subclonal copy number
    for (region in unique(seg.mat.copy$SampleID)) {
        seg.region <- seg.out[seg.out$SampleID == region,, drop = FALSE]
        prop.aber  <- 1 - length(which((seg.region$fracMaj1 == 1 | seg.region$fracMaj1 == 0) & (seg.region$fracMin1 == 1 | seg.region$fracMin1 == 0))) / nrow(seg.region)
        print(prop.aber)
    }
    # which ones work?
    seg.sorted     <- seg.out[(seg.out$fracMaj1 == seg.out$fracMin1 | seg.out$fracMaj1 == seg.out$fracMin2) | (seg.out$fracMaj1 == 1) | (seg.out$fracMaj1 == 0) | (seg.out$fracMin1 == 1) | (seg.out$fracMin1 == 0),, drop = FALSE]
    seg.problem    <- seg.out[which(!rownames(seg.out) %in% rownames(seg.sorted)),]
  
    # let's divide the problem further
    seg.out.final  <- c()

    for (a in 1:nrow(seg.problem)) {
        seg.info         <- seg.problem[a,,drop=FALSE]
        seg.info$fracA   <- seg.info$fracMaj1 * seg.info$fracMin1  
        seg.info$fracB   <- seg.info$fracMaj1 - seg.info$fracA
        seg.info$fracC   <- seg.info$fracMaj2 * seg.info$fracMin2
        seg.info$fracD   <- seg.info$fracMaj2 - seg.info$fracC

        # let's see how this works. 
        seg.info$nMaj_A <- seg.info$nMaj1
        seg.info$nMin_A <- seg.info$nMin1

        seg.info$nMaj_B <- seg.info$nMaj1
        seg.info$nMin_B <- seg.info$nMin2

        seg.info$nMaj_C <- seg.info$nMaj2
        seg.info$nMin_C <- seg.info$nMin2

        seg.info$nMaj_D <- seg.info$nMaj2
        seg.info$nMin_D <- seg.info$nMin1

        seg.out.final <- rbind(seg.out.final, seg.info)
    }
  
    # let's make the sorted easier to read
    seg.sorted$fracA            <- NA
    seg.sorted$fracB            <- NA
    seg.sorted$fracC            <- NA
    seg.sorted$fracD            <- NA
    seg.sorted$nMaj_A           <- NA
    seg.sorted$nMaj_B           <- NA
    seg.sorted$nMaj_C           <- NA
    seg.sorted$nMaj_D           <- NA
    seg.sorted$nMin_A           <- NA
    seg.sorted$nMin_B           <- NA
    seg.sorted$nMin_C           <- NA
    seg.sorted$nMin_D           <- NA

    # put the columns we do need in the right order
    cols.maj <- apply(cbind(seg.sorted$fracMaj1,seg.sorted$fracMaj2),1,which.max)
    cols.min <- apply(cbind(seg.sorted$fracMin1,seg.sorted$fracMin2),1,which.max)
    for (i in 1:nrow(seg.sorted)) {
        poss.frac.major         <- c(seg.sorted$fracMaj1[i], seg.sorted$fracMaj2[i])
        poss.frac.minor         <- c(seg.sorted$fracMin1[i], seg.sorted$fracMin2[i])
        poss.cpn.major          <- c(seg.sorted$nMaj1[i], seg.sorted$nMaj2[i])
        poss.cpn.minor          <- c(seg.sorted$nMin1[i], seg.sorted$nMin2[i])
        cols.major         <- cols.maj[i]
        cols.minor         <- cols.min[i]
        
        fracA.major         <- poss.frac.major[cols.major]
        fracA.minor         <- poss.frac.minor[cols.minor]
        fracB.major         <- poss.frac.major[which(!1:2 %in% cols.major)]
        fracB.minor         <- poss.frac.minor[which(!1:2 %in% cols.minor)]
        
        seg.sorted$nMaj_A[i]              <- poss.cpn.major[cols.major]
        seg.sorted$nMin_A[i]              <- poss.cpn.minor[cols.minor]
        seg.sorted$nMaj_B[i]              <- poss.cpn.major[which(!1:2 %in% cols.major)]
        seg.sorted$nMin_B[i]              <- poss.cpn.minor[which(!1:2 %in% cols.minor)]
            
        seg.sorted$fracA[i]               <- fracA.major
        seg.sorted$fracB[i]               <- fracB.major
        if(fracA.major == 1) {
            seg.sorted$fracA[i]  <- fracA.minor
            seg.sorted$fracB[i]  <- fracB.minor
            seg.sorted$nMaj_B[i] <- poss.cpn.major[cols.major]
        }
        
        if(fracA.minor == 1) {
            seg.sorted$fracA[i]  <- fracA.major
            seg.sorted$fracB[i]  <- fracB.major
            seg.sorted$nMin_B[i] <- poss.cpn.minor[cols.minor]
        }  
    }
    seg.final <- rbind(seg.out.final,seg.sorted) 

    #let's order this correctly
    seg.final <- seg.final[order(seg.final$SampleID, seg.final$chr,seg.final$startpos),, drop = FALSE]
    # finally, let's choose the columns we want and the order we want
    colnames(seg.final)
    col.names <- c("SampleID", "chr", "startpos", "endpos", "n.het", 
                   "cnTotal", "nMajor", "nMinor", "Ploidy", "ACF",
                   "COPY_NUMBER_A", "COPY_NUMBER_B", "fracA", "nMaj_A", "nMin_A",
                   "fracB", "nMaj_B", "nMin_B", 
                   "fracC", "nMaj_C", "nMin_C", 
                   "fracD", "nMaj_D", "nMin_D")
    seg.final  <- seg.final[, col.names]
    return(seg.final)
}


determineIndelCorrectionFactor <- function(patient, mut.table, regions.to.use, seg.mat.phylo, seg.mat.copy) {
    # First, determine which mutations are considered ubiquitous
    gender <- "male"
    head(mut.table)

    #let's select the regions of interest
    UbiqMuts   <- rownames(mut.table)[which(rowSums(ifelse(mut.table[, paste0(regions.to.use, '.VAF'), drop = FALSE] > 5, 1, 0)) %in% length(regions.to.use))]
    IndelMuts  <- intersect(UbiqMuts, rownames(mut.table[mut.table$Use.For.Plots.Indel %in% TRUE,]))
    SNVMuts    <- UbiqMuts[!UbiqMuts %in% IndelMuts]
    # Next, let's do this check each region
    regionCorrectionFactors <- c()
    if (length(IndelMuts) == 0) {
        regionCorrectionFactors <- rep(1, length(regions.to.use))
    } else {

        for (region in regions.to.use) {
            region.mut.table <- mut.table
            region.mut.table <- region.mut.table[UbiqMuts,, drop = FALSE]

            region.seg.copy  <- seg.mat.copy[seg.mat.copy$SampleID %in% region,, drop = FALSE]
            region.seg.phylo <- seg.mat.phylo[seg.mat.phylo$SampleID %in% region,, drop = FALSE]
            pyclone.table    <- data.frame(t(sapply(1:nrow(region.mut.table),identify.subclonal.mut.copy.number.ascat,region.mut.table,region.seg.phylo,region,patient)), stringsAsFactors=FALSE)
            pyclone.table    <- pyclone.table[!is.na(pyclone.table$minor_cn),]
            pyclone.table    <- pyclone.table[!is.na(pyclone.table$ref_counts),]
            pyclone.table    <- pyclone.table[!duplicated(pyclone.table$mutation_id),]

            sample.purity   <- region.seg.copy$ACF[1]

            pyclone.table              <- pyclone.table[as.numeric(pyclone.table$ref_counts)+as.numeric(pyclone.table$var_counts)>=1,,drop=FALSE]
            region.earlyLate           <- earlyORlateGender(region=region,complete.mutation.table=pyclone.table,purity=sample.purity,gender=gender)
            region.phyloCCF            <- calculate_phylo_ccf(region=region,complete.mutation.table=pyclone.table,purity=sample.purity,order.by.pos = TRUE,gender=gender)
            region.phyloCCF            <- region.phyloCCF[!is.na(region.phyloCCF$phyloCCF),,drop=FALSE]

            region.phyloCCFIndel       <- region.phyloCCF[rownames(region.phyloCCF)%in%IndelMuts,,drop=FALSE]
            region.phyloCCFSNV         <- region.phyloCCF[rownames(region.phyloCCF)%in%SNVMuts,,drop=FALSE]

            CorrectionFactorCCF        <- median(region.phyloCCFSNV$phyloCCF)/median(region.phyloCCFIndel$phyloCCF)
            regionCorrectionFactors    <- c(regionCorrectionFactors,CorrectionFactorCCF)
        }
    }
    names(regionCorrectionFactors) <- regions.to.use
    if (length(IndelMuts) == 0) {
        regionCorrectionFactors[is.na(regionCorrectionFactors)] <- 1
    }
    return(regionCorrectionFactors)
}

identify.subclonal.mut.copy.number.ascat <- function(x, sub.mat.mut, sub.mat.copy, region, sample, sex = 'male') {
  
 
  mut                  <- sub.mat.mut[x,,drop=FALSE]
  ww                   <- which(as.numeric(sub.mat.copy$chr)==as.numeric(mut$chr)
                                &as.numeric(sub.mat.copy$startpos)<=as.numeric(mut$start)
                                &as.numeric(sub.mat.copy$endpos)>=as.numeric(mut$stop))
  copy                 <- sub.mat.copy[ww,,drop=FALSE]
  
  
  mutation_id         <- paste(sample,mut$chr,mut$start,mut$ref,sep=":")
  if (TRUE%in%grepl(paste(region,"d1ex1",sep=""),colnames(mut)))
  {
    ref_counts          <- mut[,gsub("-","\\.",paste(region, "d1ex1.ref.count",sep=""))]
    var_counts          <- mut[,gsub("-","\\.",paste(region, "d1ex1.var.count",sep=""))]    
  }
  if (!TRUE%in%grepl(paste(region,"d1ex1",sep=""),colnames(mut)))
  {
    ref_counts          <- mut[,gsub("-","\\.",paste(region, ".ref_count",sep=""))]
    var_counts          <- mut[,gsub("-","\\.",paste(region, ".var_count",sep=""))]    
  }
  
  normal_cn           <- 2
  region              <- region
  Reference_Base      <- mut$ref
  Alternate_Base      <- mut$var
  
  if(nrow(copy)!=1)
  {
    minor_cn          <- NA
    major_cn          <- NA
    major_raw         <- NA
    minor_raw         <- NA
    fracA             <- NA
    nMaj_A            <- NA
    nMin_A            <- NA
    fracB             <- NA
    nMaj_B            <- NA
    nMin_B            <- NA
    fracC             <- NA
    nMaj_C            <- NA
    nMin_C            <- NA
    fracD             <- NA
    nMaj_D            <- NA
    nMin_D            <- NA  
    
    output            <-  data.frame(mutation_id
                                     , ref_counts
                                     , var_counts
                                     , normal_cn
                                     , minor_cn
                                     , major_cn
                                     , major_raw
                                     , minor_raw
                                     , region
                                     , Reference_Base
                                     , Alternate_Base
                                     , fracA            
                                     , nMaj_A           
                                     , nMin_A           
                                     , fracB            
                                     , nMaj_B           
                                     , nMin_B           
                                     , fracC            
                                     , nMaj_C           
                                     , nMin_C           
                                     , fracD            
                                     , nMaj_D           
                                     , nMin_D           
                                     , stringsAsFactors=FALSE)
    return(output)
  }
  
  minor_cn            <- min(c(copy$nMajor,copy$nMinor))
  major_cn            <- max(c(copy$nMajor,copy$nMinor))
  major_raw           <- c(copy$COPY_NUMBER_A)
  minor_raw           <- c(copy$COPY_NUMBER_B)
  fracA               <- as.numeric(copy$fracA)
  nMaj_A              <- as.numeric(copy$nMaj_A)
  nMin_A              <- as.numeric(copy$nMin_A)
  fracB               <- as.numeric(copy$fracB)
  nMaj_B              <- as.numeric(copy$nMaj_B)
  nMin_B              <- as.numeric(copy$nMin_B)
  fracC               <- as.numeric(copy$fracC)
  nMaj_C              <- as.numeric(copy$nMaj_C)
  nMin_C              <- as.numeric(copy$nMin_C)
  fracD               <- as.numeric(copy$fracD)
  nMaj_D              <- as.numeric(copy$nMaj_D)
  nMin_D              <- as.numeric(copy$nMin_D)
  
  output            <-  data.frame(mutation_id
                                   , ref_counts
                                   , var_counts
                                   , normal_cn
                                   , minor_cn
                                   , major_cn
                                   , major_raw
                                   , minor_raw
                                   , region
                                   , Reference_Base
                                   , Alternate_Base
                                   , fracA            
                                   , nMaj_A           
                                   , nMin_A           
                                   , fracB            
                                   , nMaj_B           
                                   , nMin_B           
                                   , fracC            
                                   , nMaj_C           
                                   , nMin_C           
                                   , fracD            
                                   , nMaj_D           
                                   , nMin_D           
                                   , stringsAsFactors=FALSE)
  
  return(output)
  
}

earlyORlateGender <- function(region, complete.mutation.table, purity, gender) {
  
  
  # First of all load the needed packaged
  # 
  # suppressPackageStartupMessages(library(sequenza))
  # suppressPackageStartupMessages(library(bootstrap))
  # suppressPackageStartupMessages(library(boot))
  
  # Get the table ready, with only information for specific patient
  pyClone.tsv  <- complete.mutation.table[complete.mutation.table$region==region,]
  
  # and all the other stuff
  cellularity <- as.numeric(purity)
  major.cn <- unlist(pyClone.tsv$major_cn)
  abs.cn   <- unlist(pyClone.tsv$minor_cn) + unlist(pyClone.tsv$major_cn)
  depth.t  <- unlist(pyClone.tsv$ref_counts) + unlist(pyClone.tsv$var_counts)
  max.cn   <- max(abs.cn)
  VAF      <- unlist(pyClone.tsv$var_counts)/(unlist(pyClone.tsv$var_counts)+unlist(pyClone.tsv$ref_counts))

  if(gender=='male')
  {
    CNn      <- ifelse(do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='23'|do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='24',1,2)
    pyClone.tsv$normal_cn <- ifelse(do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='23'|do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='24',1,2)
    
  }
  if(gender=='female')
  {
    CNn      <- rep(2,nrow(pyClone.tsv))
    pyClone.tsv$normal_cn <- rep(2,nrow(pyClone.tsv))
  }
    
    
                     
                     
  # Use sequenza to estimate theoretical VAFs for each type of copy number
  types    <- types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 2)
  types.xy <- types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 1)
  types    <- rbind(types, types.xy)
  types    <- types[types$Mt >= 1, ]
  types$F  <- 0
  for (i in 1:nrow(types)) {
    types$F[i] <- theoretical.mufreq(cellularity = cellularity,
                                     CNn = types$CNn[i], CNt = types$CNt[i],
                                     Mt = types$Mt[i])
  }
  
  # Let's create some functions that can estimate whether early or late
  
  get.Mt <- function(F, depth.t, types, CNt, CNn, Mt){
    types <- types[types$CNn == CNn, ]
    l <- mufreq.dpois(mufreq = F, types$F[types$CNt== CNt&types$Mt<=Mt],
                                 depth.t = depth.t)
    l <- l/sum(l)
    L <- data.frame(l = l, Mt = types$Mt[types$CNt== CNt&types$Mt<=Mt])
  }
  
  get.conf <- function(F, depth.t){
    conf.int   <- cbind(prop.test(round(F*depth.t,0),depth.t)$conf[1]
                        ,prop.test(round(F*depth.t,0),depth.t)$conf[2])
    return(conf.int)
  }
  
  bootstrap.cf <- function(Vaf, cellularity, CNn, CNt, depth.t)
  {
    #print(i)
    if(Vaf==1)
    {
      conf.int   <- cbind(prop.test(round(Vaf*depth.t,0),depth.t)$conf[1]
                          ,prop.test(round(Vaf*depth.t,0),depth.t)$conf[2])
      
      lower      <- get.mut.mult(Vaf=conf.int[1],cellularity=cellularity,CNt=CNt,CNn=CNn)
      higher     <- get.mut.mult(Vaf=conf.int[2],cellularity=cellularity,CNt=CNt,CNn=CNn)
      conf.int   <- cbind(lower,higher)
      return(conf.int)
      
    }
    
    x          <- c(rep(1,round(Vaf*depth.t,0)),rep(0,(depth.t-round(Vaf*depth.t,0))))
    theta      <- function(x,i)
    {
      data      <- x[i]
      est       <- sum(data)/length(data)
      mut.multi <- (est *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity))
      return(mut.multi)
    }
    
    bt.res      <- boot::boot(x,theta,R=1000)
    bt.ci       <- boot::boot.ci(bt.res,type='norm')
    out         <- c(bt.ci$normal[2],bt.ci$normal[3])
    
    return(out)
    
  }
  
  
  get.mut.mult <- function(CNt,Vaf,cellularity,CNn)
  {
    
    return((Vaf *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity)))
    
  }
  
  
  get.cancer.cell.fraction <- function(Max.Likelihood,mut.mult)
  {
    predicted.Mtn   <- Max.Likelihood[,'Mt']
    ccf             <- mut.mult/predicted.Mtn
    return(ccf)
  }
  
  absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number, normal.copy.number)
  {
    f.function <- function (c,purity,local.copy.number)
    {
      
      return((purity*c) / (normal.copy.number*(1-purity) + purity*local.copy.number))
      
    }
    x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number))
    if(min(x)==0)
    {
      x[length(x)] <- 1
    }
    
    names(x)       <- seq(0.01,1,length.out=100)
    sub.cint <- function(x, prob = 0.95,n.alt,depth) {
      xnorm   <- x/sum(x)
      xsort   <- sort(xnorm, decreasing = TRUE)
      xcumLik <- cumsum(xsort)
      n = sum(xcumLik < prob) + 1
      LikThresh <- xsort[n]
      cint  <- x[xnorm >= LikThresh]
      all   <- as.numeric(names(x))
      cellu <- as.numeric(names(cint))
      l.t   <- cellu[1]
      r.t   <- cellu[length(cellu)]
      m     <- cellu[which.max(cint)]
      
      prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
      prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val
      
      data.frame(left = l.t, est = m, right = r.t,prob.subclonal=prob.subclonal,prob.clonal=prob.clonal)
    }
    
    
    return(sub.cint(x,n.alt=n.alt,depth=depth))
    
    
  }
  # add an absolute estimate of the cancer cell fraction
  
  get.all.mut.info <- function(i)
  {
    #print(i)
    # First estimate the VAF confidence intervals
    obs.VAF         <- VAF[i]
    mut.conf.0.05   <- get.conf(F=VAF[i],depth.t=depth.t[i])[1]
    mut.conf.0.95   <- get.conf(F=VAF[i],depth.t=depth.t[i])[2]
    
    if(abs.cn[i]==0)
    {
      output          <- cbind(obs.VAF
                               ,mut.conf.0.05
                               ,mut.conf.0.95
                               ,mut.multi=NA
                               ,mut.multi.0.05=NA
                               ,mut.multi.bstr.0.05=NA
                               ,mut.multi.0.95=NA
                               ,mut.multi.bstr.0.95=NA
                               ,Exp.Cpn =NA
                               ,Exp.Cpn.Likelihood=NA
                               ,ccf=NA
                               ,ccf.0.05=NA
                               ,ccf.btstr.0.05=NA
                               ,ccf.0.95=NA
                               ,ccf.btstr.0.95=NA
                               ,absolute.ccf=NA
                               ,absolute.ccf.0.05=NA
                               ,absoltue.ccf.0.95=NA
                               ,prob.subclonal=NA
                               ,prob.clonal=NA
                               ,timing='Not.Poss')
      return(output)
    }
    
    # Next estimate the likelihood relating to which copy number the mutation has
    L <- get.Mt(F = VAF[i],
                depth.t = depth.t[i], CNt = abs.cn[i],
                types = types, CNn = unlist(pyClone.tsv$normal_cn[i])
                ,Mt=major.cn[i])
    
    # Next determine the mut multiplicity
    mut.multi            <- get.mut.mult(CNt=abs.cn[i],Vaf=VAF[i],cellularity=cellularity,CNn=unlist(pyClone.tsv$normal_cn[i]))
    mut.multi.0.05       <- get.mut.mult(CNt=abs.cn[i],Vaf=mut.conf.0.05,cellularity=cellularity,CNn=unlist(pyClone.tsv$normal_cn[i]))
    mut.multi.0.95       <- get.mut.mult(CNt=abs.cn[i],Vaf=mut.conf.0.95,cellularity=cellularity,CNn=unlist(pyClone.tsv$normal_cn[i]))
    mut.multi.bstr       <- bootstrap.cf(Vaf=VAF[i],cellularity=cellularity,CNn=unlist(pyClone.tsv$normal_cn[i]),CNt=abs.cn[i],depth.t=depth.t[i])
    mut.multi.bstr.0.05  <- mut.multi.bstr[1]
    mut.multi.bstr.0.95  <- mut.multi.bstr[2]
    
    if(is.na(L$l)[1])
    {
      output          <- cbind(obs.VAF
                               ,mut.conf.0.05
                               ,mut.conf.0.95
                               ,mut.multi
                               ,mut.multi.0.05
                               ,mut.multi.bstr.0.05
                               ,mut.multi.0.95
                               ,mut.multi.bstr.0.95
                               ,Exp.Cpn =NA
                               ,Exp.Cpn.Likelihood=NA
                               ,ccf=NA
                               ,ccf.0.05=NA
                               ,ccf.btstr.0.05=NA
                               ,ccf.0.95=NA
                               ,ccf.btstr.0.95=NA
                               ,absolute.ccf=NA
                               ,absolute.ccf.0.05=NA
                               ,absoltue.ccf.0.95=NA
                               ,prob.subclonal=NA
                               ,prob.clonal=NA
                               ,timing='Not.Poss')
      return(output)
    }
    
    # Now determine which likelihood should be used
    Max.Likelihood   <- L[which.max(L$l),]
    absolute.calc        <- absolute.cancer.cell.fraction(n.alt=unlist(pyClone.tsv$var_counts)[i],depth=depth.t[i],purity=cellularity,local.copy.number=abs.cn[i],normal.copy.number = unlist(pyClone.tsv$normal_cn)[i])
    absolute.ccf.0.05    <- absolute.calc[1]
    absolute.ccf.0.95    <- absolute.calc[3]
    absolute.ccf         <- absolute.calc[2]
    prob.subclonal       <- absolute.calc[4]
    prob.clonal          <- absolute.calc[5]
    
    
    
    
    # Next determine the cancer cell fraction
    ccf             <- get.cancer.cell.fraction(Max.Likelihood,mut.multi)
    ccf.0.05        <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.0.05)
    ccf.btstr.0.05  <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.bstr.0.05)
    
    ccf.0.95        <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.0.95)
    ccf.btstr.0.95  <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.bstr.0.95)
    
    # Next determine the late cancer cell fraction
    # Make sure you also output the theoretical Copy (i.e. what it's closest to using maximum likelihood)
    expected.copy     <- Max.Likelihood[2]
    
    # Finally also make a suggestion about whether the mutation is early late or not possible to tell
    
    
    if(Max.Likelihood$Mt>1)
    {
      timing    <- 'early'
    }
    
    if(Max.Likelihood$Mt<=1)
    {
      timing    <- 'late'
    }
    
    if(major.cn[i]<=1)
    {
      timing    <- 'Not.Poss'
    }
    
    
    
    
    # Let's put this all together and output it
    output          <- data.frame(obs.VAF
                                  ,mut.conf.0.05
                                  ,mut.conf.0.95
                                  ,mut.multi
                                  ,mut.multi.0.05
                                  ,mut.multi.bstr.0.05
                                  ,mut.multi.0.95
                                  ,mut.multi.bstr.0.95
                                  ,Exp.Cpn =Max.Likelihood$Mt
                                  ,Exp.Cpn.Likelihood=Max.Likelihood$l
                                  ,ccf
                                  ,ccf.0.05
                                  ,ccf.btstr.0.05
                                  ,ccf.0.95
                                  ,ccf.btstr.0.95
                                  ,absolute.ccf.0.05    
                                  ,absolute.ccf.0.95    
                                  ,absolute.ccf
                                  ,prob.subclonal
                                  ,prob.clonal
                                  ,timing
                                  ,stringsAsFactors=FALSE)
    
    #output           <- data.frame(output,stringsAsFactors=FALSE)
    return(output)
    
  }
  
  output <- t(sapply(1:nrow(pyClone.tsv),get.all.mut.info))
  output <- data.frame(output,stringsAsFactors=FALSE)
  
  colnames(output) <- c('obs.VAF'
                        ,'mut.conf.0.05'
                        ,'mut.conf.0.95'
                        ,'mut.multi'
                        ,'mut.multi.0.05'
                        ,'mut.multi.bstr.0.05'
                        ,'mut.multi.0.95'
                        ,'mut.multi.bstr.0.95'
                        ,'Exp.Cpn'
                        ,'Exp.Cpn.Likelihood'
                        ,'ccf'
                        ,'ccf.0.05'
                        ,'ccf.btstr.0.05'
                        ,'ccf.0.95'
                        ,'ccf.btstr.0.95'
                        ,'absolute.ccf.0.05'    
                        ,'absolute.ccf.0.95'    
                        ,'absolute.ccf'
                        ,'prob.subclonal'
                        ,'prob.clonal'
                        ,'timing')
  
  
  out <- cbind(pyClone.tsv,output)
  return(out)
  
}

calculate_phylo_ccf <- function(region 
                                      ,complete.mutation.table
                                      ,purity
                                      ,region.seg.copy
                                      ,order.by.pos = TRUE
                                      ,min.subclonal = 0.05
                                      ,gender
                                      ,minVAFtoEXPLAIN = 0.05) 
{
  
  complete.mutation.table <- data.frame(complete.mutation.table,stringsAsFactors=FALSE)
  mut.order               <- complete.mutation.table$mutation_id
  if(order.by.pos)
  {
    chr <- do.call(rbind,strsplit(as.character(complete.mutation.table$mutation_id),split=":"))
    complete.mutation.table  <- complete.mutation.table[order(as.numeric(as.character(chr[,2])),as.numeric(as.character(chr[,3]))),]
    mut.order                <- complete.mutation.table$mutation_id
    
  }
  
  
  # suppressPackageStartupMessages(library(sequenza))
  # suppressPackageStartupMessages(library(bootstrap))
  # suppressPackageStartupMessages(library(boot))
  
  # Get the table ready, with only information for specific patient
  pyClone.tsv  <- complete.mutation.table[complete.mutation.table$region==region,,drop=FALSE]
  row.names    <- pyClone.tsv$mutation_id
  cellularity  <- as.numeric(purity)
  
  
  f.function <- function (c,purity,local.copy.number,normal.copy.number)
  {
    
    return(min(1,c((purity*c) / (normal.copy.number*(1-purity) + purity*local.copy.number))))
    
  }
  
  get.mut.mult <- function(CNt,Vaf,cellularity,CNn)
  {
    
    return((Vaf *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity)))
    
  }
  
  
  if(gender=='male')
  {
    CNn      <- ifelse(do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='23'|do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='24',1,2)
    pyClone.tsv$normal_cn <- ifelse(do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='23'|do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='24',1,2)
    
  }
  if(gender=='female')
  {
    CNn      <- rep(2,nrow(pyClone.tsv))
    pyClone.tsv$normal_cn <- rep(2,nrow(pyClone.tsv))
  }
  
  pyClone.tsv$expected.VAF       <- f.function(rep(1,nrow(pyClone.tsv)),cellularity,unlist(pyClone.tsv$major_raw)+unlist(pyClone.tsv$minor_raw),normal.copy.number = unlist(pyClone.tsv$normal_cn))
  
  # let's put it back together
  rownames(pyClone.tsv) <- pyClone.tsv$mutation_id
  pyClone.tsv           <- pyClone.tsv[unlist(mut.order),]
  pyClone.tsv$obs.VAF   <- as.numeric(pyClone.tsv$var_counts)/(as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts))
  
  get.mutCopyNum<- function(i)
  {
    
    # let's simplify things, and assume there are max two subclonal populations
    mutCopyNumber <- get.mut.mult(CNt=(as.numeric(pyClone.tsv$major_raw[i])+as.numeric(pyClone.tsv$minor_raw[i])),Vaf=pyClone.tsv$obs.VAF[i],cellularity=cellularity,CNn=unlist(pyClone.tsv$normal_cn[i]))
    return(mutCopyNumber)
    
    
  }
  get.absolute.ccf<- function(i)
  {
    #print(i)
    absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number,normal.copy.number)
    {
      
      f.function <- function (c,purity,local.copy.number,normal.copy.number)
      {
        
        return(min(1,c((purity*c) / (normal.copy.number*(1-purity) + purity*local.copy.number))))
        
      }
      x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number,normal.copy.number))
      if(min(x)==0)
      {
        x[length(x)] <- 1
      }
      
      names(x)       <- seq(0.01,1,length.out=100)
      sub.cint <- function(x, prob = 0.95,n.alt,depth) {
        xnorm   <- x/sum(x)
        xsort   <- sort(xnorm, decreasing = TRUE)
        xcumLik <- cumsum(xsort)
        n = sum(xcumLik < prob) + 1
        LikThresh <- xsort[n]
        cint  <- x[xnorm >= LikThresh]
        all   <- as.numeric(names(x))
        cellu <- as.numeric(names(cint))
        l.t   <- cellu[1]
        r.t   <- cellu[length(cellu)]
        m     <- cellu[which.max(cint)]
        
        prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
        prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val
        
        data.frame(left = l.t, est = m, right = r.t,prob.subclonal=prob.subclonal,prob.clonal=prob.clonal)
      }
      
      
      return(sub.cint(x,n.alt=n.alt,depth=depth))
      
      
    }
    
    absolute.calc        <- absolute.cancer.cell.fraction(n.alt=unlist(pyClone.tsv$var_counts)[i]
                                                          ,depth=(as.numeric(pyClone.tsv$ref_counts[i])+as.numeric(pyClone.tsv$var_counts[i]))
                                                          ,purity=cellularity
                                                          ,local.copy.number=(as.numeric(pyClone.tsv$major_raw[i])+as.numeric(pyClone.tsv$minor_raw[i]))
                                                          ,normal.copy.number=as.numeric(unlist(pyClone.tsv$normal_cn))[i]
    )
    absolute.ccf.0.05    <- absolute.calc[1]
    absolute.ccf.0.95    <- absolute.calc[3]
    absolute.ccf         <- absolute.calc[2]
    prob.subclonal       <- absolute.calc[4]
    prob.clonal          <- absolute.calc[5]
    return(cbind(absolute.ccf,absolute.ccf.0.05,absolute.ccf.0.95,absolute.ccf,prob.subclonal,prob.clonal))
    
    
    
    
  }
  
  pyClone.tsv$mutCopyNum        <- sapply(1:nrow(pyClone.tsv),get.mutCopyNum)
  absolute.ccfs                 <- t(sapply(1:nrow(pyClone.tsv),get.absolute.ccf))
  absolute.ccfs                 <- data.frame(absolute.ccfs,stringsAsFactors=FALSE)
  
  pyClone.tsv$absolute.ccf        <- as.numeric(absolute.ccfs$est)
  pyClone.tsv$absolute.ccf.0.05   <- as.numeric(absolute.ccfs$left)
  pyClone.tsv$absolute.ccf.0.95   <- as.numeric(absolute.ccfs$right)
  pyClone.tsv$phyloCCF            <- as.numeric(pyClone.tsv$mutCopyNum)
  pyClone.tsv$phyloCCF.0.05       <- sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,CNtumor,cellularity,CNn,i){
    VAF <- prop.test(var.count[i]
                     ,depth[i]
                     ,e[i]
    )$conf.int[1]
    return(get.mut.mult(CNt = CNtumor[i],Vaf = VAF,cellularity = cellularity,CNn = CNn[i]))
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e = as.numeric(pyClone.tsv$expected.VAF)
  ,CNtumor=as.numeric(pyClone.tsv$major_raw)+as.numeric(pyClone.tsv$minor_raw)
  ,cellularity=cellularity
  ,CNn=as.numeric(unlist(pyClone.tsv$normal_cn)))
  pyClone.tsv$phyloCCF.0.95      <- sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,CNtumor,cellularity,CNn,i){
    VAF <- prop.test(var.count[i]
                     ,depth[i]
                     ,e[i]
    )$conf.int[2]
    return(get.mut.mult(CNt = CNtumor[i],Vaf = VAF,cellularity = cellularity,CNn = CNn[i]))
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e = as.numeric(pyClone.tsv$expected.VAF)
  ,CNtumor=as.numeric(pyClone.tsv$major_raw)+as.numeric(pyClone.tsv$minor_raw)
  ,cellularity=cellularity
  ,CNn=as.numeric(unlist(pyClone.tsv$normal_cn)))
  
  
  pyClone.tsv$no.chrs.bearing.mut <- 1
  
  # let's check which subclonal mutations can be explained by copy number
  #pyClone.tsv.subclonal         <- pyClone.tsv[which(pyClone.tsv$absolute.ccf.0.95<1&pyClone.tsv$mutCopyNum>0.01),]
  subclonal.mutations            <- which(pyClone.tsv$absolute.ccf.0.95<1&pyClone.tsv$mutCopyNum>0.01&pyClone.tsv$obs.VAF>=minVAFtoEXPLAIN)
  # subclonal.mutations            <- which(absolute.ccfs$prob.subclonal>0.5&pyClone.tsv$mutCopyNum>0.01&pyClone.tsv$obs.VAF>=minVAFtoEXPLAIN)
  pyClone.tsv$whichFrac          <- NA
  pyClone.tsv$CPNChange <- 0
  
  
  if(length(subclonal.mutations)>0)
  {
    
    # assume subclonal muts are on one chromosome copy
    # therefore mutation copy number must be subclonal fraction 
    # of the higher CN subclone (i.e. lost in lower CN subclone) 
    # or 1 (i.e. present in both subclones)
    for (a in subclonal.mutations)
    {
      #print(a)
      mut.info <- pyClone.tsv[a,,drop=FALSE]
      # check whether the fracA is hundred%
      if (mut.info$fracA==1)
      {
        pyClone.tsv$whichFrac[a] <- c('A,B')
        next;
      }
      
      # let's see what happens when we only have A and B
      if (is.na(mut.info$fracC))
      {
        # determine which fraction has a copy number loss
        possible.subclonal.fractions <- c(1)
        if(unlist(mut.info$nMaj_A)>unlist(mut.info$nMaj_B))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracA))
        }
        if(unlist(mut.info$nMaj_B)>unlist(mut.info$nMaj_A))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracB))
        }
        if(unlist(mut.info$nMin_B)>unlist(mut.info$nMin_A))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracB))
        }
        if(unlist(mut.info$nMin_A)>unlist(mut.info$nMin_B))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracA))
        }
        
        
        best.CN                      <- possible.subclonal.fractions[which.min(abs(mut.info$mutCopyNum/possible.subclonal.fractions - 1))]
        if(best.CN==1|best.CN>=(1-min.subclonal))
        {
          pyClone.tsv$whichFrac[a] <- c('A,B')
          next;
        }
        var.count       = as.numeric(mut.info$var_counts)
        depth.count     =  as.numeric(mut.info$var_counts)+as.numeric(mut.info$ref_counts)
        expected.prop   = pyClone.tsv$expected.VAF[a] * best.CN
        
        # check whether subclonal CN results in clonal mutation
        # otherwise subclonal CN doesn't explain subclonal MCN
        if(best.CN != 1 & prop.test(var.count,depth.count*purity,expected.prop/purity,alternative = 'less')$p.value > 0.01)
        {
          if((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN)>1)
          {
            next;
          }
          
          if ((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN)<=1)
          {
            pyClone.tsv$phyloCCF[a]  = mut.info$mutCopyNum/best.CN
            pyClone.tsv$phyloCCF.0.05[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN 
            pyClone.tsv$phyloCCF.0.95[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[2],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN 
            pyClone.tsv$no.chrs.bearing.mut[a] = best.CN
            #pyClone.tsv$whichFrac[a] = c('A')
            pyClone.tsv$expected.VAF[a] = expected.prop
            pyClone.tsv$CPNChange[a] <- 1
          }
       
        }
      }
      if (!is.na(mut.info$fracC))
      {
        possible.subclonal.fractions <- c(1)
        # there are going to be quite a few options here, let's list them
        if((unlist(mut.info$nMaj_A)+unlist(mut.info$nMaj_B))>(unlist(mut.info$nMaj_C)+unlist(mut.info$nMaj_D)))
        {
          fracAB                       <- unlist(mut.info$fracA)+unlist(mut.info$fracB)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracAB)
        }
        if((unlist(mut.info$nMaj_A)+unlist(mut.info$nMaj_B))<(unlist(mut.info$nMaj_C)+unlist(mut.info$nMaj_D)))
        {
          fracCD                       <- unlist(mut.info$fracC)+unlist(mut.info$fracD)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracCD)
        }
        if((unlist(mut.info$nMin_A)+unlist(mut.info$nMin_D))>(unlist(mut.info$nMin_C)+unlist(mut.info$nMin_B)))
        {
          fracAD                       <- unlist(mut.info$fracA)+unlist(mut.info$fracD)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracAD)
        }
        if((unlist(mut.info$nMin_A)+unlist(mut.info$nMin_D))<(unlist(mut.info$nMin_C)+unlist(mut.info$nMin_B)))
        {
          fracBC                       <- unlist(mut.info$fracB)+unlist(mut.info$fracC)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracBC)
        }
        
        best.CN                      <-  possible.subclonal.fractions[which.min(abs(mut.info$mutCopyNum/possible.subclonal.fractions - 1))]
        if(best.CN==1|best.CN>=(1-min.subclonal))
        {
          pyClone.tsv$whichFrac[a] <- c('A,B,C,D')
          next;
        }
        var.count       = as.numeric(mut.info$var_counts)
        depth.count     =  as.numeric(mut.info$var_counts)+as.numeric(mut.info$ref_counts)
        expected.prop   = pyClone.tsv$expected.VAF[a] * best.CN
        
        # check whether subclonal CN results in clonal mutation
        # otherwise subclonal CN doesn't explain subclonal MCN
        if(best.CN != 1 & prop.test(var.count,depth.count*purity,expected.prop/purity,alternative = 'less')$p.value > 0.01)
        {
          
          if((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN)>1)
          {
            next;
          }
          
          if ((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN)<=1)
          {
            pyClone.tsv$phyloCCF[a]  = mut.info$mutCopyNum/best.CN
            pyClone.tsv$phyloCCF.0.05[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN 
            pyClone.tsv$phyloCCF.0.95[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[2],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN 
            pyClone.tsv$no.chrs.bearing.mut[a] = best.CN
            #pyClone.tsv$whichFrac[a] = c('A')
            pyClone.tsv$expected.VAF[a] = expected.prop
            pyClone.tsv$CPNChange[a] <- 1
          }
                    
        }
      }
      
      
      
    }
  }
  # Next, let's deal with potentially amplified mutations
  # convert MCN to subclonal fraction - tricky for amplified mutations
  # test for mutations in more than 1 copy
  
  
  
  p.vals = sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,i){
    prop.test(var.count[i]
              ,depth[i]
              ,e[i]
              ,alternative="greater")$p.value 
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e=as.numeric(pyClone.tsv$expected.VAF))
  
  amplified.muts <- which(p.vals<=0.05&pyClone.tsv$mutCopyNum>1)
  
  
  #copy numbers of subclones can only differ by 1 or 0 (as assumed when calling subclones)
  if(length(amplified.muts)>0)
  {    
    for(a in amplified.muts)
    {
      #print(a)
      mut.info    = pyClone.tsv[a,,drop=FALSE]
      if(is.na(mut.info$fracC))
      {
        max.CN2     = 0
        max.CN1     = unlist(mut.info$nMaj_A)
        frac1.mut   = unlist(mut.info$fracA)
        frac2.mut   = 0
        
        
        #swap subclones, so that the one with the higher CN is first
        if(unlist(mut.info$nMaj_B)>max.CN1){
          max.CN2 = max.CN1
          max.CN1 = unlist(mut.info$nMaj_B)
          frac2.mut = frac1.mut
          frac1.mut = unlist(mut.info$fracB)     
        }else{
          max.CN2   = unlist(mut.info$nMaj_B)
          frac2.mut = unlist(mut.info$fracB)
        }
        
        
        best.err = mut.info$mutCopyNum - 1
        best.CN=1
        allCNs = c(1)
        for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN = j * frac1.mut + k * frac2.mut
            err = abs(mut.info$mutCopyNum/potential.CN-1)
            #UpperConf = mut.info$
            if(err<best.err){
              #print(potential.CN)
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err=err
              best.CN = potential.CN
              allCNs  = c(allCNs,best.CN)
            }
          }
        }
        
        
        
        
        # let's just make sure we haven't created a subclonal mutation
        if(as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN)<1)
        {
          if(prop.test(as.numeric(pyClone.tsv$var_counts[a])/2
                       ,round(as.numeric(pyClone.tsv$var_counts[a])/as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN))
                       ,0.5
          )$p.value<0.05&best.CN>1)
          {
            best.CN <- max(allCNs[!allCNs%in%best.CN])
          }
        }
        
        
        pyClone.tsv$phyloCCF[a] = pyClone.tsv$mutCopyNum[a] / best.CN
        pyClone.tsv$expected.VAF[a]   = pyClone.tsv$expected.VAF[a]*best.CN
        pyClone.tsv$phyloCCF.0.05[a]  = pyClone.tsv$phyloCCF.0.05[a]/best.CN
        pyClone.tsv$phyloCCF.0.95[a]  = pyClone.tsv$phyloCCF.0.95[a]/best.CN
        
      }
      if(!is.na(mut.info$fracC))
      {
        
        max.CN2         = 0
        max.CN1         = unlist(mut.info$nMaj_A)
        frac1.mut       = unlist(mut.info$fracA)+unlist(mut.info$fracB)
        frac2.mut       = 0 
        
        
        #swap subclones, so that the one with the higher CN is first
        if(mut.info$nMaj_C>max.CN1){
          max.CN2       = max.CN1
          max.CN1       = unlist(mut.info$nMaj_C)
          frac2.mut     = frac1.mut
          frac1.mut     = unlist(mut.info$fracC)+unlist(mut.info$fracD)
          
          
        }else{
          max.CN2       = unlist(mut.info$nMaj_C)
          frac2.mut     = unlist(mut.info$fracC)+unlist(mut.info$fracD)
          
        }
        
        
        best.err = mut.info$mutCopyNum - 1
        best.CN=1
        allCNs <- c(best.CN)
        for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN  = j * frac1.mut + k * frac2.mut
            err           = abs(mut.info$mutCopyNu/potential.CN-1)
            
            
            if(err<best.err){
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err= err
              best.CN = potential.CN
              allCNs  = c(allCNs,best.CN)
            }
          }
        }
        
        # next, let's see whether we should also look at the minor allele
        max.CN2         = 0
        max.CN1         = unlist(mut.info$nMin_A)
        frac1.mut       = unlist(mut.info$fracA)+unlist(mut.info$fracD)
        frac2.mut       = 0 
        
        
        #swap subclones, so that the one with the higher CN is first
        if(unlist(mut.info$nMin_B)>max.CN1){
          max.CN2       = max.CN1
          max.CN1       = unlist(mut.info$nMin_B)
          frac2.mut     = frac1.mut
          frac1.mut     = unlist(mut.info$fracB)+unlist(mut.info$fracC)
          
          
        }else{
          max.CN2       = unlist(mut.info$nMin_B)
          frac2.mut     = unlist(mut.info$fracB)+unlist(mut.info$fracC)
          
        }
        
        
        
        for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN  = j * frac1.mut + k * frac2.mut
            err           = abs(mut.info$mutCopyNu/potential.CN-1)
            
            
            if(err<best.err){
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err= err
              best.CN = potential.CN
              allCNs   = c(allCNs,best.CN)
            }
          }
        }
        # let's just make sure we haven't created a subclonal mutation
        if(as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN)<1)
        {
          if(prop.test(as.numeric(pyClone.tsv$var_counts[a])/2
                       ,round(as.numeric(pyClone.tsv$var_counts[a])/as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN))
                       ,0.5
          )$p.value<0.05&best.CN>1)
          {
            best.CN <- max(allCNs[!allCNs%in%best.CN])
          }
        }
        
        pyClone.tsv$phyloCCF[a] = pyClone.tsv$mutCopyNum[a] / best.CN
        pyClone.tsv$expected.VAF[a]   = pyClone.tsv$expected.VAF[a]*best.CN
        pyClone.tsv$phyloCCF.0.05[a]  = pyClone.tsv$phyloCCF.0.05[a]/best.CN
        pyClone.tsv$phyloCCF.0.95[a]  = pyClone.tsv$phyloCCF.0.95[a]/best.CN 
      }
    }
  }
  
  # finally, let's sort out 'missing' ones
  pyClone.tsv[pyClone.tsv$var_counts==0,'no.chrs.bearing.mut'] <- 0
  pyClone.tsv[pyClone.tsv$var_counts==0,'expected.VAF']        <- 0
  pyClone.tsv[pyClone.tsv$var_counts==0,'absolute.ccf']        <- 0
  return(pyClone.tsv)
  
}

# new version with BH correction for subclonal mutations to be corrected to clonal
calculate_phylo_ccf_withBH <- function(region
                                       ,complete.mutation.table
                                       ,purity
                                       ,region.seg.copy
                                       ,order.by.pos = TRUE
                                       ,min.subclonal = 0.05
                                       ,gender
                                       ,minVAFtoEXPLAIN = 0.05) {
  
  complete.mutation.table <- data.frame(complete.mutation.table,stringsAsFactors=FALSE)
  mut.order               <- complete.mutation.table$mutation_id
  if(order.by.pos)
  {
    chr <- do.call(rbind,strsplit(as.character(complete.mutation.table$mutation_id),split=":"))
    complete.mutation.table  <- complete.mutation.table[order(as.numeric(as.character(chr[,2])),as.numeric(as.character(chr[,3]))),]
    mut.order                <- complete.mutation.table$mutation_id
    
  }
  
  
  # suppressPackageStartupMessages(library(sequenza))
  # suppressPackageStartupMessages(library(bootstrap))
  # suppressPackageStartupMessages(library(boot))
  
  # Get the table ready, with only information for specific patient
  pyClone.tsv  <- complete.mutation.table[complete.mutation.table$region==region,,drop=FALSE]
  row.names    <- pyClone.tsv$mutation_id
  cellularity  <- as.numeric(purity)
  
  
  f.function <- function (c,purity,local.copy.number,normal.copy.number)
  {
    
    return(min(1,c((purity*c) / (normal.copy.number*(1-purity) + purity*local.copy.number))))
    
  }
  
  get.mut.mult <- function(CNt,Vaf,cellularity,CNn)
  {
    
    return((Vaf *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity)))
    
  }
  
  
  if(gender=='male')
  {
    CNn      <- ifelse(do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='23'|do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='24',1,2)
    pyClone.tsv$normal_cn <- ifelse(do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='23'|do.call(rbind,strsplit(unlist(pyClone.tsv$mutation_id),split=":"))[,2]=='24',1,2)
    
  }
  if(gender=='female')
  {
    CNn      <- rep(2,nrow(pyClone.tsv))
    pyClone.tsv$normal_cn <- rep(2,nrow(pyClone.tsv))
  }
  
  pyClone.tsv$expected.VAF       <- f.function(rep(1,nrow(pyClone.tsv)),cellularity,unlist(pyClone.tsv$major_raw)+unlist(pyClone.tsv$minor_raw),normal.copy.number = unlist(pyClone.tsv$normal_cn))
  
  # let's put it back together
  rownames(pyClone.tsv) <- pyClone.tsv$mutation_id
  pyClone.tsv           <- pyClone.tsv[unlist(mut.order),]
  pyClone.tsv$obs.VAF   <- as.numeric(pyClone.tsv$var_counts)/(as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts))
  
  get.mutCopyNum<- function(i)
  {
    
    # let's simplify things, and assume there are max two subclonal populations
    mutCopyNumber <- get.mut.mult(CNt=(as.numeric(pyClone.tsv$major_raw[i])+as.numeric(pyClone.tsv$minor_raw[i])),Vaf=pyClone.tsv$obs.VAF[i],cellularity=cellularity,CNn=unlist(pyClone.tsv$normal_cn[i]))
    return(mutCopyNumber)
    
    
  }
  get.absolute.ccf<- function(i)
  {
    #print(i)
    absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number,normal.copy.number)
    {
      
      f.function <- function (c,purity,local.copy.number,normal.copy.number)
      {
        
        return(min(1,c((purity*c) / (normal.copy.number*(1-purity) + purity*local.copy.number))))
        
      }
      x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number,normal.copy.number))
      if(min(x)==0)
      {
        x[length(x)] <- 1
      }
      
      names(x)       <- seq(0.01,1,length.out=100)
      sub.cint <- function(x, prob = 0.95,n.alt,depth) {
        xnorm   <- x/sum(x)
        xsort   <- sort(xnorm, decreasing = TRUE)
        xcumLik <- cumsum(xsort)
        n = sum(xcumLik < prob) + 1
        LikThresh <- xsort[n]
        cint  <- x[xnorm >= LikThresh]
        all   <- as.numeric(names(x))
        cellu <- as.numeric(names(cint))
        l.t   <- cellu[1]
        r.t   <- cellu[length(cellu)]
        m     <- cellu[which.max(cint)]
        
        prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
        prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val
        
        data.frame(left = l.t, est = m, right = r.t,prob.subclonal=prob.subclonal,prob.clonal=prob.clonal)
      }
      
      
      return(sub.cint(x,n.alt=n.alt,depth=depth))
      
      
    }
    
    absolute.calc        <- absolute.cancer.cell.fraction(n.alt=unlist(pyClone.tsv$var_counts)[i]
                                                          ,depth=(as.numeric(pyClone.tsv$ref_counts[i])+as.numeric(pyClone.tsv$var_counts[i]))
                                                          ,purity=cellularity
                                                          ,local.copy.number=(as.numeric(pyClone.tsv$major_raw[i])+as.numeric(pyClone.tsv$minor_raw[i]))
                                                          ,normal.copy.number=as.numeric(unlist(pyClone.tsv$normal_cn))[i]
    )
    absolute.ccf.0.05    <- absolute.calc[1]
    absolute.ccf.0.95    <- absolute.calc[3]
    absolute.ccf         <- absolute.calc[2]
    prob.subclonal       <- absolute.calc[4]
    prob.clonal          <- absolute.calc[5]
    return(cbind(absolute.ccf,absolute.ccf.0.05,absolute.ccf.0.95,absolute.ccf,prob.subclonal,prob.clonal))
    
    
    
    
  }
  
  pyClone.tsv$mutCopyNum        <- sapply(1:nrow(pyClone.tsv),get.mutCopyNum)
  absolute.ccfs                 <- t(sapply(1:nrow(pyClone.tsv),get.absolute.ccf))
  absolute.ccfs                 <- data.frame(absolute.ccfs,stringsAsFactors=FALSE)
  
  pyClone.tsv$absolute.ccf        <- as.numeric(absolute.ccfs$est)
  pyClone.tsv$absolute.ccf.0.05   <- as.numeric(absolute.ccfs$left)
  pyClone.tsv$absolute.ccf.0.95   <- as.numeric(absolute.ccfs$right)
  pyClone.tsv$phyloCCF            <- as.numeric(pyClone.tsv$mutCopyNum)
  pyClone.tsv$phyloCCF.0.05       <- sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,CNtumor,cellularity,CNn,i){
    VAF <- prop.test(var.count[i]
                     ,depth[i]
                     ,e[i]
    )$conf.int[1]
    return(get.mut.mult(CNt = CNtumor[i],Vaf = VAF,cellularity = cellularity,CNn = CNn[i]))
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e = as.numeric(pyClone.tsv$expected.VAF)
  ,CNtumor=as.numeric(pyClone.tsv$major_raw)+as.numeric(pyClone.tsv$minor_raw)
  ,cellularity=cellularity
  ,CNn=as.numeric(unlist(pyClone.tsv$normal_cn)))
  pyClone.tsv$phyloCCF.0.95      <- sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,CNtumor,cellularity,CNn,i){
    VAF <- prop.test(var.count[i]
                     ,depth[i]
                     ,e[i]
    )$conf.int[2]
    return(get.mut.mult(CNt = CNtumor[i],Vaf = VAF,cellularity = cellularity,CNn = CNn[i]))
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e = as.numeric(pyClone.tsv$expected.VAF)
  ,CNtumor=as.numeric(pyClone.tsv$major_raw)+as.numeric(pyClone.tsv$minor_raw)
  ,cellularity=cellularity
  ,CNn=as.numeric(unlist(pyClone.tsv$normal_cn)))
  
  
  pyClone.tsv$no.chrs.bearing.mut <- 1
  
  # let's check which subclonal mutations can be explained by copy number
  #pyClone.tsv.subclonal         <- pyClone.tsv[which(pyClone.tsv$absolute.ccf.0.95<1&pyClone.tsv$mutCopyNum>0.01),]
  subclonal.mutations            <- which(pyClone.tsv$absolute.ccf.0.95<1&pyClone.tsv$mutCopyNum>0.01)
  subclonal.mutations            <- which(absolute.ccfs$prob.subclonal>0.5&pyClone.tsv$mutCopyNum>0.01&pyClone.tsv$obs.VAF>=minVAFtoEXPLAIN)
  pyClone.tsv$whichFrac          <- NA
  pyClone.tsv$CPNChange <- 0
  
  
  p.vals = sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,i){
    prop.test(var.count[i]
              ,depth[i]
              ,e[i]
              ,alternative="less")$p.value 
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e=as.numeric(pyClone.tsv$expected.VAF))
  
  p.vals.adj <- p.adjust(p.vals,method = 'BH')
  subclonal.mutations <- which(pyClone.tsv$absolute.ccf.0.95<1&pyClone.tsv$mutCopyNum>0.01&p.vals.adj<0.01)
  
  if(length(subclonal.mutations)>0)
  {
    
    # assume subclonal muts are on one chromosome copy
    # therefore mutation copy number must be subclonal fraction 
    # of the higher CN subclone (i.e. lost in lower CN subclone) 
    # or 1 (i.e. present in both subclones)
    for (a in subclonal.mutations)
    {
      #print(a)
      mut.info <- pyClone.tsv[a,,drop=FALSE]
      # check whether the fracA is hundred%
      if (mut.info$fracA==1)
      {
        pyClone.tsv$whichFrac[a] <- c('A,B')
        next;
      }
      
      # let's see what happens when we only have A and B
      if (is.na(mut.info$fracC))
      {
        # determine which fraction has a copy number loss
        possible.subclonal.fractions <- c(1)
        if(unlist(mut.info$nMaj_A)>unlist(mut.info$nMaj_B))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracA))
        }
        if(unlist(mut.info$nMaj_B)>unlist(mut.info$nMaj_A))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracB))
        }
        if(unlist(mut.info$nMin_B)>unlist(mut.info$nMin_A))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracB))
        }
        if(unlist(mut.info$nMin_A)>unlist(mut.info$nMin_B))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracA))
        }
        
        
        best.CN                      <- possible.subclonal.fractions[which.min(abs(mut.info$mutCopyNum/possible.subclonal.fractions - 1))]
        if(best.CN==1|best.CN>=(1-min.subclonal))
        {
          pyClone.tsv$whichFrac[a] <- c('A,B')
          next;
        }
        var.count       = as.numeric(mut.info$var_counts)
        depth.count     =  as.numeric(mut.info$var_counts)+as.numeric(mut.info$ref_counts)
        expected.prop   = pyClone.tsv$expected.VAF[a] * best.CN
        
        # check whether subclonal CN results in clonal mutation
        # otherwise subclonal CN doesn't explain subclonal MCN
        if(best.CN != 1 & prop.test(var.count,depth.count*purity,expected.prop/purity,alternative = 'less')$p.value > 0.01)
        {
          if((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN)>1)
          {
            next;
          }
          
          if ((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN)<=1)
          {
            pyClone.tsv$phyloCCF[a]  = mut.info$mutCopyNum/best.CN
            pyClone.tsv$phyloCCF.0.05[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN 
            pyClone.tsv$phyloCCF.0.95[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[2],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN 
            pyClone.tsv$no.chrs.bearing.mut[a] = best.CN
            #pyClone.tsv$whichFrac[a] = c('A')
            pyClone.tsv$expected.VAF[a] = expected.prop
            pyClone.tsv$CPNChange[a] <- 1
          }
          
        }
      }
      if (!is.na(mut.info$fracC))
      {
        possible.subclonal.fractions <- c(1)
        # there are going to be quite a few options here, let's list them
        if((unlist(mut.info$nMaj_A)+unlist(mut.info$nMaj_B))>(unlist(mut.info$nMaj_C)+unlist(mut.info$nMaj_D)))
        {
          fracAB                       <- unlist(mut.info$fracA)+unlist(mut.info$fracB)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracAB)
        }
        if((unlist(mut.info$nMaj_A)+unlist(mut.info$nMaj_B))<(unlist(mut.info$nMaj_C)+unlist(mut.info$nMaj_D)))
        {
          fracCD                       <- unlist(mut.info$fracC)+unlist(mut.info$fracD)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracCD)
        }
        if((unlist(mut.info$nMin_A)+unlist(mut.info$nMin_D))>(unlist(mut.info$nMin_C)+unlist(mut.info$nMin_B)))
        {
          fracAD                       <- unlist(mut.info$fracA)+unlist(mut.info$fracD)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracAD)
        }
        if((unlist(mut.info$nMin_A)+unlist(mut.info$nMin_D))<(unlist(mut.info$nMin_C)+unlist(mut.info$nMin_B)))
        {
          fracBC                       <- unlist(mut.info$fracB)+unlist(mut.info$fracC)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracBC)
        }
        
        best.CN                      <-  possible.subclonal.fractions[which.min(abs(mut.info$mutCopyNum/possible.subclonal.fractions - 1))]
        if(best.CN==1|best.CN>=(1-min.subclonal))
        {
          pyClone.tsv$whichFrac[a] <- c('A,B,C,D')
          next;
        }
        var.count       = as.numeric(mut.info$var_counts)
        depth.count     =  as.numeric(mut.info$var_counts)+as.numeric(mut.info$ref_counts)
        expected.prop   = pyClone.tsv$expected.VAF[a] * best.CN
        
        # check whether subclonal CN results in clonal mutation
        # otherwise subclonal CN doesn't explain subclonal MCN
        if(best.CN != 1 & prop.test(var.count,depth.count*purity,expected.prop/purity,alternative = 'less')$p.value > 0.01)
        {
          
          if((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN)>1)
          {
            next;
          }
          
          if ((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN)<=1)
          {
            pyClone.tsv$phyloCCF[a]  = mut.info$mutCopyNum/best.CN
            pyClone.tsv$phyloCCF.0.05[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN 
            pyClone.tsv$phyloCCF.0.95[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[2],cellularity = cellularity,CNn = mut.info$normal_cn)/best.CN 
            pyClone.tsv$no.chrs.bearing.mut[a] = best.CN
            #pyClone.tsv$whichFrac[a] = c('A')
            pyClone.tsv$expected.VAF[a] = expected.prop
            pyClone.tsv$CPNChange[a] <- 1
          }
                    
        }
      }
      
      
      
    }
  }
  # Next, let's deal with potentially amplified mutations
  # convert MCN to subclonal fraction - tricky for amplified mutations
  # test for mutations in more than 1 copy
  
  
  
  p.vals = sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,i){
    prop.test(var.count[i]
              ,depth[i]
              ,e[i]
              ,alternative="greater")$p.value 
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e=as.numeric(pyClone.tsv$expected.VAF))
  
  amplified.muts <- which(p.vals<=0.05&pyClone.tsv$mutCopyNum>1)
  
  
  #copy numbers of subclones can only differ by 1 or 0 (as assumed when calling subclones)
  if(length(amplified.muts)>0)
  {    
    for(a in amplified.muts)
    {
      #print(a)
      mut.info    = pyClone.tsv[a,,drop=FALSE]
      if(is.na(mut.info$fracC))
      {
        max.CN2     = 0
        max.CN1     = unlist(mut.info$nMaj_A)
        frac1.mut   = unlist(mut.info$fracA)
        frac2.mut   = 0
        
        
        #swap subclones, so that the one with the higher CN is first
        if(unlist(mut.info$nMaj_B)>max.CN1){
          max.CN2 = max.CN1
          max.CN1 = unlist(mut.info$nMaj_B)
          frac2.mut = frac1.mut
          frac1.mut = unlist(mut.info$fracB)     
        }else{
          max.CN2   = unlist(mut.info$nMaj_B)
          frac2.mut = unlist(mut.info$fracB)
        }
        
        
        best.err = mut.info$mutCopyNum - 1
        best.CN=1
        allCNs = c(1)
        for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN = j * frac1.mut + k * frac2.mut
            err = abs(mut.info$mutCopyNum/potential.CN-1)
            #UpperConf = mut.info$
            if(err<best.err){
              #print(potential.CN)
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err=err
              best.CN = potential.CN
              allCNs  = c(allCNs,best.CN)
            }
          }
        }
        
        
        
        
        # let's just make sure we haven't created a subclonal mutation
        if(as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN)<1)
        {
          if(prop.test(as.numeric(pyClone.tsv$var_counts[a])/2
                       ,round(as.numeric(pyClone.tsv$var_counts[a])/as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN))
                       ,0.5
          )$p.value<0.05&best.CN>1)
          {
            best.CN <- max(allCNs[!allCNs%in%best.CN])
          }
        }
        
        
        pyClone.tsv$phyloCCF[a] = pyClone.tsv$mutCopyNum[a] / best.CN
        pyClone.tsv$expected.VAF[a]   = pyClone.tsv$expected.VAF[a]*best.CN
        pyClone.tsv$phyloCCF.0.05[a]  = pyClone.tsv$phyloCCF.0.05[a]/best.CN
        pyClone.tsv$phyloCCF.0.95[a]  = pyClone.tsv$phyloCCF.0.95[a]/best.CN
        
      }
      if(!is.na(mut.info$fracC))
      {
        
        max.CN2         = 0
        max.CN1         = unlist(mut.info$nMaj_A)
        frac1.mut       = unlist(mut.info$fracA)+unlist(mut.info$fracB)
        frac2.mut       = 0 
        
        
        #swap subclones, so that the one with the higher CN is first
        if(mut.info$nMaj_C>max.CN1){
          max.CN2       = max.CN1
          max.CN1       = unlist(mut.info$nMaj_C)
          frac2.mut     = frac1.mut
          frac1.mut     = unlist(mut.info$fracC)+unlist(mut.info$fracD)
          
          
        }else{
          max.CN2       = unlist(mut.info$nMaj_C)
          frac2.mut     = unlist(mut.info$fracC)+unlist(mut.info$fracD)
          
        }
        
        
        best.err = mut.info$mutCopyNum - 1
        best.CN=1
        allCNs <- c(best.CN)
        for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN  = j * frac1.mut + k * frac2.mut
            err           = abs(mut.info$mutCopyNu/potential.CN-1)
            
            
            if(err<best.err){
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err= err
              best.CN = potential.CN
              allCNs  = c(allCNs,best.CN)
            }
          }
        }
        # next, let's see whether we should also look at the minor allele
        max.CN2         = 0
        max.CN1         = unlist(mut.info$nMin_A)
        frac1.mut       = unlist(mut.info$fracA)+unlist(mut.info$fracD)
        frac2.mut       = 0 
        
        
        #swap subclones, so that the one with the higher CN is first
        if(unlist(mut.info$nMin_B)>max.CN1){
          max.CN2       = max.CN1
          max.CN1       = unlist(mut.info$nMin_B)
          frac2.mut     = frac1.mut
          frac1.mut     = unlist(mut.info$fracB)+unlist(mut.info$fracC)
          
          
        }else{
          max.CN2       = unlist(mut.info$nMin_B)
          frac2.mut     = unlist(mut.info$fracB)+unlist(mut.info$fracC)
          
        }
        
        
        #best.err = mut.info$mutCopyNum - 1
        #best.CN=1
        
        for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN  = j * frac1.mut + k * frac2.mut
            err           = abs(mut.info$mutCopyNu/potential.CN-1)
            
            
            if(err<best.err){
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err= err
              best.CN = potential.CN
              allCNs   = c(allCNs,best.CN)
            }
          }
        }
        # let's just make sure we haven't created a subclonal mutation
        if(as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN)<1)
        {
          if(prop.test(as.numeric(pyClone.tsv$var_counts[a])/2
                       ,round(as.numeric(pyClone.tsv$var_counts[a])/as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN))
                       ,0.5
          )$p.value<0.05&best.CN>1)
          {
            best.CN <- max(allCNs[!allCNs%in%best.CN])
          }
        }
        
        pyClone.tsv$phyloCCF[a] = pyClone.tsv$mutCopyNum[a] / best.CN
        pyClone.tsv$expected.VAF[a]   = pyClone.tsv$expected.VAF[a]*best.CN
        pyClone.tsv$phyloCCF.0.05[a]  = pyClone.tsv$phyloCCF.0.05[a]/best.CN
        pyClone.tsv$phyloCCF.0.95[a]  = pyClone.tsv$phyloCCF.0.95[a]/best.CN 
      }
    }
  }
  
  # finally, let's sort out 'missing' ones
  pyClone.tsv[pyClone.tsv$var_counts==0,'no.chrs.bearing.mut'] <- 0
  pyClone.tsv[pyClone.tsv$var_counts==0,'expected.VAF']        <- 0
  pyClone.tsv[pyClone.tsv$var_counts==0,'absolute.ccf']        <- 0
  return(pyClone.tsv)
  
}


fun.chrom.length <- function(seg) {
  # This function returns the length of chromosome 1,...,22
  # seg: first three columns of the minimum conistent region matrix
  
  chrom.length <- c()
  for (i in sort(as.numeric(unique(seg[,1])))) {
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.length <- c(chrom.length, sub.chrom[nrow(sub.chrom), 3])
  }
  return(chrom.length)
}


fun.add.chrom <- function(seg, chrom.length) {
  # This function adds the length of chromosome 1,...,(i - 1) to chromsomes 2,...,22
  # seg: first three columns of the minimum conistent region matrix
  # chrom.length: output of fun.chrom.length(...)
  if (1 %in% seg[, 1])
    seg.tmp <-  subset(seg, seg[, 1] == 1)
  else
    seg.tmp <- matrix(NA, nr = 0, nc = ncol(seg))  
  for (i in 2:max(max(as.numeric(seg[,1])),3)) {
    if(!(i %in% seg[, 1]))
      next()
    sub.chrom <- subset(seg, seg[, 1] == i)
    sub.chrom[, 2] <- as.numeric(sub.chrom[, 2]) + sum(as.numeric(chrom.length[1:(i - 1)]))
    sub.chrom[, 3] <- as.numeric(sub.chrom[, 3]) + sum(as.numeric(chrom.length[1:(i - 1)]))
    seg.tmp <- rbind(seg.tmp, sub.chrom)
  }
  return(seg.tmp)
} 


fun.chrom.mean <- function(seg) {
  # This function finds the middle of the positions for each chromosome.
  # seg: Matrix with three columns, containing the positions.
  chrom.mean <- c()
  for (i in sort(as.numeric(unique(seg[,1])))) {
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.mean <- c(chrom.mean, (min(sub.chrom[, 2]) + max(sub.chrom[, 3]))/2)
  }
  return(chrom.mean)
}



PasteVector <- function(v, sep = "") {
  
  vt <- v[1]
  if(length(v) > 1) {
    for(g in 2:length(v)) {
      vt <- paste(vt,v[g],sep = sep) 
    }
  }
  vt <- paste(vt," EnD",sep = "")
  out.v <- sub(" EnD","",vt)
  out.v <- sub("NA , ","",out.v)
  out.v <- sub(" , NA","",out.v)
  out.v <- sub(" , NA , "," , ",out.v)
  return(out.v)
}


copy.driven.clusterNEW      <-  function(cluster
                                      ,seg.mat.copy
                                      ,most.likely.cluster
                                      ,region.earlyLate.list
                                      ,min.prop.cens = 0
                                      ,loss.thresh = 0.25
                                      ,diff.thresh = 0.55) {
  
  cluster.muts                 <- names(most.likely.cluster[most.likely.cluster==cluster])
  cpn.pvals                    <- rep(1,length(cluster.muts))
  prop.region.loss             <- rep(0,length(cluster.muts))
  cpn.loss                     <- rep(2,length(cluster.muts))
  
  
  names(cpn.pvals)             <- cluster.muts
  names(prop.region.loss)      <- names(cpn.pvals)
  names(cpn.loss)              <- names(cpn.pvals)
  
  # create some NA matrices
  region.major.cpn             <- data.frame(matrix(NA,length(cluster.muts),ncol=length(region.earlyLate.list)))
  rownames(region.major.cpn)   <- cluster.muts
  colnames(region.major.cpn)   <- names(region.earlyLate.list)
  
  region.minor.cpn             <- region.major.cpn
  region.minor.raw             <- region.major.cpn
  region.major.raw             <- region.major.cpn
  region.ccf                   <- region.major.cpn
  region.no.chrs.bearingMuts   <- region.major.cpn
  
  for (i in 1:length(region.earlyLate.list))
  {
    region.earlyLate           <- region.earlyLate.list[[i]]
    rownames(region.earlyLate) <- region.earlyLate$mutation_id
    tmp                  <- intersect(cluster.muts,rownames(region.earlyLate))
    
    #
    region.major.cpn[tmp,i] <- round(unlist(region.earlyLate[tmp,'major_raw']))
    region.minor.cpn[tmp,i] <- round(unlist(region.earlyLate[tmp,'minor_raw']))
    
    # Let's make sure we haven't accidentally called a homozygous deletion
    tmp2 <- tmp[which(region.major.cpn[tmp,i]+region.minor.cpn[tmp,i]==0)]
    

    tmp2major <- tmp2[unlist(region.earlyLate[tmp2,'major_raw'])>=unlist(region.earlyLate[tmp2,'minor_raw'])]
    tmp2minor <- tmp2[unlist(region.earlyLate[tmp2,'major_raw'])<unlist(region.earlyLate[tmp2,'minor_raw'])]
    
    region.major.cpn[tmp2major,i] <- 1
    region.major.cpn[tmp2minor,i] <- 1
            
    region.major.raw[tmp,i] <- unlist(region.earlyLate[tmp,'major_raw'])
    region.minor.raw[tmp,i] <- unlist(region.earlyLate[tmp,'minor_raw'])
    
    
    region.no.chrs.bearingMuts[tmp,i] <- unlist(region.earlyLate[tmp,'no.chrs.bearing.mut'])
    
    region.ccf[tmp,i]       <- unlist(region.earlyLate[tmp,'absolute.ccf'])
  }
  
  region.ccf                 <- region.ccf[!is.na(rowSums(region.ccf)),,drop=FALSE]
  region.minor.cpn           <- region.minor.cpn[!is.na(rowSums(region.minor.cpn)),,drop=FALSE]
  region.major.cpn           <- region.major.cpn[!is.na(rowSums(region.major.cpn)),,drop=FALSE]
  region.no.chrs.bearingMuts <-  region.no.chrs.bearingMuts[!is.na(rowSums(region.no.chrs.bearingMuts)),,drop=FALSE]
  
  for (region in colnames(region.minor.cpn))
  {
    region.major.cpn[,region]
  }
  
  region.minor.cpn
  region.major.cpn
  
  region.ccf[region.ccf<=0.01] <- 0
  
  fun.check.con <- function(v)
  {
    return(length(unique(v))==1)
  }
  
  major.cons        <- rownames(region.major.cpn)[apply(region.major.cpn,1,fun.check.con)]
  minor.cons        <- rownames(region.minor.cpn)[apply(region.minor.cpn,1,fun.check.con)]
  cluster.muts.cens <- cluster.muts[cluster.muts%in%major.cons&cluster.muts%in%minor.cons]
  
  # is it worth looking at these further
  prop.not.cens     <- 1-(length(cluster.muts.cens)/length(cluster.muts))
  if(prop.not.cens>min.prop.cens)
  {
    # clearly there is some evidence for copy number, but how much?
    # can this be driven by minor copy number? 
    if(length(which(!rownames(region.major.cpn)%in%major.cons))>=1)
    {
      
      major.cpn.aber <- region.major.cpn[!rownames(region.major.cpn)%in%major.cons,,drop=FALSE]
      major.ccf.aber <- region.ccf[!rownames(region.major.cpn)%in%major.cons,,drop=FALSE]
      
      chrs           <- do.call(rbind,strsplit(rownames(major.cpn.aber),split=":"))[,2]
      tmp            <- major.cpn.aber
      tmp            <- rbind(apply(tmp,2,as.numeric))
      rownames(tmp)  <- rownames(major.cpn.aber)
      states         <- apply(cbind(chrs,tmp),1,PasteVector,sep=":")
      
      
      for (state in unique(states))
      {
        #print(state)
        sub.ccf    <- region.ccf[names(states)[which(states%in%state)],]
        state.cpn  <- PasteVector(unlist(strsplit(state,split=":"))[-1],sep=":")
        medianCCF  <- apply(sub.ccf,2,median)
        
        
        
        if(min(apply(sub.ccf,2,median))>loss.thresh&max(medianCCF)-min(medianCCF)<diff.thresh)
        {
          next;
        }
        
        # turn it into a vector
        cpn.states  <- c()
        ccfs        <- c()
        
        for (i in 1:ncol(sub.ccf))
        {
          
          ccf              <- sub.ccf[,i]
          names(ccf)       <- paste(colnames(sub.ccf)[i],rownames(sub.ccf),sep=":")
          cpn.state        <- rep(unlist(strsplit(state.cpn,split=":"))[i],nrow(sub.ccf))
          names(cpn.state) <- names(ccf)  
          
          ccfs          <- c(ccfs,ccf)
          cpn.states    <- c(cpn.states,cpn.state)
        }
        
        
        if(length(unique(ccfs))==1)
        {
          # we only have mutations with CCF<0.75, let's not bother
          next;
        }
        
        
        # do we only have two copy number states?
        if(length(unique(cpn.states))==2)
        {
          
          
          # Do all without loss exhibit a high ccf fraction? and vice versa
          region.state         <- as.numeric(unlist(strsplit(state.cpn,split=":")))
          names(region.state)  <- colnames(sub.ccf)
          
          # One way to check this is see if we perform a regression 
          # analysis what comes out as significant?
          region.names <- as.factor(do.call(rbind,strsplit(names(ccfs),split=":"))[,1])
          tmp <- summary(lm(ccfs~cpn.states+region.names))
          tmp <- tmp$coefficients
          tmp <- tmp[,4]
          min.pval <- min(tmp[grep('cpn',names(tmp))])
          
          # is copy number states significant
          if(min.pval>=0.99&!is.na(min.pval))
          {
            next;
          }
          
          
          # loss state
          # are any of these greater than the threshold
          # is there more than one samples that shows the loss or 
          if(length(names(which(region.state==min(region.state))))>1)
          {
            if(TRUE%in%c(apply(sub.ccf[,names(which(region.state==min(region.state))),drop=FALSE],2,median)>loss.thresh))
            {
              next;
            }
          }
          
          # gain state
          # are any of these less than the threshold
          if(length(names(which(region.state==max(region.state))))>1)
          {
            if(TRUE%in%c(apply(sub.ccf[,names(which(region.state==max(region.state)))],2,median)<=loss.thresh))
            {
              next;
            }
          }
          
          
          
          
          # let's do a wilcox test. 
          p.val <- wilcox.test(as.numeric(ccfs[order(cpn.states,decreasing=TRUE)])~as.numeric(cpn.states[order(cpn.states,decreasing=TRUE)]),alternative = 'less')$p.val
          for (mut in rownames(sub.ccf))
          {
           
            
            if(median(as.numeric(sub.ccf[mut,names(which(region.state==min(region.state)))]))>loss.thresh&(median(as.numeric(sub.ccf[mut,names(which(region.state==min(region.state)))]))-max(apply(sub.ccf,2,median)))<diff.thresh)
            {
              next;
            }
            
            
            cpn.pvals[mut] <- min(cpn.pvals[mut],p.val)
            if(p.val<0.05)
            {
              
              prop.region.loss[mut] <- length(which(apply(sub.ccf,2,median)<=loss.thresh|(max(apply(sub.ccf,2,median))-apply(sub.ccf,2,median))>diff.thresh))/ncol(sub.ccf)
              cpn.loss[mut] <- min(region.state) 
              
            }
            
          }
          next;
        }
        
        if (length(unique(cpn.states))>2)
        {
          #print(state)
          # can we combine the ccfs easily? 
          # first of all, combine ccfs with zero
          region.names <- as.factor(do.call(rbind,strsplit(names(ccfs),split=":"))[,1])
          tmp <- summary(lm(ccfs~cpn.states+region.names))
          tmp <- tmp$coefficients
          tmp <- tmp[,4]
          
          min.pval <- min(tmp[grep('cpn',names(tmp))])
          # is copy number states significant
          if(min.pval>=0.99&!is.na(min.pval))
          {
            next;
          }
          
          # Do all without loss exhibit a high ccf fraction? and vice versa
          # One question
          # Do the copy number and ccf loss samples match
          region.state         <- as.numeric(unlist(strsplit(state.cpn,split=":")))
          names(region.state)  <- colnames(sub.ccf)
          
          ccf.loss.samples  <- region.state[which(apply(sub.ccf,2,median)<=loss.thresh|abs(apply(sub.ccf,2,median)-max(apply(sub.ccf,2,median)))>=diff.thresh)]
          cpn.loss.samples  <- names(which(region.state<=max(ccf.loss.samples)))
          if(!identical(names(ccf.loss.samples),cpn.loss.samples))
          {
            next;
          }
          
          p.val <- coin::pvalue(coin::independence_test(as.numeric(ccfs[order(cpn.states)])~as.numeric(cpn.states[order(cpn.states)]),alternative = 'greater'))
          for (mut in rownames(sub.ccf))
          {
            # is this mutation actually lost? 
            if(median(as.numeric(sub.ccf[mut,cpn.loss.samples]))>loss.thresh&(median(as.numeric(sub.ccf[mut,cpn.loss.samples]))-median(as.numeric(sub.ccf[mut,cpn.loss.samples]))<diff.thresh))
            {
              next;
            }
            
            cpn.pvals[mut] <- min(cpn.pvals[mut],p.val)
            if(p.val<0.05)
            {
              # how many regions show the loss
              prop.region.loss[mut] <- length(which(apply(sub.ccf,2,median)<=loss.thresh|(max(apply(sub.ccf,2,median))-apply(sub.ccf,2,median))>diff.thresh))/ncol(sub.ccf)
              cpn.loss[mut] <- max(region.state[cpn.loss.samples]) 
              
            }
          }
          next;                
        }
        
      }
    }
    
    if(length(which(!rownames(region.minor.cpn)%in%minor.cons))>=1)
    {
      minor.cpn.aber <- region.minor.cpn[!rownames(region.minor.cpn)%in%minor.cons,,drop=FALSE]
      minor.ccf.aber <- region.ccf[!rownames(region.minor.cpn)%in%minor.cons,,drop=FALSE]
      
      chrs           <- do.call(rbind,strsplit(rownames(minor.cpn.aber),split=":"))[,2]
      tmp            <- minor.cpn.aber
      tmp            <- rbind(apply(tmp,2,as.numeric))
      rownames(tmp)  <- rownames(minor.cpn.aber)
      states         <- apply(cbind(chrs,tmp),1,PasteVector,sep=":")
      
      for (state in unique(states))
      {
        #print(state)
        sub.ccf    <- region.ccf[names(states)[which(states%in%state)],]
        state.cpn  <- PasteVector(unlist(strsplit(state,split=":"))[-1],sep=":")
        medianCCF  <- apply(sub.ccf,2,median)
        
        
        
        if(min(apply(sub.ccf,2,median))>loss.thresh&max(medianCCF)-min(medianCCF)<diff.thresh)
        {
          next;
        }
        
     
        # turn it into a vector
        cpn.states  <- c()
        ccfs        <- c()
        
        for (i in 1:ncol(sub.ccf))
        {
          
          ccf              <- sub.ccf[,i]
          names(ccf)       <- paste(colnames(sub.ccf)[i],rownames(sub.ccf),sep=":")
          cpn.state        <- rep(unlist(strsplit(state.cpn,split=":"))[i],nrow(sub.ccf))
          names(cpn.state) <- names(ccf)  
          
          ccfs          <- c(ccfs,ccf)
          cpn.states    <- c(cpn.states,cpn.state)
        }
        
        
        
        if(length(which(ccfs<=loss.thresh))==length(ccfs))
        {
          
          next;
        }
        
        if(length(unique(ccfs))==1)
        {
          next;
        }
        
        # do we only have two copy number states?
        if(length(unique(cpn.states))==2)
        {
          #print(state)
          
          # Do all without loss exhibit a high ccf fraction? and vice versa
          region.state         <- as.numeric(unlist(strsplit(state.cpn,split=":")))
          names(region.state)  <- colnames(sub.ccf)
          
          # One way to check this is see if we perform a regression 
          # analysis what comes out as significant?
          region.names <- as.factor(do.call(rbind,strsplit(names(ccfs),split=":"))[,1])
          tmp <- summary(lm(ccfs~cpn.states+region.names))
          tmp <- tmp$coefficients
          tmp <- tmp[,4]
          min.pval <- min(tmp[grep('cpn',names(tmp))])
          
          # is copy number states significant
          if(min.pval>=0.99&!is.na(min.pval))
          {
            next;
          }
          
          
          # loss state
          # are any of these greater than the threshold
          # is there more than one samples that shows the loss or 
          if(length(names(which(region.state==min(region.state))))>1)
          {
            if(TRUE%in%c(apply(sub.ccf[,names(which(region.state==min(region.state)))],2,median)>loss.thresh))
            {
              # Can this be explained by a 'non-complete loss'
              # If this is the case, we shouldn't be ignoring this sample 
              # Which of the loss ones have a CCF above the loss threshold
              HardToExplain <- names(which(apply(sub.ccf[,names(which(region.state==min(region.state)))],2,median)>loss.thresh))
              for (region in HardToExplain)
              {
                regionMutTable <- region.earlyLate.list[[region]][region.earlyLate.list[[region]]$mutation_id%in%rownames(sub.ccf),,drop=FALSE]
                PropMRC           <- length(which(regionMutTable$phyloCCF.0.95>=1&regionMutTable$absolute.ccf.0.95<1))/nrow(regionMutTable)
                PropCopyCorrected <- length(which(regionMutTable$no.chrs.bearing.mut<1))/nrow(regionMutTable)
                
                if(PropMRC<0.25&PropCopyCorrected<0.25)
                {
                  next;
                  
                }
                
              }
              
            }
          }
          
          # gain state
          # are any of these less than the threshold
          if(length(names(which(region.state==max(region.state))))>1)
          {
            if(TRUE%in%c(apply(sub.ccf[,names(which(region.state==max(region.state)))],2,median)<=loss.thresh))
            {
              next;
            }
          }
          
          
          
          
          # let's do a wilcox test. 
          p.val <- wilcox.test(as.numeric(ccfs[order(cpn.states,decreasing=TRUE)])~as.numeric(cpn.states[order(cpn.states,decreasing=TRUE)]),alternative = 'less')$p.val
          for (mut in rownames(sub.ccf))
          {
            if(median(as.numeric(sub.ccf[mut,names(which(region.state==min(region.state)))]))>loss.thresh&(median(as.numeric(sub.ccf[mut,names(which(region.state!=min(region.state)))]))-median(as.numeric(sub.ccf[mut,names(which(region.state==min(region.state)))])))<diff.thresh)
            {
              next;
            }
            
            
            cpn.pvals[mut] <- min(cpn.pvals[mut],p.val)
            if(p.val<0.05)
            {
              
              prop.region.loss[mut] <- length(which(apply(sub.ccf,2,median)<=loss.thresh))/ncol(sub.ccf)
              cpn.loss[mut] <- min(region.state) 
              
            }
            
          }
          next;
        }
        
        if (length(unique(cpn.states))>2)
        {
          #print(state)
          # can we combine the ccfs easily? 
          # first of all, combine ccfs with zero
          region.names <- as.factor(do.call(rbind,strsplit(names(ccfs),split=":"))[,1])
          tmp <- summary(lm(ccfs~cpn.states+region.names))
          tmp <- tmp$coefficients
          tmp <- tmp[,4]
          
          min.pval <- min(tmp[grep('cpn',names(tmp))])
          # is copy number states significant
          if(min.pval>=0.99&!is.na(min.pval))
          {
            next;
          }
          
          # Do all without loss exhibit a high ccf fraction? and vice versa
          # One question
          # Do the copy number and ccf loss samples match
          region.state         <- as.numeric(unlist(strsplit(state.cpn,split=":")))
          names(region.state)  <- colnames(sub.ccf)
          
         
          
          ccf.loss.samples  <- region.state[which(apply(sub.ccf,2,median)<=loss.thresh|abs(apply(sub.ccf,2,median)-max(apply(sub.ccf,2,median)))>=diff.thresh)]
          cpn.loss.samples  <- names(which(region.state<=max(ccf.loss.samples)))
          if(!identical(names(ccf.loss.samples),cpn.loss.samples))
          {
            next;
          }
          
          p.val <- coin::pvalue(coin::independence_test(as.numeric(ccfs[order(cpn.states)])~as.numeric(cpn.states[order(cpn.states)]),alternative = 'greater'))
          for (mut in rownames(sub.ccf))
          {
            # is this mutation actually lost? 
            if(median(as.numeric(sub.ccf[mut,cpn.loss.samples]))>loss.thresh)
            {
              next;
            }
            
            cpn.pvals[mut] <- min(cpn.pvals[mut],p.val)
            if(p.val<0.05)
            {
              # how many regions show the loss
              prop.region.loss[mut] <- length(which(apply(sub.ccf,2,median)<=loss.thresh))/ncol(sub.ccf)
              cpn.loss[mut] <- max(region.state[cpn.loss.samples]) 
            }
          }
          next;
          
          
        }
        
       }
    }
    
  }
    return(cbind(cpn.pvals,prop.region.loss,cpn.loss))
  
}



getMinCPN <- function(mutation_id, phylo.region.list) {
  MinCPN      <- 10
  RegionLoss  <- c()
  for (region in names(phylo.region.list))
  {
    RegionMut <- phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%mutation_id,]
    if(nrow(RegionMut)==0)
    {
      next;
    }
    if(unlist(RegionMut$minor_cn)==MinCPN)
    {
      RegionLoss <- c(RegionLoss,region)
    }
    if(unlist(RegionMut$minor_cn)<MinCPN)
    {
      RegionLoss <- region
      MinCPN     <- unlist(RegionMut$minor_cn)
    }
    
  }
  RegionLoss  <- PasteVector(RegionLoss,sep=":")
  out  <- cbind(MinCPN,RegionLoss)
  rownames(out) <- mutation_id
  return(out)
  
}


correct.clusters.from.table <- function(full_table_out 
                                        ,max_cpn_proportion = 0.75
                                        ,min_muts_present_in_region = 0.25
                                        ,p_val_diff = 0.001
                                        ,copy_corrected_pval = 0.01) {
  
  #full_table_input <- full_table_out
  full_table_out <- full_table_out[!is.na(full_table_out$CLUSTER),,drop=FALSE]
  full_table_out$mutation_id <- paste(full_table_out$PATIENT,full_table_out$CHR,full_table_out$POS,full_table_out$REF,full_table_out$ALT,sep=":")
  
  most_likely_clusters <- as.character(full_table_out[!duplicated(full_table_out$mutation_id),]$CLUSTER)
  names(most_likely_clusters) <- full_table_out[!duplicated(full_table_out$mutation_id),]$mutation_id
  all_tumour_regions   <- unique(as.character(full_table_out$SAMPLE))
  
  # first check for each cluster the proportion of mutations that are copy# corrected
  input_most_likely_cluster <- most_likely_clusters
  cluster_prop <- c()
  cluster_presence_absence <- c()
  for (cluster in unique(most_likely_clusters))
  {
    cluster_muts <- most_likely_clusters[most_likely_clusters%in%cluster]
    cluster_copy_num <- cluster_muts
    cluster_copy_num[names(cluster_muts)] <- 0
    spec_cluster_presence_absence <- rep(0,length(unique(as.character(full_table_out$SAMPLE))))
    names(spec_cluster_presence_absence)  <- unique(as.character(full_table_out$SAMPLE))
    
    for (region in all_tumour_regions)
    {
      region_table <- full_table_out[full_table_out$SAMPLE%in%region,,drop=FALSE]
      region_table <- region_table[region_table$mutation_id%in%names(cluster_muts),,drop=FALSE]
      if(nrow(region_table[region_table$CCF_PHYLO>0.01,])/length(cluster_muts)>min_muts_present_in_region)
      {
        spec_cluster_presence_absence[region] <- 1
      }
      
      if(length(which(region_table$no.chrs.bearing.mut<1&region_table$no.chrs.bearing.mut>0))>10)
      {
        # print(region)
      }
      cluster_copy_num[which(region_table$MUT_COPY<1&region_table$CCF_PHYLO>region_table$MUT_COPY)] <- 1
      
    }
    cluster_prop<- c(cluster_prop, sum(as.numeric(cluster_copy_num))/length(cluster_copy_num))
    cluster_presence_absence <- rbind(cluster_presence_absence,spec_cluster_presence_absence)
  }
  
  rownames(cluster_presence_absence) <- unique(most_likely_clusters)
  names(cluster_prop) <- unique(most_likely_clusters)
  
  authentic_clusters  <- names(cluster_prop)[cluster_prop<max_cpn_proportion]
  suspicious_clusters <- names(cluster_prop)[cluster_prop>=max_cpn_proportion]
  
  cluster_vec_presence_absence <- apply(cluster_presence_absence,1,PasteVector,sep=":")
  if(length(suspicious_clusters)==0)
  {
    return(full_table_out)
    
  }
  
  
  # attempt to merge into another cluster
  
  merge_table <- c()
  changed_muts <- c()
  for (cluster in suspicious_clusters)
  {
    cluster_muts                      <- most_likely_clusters[most_likely_clusters%in%cluster]
    spec_cluster_presence_absence     <- cluster_presence_absence[cluster,,drop=FALSE]
    spec_cluster_vec_presence_absence <- cluster_vec_presence_absence[cluster]
    other_clusters_to_consider        <- authentic_clusters[cluster_vec_presence_absence[authentic_clusters]==spec_cluster_vec_presence_absence]
    
    
    regions_to_look_at <- colnames(spec_cluster_presence_absence)[spec_cluster_presence_absence%in%1]
    potential_merges   <- c()
    
    for (region in regions_to_look_at )
    {
      
      region_table <- full_table_out[full_table_out$SAMPLE%in%region,,drop=FALSE]
      region_table_cluster <- region_table[region_table$mutation_id%in%names(cluster_muts),,drop=FALSE]
      
      
      for (other_cluster in other_clusters_to_consider)
      {
        other_cluster_muts                  <- most_likely_clusters[most_likely_clusters%in%other_cluster]
        region_table_other_cluster <- region_table[region_table$mutation_id%in%names(other_cluster_muts),,drop=FALSE]
        
        if(wilcox.test(region_table_other_cluster$CCF_PHYLO,region_table_cluster$MUT_COPY)$p.val>p_val_diff|wilcox.test(region_table_other_cluster$MUT_COPY,region_table_cluster$MUT_COPY)$p.val>p_val_diff)
        {
          potential_merges <- c(potential_merges,other_cluster)
          
        }
        
      }
      
    }
    
    potential_merges <- names(table(potential_merges))[table(potential_merges)==length(regions_to_look_at)]
    
    if(length(potential_merges)>=1)
    {
      cluster_to_merge_into <- as.character(sort(as.numeric(potential_merges),decreasing = F)[1])
      merge_table <- rbind(merge_table,c(cluster,cluster_to_merge_into))
    }
    
    
    
    
    
  }
  
  #finally do the merges
  if(length(merge_table)==0)
  {
    return(full_table_out)
  }
  
  for (row in 1:nrow(merge_table))
  {
    
    for (region in all_tumour_regions)
    {
      region_table <- full_table_out[full_table_out$SAMPLE%in%region,,drop=FALSE]
      region_ccf_cluster <- region_table[region_table$mutation_id%in%names(most_likely_clusters)[most_likely_clusters%in%merge_table[row,1]],]
      
      if(identical(region_ccf_cluster$CCF_PHYLO,region_ccf_cluster$MUT_COPY))
      {
        next;
      }
      
      muts_to_change <- region_ccf_cluster[region_ccf_cluster$CCF_PHYLO>region_ccf_cluster$MUT_COPY&region_ccf_cluster$MUT_COPY<1,]$mutation_id
      changed_muts <- c(changed_muts,muts_to_change)
      
      full_table_out[full_table_out$SAMPLE%in%region&full_table_out$mutation_id%in%muts_to_change,]$CCF_PHYLO <-full_table_out[full_table_out$SAMPLE%in%region&full_table_out$mutation_id%in%muts_to_change,]$MUT_COPY 
    }
    
    most_likely_clusters[most_likely_clusters%in%merge_table[row,1]] <- merge_table[row,2]
  }
  
  
  
  
  current_names <- names(sort(table(most_likely_clusters),decreasing = T))
  new_names     <- as.character(1:length(current_names))
  new_most_likely_cluster <- most_likely_clusters
  
  for (i in 1:length(most_likely_clusters))
  {
    old_cluster_name <- current_names[i]
    new_cluster_name <- new_names[i]
    new_most_likely_cluster[most_likely_clusters%in%old_cluster_name] <- new_cluster_name
  }
  
  for (cluster in unique(new_most_likely_cluster))
  {
    
    full_table_out[full_table_out$mutation_id%in%names(new_most_likely_cluster[new_most_likely_cluster%in%cluster]),]$CLUSTER <- cluster
    
  }
  
  
  
  most_likely_clusters <- new_most_likely_cluster
  cluster_prop <- c()
  cluster_presence_absence <- c()
  for (cluster in unique(most_likely_clusters))
  {
    cluster_muts <- most_likely_clusters[most_likely_clusters%in%cluster]
    cluster_copy_num <- cluster_muts
    cluster_copy_num[names(cluster_muts)] <- 0
    spec_cluster_presence_absence <- rep(0,length(unique(as.character(full_table_out$SAMPLE))))
    names(spec_cluster_presence_absence)  <- unique(as.character(full_table_out$SAMPLE))
    
    for (region in all_tumour_regions)
    {
      region_table <- full_table_out[full_table_out$SAMPLE%in%region,,drop=FALSE]
      region_table <- region_table[region_table$mutation_id%in%names(cluster_muts),,drop=FALSE]
      if(nrow(region_table[region_table$CCF_PHYLO>0.01,])/length(cluster_muts)>min_muts_present_in_region)
      {
        spec_cluster_presence_absence[region] <- 1
      }
      
      if(length(which(region_table$no.chrs.bearing.mut<1&region_table$no.chrs.bearing.mut>0))>10)
      {
        # print(region)
      }
      cluster_copy_num[which(region_table$MUT_COPY<1&region_table$CCF_PHYLO>region_table$MUT_COPY)] <- 1
      
    }
    cluster_prop<- c(cluster_prop, sum(as.numeric(cluster_copy_num))/length(cluster_copy_num))
    cluster_presence_absence <- rbind(cluster_presence_absence,spec_cluster_presence_absence)
  }
  
  rownames(cluster_presence_absence) <- unique(most_likely_clusters)
  names(cluster_prop) <- unique(most_likely_clusters)
  cluster_vec_presence_absence <- apply(cluster_presence_absence,1,PasteVector,sep=":")
  
  # return(full_table_out)
  # next look at each individual cluster and check whether copy number correction is necessary
  reversing_cluster <- c()
  for (cluster in unique(new_most_likely_cluster))
  {
    cluster_table_out <- full_table_out[full_table_out$CLUSTER%in%cluster,,drop=FALSE]
    cluster_table_out <- cluster_table_out[cluster_table_out$SAMPLE%in%colnames(cluster_presence_absence)[cluster_presence_absence[cluster,]==1],,drop=FALSE]
    
    copy_corrected_table_out <- cluster_table_out[cluster_table_out$MUT_COPY<cluster_table_out$CCF_PHYLO,,drop=FALSE]
    ccf_corrected_table_out  <- cluster_table_out[cluster_table_out$MUT_COPY>=cluster_table_out$CCF_PHYLO,,drop=FALSE]
    
    if(nrow(copy_corrected_table_out)==0)
    {
      next;
    }
    
    
    # reverse_copy_correction <- c()
    for (region in unique(copy_corrected_table_out$SAMPLE))
    {
      region_copy_corrected_table_out <- copy_corrected_table_out[copy_corrected_table_out$SAMPLE%in%region,,drop=FALSE]
      region_ccf_corrected_table_out  <- ccf_corrected_table_out[ccf_corrected_table_out$SAMPLE%in%region,,drop=FALSE]
      
      if(nrow(region_ccf_corrected_table_out)==0|nrow(region_copy_corrected_table_out==0))
      {
        next;
      }
      
      
      if(nrow(region_copy_corrected_table_out)<=2|nrow(region_ccf_corrected_table_out)<=2)
      {
        dist_copy_corrected    <- abs(median(region_copy_corrected_table_out$CCF_PHYLO) - median(region_ccf_corrected_table_out$CCF_PHYLO))
        dist_noncopy_corrected <- abs(median(region_copy_corrected_table_out$MUT_COPY) - median(region_ccf_corrected_table_out$CCF_PHYLO))
        
        if(dist_copy_corrected>dist_noncopy_corrected)
        {
          # reverse_copy_correction <- c(reverse_copy_correction,region)
          full_table_out[full_table_out$mutation_id%in%region_copy_corrected_table_out$mutation_id&full_table_out$SAMPLE%in%region,]$CCF_PHYLO <- full_table_out[full_table_out$mutation_id%in%region_copy_corrected_table_out$mutation_id&full_table_out$SAMPLE%in%region,]$MUT_COPY
          reversing_cluster <- c(reversing_cluster,cluster)
          next;
        }
        
      }
      
      region_p_val_copy_corrected    <- wilcox.test(region_copy_corrected_table_out$CCF_PHYLO,region_ccf_corrected_table_out$CCF_PHYLO,alternative = 'greater')$p.val
      region_p_val_noncopy_corrected <- wilcox.test(region_copy_corrected_table_out$MUT_COPY,region_ccf_corrected_table_out$CCF_PHYLO)$p.val
      
      dist_copy_corrected    <- abs(median(region_copy_corrected_table_out$CCF_PHYLO) - median(region_ccf_corrected_table_out$CCF_PHYLO))
      dist_noncopy_corrected <- abs(median(region_copy_corrected_table_out$MUT_COPY) - median(region_ccf_corrected_table_out$CCF_PHYLO))
      
      if(region_p_val_copy_corrected<region_p_val_noncopy_corrected&dist_copy_corrected>dist_noncopy_corrected)
      {
        # reverse_copy_correction <- c(reverse_copy_correction,region)
        full_table_out[full_table_out$mutation_id%in%region_copy_corrected_table_out$mutation_id&full_table_out$SAMPLE%in%region,]$CCF_PHYLO <- full_table_out[full_table_out$mutation_id%in%region_copy_corrected_table_out$mutation_id&full_table_out$SAMPLE%in%region,]$MUT_COPY
        reversing_cluster <- c(reversing_cluster,cluster)
      }
            
    }
    
  }
  
  return(full_table_out)
  
} 


plot.pyclone.clusters <- function(patient
                                  ,regionList
                                  ,mut.table
                                  ,regions.to.use
                                  ,mostLikelyClusters
                                  ,driverCat
                                  ,ccf = 'sanger'
                                  ,colour.choice = NA
                                  ,plot.median = TRUE) {
  # suppressPackageStartupMessages(library(plyr))
  if (ccf=='absolute')
  {
    
  #its hard to distinguish more than 8 different colours
  no.optima = length(unique(mostLikelyClusters))
  max.cols = 12
  # require(RColorBrewer)
  cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
  cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
  cols.opac      = paste(cols,'99',sep="")
  
  
  
  pch.vals       = c()
  k <- 21
  m <- 1  
  while(m<=no.optima)
  {
    pch.vals    = c(pch.vals,k)
    if (m==max.cols)
    {
      k <- k+1
    }
    m <- m+1
  }
  
  pch.vals     = pch.vals
  
  
  for(i in 1:(length(regionList)-1)){
    for(j in (i+1):length(regionList)){
      regionX = regionList[[i]]
      regionY = regionList[[j]]
      #regionX$region <- regionList[i]
      #regionY$region <- regionList[j]
      rownames(regionX) <- regionX$mutation_id
      rownames(regionY) <- regionY$mutation_id
      # let's make sure we're looking at the same thing
      tmp     <- intersect(rownames(regionX),rownames(regionY))
      tmp     <- intersect(names(mostLikelyClusters),tmp)
      regionX <- regionX[tmp,,drop=FALSE]
      regionY <- regionY[tmp,,drop=FALSE]
      mostLikelyClusters <- mostLikelyClusters[tmp]
      
      layout(cbind(rbind(c(1,1,1,1,1)),2))
      par(mar=c(5,5,7,1))
      plot(regionX$absolute.ccf,regionY$absolute.ccf
           ,type = "n"
           , main = patient
           , xlab = paste(names(regionList)[i]," CCF",sep="")
           , ylab = paste(names(regionList)[j]," CCF",sep="")
           , xlim = c(-0.05,1.25)
           , ylim = c(-0.05,1.25)
           ,xaxt='n'
           ,yaxt='n'
           ,cex.axis=1.3)
      
      abline(h=1,lty='dashed')
      abline(v=1,lty='dashed')
      axis(side = 1,at = seq(0,1,by=0.2))
      axis(side = 2,at = seq(0,1,by=0.2),las=2) 
      # c(0,max(plot.data[,i])*1.25))
      for(n in 1:no.optima){
        
        points(unlist(regionX$absolute.ccf)[mostLikelyClusters==n],unlist(regionY$absolute.ccf)[mostLikelyClusters==n]
               ,bg=cols.opac[n]
               ,col = cols[n]
               ,pch=pch.vals[n]
               ,cex=1.2)
        
      }
      for(n in 1:no.optima){
        
        points(unlist(regionX$absolute.ccf)[mostLikelyClusters==n],unlist(regionY$absolute.ccf)[mostLikelyClusters==n]
               ,bg=cols.opac[n]
               ,col = cols[n]
               ,pch=pch.vals[n]
               ,cex=1.2)
        
      }
      
      
      x.vals <- c()
      y.vals <- c()
      colors <- c()
      labs   <- c()
      for (n in 1:no.optima)
      {
        
        x.vals <- c(x.vals,median(unlist(regionX$absolute.ccf)[mostLikelyClusters==n]))
        y.vals <- c(y.vals,median(unlist(regionY$absolute.ccf)[mostLikelyClusters==n]))
        labs   <- c(labs,n)
        points(median(unlist(regionX$absolute.ccf)[mostLikelyClusters==n])
               ,median(unlist(regionY$absolute.ccf)[mostLikelyClusters==n])
               ,pch=21
               ,col=cols[n]
               ,bg=cols[n]
               ,cex=1.2*2.5)
        
      }
      
      # require(wordcloud)
      
      if(length(x.vals)>1)
      {
        wordcloud::textplot(x.vals,y.vals,words = labs,new=FALSE
                 ,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white')
                 ,show.lines=TRUE) 
      }

      if(length(x.vals)==1)
      {
        text(x.vals,y.vals,labels = labs)
      }
      
      
      
      
      
      # let's add the drivers
      driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
      driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
      driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
      driver.table$driverCategoryNumeric<-as.numeric(plyr::revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2, "2A"=3, "2"=4, "3"=5))) 
      if(nrow(driver.table)>1)
      {
         
        wordcloud::textplot(unlist(regionX[rownames(driver.table),'absolute.ccf'])
                 ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
                 ,words=c(driver.table$Gene.refGene)
                 ,new=FALSE
                 ,col=gray(driver.table$driverCategoryNumeric/7)
                 ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                 ,cex=0.7
                 ,show.lines=FALSE
                 ,pos=4)
      }
      if(nrow(driver.table)==1)
      {
        
        text(x = unlist(regionX[rownames(driver.table),'absolute.ccf']),
                 ,y=unlist(regionY[rownames(driver.table),'absolute.ccf'])
                 ,labels=c(driver.table$Gene.refGene)
                 #,new=FALSE
                 ,col=gray(driver.table$driverCategoryNumeric/7)
                 ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                 ,cex=0.7
                 #,show.lines=FALSE
                 ,pos=4)
      }
      
      
      
      
      #plot.new()
      par(mar=c(3,0.5,7,0))
      plot(1, type="n"
           , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
      for (k in 1:no.optima)
      {
        points(0
               ,k
               ,col=cols.opac[k]
               ,bg = cols[k]
               ,pch=pch.vals[k]
               ,cex=1.5)
        text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
             ,pos=4)
      }
      
      text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
      
    }
  }
  }
  
  if (ccf=='sanger')
  {
    #its hard to distinguish more than 8 different colours
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    # require(RColorBrewer)
    cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    
    
    for(i in 1:(length(regionList)-1)){
      for(j in (i+1):length(regionList)){
        regionX = regionList[[i]]
        regionY = regionList[[j]]
        regionX$region <- regions.to.use[i]
        regionY$region <- regions.to.use[j]
        rownames(regionX) <- regionX$mutation_id
        rownames(regionY) <- regionY$mutation_id
        # let's make sure we're looking at the same thing
        tmp     <- intersect(rownames(regionX),rownames(regionY))
        tmp     <- intersect(names(mostLikelyClusters),tmp)
        regionX <- regionX[tmp,,drop=FALSE]
        regionY <- regionY[tmp,,drop=FALSE]
        mostLikelyClusters <- mostLikelyClusters[tmp]
        
        layout(cbind(rbind(c(1,1,1,1,1)),2))
        par(mar=c(5,5,7,1))
        plot(regionX$ccf,regionY$ccf
             ,type = "n"
             , main = patient
             , xlab = paste(names(regionList)[i]," CCF",sep="")
             , ylab = paste(names(regionList)[j]," CCF",sep="")
             , xlim = c(-0.05,1.5)
             , ylim = c(-0.05,1.5)
             ,xaxt='n'
             ,yaxt='n'
             ,cex.axis=1.3)
        
        abline(h=1,lty='dashed')
        abline(v=1,lty='dashed')
        axis(side = 1,at = seq(0,1,by=0.2))
        axis(side = 2,at = seq(0,1,by=0.2),las=2) 
        # c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.optima){
          
          points(unlist(regionX$ccf)[mostLikelyClusters==n],unlist(regionY$ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        for(n in 1:no.optima){
          
          points(unlist(regionX$ccf)[mostLikelyClusters==n],unlist(regionY$ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        
        
        x.vals <- c()
        y.vals <- c()
        colors <- c()
        labs   <- c()
        for (n in 1:no.optima)
        {
          
          x.vals <- c(x.vals,median(unlist(regionX$ccf)[mostLikelyClusters==n]))
          y.vals <- c(y.vals,median(unlist(regionY$ccf)[mostLikelyClusters==n]))
          labs   <- c(labs,n)
          points(median(unlist(regionX$ccf)[mostLikelyClusters==n])
                 ,median(unlist(regionY$ccf)[mostLikelyClusters==n])
                 ,pch=21
                 ,col=cols[n]
                 ,bg=cols[n]
                 ,cex=1.2*2.5)
          
        }
        
        # require(wordcloud)
        wordcloud::textplot(x.vals,y.vals,words = labs,new=FALSE
                 ,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white')
                 ,show.lines=TRUE)
        
        # let's add the drivers
        driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
        driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
        driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
        driver.table$driverCategoryNumeric<-as.numeric(plyr::revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2, "2A"=3, "2"=4, "3"=5)))
        
        
        if(nrow(driver.table)>=1)
        {
          
          wordcloud::textplot(unlist(regionX[rownames(driver.table),'absolute.ccf'])
                   ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
                   ,words=c(driver.table$Gene.refGene)
                   ,new=FALSE
                   ,col=gray(driver.table$driverCategoryNumeric/7)
                   ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                   ,cex=0.7
                   ,show.lines=FALSE
                   ,pos=4)
        }
        
        
        
        
        #plot.new()
        par(mar=c(3,0.5,7,0))
        plot(1, type="n"
             , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
        for (k in 1:no.optima)
        {
          points(0
                 ,k
                 ,col=cols.opac[k]
                 ,bg = cols[k]
                 ,pch=pch.vals[k]
                 ,cex=1.5)
          text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
               ,pos=4)
        }
        
        text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
        
      }
    }
  }
  
  if (ccf=='phylo')
  {
    #its hard to distinguish more than 8 different colours
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    # require(RColorBrewer)
    cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    if(!is.na(colour.choice))
    {
      cols        = rep(colour.choice,length(cols)) 
      cols.opac   = rep(colour.choice,length(cols)) 
    }
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    
    
    for(i in 1:(length(regionList)-1)){
      for(j in (i+1):length(regionList)){
        regionX = regionList[[i]]
        regionY = regionList[[j]]
        regionX$region <- regions.to.use[i]
        regionY$region <- regions.to.use[j]
        rownames(regionX) <- regionX$mutation_id
        rownames(regionY) <- regionY$mutation_id
        # let's make sure we're looking at the same thing
        tmp     <- intersect(rownames(regionX),rownames(regionY))
        tmp     <- intersect(names(mostLikelyClusters),tmp)
        regionX <- regionX[tmp,,drop=FALSE]
        regionY <- regionY[tmp,,drop=FALSE]
        mostLikelyClusters <- mostLikelyClusters[tmp]
        
        layout(cbind(rbind(c(1,1,1,1,1)),2))
        par(mar=c(5,5,7,1))
        plot(regionX$phyloCCF,regionY$phyloCCF
             ,type = "n"
             , main = patient
             , xlab = paste(names(regionList)[i]," CCF",sep="")
             , ylab = paste(names(regionList)[j]," CCF",sep="")
             , xlim = c(-0.05,1.5)
             , ylim = c(-0.05,1.5)
             ,xaxt='n'
             ,yaxt='n'
             ,cex.axis=1.3)
        
        abline(h=1,lty='dashed')
        abline(v=1,lty='dashed')
        axis(side = 1,at = seq(0,1,by=0.2))
        axis(side = 2,at = seq(0,1,by=0.2),las=2) 
        # c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.optima){
          
          points(unlist(regionX$phyloCCF)[mostLikelyClusters==n],unlist(regionY$phyloCCF)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        for(n in 1:no.optima){
          
          points(unlist(regionX$phyloCCF)[mostLikelyClusters==n],unlist(regionY$phyloCCF)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        
        if(plot.median)
        {
          
        
        x.vals <- c()
        y.vals <- c()
        colors <- c()
        labs   <- c()
        for (n in 1:no.optima)
        {
          
          x.vals <- c(x.vals,median(unlist(regionX$phyloCCF)[mostLikelyClusters==n]))
          y.vals <- c(y.vals,median(unlist(regionY$phyloCCF)[mostLikelyClusters==n]))
          labs   <- c(labs,n)
          points(median(unlist(regionX$phyloCCF)[mostLikelyClusters==n])
                 ,median(unlist(regionY$phyloCCF)[mostLikelyClusters==n])
                 ,pch=21
                 ,col=cols[n]
                 ,bg=cols[n]
                 ,cex=1.2*2.5)
          
        }
        
        # require(wordcloud)
        if(length(x.vals)>1)
        {        wordcloud::textplot(x.vals,y.vals,words = labs,new=FALSE
                          ,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white')
                          ,show.lines=TRUE)
          
        }
        if(length(x.vals)>1)
        {        
          text(x.vals,y.vals,labels = labs,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white'))
                 
        }

        
        }
        
        
        
        # let's add the drivers
        driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
        driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
        driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
        driver.table$driverCategoryNumeric<-as.numeric(plyr::revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2,"2A"=3, "2"=4, "3"=5)))
        
  if(nrow(driver.table)>1)
        {
          
          wordcloud::textplot(unlist(regionX[rownames(driver.table),'phyloCCF'])
                   ,unlist(regionY[rownames(driver.table),'phyloCCF'])
                   ,words=c(driver.table$Gene.refGene)
                   ,new=FALSE
                   ,col=gray(driver.table$driverCategoryNumeric/7)
                   ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                   ,cex=0.7
                   ,show.lines=FALSE
                   ,pos=4)
        }
        
        if(nrow(driver.table)==1)
        {
          
          text(x = unlist(regionX[rownames(driver.table),'phyloCCF']),
               ,y=unlist(regionY[rownames(driver.table),'phyloCCF'])
               ,labels=c(driver.table$Gene.refGene)
               #,new=FALSE
               ,col=gray(driver.table$driverCategoryNumeric/7)
               ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
               ,cex=0.7
               #,show.lines=FALSE
               ,pos=4)
        }
        
        
        
        #plot.new()
        par(mar=c(3,0.5,7,0))
        plot(1, type="n"
             , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
        for (k in 1:no.optima)
        {
          points(0
                 ,k
                 ,col=cols.opac[k]
                 ,bg = cols[k]
                 ,pch=pch.vals[k]
                 ,cex=1.5)
          text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
               ,pos=4)
        }
        
        text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
        
      }
    }
  }

}

plot.pycloneMutCpn.clusters <- function(patient
                                  ,regionList
                                  ,mut.table
                                  ,regions.to.use
                                  ,mostLikelyClusters
                                  ,driverCat
                                  ,ccf = 'sanger'
                                  ,colour.choice = NA
                                  ,plot.median = TRUE) {
  # suppressPackageStartupMessages(library(plyr))
  if (ccf=='absolute')
  {
    
    #its hard to distinguish more than 8 different colours
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    # require(RColorBrewer)
    cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    
    
    for(i in 1:(length(regionList)-1)){
      for(j in (i+1):length(regionList)){
        regionX = regionList[[i]]
        regionY = regionList[[j]]
        rownames(regionX) <- regionX$mutation_id
        rownames(regionY) <- regionY$mutation_id
        # let's make sure we're looking at the same thing
        tmp     <- intersect(rownames(regionX),rownames(regionY))
        tmp     <- intersect(names(mostLikelyClusters),tmp)
        regionX <- regionX[tmp,,drop=FALSE]
        regionY <- regionY[tmp,,drop=FALSE]
        mostLikelyClusters <- mostLikelyClusters[tmp]
        
        layout(cbind(rbind(c(1,1,1,1,1)),2))
        par(mar=c(5,5,7,1))
        plot(regionX$absolute.ccf,regionY$absolute.ccf
             ,type = "n"
             , main = patient
             , xlab = paste(names(regionList)[i]," CCF",sep="")
             , ylab = paste(names(regionList)[j]," CCF",sep="")
             , xlim = c(-0.05,1.25)
             , ylim = c(-0.05,1.25)
             ,xaxt='n'
             ,yaxt='n'
             ,cex.axis=1.3)
        
        abline(h=1,lty='dashed')
        abline(v=1,lty='dashed')
        axis(side = 1,at = seq(0,1,by=0.2))
        axis(side = 2,at = seq(0,1,by=0.2),las=2) 
        # c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.optima){
          
          points(unlist(regionX$absolute.ccf)[mostLikelyClusters==n],unlist(regionY$absolute.ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        for(n in 1:no.optima){
          
          points(unlist(regionX$absolute.ccf)[mostLikelyClusters==n],unlist(regionY$absolute.ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        
        
        x.vals <- c()
        y.vals <- c()
        colors <- c()
        labs   <- c()
        for (n in 1:no.optima)
        {
          
          x.vals <- c(x.vals,median(unlist(regionX$absolute.ccf)[mostLikelyClusters==n]))
          y.vals <- c(y.vals,median(unlist(regionY$absolute.ccf)[mostLikelyClusters==n]))
          labs   <- c(labs,n)
          points(median(unlist(regionX$absolute.ccf)[mostLikelyClusters==n])
                 ,median(unlist(regionY$absolute.ccf)[mostLikelyClusters==n])
                 ,pch=21
                 ,col=cols[n]
                 ,bg=cols[n]
                 ,cex=1.2*2.5)
          
        }
        
        # require(wordcloud)
        
        if(length(x.vals)>1)
        {
          wordcloud::textplot(x.vals,y.vals,words = labs,new=FALSE
                   ,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white')
                   ,show.lines=TRUE) 
        }
        
        if(length(x.vals)==1)
        {
          text(x.vals,y.vals,labels = labs)
        }
        
        
        
        
        
        # let's add the drivers
        driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
        driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
        driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
        driver.table$driverCategoryNumeric<-as.numeric(plyr::revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2,"2A"=3,"2"=4, "3"=5))) 
        if(nrow(driver.table)>1)
        {
          
          wordcloud::textplot(unlist(regionX[rownames(driver.table),'absolute.ccf'])
                   ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
                   ,words=c(driver.table$Gene.refGene)
                   ,new=FALSE
                   ,col=gray(driver.table$driverCategoryNumeric/7)
                   ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                   ,cex=0.7
                   ,show.lines=FALSE
                   ,pos=4)
        }
        if(nrow(driver.table)==1)
        {
          
          text(x = unlist(regionX[rownames(driver.table),'absolute.ccf']),
               ,y=unlist(regionY[rownames(driver.table),'absolute.ccf'])
               ,labels=c(driver.table$Gene.refGene)
               #,new=FALSE
               ,col=gray(driver.table$driverCategoryNumeric/7)
               ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
               ,cex=0.7
               #,show.lines=FALSE
               ,pos=4)
        }
        
        
        
        
        #plot.new()
        par(mar=c(3,0.5,7,0))
        plot(1, type="n"
             , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
        for (k in 1:no.optima)
        {
          points(0
                 ,k
                 ,col=cols.opac[k]
                 ,bg = cols[k]
                 ,pch=pch.vals[k]
                 ,cex=1.5)
          text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
               ,pos=4)
        }
        
        text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
        
      }
    }
  }
  
  if (ccf=='sanger')
  {
    #its hard to distinguish more than 8 different colours
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    # require(RColorBrewer)
    cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    
    
    for(i in 1:(length(regionList)-1)){
      for(j in (i+1):length(regionList)){
        regionX = regionList[[i]]
        regionY = regionList[[j]]
        regionX$region <- regions.to.use[i]
        regionY$region <- regions.to.use[j]
        rownames(regionX) <- regionX$mutation_id
        rownames(regionY) <- regionY$mutation_id
        # let's make sure we're looking at the same thing
        tmp     <- intersect(rownames(regionX),rownames(regionY))
        tmp     <- intersect(names(mostLikelyClusters),tmp)
        regionX <- regionX[tmp,,drop=FALSE]
        regionY <- regionY[tmp,,drop=FALSE]
        mostLikelyClusters <- mostLikelyClusters[tmp]
        
        layout(cbind(rbind(c(1,1,1,1,1)),2))
        par(mar=c(5,5,7,1))
        plot(regionX$ccf,regionY$ccf
             ,type = "n"
             , main = patient
             , xlab = paste(names(regionList)[i]," CCF",sep="")
             , ylab = paste(names(regionList)[j]," CCF",sep="")
             , xlim = c(-0.05,1.5)
             , ylim = c(-0.05,1.5)
             ,xaxt='n'
             ,yaxt='n'
             ,cex.axis=1.3)
        
        abline(h=1,lty='dashed')
        abline(v=1,lty='dashed')
        axis(side = 1,at = seq(0,1,by=0.2))
        axis(side = 2,at = seq(0,1,by=0.2),las=2) 
        # c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.optima){
          
          points(unlist(regionX$ccf)[mostLikelyClusters==n],unlist(regionY$ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        for(n in 1:no.optima){
          
          points(unlist(regionX$ccf)[mostLikelyClusters==n],unlist(regionY$ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        
        
        x.vals <- c()
        y.vals <- c()
        colors <- c()
        labs   <- c()
        for (n in 1:no.optima)
        {
          
          x.vals <- c(x.vals,median(unlist(regionX$ccf)[mostLikelyClusters==n]))
          y.vals <- c(y.vals,median(unlist(regionY$ccf)[mostLikelyClusters==n]))
          labs   <- c(labs,n)
          points(median(unlist(regionX$ccf)[mostLikelyClusters==n])
                 ,median(unlist(regionY$ccf)[mostLikelyClusters==n])
                 ,pch=21
                 ,col=cols[n]
                 ,bg=cols[n]
                 ,cex=1.2*2.5)
        }
        
        # require(wordcloud)
        wordcloud::textplot(x.vals,y.vals,words = labs,new=FALSE
                 ,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white')
                 ,show.lines=TRUE)
        
        
        
        
        
        # let's add the drivers
        driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
        driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
        driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
        driver.table$driverCategoryNumeric<-as.numeric(plyr::revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2,"2A"=3,"2"=4, "3"=5)))
        
        
        if(nrow(driver.table)>=1)
        {
          
          wordcloud::textplot(unlist(regionX[rownames(driver.table),'absolute.ccf'])
                   ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
                   ,words=c(driver.table$Gene.refGene)
                   ,new=FALSE
                   ,col=gray(driver.table$driverCategoryNumeric/7)
                   ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                   ,cex=0.7
                   ,show.lines=FALSE
                   ,pos=4)
        }
        
        
        
        
        #plot.new()
        par(mar=c(3,0.5,7,0))
        plot(1, type="n"
             , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
        for (k in 1:no.optima)
        {
          points(0
                 ,k
                 ,col=cols.opac[k]
                 ,bg = cols[k]
                 ,pch=pch.vals[k]
                 ,cex=1.5)
          text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
               ,pos=4)
        }
        
        text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
        
      }
    }
  }
  
  if (ccf=='phylo')
  {
    #its hard to distinguish more than 8 different colours
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    # require(RColorBrewer)
    cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    if(!is.na(colour.choice))
    {
      cols        = rep(colour.choice,length(cols)) 
      cols.opac   = rep(colour.choice,length(cols)) 
    }
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    
    
    for(i in 1:(length(regionList)-1)){
      for(j in (i+1):length(regionList)){
        regionX = regionList[[i]]
        regionY = regionList[[j]]
        regionX$region <- regions.to.use[i]
        regionY$region <- regions.to.use[j]
        rownames(regionX) <- regionX$mutation_id
        rownames(regionY) <- regionY$mutation_id
        # let's make sure we're looking at the same thing
        tmp     <- intersect(rownames(regionX),rownames(regionY))
        tmp     <- intersect(names(mostLikelyClusters),tmp)
        regionX <- regionX[tmp,,drop=FALSE]
        regionY <- regionY[tmp,,drop=FALSE]
        mostLikelyClusters <- mostLikelyClusters[tmp]
        
        layout(cbind(rbind(c(1,1,1,1,1)),2))
        par(mar=c(5,5,7,1))
        plot(regionX$mutCopyNum,regionY$mutCopyNum
             ,type = "n"
             , main = patient
             , xlab = paste(names(regionList)[i]," CCF",sep="")
             , ylab = paste(names(regionList)[j]," CCF",sep="")
             , xlim = c(-0.05,5.5)
             , ylim = c(-0.05,5.5)
             ,xaxt='n'
             ,yaxt='n'
             ,cex.axis=1.3)
        
        abline(h=1:5,lty='dashed')
        abline(v=1:5,lty='dashed')
        axis(side = 1,at = seq(0,5,by=0.2))
        axis(side = 2,at = seq(0,5,by=0.2),las=2) 
        # c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.optima){
          
          points(unlist(regionX$mutCopyNum)[mostLikelyClusters==n],unlist(regionY$mutCopyNum)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        for(n in 1:no.optima){
          
          points(unlist(regionX$mutCopyNum)[mostLikelyClusters==n],unlist(regionY$mutCopyNum)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        
        
        
        # let's add the drivers
        driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
        driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
        driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
        driver.table$driverCategoryNumeric<-as.numeric(plyr::revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2,"2A"=3,"2"=4, "3"=5)))
        
        
        if(nrow(driver.table)>1)
        {
          
          wordcloud::textplot(unlist(regionX[rownames(driver.table),'mutCopyNum'])
                   ,unlist(regionY[rownames(driver.table),'mutCopyNum'])
                   ,words=c(driver.table$Gene.refGene)
                   ,new=FALSE
                   ,col=gray(driver.table$driverCategoryNumeric/7)
                   ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                   ,cex=0.7
                   ,show.lines=FALSE
                   ,pos=4)
        }
        
        if(nrow(driver.table)==1)
        {
          
          text(x = unlist(regionX[rownames(driver.table),'mutCopyNum']),
               ,y=unlist(regionY[rownames(driver.table),'mutCopyNum'])
               ,labels=c(driver.table$Gene.refGene)
               #,new=FALSE
               ,col=gray(driver.table$driverCategoryNumeric/7)
               ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
               ,cex=0.7
               #,show.lines=FALSE
               ,pos=4)
        }
        
        
        
        #plot.new()
        par(mar=c(3,0.5,7,0))
        plot(1, type="n"
             , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
        for (k in 1:no.optima)
        {
          points(0
                 ,k
                 ,col=cols.opac[k]
                 ,bg = cols[k]
                 ,pch=pch.vals[k]
                 ,cex=1.5)
          text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
               ,pos=4)
        }
        
        text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
        
      }
    }
  }
  
}


plot.region.mutCopyNum  <- function(phylo.region.list
                                    , seg.mat.copy
                                    , mostLikelyClusters = most.likely.cluster
                                    , plot.separate.clusters = FALSE) {
  if(!plot.separate.clusters)
  {
    regions.to.use <- names(phylo.region.list)
    lyout <- c()
    for (i in 1:length(regions.to.use))
    {
      lyout <- rbind(lyout,rbind(c(rep(i,9),length(regions.to.use)+1)))
    }
    layout(lyout)
    for (region in regions.to.use)
    {
      
      # Failsafe in case column SAMPLE not present
      if ((! 'SAMPLE' %in% colnames(seg.mat.copy)) & ('SampleID' %in% colnames(seg.mat.copy))){
        seg.mat.copy$SAMPLE <- seg.mat.copy$SampleID
      }

      # rename columns
      sub.mat.copy               <- seg.mat.copy[seg.mat.copy$SAMPLE%in%region,,drop=FALSE]
      colnames(sub.mat.copy)[2]  <- 'Chromosome'
      colnames(sub.mat.copy)[3]  <- 'StartPosition'
      colnames(sub.mat.copy)[4]  <- 'EndPosition'
      region.phyloCCF            <- phylo.region.list[[region]]
      
      #pdf(early.late.pdf)
      par(mar=c(0.5,4,0.5,0.2))
      if(nrow(lyout)==1)
      {
        par(mar=c(14,4,14,0.2))
      }
      
      par(lend=1)
      plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                              ,TCGA.earlyLate=region.phyloCCF
                              ,sub.clonal=1
                              ,most.likely.cluster=mostLikelyClusters
      )
      mtext(unlist(strsplit(region,split="_"))[length(unlist(strsplit(region,split="_")))],side=2,cex=0.6,line=2,las=2)
      #its hard to distinguish more than 8 different colours
    }
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    # require(RColorBrewer)
    cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    par(mar=c(2,0.5,2,0))
    if(nrow(lyout)==1)
    {
      par(mar=c(14,0.5,14,0))
    }
    plot(1, type="n"
         , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
    for (k in 1:no.optima)
    {
      points(0
             ,k
             ,col=cols.opac[k]
             ,bg = cols[k]
             ,pch=pch.vals[k]
             ,cex=1.5)
      text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
           ,pos=4)
    }
    
    text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
  }
  if(plot.separate.clusters)
  {
    regions.to.use <- names(phylo.region.list)
    lyout <- c()
    for (i in 1:length(regions.to.use))
    {
      lyout <- rbind(lyout,rbind(c(rep(i,9),length(regions.to.use)+1)))
    }
    layout(lyout)
    for (region in regions.to.use)
    {
      # Failsafe in case column SAMPLE not present
      if ((! 'SAMPLE' %in% colnames(seg.mat.copy)) & ('SampleID' %in% colnames(seg.mat.copy))){
        seg.mat.copy$SAMPLE <- seg.mat.copy$SampleID
      }
      
      # rename columns
      sub.mat.copy               <- seg.mat.copy[seg.mat.copy$SAMPLE%in%region,,drop=FALSE]
      colnames(sub.mat.copy)[2]  <- 'Chromosome'
      colnames(sub.mat.copy)[3]  <- 'StartPosition'
      colnames(sub.mat.copy)[4]  <- 'EndPosition'
      region.phyloCCF            <- phylo.region.list[[region]]
      
      #pdf(early.late.pdf)
      par(mar=c(0.5,4,0.5,0.2))
      #par(mar=c(0.5,4,0.5,0.2))
      if(nrow(lyout)==1)
      {
        par(mar=c(14,4,14,0.2))
      }
      par(lend=1)
      plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                              ,TCGA.earlyLate=region.phyloCCF
                              ,sub.clonal=1
                              ,most.likely.cluster=mostLikelyClusters
      )
      mtext(unlist(strsplit(region,split="_"))[length(unlist(strsplit(region,split="_")))],side=2,cex=0.6,line=2,las=2)
    }
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    # require(RColorBrewer)
    cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    par(mar=c(2,0.5,2,0))
    if(nrow(lyout)==1)
    {
      par(mar=c(14,0.5,14,0))
    }
    
    plot(1, type="n"
         , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
    for (k in 1:no.optima)
    {
      points(0
             ,k
             ,col=cols.opac[k]
             ,bg = cols[k]
             ,pch=pch.vals[k]
             ,cex=1.5)
      text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
           ,pos=4)
    }
    
    text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
    
    lyout <- c()
    for (i in seq(1,length(regions.to.use)*2,by=2))
    {
      lyout <- rbind(lyout,rbind(c(rep(i,9),i+1)))
    }
    
    layout(lyout)
    
    # next, let's plot each cluster separately
    for (cluster in 1:max(mostLikelyClusters))
    {
      for (region in regions.to.use)
      {
        # Failsafe in case column SAMPLE not present
        if ((! 'SAMPLE' %in% colnames(seg.mat.copy)) & ('SampleID' %in% colnames(seg.mat.copy))){
          seg.mat.copy$SAMPLE <- seg.mat.copy$SampleID
        } 
        # rename columns
        sub.mat.copy               <- seg.mat.copy[seg.mat.copy$SAMPLE%in%region,,drop=FALSE]
        colnames(sub.mat.copy)[2]  <- 'Chromosome'
        colnames(sub.mat.copy)[3]  <- 'StartPosition'
        colnames(sub.mat.copy)[4]  <- 'EndPosition'
        region.phyloCCF            <- phylo.region.list[[region]]
        region.phyloCCF.cluster    <- region.phyloCCF[names(which(mostLikelyClusters==cluster)),,drop=FALSE]
        #pdf(early.late.pdf)
        par(mar=c(0.5,4,0.5,0.2))
        if(nrow(lyout)==1)
        {
          par(mar=c(14,4,14,0.2))
        }
        par(lend=1)
        plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                                ,TCGA.earlyLate=region.phyloCCF
                                ,sub.clonal=1
                                ,most.likely.cluster=mostLikelyClusters
                                ,cluster=cluster
        )
        mtext(unlist(strsplit(region,split="_"))[length(unlist(strsplit(region,split="_")))],side=2,cex=0.6,line=2,las=2)
        if(nrow(region.phyloCCF.cluster)==1)
        {
          plot(NULL, type = "n"
               , xlim = c(0, max(A$density))
               , ylim=c(-0.25,6)
               ,bty='n'
               ,xaxs='i'
               ,xaxt='n'
               ,yaxt='n'
               ,yaxs='i'
               ,xlab=""
               ,main=""
               ,ylab=""
          )
          next;
        }
        ds <- density(ifelse(as.numeric(region.phyloCCF.cluster$mutCopyNum)>5,5,as.numeric(region.phyloCCF.cluster$mutCopyNum)))
        ds1 <- ds
        ds1$x <- ds$y
        ds1$y <- ds$x
        par(mar=c(0.5,0,0.5,4))
        if(nrow(lyout)==1)
        {
          par(mar=c(14,0,14,4))
        }
        A <- hist(ifelse(as.numeric(region.phyloCCF.cluster$mutCopyNum)>5,5,as.numeric(region.phyloCCF.cluster$mutCopyNum)),breaks=seq(-0.25,6,by=0.1),plot=FALSE)
        plot(NULL, type = "n"
             , xlim = c(0, max(A$density))
             , ylim=c(-0.25,6)
             ,bty='n'
             ,xaxs='i'
             ,xaxt='n'
             ,yaxt='n'
             ,yaxs='i'
             ,xlab=""
             ,main=""
             ,ylab=""
        )
        rect(0, A$breaks[1:(length(A$breaks) - 1)], A$density, A$breaks[2:length(A$breaks)]
             ,border=TRUE,col="#CC6666")
        
        lines(ds1)       
      }
    }
  }
}


plot.simpleClusters.raw    <- function( seg.mat.patient
                                    , TCGA.earlyLate
                                    , TCGA.barcode=""
                                    , max.cpn = 5
                                    , min.probes = 10
                                    , sub.clonal = 1
                                    , min.vaf.present = 0.01
                                    , most.likely.cluster
                                    ,cluster=NA) {
  
  
  seg.mat.patient$Chromosome    <- as.numeric(as.character(seg.mat.patient$Chromosome))
  seg.mat.patient$StartPosition <- as.numeric(as.character(seg.mat.patient$StartPosition))
  seg.mat.patient$EndPosition   <- as.numeric(as.character(seg.mat.patient$EndPosition))
  seg.mat.patient$major         <- seg.mat.patient$COPY_NUMBER_A
  seg.mat.patient$min.allele    <- seg.mat.patient$COPY_NUMBER_B
  
  
  chrom.length.copy   <- fun.chrom.length(seg.mat.patient[,2:4])
  chrom.segs          <- fun.add.chrom(seg.mat.patient[,2:4],chrom.length.copy)
  seg.mat.plot        <- seg.mat.patient
  seg.mat.plot[,2:4]  <- chrom.segs
  major               <- seg.mat.plot$COPY_NUMBER_A
  seg.mat.plot        <- cbind(seg.mat.plot,major)
  
  
  # Order correctly
  TCGA.plot          <- TCGA.earlyLate
  Chromosome         <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,2])
  Start_pos          <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,3])
  TCGA.plot          <- cbind(TCGA.plot,Chromosome,Start_pos)
  TCGA.plot          <- data.frame(apply(TCGA.plot,2,unlist),stringsAsFactors=FALSE)
  
  min.x               <- 0
  min.y               <- -0.25
  max.y               <- max.cpn +1
  max.x               <- as.numeric(max(seg.mat.plot$EndPosition))
  
  # Make sure any copy numbers that are too high are still included
  seg.mat.plot$major       <- ifelse(seg.mat.plot$major>max.cpn,max.cpn,seg.mat.plot$major)
  seg.mat.plot$min.allele  <- ifelse(seg.mat.plot$min.allele>max.cpn,max.cpn,seg.mat.plot$min.allele)
  TCGA.plot$mutCopyNum      <- ifelse(TCGA.plot$mutCopyNum>max.cpn,max.cpn,TCGA.plot$mutCopyNum)
  
  plot(1
       ,xlim=c(min.x,max.x)
       ,ylim=c(min.y,max.y)
       ,type='n'
       ,xaxt='n'
       ,yaxs='i'
       ,yaxt='n'
       ,xaxs='i'
       ,ylab=""
       ,xlab=""
       ,lwd=2
       ,bty="n")
  
  box(lwd=0.5) 
  
  axis(side=2
       ,at=seq(0,max.cpn,by=1)
       ,labels = c(seq(0,max.cpn-1,by=1),paste('>',max.cpn,sep=""))
       ,las=2
       ,cex.axis=0.7
       ,lwd=0.5)
  
  mtext(text=TCGA.barcode
        ,side=3
        ,line=2
  )
  
  
  # Now plot each segment
  fun.plot.segment <- function(x,seg.mat.plot)
  {
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- seg.mat.plot[x,'major']+0.05
    y1 <- y0
    
    segments(x0,y0,x1,y1,col='black',lwd=2.4)
    
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- as.numeric(seg.mat.plot[x,'min.allele'])-0.05
    y1 <- y0
    
    segments(x0,y0,x1,y1,col='#009E73',lwd=2.4)
    
  }
  
  sapply(1:nrow(seg.mat.plot),fun.plot.segment,seg.mat.plot)
  # Draw lines to separate the different chromosomes
  for (i in sort(unique(as.numeric(seg.mat.plot$Chromosome))))
  {
    abline(v=max(as.numeric(seg.mat.plot[seg.mat.plot$Chromosome==i,'EndPosition']))
           ,lwd=0.5)
  }
  
  # We need to update the 
  TCGA.plot          <- TCGA.plot[order(as.numeric(TCGA.plot$Chromosome),as.numeric(TCGA.plot$Start_pos)),]
  TCGA.plot$Start_pos <- as.numeric(fun.add.chrom(cbind(TCGA.plot$Chromosome,TCGA.plot$Start_pos,TCGA.plot$Start_pos),chrom.length.copy)[,2])
  
  
  max.cols = 12
  # require(RColorBrewer)
  no.optima      = length(unique(most.likely.cluster))
  cols           = paste(RColorBrewer::brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
  cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
  cols.opac      = paste(cols,'99',sep="")
  
  pch.vals       = c()
  k <- 21
  m <- 1  
  while(m<=no.optima)
  {
    pch.vals    = c(pch.vals,k)
    if (m==max.cols)
    {
      k <- k+1
    }
    m <- m+1
  }
  
  pch.vals     = rev(pch.vals)
  
  
  if (is.na(cluster[1]))
  {
    for (i in 1:length(unique(most.likely.cluster)))
    {
      
      
      cluster.muts <- TCGA.plot[TCGA.plot$mutation_id%in%names(which(most.likely.cluster==unique(most.likely.cluster)[i])),,drop=FALSE]
      points(cluster.muts$Start_pos
             ,cluster.muts$mutCopyNum
             ,cex=0.7
             ,pch=pch.vals[i]
             ,col=cols[i]
             ,bg =cols.opac[i] ) # ca
      
      
      
    }
    
  }
  if (!is.na(cluster[1]))
  {
    for (i in cluster)
    {
      
      cluster.muts <- TCGA.plot[TCGA.plot$mutation_id%in%names(which(most.likely.cluster==unique(most.likely.cluster)[i])),,drop=FALSE]
      points(cluster.muts$Start_pos
             ,cluster.muts$mutCopyNum
             ,cex=0.7
             ,pch=pch.vals[as.numeric(i)]
             ,col=cols[i]
             ,bg =cols.opac[as.numeric(i)] ) # ca
      
      
      
    }
    
  }
  

  mtext(side=3
        ,at=fun.chrom.mean(seg.mat.plot[,2:4])
        ,text=sort(unique(seg.mat.plot[,2]))
        ,cex=seq(0.6,0.4,length.out=length(unique(seg.mat.plot[,2])))
        ,line=-1
        #,lwd=0.5
  )
  
  
}

mtexti <- function(text, side, off = 0.25,
                   srt = if(side == 2) 90  else
                     if(side == 4) 270 else 0, ...) {
  # dimensions of plotting region in user units
  usr <- par('usr')
  # dimensions of plotting region in inches
  pin <- par('pin')
  # user units per inch
  upi <- c(usr[2]-usr[1],
           usr[4]-usr[3]) / pin
  # default x and y positions
  xpos <- (usr[1] + usr[2])/2
  ypos <- (usr[3] + usr[4])/2
  if(1 == side)
    ypos <- usr[3] - upi[2] * off
  if(2 == side)
    xpos <- usr[1] - upi[1] * off
  if(3 == side)
    ypos <- usr[4] + upi[2] * off
  if(4 == side)
    xpos <- usr[2] + upi[1] * off
  text(x=xpos, y=ypos, text, xpd=NA, srt=srt, ...)
}