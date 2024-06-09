########################
loadassays <- function(){
  bVals <- readRDS("data/bVals081820.rds")
  mVals <- readRDS("data/mVals081820.rds")
  ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  baseDir <- file.path(getwd(),"meth_081720")
  mVals1 <- readRDS("data/mVals.rds")
  bVals1 <- readRDS("data/bVals.rds")
  targets1 <- readRDS("data/targets.rds")
  man1 <- readRDS("data/methylationmanifest.rds")
  targets <- createtargets(baseDir)
  assays <- list(bVals,mVals,ann850k,baseDir,mVals1,bVals1,targets1,targets,man1)
  names(assays) <- c("bVals","mVals","ann850k","baseDir","mVals1","bVals1","targets1",
                     "targets","man1")
  assays
}

###########################
mergetargets <- function(targets,targets1){
  targets1$Sample_subtype[targets1$Sample_subtype == "colorectal"] <- "Colorectal mucinous"
  targets1$Sample_subtype[targets1$Sample_subtype == "ovarian"] <- "Ovarian mucinous"
  targets1$Sample_subtype[targets1$Sample_subtype == "pancreas"] <- "Pancreas mucinous"
  targets1$Sample_subtype[targets1$Sample_subtype == "stomach"] <- "Stomach mucinous"
  colnames(targets1)[9] <- "T.N"
  colnames(targets1)[10] <- "Diagnosis"
  colnames(targets1)[11] <- "Sample"
  colnames(targets1)[12] <- "Tissue"
  targets1$Other.ID <- rep(NA, length(targets1$Sample_Name))
  targets1$sampletype <- paste0(targets1$Diagnosis, "_", targets1$T.N)
  targets1$batch <- rep(1,length(targets1$Sample_Name))
  targets$batch <- rep(2,length(targets$Sample_Name))
  targets1 <- targets1[match(colnames(targets),colnames(targets1))]
  targets <- rbind(targets,targets1)
}

########################
mergemVals <- function(mVals,mVals1, targets){
  t1 <- mVals
  t2 <- mVals1
  t1 <- as.data.frame(t1)
  t2 <- as.data.frame(t2)
  t1$row.names <- row.names(t1)
  t2$row.names <- row.names(t2)
  t3 <- left_join(t1, t2, by = "row.names")
  rownames(t3) <- t3$row.names
  t3 <- t3[!is.na(t3$CGOV474N),]
  t3 <- t3[!is.na(t3$CGOV169N),]
  t3 <- t3[,match(targets$Sample_Name,colnames(t3))]
  mVals <- as.matrix(t3)
}

#######################
##left_join
mergebVals <- function(bVals,bVals1,targets){
  t1 <- bVals
  t2 <- bVals1
  t1 <- as.data.frame(t1)
  t2 <- as.data.frame(t2)
  t1$row.names <- row.names(t1)
  t2$row.names <- row.names(t2)
  t3 <- left_join(t1, t2, by = "row.names")
  rownames(t3) <- t3$row.names
  t3 <- t3[!is.na(t3$CGOV474N),]
  t3 <- t3[!is.na(t3$CGOV169N),]
  t3 <- t3[,match(targets$Sample_Name,colnames(t3))]
  bVals <- as.matrix(t3)
}

################
createtbVals <- function(bVals, targets){
  tbVals <- bVals[,colnames(bVals) %in%
                    targets$Sample_Name[targets$T.N == "T"]]
}

###############
createtmVals <- function(mVals, targets){
  tmVals <- mVals[,colnames(mVals) %in% targets$Sample_Name[targets$T.N == "T"]]
}

#############
createbitbVals <- function(tbVals, cutoff.beta=0.2){
  bintbVals <- tbVals
  bintbVals[bintbVals < cutoff.beta] <- 0L
  bintbVals[bintbVals > cutoff.beta] <- 1L
  ##print(bintbVals)
  bintbVals
}

###############
createbitmVals <- function(tmVals, cutoff.M=-1.5){
  bintmVals <- tmVals
  bintmVals[bintmVals < cutoff.M] <- 0L
  bintmVals[bintmVals > cutoff.M] <- 1L
  ##print(bintmVals)
  bintmVals
}

############
creategann850ksub <- function(ann850k, bVals){
  ann850kSub <- ann850k[rownames(ann850k) %in% rownames(bVals),]
  ann850kSub <- ann850kSub[match(rownames(bVals), ann850kSub$Name),]
  gann850ksub <- makeGRangesFromDataFrame(ann850kSub,
                                          keep.extra.columns = T,
                                          seqinfo = NULL,
                                          seqnames.field = "chr",
                                          start.field = "pos",
                                          end.field = "pos",
                                          strand.field = "strand")
}


#################################
createtargets <- function(baseDir){
  targets <- read.metharray.sheet(baseDir)
  samp <- targets$Sample_Name
  sampdat <- read.csv("extdata/methdat082620.csv")
  sampdat <- sampdat[,c(1,3:7)]
  colnames(sampdat)[1] <- "Sample_Name"
  targets <- merge(targets, sampdat, by = "Sample_Name")
  targets$sampletype <- paste0(targets$Diagnosis,"_",targets$T.N)
  return(targets)
}

##################
createttargets <- function(targets){
  ttargets <- targets[targets$T.N == "T",]
}

########################
clean.targets <- function(targets){
  samp <- unlist(lapply(strsplit(targets$Sample_Name,"_"),function(x) tail(strsplit(x,split=" ")[[2]],1)))
  samp[31] <- "CGPA365T"
  samp[38] <- "CGST1N"
  samp[3] <- "CGCRC330T1"
  targets$Sample_Name <- samp
  targets$Sample_Name[12] <- "CGOV179T"
  return(targets)
}

#####################################
clean.manifest <- function(man3){
  colnames(man3) <- c("Sample_type","Sample_subtype","Sample_Source","Sample_Name","Sample_Origin")
  man3$Sample_Name <- as.character(man3$Sample_Name)
  man3$Sample_Name[3] <- "CGCRC330T1"
  man3$Sample_Name[12] <- "CGOV179T"
  man3$Sample_Name[31] <- "CGPA365T"
  man3$Sample_Name[38] <- "CGST1N"
  man3$Sample_Source <- as.character(man3$Sample_Source)
  man3$Sample_Source[12] <- "CGOV179"
  man3$Sample_Source[3] <- "CGCRC330"
  man3$Sample_Source[38] <- "CGST1"
  man3$Sample_Source[31] <- "CGPA365"
  return(man3)
}

###############################
clean.targets <- function(targets){
  samp <- unlist(lapply(strsplit(targets$Sample_Name,"_"),function(x) tail(strsplit(x,split=" ")[[2]],1)))
  samp[31] <- "CGPA365T"
  samp[38] <- "CGST1N"
  samp[3] <- "CGCRC330T1"
  targets$Sample_Name <- samp
  targets$Sample_Name[12] <- "CGOV179T"
  return(targets)
}

###################
createseobject <- function(assays, cutoff.beta=0.2, cutoff.M=-1.5){
  assays$targets1 <- clean.targets(assays$targets1)
  assays$man1 <- clean.manifest(assays$man1)
  assays$targets1 <- merge(assays$targets1,assays$man1,by = "Sample_Name")
  assays$targets <- mergetargets(assays$targets,assays$targets1)
  assays$mVals <- mergemVals(assays$mVals,assays$mVals1, assays$targets)
  assays$bVals <- mergebVals(assays$bVals, assays$bVals1, assays$targets)

  ## tumor beta values
  assays$tbVals <- createtbVals(assays$bVals, assays$targets)
  ## tumor M values
  assays$tmVals <- createtmVals(assays$mVals, assays$targets)
  ## binary tumor beta Values
  assays$bintbVals <- createbitbVals(assays$tbVals, cutoff.beta)
  ## binary tumor M values
  assays$bintmVals <- createbitmVals(assays$tmVals, cutoff.M)
  assays$gann850ksub <- creategann850ksub(assays$ann850k, assays$tbVals)
  assays$gann850ksuborig <- creategann850ksub(assays$ann850k, assays$bVals)
  assays$ttargets <- createttargets(assays$targets)
  assays$bibVals <- createbitbVals(assays$bVals)
  assays$bimVals <- createbitmVals(assays$mVals)
  assays
}


#############################################
createoriginalSE <- function(assays){
  se = SummarizedExperiment(assays=SimpleList(beta= assays$bVals, M = assays$mVals,
                                              bibeta = assays$bibVals,
                                              biM = assays$bimVals),
                            rowRanges = assays$gann850ksuborig,
                            colData = assays$targets)
}


##################################################
createSummarizedExperiment <- function(assays){
    sl <- SimpleList(beta= assays$tbVals,
                     M= assays$tmVals,
                     bibeta = assays$bintbVals,
                     bimVals = assays$bintmVals)
    se = SummarizedExperiment(assays= sl,
                            rowRanges = assays$gann850ksub,
                            colData = assays$ttargets)
}

M <- function(se) assays(se)[["M"]]
beta <- function(se) assays(se)[["beta"]]
binary <- function(se) assays(se)[["bibeta"]]

isEndometrioid <- function(se) colData(se)$Diagnosis == "Ovarian endometrioid"

prop.test2 <- function(x, y) tidy(prop.test(x, y))[c("statistic", "p.value")]

proptest.methyl <- function(se){
    numMeth1 <- rowSums(binary(se)[,isEndometrioid(se)])
    numMeth2 <- rowSums(binary(se)[, !isEndometrioid(se)])
    n1 <- sum(isEndometrioid(se))
    n2 <- ncol(se)-n1
    x <- cbind(numMeth1, numMeth2)
    y <- c(n1, n2)
    results <- apply(x, 1, prop.test2, y=y)
    results <- rbindlist(results)
    colnames(results) <- c("prop.test","prop.p.value")
    results
}

wilcox.test2 <- function(x,y) tidy(wilcox.test(x~y))[c("statistic", "p.value")]

wilcoxtest.methyl <- function(se){
    x <- beta(se)
    y <- colData(se)$sampletype
    results <- apply(x, 1 ,wilcox.test2, y = y)
    results <- rbindlist(results)
    colnames(results) <- c("wilcox.test","wilcox.p.value")
    results
}
