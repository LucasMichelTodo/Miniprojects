  ########################################################################################
  ############################## Modify this part ########################################
  ########################################################################################

  datadir <- "/media/lucas/Disc4T/Projects/Miniprojects/Microarray_R_analysis/Raw_Data/"
  outdir <- "/media/lucas/Disc4T/Projects/Miniprojects/Microarray_R_analysis/R_results_NewTry"
  annot <- read.csv("/media/lucas/Disc4T/Projects/Miniprojects/Microarray_R_analysis/Files/array_anotation.csv",
                    sep="\t", header = F)
  infiles <- "/media/lucas/Disc4T/Projects/Miniprojects/Microarray_R_analysis/Files/"

  nsamples <- 9
  sample_names <- c("Ctl_1",
                    "Ctl_2",
                    "Ctl_3",
                    "A_1",
                    "A_2",
                    "A_3",
                    "B_1",
                    "B_2",
                    "B_3")

  times <- c(1,2,3, 1,2,3, 1,2,3)
  types <- c("ctl", "ctl", "ctl", "treatA", "treatA", "treatA", "treatB", "treatB", "treatB")

  Ctl_1 <- read.table(paste0(datadir, "ctl_tp1.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
  Ctl_2 <- read.table(paste0(datadir, "ctl_tp2.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
  Ctl_3 <- read.table(paste0(datadir, "ctl_tp3.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)


  A_1 <- read.table(paste0(datadir, "treatA_tp1.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
  A_2 <- read.table(paste0(datadir, "treatA_tp2.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
  A_3 <- read.table(paste0(datadir, "treatA_tp3.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

  B_1 <- read.table(paste0(datadir, "treatB_tp1.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
  B_2 <- read.table(paste0(datadir, "treatB_tp2.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
  B_3 <- read.table(paste0(datadir, "treatB_tp3.txt"), sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

  array_list <- list(Ctl_1,
                     Ctl_2,
                     Ctl_3,
                     A_1,
                     A_2,
                     A_3,
                     B_1,
                     B_2,
                     B_3)

run_all_plots <- "no"

  ########################################################################################
  ############################ End of modifiable part ####################################
  ########################################################################################

## Import libraries

  print("Importing Libraries...")
  list.of.packages <- c("reshape2", "ggfortify", "tidyverse", "RColorBrewer", "sp")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  library(reshape2)
  library(ggfortify)
  library(tidyverse)
  library(RColorBrewer)
  library(sp)

  if (!("Biobase" %in% installed.packages()[,"Package"])){
      if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
      BiocManager::install()
  }

  if (!("org.Pf.plasmo.db" %in% installed.packages()[,"Package"])){
      if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
      BiocManager::install("org.Pf.plasmo.db")
  }

  library(org.Pf.plasmo.db)
  library(Biobase)

## Create folders for output

  print("Creating folders for Output...")
  dir.create(paste0(outdir))
  dir.create(paste0(outdir, "/Plots"))
  dir.create(paste0(outdir, "/Plots/Array_Plots"))
  dir.create(paste0(outdir, "/Plots/MA_Plots"))
  dir.create(paste0(outdir, "/Plots/Time_estim/"))
  dir.create(paste0(outdir, "/Plots/Ratio"))
  dir.create(paste0(outdir, "/Plots/Ratio/Gene_Level"))
  dir.create(paste0(outdir, "/Plots/Ratio/Probe_Level"))
  dir.create(paste0(outdir, "/Plots/Ratio/Gene_Level"))
  dir.create(paste0(outdir, "/Plots/Red_Signal"))
  dir.create(paste0(outdir, "/Plots/Red_Signal/Gene_Level"))
  dir.create(paste0(outdir, "/Plots/Red_Signal/Probe_Level"))
  figPath <- paste0(outdir, "/Plots/")

## Load Array annotation, Gene-list and Variant Genes.
## Load ~array_anotation.csv~ a table with the array annotation.
## Load ~gene_list.txt~ a list of genes present in the array.
## This will be used for calculating the mean intensity in the array
## for excluding low expressed genes.

  print("Loading Array Annotation...")
  gene_list <- readLines(paste0(infiles, "gene_list.txt"))
  variant <- read.csv(paste0(infiles, "Gens_variants_extended.txt"), header = TRUE, sep = "\t")

## Create Probe-DF

  print("Creating Probe DF...")
  probe_df <- array_list[[1]][,c(7,11,14,15)]

  getCols <- function(df){
      return(df[,c(11,14,15)])
  }

  goodCols <- lapply(array_list[2:nsamples], function(x) getCols(x))
  df <- do.call("cbind", goodCols)

  probe_df <- cbind(probe_df, df)

  probe_df["Gene_id"] <- annot$V2
  probe_df["name"] <- annot$V4
  probe_df["Annot"] <- annot$V5

  probe_df["Annot"] <- gsub("Plasmodium", "Pl.", probe_df$Annot)
  probe_df["Annot"] <- gsub("protein", "prot.", probe_df$Annot)
  probe_df["Annot"] <- gsub("membrane", "memb.", probe_df$Annot)
  probe_df["Annot"] <- gsub("conserved", "cvd.", probe_df$Annot)
  probe_df["Annot"] <- gsub("function", "func.", probe_df$Annot)
  probe_df["Annot"] <- gsub("unknown", "ukwn.", probe_df$Annot)
  probe_df["Annot"] <- gsub("exported", "xptd.", probe_df$Annot)
  probe_df["Annot"] <- gsub("pseudogene", "pseudo", probe_df$Annot)
  probe_df["Annot"] <- gsub("putative", "put.", probe_df$Annot)
  probe_df["Annot"] <- gsub("%2C", "", probe_df$Annot)

  # Remove probes that map tu multiple genes.
  probe_df <- probe_df[annot$V3 != "drop" | is.na(annot$V3),]

  # Add Variant Genes information
  probe_df["Variant"] <- probe_df$Gene_id %in% variant$ID

## Group Columns

  signalCols <- nsamples*3+1
  allcols <- dim(probe_df)[2]

  ratioCols <- seq(2, signalCols, 3)
  redCols <- seq(4, signalCols, 3)
  greCols <- seq(3, signalCols, 3)

  infoCols <- c(1, (signalCols+1):allcols)

## Remove low expression probes

  print("Removing low expression probes...")
  medians <- list()
  for (i in 1:nsamples){
      redm <- median(sort(probe_df[probe_df$Gene_id %in% gene_list, redCols][,i])[1:100])
      grenm <- median(sort(probe_df[probe_df$Gene_id %in% gene_list, greCols][,i])[1:100])
      medians[[i]] <- c(grenm, redm)
  }

  passTest <- list()
  for(i in 1:nsamples){
      g <- probe_df[,greCols][i] < 3*medians[[i]][1]
      r <- probe_df[,redCols][i] < 3*medians[[i]][2]
      all <- g & r
      passTest[[i]] <- !all
  }

  testDF <- as.data.frame(passTest)
  pass <- rowSums(testDF) > 0
  write.csv(table(!pass), paste0(outdir, "/NA_probes.csv"))
  probe_df[!pass,c(ratioCols)] <- NA

## Array Plots

  print("Plotting Arrays...")
  cols <- rev(brewer.pal(11, 'Spectral'))

  arrayPlot <- function(df) {
      df_name <- sample_names[i]
      p1 <- qplot(Col, Row, data=df, color=log2(rMedianSignal)<7) + scale_color_manual(values=c("aliceblue", "black")) + ggtitle(df_name)
      ggsave(p1, filename = paste0(figPath, "Array_Plots/sample_", df_name, "_boolean.jpeg"), device = "jpeg")

      p2 <- qplot(Col, Row, data=df, color=log2(rMedianSignal)) + scale_colour_gradientn(colours = cols) + ggtitle(df_name)
      ggsave(p2, filename = paste0(figPath, "Array_Plots/sample_", df_name, ".jpeg"), device = "jpeg")

      p3 <- qplot(Col, Row, data=df, color=is.na(LogRatio)) + scale_color_manual(values=c("aliceblue", "red")) + ggtitle(df_name)
      ggsave(p3, filename = paste0(figPath, "Array_Plots/sample_", df_name, "_NAs.jpeg"), device = "jpeg")
  }

  for (i in 1:length(array_list)) {arrayPlot(array_list[[i]])}

## MA Plots

  print("Plotting MA Plots...")
  myMAplot  <- function(mray){
        df_name <- sample_names[i]
        m_vals <- log2(mray$gProcessedSignal) - log2(mray$rProcessedSignal)
        a_vals <- (log2(mray$gProcessedSignal) + log2(mray$rProcessedSignal))/2
        ma_df <- cbind(a_vals, m_vals)
        p <- ggplot(ma_df, aes(x=a_vals, y=m_vals))
        p <- p + geom_point()
        p <- p + geom_smooth(method = "lm", se=F, color= "red")
        p <- p + geom_hline(yintercept=0, color = "blue", size = 1)
        ggsave(p, filename = paste0(figPath, "MA_Plots/sample_", df_name, "_MA.jpeg"), device = "jpeg")
    }

    for (i in 1:length(array_list)) {myMAplot(array_list[[i]])}

## Change to Log2 (Log Ratio cols, originally log10)

  print("Creating Eset...")
  probe_df[,ratioCols] <- log2(10**probe_df[,ratioCols])

## Change to Log2 (Raw signal Cols, originally unlogged)

  probe_df[,redCols]  <- log2(probe_df[,redCols])

## Create eSet: xprobe

  exprsx <- as.matrix(probe_df[,ratioCols])
  colnames(exprsx) <- sample_names
  fdata <- new("AnnotatedDataFrame", probe_df[,infoCols])
  teor_time <- times
  type <- types
  pdata <- data.frame(type=type, teor_time=teor_time); rownames(pdata) <- sample_names
  pdata <- new("AnnotatedDataFrame", pdata)
  xprobe <- new("ExpressionSet", exprs=exprsx, featureData=fdata, phenoData=pdata)
  save(xprobe,file=paste0(outdir, '/probeLevel.RData'))

  exprsx <- as.matrix(probe_df[,redCols])
  colnames(exprsx) <- sample_names
  xprobe_red <- new("ExpressionSet", exprs=exprsx, featureData=fdata, phenoData=pdata)

  write.csv(cbind(fdata@data, exprs(xprobe)), paste0(outdir, "/probeLevel_exp.csv"), row.names=F)
  write.csv(cbind(fdata@data, exprs(xprobe_red)), paste0(outdir, "/probeLevel_redSignalexp.csv"), row.names=F)

## Rename and Summarize
## Summarize all probes for a same gene, using median polish.
## Here we are considering the ratio between green and red signal as the expression
## value.

  myRma <- function(x) {
      if (class(x)=='numeric') {
          ans <- x
      } else {
          ans <- medpolish(x,trace.iter=FALSE,na.rm=TRUE)
          ans <- ans$overall + ans$col
      }
      return(ans)
  }

  renameGenesAndSummarize <- function(genesToRename.sd,exprsx,geneid,summaryMethod=myRma,type) {

      if (type == "ratio"){
          xgene <- by(exprsx[,ratioCols],geneid,myRma)
      } else if (type == "red"){
          xgene <- by(exprsx[,redCols],geneid,myRma)
      }

      xgene <- do.call('rbind',xgene)

      mysd <- function(x) { ans <- ifelse(sum(!is.na(x))==1,0,sd(x,na.rm=TRUE)); return(ans) }
      sdgene <- aggregate(exprsx[, ratioCols],by=list(geneid),FUN=mysd)

      names(sdgene)[1] <- 'geneid'
      xgene <- data.frame(geneid=rownames(xgene),xgene); rownames(xgene) <- NULL

      fdata <- by(exprsx[,(signalCols+1):allcols],geneid,unique)

      genenames <- names(fdata)
      fdata <- do.call('rbind',fdata)

      fdata <- new("AnnotatedDataFrame", data.frame(fdata))
      rownames(fdata) <- as.character(xgene$geneid)

      exprsxgene <- as.matrix(xgene[,-1])
      rownames(exprsxgene) <- as.character(xgene$geneid);
      colnames(exprsxgene) <- sample_names
      eset <- new("ExpressionSet",exprs=exprsxgene, featureData=fdata, phenoData=pdata)
      return(list(eset=eset,sdgene=sdgene,fdata=fdata,geneid=geneid))
  }

  geneid <- probe_df$Gene_id
  geneid <- as.character(geneid)
  genesToRename.sd <- NA

  tmp <- renameGenesAndSummarize(genesToRename.sd=genesToRename.sd,exprsx=probe_df,geneid=geneid,summaryMethod=myRma, type="ratio")
  xgene <- tmp[['eset']]; sdgene <- tmp[['sdgene']]; fdata <- tmp[['fdata']]; geneid <- tmp[['geneid']]

  tmp2 <- renameGenesAndSummarize(genesToRename.sd=genesToRename.sd,exprsx=probe_df,geneid=geneid,summaryMethod=myRma, type="red")
  xgene_red <- tmp2[['eset']]; sdgene <- tmp2[['sdgene']]; fdata <- tmp2[['fdata']]; geneid <- tmp2[['geneid']]

## Estimate times

  print("Estimating times...")
  bozdechPath <- paste0(infiles, 'bozdech_Hb3_clean2.csv')
  LemieuxFunctionsPath <- paste0(infiles, 'lemieux_et_al_pipeline_functions.r')

  getTimeEstimation <- function(x,dataPath,functionsPath,figuresPath,B=100) {
                                          #  x: the expressionSet for which we want to estimate times (our data).
                                          #  dataPath: path to data that will be used to estimate timepoints (from Bozdech et al)
                                          #  functionsPath: path to the script containing the functions from Lemieux's paper.
                                          #  figuresPath: where we want to save the output plots.
      source(functionsPath)
                                          #  z <- read.csv(dataPath, as.is = T,sep='\t')
      z <- read.csv(dataPath, as.is = T)
                                          #  colnames(z)[1] <- 'Name'
                                          #  oldTime <- as.numeric(as.character(pData(x)$time))
      oldTime <- as.numeric(teor_time)
      x <- exprs(x)
      x <- data.frame(Name=as.character(rownames(x)),x,stringsAsFactors=FALSE); rownames(x) <- NULL
      data <- sync_data(x, z)
      x <- data[[1]]
      z <- data[[2]]
      x <- ordinal(x, use.name = T)
      z <- ordinal(z, use.name = T)
                                          #  z.na <- cbind(z[,1:22], rep(NA, nrow(z)), z[,23:27], rep(NA, nrow(z)), z[,28:56])
      z.na <- cbind(z[,1:22], rep(NA, nrow(z)), z[,23:27], rep(NA, nrow(z)), z[,28:46])
      z <- t(apply(z.na, 1, smooth.missing))
      sigma.epsilon <- 789.9056
      z.smooth <- smooth.ref(z, method = "spline", spar = 0.5)
      z.smooth.hourly <- z.smooth[,ll.par$hourly]
                                          #  sigma.eta <- mean(sd(z[,11:ncol(z)] - z.smooth.hourly, na.rm = T), na.rm=T)
      sigma.eta <- mean(sd(z - z.smooth.hourly, na.rm = T), na.rm=T)
      new.sigma <- sqrt(sigma.eta^2 + sigma.epsilon^2)
      ll <- compute.ll(x = x, z = z.smooth, sigma = new.sigma, bootstrap = T, B = B, sample.rate = 0.50)
      myTimes <- mle(ll)
      png(file.path(figuresPath,'/Time_estim/defaultPlots1.png'))
      plot.ll(ll)
      dev.off()
      png(file.path(figuresPath,'/Time_estim/defaultPlots2.png'))
      plot.mle(ll)
      dev.off()
      png(file.path(figuresPath,'/Time_estim/ownPlots1.png'))
      plot(density(myTimes),main='Estimated times density')
      dev.off()
      png(file.path(figuresPath,'/Time_estim/ownPlots2.png'))
      plot(oldTime, as.numeric(myTimes),xlab='Old times',ylab='Estimated times',xlim=c(-5,50),ylim=c(-5,50))
      abline(0,1,col=2,lwd=2)
      abline(v=oldTime,lwd=0.5,lty=3)
      dev.off()
      return(myTimes)
  }

  ascendingTime <- function(x){
      current <- 0
      ncycle <- 0
      for (i in 1:length(x)){
          val  <- x[i]+(48*ncycle)
          if (val < current){
              current <- val
              val <- val+48
              ncycle <- ncycle+1
          }
          current <- val
          x[i] <- val
      }
      return(x)
  }


  estimatedTimes <- getTimeEstimation(xgene,bozdechPath,LemieuxFunctionsPath,file.path(figPath),B=100)
  estimatedTimes[estimatedTimes < 0] <- 0
  hpi <- estimatedTimes

  for (type in pData(xgene)$type){
      sel <- pData(xgene)$type == type
      typetime <- estimatedTimes[sel]
      time <- ascendingTime(typetime)
      estimatedTimes[sel] <- time
  }

  write.csv(estimatedTimes, paste0(outdir, "/Estimated_Times.csv"))
  pData(xgene)$time <- estimatedTimes
  pData(xgene_red)$time <- estimatedTimes
  pData(xgene)$hpi <- hpi
  pData(xgene_red)$hpi <- hpi

  # Save ExpressionSet at gene level
  save(xgene,file=paste0(outdir, '/geneLevel.RData'))
  save(xgene_red,file=paste0(outdir, '/geneLevel_redSignal.RData'))

  # Boxplot after summarization
  pdf(file.path(figPath,'boxplot_afterSummarization.pdf'))
  boxplot(exprs(xgene),main='summarization method: median poslish')
  dev.off()

## Save results in CSVs

  write.csv(xgene@phenoData@data, file = paste0(outdir, "/experiment_data.csv"))

  write.csv(cbind(xgene@featureData@data, exprs(xgene)), paste0(outdir, "/geneLevel_exp.csv"), row.names=F)
  write.csv(cbind(xgene_red@featureData@data, exprs(xgene_red)), paste0(outdir, "/geneLevel_redSignal_exp.csv"), row.names=F)

## Probe Level

  print("Calculating Fold-Changes...")
  filter <- function(x,y) {x==y}
  combs <- cross2(sample_names, sample_names, .filter = filter)

  fc <- list()
  for (i in combs) {
      fc[[paste0(i[[1]], "_",  i[[2]])]] <- exprs(xprobe)[,i[[1]]] - exprs(xprobe)[,i[[2]]]
  }

  probe_fc <- as.data.frame(fc)
  probe_fc <- cbind(fData(xprobe), probe_fc)

  write.csv(probe_fc, file = paste0(outdir, "/proveLevel_FC.csv"), row.names = F)

## Gene Level

  filter <- function(x,y) {x==y}
  combs <- cross2(sample_names, sample_names, .filter = filter)

  fc <- list()
  for (i in combs) {
      fc[[paste0(i[[1]], "_",  i[[2]])]] <- exprs(xgene)[,i[[1]]] - exprs(xgene)[,i[[2]]]
  }

  gene_fc <- as.data.frame(fc)
  gene_fc <- cbind(fData(xgene), gene_fc)

  write.csv(gene_fc, file = paste0(outdir, "/geneLevel_FC.csv"), row.names = F)

## Functions

  imputePoint <- function(xs, ys, tp){

    ## "xs" and "ys" must be two vectors of equal length
    ## with the corresponding y(expression) and x(timepoint)
    ## values that form the expression plot of interest (one gene).
    ## "tp" must be the timepoint to impute.
    ## If the timepoint to be imputed is already present, leave it as is.
    ## Returns NA if missing the previous or next tp

    if (tp %in% xs){

        idx <- which(xs == tp)
        imputed <- list(x=xs[idx], y=ys[idx])

    } else {

        before <- which(xs == max(xs[xs < tp]))
        after <- which(xs == min(xs[xs > tp]))

        if (is.na(ys[before]) | is.na(ys[after])){

            imputed <- list(x=tp, y=NA)

        } else {

            x <- c(xs[c(before, after)])
            y <- c(ys[c(before, after)])

            imputed <- approx(x, y, xout=tp)
        }
    }
    return(imputed)
}

computeArea <- function(eset){

    ## Takes an eset and computes areas.
    ## pData(eset) must contain a field named "time" with the time-points.
    ## pData(eset) must have a field named "type" with the grouping variable.

    ## Set needed variables
    types <- unique(phenoData(eset)$type)
    type <- phenoData(eset)$type
    times <- phenoData(eset)$time

    maxminTP <- max(sapply(types,function(x) min(times[type==x])))
    minmaxTP <- min(sapply(types,function(x) max(times[type==x])))
    mybreaks <- seq(maxminTP, minmaxTP, length.out=5)
    tp1 <- mybreaks[2]
    tp2 <- mybreaks[3]
    tp3 <- mybreaks[4]

    xsList <- c()
    for (type in types){
        xsList <- c(xsList, list(pData(eset)$time[phenoData(eset)$type == type]))
    }

    ## Main loop
    all_areas <- c()
    for (i in 1:dim(eset)[1]){

        gene <- fData(eset)$geneID[i]

        ysList <- c()
        for (type in types){
            ysList <- c(ysList, list(exprs(eset)[i, phenoData(eset)$type == type]))
        }

        ## Estimate points where needed
        dfs <- list()
        for (i in 1:length(xsList)){

            x <- unlist(xsList[[i]])
            y <- unlist(ysList[[i]])

            points <- as.data.frame(cbind(x, y))
            midpoints <- points[points$x > maxminTP &
                                points$x < minmaxTP, ]

            first <- imputePoint(x, y, maxminTP)
            last <- imputePoint(x, y, minmaxTP)
            p1 <- imputePoint(x, y, tp1)
            p2 <- imputePoint(x, y, tp2)
            p3 <- imputePoint(x, y, tp3)

            impPoints <- rbind(first, last, p1, p2, p3)
            allpoints <- rbind(midpoints, impPoints)
            allpoints$x <- as.numeric(allpoints$x)
            allpoints$y <- as.numeric(allpoints$y)

            ordered <- arrange(allpoints, allpoints$x)

            dfs[[i]] <- ordered
        }

        ## Calculate minY on estimated DFs
        minY <- min(sapply(dfs, function(df) min(df$y, na.rm = T)))

        rowareas <- c()
        for (df in dfs){

            df$y <- df$y - minY

            ## Whole polygon from expression data
            polDF <- rbind(df,
                           c(minmaxTP, 0),
                           c(maxminTP, 0),
                           c(df[1,]))

            ## Create Polygons
            leftHalf <- rbind(df[which(df$x <= tp2),],
                              c(tp2, 0),
                              c(maxminTP, 0),
                              c(df[1,]))

            rightHalf <- rbind(df[which(df$x >= tp2),],
                               c(minmaxTP, 0),
                               c(tp2, 0),
                               c(df[df$x == tp2,]))

            mid <- rbind(df[which(df$x >= tp1 & df$x <= tp3),],
                         c(tp3, 0),
                         c(tp1, 0),
                         c(df[df$x == tp1,]))

            sides <- rbind(df[which(df$x <= tp1),],
                           c(tp1, 0),
                           c(tp3, 0),
                           df[which(df$x >= tp3),],
                           c(minmaxTP, 0),
                           c(maxminTP, 0),
                           df[1,])

            pols <- list(leftHalf, rightHalf, mid, sides)

            calcArea <- function(x) {ifelse(any(is.na(x)), NA, Polygon(x)@area)}
            areas <- unlist(lapply(pols, function(x) calcArea(x)))

            rowareas <- c(rowareas, areas)

            ## Plot polygons (for debugging purposes)
            ##pol <- Polygon(polDF)
            ##ps = Polygons(list(pol),1)
            ##sps = SpatialPolygons(list(ps))
            ##plot(sps)
        }
        all_areas <- c(all_areas, list(rowareas))
    }

    ## Set row and col names for output
    areaDF <- do.call(rbind, all_areas)
    titles <- c("Left", "Right", "Middle", "Sides")

    cols <- c()
    for (i in types){
        for (t in titles){
            name <- paste0(i, "_", t)
            cols <- c(cols, name)
        }
    }
    colnames(areaDF)  <- cols
    rownames(areaDF) <- rownames(exprs(eset))
    return(areaDF)
}

## Calls: Compute Areas

  print("Computing Areas...")
  areasDF <- computeArea(xgene)
  xout <- cbind(fData(xgene), areasDF)

  write.csv(xout, paste0(outdir, "/area_geneLevel.csv"), row.names=F)

## Custom Max

  myMax <- function(x){
      y <- abs(x)
      if (all(is.na(y))){
          return(list(NA, NA))
      } else {
          pos <- which.max(y)
          times <- c("Left", "Right", "Mid", "Sides")
          return(list(maxVal = x[pos],
                      maxTime = times[pos]))
      }
  }

## Calculate Differences and Max

  ## Get number of categories.
  n <- length(levels(pData(xgene)$type))

  ns <- list()
  i <- 1
  while (i < n+1){
      ns[[i]] <- 1:4+(4*(i-1))
      i <- i+1
  }

  ## Convert de areasDF into a list of DFs separated by types.
  areas <- list()
  for (i in 1:length(ns)) {
      areas[[i]] <- areasDF[,ns[[i]]]
  }


  ## Calculate area differences

  types <- unique(phenoData(xgene)$type)
  type <- phenoData(xgene)$type
  times <- phenoData(xgene)$time

  maxminTP <- max(sapply(types,function(x) min(times[type==x])))
  minmaxTP <- min(sapply(types,function(x) max(times[type==x])))
  mybreaks <- seq(maxminTP, minmaxTP, length.out=5)

  sets <- 1:length(areas)
  combs <- combn(sets, 2)
  span <- mybreaks[3] - mybreaks[1]
  titles <- c("Left", "Right", "Middle", "Sides")


  areaDifs <- list()
  for (i in 1:dim(combs)[2]){

      one <- combs[1,i]
      two <- combs[2,i]

      dif1 <- as.data.frame((areas[[one]] - areas[[two]])/span)
      names  <- paste0(colnames(areas[[one]]),
                       "_minus_",
                       colnames(areas[[two]]))

      prefix <- paste0(strsplit(colnames(areas[[one]])[1], "_")[[1]][1],
                       "-",
                       strsplit(colnames(areas[[two]])[1], "_")[[1]][1])
      names <- sapply(titles, function(x) paste(prefix, x, sep = "_"))

      colnames(dif1) <- names

      maxval <- apply(dif1, 1, function(x) myMax(x)[[1]])
      maxtime <- apply(dif1, 1, function(x) myMax(x)[[2]])

      mv <- paste0(prefix, "_MaxVal")
      mt <- paste0(prefix, "_MaxTime")

      dif1[mv] <- maxval
      dif1[mt] <- maxtime

      #dif2 <- -dif1

      ## prefix <- paste0(strsplit(colnames(areas[[two]])[1], "_")[[1]][1],
      ##                  "-",
      ##                  strsplit(colnames(areas[[one]])[1], "_")[[1]][1])
      ## names <- sapply(titles, function(x) paste(prefix, x, sep = "_"))

      ## colnames(dif2) <- names

      areaDifs <- c(areaDifs, list(dif1))#, list(dif2))
  }

  allDifs <- do.call(cbind, areaDifs)
  head(allDifs)

## Convert Partitions into life stages
## We consider the the timepoint in the middle of a partition (left,right,mid,sides) as it's timepoint.

  timetostage <- function(tp){
    while(tp > 48){
      tp = tp-48
    }
    return(tp)
  }

  getStage <- function(tp) {
      if ((tp >= 0) & (tp < 26)) {stg = "ring"}
      else if ((tp >= 26) & (tp < 38)) {stg = "troph"}
      else if ((tp >= 38) & (tp <= 48)) {stg = "schizont"}
      return(stg)
  }

  fracToStage <- function(frac, tps){
      if (is.na(frac)){
          return(NA)
      } else if (frac == "Left"){
          tp = tps[2]
          getStage(tp)
      } else if (frac == "Right"){
          tp = tps[4]
          getStage(tp)
      } else if (frac == "Mid"){
           tp = tps[3]
           getStage(tp)
      } else if (frac == "Sides"){
          tp1 = tps[1]+(span/4)
          tp2 = tps[4]+(span/4)
          stg1 <- getStage(tp1)
          stg2 <- getStage(tp2)
          return(paste0(stg1,"-",stg2))
      }
  }


  ncombs <- dim(combs)[2]

  maxValCols <- seq(5, ncombs*6, 6)
  maxTimeCols <- seq(6, ncombs*6, 6)

  aMAFC <- cbind(allDifs[,c(maxValCols, maxTimeCols)])

  timeCols <- (ncombs+1):(ncombs*2)

  tps <- sapply(mybreaks, function(x) timetostage(x))

  for (i in timeCols){
      aMAFC[,i] <- sapply(aMAFC[,i], function(x) fracToStage(x, tps))
  }

  write.csv(allDifs, paste0(outdir, "/areaDiferences_geneLevel.csv"))
  write.csv(aMAFC, paste0(outdir, "/aMAFC_geneLevel.csv"))

## PCAs

  print("Plotting PCA..")
  noNA <- xgene[complete.cases(exprs(xgene))]
  df <- t(exprs(noNA))
  df <- as.data.frame(df)

  pca <- prcomp(df)
  cmp1 <- format(round(summary(pca)$importance[2,1]*100, 2), nsmall = 2)
  cmp2 <- format(round(summary(pca)$importance[2,2]*100, 2), nsmall = 2)

  df_pca <- as.data.frame(pca$x)
  df_pca$Type <- noNA@phenoData@data$type
  df_pca$Time <- noNA@phenoData@data$time

  p <- ggplot(df_pca, aes(x=PC1,y=PC2, col = Type, group = Type))
  p <- p + geom_point(aes(size= Time))
  p <- p + geom_path()
  p <- p + scale_x_continuous(name=paste0("PC1: ", cmp1, "%"))
  p <- p + scale_y_continuous(name=paste0("PC2: ", cmp2, "%"))

  ggsave(p, filename = paste0(figPath, "PCA.png"), device = "png", dpi = "retina")

## Expression Plots

print("Plotting Expression Plots...")
expressionPlot <- function(type){

  ## Set df and lims depending on what we are plotting.
  if (type == "gene"){
    df = xgene
    path = "/Ratio/Gene_Level/"

  } else if (type == "gene_red"){
    df = xgene_red
    path = "/Red_Signal/Gene_Level/"

  } else if (type == "probe"){
    df = xprobe
    path = "/Ratio/Probe_Level/"

  } else if (type == "probe_red"){
    df = xprobe_red
    path = "/Red_Signal/Probe_Level/"

  }

  ## Set ylims
  ylim = c(min(exprs(df), na.rm = T), max(exprs(df), na.rm = T))

  ## Set number of plots
  if (run_all_plots == "yes"){
    nplots = dim(df)[1]
  } else {
    nplots = 20
  }

  ## Main Loop
  for (i in 1:nplots){

    ## Set gene for title or gene and probe for probe-level plots.
    gn <- gsub("[/:;.]", "_", fData(df)$Gene_id[i])
    if (type %in% c("probe", "probe_red")){
      prb <- paste0("_", gsub("[/:;.]", "_" , fData(xprobe)$ProbeName[i]))
    } else {
      prb <- ""
    }
    title <- paste0(gn, prb)

    ## Plot
    graf <- melt(df[i,])
    graf["Type"] <- xgene@phenoData@data$type
    graf["Time"] <- xgene@phenoData@data$time
    p <- ggplot(graf, aes(x = Time, y = value, col = Type, group = Type))
    p <- p + geom_point(aes(color = Type, shape = Type)) + geom_line()
    p <- p + coord_cartesian(ylim = ylim)
    p <- p + ggtitle(title)
    ggsave(p, file=paste0(figPath, path, gn, prb, ".jpeg"),
           device = "jpeg", width = 14, height = 10, units = "cm")

  }

}

expressionPlot("gene")
expressionPlot("gene_red")
expressionPlot("probe")
expressionPlot("probe_red")
