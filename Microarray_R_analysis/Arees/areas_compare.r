library(tidyverse)
library(sp)

load("/media/lucas/Disc4T/Projects/Anastasia/Array_EpiReset/Arrays_EpiReset_rep1/R_results/geneLevel.RData")

######## Old Areas ########

 # 5.3 AREAS

timeCorrectedEpxrs <- function(exprsx,maxMin.time.est,minMax.time.est,time,soque) {
  selClose <- function(y,x,time.est,pos='before') {
    myTime <- time.est[soque==y]
    if (any(myTime==x)) {
      ans <- rep(FALSE,length(myTime))
    } else {
      if (pos=='previous') ans <- myTime==max(myTime[myTime<x]) else ans <- myTime==min(myTime[myTime>x])
    }
    return(ans)
  }
  getNewVals <- function(x,sel1,sel2) {
    y1 <- exprsx[,sel1]
    y2 <- exprsx[,sel2]
    x1 <- time.est[sel1]
    x2 <- time.est[sel2]  
    a <- matrix((x-x1)/(x2-x1))
    b <- y2-y1
    c <- t(apply(b,1,function(x) x * a))
    ans <- c + y1
    return(ans)
  }
  sel1 <- as.logical(sapply(unique(soque),function(y) selClose(y,maxMin.time.est,time.est,'previous')))
  sel2 <- as.logical(sapply(unique(soque),function(y) selClose(y,maxMin.time.est,time.est,'next')))
  newValsBefore <- getNewVals(maxMin.time.est,sel1,sel2)  
  exprsx[,sel1] <- newValsBefore
  sel1 <- as.logical(sapply(unique(soque),function(y) selClose(y,minMax.time.est,time.est,'previous')))
  sel2 <- as.logical(sapply(unique(soque),function(y) selClose(y,minMax.time.est,time.est,'next')))
  newValsAfter <- getNewVals(minMax.time.est,sel1,sel2)  
  exprsx[,sel2] <- newValsAfter
  return(exprsx)
}

getMaxDif <- function(x) {
  if (sum(!is.na(x))>1) abs(max(x[!is.na(x)])-min(x[!is.na(x)])) else NA
}

compArea <- function(soqueName,soque,exprsx,time.est,from,to) {
  if (missing(from) & missing(to)) {
    sel <- soque == soqueName
    exprsx <- exprsx[,sel]
    x <- time.est[sel]
  } else if (!missing(from) & missing(to)) {
    uniqueTimesInSoque <- as.logical(sapply(unique(soque),function(x) !duplicated(time.est[soque==x]))) #if we have repeated times use only first
    sel <- soque == soqueName & time.est>=from & uniqueTimesInSoque
    prevTime <- c(max(time.est[soque == soqueName & time.est<from]),min(time.est[soque == soqueName & time.est>from]))
    prevSel <- soque == soqueName & time.est %in% prevTime & uniqueTimesInSoque
    firstexprs <- apply(exprsx[,prevSel],1,function(x) ifelse(any(is.na(x)),NA,approxfun(x=prevTime,y=x)(from)))
    exprsx <- cbind(firstexprs,exprsx[,sel])
    x <- c(from,time.est[sel])
  } else if (!missing(from) & !missing(to)) {
    uniqueTimesInSoque <- as.logical(sapply(unique(soque),function(x) !duplicated(time.est[soque==x]))) #if we have repeated times use only first
    sel <- soque == soqueName & (time.est>=from & time.est<=to) & uniqueTimesInSoque
    prevTime <- c(max(time.est[soque == soqueName & time.est<from]),min(time.est[soque == soqueName & time.est>from]))
    prevSel <- soque == soqueName & time.est %in% prevTime & uniqueTimesInSoque
    nextTime <- c(max(time.est[soque == soqueName & time.est<to]),min(time.est[soque == soqueName & time.est>to]))
    nextSel <- soque == soqueName & time.est %in% nextTime & uniqueTimesInSoque
    firstexprs <- apply(exprsx[,prevSel],1,function(x) ifelse(any(is.na(x)),NA,approxfun(x=prevTime,y=x)(from)))
    lastexprs <- apply(exprsx[,nextSel],1,function(x) ifelse(any(is.na(x)),NA,approxfun(x=nextTime,y=x)(to)))
    exprsx <- cbind(firstexprs,exprsx[,sel],lastexprs)
    x <- c(from,time.est[sel],to)
  } else if (missing(from) & !missing(to)) {
    uniqueTimesInSoque <- as.logical(sapply(unique(soque),function(x) !duplicated(time.est[soque==x]))) #if we have repeated times use only first
    sel <- soque == soqueName & time.est<=to & uniqueTimesInSoque
    nextTime <- c(max(time.est[soque == soqueName & time.est<to]),min(time.est[soque == soqueName & time.est>to]))
    nextSel <- soque == soqueName & time.est %in% nextTime & uniqueTimesInSoque
    lastexprs <- apply(exprsx[,nextSel],1,function(x) ifelse(any(is.na(x)),NA,approxfun(x=nextTime,y=x)(to)))        
    exprsx <- cbind(exprsx[,sel],lastexprs)
    x <- c(time.est[sel],to)
  }
  ans <- matrix(NA,nrow(exprsx))
  notNaGenes <- !apply(exprsx,1,function(x) any(is.na(x)))
  y <- exprsx[notNaGenes,]
  idx <- 1:(ncol(exprsx)-1)
  base <- x[idx+1] - x[idx]
  high.square <- sapply(idx,function(x) rowMin(y[,c(x,x+1)]))
  high.triangle <- sapply(idx,function(x) rowMax(y[,c(x,x+1)]))
  square <- rowSums(t(apply(high.square,1,function(x) x * base)))
  triangle <- high.triangle - high.square
  triangle <- rowSums(t(apply(triangle,1,function(x) x * base / 2)))
  area <- square + triangle
  ans[notNaGenes] <- area
  return(ans)
}

getArea <- function(mybreaks,soques,exprsx,time.est,soque) {
  area1 <- sapply(soques,function(x) compArea(x,soque,exprsx,time.est,to=mybreaks[3]))
  area2 <- sapply(soques,function(x) compArea(x,soque,exprsx,time.est,from=mybreaks[2],to=mybreaks[4]))
  area3 <- sapply(soques,function(x) compArea(x,soque,exprsx,time.est,from=mybreaks[3]))
  area4 <- area1 + area3 - area2
  colnames(area1) <- paste('left',colnames(area1),sep='.')
  colnames(area2) <- paste('mid',colnames(area2),sep='.')
  colnames(area3) <- paste('right',colnames(area3),sep='.')
  colnames(area4) <- paste('sides',colnames(area4),sep='.')
  area <- cbind(area1,area2,area3,area4)
  area1.maxDif <- apply(area1,1,getMaxDif)
  area2.maxDif <- apply(area2,1,getMaxDif)
  area3.maxDif <- apply(area3,1,getMaxDif)
  area4.maxDif <- apply(area4,1,getMaxDif)
  area.maxDif <- cbind(area1.maxDif,area2.maxDif,area3.maxDif,area4.maxDif)
  colnames(area.maxDif) <- cbind('area.left','area.mid','area.right','area.sides')
  ans <- list(area=area,area.maxDif=area.maxDif)
  return(ans)
}

getSampledSoque <- function(x,soque,time) {
  ans <- vector('character',length=length(soque))
  for (i in 1:length(unique(time))) {
    if (i==1) {
      ans[time==unique(time)[i]] <- as.character(sample(soque[time==unique(time[i])])) #sample soques (on each time)
    } else {
      tmpsoque <- soque[time==unique(time[i])]
      numSoques <- length(time[time==unique(time)[i]])
      prvTimeLast <- ans[time==unique(time[i-1])][numSoques]
      ans[time==unique(time)[i]][sample(numSoques-1)[1]] <- prvTimeLast
      tmpsoque <- tmpsoque[tmpsoque!=prvTimeLast]
      for (j in 1:numSoques) {
        if (ans[time==unique(time)[i]][j]=='') {
          tmpsoqueIn <- as.character(sample(tmpsoque[tmpsoque!=ans[time==unique(time)[i-1]][j]]))[1]        
          ans[time==unique(time)[i]][j] <- tmpsoqueIn
          tmpsoque <- tmpsoque[tmpsoque!=tmpsoqueIn]
        }
      }
    }
  }
  return(ans)
}

getPval <- function(soqueNames,sampledSoque,exprsx,time,time.est,maxDif,b,mybreaks) {
  myFun <- function(soqueNames,soque,exprsx,time,time.est,mybreaks) {
    myOrder <- order(soque,time.est)
    area.maxDif <- getArea(mybreaks,soqueNames,exprsx[,myOrder],time.est[myOrder],soque[myOrder])[['area.maxDif']]
    maxDif <- apply(area.maxDif,1,function(x) ifelse(any(!is.na(x)),max(x,na.rm=TRUE),NA))
    return(maxDif)
  }
  if (mc.cores==1) {
    maxDif.perm <- lapply(sampledSoque,function(x) myFun(soqueNames,x,exprsx,time,time.est,mybreaks))
  } else {
    maxDif.perm <- mclapply(sampledSoque,function(x) myFun(soqueNames,x,exprsx,time,time.est,mybreaks),mc.preschedule=TRUE,mc.cores=mc.cores)
  }
  maxDif.perm <- do.call(cbind,maxDif.perm)
  maxDif.perm <<- maxDif.perm
  ans <- sapply(1:nrow(maxDif.perm),function(i) sum(maxDif.perm[i,][!is.na(maxDif.perm[i,])] > maxDif[i]) / sum(!is.na(maxDif.perm[i,])) )
  return(ans)
}

#params
set.seed('20101028')
b <- 10000 #number of permutations

#load and order
#load(file.path(dir,'xgene_estimated.RData'))
myOrder <- order(pData(xgene)$type,pData(xgene)$time)
xgene <- xgene[,myOrder]
#estim <- estim[,myOrder]
geneDesc <- fData(xgene)$name

#preprocess
soques <- levels(pData(xgene)$type)
#if (!file.exists(file.path(figuresPath,'soques_genelevel_estimated/'))) dir.create(file.path(figuresPath,'soques_genelevel_estimated/'))
#outDir <- file.path(figuresPath,'soques_genelevel_estimated/')
time.est <- pData(xgene)$time
time <- pData(xgene)$teor_time
soque <- pData(xgene)$type
exprsx <- exprs(xgene)

#set same min and max time
maxMin.time.est <- max(sapply(unique(soque),function(x) min(time.est[soque==x])))
minMax.time.est <- min(sapply(unique(soque),function(x) max(time.est[soque==x])))
exprsx <- timeCorrectedEpxrs(exprsx,maxMin.time.est,minMax.time.est,time,soque)
sel.previous2first <- as.logical(sapply(unique(soque),function(x) time.est[soque==x]==max(time.est[soque==x][time.est[soque==x]<=maxMin.time.est])))
sel.next2last <- as.logical(sapply(unique(soque),function(x) time.est[soque==x]==min(time.est[soque==x][time.est[soque==x]>=minMax.time.est])))
pData(xgene)$time[sel.previous2first] <- maxMin.time.est
pData(xgene)$time[sel.next2last] <- minMax.time.est
time.est <- pData(xgene)$time
#solve issue where two points are before first or after last
myFun <- function(x,gb) {
  mytime <- time.est[soque==x]
  afterLast <- mytime[mytime>=minMax.time.est]
  if (length(afterLast)>1) {
    if (gb=='good') ans <- mytime==minMax.time.est else ans <- mytime>minMax.time.est
  } else {
    ans <- rep(FALSE,length(mytime))
  }
  return(ans)
}
good <- as.logical(sapply(unique(soque),function(x) myFun(x,'good')))
wrong <- as.logical(sapply(unique(soque),function(x) myFun(x,'wrong')))
exprsx[,wrong] <- exprsx[,good]
time.est[wrong] <- time.est[good]

#set the minimum value of each gene equal to 0
notAllNa <- apply(exprsx,1,function(x) any(!is.na(x)))
exprsx[notAllNa,] <- t(apply(exprsx[notAllNa,],1,function(x) x - min(x,na.rm=TRUE))) #

#compute areas
mybreaks <- seq(maxMin.time.est,minMax.time.est,length.out=5)
tmp <- getArea(mybreaks,soques,exprsx,time.est,soque)
oldArea <- tmp[['area']]; area.maxDif <- tmp[['area.maxDif']]

######## New Areas ##########
load("/media/lucas/Disc4T/Projects/Anastasia/Array_EpiReset/Arrays_EpiReset_rep1/R_results/geneLevel.RData")

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
computeArea <- function(eset, types, type, mybreaks){
  
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
  
  maxminTP <- mybreaks[1]
  tp1 <- mybreaks[2]
  tp2 <- mybreaks[3]
  tp3 <- mybreaks[4]
  minmaxTP <- mybreaks[5]
  
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

start_time <- Sys.time()
areasDF <- computeArea(xgene)
end_time <- Sys.time()
end_time - start_time


newArea <- cbind(fData(xgene), areasDF)

######## New Areas 2 ########

## Necessary variables for the functions.

load("/media/lucas/Disc4T/Projects/Anastasia/Array_EpiReset/Arrays_EpiReset_rep1/R_results/geneLevel.RData")

eset <- xgene

types <- unique(phenoData(eset)$type)
type <- phenoData(eset)$type
times <- phenoData(eset)$time

maxminTP <- max(sapply(types,function(x) min(times[type==x])))
minmaxTP <- min(sapply(types,function(x) max(times[type==x])))
mybreaks <- seq(maxminTP, minmaxTP, length.out=5)

subtimes <- lapply(types, function(x) unique(sort(c(times[type==x], mybreaks))))
subtimes <- lapply(subtimes, function(x) x[x >= mybreaks[1] & x <= mybreaks[5]])
names(subtimes) <- types


## Impute all necessary points. Subset must be a subset of an eset
## where only samples of one of the conditions are present (i.e. controls)

imputePoints <- function(tp, subset){
  xs <- pData(subset)$time
  if (tp %in% xs){
    return(exprs(subset[,xs == tp]))
  } else {
    x1 <- max(xs[xs < tp])
    x2 <- min(xs[xs > tp])
    y1 <- exprs(subset[,xs == x1])
    y2 <- exprs(subset[,xs == x2])
    
    ms <- (y2-y1)/(x2-x1)
    bs <- y1 - (ms*x1)
    yimp <- (ms*tp)+bs
    
    return(yimp)
  }
}

imputeSubSet <- function(type){
  subset <- eset[, pData(eset)$type == type]
  cols <- lapply(subtimes[[type]], function(x) imputePoints(x, subset))
  impDF <- do.call(cbind, cols)
}

getArea <- function(type){
  
  df <- impDFs[[type]]
  
  allAreas <- c()
  for (i in 1:dim(eset)[1]){
    
    #type <- "treat"
    #i <- 1
    #df <- impDFs[[type]]
    
    min <- mins[i]
    points <- as.data.frame(cbind(subtimes[[type]], df[i,]))
    points[,2] <- points[,2]-min
    colnames(points) <- c("x", "y")
    
    left <- rbind(points[which(points$x <= mybreaks[3]),],
                  c(mybreaks[3], 0),
                  c(mybreaks[1], 0),
                  points[1,])
    
    right <- rbind(points[which(points$x >= mybreaks[3]),],
                   c(mybreaks[5], 0),
                   c(mybreaks[3], 0),
                   points[which(points$x == mybreaks[3]),])
    
    mid <- rbind(points[which(points$x >= mybreaks[2] & points$x <= mybreaks[4]),],
                 c(mybreaks[4], 0),
                 c(mybreaks[2], 0),
                 points[which(points$x == mybreaks[2]),])
    
    sides <- rbind(points[which(points$x <= mybreaks[2]),],
                   c(mybreaks[2], 0),
                   c(mybreaks[4], 0),
                   points[which(points$x >= mybreaks[4]),],
                   c(mybreaks[5], 0),
                   c(mybreaks[1], 0),
                   points[1,])
    
    pols <- list(left, right, mid, sides)
    
    calcArea <- function(x) {ifelse(any(is.na(x)), NA, Polygon(x)@area)}
    areas <- unlist(lapply(pols, function(x) calcArea(x)))
    allAreas <- c(allAreas, list(areas))
  }
  areaDF <- do.call(rbind, allAreas)
  
  titles <- c("Left", "Right", "Middle", "Sides")
  names <- sapply(titles, function(x) paste(type, x, sep = "_"))
  colnames(areaDF)  <- names
  rownames(areaDF) <- rownames(exprs(eset))
  return(areaDF)
}

start_time <- Sys.time()

impDFs <- lapply(types, function(x) imputeSubSet(x))
names(impDFs) <- types

allImputed <- do.call(cbind, impDFs)
mins <- apply(allImputed, 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm=T)))

areas <- lapply(types, function(x) getArea(x))
names(areas) <- types

allAreas <- do.call(cbind, areas)

end_time <- Sys.time()
end_time - start_time



#write.csv(cbind(fData(xgene), allAreas),
#          file = paste0(wd, 'R_results/areas_geneLevel.csv'),
#          row.names=F)

sets <- 1:length(areas)
combs <- combn(sets, 2)
span <- mybreaks[3] - mybreaks[1]

areaDifs <- c()
for (i in 1:dim(combs)[2]){
  
  one <- combs[1,i]
  two <- combs[2,i]
  
  dif1 <- (areas[[one]] - areas[[two]])/span
  prefix  <- paste0(names(areas)[[one]],
                    "_minus_",
                    names(areas)[[two]])
  
  names <- sapply(titles, function(x) paste(prefix, x, sep = "_"))
  colnames(dif1) <- names
  
  dif2 <- -dif1
  prefix  <- paste0(names(areas)[[two]],
                    "_minus_",
                    names(areas)[[one]])
  
  names <- sapply(titles, function(x) paste(prefix, x, sep = "_"))
  colnames(dif2) <- names
  
  areaDifs <- c(areaDifs, list(dif1), list(dif2))
}

allDifs <- do.call(cbind, areaDifs)

#write.csv(cbind(fData(xgene), allDifs),
#          file = paste0(wd, 'R_results/areas_geneLevel_Difs.csv'),
#          row.names=F)

#### Comparisons #######

newOrder <- c(1,7,4,10)
old <- oldArea[, c(newOrder, newOrder+2, newOrder+1)]
colnames(old) <- colnames(allAreas)
rownames(old) <- rownames(allAreas)
new1 <- newArea[,c(-1,-2,-3,-4)]

head(allAreas)
head(old)
head(new1)

table(allAreas == old)
table(new1 == old)

new1[new1 != old,][1:10,]
old[rownames(old) == "1396.pre-tRNA-Pro-1",]
