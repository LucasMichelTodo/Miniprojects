library(sp)
library(Biobase)
library(reshape2)
library(ggfortify)
library(tidyverse)

#### Load Data ####
load('~/ISGlobal/Arrays/AnalisisR/Original_Results/xgene_estimated.RData')
colnames(pData(xgene)) <- c("type", "time")

## Functions

## Imputing function

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

df <- computeArea(xgene)
head(df)
write.csv(df, "/media/lucas/Disc4T/Projects/Miniprojects/Microarray_R_analysis/Arees/newMethod_NAs.csv")

