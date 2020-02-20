#### Fake Data (microarrays) ####

exprsx <- rbind(c(1,0,1,0, 0.25,0.75,1.5,2), c(1,1,1,1, 1,1,1,1))
rownames(exprsx) <- c("gene1", "gene2")
colnames(exprsx) <- c("Ctl_1", "Ctl_2", "Ctl_3", "Ctl_4",
                      "Treat_1", "Treat_2", "Treat_3", "Treat_4")

pdata <- data.frame(type=c(rep("Ctl",4), rep("Treat",4)), time=(c(3,5,7,9, 2,4,7,9)))
rownames(pdata) <- colnames(exprsx)
pdata <- new("AnnotatedDataFrame", pdata)
fdata <- data.frame(anot=c("Fake gene 1", "Fake gene 2"))
rownames(fdata) <- rownames(exprsx)
fdata <- new("AnnotatedDataFrame", fdata)
xtoy <- new("ExpressionSet", exprs=exprsx, featureData=fdata, phenoData=pdata)

#### Fake Data 2 (phd project) ####

exprsx <- rbind(c(-5,-3,0,3,3,3,3,3, 5,5,5,5,5,5,5), c(7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7))
rownames(exprsx) <- c("gene1", "gene2")
colnames(exprsx) <- c("Ctl_1", "Ctl_2", "Ctl_3", "Ctl_4", "Ctl_5", "Ctl_6", "Ctl_7", "Ctl_8",
                      "Treat_1", "Treat_2", "Treat_3", "Treat_4", "Treat_5", "Treat_6", "Treat_7" )

pdata <- data.frame(type=c(rep("Ctl",8), rep("Treat",7)), time=(c(0,12,24,32,36,40,44,48,
                                                                  7,15,28,34,41,44,47)))
rownames(pdata) <- colnames(exprsx)
pdata <- new("AnnotatedDataFrame", pdata)
fdata <- data.frame(anot=c("Fake gene 1", "Fake gene 2"))
rownames(fdata) <- rownames(exprsx)
fdata <- new("AnnotatedDataFrame", fdata)
xtoy <- new("ExpressionSet", exprs=exprsx, featureData=fdata, phenoData=pdata)
