#### Standarize whole experiment, plots keep shape (not if we stz columns) ####

scaled.dat <- (exp-mean(exp))/sd(exp)
colMeans(scaled.dat, na.rm = T)
apply(scaled.dat, 2, function(x) sd(x, na.rm = T))

figPath <- "/media/lucas/Disc4T/Projects/Miniprojects/Arees/"

for (i in 1:dim(exprs(toyset))[2]){
  graf <- melt(exprs(toyset)[i,])
  graf["Type"] <- toyset@phenoData@data$type
  graf["Time"] <- toyset@phenoData@data$time
  p <- ggplot(graf, aes(x = Time, y = value, col = Type))
  p <- p + geom_point(aes(color = Type, shape = Type)) + geom_line()
  #p <- p + coord_cartesian(ylim = c(minRatioG, maxRatioG))
  p <- p + ggtitle(toyset@featureData@data$Gene_id[i])
  ggsave(p, file=paste0(figPath, gsub("/", ":", toyset@featureData@data$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
}

for (i in 1:dim(exprs(toyset))[2]){
  graf <- melt(scaled.dat[i,])
  graf["Type"] <- toyset@phenoData@data$type
  graf["Time"] <- toyset@phenoData@data$time
  p <- ggplot(graf, aes(x = Time, y = value, col = Type))
  p <- p + geom_point(aes(color = Type, shape = Type)) + geom_line()
  #p <- p + coord_cartesian(ylim = c(minRatioG, maxRatioG))
  p <- p + ggtitle(toyset@featureData@data$Gene_id[i])
  ggsave(p, file=paste0(figPath, gsub("/", ":", toyset@featureData@data$Gene_id[i]), "_Z.jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
}

mean(scaled.dat)
sd(scaled.dat)
