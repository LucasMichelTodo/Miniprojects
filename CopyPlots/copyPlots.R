###########################################################################
## PROJECT: ACORTES_ARRAY
## DESCRIPTION OF FILE: Copy files from a list to a specified path
## AUTHOR: Evarist Planet
## DATE: Sep 2010
###########################################################################

#Paths
pathGeneList <- 'C:/Arxius Alfred/MicroarrayData/Programes/ProgramesEvarist&Rosell/CopyPlots/geneList.txt'
pathLog <- 'C:/Arxius Alfred/MicroarrayData/Programes/ProgramesEvarist&Rosell/CopyPlots'
#pathOrigin <- 'C:/Arxius Alfred/MicroarraysHeatshock/Analysis/GoodResults/soques_genelevel_estimated'
pathOrigin <- 'C:/Arxius Alfred/MicroarraysHeatshock/Analysis/GoodResults/soques_genelevel_noRatio'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/3D7/Analysis/soques_genelevel_estimated'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/3D7/Analysis/soques_genelevel_noRatio'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/D10/Analysis/soques_genelevel_noRatio'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/D10/Analysis/soques_genelevel_estimated'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/7G8/Analysis/soques_genelevel_estimated'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/7G8/Analysis/soques_genelevel_noRatio'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/HB3/Analysis/soques_genelevel_estimated'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/HB3/Analysis/soques_genelevel_noRatio'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/ParentalLines/Analysis/soques_genelevel_estimated_wo3D7'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/ParentalLines/Analysis/soques_genelevel_noRatio'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/HB3/Analysis/soques_probelevel'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/HB3/Analysis/soques_probelevel_noRatio'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/3D7/Analysis/soques_probelevel'
#pathOrigin <- 'C:/Arxius Alfred/MicroarrayData/Organized_by_strain/3D7/Analysis/soques_probelevel_noRatio'

pathDestiny <- 'C:/Arxius Alfred/MicroarrayData/Programes/ProgramesEvarist&Rosell/CopyPlots/tmp'

#start log file
logFile <- file.path(pathLog,paste('Log ',gsub(":","_",Sys.time()),'.txt',sep=''))
cat('LOG FILE\n',file=logFile)

#check paths
if (!file.exists(pathGeneList)) {cat(paste('Gene list has not been found on the given path: ',pathGeneList,'\n'),file=logFile,append=TRUE) ; stop('')}
if (!file.exists(pathOrigin)) {cat(paste('Given path for plots does not exist: ',pathOrigin,'\n'),file=logFile,append=TRUE) ; stop('')}
if (!file.exists(pathDestiny)) {cat(paste('Path where the plots are suposed to be copied does not exist: ',Destiny,'\n'),file=logFile,append=TRUE) ; stop('')}

#read geneList
val <- try(geneList <- scan(pathGeneList,'character'),silent=TRUE)
if (inherits(val, "try-error")) {
  cat(paste('Error reading geneList. File ',geneList,' not found or is not a plain text file containing one gene per row.\n',sep=''),file=logFile,append=TRUE) ; stop('')
}
geneList <- gsub(':','#',geneList) #colon is not supported under windows

#copy files
numErr <- 0
for (i in 1:length(geneList)) {
  filePath <- file.path(pathOrigin,paste(geneList[i],'.png',sep=''))
  if (file.exists(filePath)) {
    file.copy(filePath,file.path(pathDestiny,paste(geneList[i],'.png',sep='')))
  } else {
    numErr <- numErr +1 
    cat(paste('The following gene has no file to be copied: ',geneList[i],'\n',sep=''),file=logFile,append=TRUE)
  }
}

#Final report
if (numErr==0) {
  cat('Process completed without mistakes. All genes in the gene list have been copied.\nProcess OK',file=logFile,append=TRUE)
} else {
  cat('\n',paste(length(geneList)-numErr,' out of ',length(geneList),' files have been copied. ',numErr,' files that were in the gene list have not been copied.',sep=''),file=logFile,append=TRUE)
}
