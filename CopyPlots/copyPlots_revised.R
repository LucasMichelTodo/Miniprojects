###########################################################################
## PROJECT: ACORTES_ARRAY
## DESCRIPTION OF FILE: Copy files from a list to a specified path
## AUTHOR: Evarist Planet
## DATE: Sep 2010
###########################################################################


##### PATHS: Modify this bit to fit your personal folders #################

# Path to GeneList file
pathGeneList <- '/media/lucas/Disc4T/Projects/Miniprojects/CopyPlots/geneList.txt'

# Path to Log file
pathLog <- '/media/lucas/Disc4T/Projects/Miniprojects/CopyPlots/'

# Path to the folder where the plots are
pathOrigin <- '/home/lucas/ISGlobal/Arrays/Eli_Arrays/Plots/Imputed/10E'

# Path were you wish to move only the selected plots
pathDestiny <- '/media/lucas/Disc4T/Projects/Miniprojects/CopyPlots/Selected_Plots/'

##########################################################################

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
geneList <- gsub('[.:/]','_',geneList) #colon is not supported under windows

# Generate a list of plot-files
plots <- list.files(pathOrigin)

# #copy files
# numErr <- 0
# for (i in 1:length(geneList)) {
#   filePath <- file.path(pathOrigin,paste(geneList[i],'.png',sep=''))
#   if (file.exists(filePath)) {
#     file.copy(filePath,file.path(pathDestiny,paste(geneList[i],'.png',sep='')))
#   } else {
#     numErr <- numErr +1 
#     cat(paste('The following gene has no file to be copied: ',geneList[i],'\n',sep=''),file=logFile,append=TRUE)
#   }
# }

#copy files
numErr <- 0
for (i in 1:length(geneList)) {
  hit_list <- plots[grepl(geneList[i], plots)]
  for (x in hit_list) {
    filePath <- paste(pathOrigin, x, sep = "/")
    file.copy(filePath, paste(pathDestiny, x, sep = "/"))
  } 
  if (length(hit_list) < 1) {
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

