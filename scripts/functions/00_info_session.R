# Description  ------------------------------------------------------------
#### ### ### ## #### ### ### ## #### ### ### ## 
# Info session
# Created by Marc-Olivier Beausoleil
# Thursday, May 12, 2022
# Why:
  # Prints R version and Collect Information About the Current R Session to a file 
# Requires:
# NOTES:
#### ### ### ## #### ### ### ## #### ### ### ## 

info.session <- function(path = "output/session_info", 
                         file = "session_information.txt") {
  # Create session info folder 
  dir.create(path, showWarnings = FALSE)
  
  # Will add session info to this file 
  sink(paste(path, file, sep="/"),append = FALSE)
  
  # sink("~/Desktop/session_information.txt",append = FALSE)
  cat("##### R Version Information ############################################################\n")
  print(version)
  
  cat("\n\n##### Collect Information About the Current R Session ############################################################\n\n")
  
  # Reorder session information 
  # Source: https://stackoverflow.com/questions/41794350/how-do-i-sort-the-output-of-sessioninfo-in-r
  si=sessionInfo()
  si[] <- lapply(si, function(x) if (is.list(x)) x[sort(names(x))] else sort(x))
  si
  print(si)
  
  cat("\n\n##### Other programs ############################################################\n\n")
  
  cat("\n\n##### ImageMagick ##############################\n\n")
  im.v = system("magick -version", intern = TRUE)
  sapply(X = im.v, FUN = function(x) cat(x, '\n'))
  
  cat("\n\n##### FFMPEG ##############################\n\n")
  ff.v = system("ffmpeg -version", intern = TRUE)
  sapply(X = ff.v, FUN = function(x) cat(x, '\n'))
  
  sink()
  
}

info.session(path = "output/session_info", 
             file = "session_information.txt")

# cat("Done!\n")
