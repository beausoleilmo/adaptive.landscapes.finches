# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created Tuesday, February 12, 2020 
# Why: 
  # Script that will 
  #  - load necessary packages, 
  #  - functions and
  #  - run information session 
# Requires 
  #  - refer to the folder structure and all scripts 
  #  - Note that to run all scripts you'll require 
  #     - ImageMagick and 
  #     - FFMPEG
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


# Install libraries only if needed ----------------------------------------------------------------------------
# Source: https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c("adehabitatHR","cowplot","data.table","faux","fields",
                      "GA","GGally","ggnewscale","ggpubr","ggrepel","ggsignif","ggExtra",
                      "ggspatial","grid","gridExtra","gss","kableExtra",
                      "lubridate","magick","mapview","marked","MASS","metR",
                      "mgcv","mixtools","parallel","plotly","readxl","rgdal",
                      "rgl","rJava","rptR","scales","scatterplot3d","sf",
                      "tibble","tidyverse","vegan","viridis","xlsx")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries ----------------------------------------------------------
library(data.table)   # check this https://mac.r-project.org/openmp/
library(tidyverse)    # %>% 
library(cowplot)      # plot_grid in climate.analysis.R 
library(gridExtra)    # for grid.arrange 
library(lubridate)    # for data management 
library(MASS)         # LDA
library(mgcv)         # gam and bam
library(fields)       # Tps
library(parallel)     # for detectCores()
library(gss)          # gssanova
library(rgl)          # 3D plots
library(vegan)        # rda
library(mixtools)     # mixtools::rmvnorm()
library(kableExtra)   # for report 
library(plotly)       # ggplotly() and plot_ly
library(scales)       # pretty_breaks, alpha
library(GGally)       # for ggpairs
library(GA)           # for persp3D

library(rptR)         # For repeatability
library(tibble)       # add_column

# For maps 
library(mapview)      # Quickly show the maps from SF objects 
library(adehabitatHR) # For kerne UD and MCP functions 
library(rJava)        # Install JDK (Java SE Development Kit; Arm 64 DMG Installer)
library(xlsx)         # Deal with Excel files 
library(readxl)       # Deal with Excel files 
library(sf)           # 
library(rgdal)        # 
library(ggspatial)    # 
library(ggpubr)       # 
library(ggsignif)     # 
library(ggrepel)      # 
library(ggExtra)      # 
library(grid)         # 


# For Gif generation ------------------------------------------------------
library(viridis)      # 
library(magick)       # For making the gif 

# Used for adaptive landscape
library(metR)          # geom_text_contour
library(ggnewscale)    # new_scale_color
library(faux)          # rnorm_multi

# For plot_fitness_no_model.R
library(scatterplot3d) #

# fitness landscape in 3D
library(plot3D) #

# For CMR 
library(marked)        #

# Source scripts ----------------------------------------------------------
source('scripts/functions/0.1_misc.R') # Helper functions 
# Load functions 
source('scripts/functions/plot_gams.R') # retrieve data from the spline. 
source("scripts/functions/cor_fun.R")

# Preparation of variables and data  --------------------------------------
source("scripts/functions/isolate_surface.R") # for isolate.surf (find the maximum in a certain region)

source("scripts/functions/euclidean_dist.R")

# For prospective selection 
source("scripts/functions/intersect.lines.R")


# Colour palette ----------------------------------------------------------
pal  = c(`fortis large` = "black", 
         `fortis small` = "#808080FF", 
         fuliginosa = "orange", 
         magnirostris = "darkgreen", 
         scandens = "pink")

# Fitness landscapes ------------------------------------------------------

## Coefficient of variation ------------------------------------------------
cv.var = function(mean, sd) {
  sd/mean*100
}

## Euclidean distance pairwise ---------------------------------------------
dist.mat.all.sp <- function(data, # Data where the POSITION (phenotypic space) or the average phenotypes or fitness peaks are found.
                            trait1, trait2, # will subset the data for the traits you want to look at. 
                            sp = c("fortis large", "fortis small", "fuliginosa", # species names in order of the data
                                   "magnirostris", "scandens"),
                            round = NULL) {
  # Empty vector 
  comparison.eucl.between.sp = NULL
  # Loop for all the different speices 
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(data)) {
      
      # Get euclidean distance 
      ecl.dat = eucli.d(x1 = data[i,trait1], x2 = data[j,trait1], 
                        y1 = data[i,trait2], y2 = data[j,trait2])
      
      # Get species names
      dt.tmp = cbind(sp.comp1 = as.character(sp[i]),
                     sp.comp2 = as.character(sp[j]), 
                     eucl = ecl.dat$eulc, # Euclidean distance s
                     x = ecl.dat$d.x, # difference on X
                     y = ecl.dat$d.y) |> as.data.frame()
      # Record the data 
      comparison.eucl.between.sp = rbind(comparison.eucl.between.sp, dt.tmp)
    }
  }
  # Make matrix 
  eucl.dist.all.comp = matrix(data = as.numeric(comparison.eucl.between.sp$eucl),
                              nrow = nrow(data),
                              ncol = nrow(data))
  
  # Get names of the matrix
  rownames(eucl.dist.all.comp) <- sp
  colnames(eucl.dist.all.comp) <- sp
  
  # Number of comparisons 
  n.comp = length(eucl.dist.all.comp[lower.tri(eucl.dist.all.comp)])
  
  # Clean the output to get rounded numbers 
  if (!is.null(round)) {
  eucl.dist.all.comp.rd = round(eucl.dist.all.comp, digits = round)
  avg.dist = round(mean(eucl.dist.all.comp[lower.tri(eucl.dist.all.comp)]), digits = round)
  sd.dist = round(sd(eucl.dist.all.comp[lower.tri(eucl.dist.all.comp)]), digits = round)
  ran.dist = round(range(eucl.dist.all.comp[lower.tri(eucl.dist.all.comp)]), digits = round)
  } else {
    eucl.dist.all.comp.rd = eucl.dist.all.comp
    avg.dist = mean(eucl.dist.all.comp[lower.tri(eucl.dist.all.comp)])
    sd.dist = sd(eucl.dist.all.comp[lower.tri(eucl.dist.all.comp)])
    ran.dist = range(eucl.dist.all.comp[lower.tri(eucl.dist.all.comp)])
    
  }
  # Make upper matrix to NA (since it's the same numbers )
  eucl.dist.all.comp.rd[upper.tri(eucl.dist.all.comp.rd)] <- NA
  # Report the values 
  return(list(eucl.dist.all.comp.rd = eucl.dist.all.comp.rd, 
              avg.dist = avg.dist, 
              sd.dist = sd.dist, 
              n.comp = n.comp,
              range.dist = ran.dist))
}


bootstrap.eucl <- function(iter) {
  rec.eucl.perm = NULL
  for (i in 1:iter) {
    bd.rdm = bird.data %>% 
      # Randomize species, trait values 
      mutate(sp2.rdm = sample(sp2)) %>% 
      mutate(mbl.rdm = sample(avg.mbl)) %>% 
      mutate(mbd.rdm = sample(avg.mbd)) %>% 
      # Group by the new randomized variable 
      group_by(sp2.rdm) %>% 
      dplyr::summarize(avg.trait1 = mean(get(traits[1])),
                       avg.trait2 = mean(get(traits[2])),
                       avg.mbl.rdm= mean(mbl.rdm),
                       avg.mbd.rdm= mean(mbd.rdm),)
    
    # Since the group is now sp2.rdm, it's ok to use the "avg.traits" for the permutation 
    ecl.dat.rdm = eucli.d(x1 = bd.rdm$avg.trait1, x2 = all.local.max$x, 
                          y1 = bd.rdm$avg.trait2, y2 = all.local.max$y)
    # ecl.dat$eulc
    rec.eucl.perm = cbind(rec.eucl.perm,ecl.dat.rdm$eulc)
  }
  rec.eucl.perm = as.data.frame(rec.eucl.perm)
  row.names(rec.eucl.perm) = bd.rdm$sp2.rdm
  return(rec.eucl.perm)
}


# Adaptive landscape ------------------------------------------------------
pt.avg <- function(adland.fun.out, add.points = FALSE, 
                   col = "red", bg = "red", pch = 19, cex = 3) {
  avg = apply(adland.fun.out$traits, 2, mean)
  points(avg[1],
         avg[2],
         col=col,
         bg = bg,
         pch=pch, cex = cex)
  if (add.points) {
    points(adland.fun.out$traits[,1],
           adland.fun.out$traits[,2], 
           col = scales::alpha("black",alpha = .5),
           pch=19, cex = .1)
  }
}

# Adaptive landscape function  --------------------------------------------
adapt.land.fun <- function(sp.check = c("fortis large", 
                                        "fortis small", 
                                        "fuliginosa", 
                                        "magnirostris", 
                                        "scandens"),   # Species to include when calculating the landscape. Can be 1 or more. 
                           data = bird.d,              # Data to calculate the phenotypes or use the data for empirical measurement
                           n.simul = 50,               # Number of iterations and individuals to be created
                           n.iter = 50,                # n.iter is the 'resolution' from which we calculate the landscape
                           uppy = NULL, upx = NULL,    # Frame to look into (if you want to extend the predictions)
                           downy = NULL, downx = NULL,    # Allows to fine-tune the predicted landscape for visualisation 
                           empirical = NULL,           # Focal species for which you want to empirically calculate the landscape
                           seed = NULL) {              # Replication of the simulated population 
  
  if (!is.null(seed)) {
    set.seed(seed)  
  }
  
  # Number of iterations and individuals to be created 
  n <- n.simul
  
  # setting colour 
  b.5 = scales::alpha("black",alpha = .5)
  
  # Make a summary of the bird data (used to simulate more points)
  summary.data = data %>% 
    group_by(sp2) %>% 
    filter(sp2 %in% sp.check) %>% 
    summarise(m.bl = mean(MedianBeakLength),
              m.bd = mean(MedianBeakDepth),
              m.bw = mean(MedianBeakWidth),
              sd.bl = sd(MedianBeakLength),
              sd.bd = sd(MedianBeakDepth),
              sd.bw = sd(MedianBeakWidth),
              cor.bl.bd = cor(MedianBeakLength,MedianBeakDepth),
              cor.bl.bw = cor(MedianBeakLength,MedianBeakWidth),
              cor.bd.bw = cor(MedianBeakDepth,MedianBeakWidth))
  
  # extract relevant summary and subset for a certain species 
  mu.sp = summary.data[ , c("m.bl", "m.bd", "m.bw")]
  sd.sp = summary.data[ , c("sd.bl", "sd.bd", "sd.bw" )]
  corr.sp = summary.data[ , c("cor.bl.bd", "cor.bl.bw", "cor.bd.bw")]
  
  # Mean across all columns 
  # apply(summary.data[,-which(names(summary.data)=="sp2")], 2, mean)
  
  # If more than 1 species, take the mean of all 
  if (length(sp.check) > 1) {
    mu.sp = apply(mu.sp, 2, mean)
    sd.sp = apply(sd.sp, 2, mean)
    corr.sp = apply(corr.sp, 2, mean)
  }
  
  # Coefficient of variation. 
  # Schluter (2000) used 5% on page 105 in
  coef.variation = sd.sp/mu.sp*100
  names(coef.variation) <- c("bl","bd","bw")
  round(coef.variation, 2)
  
  # If a species name is given in empirical variable, the species data will be used to calculate the adaptive landscape 
  if (!is.null(empirical)) { 
    traits = data %>% 
      filter(sp2 %in% empirical) %>% 
      dplyr::select(avg.mbl, avg.mbd, avg.mbw)
    # Rename traits 
    names(traits) <- c("Beak length", "Beak depth", "Beak width")
  } else {
    # Simulation of traits for a target species 
    traits = faux::rnorm_multi(n = n, # ??rnorm_multi
                               mu = unlist(mu.sp), # Mean of 3 traits (variables)
                               sd = unlist(sd.sp), # SD of variables 
                               r  = unlist(corr.sp), # correlation between variables 
                               varnames = c("Beak length", "Beak depth", "Beak width"),
                               empirical = TRUE)
  }
  
  round(apply(traits, 2, mean), 2)
  round(apply(traits, 2, sd), 2)
  # coefficient of variation. 
  round(apply(traits, 2, sd)/apply(traits, 2, mean)*100, 2)

  round(cor(traits)[lower.tri(cor(traits))], 2)
  

  # Getting the range of the traits to navigate the space from which the adaptive landscape will be calculated from 
  r.bl = range(traits[,"Beak length"])
  r.bd = range(traits[,"Beak depth"])
  
  # overwriting the default range if you want to force the adaptive landscape to be explored beyond the data 
  if(!is.null(uppy)){  r.bd[2] = max(uppy, r.bd[2])}
  if(!is.null(downy)){ r.bd[1] = min(downy, r.bd[1])}
  if(!is.null(downx)){ r.bl[1] = min(downx, r.bl[1])}
  if(!is.null(upx)){   r.bl[2] = max(upx, r.bl[2])}
  
  # make a sequence from the min and max of each trait 
  # n.iter is the 'resolution' from which we calculate the landscape
  seq.bl = seq(r.bl[1], 
               r.bl[2], 
               length.out = n.iter)
  
  seq.bd = seq(r.bd[1], 
               r.bd[2], 
               length.out = n.iter)
  
  # Blank 2-D object in which the mean at all phenotypic values will be calculated 
  tmp2 = matrix(NA, 
                nrow = length(seq.bl), 
                ncol = length(seq.bd))
  
  # Move the phenotypic distribution at all positions of the phenotypic space in order to find the average 
  # expected fitness for all average phenotypic values
  # Move the phenotypic distribution at a certain Z 
  avg.traits = apply(traits, 2, mean)
  avg.bl = avg.traits[1]
  avg.bd = avg.traits[2]
  
  # Make a matrix of the traits 
  tmp.mat.bl = matrix(traits[,c("Beak length")], nrow = nrow(traits), ncol = length(seq.bl))
  tmp.mat.bd = matrix(traits[,c("Beak depth")], nrow = nrow(traits), ncol = length(seq.bd))
  # Explanation: 
  # Each row is a unique simulated trait and each column is a copy of the traits 
  # Each column will be modified in order to move the phenotypes in phenotypic space 
  
  # For all columns, move the phenotypic distribution according to the sequence specified 
  tmp.mat.moved.bl = sweep(x = tmp.mat.bl, MARGIN = 2, FUN = "-", STATS = (avg.bl-seq.bl))
  # Explanation: 
  # This will move the phenotypes in the X direction only. Will need to do the same thing for the Y axis
  
  # Make an empty matrix that will be filled with the average at each point in the sequence 
  mat.mat = matrix(NA,  nrow = length(seq.bd), ncol = length(seq.bl))
  mat.mat2 = matrix(NA, nrow = length(seq.bd), ncol = length(seq.bl))
  mat.mat.link = matrix(NA, nrow = length(seq.bd), ncol = length(seq.bl))
  # Explanation:
  # I will record the AVERAGE of the expected fitness for the simulated phenotypes that have been moved to a specific location
  
  # Find the inverse link function from the model 
  ginv2 = gam3.p$family$linkinv # check https://stackoverflow.com/questions/40985366/prediction-of-poisson-regression
  
  for (j in 1:length(seq.bd)) {
    # Need to iterate for each y-value
    tmp.mat.moved.bd = sweep(x = tmp.mat.bd, 
                             MARGIN = 2, FUN = "-", 
                             STATS = (avg.bd-seq.bd[j]))
    
    tmp.mat.2d = matrix(NA, nrow = nrow(tmp.mat.moved.bl), ncol = length(seq.bl))
    tmp.mat.2d.link = matrix(NA, nrow = nrow(tmp.mat.moved.bl), ncol = length(seq.bl))
    
    for (i in 1:ncol(tmp.mat.moved.bl)) {
      
      tmp.mat.2d[,i] <- predict(object = gam3.p, 
                                newdata = data.frame(avg.mbl = tmp.mat.moved.bl[,i],
                                                     avg.mbd = tmp.mat.moved.bd[,i]), 
                                type = "response") 
      # I think that the better approach to this is to calculate everything on the link scale, 
      # take the mean and then transfer to the response scale
      tmp.mat.2d.link[,i] <- predict(object = gam3.p, 
                                     newdata = data.frame(avg.mbl = tmp.mat.moved.bl[,i],
                                                          avg.mbd = tmp.mat.moved.bd[,i]), 
                                     type = "link", se.fit = F)
      
    }
    # Calculate the AVERAGE fitness (W^{bar}) for each column on the 
    # Response scale 
    W.bar  = apply(tmp.mat.2d,      2, function(x) mean(x, na.rm = TRUE))
    # Would be better to calculate the geometric mean of the RAW expected fitness... Otherwise, not scaled properly 
    # The link data 
    W.bar.link = apply(tmp.mat.2d.link, 2, function(x) mean(x, na.rm = TRUE))
    
    # Store that in the matrix of average fitness for 
    mat.mat[j,]  <- W.bar
    # mat.mat[j,]  <- ginv2(W.bar.t) # same as below 
    mat.mat2[j,] <- ginv2(W.bar.link)
    # link scale without transformation (for furture manipulation)
    mat.mat.link[j,] <- (W.bar.link)
    print(j)
  }
  
  # Export the data 
  return(list(x = seq.bl, 
              y = seq.bd, 
              z = t(mat.mat), 
              z2 = t(mat.mat2),
              z.link = t(mat.mat.link),
              traits = traits,
              coef.variation = coef.variation,
              sd.sp = sd.sp,
              mu.sp = mu.sp))
}

# function to draw the adaptive landscapes on the main plot (ADDING with add = TRUE)
image.adapt.land <- function(data, 
                             xlim = range(my.persp.data$m1), 
                             ylim = range(my.persp.data$m2), 
                             col = color_plate.dens2, 
                             lwd = 3, 
                             brks, n.contour = 5, pts = FALSE,
                             xadj = 0, # Adjust the linesaround the images 
                             yadj = 0) {
  # Draw the adaptive surface 
  image(data$x,data$y,data$z2, 
        add = TRUE,
        xlim = xlim, ylim = ylim,
        xlab = "", ylab = "",
        main = "",
        zlim = zlim, axes=F,
        col = col, breaks=brks)
  
  # Draw contour on the adaptive surface 
  contour(x = data$x, y = data$y, z = data$z2, 
          labcex = 1,
          drawlabels = TRUE,
          nlevels = n.contour, 
          add = T, 
          col = scales::alpha("black",.9), lwd = 2) 
  # Draw box around the adaptive surface 
  segments(x0 = min(data$x)-xadj, y0 = min(data$y)-yadj, x1 = min(data$x)- xadj, y1 = max(data$y) + yadj, lwd = lwd, lty = 2) # left
  segments(x0 = min(data$x)-xadj, y0 = max(data$y)+yadj, x1 = max(data$x)+ xadj, y1 = max(data$y) + yadj, lwd = lwd, lty = 2) # top
  segments(x0 = max(data$x)+xadj, y0 = max(data$y)+yadj, x1 = max(data$x)+ xadj, y1 = min(data$y) - yadj, lwd = lwd, lty = 2) # right
  segments(x0 = max(data$x)+xadj, y0 = min(data$y)-yadj, x1 = min(data$x)- xadj, y1 = min(data$y) - yadj, lwd = lwd, lty = 2) # Bottom
  
  # Adding the points of the simulated individuals conditionally
  if (pts) {
    points(data$traits[,c("Beak length", 
                          "Beak depth")],col=b.5,pch=19)
  }
}



#computation of the standard error of the mean
# Keep attributes of the original object 
sem <- function(x) {
  att.fct = attributes(x)
  se.tmp=sd(x, na.rm = FALSE)/sqrt(length(na.omit(x)))
  attributes(se.tmp) <- att.fct
  return(se.tmp)
} # Careful, length is counting the number of NA




# Plot fitness no model  --------------------------------------------------


plot.fit.nomodel <- function(x,y,z, 
                             ncol = 100, 
                             res = 0.1, # Increment by which the min/max values will be reached, the smaller, the better the resolution 
                             alpha.v = .8, 
                             seed = 12345, 
                             asp = 0, 
                             points = TRUE,
                             alpha.pt = 0.3, 
                             grid.col = scales::alpha("black",.1),
                             ln = FALSE, # Transform z into ln 
                             trans.fun = NULL, # Further transform z into something you want
                             main = "",
                             xlab = "",
                             ylab = "") {
  set.seed(seed)
  # Transform to ln scale 
  if (ln) {
    z = log(z)
  }
  
  # Take function as input to modify your data 
  if (!is.null(trans.fun)) {
    # Source: https://stackoverflow.com/questions/14046195/r-pass-function-in-as-variable
    transform.dat <- function(f) f(z)
    z = transform.dat(trans.fun)  
  }
  
  
  df.cut = data.frame(x,y,z#,xc,yc
  )
  
  # Empty plot
  plot(x = df.cut$x,
       y = df.cut$y, 
       main = main,
       xlab = xlab,
       ylab = ylab,
       cex = df.cut$z/max(df.cut$z)*2, 
       type = "n", 
       asp = asp)
  
  pal.1 = colorRampPalette(colors = c("bisque1","red"), space="rgb")
  breaks <- seq(min(df.cut$z, na.rm = TRUE), 
                max(df.cut$z, na.rm = TRUE),
                length.out=ncol)
  color.gradient <- function(x, xmin = 0, xmax = 1, colors=c("bisque1","red"), colsteps=100) {
    return( colorRampPalette(colors, space="rgb") (colsteps) [ findInterval(x, seq(xmin ,xmax, length.out=colsteps)) ] )
  }
  
  # Add grid
  i.seq = seq(min(x), max(x), by = res)
  j.seq = seq(min(y), max(y), by = res)
  mid.x = NULL
  mid.y = NULL
  z.val = NULL
  pol.valx = NULL
  pol.valy = NULL
  pol.valc = NULL
  col.vec = NULL
  
  # Add the GRID 
  for (i in 1:(length(i.seq)-1)) {
    for (j in 1:(length(j.seq)-1)) {
      tmpdat = df.cut[df.cut$x >= i.seq[i] & df.cut$x < i.seq[i+1] & df.cut$y >= j.seq[j] & df.cut$y < j.seq[j+1], ]
      
      # Calculate average fitness in the bin 
      tmp.z = mean(tmpdat$z, na.rm = TRUE)

      col.gd = ifelse(is.na(tmp.z), "white", color.gradient(tmp.z, xmin = min(df.cut$z, na.rm = TRUE),xmax = max(df.cut$z, na.rm = TRUE)))
      
      polygon(x = c(i.seq[i], i.seq[i+1], i.seq[i+1], i.seq[i]), 
              y = c(j.seq[j], j.seq[j], j.seq[j+1], j.seq[j+1]), 
              col = scales::alpha(col.gd, alpha.v), # selecting the colour based on the z value 
              border = grid.col)
      # Middle of square 
      mid.x = c(mid.x, c(i.seq[i] + i.seq[i+1])/2)
      mid.y = c(mid.y, c(j.seq[j] + j.seq[j+1])/2)
      z.val = c(z.val, tmp.z)
      pol.valx = c(pol.valx, c(i.seq[i], i.seq[i+1], i.seq[i+1], i.seq[i]))
      pol.valy = c(pol.valy, c(j.seq[j], j.seq[j], j.seq[j+1], j.seq[j+1]))
      pol.valc = c(pol.valc, scales::alpha(pal.1(length(breaks))[ceiling(tmp.z)], alpha.v))
      col.vec = c(col.vec, scales::alpha(pal.1(length(breaks))[ceiling(tmp.z)], alpha.v))
    }
  }
  if (points) {
    points(x = df.cut$x,
           y = df.cut$y, 
           cex = df.cut$z/max(df.cut$z)*2, 
           col = scales::alpha("black",alpha.pt), 
           pch = 19)
  }
  
  return(list(mid.x = mid.x, 
              mid.y = mid.y, 
              z.val = z.val, 
              i.seq = i.seq, 
              j.seq = j.seq,
              pol.valx = pol.valx,
              pol.valy = pol.valy,
              pol.valc = pol.valc,
              col.vec = col.vec,
              breaks = breaks))
}

# Info session ------------------------------------------------------------
# Prints information of the session  
source('scripts/functions/00_info_session.R')
