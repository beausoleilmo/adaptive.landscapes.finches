# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created November 29, 2022
# Why: 
  # Draw perpendicular lines in a plot 
# Requires 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# NOTE: it only works when the plot is in asp = 1 
intersect.lines <- function(p, # position of the intersecting lines 
                            xlength, ylength, # Length of each line 
                            angl = 0, # if 0, a cross will be drawn 
                            lty = 2,# Type of line
                            col.line  = "black", 
                            text1 = c(0,0), text2 = c(0,0), # Adjust the text position 
                            txt1 = "Δβ (Beak size)", # labels of the text 
                            txt2 = "Δβ (Beak shape)") {
  
  ang2rad = function(angl) {
    rad = angl * pi/180 
    return(rad)
  }
  
  x1x1 = p[1] + cos(ang2rad(angl = angl)) * xlength/2
  y1y1 = p[2] + sin(ang2rad(angl = angl)) * xlength/2
  x1x2 = p[1] - cos(ang2rad(angl = angl)) * xlength/2
  y1y2 = p[2] - sin(ang2rad(angl = angl)) * xlength/2
  
  x2x1 = p[1] - cos(ang2rad(angl = (angl + 90))) * ylength/2
  y2y1 = p[2] - sin(ang2rad(angl = (angl + 90))) * ylength/2
  x2x2 = p[1] + cos(ang2rad(angl = (angl + 90))) * ylength/2
  y2y2 = p[2] + sin(ang2rad(angl = (angl + 90))) * ylength/2
  
  
  segments(x1x1, y1y1,
           x1x2, y1y2, lty = lty, col = col.line)
  segments(x2x1, y2y1,
           x2x2, y2y2, lty = lty, col = col.line)
  
  if (!is.null(text1)) {
    text(x = text1[1]+p[1], y = text1[2]+p[2], labels = txt1, srt=angl, cex = .7)
    text(x = text2[1]+p[1], y = text2[2]+p[2], labels = txt2, srt=angl-90, cex = .7)
  }
}

# Example
# var = c(8:18)
# plot(x = var, y = var,
#      pch = c(21,22,23,24,25),
#      # log = "xy",
#      bg = "black",
#      col = "white", asp = 1,
#      cex = 1);abline(h =0, v =0 , lty = 3)
# 
# intersect.lines(p = (c(10,10)), 
#                 xlength = 2, ylength = 3, angl = 10, lty = 2, 
#                 col.line  = "black", text1 = c(0.2, 0.8), text2 = c(-0.9, 0.2))
