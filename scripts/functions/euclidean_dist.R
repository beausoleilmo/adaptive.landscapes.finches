# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Euclidean distance
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created April 13, 2022
# Why: 
  # calculates the Euclidean distance for 2 points in a 2D graph 
# Requires 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


# Converters --------------------------------------------------------------
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}


# Euclidean distance
eucli.d <- function(x1,x2,y1,y2, plot=FALSE) {
  
  # Difference in the coordinate
  diff.x = x2-x1
  diff.y = y2-y1
  # Calculate the Euclidean distance
  ed = sqrt((diff.x)^2+(diff.y)^2)
  
  # Calculate the angle if there is only ONE input value 
  if (length(diff.y) == 1) {
    # https://math.stackexchange.com/questions/1910825/how-do-i-find-the-angle-a-vector-makes-to-the-x-axis
    if (sign(diff.y) < 0) { # If "diff.y" is negative, will need to take the angle value and subtract it from 360 
      angle <- 360 - (acos(diff.x/ed) * 180) / (pi)#(theta * 180) / (pi) + add.180
    } else {
      angle <- (acos(diff.x/ed) * 180) / (pi)
    } # end else
  } else {
    # Calculate the angle if there is MORE THAN ONE input value 
    tmp.v = numeric(length(diff.y)) # create empty vector where the data will be recorded 
    for (i in 1:length(diff.y)) { # loop for all the values 
      if (sign(diff.y[i]) < 0) {
        tmp.v[i] <- 360 - (acos(diff.x[i]/ed[i]) * 180) / (pi)
      } else {
        tmp.v[i] <- (acos(diff.x[i]/ed[i]) * 180) / (pi)
      } # end else
    } # end for
    angle = tmp.v
  }
  
  if (plot) {
    plot(data.frame(x = c(x1,x2),y = c(y1,y2)), 
         col = c("red","blue"), cex =c(2,1),
         pch = 19, asp = 1, 
         xlim = range(x1,x2), 
         ylim = range(y1,y2))
    segments(x1,y1,x2,y2)
    segments(x1,y1,x2,y1)
    segments(x2,y1,x2,y2)
    
  }
  
  # Return the values calculated 
  return(list(d.x = diff.x,
              d.y = diff.y,
              eulc = ed, 
              angle = angle))
}

# Example: 3-4-5 triangle 
# eucli.d(x1 = 0, x2 = 3, # euclidean
#         y1 = 0, y2 = 4, plot = TRUE)
# eucli.d(x1 = 3, x2 = 0, # euclidean
#         y1 = 4, y2 = 0, plot = TRUE)
