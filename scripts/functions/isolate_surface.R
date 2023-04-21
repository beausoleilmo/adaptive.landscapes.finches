# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created May 12, 2022
# Why: 
  # Isolate par of a 3D surface and find the maximum in the region
# Requires 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

isolate.surf = function(dat, #list from the plot.gam.cust() custom function 
                        mod = NULL,
                        # mod.var = c("x1","x2"),
                        mod.var = c("MedianBeakLength","MedianBeakDepth"),
                        variables = c("m1","m2","z"),
                        locator = FALSE,
                        xmin = NULL, # Length
                        xmax = NULL, 
                        ymin = NULL, # Depth
                        ymax = NULL, 
                        points.id = NULL, # Add x, y (for position), and species (for colour) the help find the peaks
                        plot=FALSE) {
  
  # Extract the columns from the list 
  get.df.from.list = data.frame(x = dat[[variables[1]]], y = dat[[variables[2]]])
  get.z = dat[[variables[3]]]
  
  # Instead of providing the values, you can locate them by clicking on the graph 
  if(locator){
    cat("Please, provide 2 points starting from the bottom left corner and top right:")
    # var = readline()
    # var = as.integer(var)
    var = 2
    par(mfrow = c(1,1))
    image(x = get.df.from.list$x, 
          y = get.df.from.list$y, 
          z = get.z)
    
    if (!is.null(points.id)) {
      points(points.id[,1], points.id[,2], 
             col = scales::alpha(pal[points.id[,3]], .4), pch = 19)
    }
    
    my.pts = locator(n = var)
    xmin = my.pts$x[1]
    xmax = my.pts$x[2]
    ymin = my.pts$y[1]
    ymax = my.pts$y[2]
  }
  
  # Isolate the location with the user specified x- and y-axis limits
  find.x = which(get.df.from.list$x>xmin & get.df.from.list$x<xmax)
  find.y = which(get.df.from.list$y>ymin & get.df.from.list$y<ymax)
  
  # Subset the vector for each 
  x = get.df.from.list$x[find.x]
  y = get.df.from.list$y[find.y]
  z = get.z[find.x,find.y] # rows are y, column are x in this type of matrix
  
  # Find the maximum value in the subset 
  max.z = max(z, na.rm = TRUE)
  # which(z == max.z)
  
  # Find the position of the maximum Z value 
  loc.z = which(z == max.z, arr.ind = TRUE)
  # Find the X and Y coordinate in the subsetted dataset 
  x.y.z = c(x[loc.z[1]], y[loc.z[2]])
  
  # More robust: find the location of the max z from the ORIGINAL indicies
  get.df.from.list$x[find.x[loc.z[1]]]
  get.df.from.list$y[find.y[loc.z[2]]]
  get.z[find.x[loc.z[1]],find.y[loc.z[2]]]
  # pt.max.all = which(get.z == max.z, arr.ind = TRUE)
  
  # Get the coordinate from the original data 
  find.local.max = data.frame(
    x = get.df.from.list$x[find.x[loc.z[1]]],
    y = get.df.from.list$y[find.y[loc.z[2]]],
    z = get.z[find.x[loc.z[1]],find.y[loc.z[2]]])
  
  # you can Predict the Z value based on the phenotypes and the model 
  if(!is.null(mod)){
    tmp.df = data.frame(x = x.y.z[1], y = x.y.z[2])
    colnames(tmp.df) <- mod.var
    mod.pred.z = predict(mod, newdata = tmp.df, type = "response")
  } else {mod.pred.z = NULL}
  
  
  # Make a visual verification 
  if(plot){
    par(mfrow = c(1,2))
    image(x = x, y = y, z = z
          # ,xlim = range(get.df.from.list$x),
          # ylim = range(get.df.from.list$y)
          )
    points(x.y.z[1],x.y.z[2], cex = 3,pch = 19)
    image(x = get.df.from.list$x, y = get.df.from.list$y, z = get.z)
    points(x.y.z[1],x.y.z[2], cex = 1,pch = 19)
  }
  
  # Get the data out of the function 
  return(list(x = x, # This is the data from the subset
              y = y, # This is the data from the subset
              z = z, # This is the data from the subset
              max.z = max.z, # Maximum value of Z found
              loc.z = loc.z, # Location in the subset 
              x.y.z = x.y.z, # X-Y coordinate
              local.max = find.local.max,
              mod.pred.z = mod.pred.z) # X-Y-Z coordinate in a dataframe 
         )
}
