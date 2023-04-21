# Description -------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created March 1, 2023 
# Why:
    # Animation that shows how an adaptive landscape is calculated based on the fitness function (landscape)
# Requires 
    # ImageMagick and FFMPEG
    # All required functions are within the script 
# NOTES: 
    # You can play with a mathematical equation to see what you to see 
    # https://www.desmos.com/calculator
    # You can use the image_animate() function to lengthen or shorten the time to see a particular frame 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
library(viridis)
library(magick)

wide = FALSE

# Folder path 
export.path = "output/adaptive_land.animation/_fitplot.gif"
dir.create(export.path, showWarnings = FALSE, recursive = TRUE)
# Name of gif 
gif.name = "animated.adaptive.landscape.example.gif"

# You can print an adaptive landscape for a population with more variance in their phenotypes
if(wide){
  export.path = "output/adaptive_land.animation/_fitplot_wide.gif"
  gif.name = "animated.adaptive.landscape.example_wide.gif"
}

set.seed(1234)
# Parameter for plotting
# Set x variable for the fitness function 
x = seq(-1,7,by = .01)

# Math parameters for Fitness function 
s = .5
m = 1.5
t = 1.7
u = 3.4
div1 = 2.4
div2 = 3.5

# Set population parameters  
n = 1000 # NB individuals, try with 10, and 1000 
# Get fake phenotypes 
pheno = rnorm(n = n, mean = 3.2, sd = 0.5) # Phenotypes for the population
# Make y height just for vizualisation purpose 
fit = rnorm(length(pheno), mean = .12, sd =.01)


# Load functions ----------------------------------------------------------
# Will redraw the plot to add points for gif 
reset.plot = function(mar = c(4,4,.25,.25), col = "red") {
  # Set parameter of plotting area 
  par(mar = mar)
  # Plot fitness function 
  y = fit.function(x = x,s = s, m = m,t = t,u = u, div1 = div1, div2 = div2)
  plot(y~x, 
       xlim = c(-0.5,6),
       ylim = c(0,.8),
       type = "l", lwd = 5, col = col, 
       ylab = "Fitness", xlab = "Phenotypes")
  
}
# Draw arrow to specific location 
add.arrow.text = function(text,x, y, code = 2, x.shift = 0.1, y.shift = 0, col = "red", ...) {
  text(x +x.shift, y = y + y.shift, labels = text, pos = 4, ...)
  arrows(x1 = x[1],    y1 = y[1], 
         x0 = x[1]+x.shift+0.05, y0 = y[1] + y.shift,
         angle = 26, length = .1, code = code, lwd = 3, 
         col = col)
}

# Make a fake population 
add.fake.pop =  function(x = pheno, y = fit,  alp = 1) {
  # Add pop to plot 
  points(x = x, y = y, col = scales::alpha("black", alpha = alp), cex = .2) 
  # Add population mean 
  points(mean(x), y = mean(y), pch = 21, bg = scales::alpha("grey50", .9), col = "black", cex = 2.5)
  # Line at the mean 
  abline(v = mean(x), lty = 3)
  # Calculate the density of the population 
  dens.xy = density(x = x, n = 512, adjust = 2) # Adjust to make it smoother 
  # Add density 
  lines(x = dens.xy$x, y = dens.xy$y/3, lty = 2)
}

# Mathematical fitness function (2 peaks )
# This is a speculative fitness function, you could design your own.
fit.function  <- function(x,s,m,t,u, div1 = 2, div2 = 2.5) (
  exp(-s^(-2)*(x-m)^(2))/div1+exp(-t^(-2)*(x-u)^(2))/div2
)


# Adds a legend to the plot 
add.legend <- function(variables) {
  legend("topright",legend = c("Population", 
                               "Fitness function",
                               "Phenotypic distribution",
                               "Recentered population",
                               "Adaptive landscape"
  ),
  col = c("black",
          "red",
          # viridis(1),
          # viridis(3)[2],
          "black",
          viridis(5, alpha = .5)[2],
          "black"),
  pt.bg = c("black",
            "red",
            # viridis(1),
            viridis(3)[2],
            viridis(5, alpha = .5)[2],
            viridis(3)[3]),
  bg = scales::alpha("white",.5),
  lty = c(NA,1,#1,
          2,NA,NA,NA,NA),
  lwd = c(NA,3,#3,
          2,NA,NA,NA,NA),
  pch = c(21,NA,#NA,
          NA,21,21,21,21))
}

add.pop.moved = function(index, y = fit, col = viridis(5, alpha = .5)[2], seed = 123456) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # Add fake points of population traits
  points(list.fun[[index]]$pop.traits, y = y-hpop,                           pch = 21, col = col,cex = .2) # y pos before rnorm(length(list.fun[[5]]$pop.traits), mean = hpop, sd =.01)
  # Mean of population 
  points(list.fun[[index]]$mean.pheno, y = mean(y)-hpop,                     pch = 21, bg = col, cex = 2.5)
  # Point on the fitness landspcae 
  points(list.fun[[index]]$mean.pheno, y = list.fun[[index]]$adapt.mean.fit, pch = 21, bg = col, cex = 1)
  # Point on adaptive landscape 
  points(list.fun[[index]]$mean.pheno, y = list.fun[[index]]$mean.fit,       pch = 21, bg = col, cex = 1)
}

add.segments.to.fit = function(x, n, alpha.lvl = .1, add.points = FALSE) {
  segments(x0 = x, 
           x1 = x, 
           y0 = 0, 
           y1 = fit.function(x = x, s,m = m,t = t,u = u, div1 = div1, div2 = div2),
           col = scales::alpha("black",alpha.lvl))
  if (add.points) {
    points(x = pheno, 
           y = fit.function(x = pheno, s,m = m,t = t,u = u, div1 = div1, div2 = div2), pch = 19, cex = 1.1)  
  }
}



# Draw plots --------------------------------------------------------------
png(filename = paste0(export.path, "/adapt.land.example%03d.png"),
    width = 7,height = 4.5,units = "in", res = 300, pointsize = 12, bg = "white")

# Plot fitness function 
# Frame 1
reset.plot(); add.legend()
text.explain = "\n\nThis is an example of a fitness landscape (solid red line)."
legend(#"topleft", 
  x = -.75, y = .850,
  legend = paste0(strwrap(text.explain, 35), collapse = "\n"), bty = "n")

# Frame 2
reset.plot(); add.legend()
text.explain = "Imagine we sample individuals (black dots) and measure a phenotype (z) for each (we are going to use the individual fitness landscape to estimate their fitness)."
legend(
  x = -.75, y = .860,
  legend = paste0(strwrap(text.explain, 35), collapse = "\n"), bty = "n")
add.fake.pop(y = rep(0,length(pheno))); add.arrow.text(text = bquote(bar(z)), x = mean(pheno), y = mean(0), y.shift = 0.05, col = "grey70")

# Frame 3
reset.plot(); add.legend()
text.explain = "\n\nWe can get the expected fitness value for all these individuals from the fitness landscape."
# Add segments to show what would be expected fitness 
add.segments.to.fit(x = pheno, n = n, alpha.lvl = .05)
legend(
  x = -.75, y = .850,
  legend = paste0(strwrap(text.explain, 35), collapse = "\n"), bty = "n")
add.fake.pop(y = rep(0,length(pheno))); add.arrow.text(text = bquote(bar(z)), x = mean(pheno), y = mean(0), y.shift = 0.05, col = "grey70")


# Frame 4
reset.plot(); add.legend()
# Add segments to show what would be expected fitness 
add.segments.to.fit(x = pheno, n = n, alpha.lvl = .05)
text.explain = "\n\nLet's move the phenotypes of each individual up so we can see them better."
legend(
  x = -.75, y = .850,
  legend = paste0(strwrap(text.explain, 35), collapse = "\n"), bty = "n")
add.fake.pop(); add.arrow.text(text = bquote(bar(z)), x = mean(pheno), y = mean(fit), y.shift = 0.05, col = "grey70")

# Frame 5
reset.plot(); add.legend()
# Add segments to show what would be expected fitness 
add.segments.to.fit(x = pheno, n = n, alpha.lvl = .05)
text.explain = "\nThe black point (pointed to by the black arrow) is the mean phenotype and mean expected fitness (which is one point of the adaptive landscape)."
legend(
  x = -.75, y = .850,
  legend = paste0(strwrap(text.explain, 35), collapse = "\n"), bty = "n")
add.fake.pop(); add.arrow.text(text = bquote(bar(z)), x = mean(pheno), y = mean(fit), y.shift = 0.05, col = "grey70")

points(mean(pheno), 
       mean(fit.function(x = pheno, s,m = m,t = t,u = u, div1 = div1, div2 = div2)), pch = 19)
# Change text to match symbols
add.arrow.text(text = (c(expression(Mean~of~phenotypes~(bar(z))~and~phantom(0)), 
                         bquote("mean of fitness" ~(bar(W))~"for"),
                         bquote("each individual")
                         )), 
               x = c(mean(pheno),
                     mean(pheno),
                     mean(pheno)+.51), 
               y = c(mean(fit.function(x = pheno, s,m = m,t = t,u = u, div1 = div1, div2 = div2))+0.02,
                     mean(fit.function(x = pheno, s,m = m,t = t,u = u, div1 = div1, div2 = div2))-0.02,
                     mean(fit.function(x = pheno, s,m = m,t = t,u = u, div1 = div1, div2 = div2))-0.055),
               adj=0,
               y.shift = .1,col = "black")

# Frame 6+++
# setting variables to record values for the adaptive landscape 
adapt.line = NULL
adapt.x = NULL
list.fun = list()
# For loop to get the adaptive landscape values 
extend = 2 # How far from the mean phenotype should you go 
for (i in seq(mean(pheno)*extend,-mean(pheno)*extend, by = -0.01)) {
  pheno.moved = pheno - i # Move the population mean (by moving each individual points)
  fit.out = fit.function(x = pheno.moved, s = s, m = m,t = t,u = u, div1 = div1, div2 = div2) # get the fitness of each individual in the simulated population 
  adapt.x = c(adapt.x, mean(pheno.moved)) # record the position of phenotype for the mean 
  adapt.line = c(adapt.line, mean(fit.out)) # record mean fitness of population
  
  # Iterate points 
  adapt.mean.fit = fit.function(x = mean(pheno.moved), s = s, m = m,t = t,u = u, div1 = div1, div2 = div2)
  list.fun[[length(list.fun) + 1]] =  list(mean.fit = mean(fit.out), 
                                           mean.pheno = mean(pheno.moved), 
                                           adapt.mean.fit = adapt.mean.fit,
                                           pop.traits = pheno.moved,
                                           fit.out = fit.out)
}

length(list.fun)
hpop = 0
index1 = 350
index2 = 580
index3 = 850

# Will record the position of each points 
x.adapt.all = NULL
y.adapt.all = NULL

for (j in seq(250, 1000, by = 10)) {
  # Make the same plot as before but change certain plotting information 
  reset.plot();add.fake.pop(alp = .1); add.arrow.text(text = bquote(bar(z)), x = mean(pheno), y = mean(fit), y.shift = 0.05, col = "grey70")
  points(mean(pheno), 
         mean(fit.function(x = pheno, s,m = m,t = t,u = u, div1 = div1, div2 = div2)), pch = 19)
  
  # Add points of the moved population (moved based on the averge of the traits)
  add.pop.moved(index = j, col = viridis(5, alpha = .2)[2])
  # Add segments to show what would be expected fitness 
  add.segments.to.fit(x = list.fun[[j]]$pop.traits, n = n, alpha.lvl = .01)
  
  # Get the x position
  x.adapt = list.fun[[j]]$mean.pheno
  # Get the y position
  y.adapt = list.fun[[j]]$mean.fit
  # Keep track of the positions for each iteration
  x.adapt.all = c(x.adapt.all, x.adapt)
  y.adapt.all = c(y.adapt.all, y.adapt)
  points(x.adapt.all, y.adapt.all, pch = 21, bg = viridis(3)[3], cex = 1)
  
  # Adding the zbar information 
  add.arrow.text(text = bquote(bar(z)[i]), x = list.fun[[j]]$mean.pheno, y = mean(fit)-hpop, x.shift = .25, y.shift = .03, col = "blue4")
  # Adding the adaptive landscape information 
  add.arrow.text(text = bquote(bar(W)[i]), x = x.adapt, y = y.adapt, x.shift = .25, y.shift = .04, col = "goldenrod3")
  # Adding the f(z) landscape information 
  add.arrow.text(text = bquote(f(z)==hat(W)), x = list.fun[[j]]$mean.pheno, y = list.fun[[j]]$adapt.mean.fit, x.shift = .25, y.shift = .09, col = "red")
  
  # Add the same legend each time 
  add.legend()
  
  # Text explanation
  title = list("We 'move' the original population's",
               bquote("mean"~bar(z)~"to a new"~bar(z)[i]~"and calculate the"),
               "average fitness at that new mean",
               bquote("phenotype"~(bar(W)[i])~"of the population to get"),
               "the adaptive landscape.")
  leg = do.call(expression, title)
  text.explain = "We 'move' the original population's mean to a new z bar_i and calculate the average fitness at that new mean phenotype (W bar_i) of the population to get the adaptive landscape. "
  legend(x = -.75, 
         y = .770,
         y.intersp = 1.0,
         legend = leg, 
         bty = "n")
}

reset.plot();add.fake.pop(alp = .1); add.arrow.text(text = bquote(bar(z)), x = mean(pheno), y = mean(fit), y.shift = 0.05, col = "grey70")
points(mean(pheno), 
       mean(fit.function(x = pheno, s,m = m,t = t,u = u, div1 = div1, div2 = div2)), pch = 19)

# Get the x position
x.adapt = list.fun[[j]]$mean.pheno
# Get the y position
y.adapt = list.fun[[j]]$mean.fit
# Keep track of the positions for each iteration
x.adapt.all = c(x.adapt.all, x.adapt)
y.adapt.all = c(y.adapt.all, y.adapt)
points(x.adapt.all, y.adapt.all, pch = 21, bg = viridis(3)[3], cex = 1)

# Adding the zbar information 
add.arrow.text(text = bquote(bar(z)[i]), x = list.fun[[j]]$mean.pheno, y = mean(fit)-hpop, x.shift = .25, y.shift = .03, col = "blue4")
# Adding the adaptive landscape information 
add.arrow.text(text = bquote(bar(W)[i]), x = x.adapt, y = y.adapt, x.shift = .25, y.shift = .04, col = "goldenrod3")
# Adding the f(z) landscape information 
add.arrow.text(text = bquote(f(z)==hat(W)), x = list.fun[[j]]$mean.pheno, y = list.fun[[j]]$adapt.mean.fit, x.shift = .25, y.shift = .09, col = "red")

# Add the same legend each time 
add.legend()

# Text explanation
legend(x = -.75, 
       y = .770,
       y.intersp = 1.0,
       legend = leg, 
       bty = "n")

# When plot is done, add a line for the adaptive landscape 
# Show Adaptive landscape from the population sliding mean
points(adapt.x, adapt.line, type = "l", lwd = 5, col = viridis(3)[2])

# Add arrow of adaptive function 
add.arrow.text(text = bquote(f(bar(z))==bar(W)), x = 2.4, y = .3, x.shift = .25, y.shift = .03, col = viridis(3)[2])


dev.off()


# Make GIF animation ------------------------------------------------------
## list file names and read in
imgs <- list.files(path = export.path,
                   pattern = glob2rx("*.png"), 
                   full.names = TRUE)
nb.frames = length(imgs)

delay.frames = nb.frames
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)
## animate at DELAY X*1/100 seconds OR (FPS: 20 frames per second)
nb.introframes = 5
nb.exitframes = 1
img_animated <- image_animate(img_joined, 
                              # fps = 20,
                              delay = c(rep(500, 1),
                                        rep(700, 4),
                                        rep(20, 76),
                                        300))
## save to disk
image_write(image = img_animated,
            path = paste(export.path, gif.name, sep = "/"))

# Convert GIF to MP4 ------------------------------------------------------
# -y overwrites
system(paste("ffmpeg -y -i", paste(export.path, gif.name, sep = "/"), 
             "-c:v libx264 -pix_fmt yuv420p -movflags +faststart",
             paste("output/adaptive_land.animation", 
                   paste0(sub('\\.gif$', '', gif.name) , 
                          ".mp4"), 
                   sep = "/"), 
             sep = " "))



# Make gif smaller in size 
# Source: https://legacy.imagemagick.org/discourse-server/viewtopic.php?t=29205

cmd.smaller.gif = paste("convert", 
      paste(export.path, gif.name, sep = "/"),
      "-coalesce -scale 700x525 -fuzz 2% +dither -layers Optimize ",
      paste(export.path, 
            paste0(sub('\\.gif$', '', gif.name) , "_small", 
                   ".gif"), 
            sep = "/"))

system(cmd.smaller.gif)
