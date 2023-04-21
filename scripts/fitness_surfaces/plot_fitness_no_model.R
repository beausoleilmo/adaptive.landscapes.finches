# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Raw Darwin's finches' fitness landscapes without model
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created October 18, 2022
# Why: 
  # Mean fitness for a phenotypic area without model 
# Requires 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

plot.extact.pheno.fit.no.land = FALSE

# Source scripts ----------------------------------------------------------
source('scripts/0.0_initialize.R')

# read data ---------------------------------------------------------------
load('data/bird.data.RData', verbose = TRUE)
text.size = 16

cat(paste("Plotting to file", 
          paste("output/images/landscape_plots/fitness_lands.nomodel",ext.file,".png", sep = "")), 
    fill = TRUE)

# PLOT fintess without model so you can see a raw fitness landscape --------
# Needed to get the data for ggplot 
gridxy = plot.fit.nomodel(x = bird.data$avg.mbl,
                          y = bird.data$avg.mbd, 
                          z = bird.data$mxcpois+1, 
                          main = "Fitness average",
                          xlab = "Beak length (mm)",
                          ylab = "Beak depth (mm)",
                          ncol = 8, 
                          res = .20, # Size of grid for ggplot (reusing the same data)
                          alpha.v = 1, asp = 1,
                          points = FALSE, alpha.pt = .1, 
                          ln = F, # transform to ln
                          grid.col = scales::alpha("black",.05))

bird.data$mxcpois.cex = c(bird.data$mxcpois+1)/max(c(bird.data$mxcpois+1))*2

# Polygon data. Each polygon is defined as a quadruplet of points (row 1:4, 5:8, etc.)
gd.pt = data.frame(i = gridxy$pol.valx, 
                   j = gridxy$pol.valy,
                   c = rep(gridxy$pol.valc,each =4),
                   z = rep(gridxy$z.val,each =4),
                   col = rep(gridxy$col.vec,each =4),
                   ids = factor(rep(1:length(gridxy$pol.valc),each =4)))

# Add grid to the plot 
new.gridx = data.frame(x = seq(5,16, by = 2.5))
new.gridy = data.frame(y = seq(6,18, by = 2*2))

# GGplot Fitness landscape no model ----------------------------------------------------------
gd.pt$z2 = log(x = gd.pt$z,base = exp(1))
ggp.fit.no.model = ggplot(data = gd.pt, 
                          mapping = aes(i, j)) +
  geom_polygon(aes(fill = (z2), 
                   group = ids), 
               color = "grey", # Colour of grid 
               linewidth = 0.07) + 
  scale_fill_gradient(low = "bisque1",  # add colour to
                      high = "red",
                      na.value = "#FFFFFF",
                      name = "Ln mean\napparent\nsurvival\n(years)") +
  coord_cartesian(xlim = c(6.96, 16.33), # Make the coordinate system equivalent to the other plot 
                  ylim = c(5.02, 18.67)) +
  # add custom grid to plot 
  geom_hline(data = new.gridy, aes(yintercept = y), colour = 'grey50', linetype = "dotted", linewidth = .25) +
  geom_vline(data = new.gridx, aes(xintercept = x), colour = 'grey50', linetype = "dotted", linewidth = .25) +
  # Theme stuff
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "grey100", linetype = "dotted", linewidth = .25),
        panel.grid.minor = element_line(colour = "grey100", linetype = "dotted", linewidth = .25),
        axis.title = element_text(size = text.size),
        axis.text = element_text(size = text.size, colour = "black", vjust = 0.25), 
        axis.text.x = element_text(size = text.size), 
        axis.text.y = element_text(size = text.size),
        plot.title = element_text(size = text.size),
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size),
        legend.key = element_rect(fill = "white", color = "grey80"), 
        plot.tag =element_text(face = "bold", size = 18),
        strip.background =element_rect(fill="white", colour = "white"),
        strip.text = element_text(colour = 'black', size = 13),
        legend.background = element_rect(fill = "white", color = NA)) +
  labs(tag = "A",
       x = "Beak length (mm)",
       y = "Beak depth (mm)", 
       size = 12, alpha = 1); ggp.fit.no.model

cat(paste("Saving ggplot to file", paste("output/data.out/ggplot_data/ggp.fit.no.model.RDSgp.fit.no.model.RDS")), 
    fill = TRUE)

saveRDS(object = ggp.fit.no.model, 
        file = "output/data.out/ggplot_data/ggp.fit.no.model.RDS")