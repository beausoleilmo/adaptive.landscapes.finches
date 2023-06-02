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
        # legend.position = "bottom",
        plot.tag =element_text(face = "bold", size = 18),
        strip.background =element_rect(fill="white", colour = "white"),
        strip.text = element_text(colour = 'black', size = 13),
        legend.background = element_rect(fill = "white", color = NA)) +
  labs(tag = "A",
       x = "Beak length (mm)",
       y = "Beak depth (mm)", 
       size = 12, alpha = 1); ggp.fit.no.model


# Marginal rug plot ------------------------------------------------------------------------------------------
mul.cex = 1
ggp.fit.no.model.points = ggp.fit.no.model + 
  geom_point(data = bird.data, 
             mapping = aes(x = avg.mbl,
                           y = avg.mbd, 
                           col = sp2.fct), 
             alpha = 0.05,
             shape = 19, #".", # if use ".", size you should increase the size of points
             show.legend = FALSE,
             size = ifelse(bird.data$mxcpois == 0, 
                           yes = .7, 
                           no = bird.data$mxcpois/max(bird.data$mxcpois)*mul.cex)) +
  # Add tick marks called 'rugs' to show the density of 
  # points compact visualisation of marginal distributions 
  geom_rug(data = bird.data, 
           mapping = aes(x = avg.mbl,
                         y = avg.mbd, 
                         col = sp2.fct), 
           show.legend = FALSE,
           alpha = .1) +
  guides(#fill = guide_legend(order = 1), 
         colour = guide_legend(order = 2, ncol=1,
                               override.aes = list(size = 3, 
                                                   alpha = 1))) +
  scale_color_manual(values = pal, 
                     name = "Species",
                     labels = c(bquote(italic("G. fuliginosa")),
                                bquote(italic("G. fortis")*" small"),
                                bquote(italic("G. fortis")*" large"),
                                bquote(italic("G. magnirostris")),
                                bquote(italic("G. scandens"))
                     ));ggp.fit.no.model.points


# Marginal plots manual density ------------------------------------------------------------------------------
# Add density plot directly to ggplot 
ggp.fit.no.model.points + 
  # Making more space for the density plots 
  coord_cartesian(xlim = c(min(bird.data$avg.mbl), ceiling(max(bird.data$avg.mbl))+1),
                  ylim = c(min(bird.data$avg.mbd), ceiling(max(bird.data$avg.mbd))+1)) +
  # Adding and moving the density plot 
  geom_density(data = bird.data, 
               mapping = aes(y = avg.mbd, 
                             col = sp2.fct),
               position = position_nudge(x = max(bird.data$avg.mbl)+.1),
               inherit.aes = FALSE) +
  geom_density(data = bird.data,
               mapping = aes(x = avg.mbl,
                             col = sp2.fct),
               position = position_nudge(y = max(bird.data$avg.mbd)+.1),
               inherit.aes = FALSE)
  
# Marginal plots with cowplot -------------------------------------------------------------------------------
# Allows for more control since plotting full GGPLOTS in a grid 

# Cleaning the main plots by removing the legend 
sp <- ggp.fit.no.model.points + rremove("legend") 

# Making density plot for x axis
xplot <- ggdensity(data = bird.data, x = "avg.mbl", fill = "sp2.fct", palette = pal) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Making density plot for y axis. Note that we are providing X values, but ROTATING the plot 
yplot <- ggdensity(data = bird.data, x = "avg.mbd", fill = "sp2.fct", palette = pal) + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  rotate()

# Remove legend and rename axes. In both cases, don't print x axis
# labels as they are already present in the main plot 
xplot <- xplot + rremove("legend") + labs(y = "Density", x = "") #+ clean_theme() 
yplot <- yplot + rremove("legend") + labs(y = "Density", x = "") #+ clean_theme() 

# Arranging the plot using cowplot
cowplot::plot_grid(xplot, NULL, sp,
          yplot, ncol = 2, align = "hv", 
          rel_widths = c(5, 1), rel_heights = c(1, 5))

# Marginal plots with ggExtra --------------------------------------------------------------------------------
# Add histogram marginal plots with ggExtra
ggp.fit.no.model.points.marginals = ggExtra::ggMarginal(ggp.fit.no.model.points, 
                                                        type = "histogram", 
                                                        groupColour = TRUE, groupFill = TRUE); ggp.fit.no.model.points.marginals


# SAVING GGplots ---------------------------------------------------------------------------------------------
cat(paste("Saving ggplot to file", 
          paste("output/data.out/ggplot_data/ggp.fit.no.model.RDSgp.fit.no.model.RDS")), 
    fill = TRUE)

ggsave(filename = paste("output/images/landscape_plots/fitplot.no.model_ggplt",ext.file,".png", sep = ""),
       plot = ggp.fit.no.model.points, device = "png", units = "in", width = 12, height = 6)

saveRDS(object = ggp.fit.no.model.points, 
        file = "output/data.out/ggplot_data/ggp.fit.no.model.RDS")
