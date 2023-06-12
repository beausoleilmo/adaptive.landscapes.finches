# Description -------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created July 26, 2022
# Why: 
  #  simulation of an adaptive landscape based on the models that is calculated from the finches data 
# Requires
  # the model file, the original dataset and the fitness surface estimation from the model 

# NOTES:
  # The adaptive landscape was calculated based on the all the phenotypes and the fitness surface.
  # Basically, you need to MOVE the phenotypes of your individuals to a specific location ()

save.adaptive.landscape.data = FALSE # If TRUE, the adaptive landscape will be computed. 
                                     # If you already computed the adaptive landscape, 
                                     # you can set it to false (run faster)
export.png = TRUE
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# SETUP -------------------------------------------------------------------
## Preparation of variables and data  --------------------------------------
source("scripts/0.0_initialize.R")

## Load the data -----------------------------------------------------------
load('data/bird.data.RData', verbose = TRUE)
bird.d = bird.data

my.persp.data = readRDS(file = "output/data.out/landscape_data/my.persp.data.RDS")
gam3.p = readRDS(file = "output/model.out/poisson.spline.model_EG.RDS")
all.local.max = readRDS(file = "output/data.out/fitness.landscape.data/all.local.max.RDS")
w.data = readRDS(file = "output/data.out/fitness.landscape.data/w.data.RDS")
summ.mean.traits = readRDS(file = "output/data.out/fitness.landscape.data/summ.mean.traits.RDS")

# Summary data ------------------------------------------------------------
# Make a summary of the bird data (used to simulate more points)
summary.data = bird.d %>% 
  group_by(sp2) %>% 
  summarise(m.bl = mean(MedianBeakLength),
            m.bd = mean(MedianBeakDepth),
            m.bw = mean(MedianBeakWidth),
            sd.bl = sd(MedianBeakLength),
            sd.bd = sd(MedianBeakDepth),
            sd.bw = sd(MedianBeakWidth),
            cor.bl.bd = cor(MedianBeakLength,MedianBeakDepth),
            cor.bl.bw = cor(MedianBeakLength,MedianBeakWidth),
            cor.bd.bw = cor(MedianBeakDepth,MedianBeakWidth))


# Adaptive landscape function  --------------------------------------------
# Load the bird data 
# load(file = "data/bird.data.RData", verbose = TRUE) # was bird.d = readRDS(file = "output/data.out/bird.data.RDS")

if(save.adaptive.landscape.data){
  # Since calculating the adaptive landscape takes some time, you don't need to recalculate each time
  # Do it once and use the data after for plots. 
  
  ## Simulate and calculate adaptive landscape -------------------------------
  n.simul = 1000 # number of individuals (traits) to simulate 
  n.iter = 30;n.iter*n.iter # Resolution (exponential) of the adaptive landscape 40 seems fine for a 'per sp' adaptive landscape 
  empirical = NULL
  # empirical = "fortis small" # Focal species to verify the effect of simulated vs real population data 
  # specify desired contour levels:
  prob <- c(0.95,0.90,0.5)
  
  bird.data %>% 
    filter(sp2 == "fortis small") %>% 
  ggplot(data =., 
         mapping = aes(x = avg.mbl, y = avg.mbd)) +
    # geom_contour_filled() +
    geom_density_2d_filled(binwidth = .1) +
    scale_fill_manual(values =  alpha(hcl.colors(4, palette = hcl.pals()[65], rev = TRUE), alpha =  .99))+
    # geom_contour() +
    geom_point(alpha = .1) +
    # multivariate t-distribution
    stat_ellipse(type = "t", level = 0.99, color = "grey0",  linewidth = 1, alpha = .8, linetype = 2) +
    stat_ellipse(type = "t", level = 0.95, color = "grey20",  linewidth = 1, alpha = .8, linetype = 2) +
    # stat_ellipse(type = "t", level = 0.90, color = "grey20", linewidth = 1, alpha = .8, linetype = 2) +
    stat_ellipse(type = "t", level = 0.50, color = "grey90", linewidth = 1, alpha = .8, linetype = 2) +
    labs(#tag = "B",
         x = "Beak length (mm)",
         y = "Beak depth (mm)", fill = "Levels", 
         # col = "Species",
         size = 12, alpha = 1) +
    theme_bw()
  
  # More iterations since corvering a larger area, needs a larger resolution 
  all.sp =    adapt.land.fun(sp.check = levels(bird.data$sp2), 
                             data = bird.data, 
                             empirical = empirical,
                             n.simul = n.simul, 
                             n.iter = 90, # 90X90
                             uppy = 20, upx = 17, downy = 5, downx = 5)
  
  # Approximative size of the grid 
  round(all.sp$x[2] - all.sp$x[1], digits = 2)
  round(all.sp$y[2] - all.sp$y[1], digits = 2)
  
  # Saving the data ---------------------------------------------------------
  dir.create(path = "output/data.out/adpt_land.sp", showWarnings = FALSE)
  
  # Make sure empirical data is exported in new object 
  if (is.null(empirical)) {
    ext.file.emp = ext.file
  } else {
    ext.file.emp = paste(ext.file, gsub(pattern = " ", replacement = "_", x = empirical), sep = "_")
  }
  
  saveRDS(all.sp, file = paste("output/data.out/adpt_land.sp/all.sp",ext.file.emp,".RDS", sep = ""))

}

# Loading saved estiamtes of adaptive landscapes --------------------------
# If ran steps before, this will be much faster 
all.sp = readRDS(file = paste("output/data.out/adpt_land.sp/all.sp",ext.file,".RDS", sep = ""))

if(save.adaptive.landscape.data){
  all.sp = readRDS(file = paste("output/data.out/adpt_land.sp/all.sp",ext.file.emp,".RDS", sep = ""))
}

all.sp$coef.variation

# PLOT: adaptive landscape ---------------------------------------------
# Number of colours (not too many as it is difficult to read )
n.col = 10 # 10 for nice scale, 13 for better resolution betwee 0.3 and 0.4 
# Make some palettes 
color_plate.dens2 = hcl.colors(n.col, "YlOrRd", rev = TRUE, alpha = 1)
color_plate.dens2.1 = hcl.colors(n.col, "YlOrRd", rev = TRUE, alpha = .4)
color_plate.dens = viridis::plasma(n.col)
color_plate.dens3 = viridis::plasma(n.col,alpha = .1)
color_plate=viridis::viridis(n.col)
# define max z axis to be consistent between the graphs 
zlim <- c(0,1)
# Create breaks based on the number of colours specified 
brks1 = seq(zlim[1],zlim[2],length.out=length(color_plate.dens2)+1)

# Manual breaks to highlight the 5 peaks 
brks = c(0,.1, .2, 0.3, .35, .36, .37,  .4, .5,.6, .7, 0.8, 0.9, 1)

# Find overall range to set up the colour scale 
overall.range.fit.land = range(c(all.sp$z2,my.persp.data$z),na.rm = TRUE)

# GGplot version ----------------------------------------------------------

# Get the fitness landscape data 
fit.data = cbind(expand.grid(my.persp.data$m1, my.persp.data$m2), 
                 as.vector(my.persp.data$z))
# Rename columns 
colnames(fit.data) <- c("x", "y", "z")

# Mean of each trait 
mean.beak.per.sp = bird.d %>% 
  group_by(sp2) %>% 
  summarise(mean.bl = mean(avg.mbl),
            mean.bd = mean(avg.mbd))


# plotting parameters 
size.pheno = 5 # Size of points 
size.peaks = 5 # Size of triangles
text.size = 16

# Data to draw segments 
peak.pheno = cbind(w.data[,c("sp2", "avg.trait1", "avg.trait2")], 
                   all.local.max[,c("x","y","sp")])


## GGplot: FITNESS LANDSCAPE  --------------------------------------------

ggp.fit.land = ggplot(data = fit.data, 
                      mapping = aes(x = x, y = y, z = z)) + 
  # geom_contour_filled(breaks = brks) +
  geom_contour_filled(breaks = brks1) +
  geom_contour(col = alpha("black",.8), 
               linewidth = .2, breaks = brks) +
  # Add points of all individuals and draw contour to the points 
  geom_point(data = bird.d, mapping = aes(x = avg.mbl, y = avg.mbd, col = sp2.fct), color = "grey30", alpha = .02, size = 1.05, inherit.aes = FALSE) +
  geom_point(data = bird.d, mapping = aes(x = avg.mbl, y = avg.mbd, col = sp2.fct), alpha = .02, size = 1, inherit.aes = FALSE) + 
  # Average of each trait 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), color = "black", inherit.aes = FALSE, size = size.pheno+.2) + 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), inherit.aes = FALSE, size = size.pheno) + 
  # Add peaks 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), color = "black", shape = 17, inherit.aes = FALSE, size = size.peaks+.6, show.legend = FALSE) + 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17, inherit.aes = FALSE, size = size.peaks, show.legend = FALSE) + 
  geom_segment(data = peak.pheno, 
               mapping = aes(x = avg.trait1, y = avg.trait2, xend = x, yend = y), inherit.aes = FALSE) +
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "grey50", linetype = "dotted", linewidth = .25),
        panel.grid.minor = element_line(colour = "grey50", linetype = "dotted", linewidth = .25),
        axis.title = element_text(size = text.size),
        axis.text = element_text(size = text.size, colour = "black", vjust = 0.25), 
        axis.text.x = element_text(size = text.size), 
        axis.text.y = element_text(size = text.size),
        plot.title = element_text(size = text.size),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = text.size),
        legend.key = element_rect(fill = "white", color = "grey80"), 
        plot.tag =element_text(face = "bold", size = 18),
        strip.background =element_rect(fill="white", colour = "white"),
        strip.text = element_text(colour = 'black', size = 13),
        legend.background = element_rect(fill = "white")) +
  labs(tag = "B",
       x = "Beak length (mm)",
       y = "Beak depth (mm)", fill = "Levels", col = "Species",
       size = 12, alpha = 1) +
  scale_fill_manual(values =  alpha(color_plate.dens2, .99))+
  scale_color_manual(values = alpha(pal,1), 
                     name = "Species",
                     labels = c(bquote(italic("G. fuliginosa")),
                                bquote(italic("G. fortis")*" small"),
                                bquote(italic("G. fortis")*" large"),
                                bquote(italic("G. magnirostris")),
                                bquote(italic("G. scandens"))
                     ))+
  coord_cartesian(xlim = range(fit.data$x), # Make the coordinate system equivalent to the other plot 
                  ylim = range(fit.data$y)) +
  guides(fill = guide_legend(order = 1), 
         colour = guide_legend(order = 2)) +
  metR::geom_text_contour(aes(z = z), size = 2.5, 
                          label.placer = label_placer_flattest()
  ); ggp.fit.land

# Add grid to the plot 
new.gridx = data.frame(x = seq(5,16, by = 2.5))
new.gridy = data.frame(y = seq(6,18, by = 2*2))
ggp.fit.land = ggp.fit.land + 
  geom_hline(data = new.gridy, aes(yintercept = y), colour = 'grey50', linetype = "dotted", linewidth = .25) +
  geom_vline(data = new.gridx, aes(xintercept = x), colour = 'grey50', linetype = "dotted", linewidth = .25)

name.adap.file.rev1 = "fit.surf.rev1_ggpt"
ggsave(filename = paste("output/images/landscape_plots/",name.adap.file.rev1,ext.file,".png", sep = ""),
       device = "png",
       plot = ggp.fit.land, units = "in", width = 7, height = 5)


## GGplot: ADAPTIVE LANDSCAPE  --------------------------------------------
ada.data = cbind(expand.grid(all.sp$x, 
                             all.sp$y), 
                 as.vector(all.sp$z2)) 
colnames(ada.data) <- c("x", "y", "z")

ggp.adap.fit = ggplot(data = ada.data, mapping = aes(x = x, y = y, z = z)) + 
  geom_contour_filled(breaks = brks1) +
  geom_contour(col = alpha("black",.8), linewidth = .2, breaks = brks1) +
  # Add points of all individuals and draw contour to the points 
  geom_point(data = bird.d, mapping = aes(x = avg.mbl, y = avg.mbd, col = sp2.fct), color = "grey30", alpha = .02, size = 1.05, inherit.aes = FALSE) +
  geom_point(data = bird.d, mapping = aes(x = avg.mbl, y = avg.mbd, col = sp2.fct), alpha = .02, size = 1, inherit.aes = FALSE) +
  # Add the trait means 
  geom_point(data = mean.beak.per.sp,
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), color = "black", inherit.aes = FALSE, size = size.pheno+.2) + 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), inherit.aes = FALSE, size = size.pheno) + 
  # Add the peaks on the adaptive landscape 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp),color = "grey10", alpha = .5, shape = 17, inherit.aes = FALSE, size = size.peaks+.6, show.legend = FALSE) + 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17, inherit.aes = FALSE, size = size.peaks, show.legend = FALSE) + 
  # Add segments 
  geom_segment(data = peak.pheno, 
               mapping = aes(x = avg.trait1, y = avg.trait2, xend = x, yend = y), inherit.aes = FALSE) +
  theme_classic() +
  theme(axis.ticks = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "grey50", linetype = "dotted", linewidth = .25),
        panel.grid.minor = element_line(colour = "grey50", linetype = "dotted", linewidth = .25),
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
  labs(tag = "C",
       title = "Adaptive landscape based simulated bivariate normal",
       # title = "Adaptive landscape based on small G. fortis",
       x = "Mean beak length (mm)",
       y = "Mean beak depth (mm)", 
       fill = "Levels", col = "Levels",
       size = 12, alpha = 1) +
  scale_fill_manual(values =  alpha(color_plate.dens2, .99),
                    guide = guide_legend(
                      direction = "vertical",
                      title.position = "top",
                      label.position = "bottom",
                      label.hjust = 0.5,
                      label.vjust = 1,
                      label.theme = element_text(angle = 90)
                    )) +
  scale_color_manual(values = alpha(pal,1), 
                     name = "Species",
                     labels = c(bquote(italic("G. fuliginosa")),
                                bquote(italic("G. fortis")*" small"),
                                bquote(italic("G. fortis")*" large"),
                                bquote(italic("G. magnirostris")),
                                bquote(italic("G. scandens"))
                     ))+
  coord_cartesian(xlim = range(fit.data$x), # Make the coordinate system equivalent to the other plot 
                  ylim = range(fit.data$y)) +
  guides(fill = "none",
         colour = "none"
  ) +
  metR::geom_text_contour(aes(z = z), size = 2.5, 
                          label.placer = label_placer_flattest()
  ); ggp.adap.fit

gg.adapt.sim.emp = ggpubr::ggarrange(ggp.adap.fit,
                                 ggp.adap.fit2,
                                 align = "v",
                                 ncol=2, common.legend = TRUE, legend = "right");gg.adapt.sim.emp

ggsave(filename = paste("~/Desktop/adaptiveland_sim_emp",ext.file,".png", sep = ""),
# ggsave(filename = paste("~/Desktop/adaptiveland_empirical",ext.file,".png", sep = ""),
       device = "png",
       plot = gg.adapt.sim.emp, 
units = "in", width = 14, height = 7)

# Add grid 
ggp.adap.fit = ggp.adap.fit + 
  geom_hline(data = new.gridy, aes(yintercept = y), colour = 'grey50', linetype = "dotted", linewidth = .25) +
  geom_vline(data = new.gridx, aes(xintercept = x), colour = 'grey50', linetype = "dotted", linewidth = .25)


# from plot_fitness_no_model.R
ggp.fit.no.model = readRDS(file = "output/data.out/ggplot_data/ggp.fit.no.model.RDS")

## GGplot: Combine panels together -------------------------------------------------
# The grid
gg.fit.adapt = ggpubr::ggarrange(ggp.fit.no.model,
                                 ggp.fit.land, ggp.adap.fit, 
                                 align = "v",
                                 legend.grob = rbind(get_legend(ggp.fit.no.model),
                                                     get_legend(ggp.fit.land)),
                                 ncol=1, common.legend = TRUE, legend = "right");gg.fit.adapt

name.adap.file = "fit.surf.adapt.land.simulations_complete_single_ggpt"

ggsave(filename = paste("output/images/landscape_plots/", name.adap.file, ext.file,".pdf", sep = ""),
       device = "pdf",
       plot = gg.fit.adapt, units = "in", width = 7, height = 13)

# GGplot: year changing  ----------------------------------------------------------
ggp.fit.land.yr.chg = ggplot(data = fit.data, mapping = aes(x = x, y = y, z = z)) + 
  geom_contour_filled(breaks = brks1)+
  geom_contour(col = alpha("black",.8), linewidth = .2, breaks = brks1) +
  # Add the population points 
  geom_point(data = bird.d, mapping = aes(x = avg.mbl, y = avg.mbd, col = sp2.fct), color = "grey30", alpha = .02, size = 1.05, inherit.aes = FALSE) +
  geom_point(data = bird.d, mapping = aes(x = avg.mbl, y = avg.mbd, col = sp2.fct), alpha = .02, size = 1, inherit.aes = FALSE) + 
  # Aad the mean beak traits
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), color = "black", inherit.aes = FALSE, size = size.pheno+.2) + 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), inherit.aes = FALSE, size = size.pheno) + 
  # Add peaks 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), color = "black", shape = 17, inherit.aes = FALSE, size = size.peaks+.6, show.legend = FALSE) + 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17, inherit.aes = FALSE, size = size.peaks, show.legend = FALSE) + 
  # Add segment between peaks and mean phenotype
  geom_segment(data = peak.pheno, 
               mapping = aes(x = avg.trait1, y = avg.trait2, xend = x, yend = y), inherit.aes = FALSE) +
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "grey50", linetype = "dotted", linewidth = .25),
        panel.grid.minor = element_line(colour = "grey50", linetype = "dotted", linewidth = .25),
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
        legend.background = element_rect(fill = "white")) +
  labs(tag = "",
       x = "Beak length (mm)",
       y = "Beak depth (mm)", fill = "Levels", col = "Species",
       size = 12, alpha = 1) +
  scale_fill_manual(values =  alpha(color_plate.dens2, .99))+
  scale_color_manual(values = alpha(pal,1), 
                     name = "Species",
                     labels = c(bquote(italic("G. fuliginosa")),
                                bquote(italic("G. fortis")*" small"),
                                bquote(italic("G. fortis")*" large"),
                                bquote(italic("G. magnirostris")),
                                bquote(italic("G. scandens"))
                     ))+
  coord_cartesian(xlim = range(fit.data$x), # Make the coordinate system equivalent to the other plot 
                  ylim = range(fit.data$y)) +
  guides(fill = guide_legend(order = 1), 
         colour = guide_legend(order = 2)) +
  metR::geom_text_contour(aes(z = z), size = 2.5, 
                          label.placer = label_placer_flattest()
  ) + 
  new_scale_color() +
  geom_point(data = summ.mean.traits, 
             mapping = aes(x = mean.bl, y = mean.bd, 
                           color = yr, group = sp2), 
             shape = 19, 
             inherit.aes = FALSE, 
             show.legend = c(colour = TRUE, fill = FALSE)) + 
  geom_path(data = summ.mean.traits[order(summ.mean.traits$yr),], 
            mapping = aes(x = mean.bl, y = mean.bd, color = yr, 
                          group = as.factor(sp2)), inherit.aes = FALSE, 
            show.legend = FALSE) +
  scale_color_viridis_d(name = "Years", guide = guide_legend(ncol = 2)); ggp.fit.land.yr.chg

# Add grid to the plot 
new.gridx = data.frame(x = seq(5,16, by = 2.5))
new.gridy = data.frame(y = seq(6,18, by = 2*2))
ggp.fit.land.yr.chg = ggp.fit.land.yr.chg + 
  geom_hline(data = new.gridy, aes(yintercept = y), colour = 'grey50', linetype = "dotted", linewidth = .25) +
  geom_vline(data = new.gridx, aes(xintercept = x), colour = 'grey50', linetype = "dotted", linewidth = .25)

ggp.fit.land.yr.chg

ggsave(filename = paste("output/images/landscape_plots/fit.surf.yr.chg_ggpt",ext.file,".png", sep = ""),
       device = "png",
       plot = ggp.fit.land.yr.chg, units = "in", width = 9, height = 8)
