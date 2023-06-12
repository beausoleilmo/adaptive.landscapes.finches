# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created January 11, 2023 
# Why:
  # Walking the fitness landscape (finding minimum between each peaks of all 
  # species and within species peaks to average phenotype)
# Requires :
  # data generated from fit.land.analysis_3.R 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Preparation of variables and data  --------------------------------------
source('scripts/0.0_initialize.R')
load('data/bird.data.RData', verbose = TRUE)

# Above to get pal 
my.persp.data = readRDS(file = "output/data.out/landscape_data/my.persp.data.RDS")
all.local.max = readRDS(file = "output/data.out/fitness.landscape.data/all.local.max.RDS")
w.data = readRDS(file = "output/data.out/fitness.landscape.data/w.data.RDS")
# Load the final model 
gam3.p = readRDS(file = "output/model.out/poisson.spline.model_EG.RDS")
mod = gam3.p

# Number of colours (not too many as it is difficult to read )
n.col = 10
# Make some palettes 
color_plate.dens2 = hcl.colors(n.col, "YlOrRd", rev = TRUE, alpha = 1)
color_plate=viridis::viridis(n.col)

# Check out range(my.persp.data$z, na.rm = TRUE) to make sure this makes sense 
zlim <- c(0,1)

# Create breaks based on the number of colours specified 
brks = seq(zlim[1],zlim[2],length.out=length(color_plate)+1)

fit.land.data = cbind(expand.grid(my.persp.data$m1, my.persp.data$m2), as.vector(my.persp.data$z))
colnames(fit.land.data) <- c("x", "y", "z")

mean.beak.per.sp = bird.data %>% 
  group_by(sp2) %>% 
  summarise(mean.bl = mean(MedianBeakLength),
            mean.bd = mean(MedianBeakDepth))
mean.beak.per.sp$sp2 = factor(mean.beak.per.sp$sp2, levels = levels(bird.data$sp2.fct))

size.pheno = 5
size.peaks = 5
text.size = 16
peak.pheno = cbind(w.data[,c("sp2", "avg.trait1", "avg.trait2")], all.local.max[,c("x","y","sp")])


# make sequence definition (higher = more numbers between 2 points)
length.out.seq = 250
min.fit.sca.comp = NULL
pred.walk.all = NULL

# plot the expected fitness for the distance between the peaks 
# The focus here is ALL species to ALL the other species 
for (target in c("scandens", "fortis large", "fortis small", "fuliginosa", "magnirostris")) {
  for (i in c("scandens", "fortis large", "fortis small", "fuliginosa", "magnirostris")) {
    # new dataframe with sequence values of x-y axes
    pred.walk = data.frame(spcomp1 = target,    # Record which species analysed 
                           spcomp2 = i,
                           avg.mbl = seq(all.local.max[target,c("beak.l")], # Make phenotype sequence to reach 
                                         all.local.max[i,c("beak.l")], 
                                         length.out = length.out.seq),
                           avg.mbd = seq(all.local.max[target,c("beak.d")], 
                                         all.local.max[i,c("beak.d")], 
                                         length.out = length.out.seq),
                           spcomp1.peak = all.local.max[target,c("z")], # Keep peak height (fitness) info 
                           spcomp2.peak = all.local.max[i,c("z")],      # Keep peak height (fitness) info 
                           res = length.out.seq,
                           sp1.sp2 = paste(target, i, sep = "-"))
    
    pred.walk.pheno = data.frame(avg.mbl = seq(all.local.max[target,c("beak.l")],
                                               w.data[w.data$sp2 %in% target, c("avg.trait1")], 
                                               length.out = length.out.seq),
                                 avg.mbd = seq(all.local.max[target,c("beak.d")], 
                                               w.data[w.data$sp2 %in% target, c("avg.trait2")], 
                                               length.out = length.out.seq))
    # generate new predicted fitness values (to find the minimum)
    pred.walk$min.fitval.pheno = predict(mod, newdata = pred.walk.pheno, type = "response")
    
    # generate new predicted fitness values (to find the minimum)
    pred.walk$min.fitval = predict(mod, newdata = pred.walk, type = "response")
    pred.walk$avg.mbl.pheno = pred.walk.pheno$avg.mbl
    pred.walk$avg.mbd.pheno = pred.walk.pheno$avg.mbd
    pred.walk.min.only = pred.walk[which.min(pred.walk$min.fitval),]
    pred.walk.min.only[,c("min.fitval.pheno",
                          "avg.mbl.pheno", 
                          "avg.mbd.pheno")] <- pred.walk[which.min(pred.walk$min.fitval.pheno), 
                                                         c("min.fitval.pheno",
                                                           "avg.mbl.pheno",
                                                           "avg.mbd.pheno")]
    min.fit.sca.comp = rbind(min.fit.sca.comp, 
                             pred.walk.min.only)
    pred.walk.all = rbind(pred.walk.all, 
                          pred.walk)
  }
}


# GGplot : between species walk -------------------------------------------

ggp.walk.fit.land.bridge = ggplot(data = fit.land.data, 
                                  mapping = aes(x = x, y = y, z = z)) + 
  geom_contour_filled(breaks = brks)+
  geom_contour(col = alpha("black",.8), linewidth = .2, breaks = brks) +
  # Add path relative to scandens
  geom_point(data = pred.walk.all[pred.walk.all$spcomp1=="scandens",], 
             mapping = aes(x = avg.mbl, y = avg.mbd, size = min.fitval), 
             color = alpha("darkolivegreen3",.09), show.legend = FALSE,
             inherit.aes = FALSE) +
  # Add minimum on path 
  geom_point(data = min.fit.sca.comp[min.fit.sca.comp$spcomp1=="scandens" & 
                                       min.fit.sca.comp$spcomp2!="scandens",], 
             mapping = aes(x = avg.mbl, y = avg.mbd), 
             shape = 25, color = alpha("purple",.5), fill = alpha("purple",.5), 
             size = 4,show.legend = FALSE,
             inherit.aes = FALSE) +
  # Add path in ground finches except scandens 
  geom_point(data = pred.walk.all[pred.walk.all$sp1.sp2 %in% c("fortis large-fortis small", 
                                                               "fortis small-fuliginosa", 
                                                               "magnirostris-fortis large"),], 
             mapping = aes(x = avg.mbl, y = avg.mbd, size = min.fitval), 
             color = alpha("darkolivegreen3",.09), show.legend = FALSE,
             inherit.aes = FALSE) +
  # Add minimum on path 
  geom_point(data = min.fit.sca.comp[min.fit.sca.comp$sp1.sp2 %in% c("fortis large-fortis small", 
                                                                     "fortis small-fuliginosa", 
                                                                     "magnirostris-fortis large"),], 
             mapping = aes(x = avg.mbl, y = avg.mbd), 
             shape = 25, color = alpha("purple",.5), fill = alpha("purple",.5), 
             size = 4,show.legend = FALSE,
             inherit.aes = FALSE) +
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), color = "black", 
             inherit.aes = FALSE, size = size.pheno+.6) + 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), 
             inherit.aes = FALSE, size = size.pheno) + 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17, color = "black",
             inherit.aes = FALSE, size = size.peaks+.6, show.legend = FALSE) + 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17, 
             inherit.aes = FALSE, size = size.peaks, show.legend = FALSE) + 
  theme_classic() + # https://ggplot2.tidyverse.org/reference/theme.html
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
  labs(#title = "A", 
    tag = "A",
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
  coord_cartesian(xlim = range(fit.land.data$x), # Make the coordinate system equivalent to the other plot 
                  ylim = range(fit.land.data$y)) +
  guides(fill = guide_legend(order = 1), 
         colour = guide_legend(order = 2)) +
  metR::geom_text_contour(aes(z = z), size = 2.5, 
                          label.placer = label_placer_flattest()
  ); ggp.walk.fit.land.bridge

###


# GGplot : within species walk --------------------------------------------

ggp.walk.fit.land.pheno = ggplot(data = fit.land.data, mapping = aes(x = x, y = y, z = z)) + 
  geom_contour_filled(breaks = brks)+
  geom_contour(col = alpha("black",.8), linewidth = .2, breaks = brks) +
  # add the path 
  geom_point(data = pred.walk.all[pred.walk.all$sp1.sp2 %in%  c("scandens-scandens","fortis large-fortis large", 
                                                                "fortis small-fortis small",
                                                                "fuliginosa-fuliginosa", 
                                                                "magnirostris-magnirostris"),], 
             mapping = aes(x = avg.mbl.pheno, y = avg.mbd.pheno, size = min.fitval.pheno), 
             color = alpha("darkolivegreen3",.09), show.legend = FALSE,
             inherit.aes = FALSE) +
  # Add phenotypic mean 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), color = "black",
             inherit.aes = FALSE, size = size.pheno+.6) + 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), 
             inherit.aes = FALSE, size = size.pheno) + 
  
  geom_point(data = min.fit.sca.comp[min.fit.sca.comp$sp1.sp2 %in%  c("scandens-scandens",
                                                                      "fortis large-fortis large", 
                                                                      "fortis small-fortis small",
                                                                      "fuliginosa-fuliginosa", 
                                                                      "magnirostris-magnirostris"),], 
             mapping = aes(x = avg.mbl.pheno, y = avg.mbd.pheno), 
             shape = 25, color = alpha("purple",.5), fill = alpha("purple",.5), 
             size = 4,show.legend = FALSE,
             inherit.aes = FALSE) +
  # Add peak 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17, color = "black",
             inherit.aes = FALSE, size = size.peaks+.6, show.legend = FALSE) + 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17, 
             inherit.aes = FALSE, size = size.peaks, show.legend = FALSE) + 
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
  coord_cartesian(xlim = range(fit.land.data$x), # Make the coordinate system equivalent to the other plot 
                  ylim = range(fit.land.data$y)) +
  guides(fill = guide_legend(order = 1), 
         colour = guide_legend(order = 2)) +
  metR::geom_text_contour(aes(z = z), size = 2.5, 
                          label.placer = label_placer_flattest()
  );ggp.walk.fit.land.pheno

# The grid
gg.fit.walk = ggpubr::ggarrange(ggp.walk.fit.land.bridge, 
                                ggp.walk.fit.land.pheno, 
                                 ncol=2, common.legend = TRUE, legend="right")

ggsave(filename = paste("output/images/landscape_plots/fit.surf.walk_ggpt",ext.file,".png", sep = ""),
       plot = gg.fit.walk, device = "png", units = "in", width = 12, height = 6)


