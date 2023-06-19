# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Analysis of the dynamics of Darwin's finches' fitness landscapes  
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created June 1, 2023

# Why: 
  # projection pursuit searches for low-dimensional linear projections in high-dimensional data structures
  # Can be done on SCALED fitness (mxcpois.scale) OR none scaled (mxcpois) fitness. 
  # Try using the 2 (beak length and depth)and 3 traits (adding width). 

# Requires

# NOTES: 
  # See “Interesting” Projections — Where PCA Fails: 
    # https://towardsdatascience.com/interesting-projections-where-pca-fails-fe64ddca73e6 
    # Basically, PPR has a hard time for: 
      # - Unbalanced classes 
      # - The number of classes is not a power of 2 (because of the way dimensional reduction is done)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Preparation of variables and data  --------------------------------------
source('scripts/0.0_initialize.R')

load('data/bird.data.RData', verbose = TRUE)
# bird.data$avg.mbl = c(scale(bird.data$avg.mbl, center = TRUE, scale = TRUE))
# bird.data$avg.mbd = c(scale(bird.data$avg.mbd, center = TRUE, scale = TRUE))
# bird.data$avg.mbw = c(scale(bird.data$avg.mbw, center = TRUE, scale = TRUE))
bird.data$mxcpois.scale = scale(bird.data$mxcpois,center = TRUE, scale = FALSE)

table(bird.data$sp2)

bird.data.sub = bird.data #%>% 
  # filter(!c(sp2 %in% c("scandens", 
  #                      "magnirostris")))
nrow(bird.data)
nrow(bird.data.sub)

# Create directory where the models will be recorded 
dir.create("output/model.out", showWarnings = FALSE)

# MODELS ------------------------------------------------------------------
## PPR (traits) ---------------------------------------------------
### EG ----------------------------------------------------------------------
# if (site.check == "El Garrapatero") {
# finch.ppr <- ppr(mxcpois ~ avg.mbl  + avg.mbw,
#                 data = bird.data,
#                                 nterms = 2,      # number of terms to include in the final model
#                                 max.terms = 5, # maximum number of terms to choose from when building the model
#                                 sm.method = c("spline"), df = 9, trace = TRUE,
#                                 # sm.method = c("gcvspline"), gcvpen = 1.2,
#                                 # sm.method = c("supsmu"), bass = 2, span = 0,
# )

# Check if the fitness is SCALED or NOT 
# finch.ppr <- ppr(mxcpois.scale ~ avg.mbl * avg.mbd + avg.mbw + sp2,
# finch.ppr <- ppr(mxcpois.scale ~ avg.mbl * avg.mbd + avg.mbw,
# finch.ppr <- ppr(mxcpois ~ avg.mbl * avg.mbd + sp2,
finch.ppr <- ppr(mxcpois ~ avg.mbl * avg.mbd ,
                data = bird.data.sub,
                nterms = 2,      # number of terms to include in the final model
                max.terms = 2, # maximum number of terms to choose from when building the model
                sm.method = c("spline"), df = 8.8, trace = TRUE
                # sm.method = c("gcvspline"), gcvpen = 1.1
                # sm.method = c("supsmu"), bass = 1, span = 0
                ); summary(finch.ppr)

# }
par(mfrow = c(1,1))   # maybe: , pty = "s")
# Nb of new points ('resolution')
leng.tot = 50
# Generate increasing x and y 
mbl.x = seq(min(bird.data$avg.mbl), max(bird.data$avg.mbl), length.out = leng.tot)
mbd.x = seq(min(bird.data$avg.mbd), max(bird.data$avg.mbd), length.out = leng.tot)
# mbw.x = seq(min(bird.data$avg.mbw), max(bird.data$avg.mbw), length.out = leng.tot)
mbw.x = bird.data.sub %>%
  group_by(sp2.fct) %>%
  summarize(mean.bw = mean(avg.mbw)) %>% 
  as.data.frame()
# mbw.x = 6.83
mbw.x = mean(bird.data.sub$avg.mbw)
export.path = "~/Desktop/test_ppr/"

# png(filename = paste0(export.path, "/ppr.land.example%03d.png"),
# png(filename = paste0(export.path, "/ppr.land.all.sp.png"),
#     width = 10, # was 7
#     height = 4.5,units = "in", res = 300, pointsize = 12, bg = "white")
# par(mfrow = c(2,3), mar = c(4,4,2,6))
par(mfrow = c(1,1))

# for (i in 1:nrow(mbw.x)) {
# Combine to have a grid of points 
new.df.mbx = expand.grid(avg.mbl = mbl.x, 
                         # avg.mbw = i,
                         avg.mbw = mbw.x,#[i,"mean.bw"],
                         avg.mbd = mbd.x#, 
                         # sp2 = mbw.x[i,"sp2.fct"]
                         )
nrow(new.df.mbx)

# Predict using EACH POINT PAIR 
z = predict(finch.ppr, 
            newdata = new.df.mbx, 
            type = c("response"))

new.df.mbx.z = cbind(new.df.mbx, z)

# Put into a matrix 
z.mat = matrix(z, nrow = leng.tot, ncol = leng.tot)
# Draw the plot 
# image(x = mbl.x, y = mbd.x, # New x and y
#       z = z.mat,            # Model prediction 
#       col = hcl.colors(n = 12, palette = "viridis", rev = TRUE),)  # YlOrRd


# Make the breaks for the plot 
z.rge = range(z.mat)
# brk <- seq(0,z.rge[2],by = .01)
brk <- seq(0,.7,by = .01)

library(fields) # for image.plot 
# par(mfrow = c(8,8), mar = c(0,0,2,0)) 
# for (i in 1:length(hcl.pals())) {
fields::image.plot(x = mbl.x, 
                   y = mbd.x, 
                   z = z.mat, 
                   nlevel = length(brk)-1,
                   breaks = brk,
                   col=hcl.colors(length(brk)-1, 
                                  palette = "Zissou 1",  # "YlOrBr", hcl.pals()
                                  alpha = NULL, rev = FALSE),
                   # main = "",
                   # main = "PPR all species, mbl*mbd",
                   main = "PPR all species, mbl*mbd+mbw (at mean mbw)",
                   # main = mbw.x[i,"sp2.fct"],
                   ylab ="Beak depth (mm)",
                   xlab ="Beak length (mm)")
# }
# Add contour lines 
contour(x = mbl.x, y = mbd.x,z = z.mat,
        add = T, 
        levels = seq(0, z.rge[2], by = .1),
        nlevels = 10)

# Multiplication for CEX 
mul.cex = 5 # was 10 
# Check min value for cex 
min(bird.data[bird.data$mxcpois>0, "mxcpois"]/max(bird.data[bird.data$mxcpois>0, "mxcpois"])*mul.cex)

# Add each species on top 
points(bird.data$avg.mbl, 
       bird.data$avg.mbd, 
       pch = ifelse(bird.data$mxcpois == 0, 
                    yes = 25, 
                    no = 21), 
       col = scales::alpha("black", alpha = .2),
       bg = scales::alpha(pal, alpha = .2)[bird.data$sp2],  # was .5 alpha 
       cex = ifelse(bird.data$mxcpois == 0, 
                    yes = .7, 
                    no = bird.data$mxcpois/max(bird.data$mxcpois)*mul.cex))
# }
# dev.off()


plot(finch.ppr)
class(finch.ppr)
par(mfrow = c(3,2))   # maybe: , pty = "s")
plot(finch.ppr, main = "ppr(mxcpois ~ avg.mbl + avg.mbd + avg.mbw)")
plot(update(finch.ppr, bass = 5), main = "update(..., bass = 5)")
plot(update(finch.ppr, sm.method = "gcv", gcvpen = 2),
     main = "update(..., sm.method=\"gcv\", gcvpen=2)")
cbind(perm = finch.ppr$perm, prediction = round(exp(predict(finch.ppr)), 1))


# GGplot -----------------------------------------------------------------------------------------------------



size.pheno = 3 
text.size = 12

ggp.ppr = ggplot(data = new.df.mbx.z, mapping = aes(x = avg.mbl, y = avg.mbd, z = z)) + 
  geom_contour_filled(breaks = brk, show.legend = FALSE)+
  geom_contour(col = alpha("black",.8), linewidth = .2, breaks = brk) +
  geom_point(data = bird.data, 
             mapping = aes(x = avg.mbl, y = avg.mbd, color = sp2), 
             inherit.aes = FALSE, size = size.pheno, 
             show.legend = FALSE) + 
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
  labs(tag = "A",
       x = "Beak length (mm)",
       y = "Beak depth (mm)", fill = "Levels", col = "Species",
       size = 12, alpha = 1) +
  # scale_fill_manual(values =  alpha(color_plate.gray, .99)) +
  # scale_color_manual(values = alpha(pal,1), 
  #                    name = "Species",
  #                    labels = c(bquote(italic("G. fuliginosa")),
  #                               bquote(italic("G. fortis")*" small"),
  #                               bquote(italic("G. fortis")*" large"),
  #                               bquote(italic("G. magnirostris")),
  #                               bquote(italic("G. scandens"))
  #                    ))+
  coord_cartesian(xlim = range(new.df.mbx.z$avg.mbl), # Make the coordinate system equivalent to the other plot 
                  ylim = range(new.df.mbx.z$avg.mbd)) +
  guides(fill = guide_legend(order = 1), 
         colour = guide_legend(order = 2)) +
  metR::geom_text_contour(aes(z = z), size = 2.5, 
                          label.placer = label_placer_flattest()
  ); ggp.ppr

# Number of colours (not too many as it is difficult to read )
brk <- seq(0,round(z.rge[2],1),by = .1)
n.col = length(brk)-1
# Make some palettes 
color_plate.dens2 = hcl.colors(n.col, "YlOrRd", rev = TRUE, alpha = 1)

ggp.ppr = ggplot(new.df.mbx.z, mapping = aes(x = avg.mbl, y = avg.mbd, z = z)) + 
  geom_contour_filled(breaks = brk) +
  geom_contour(col = alpha("black",.8), 
               linewidth = .2, 
               breaks = brk) +
  # Add points of all individuals and draw contour to the points 
  geom_point(data = bird.data, mapping = aes(x = avg.mbl, y = avg.mbd, col = sp2.fct), color = "grey30", alpha = .02, size = 3.05, inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(data = bird.data, mapping = aes(x = avg.mbl, y = avg.mbd, col = sp2.fct), alpha = .2, size = 3, inherit.aes = FALSE) + 
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
  labs(tag = "D",
       x = "Beak length (mm)",
       y = "Beak depth (mm)", fill = "Levels", col = "Species",
       size = 12, alpha = 1) +
  scale_fill_manual(values =  alpha(color_plate.dens2, .99))+
  scale_color_manual(values = pal, 
                     name = "Species",
                     labels = c(bquote(italic("G. fuliginosa")),
                                bquote(italic("G. fortis")*" small"),
                                bquote(italic("G. fortis")*" large"),
                                bquote(italic("G. magnirostris")),
                                bquote(italic("G. scandens"))
                     ))+
  coord_cartesian(xlim = range(bird.data$avg.mbl), # Make the coordinate system equivalent to the other plot 
                  ylim = range(bird.data$avg.mbd)) +
  guides(fill = guide_legend(order = 1), 
         colour = guide_legend(order = 2, 
                               override.aes = list(size = 3, 
                                                   alpha = 1))) +
  # Add contour text with only positive values 
  metR::geom_text_contour(data = new.df.mbx.z[new.df.mbx.z$z >0,], 
                          mapping = aes(z = z), size = 2.5, 
                          label.placer = label_placer_flattest()
  ); ggp.ppr

ggsave(filename = paste("output/images/landscape_plots/ppr_ggpt",ext.file,".png", sep = ""),
       device = "png",
       plot = ggp.ppr, units = "in", width = 9, height = 8)
