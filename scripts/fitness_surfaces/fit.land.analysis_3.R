# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Analysis of the dynamics of Darwin's finches' fitness landscapes  
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created April 13, 2022
# Why: 
  # Calculation of tables, plots of fitness, euclidean distance, walking the fitness landscape, prospective selection 
# Requires 
  # You need the data, the models,  
  # the summary of the data 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# SETUP -------------------------------------------------------------------
## Preparation of variables and data  --------------------------------------
source('scripts/0.0_initialize.R')
load('data/bird.data.RData', verbose = TRUE)

# Load the final model 
gam3.p = readRDS(file = paste("output/model.out/poisson.spline.model",ext.file,".RDS", sep = ""))
gam3.pc.s = readRDS(file = "output/model.out/poisson.spline.model_pca.RDS")
gam.pca = summary(gam3.pc.s)

mod = gam3.p
summ.mean.traits = readRDS(file = "output/data.out/fitness.landscape.data/summary.traits.finches.RDS")

## Trait summary -----------------------------------------------------------
stats.traits = bird.data %>% 
  group_by(sp2) %>% 
  summarize(mean.avg.mbl = mean(avg.mbl),
            min.avg.mbl = min(avg.mbl),
            max.avg.mbl = max(avg.mbl),
            qut.avg.mbl.15 = quantile(avg.mbl, probs = 0.15),
            qut.avg.mbl.85 = quantile(avg.mbl, probs = 0.85),
            mean.avg.mbd = mean(avg.mbd),
            mean.avg.mbw = mean(avg.mbw),
            sd.avg.mbl = sd(avg.mbl),
            sd.avg.mbd = sd(avg.mbd),
            sd.avg.mbw = sd(avg.mbw),
            cor.mbl.mbd = cor(avg.mbl,avg.mbd),
            cor.mbl.mbw = cor(avg.mbl,avg.mbw),
            cor.mbd.mbw = cor(avg.mbd,avg.mbw)) %>% 
  as.data.frame()

# See how to transform from raw trait to scaled 
sp.test = "fortis large"
sd.test = stats.traits[stats.traits$sp2==sp.test,"sd.avg.mbl"]
mn.test = stats.traits[stats.traits$sp2==sp.test,"mean.avg.mbl"]
as.vector(head(bird.data[bird.data$sp2 == sp.test, "scaled.mbl"]))*sd.test+mn.test
head(bird.data[bird.data$sp2 == sp.test, "avg.mbl"])
(head(bird.data[bird.data$sp2 == sp.test, "avg.mbl"])-mn.test)/sd.test
as.vector(head(bird.data[bird.data$sp2 == sp.test, "scaled.mbl"]))

scale(bird.data[bird.data$sp2 == "magnirostris", "avg.mbl"])
mean(bird.data[bird.data$sp2 == "magnirostris",  "avg.mbl"])
sd(bird.data[bird.data$sp2 == "magnirostris",    "avg.mbl"])

# Similar to "summ.mean.traits" in the fit.land.analysis.R script
summary.bird.dat = bird.data %>%   
  group_by(sp2, Year) %>% # group by year, since I want the YEARLY average for each species and morphotype 
  summarise(mean.bl  = mean(MedianBeakLength, na.rm = TRUE), # calculate the mean
            mean.bd  = mean(MedianBeakDepth, na.rm = TRUE),
            mean.bw  = mean(MedianBeakWidth, na.rm = TRUE),
          mean.s.bl  = mean(scaled.mbl, na.rm = TRUE),
            mean.s.bd  = mean(scaled.mbd, na.rm = TRUE),
            sp.mean.uni.bl= unique(sp.mean.bl),
            sp.mean.uni.bd= unique(sp.mean.bd),
            sp.sd.uni.bl= unique(sp.sd.bl),
            sp.sd.uni.bd= unique(sp.sd.bd),
            mean.mas = mean(Mass, na.rm = TRUE),
            mean.wcd = mean(Wing.Chord, na.rm = TRUE),
            mean.tar = mean(Tarsus, na.rm = TRUE),
            sd.bl  = sd(MedianBeakLength), # calculate the sd
            sd.bd  = sd(MedianBeakDepth), # calculate the sd
            se.bl  = sem(MedianBeakLength), # calculate the SE
            se.bd  = sem(MedianBeakDepth),
            se.bw  = sem(MedianBeakWidth),
            se.mas = sem(Mass),
            se.wcd = sem(Wing.Chord),
            se.tar = sem(Tarsus),
            nb = n(),# get the number of individuals that we calculated the mean from 
  )

testline=4;(summary.bird.dat$mean.bl[testline]-summary.bird.dat$sp.mean.uni.bl[testline])/summary.bird.dat$sp.sd.uni.bl[testline]


## Set other parameters ----------------------------------------------------
bss = c("tp")
print.gif = FALSE
traits = c("avg.mbl","avg.mbd","avg.mbw")
view.traits = c("avg.mbl","avg.mbd")

## TABLE:  correlation and covariance ---------------------------------------
varcov = cov(bird.data[,c(traits[1],traits[2],traits[3])])
cor.mat = cor(bird.data[,c(traits[1],traits[2],traits[3])])
corvarcov = varcov
corvarcov[lower.tri(corvarcov)] <- cor.mat[lower.tri(cor.mat)]

## TABLE: model smooth parameters ----------------------------------------------------
mod.summ=summary(gam3.p)
rownames(mod.summ$s.table) <- c("mbl","mbd","mbl, mbd")
colnames(mod.summ$s.table) <- c("Estimated degrees of freedom","Reference degrees of freedom","Chi square", "P-value")
# Save parameters from model 
write.csv(x = round(mod.summ$s.table[,c(1,3:4)],2), 
          file = "output/model.out/spline.model.parameters.csv")

# TABLE /PUB\ MODEL Diagnostic --------------------------------------------------------
summary.mod = summary(mod) # regression coefficients
# Formule 
summary.mod$formula
# Parametric coefficients terms
summary.mod$p.table
# Smooth terms
summary.mod$s.table
summary.mod$chi.sq

log(mod$sp)    # "best" value for the smoothing parameter
mod$gcv.ubre   # the corresponding (minimum) GCV score
sum(mod$hat)   # effective number of parameters of fitted curve

## PLOT traits my.persp.data -----------------------------------------------
cat("Extracting model information to plot landscape...",fill = TRUE)
par(mfrow = c(1,2))
kk = 6
cat("   For traits.",fill = TRUE)
my.persp.data = plot.gam.cust(
  mod = mod,
  plots = c(FALSE,FALSE),
  data = bird.data,
  sp="sp2", 
  survival = "mxcpois",
  n.grid = 300, 
  nCol = 200,
  col.sch = "terrain", bdr = TRUE, 
  too.far = .15, 
  type = "response",
  view = view.traits,
  zlim = c(0,2),
  bss = bss, kk = kk,
  se = 2, # NULL
  theta=-120, phi=30,
  my.col = c("black","orange","darkgreen","pink"),
  title = "GAM PC1-2,")

## PLOT PCA my.persp.data -----------------------------------------------
cat("   For PCA.",fill = TRUE)
# For PCA 
my.persp.data_pca = plot.gam.cust(
  mod = gam3.pc.s,
  plots = c(FALSE,FALSE),
  data = bird.data,
  sp="sp2", 
  survival = "mxcpois",
  n.grid = 300, 
  nCol = 200,
  col.sch = "terrain", bdr = TRUE, 
  too.far = .15,
  type = "response",
  view = c("pc1.new","pc2.new"),
  zlim = c(0,2),
  bss = bss, kk = kk,
  se = 2,
  theta=-120, phi=30,
  my.col = c("black","orange","darkgreen","pink"),
  title = "GAM PC1-2,")

## SAVE plotting data ------------------------------------------------------
# Especially useful when plotting the SE for the fitness surface 
saveRDS(object = my.persp.data, file = "output/data.out/landscape_data/my.persp.data.RDS")
saveRDS(object = my.persp.data_pca, file = "output/data.out/landscape_data/my.persp.data_pca.RDS")

# w.data -------------------
## Mean fitness (w.data) ------------------------------------------------------------
# Make w.data which contains the average, sd of different traits the average 
# fitness and the predicted fitness at the trait average 

# Add the predicted fitness to the original data 
pred.spline.mod = predict(mod, 
                      newdata = bird.data[, view.traits], 
                      type = "response")
bird.data$pred.spline = pred.spline.mod

# The average of the WHOLE population regardless of the year
w.data.nopal = bird.data %>% 
  group_by(sp2) %>% 
  dplyr::summarize(avg.trait1 = mean(get(traits[1])),
                   avg.trait2 = mean(get(traits[2])),
                   sd.trait1 = sd(get(traits[1])),
                   sd.trait2 = sd(get(traits[2])),
                   avg.pc1 = mean(pc1.new),
                   avg.pc2 = mean(pc2.new),
                   w.bar = mean(pred.spline)  ) %>% # Average fitness for all the individuals for each species 
  mutate(pd.z = predict(mod, # Predicted fitness based on the average traits ONLY for each species 
                        newdata = data.frame(avg.mbl = avg.trait1,
                                             avg.mbd = avg.trait2),
                        type = "response"))

# Combine the palette to get the colours 
pal.df = as.data.frame(pal)
pal.df$sp2=rownames(pal.df)
w.data = merge(w.data.nopal, 
               pal.df, by.x = "sp2",by.y = "sp2")

# All local max found on Fit_land ----------------------------------------------------------
cat("Extracting local maxima from landscapes...",fill = TRUE)
# It is suggested to check manually first (locator = TRUE) and then add the coordinate so 
# that it can reach the values automatically (locator = FALSE). The values selected below 
# were already manually found. 

# find local fitness maxima on landscape 
if (site.check == "El Garrapatero") {
  # Beak traits model 
  fit.surf.iso.for.l = isolate.surf(dat = my.persp.data, xmin = 11.2, xmax = 13.5, ymin = 13,ymax = 15,   plot = T, locator = FALSE)
  fit.surf.iso.for.s = isolate.surf(dat = my.persp.data, xmin = 9.5,  xmax = 11.5, ymin = 9, ymax = 12.5, plot = T, locator = FALSE)
  fit.surf.iso.ful   = isolate.surf(dat = my.persp.data, xmin = 7.5,  xmax = 9.5,  ymin = 6, ymax = 9.5 , plot = T, locator = FALSE)
  fit.surf.iso.mag   = isolate.surf(dat = my.persp.data, xmin = 13,   xmax = 17,   ymin = 14,ymax = 17,   plot = T, locator = FALSE)
  fit.surf.iso.sca   = isolate.surf(dat = my.persp.data, xmin = 12,   xmax = 14,   ymin = 6, ymax = 22,   plot = T, locator = FALSE)
  
  # PCA model 
  fit.surf.iso.pca.for.l = isolate.surf(dat = my.persp.data_pca, xmin = .25, xmax = .32,  ymin = -.29, ymax = 0,  plot = T, points.id = bird.data[,c("pc1.new","pc2.new", "sp2")] ,locator = FALSE)
  fit.surf.iso.pca.for.s = isolate.surf(dat = my.persp.data_pca, xmin = -.1, xmax = .2,   ymin = -.5,  ymax = 0,  plot = T, points.id = bird.data[,c("pc1.new","pc2.new", "sp2")] ,locator = FALSE)
  fit.surf.iso.pca.ful   = isolate.surf(dat = my.persp.data_pca, xmin = -.4, xmax = -.25, ymin = -.5,  ymax = 0 , plot = T, points.id = bird.data[,c("pc1.new","pc2.new", "sp2")] ,locator = FALSE)
  fit.surf.iso.pca.mag   = isolate.surf(dat = my.persp.data_pca, xmin = .4,  xmax = .8,   ymin = -.4,  ymax = 0,  plot = T, points.id = bird.data[,c("pc1.new","pc2.new", "sp2")] ,locator = FALSE)
  fit.surf.iso.pca.sca   = isolate.surf(dat = my.persp.data_pca, xmin = -.2, xmax = .6,   ymin = .3,   ymax = 1,  plot = T, points.id = bird.data[,c("pc1.new","pc2.new", "sp2")] ,locator = FALSE)
}

# Combine local maxima of fitness landscape for each species 
# add in order of levels(bird.data$sp2)
all.local.max = rbind(
  fit.surf.iso.for.l$local.max,
  fit.surf.iso.for.s$local.max, 
  fit.surf.iso.ful$local.max,
  fit.surf.iso.mag$local.max,
  fit.surf.iso.sca$local.max)
# Add species 
all.local.max$sp = levels(bird.data$sp2)
# Rename column with trait 
all.local.max$beak.l = all.local.max$x
all.local.max$beak.d = all.local.max$y
row.names(all.local.max) = all.local.max$sp

# Same for PCA 
all.local.pca.max = rbind(
  fit.surf.iso.pca.for.l$local.max,
  fit.surf.iso.pca.for.s$local.max, 
  fit.surf.iso.pca.ful$local.max,
  fit.surf.iso.pca.mag$local.max,
  fit.surf.iso.pca.sca$local.max)
all.local.pca.max$sp = levels(bird.data$sp2)
all.local.pca.max$pca1 = all.local.pca.max$x
all.local.pca.max$pca2 = all.local.pca.max$y
row.names(all.local.pca.max) = all.local.pca.max$sp

# Euclidean distance ------------------------------------------------------
# Calculate Euclidean distance between features of the landscapes (mean populaiton vs fitness peak)

# Calculate standardized traits using the summary of trait 
all.local.max$scaled.all.loc.x = (all.local.max$x - stats.traits$mean.avg.mbl) / stats.traits$sd.avg.mbl
all.local.max$scaled.all.loc.y = (all.local.max$y - stats.traits$mean.avg.mbd) / stats.traits$sd.avg.mbd

## EUC: traits VS local max fit ------------------------------------------------
# Euclidean distance for each species: the local fitness maximum VS the average trait
ecl.dat = eucli.d(x1 = w.data$avg.trait1, x2 = all.local.max$x, 
                  y1 = w.data$avg.trait2, y2 = all.local.max$y)

# Same but for scaled values 
ecl.dat.scaled = eucli.d(x1 = all.local.max$scaled.all.loc.x, x2 = 0, 
                         y1 = all.local.max$scaled.all.loc.y, y2 = 0)

## EUC: trait Paiwize diff. between sp -------------------------------------
levels(bird.data$sp2.fct)
sp.ord = c("fuliginosa", "fortis small", "fortis large", "magnirostris", "scandens")

(between.species.pheno.eucli = dist.mat.all.sp(w.data[c(3,2,1,4,5),], 
                                               trait1 = "avg.trait1", 
                                               trait2 = "avg.trait2", 
                                               sp = sp.ord, 
                                               round = 2))
num.pairwize.pheno = between.species.pheno.eucli
dist.mat.all.sp(all.local.max, trait1 = "x", trait2 = "y")

# Change row and col names to fet the species names 
rownames(between.species.pheno.eucli$eucl.dist.all.comp.rd) <- paste("G.",rownames(between.species.pheno.eucli$eucl.dist.all.comp.rd))
colnames(between.species.pheno.eucli$eucl.dist.all.comp.rd) <- paste("G.",colnames(between.species.pheno.eucli$eucl.dist.all.comp.rd))
# Empty cells with NAs 
between.species.pheno.eucli$eucl.dist.all.comp.rd[which(is.na(between.species.pheno.eucli$eucl.dist.all.comp.rd))] <- ""

# raw traits 
summary(ecl.dat$eulc)
round(mean(ecl.dat$eulc),2)
round(mean(ecl.dat$angle),2)
round(range(ecl.dat$eulc),2)

# Scaled 
summary(ecl.dat.scaled$eulc)
round(mean(ecl.dat.scaled$eulc),2)
round(range(ecl.dat.scaled$eulc),2)


# MERGE: fitness.summary (w.data + all.local.max + EUC) ---------------------------------------------------------
# Merge the predicted fitness with the local maxima information 
fitness.summary = merge(w.data, # includes the average predicted fitness for all individuals of each species (w.bar) and the predicted fitness at the MEAN phenotype of all individuals of a species 
                        all.local.max, # Position and height of the local maximum of the fitness surface
                        by.x = "sp2",by.y = "sp")
# Absolute distance 
fitness.summary$d.x = abs(ecl.dat$d.x)
fitness.summary$d.y = abs(ecl.dat$d.y)
fitness.summary$eulc = ecl.dat$eulc # Euclidean distance between trait mean for each mode and their fitness peaks 
fitness.summary$eulc.angle = ecl.dat$angle
# Same with scaled data 
fitness.summary$eulc.scaled = ecl.dat.scaled$eulc
fitness.summary$eulc.sca.angle = ecl.dat.scaled$angle

# Coefficient of variation (CV)
fitness.summary$cv.trait1 = cv.var(mean = fitness.summary$avg.trait1, sd = fitness.summary$sd.trait1)
fitness.summary$cv.trait2 = cv.var(mean = fitness.summary$avg.trait2, sd = fitness.summary$sd.trait2)

# Get mean coefficient of variation (CV)
mean(fitness.summary$cv.trait1)
mean(fitness.summary$cv.trait2)

# Prepare data for pretty representation in TABLE 
avg.mead.trait1 = paste(format(round(fitness.summary$avg.trait1,2),digits = 3), 
                        " [", 
                        format(round(fitness.summary$sd.trait1,2),digits = 3),";", 
                        format(round(fitness.summary$cv.trait1,2),digits = 3),
                        "]", 
                        sep = "" )
avg.mead.trait2 = paste(format(round(fitness.summary$avg.trait2,2),digits = 3), 
                        " [", 
                        format(round(fitness.summary$sd.trait2,2),digits = 3),";", 
                        format(round(fitness.summary$cv.trait2,2),digits = 3),
                        "]", 
                        sep = "" )

delta.l.d = paste(format(round(fitness.summary$d.x,2),digits = 3), 
                  "; ", 
                  format(round(fitness.summary$d.y,2),digits = 3), 
                  sep = "" )

# GGplot: Fitland SE (B/W) --------------------------------------------------------------
my.persp.dat.frame =cbind(expand.grid(my.persp.data$m1, my.persp.data$m2), 
                          as.vector(my.persp.data$z), as.vector(my.persp.data$z.se.l), as.vector(my.persp.data$z.se.u))
colnames(my.persp.dat.frame) <- c("x", "y", "z","z.l", "z.u")
n.col = 10
# Make some palettes 
color_plate.gray = hcl.colors(n.col, "grays", rev = TRUE, alpha = 1)

zlim <- c(0,1)  

# Create breaks based on the number of colours specified 
brks = seq(zlim[1],zlim[2],length.out=length(color_plate.gray)+1)

size.pheno = 5
size.peaks = 5
text.size = 16

mean.beak.per.sp = bird.data %>% 
  group_by(sp2) %>% 
  summarise(mean.bl = mean(avg.mbl),
            mean.bd = mean(avg.mbd),
            mean.pc1 = mean(pc1.new),
            mean.pc2 = mean(pc2.new))
fit.data =cbind(expand.grid(my.persp.data$m1, my.persp.data$m2), as.vector(my.persp.data$z))
colnames(fit.data) <- c("x", "y", "z")
fit.data.pca =cbind(expand.grid(my.persp.data_pca$m1, my.persp.data_pca$m2), as.vector(my.persp.data_pca$z))
colnames(fit.data.pca) <- c("x", "y", "z")
peak.pheno = cbind(w.data[,c("sp2", "avg.trait1", "avg.trait2","avg.pc1","avg.pc2")], all.local.pca.max[, c("x","y","sp")]#, all.local.max[,c("x","y","sp")]
)

ggp.se.low = ggplot(data = my.persp.dat.frame, mapping = aes(x = x, y = y, z = z.l)) + 
  geom_contour_filled(breaks = brks, show.legend = FALSE)+
  geom_contour(col = alpha("black",.8), linewidth = .2, breaks = brks) +
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), color = "black",
             inherit.aes = FALSE, size = size.pheno+.6, 
             show.legend = FALSE) + 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), 
             inherit.aes = FALSE, size = size.pheno, 
             show.legend = FALSE) + 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17,  color = "black",
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
        legend.background = element_rect(fill = "white", color = NA)) +
  labs(tag = "A",
       x = "Beak length (mm)",
       y = "Beak depth (mm)", fill = "Levels", col = "Species",
       size = 12, alpha = 1) +
  scale_fill_manual(values =  alpha(color_plate.gray, .99)) +
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
  metR::geom_text_contour(aes(z = z.l), size = 2.5, 
                          label.placer = label_placer_flattest()
  ); ggp.se.low


mean.beak.per.sp$sp2 = factor(mean.beak.per.sp$sp2, levels = levels(bird.data$sp2.fct))
ggp.se.up = ggplot(data = my.persp.dat.frame, 
                   mapping = aes(x = x, y = y, z = z.u)) + 
  geom_contour_filled(breaks = brks)+
  geom_contour(col = alpha("black",.8), linewidth = .2, breaks = brks) +
  geom_point(data = mean.beak.per.sp,   color = "black",
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), inherit.aes = FALSE, size = size.pheno+.6) + 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.bl, y = mean.bd, color = sp2), inherit.aes = FALSE, size = size.pheno) + 
  geom_point(data = all.local.max,   color = "black", 
             mapping = aes(x = x, y = y, col = sp), shape = 17, inherit.aes = FALSE, size = size.peaks+.6, show.legend = FALSE) + 
  geom_point(data = all.local.max, 
             mapping = aes(x = x, y = y, col = sp), shape = 17, inherit.aes = FALSE, size = size.peaks, show.legend = FALSE) + 
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
  labs(tag = "B",
       x = "Beak length (mm)",
       y = "Beak depth (mm)", fill = "Levels", col = "Species",
       size = 12, alpha = 1) +
  scale_fill_manual(values =  alpha(color_plate.gray, .99), drop = FALSE) +
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
  metR::geom_text_contour(aes(z = z.u), size = 2.5, 
                          label.placer = label_placer_flattest()
  ); ggp.se.up

# The grid
gg.fit.se = ggpubr::ggarrange(ggp.se.low, 
                              ggp.se.up, 
                                 ncol=2, common.legend = TRUE, legend="right")

ggsave(filename = paste("output/images/landscape_plots/fitplot.average.change_se_ggplt",ext.file,".png", sep = ""),
       plot = gg.fit.se, device = "png", units = "in", width = 12, height = 6)

# Prospective selection ---------------------------------------------------
traits.fc = w.data[,c("avg.trait2","avg.trait1")]
row.names(traits.fc) <- w.data$sp2

# Population parameters (h2, sd)
par = structure(list(trait = c("blg", "bdt"), h2 = c(0.65, 0.79), 
                     table6_rp_diag = c(2.79, 3.64), 
                     sd_calc = c(0.0642421240945339, 0.0838140973849833), 
                     sd = c(0.06424, 0.08381),  cv = c(6.64, 9.29), 
                     avg = c(10.74, 9.56), sd_fromcv = c(0.713136, 0.888124), 
                     va = c(0.00268240544, 0.005549051719)), 
                row.names = 3:4, class = "data.frame")

# Find the additive genetic variation (h^2 = V_A/V_P, h^2*V_P = V_A)
(par$va = par$h2*(par$sd)^2)

# Genetic correlation matrix (from Boag 1983)
gen_cor <- structure(list(traits = c("wng", "trs", "blg", "bdt", "bwd"), 
                          wng = c(1, 0.68, 0.95, 0.87, 0.78), 
                          trs = c(0.68, 1, 0.71, 0.75, 0.61), 
                          blg = c(0.95, 0.71, 1, 0.9, 0.89), 
                          bdt = c(0.87, 0.75, 0.9, 1, 0.93), 
                          bwd = c(0.78, 0.61, 0.89, 0.93, 1)), 
                     row.names = c(NA, -5L), 
                     class = "data.frame")
row.names(gen_cor) <- gen_cor$traits
gen_cor = as.matrix(gen_cor[,-which(names(gen_cor) == "traits")])


# get the genetic covariance matrix from V_A and genetic correlation, equation 9.4, p. 218
gcovmat = outer(sqrt(par$va), sqrt(par$va)) * gen_cor[c("bdt","blg"),c("bdt","blg")]

# Computing equation 9.7 p.232 
(gtrans = as.matrix(log(traits.fc)) %*% solve(gcovmat[c("bdt","blg"),c("bdt","blg")]))

### ALTERNATIVE
# The values in this script were extracted from Boag 1983 and from Schluter 2000
# Extract correlation and heritability values 
table6.gen.corr.fortis.bl.bd = 0.9 # Genetic correlation beak length and depth 
table3.gen.h2.fortis.bl = 0.65 # h^2 for midparrent beak length 
table3.gen.h2.fortis.bd = 0.79 # h^2 for midparrent beak depth 

# See also table 6 for h^2 (seems there is an error in the log_10(h^2), not in log10 scale)
# Data from table 4 adults X_bar Adults in 1976
# Beak length 
fortis.blg.mean =  10.74
fortis.blg.cv = 6.64
# Extracting sd from C.V.
sd.fortis.bl = (fortis.blg.cv*fortis.blg.mean)/100
# Phenotypic variance
var.bl = sd.fortis.bl^2
# Calculating additive genetic variation (heritability * phenotypic variance)
V.A.bl = table3.gen.h2.fortis.bl*var.bl

# Beak depth 
fortis.bdt.mean = 9.56
fortis.bdt.cv = 9.29
# Extracting sd from C.V.
sd.fortis.dt = (fortis.bdt.cv*fortis.bdt.mean)/100
# Phenotypic variance
var.dt = sd.fortis.dt^2
# Calculating additive genetic variation (heritability * phenotypic variance)
V.A.bd = table3.gen.h2.fortis.bd*var.dt

# Genetic Covariance (equation 9.4, p. 218; cov = r_A*genetic standard deviation traits1 * genetic standard deviation trait2)
cov.for.bl.bd = table6.gen.corr.fortis.bl.bd*sqrt(V.A.bl)*sqrt(V.A.bd)

# This is the (additive genetic) variance and covariance
Gmat = matrix(c(V.A.bd, cov.for.bl.bd,
                cov.for.bl.bd, V.A.bl), nrow = 2, ncol = 2)
# Similar to equation 9.7 on page 232, Schluter 2000
# solve calculates the inverse of the g-matrix
(gtransformed = as.matrix(log(traits.fc[,c(1,2)])) %*% solve(Gmat))
(gtransformed.peaks = as.matrix(log(all.local.max[,c(2,1)])) %*% solve(Gmat))

# Calculate euclidean distance of the prospective selection (Beta_p)
# between species means and the nearest fitness peak
# Distance between FITNESS PEAKS and MEAN of TRAITS 
eucl.gtransfor = eucli.d(x1 = gtransformed.peaks[,1],x2 = gtransformed[,1],
                        y1 = gtransformed.peaks[,2],y2 = gtransformed[,2])

length(eucl.gtransfor$eulc)
round(eucl.gtransfor$eulc, digits = 2)
round(mean(eucl.gtransfor$eulc), digits = 2)
round(median(eucl.gtransfor$eulc), digits = 2)
round(sd(eucl.gtransfor$eulc), digits = 2)
round(range(eucl.gtransfor$eulc), digits = 2)

eucl.logtrait = eucli.d(x1 = log(traits.fc[,1]),x2 = log(all.local.max[,1]),
                        y1 = log(traits.fc[,2]),y2 = log(all.local.max[,2]))
eucl.logtrait
delta.beta = round(gtransformed.peaks - gtransformed,2)

# Add delta beta for G-beak depth and length 
paste("[",delta.beta[,2],", ", delta.beta[,1], "]", sep ="")


# TABLE /PUB\ EUC: export euclidean table ---------------------------------------------
#### Get the table of mean traits, position of peaks, euclidean distance, etc. 
df.euclid = data.frame(Species = paste("G.",as.character(fitness.summary$sp2)), 
                       avg.mead.trait1, 
                       avg.mead.trait2, 
                       peak.pos.l = round(fitness.summary$x,2),
                       peak.pos.d = round(fitness.summary$y,2),
                       delta.l.d, 
                       euc = round(fitness.summary$eulc, 2),
                       euc.sca = round(fitness.summary$eulc.scaled, 2),
                       euc.ang = round(fitness.summary$eulc.angle, 2),
                       euc.gtransformed = round(eucl.gtransfor$eulc, digits = 2),
                       delta.beta = paste("[",delta.beta[,2],", ", delta.beta[,1], "]", sep =""))

# new order of species in the data 
df.euclid = df.euclid[match(levels(bird.data$sp2.fct), rownames(df.euclid)),]

avg.summ = apply(df.euclid[,c("euc", "euc.sca", "euc.ang", "euc.gtransformed")],2, mean)
df.euclid.exp = rbind(df.euclid,c("Mean", rep("",5), round(avg.summ,2), ""))
write.csv(x = df.euclid.exp, 
          file = "output/data.out/tables/euclidean.csv")


## Distance matrix (/PUB\) ---------------------------------------------------------

sp.ord = c("fuliginosa", "fortis small", "fortis large", "magnirostris", "scandens")
# Species comparisons of the euclidean distance between mean traits and peaks 
dist.mat.all.sp(data = gtransformed.peaks[c(3,2,1,4,5),], 
                trait1 = 1, trait2 = 2, 
                sp = sp.ord)

# Euclidean distance between all species: take the comparisons and calculate the average 
# From the g-transformed species means (between the different species means)
(between.species.gtrans.pheno.eucli = dist.mat.all.sp(data = gtransformed[c(3,2,1,4,5),], 
                                                      trait1 = 1, trait2 = 2, 
                                                      sp = sp.ord, round = 2))
# Summary statistics of table 
low.tri.mat = as.numeric(num.pairwize.pheno$eucl.dist.all.comp.rd[lower.tri(num.pairwize.pheno$eucl.dist.all.comp.rd)])
mean(low.tri.mat)
range(low.tri.mat)
up.tri.mat = t(between.species.gtrans.pheno.eucli$eucl.dist.all.comp.rd)[upper.tri(t(between.species.gtrans.pheno.eucli$eucl.dist.all.comp.rd))]
mean(up.tri.mat)
range(up.tri.mat)

between.species.pheno.eucli$eucl.dist.all.comp.rd[upper.tri(between.species.pheno.eucli$eucl.dist.all.comp.rd)] <- t(between.species.gtrans.pheno.eucli$eucl.dist.all.comp.rd)[upper.tri(t(between.species.gtrans.pheno.eucli$eucl.dist.all.comp.rd))]
diag(between.species.pheno.eucli$eucl.dist.all.comp.rd) <- "-"
between.species.pheno.eucli$eucl.dist.all.comp.rd
write.csv(x = between.species.pheno.eucli$eucl.dist.all.comp.rd, 
          file = "output/data.out/tables/between.species.pheno.eucli.csv", 
          row.names = TRUE)


## PLOT: prospective selection Export for publication --------------------------------------------------
png(filename = paste("output/images/prospective_selection",ext.file,".png", sep = ""),
    width = 8, height = 4.5,
    units = "in", res = 300, pointsize = 12, bg = "white")
par(mfrow = c(1,2), mar = c(3.25,3.5,1.2,.5))
pos.tag.letter = -0.2
cex.all = 2
#### 
range.x.traits = range(all.local.max$beak.d,w.data$avg.trait2)
range.y.traits = range(all.local.max$beak.l,w.data$avg.trait1)

plot(w.data$avg.trait1~w.data$avg.trait2,
     pch = c(21,22,23,24,25), 
     xlab = "",
     ylab = "",
     xlim = range.x.traits*c(.9,1.1),
     ylim = range.y.traits,
     bg = pal[as.factor(all.local.max$sp)],
     col = "white", 
     asp = 1,
     cex = cex.all)

segments(x0 = w.data$avg.trait2-w.data$sd.trait2, y0 = w.data$avg.trait1, 
         x1 = w.data$avg.trait2+w.data$sd.trait2, y1 = w.data$avg.trait1)
segments(x0 = w.data$avg.trait2, y0 = w.data$avg.trait1-w.data$sd.trait1, 
         x1 = w.data$avg.trait2, y1 = w.data$avg.trait1+w.data$sd.trait1)
points(x = all.local.max$beak.d, 
       y = all.local.max$beak.l, 
       pch = c(21,22,23,24,25), 
       bg = "gray",
       cex = cex.all, 
       col = "grey70")
mtext(text = "A",side = 3, 
      line = 0.0,
      adj = pos.tag.letter,
      cex = 1.5, 
      font=2, col="black" )
mtext("Beak length (mm)", side=2, line=2, cex=1)
mtext("Beak depth (mm)", side=1, line=2,  cex=1)

intersect.lines(p = c(15, 9.0), 
                xlength = 5, ylength = 2, 
                angl = 45, lty = 2, col.line  = "black", 
                text1 = c(0.80, 1.40),
                text2 = c(-0.30, -0.20),
                txt1 = "ΔBeak size", 
                txt2 = "ΔBeak shape")
### Add segments between the 2 traits 
segments(x0 = all.local.max$beak.d, 
         y0 = all.local.max$beak.l,
         x1 = w.data$avg.trait2,
         y1 = w.data$avg.trait1, col = "gray50", lty = 2)


#### 
# B panel 
range.xy = apply(rbind(gtransformed,gtransformed.peaks),2,range)

plot(gtransformed[,2]~gtransformed[,1],
     pch = c(21,22,23,24,25), 
     xlab = "",
     ylab = "",
     xlim = c(range.xy[1,1],range.xy[2,1]),
     ylim = c(range.xy[1,2],range.xy[2,2]),
     asp = 1,
     bg = pal[as.factor(all.local.max$sp)], 
     col = "white", 
     cex = cex.all)
#### 
mtext(expression(bolditalic(G)*"-transformed beak length"), side=2, line=2.3, cex=1)
mtext(expression(bolditalic(G)*"-transformed beak depth"), side=1, line=2.3,  cex=1)

points(gtransformed.peaks[,1], gtransformed.peaks[,2], 
       pch = c(21,22,23,24,25), 
       bg = "gray", 
       col = "white", 
       cex = cex.all)
l1 = c(as.expression(bquote(italic("G. fuliginosa"))),
       as.expression(bquote(italic("G. fortis")*" small")),
       as.expression(bquote(italic("G. fortis")*" large")),
       as.expression(bquote(italic("G. magnirostris"))),
       as.expression(bquote(italic("G. scandens")))
)
legend("topright",legend = l1, 
       col = c(rep("black",5)),
       pt.cex = c(rep(1.3, 5)), #size points in legend 
       pch = c(23,22,21,24,25),
       bg = scales::alpha("white",.2),
       pt.bg = pal[c(3,2,1,4,5)], #"black", 
       y.intersp = 1.0, # space between each legend entry 
       cex = 1)

mtext(text = "B",side = 3, 
      line = 0.0,
      adj = pos.tag.letter,
      cex = 1.5, 
      font=2, col="black" )

### Add segments between the 2 traits 
segments(x0 = gtransformed[,1], 
         y0 = gtransformed[,2],
         x1 = gtransformed.peaks[,1], 
         y1 = gtransformed.peaks[,2], col = "gray50", lty = 2)

# draw line 90 degrees 
pts <- matrix(c(-7, 11, -5, 9, -5.5, 10.3), ncol = 2, byrow = TRUE)
intersect.lines(p = c(-8.5, 13.05), 
                xlength = 2, ylength = 4.5, 
                angl = 45, lty = 2, col.line  = "black", 
                text1 = c(0.56, 0.28), 
                text2 = c(-1.20, 0.80), 
                txt1 = "Δβ (Beak size)", 
                txt2 = "Δβ (Beak shape)")
dev.off()

# BETWEEN SPECIES ----
## Walking fitness surface fitness-peak between species ---------------------------------------------
# Peak to peak measurements 
# Makes a sequence of numbers to walk the fitness landscape 
# The image will also generate the data needed to get the species 
# make sequence definition (higher = more numbers between 2 points)
length.out.seq = 250
min.fit.sca.comp = NULL
# plot the expected fitness for the distance between the peaks 
# The focus here is ALL species to ALL the other species 
for (target in c("scandens", "fortis large", "fortis small", "fuliginosa", "magnirostris")) {
  for (i in c("fortis large", "fortis small", "fuliginosa", "magnirostris")) {
    if (target == i) {
      next
    }
    # new dataframe with sequence values of x-y axes
    pred.walk = data.frame(spcomp1 = target,    # Record which species analysed 
                           spcomp2 = i,
                           sp1.sp2 = paste(target, i, sep = "-"),
                           avg.mbl = seq(all.local.max[target,c("beak.l")], # Make a sequence of phenotypes to reach 
                                         all.local.max[i,c("beak.l")], 
                                         length.out = length.out.seq),
                           avg.mbd = seq(all.local.max[target,c("beak.d")], 
                                         all.local.max[i,c("beak.d")], 
                                         length.out = length.out.seq),
                           spcomp1.peak = all.local.max[target,c("z")], # Keep the peak height (fitness) information 
                           spcomp2.peak = all.local.max[i,c("z")],      # Keep the peak height (fitness) information 
                           res = length.out.seq)
    # generate new predicted fitness values (to find the minimum)
    pred.walk$min.fitval = predict(mod, newdata = pred.walk, type = "response")
    min.fit.sca.comp = rbind(min.fit.sca.comp, 
                             pred.walk[which.min(pred.walk$min.fitval),])
  }
}
# Subset for scandens only 
min.fit.sca.comp[min.fit.sca.comp$spcomp1=="scandens",]

## TABLE /PUB\ Different species-peaks min fitness (after the walk, scandens-all) -------------------------------------
# Calculate the difference between the fitness peak and the fitness MINIMUM value 
# between 2 species' fitness peak (i.e., maximum of fitness for a sp MINUS the minimum between the 2 species (from their 2 fitness peaks))
min.fit.sca.comp$diff.peak.min1  = min.fit.sca.comp$spcomp1.peak - min.fit.sca.comp$min.fitval
min.fit.sca.comp$diff.peak.min2  = min.fit.sca.comp$spcomp2.peak - min.fit.sca.comp$min.fitval
# Calcualte PERCENTAGE the peak difference to the MINIMUM value 
min.fit.sca.comp$percent.peak.min1  = min.fit.sca.comp$diff.peak.min1/min.fit.sca.comp$spcomp1.peak*100
min.fit.sca.comp$percent.peak.min2  = min.fit.sca.comp$diff.peak.min2/min.fit.sca.comp$spcomp2.peak*100

# Make a table of the comparison
cbind(species1 = min.fit.sca.comp$spcomp1,
      species2 = min.fit.sca.comp$spcomp2, 
      round(min.fit.sca.comp[,c("percent.peak.min1","percent.peak.min2")], 1),
      round(min.fit.sca.comp[,c("diff.peak.min1","diff.peak.min2")], 2))

min.fit.comp.exp = min.fit.sca.comp
min.fit.comp.exp [,-c(1:3)]= apply(min.fit.comp.exp[,-c(1:3)], 2, function(x) round(x, 2))

# Export all comparisons 
min.fit.comp.exp.csv = min.fit.comp.exp

min.fit.comp.exp.csv$sp2.1 = factor(min.fit.comp.exp.csv$spcomp1, levels = c("fuliginosa", "fortis small", "fortis large", "magnirostris", "scandens"))
min.fit.comp.exp.csv$sp2.2 = factor(min.fit.comp.exp.csv$spcomp2, levels = c("fuliginosa", "fortis small", "fortis large", "magnirostris", "scandens"))
min.fit.comp.exp.csv$spcomp1 = paste("G.", min.fit.comp.exp.csv$spcomp1)
min.fit.comp.exp.csv$spcomp2 = paste("G.", min.fit.comp.exp.csv$spcomp2)
names(min.fit.comp.exp.csv) <- c("Species 1", "Species 2", "sp1.sp2",
                                 "Beak length (mm)", "Beak depth (mm)", 
                                 "Fitness peak for species 1", 
                                 "Fitness peak for species 2", 
                                 "Resolution", 
                                 "Minimum fitness value", 
                                 "∆ fitness from peak to minimum sp1", # Difference in ... for species 1
                                 "∆ fitness from peak to minimum sp2", # Difference in ... for species 2
                                 "Percentage of fitness drop from the peak sp1 (%)", #of species 
                                 "Percentage of fitness drop from the peak sp2 (%)", #of species 
                                 "sp2.1", "sp2.2")
row.names(min.fit.comp.exp.csv) <- min.fit.comp.exp.csv$sp1.sp2
keep.sp.pairs = c("scandens-fortis large", "scandens-fortis small", "scandens-fuliginosa", 
  "scandens-magnirostris", "fortis large-fortis small", "fortis large-fuliginosa", 
  "fortis large-magnirostris", 
  # "fortis small-fortis large", 
  "fortis small-fuliginosa", "fortis small-magnirostris", 
  # "fuliginosa-fortis large", "fuliginosa-fortis small", 
  "fuliginosa-magnirostris"#, 
  # "magnirostris-fortis large", "magnirostris-fortis small", "magnirostris-fuliginosa"
  )
min.fit.comp.exp.csv[min.fit.comp.exp.csv$sp1.sp2 %in% keep.sp.pairs, c(1:2,c(ncol(min.fit.comp.exp.csv)-1):ncol(min.fit.comp.exp.csv))]

min.fit.comp.exp.csv.keep = min.fit.comp.exp.csv %>% 
  filter(sp1.sp2 %in% keep.sp.pairs) %>% # Keep only certain pairs 
  arrange(sp2.1, sp2.2) %>% # Reorder the data based on the factors of species names 
  dplyr::select(-c(sp1.sp2, sp2.1, sp2.2)) # Drop certain columns  
min.fit.comp.exp.csv.keep

min.fit.comp.exp.csv.keep[order(min.fit.comp.exp.csv.keep$`Species 1`, decreasing = T),]

write.csv(x = min.fit.comp.exp.csv.keep, 
          file = "output/data.out/tables/fitness.landscape.features.csv", 
          row.names = FALSE)

# WITHIN SPECIES -----
## Walking fitness surface fitness-peak-pheno-peak within species ---------------------------------------------
length.out.seq = 250
min.fit.all.sp.comp = NULL
# plot the expected fitness for the distance between the peaks 
for (target in c("scandens","fortis large", "fortis small", "fuliginosa", "magnirostris")) {
  
  # new dataframe with sequence values of x-y axes
  pred.walk = data.frame(avg.mbl = seq(all.local.max[target,c("beak.l")],
                                       w.data[w.data$sp2 %in% target, c("avg.trait1")], length.out = length.out.seq),
                         avg.mbd = seq(all.local.max[target,c("beak.d")], 
                                       w.data[w.data$sp2 %in% target, c("avg.trait2")], length.out = length.out.seq))
  pred.walk.at.phenotype.avg = data.frame(avg.mbl = w.data[w.data$sp2 %in% target, c("avg.trait1")],
                         avg.mbd = w.data[w.data$sp2 %in% target, c("avg.trait2")])
  # generate new predicted fitness values (to find the minimum)
  pred.walk$min.fitval = predict(mod, newdata = pred.walk, type = "response")
  pred.walk$pheno.avg.fitval = predict(mod, newdata = pred.walk.at.phenotype.avg, type = "response")
  # Record which species analysed 
  pred.walk$spcomp1 = target
  # Keep the peak height (fitness) information (from the model, just the peak *not necessarily* at the mean phenotypes)
  pred.walk$spcomp1.peak = all.local.max[target,c("z")]
  # Keep only the minimum value 
  min.fit.all.sp.comp = rbind(min.fit.all.sp.comp, 
                           pred.walk[which.min(pred.walk$min.fitval),])
}
## TABLE /PUB\ Different species-peaks min fitness (after the walk) -------------------------------------
# Calculate the difference between the fitness peak and the fitness MINIMUM value between 2 species' fitness peak (i.e., maximum of fitness for a sp MINUS the minimum between the 2 species (from their 2 fitness peaks))
# In MY case, these are the same 
min.fit.all.sp.comp$diff.peak.min  = min.fit.all.sp.comp$spcomp1.peak - min.fit.all.sp.comp$min.fitval
min.fit.all.sp.comp$diff.peak.pheno.avg  = min.fit.all.sp.comp$spcomp1.peak - min.fit.all.sp.comp$pheno.avg.fitval
# Calcualte PERCENTAGE the peak difference to the MINIMUM value 
min.fit.all.sp.comp$percent.peak.min  = min.fit.all.sp.comp$diff.peak.min/min.fit.all.sp.comp$spcomp1.peak

# Make a table of the comparison
fitness.diff.within.species = cbind(Species = paste("G.", min.fit.all.sp.comp$spcomp1),
      Percentage = round(min.fit.all.sp.comp$percent.peak.min*100, 1),
      Fitness.at.peak = round(min.fit.all.sp.comp$spcomp1.peak, 2),
      Minimum.fitness.value = round(min.fit.all.sp.comp$min.fitval, 2),
      Difference.peak.min = round(min.fit.all.sp.comp$diff.peak.min, 2),
      Difference.peak.pheno = round(min.fit.all.sp.comp$diff.peak.pheno.avg, 2)
      )

fitness.diff.within.species[c(4,3,2,5,1),-1] |> 
  as.data.frame() |> 
  apply(2,as.numeric) |> 
  addmargins(margin = seq_along(ncol(fitness.diff.within.species)-1), FUN = mean)
fitness.diff.within.species[,2] |> as.numeric() |> range()

write.csv(x = fitness.diff.within.species[c(4,3,2,5,1),], 
          file = "output/data.out/tables/fitness.diff.within.species.adaptive.walk.csv", row.names = FALSE)

# Number of colours (not too many as it is difficult to read )
n.col = 10
# Make some palettes 
color_plate.dens2 = hcl.colors(n.col, "YlOrRd", rev = TRUE, alpha = 1)

# GGplot: fitness landscape PCA  --------------------------------------------------
ggp.fit.land.pca = ggplot(data = fit.data.pca, 
                          mapping = aes(x = x, y = y, z = z)) + 
  geom_contour_filled(breaks = brks)+
  geom_contour(col = alpha("black",.8), linewidth = .2, breaks = brks) +
  geom_point(data = bird.data, mapping = aes(x = pc1.new, y = pc2.new, col = sp2.fct), color = "grey30", alpha = .02, size = 1.05, inherit.aes = FALSE) +
  geom_point(data = bird.data, mapping = aes(x = pc1.new, y = pc2.new, col = sp2.fct), alpha = .02, size = 1, inherit.aes = FALSE) +
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.pc1, y = mean.pc2, color = sp2), color = "black", inherit.aes = FALSE, size = size.pheno+.2) + 
  geom_point(data = mean.beak.per.sp, 
             mapping = aes(x = mean.pc1, y = mean.pc2, color = sp2), inherit.aes = FALSE, size = size.pheno) + 
  geom_vline(xintercept = 0, linetype = "dotted", color = alpha("grey50",.5)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = alpha("grey50",.5)) +
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
  labs(
    x = "PC1 - Beak size",
    y = "PC2 - Beak shape", fill = "Levels", col = "Species",
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
  coord_cartesian(xlim = range(fit.data.pca$x), # Make the coordinate system equivalent to the other plot 
                  ylim = range(fit.data.pca$y)) +
  guides(fill = guide_legend(order = 1), 
         colour = guide_legend(order = 2)) +
  metR::geom_text_contour(aes(z = z), size = 2.5, 
                          label.placer = label_placer_flattest()
  ); ggp.fit.land.pca

ggp.fit.land.pca + 
  geom_point(data = bird.data, mapping = aes(x = pc1.new, y = pc2.new, col = sp2), alpha = .1, size = 1, inherit.aes = FALSE)

ggsave(filename = paste("output/images/landscape_plots/fit.surf.PCA_ggpt",ext.file,".png", sep = ""),
       plot = ggp.fit.land.pca, device = "png", units = "in", width = 9, height = 5)

# SAVE data ---------------------------------------------------------------
saveRDS(object = w.data, file = "output/data.out/fitness.landscape.data/w.data.RDS")
saveRDS(object = all.local.max, file = "output/data.out/fitness.landscape.data/all.local.max.RDS")
saveRDS(object = all.local.pca.max, file = "output/data.out/fitness.landscape.data/all.local.max_pca.RDS")
saveRDS(object = fitness.summary, file = "output/data.out/fitness.landscape.data/fitness.summary.RDS")
saveRDS(object = min.fit.sca.comp, file = "output/data.out/fitness.landscape.data/min.fit.sca.comp.RDS")

