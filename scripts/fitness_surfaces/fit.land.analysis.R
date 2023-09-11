# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Analysis of the dynamics of Darwin's finches' fitness landscapes  
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created April 13, 2022
# Why: 
  # Prepare the data and generate diagnostic plots and tables 
# Requires 
  # data and scripts 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Preparation of variables  --------------------------------------
source('scripts/0.0_initialize.R')

# Data --------------------------------------------------------------------
load('data/bird.data.RData', verbose = TRUE) 
write.csv(x = bird.data, file = "data/bird.data.csv")

# Check fitness values 
gg.bar.fitness = ggplot(bird.data, 
                        mapping = aes(x = mxcpois, 
                                      col = sp2,
                                      fill = sp2)) +
  geom_bar(position = "dodge") + 
  theme_bw() + 
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12,
                                    face = "bold")) +
  labs(x = "Apparent survival (Fitness )",
       y = "Count", 
       colour = "Species", 
       fill = "Species");gg.bar.fitness

col.check = c("BANDFINAL", "Date", "Species1", "Sex0", 
              "first", "last",  "maxseen", "ch",  "maxseen.corr", "ch.corr",
              "MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth", "PC1", "PC2","Site", "By.")

# Calculate correlation between variables 
cor(bird.data[,"MedianBeakDepth"],bird.data[,"MedianBeakLength"])
cor(bird.data[,"avg.mbl"],bird.data[,"avg.mbd"])
cor(bird.data[bird.data$Species1 %in% c("fortis", "fuliginosa", "magnirostris"),"MedianBeakDepth"],
    bird.data[bird.data$Species1 %in% c("fortis", "fuliginosa", "magnirostris"),"MedianBeakLength"])


## Check age classes  -----------------------------------------------------
band.juv = bird.data %>%
  filter((age %in% c('',"fldg","j")), mxcpois > 1) %>%
  select(BANDFINAL)

# TABLE: ----------------------------------

## TABLE:  PCA loadings ---------------------------------------
pca.out = rda(bird.data[,c("MedianBeakLength","MedianBeakDepth","MedianBeakWidth")])
sum.pca.out = summary(pca.out)
loadings = scores(sum.pca.out, choices = 1:3, display = "species", scaling = 0)

# Write loadings to ".csv"
write.csv(loadings,
          file = paste("output/data.out/tables/pca.loadings",ext.file,".csv", sep = ""),
          row.names = TRUE)
sort (abs (loadings[,1]), decreasing = TRUE)
sort (abs (loadings[,2]), decreasing = TRUE)

# How much variance explained by 2 first axes
(var.pc1 = round((sum.pca.out$cont$importance[,"PC1"][2])*100, 2))
(var.pc2 = round((sum.pca.out$cont$importance[,"PC2"][2])*100, 2))
(var.pc3 = round((sum.pca.out$cont$importance[,"PC3"][2])*100, 2))
(var.pc1 + var.pc2)

## Avg each year per sp ----------------------------------
# Column names that will be keeped or removed when changing the format of the data 
names.check = c("BANDFINAL", "Tarsus", "Wing.Chord", "Mass",
                "MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth", 
                "pc1.new", "pc2.new", "sp2",
                "Year", "By.","Sex0", "Species1", 
                "maxseen", "ch", "ch.corr", "maxseen.corr", "date2"
                )

# Combine the dataset with only the CORRECTED capture history 
dat.mean.summary = cbind(bird.data[,names.check],
                         as.data.frame(bird.data$X.corr))

# Summary of the beak and body traits YEARLY
summ.mean.traits = dat.mean.summary %>% 
  pivot_longer(-c(names.check), # make a long data frame to ease the calculation (keeps only the 
                                # recapture matrix since I want to calculate the mean phenotypes 
                                # for all individuals whitin each year)
               names_to = "year",
               values_to = "obs") %>% 
  group_by(sp2) %>% # scale traits per species for all years (will be able 
                    # to see how far away each year-species is far from the overall mean)
  # New columns (used by the summary command below)
  mutate(scaled.mbl = scale(MedianBeakLength), 
         scaled.mbd = scale(MedianBeakDepth),
         scaled.mbw = scale(MedianBeakWidth),
         sp.mean.bl = mean(MedianBeakLength, na.rm = TRUE),
         sp.mean.bd = mean(MedianBeakDepth, na.rm = TRUE),
         sp.mean.bw = mean(MedianBeakWidth, na.rm = TRUE),
         sp.sd.bl   = sd(MedianBeakLength, na.rm = TRUE),
         sp.sd.bd   = sd(MedianBeakDepth, na.rm = TRUE),
         sp.sd.bw   = sd(MedianBeakWidth, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(sp2, year) %>% # group by year, since I want the YEARLY average for each species and morphotype 
  filter(obs ==1) %>%  # find the individuals that were present in a particular year only 
  # Summary 
  summarise(mean.bl   = mean(MedianBeakLength, na.rm = TRUE),       # calculate the mean
            mean.bd   = mean(MedianBeakDepth, na.rm = TRUE),        #          "
            mean.bw   = mean(MedianBeakWidth, na.rm = TRUE),        #          "
            mean.pc1  = mean(pc1.new, na.rm = TRUE),                #          "
            mean.pc2  = mean(pc2.new, na.rm = TRUE),                #          "
            median.bl = median(MedianBeakLength, na.rm = TRUE),     # calculate the median
            median.bd = median(MedianBeakDepth, na.rm = TRUE),      #          "
            median.bw = median(MedianBeakWidth, na.rm = TRUE),      #          "
            mean.s.bl = mean(scaled.mbl, na.rm = TRUE),             # calculate the mean
            mean.s.bd = mean(scaled.mbd, na.rm = TRUE),             #          "
            mean.mas  = mean(Mass, na.rm = TRUE),                   #          "
            mean.wcd  = mean(Wing.Chord, na.rm = TRUE),             #          "
            mean.tar  = mean(Tarsus, na.rm = TRUE),                 #          "
            sp.mean.uni.bl = unique(sp.mean.bl),         # Unique vals 
            sp.mean.uni.bd = unique(sp.mean.bd),         #
            sp.sd.uni.bl   = unique(sp.sd.bl),           #
            sp.sd.uni.bd   = unique(sp.sd.bd),           #
            sd.bl  = sd(MedianBeakLength),              # calculate the sd
            sd.bd  = sd(MedianBeakDepth),               #         "
            se.bl  = sem(MedianBeakLength),             # calculate the SE
            se.bd  = sem(MedianBeakDepth),              #         "
            se.bw  = sem(MedianBeakWidth),              #         "
            se.mas = sem(Mass),                         #         "
            se.wcd = sem(Wing.Chord),                   #         "
            se.tar = sem(Tarsus),                       #         "
            nb = n() # get NB of individuals that we calculated the mean from 
  ) %>% 
  # Add year column 
  mutate(yr = as.factor(as.numeric(substr(x = year, # Changes the year columns "y.XXXX" to "XXXX" (that were 
                                                    # the column names in X.corr)
                                          start = 3, 
                                          stop = 6))),
         year.fct = as.factor(year)) 

## Summary mean traits per sp -----------------------------------------------------
# Make table of species and average beak traits 
mean.across.yr.per.sp = dat.mean.summary %>% 
  pivot_longer(-c(names.check), # make a long data frame to ease the calculation 
                                # (keeps only the recapture matrix since I want 
                                # to calculate the mean phenotypes for all individuals whitin each year)
               names_to = "year", 
               values_to = "obs") %>% 
  group_by(sp2) %>%      # group by year, since I want the YEARLY average for each species and morphotype 
  filter(obs ==1) %>%    # find the individuals that were present in a particular year only 
  summarise(mean.bl   = mean(MedianBeakLength,   na.rm = TRUE), # calculate the mean
            mean.bd   = mean(MedianBeakDepth,    na.rm = TRUE),
            mean.bw   = mean(MedianBeakWidth,    na.rm = TRUE),
            median.bl = median(MedianBeakLength, na.rm = TRUE),
            median.bd = median(MedianBeakDepth,  na.rm = TRUE),
            median.bw = median(MedianBeakWidth,  na.rm = TRUE),
            nb = n()) # Nb of individuals in calculation 

# Mean standard error of beak length 
sd.mean.traits = summ.mean.traits %>% 
  group_by(sp2) %>% 
  summarise(sd.bl = sd(mean.bl),
            se.bl = sem(mean.bl),
            mean.se.bl = mean(se.bl)) 
range(sd.mean.traits$sd.bl)

# Effort plot per species -------------------------------------------------
df.bd = cbind(bird.data$X.corr, 
              sp = bird.data$sp2)

# Make dataframe with sp-year-presence 
ch.recp = pivot_longer(data = as.data.frame(df.bd),
                       names_prefix = "y.",
                       names_to = "year", 
                       values_to = "pres",
                       cols = -sp)

# Keep only the birds present 
ch.recp.1 = ch.recp[ch.recp$pres == 1,]

# relevel factors 
ch.recp.1$sp2 = levels(bird.data$sp2)[ch.recp.1$sp]

# Rename species for plotting (see as_labeller)
finch_names <- c(
  `fortis large` = "Geospiza fortis large",
  `fortis small` = "Geospiza fortis small",
  `fuliginosa` = "Geospiza fuliginosa",
  `magnirostris` = "Geospiza magnirostris",
  `scandens` = "Geospiza scandens"
)

ch.recp.1$sp3 = factor(ch.recp.1$sp2, 
                       levels = c("fuliginosa", 
                                  "fortis small", 
                                  "fortis large", 
                                  "magnirostris", 
                                  "scandens"))
levels(ch.recp.1$sp3) <- c("italic('Geospiza fuliginosa')", 
                           "italic('Geospiza fortis')*' small morphotype'", 
                           "italic('Geospiza fortis')*' large morphotype'", 
                           "italic('Geospiza magnirostris')", 
                           "italic('Geospiza scandens')" )

## GGPLOT effort per sp ----------------------------------------------------
plot.effort = ch.recp.1 %>% 
  group_by(sp3, year) %>% 
  summarise(sum = sum(pres)) %>% 
  ggplot(aes(x = year, y = sum, group = 1)) + 
  geom_point(size = 3.5) + 
  geom_line(linewidth = 1.5) +
  facet_wrap(vars(sp3), scales = "free_y", ncol = 1,
             labeller=label_parsed) + 
  theme_bw()+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.text = element_text(colour = 'black', size =20),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.ticks = element_line(colour = "black"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20, colour = "black"),
    plot.title = element_text(size = 15)) +
  scale_y_continuous(breaks= pretty_breaks()) +
  labs(x = "Year", 
       y = "Number of individuals captured", 
       title = "");plot.effort

ggsave(paste("output/images/data_exploration/effort.plot.per.species",ext.file,".png", sep = ""),
       plot = plot.effort, device = "png",
       width = 16, height = 12,units = "in")

# look at the number of species in each year 
cap.sp.yr = table(ch.recp.1$sp2, ch.recp.1$year)

# Get rownames 
rnames = rownames(t(cap.sp.yr))
# Dataframe yr = rows, species=col
cap.t = as.data.frame.matrix(t(cap.sp.yr))

cap.t$yr = rnames

# TABLE /PUB\ : ------------------------------------------------------------------
## sample size and sampling date range ---------------------------------
date.range.sampling.effort = merge(find.date.range, 
                                   cap.t, 
                                   by.x = "Year",
                                   by.y = "yr")
drg.eff = date.range.sampling.effort[,c("Year", "mintime", "maxtime", 
                                        "fortis large", "fortis small", 
                                        "fuliginosa", "magnirostris", "scandens")]
sum(drg.eff$`fortis large`)
sum(drg.eff$`fortis small`)
sum(drg.eff$fuliginosa)
margg = addmargins(as.matrix(drg.eff[,c("fortis large", "fortis small", 
                                        "fuliginosa", "magnirostris", "scandens")]))
drg.eff$min.mth.day = substr(drg.eff$mintime, 6, 10)
drg.eff$max.mth.day = substr(drg.eff$maxtime, 6, 10)
drg.eff.r = rbind(drg.eff,c("","","",margg[nrow(margg),-ncol(margg)],"",""))
drg.eff.r.c = cbind(drg.eff.r,Sum = margg[,ncol(margg)])
Percentage.sampling=as.numeric(drg.eff.r.c[nrow(drg.eff.r.c),])/drg.eff.r.c[dim(drg.eff.r.c)[1],
                                                                            dim(drg.eff.r.c)[2]]*100
drg.eff.r.c = rbind(drg.eff.r.c, Percent = round(Percentage.sampling,2))
write.csv(drg.eff.r.c[,c("Year", "min.mth.day", "max.mth.day", 
                         "fuliginosa", "fortis small", 
                         "fortis large", "magnirostris", "scandens", 
                         "Sum")], 
          file = paste("output/data.out/tables/sampling.effort.sampling.range",ext.file,".csv", sep = ""),
          row.names = FALSE)

# Get column indices
col.index = which(names(bird.data) %in% c("MedianBeakLength","MedianBeakDepth","MedianBeakWidth",
                            "Mass","Tarsus",
                            "Wing.Chord","pc1.new","pc2.new","pc1.bod","pc2.bod"))

# Find the ones that are most likely to give you the output wished
# Correlation 
cor(bird.data[,c("MedianBeakLength","MedianBeakDepth")])

# Testing the significance Correlation 
cor.all.sp.DW = cor.test(x = bird.data$MedianBeakWidth,bird.data$MedianBeakDepth,method = "pearson")
cor.all.sp = cor.test(x = bird.data$MedianBeakLength,bird.data$MedianBeakDepth,method = "pearson")
cor.per.sp = bird.data %>% 
  group_by(sp2) %>% 
  summarise(cor.sp = cor(MedianBeakLength, 
                         MedianBeakDepth),
            cor.pval = ifelse(test = cor.test(MedianBeakLength, 
                                              MedianBeakDepth)$p.value<2.2e-16,
                              yes = 2.2e-16,
                              no = cor.test(MedianBeakLength, 
                                            MedianBeakDepth)$p.value))

# GGplot pairs ------------------------------------------------------------

## Legend for ggpairs ------------------------------------------------------
# Set up the legend for ggpairs
sp.order = c(3,2,1,4,5)
plot.leg = ggplot(data = data.frame(x=1:5, 
                                    y = 1:5, 
                                    sp = factor(names(pal)[sp.order], 
                                                levels = c("fuliginosa", "fortis small", "fortis large", 
                                                           "magnirostris", "scandens"))), # reorder species names 
                  aes(x,y, col = sp)) + 
  geom_point() +   
  theme_classic() +   
  theme(legend.text.align = 0,
        axis.line = element_line(linewidth = 0.8),
        axis.ticks = element_line(colour = "black",
                                  linewidth = 0.8), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "bottom",
        strip.text.x = element_text(size = 10, colour = "black", angle = 0)) +
  scale_color_manual(values = pal, 
                     name = "", # Title of the legend
                     labels = c(bquote(italic("G. fuliginosa")),
                                bquote(italic("G. fortis")*" small"),
                                bquote(italic("G. fortis")*" large"),
                                bquote(italic("G. magnirostris")),
                                bquote(italic("G. scandens"))
                     )) + 
  guides(colour = guide_legend(override.aes = list(size=5))); plot.leg

## GGpairs -----------------------------------------------------------------
traits.pairs = GGally::ggpairs(data = bird.data, 
                               mapping = ggplot2::aes(colour = sp2.fct, 
                                                      alpha = 0.4),
                               # Relabel the tags 
                               columnLabels = c("Beak length","Beak depth","Beak width",
                                                "Mass","Tarsus", "Wing chord",
                                                "PC1","PC2","PC1 body","PC2 body"),
                               # Order the columns 
                               columns = col.index[c(4,6,5,3,1,2,7,8,9,10)],
                               # Get the legend from another place 
                               legend = grab_legend(plot.leg),
                               # from source("scripts/functions/cor_fun.R")
                      upper = list(continuous = wrap(funcVal = cor_fun, 
                                                     stars = TRUE))) + 
  theme_classic() +   
  theme(legend.text.align = 0,
        axis.line = element_line(linewidth = 0.8),
        axis.ticks = element_line(colour = "black", linewidth = 0.8), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(size = 10), # Adjusted the size of x and y axes numbers to fit on the image
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "bottom",
        strip.text.x = element_text(size = 10, colour = "black", angle = 0)) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal)

ggsave(paste("output/images/data_exploration/pairs.traits",ext.file,".png", sep = ""),
       plot = traits.pairs, device = "png",
       width = 15,height = 10,units = "in")

# GGplot density ----------------------------------------------------------
## Density on beak traits (not survival) ------------------------------------------
ft.no.surv = ggplot(bird.data, aes(x=MedianBeakLength, y=MedianBeakDepth) ) +
  geom_point(data = bird.data, aes(x=MedianBeakLength, y=MedianBeakDepth, 
                                   colour = sp2.fct),
             inherit.aes = FALSE, show.legend = TRUE)+
  scale_color_manual(values = alpha(pal,1), 
                     name = "Species",
                     labels = c(bquote(italic("G. fuliginosa")),
                                bquote(italic("G. fortis")*" small"),
                                bquote(italic("G. fortis")*" large"),
                                bquote(italic("G. magnirostris")),
                                bquote(italic("G. scandens"))
                                ))+
  stat_density_2d(geom = "polygon", 
                  contour = TRUE,
                  breaks = seq(0, .3,length.out = 20),
                  aes(fill = after_stat(level)), colour = alpha("black",.7),
                  bins = 5)+
  scale_fill_viridis_c(alpha = .6, name = "Density level", guide = "none")+
  theme_bw() + 
  theme(legend.text.align = 0,
        axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.8), 
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 16),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = "topleft",
    legend.position="none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)) +
  labs(x = "Beak length (mm)", 
       tag = "A",
       y = "Beak depth (mm)"); ft.no.surv 

## Density on PC traits (not survival) ------------------------------------------
ft.no.surv.pca = ggplot(bird.data, aes(x=pc1.new, y=pc2.new) ) +
  geom_point(data = bird.data, aes(x=pc1.new, y=pc2.new, colour = sp2.fct),
             inherit.aes = FALSE, show.legend = TRUE)+
  scale_color_manual(values = alpha(pal,1), 
                     name = "Species",
                     labels = c(bquote(italic("G. fuliginosa")),
                                bquote(italic("G. fortis")*" small"),
                                bquote(italic("G. fortis")*" large"),
                                bquote(italic("G. magnirostris")),
                                bquote(italic("G. scandens"))
                                ))+
  stat_density_2d(geom = "polygon", 
                  contour = TRUE,
                  breaks = seq(0, 9,length.out = 9),
                  aes(fill = after_stat(level)), colour = alpha("black",.7),
                  bins = 5)+
  scale_fill_viridis_c(alpha = .6, name = "Density level",guide = "none", option = "plasma")+
  theme_bw() + 
  theme(legend.text.align = 0,
        axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.8), 
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 16),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = "topleft",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)) +
  labs(x = paste0("Beak size - PC1 ", var.pc1,"%"), 
       tag = "B",
       y = paste0("Beak shape - PC2 ", var.pc2,"%")); ft.no.surv.pca 

# # GGplot density grid -------------------------------------------------------------
prow <- plot_grid( ft.no.surv + theme(legend.position="none"),
                   ft.no.surv.pca + theme(legend.position="none"),
                   label_size = 18,
                   hjust = -1,
                   vjust = 1.5,
                   nrow = 1
)
# extract the legend from one of the plots
# (Only makes sense if all plots
# have the same legend, so we can arbitrarily pick one.)
legend_b <- get_legend(ft.no.surv + 
                         theme(legend.position="bottom") + 
                         guides(colour = guide_legend(override.aes = list(size=5))) )
# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))

ggsave(paste("output/images/data_exploration/density.traits_and_pca",ext.file,".png", sep = ""),
       plot = p, device = "png",
       width = 9,height = 5,units = "in")

# Poisson data ------------------------------------------------------------
# Generating count data from the "inferred" (corrected) capture history
nmax = c(2020-2003)
bird.data %>% 
  group_by(sp2) %>% 
  mutate(rate.id = maxseen/(nmax)) %>% # proportion of the total study time for which an individual was seen
  summarise(m.rate = mean(rate.id)*100, # Average proportion for each species. 
            exp.yr = m.rate/100*nmax)   # Expected number of year each individual was seen 
                                        # in the data out of the total number of years in the study

mean(bird.data$maxseen/(nmax))*100
mean(bird.data$maxseen.corr/(nmax))*100
mean(apply(bird.data$X[,7:ncol(bird.data$X)],1,sum)/(nmax))*100


# TABLE: -------------------------------------------------------------------

## TABLE /PUB\ : recapnumber total percentage -------------------------------------
recapnumber = table(bird.data$sp2.fct, 
                    bird.data$mxcpois+1) # this is what was modeled with only mxcpois
tab.recap  = t(recapnumber) # Transposed version 
# Add totals 
recapnumber.marg = addmargins(recapnumber)
recapnumber.mat = as.data.frame.matrix(recapnumber.marg) 

# Add missing year (no bird was seen only 10 years in a row)
recapnumber.all = tibble::add_column(recapnumber.mat, 
                             `10`=rep(0, nrow(recapnumber.mat)), 
                             .after = 9)
# Add percentage in columns  (as a new row)
rcpnb.dim = dim(recapnumber.all)
Percentage.recp = recapnumber.all[rcpnb.dim[1],]/recapnumber.all[rcpnb.dim[1],rcpnb.dim[2]]*100
recapnumber.all.p = rbind(recapnumber.all, 
                          Percent = round(Percentage.recp,1))

fq.obs = c(recapnumber.all.p[nrow(recapnumber.all.p)-1,-ncol(recapnumber.all.p)]/max(recapnumber.all.p))
# sd(bird.data$MedianBeakLength)
recap.num = as.numeric(names(recapnumber.all.p)[-length(names(recapnumber.all.p))])
recap.tot = recapnumber.all.p[nrow(recapnumber.all.p)-1,-ncol(recapnumber.all.p)]

weights.mean = as.numeric(names(recapnumber.all.p)[-length(names(recapnumber.all.p))])*unlist(fq.obs)

# Add percentage in rows (as a new column)
percent.last = apply(X = recapnumber.all.p[,-which(names(recapnumber.all.p)=="Sum")],
                     MARGIN = 1,
                     FUN = function(x) round(sum(x)/recapnumber.all.p["Sum","Sum"]*100,1))
recapnumber.all.p = cbind(recapnumber.all.p, 
                          Percent = percent.last)
recapnumber.all.p["Percent","Percent"] <- "" # Empty cell 

# Change row and col names to "total"
rownames(recapnumber.all.p)[which(rownames(recapnumber.all.p) == "Sum")] <- "Total"
colnames(recapnumber.all.p)[which(colnames(recapnumber.all.p) == "Sum")] <- "Total"

# Add Geospiza at the front of species names 
rownames(recapnumber.all.p)[1:5] <- paste("Geospiza", rownames(recapnumber.all.p)[1:5])

# Export table 
write.csv(recapnumber.all.p, 
          file = paste("output/data.out/tables/recapnumber",ext.file,".csv", sep = ""),)

# Save data ---------------------------------------------------------------
saveRDS(summ.mean.traits, file = "output/data.out/fitness.landscape.data/summary.traits.finches.RDS")

# exploring fitness without model -----------------------------------------
# Please keep this script at this position (not at the top!)
source("scripts/fitness_surfaces/plot_fitness_no_model.R")

