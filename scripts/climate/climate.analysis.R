# Description -------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created March 11, 2023
# Why:
# Requires 
# NOTES: 
  # Source: 
    # https://www.darwinfoundation.org/en/datazone/climate/puerto-ayora
    # https://www.darwinfoundation.org/images/climate/climate_puerto-ayora.csv
    # https://www.ncei.noaa.gov/access/monitoring/enso/sst
    # https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt (also https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php)
    # See also links in the code 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# SETUP -------------------------------------------------------------------
## Preparation of variables and data  --------------------------------------
source("scripts/0.0_initialize.R")

dir.create(path = "data/climate", showWarnings = FALSE)

# Download the climate data  ----------------------------------------------
cat("Downloading and read the climate data", fill = TRUE)
# Record the version number of the downloaded data 
version = paste0("_V.",Sys.Date()) # Or version = NULL
###
# Check if climate data is already downloaded 
list.f.clim = list.files(path = "data/climate/")
if (length(grep(pattern = c("climate_puerto-ayora"), x = list.f.clim)) > 0) {
  # Find newest file 
  last.file = tail(sort(grep(pattern = c("climate_puerto-ayora"), x = list.f.clim)), 1)
  climate.pa <- read.csv(paste("data/climate/",list.f.clim[last.file], sep = ""))
} else {
  download.file("https://www.darwinfoundation.org/images/climate/climate_puerto-ayora.csv",
                destfile = paste("data/climate/climate_puerto-ayora",version,".csv", sep = ""))
  climate.pa = read.csv(file = paste("data/climate/climate_puerto-ayora",version,".csv", sep = ""))
  
}

list.f.enso = list.files(path = "data/climate/ENSO/")
if (length(grep(pattern = c("oni.ascii"), x = list.f.enso)) > 0) {
  # Find newest file 
  last.file = tail(sort(grep(pattern = c("oni.ascii"), x = list.f.enso)), 1)
  ONI <- read.table(file = paste("data/climate/ENSO/",list.f.enso[last.file], sep = ""), header = TRUE)
} else {
  download.file("https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
                destfile = paste("data/climate/ENSO/oni.ascii",version,".txt", sep = ""))
  ONI = read.table(file = paste("data/climate/ENSO/oni.ascii",version,".txt", sep = ""), header = TRUE)
}
###

# TABLE: climate data -----------------------------------------------------

# Date column 
climate.pa$date = climate.pa$observation_date = as.Date(climate.pa$observation_date, "%Y-%m-%d")
climate.pa$yr = lubridate::year(climate.pa$date)
climate.pa$month = lubridate::month(climate.pa$date)

# Get sum precipitation per month in every year. 
yr.range = 2000:2020
sum.permonth = climate.pa %>% 
  filter(yr %in% yr.range) %>% 
  group_by(yr, month) %>% 
  summarise(sum.per.month = sum(precipitation)) %>% 
  mutate(day = 1) # Fake day 

sum.permonth$date = as.Date(paste(sum.permonth$yr,
                                  sum.permonth$month, 
                                  sum.permonth$day, sep = "-"), "%Y-%m-%d")
# Get sum precipitation per year. 
sum.peryr = climate.pa %>% 
  filter(yr %in% yr.range) %>% 
  group_by(yr) %>% 
  summarise(sum.per.yr = sum(precipitation)) %>% 
  mutate(day = 1) # Fake day 


# PLOT: precipitation pattern ---------------------------------------------
# Plot per month
plot(sum.permonth$sum.per.month~sum.permonth$date, 
     xlab = "Year",
     ylab = "Sum precipiration per month (mm)",
     type ="l", las = 2)
abline(h = seq(0, 400, by = 25), lty = 3, col = scales::alpha("grey",.5))
abline(v = seq(as.Date("2000/01/01"), as.Date("2022/01/01"), by = "year"), lty = 3, col = scales::alpha("grey",.5))

# Plot per year
plot(sum.peryr$sum.per.yr~sum.peryr$yr, 
     xlab = "Year",
     ylab = "Sum precipiration per year (mm)",
     type ="l", las = 2)

# ONI El nino La nina  ----------------------------------------------------
# Add El Nino and La Nina data 
# Need to modify the format of the data to plot easily 
onisoiquimalypense = ONI %>% 
  filter(YR %in% 2000:2022) %>% 
  mutate(month.mean = plyr::revalue(SEAS, c("DJF" = 1, "JFM" = 2 , "FMA" = 3, "MAM" = 4, 
                                            "AMJ" = 5, "MJJ" = 6, "JJA" = 7, "JAS" = 8, "ASO" = 9, 
                                            "SON" = 10, "OND" = 11, "NDJ" = 12)),
         date = as.Date(paste(YR, month.mean, "01", sep = "-")))

onisoiquimalypense$elnlan = ifelse(onisoiquimalypense[,"ANOM"] > 0.5, 
                                   yes = "El Nino", 
                                   no = ifelse(onisoiquimalypense[,"ANOM"] < -0.5, 
                                               yes = "La Nina", 
                                               no = 'none'))

oni.sumpre = merge(sum.permonth, onisoiquimalypense, by = "date")

oni.sumpre$elnlan = as.factor(oni.sumpre$elnlan)
oni.sumpre.filt = oni.sumpre %>% 
  filter(yr >= 2002)


# GGplot precipitation pattern with ONI ------------------------------------------------------------------
el.nino.la.nina.precipitation = ggplot(data = oni.sumpre.filt, 
                                       mapping = aes(x = date, 
                                                     y = sum.per.month))+
  geom_path(aes(colour = elnlan, group = 1), 
            linewidth = .5) + 
  # Revalue the tick marks on plot 
  scale_x_date(labels = date_format("%Y"), 
               breaks = as.Date(seq(min(oni.sumpre$date), 
                                    max(oni.sumpre$date), by = "year"), 
                                format = "%Y-%m-%d"))+
  scale_y_continuous(breaks = round(seq(min(oni.sumpre$sum.per.month), 
                                        max(oni.sumpre$sum.per.month), by = 50),1)) +
  theme(axis.line = element_line(linetype = "solid"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major = element_line(colour = "gray80",linetype = "dotted"), 
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 14,
        vjust = 0.5, angle = 90), 
    axis.text.y = element_text(size = 14),
    legend.key = element_rect(fill = NA),
    legend.position="top", 
    panel.background = element_rect(fill = NA)) + 
  # ylim(c(0,300))+
  labs(x = "Year", y = "Sum precipitation per month (mm)") +
  # Change the size of lineinside legend key only 
  guides(colour = guide_legend(override.aes = list(linewidth=1.2))) + 
  scale_color_manual(name = "",
                     values = c("none" = "black", 
                                "El Nino" = "red",
                                "La Nina" = "blue")); el.nino.la.nina.precipitation

ggsave(filename = "output/images/climate/precipitation.Puerto.ayora.el.nino.la.nina.png", 
       plot = el.nino.la.nina.precipitation, device = 'png', 
       width = 9, height = 4)

# Get fuliginosa (low- highland) climate data ---------------------------------------------
load('data/bird.data.RData', verbose = TRUE)
nrow(bird.data)

mean.fuli = bird.data %>% 
  filter(sp2 =="fuliginosa", Year >= 2000) %>% 
  group_by(sp2, Year) %>% 
  summarize(mean.l = mean(MedianBeakLength),
            mean.d = mean(MedianBeakDepth),
            mean.w = mean(MedianBeakWidth),
            sem.bl = sem(MedianBeakLength),
            sem.bd = sem(MedianBeakDepth),
            sem.bw = sem(MedianBeakWidth),
            mean.fit =mean(mxcpois))
range(mean.fuli$mean.l)

mean.fuli$land = "Lowlands EG"
mean.fuli$low = mean.fuli$mean.l-mean.fuli$sem.bl
mean.fuli$upp = mean.fuli$mean.l+mean.fuli$sem.bl

kleindorferetal2006.fig.3 = read.csv("data/climate/Kleindorfer.et.al._2006/H.L.lands.fig3.Kleindorfer.2006.csv", 
                                     header = FALSE)
colnames(kleindorferetal2006.fig.3) <- c("yr", "beak.length","gr","type","land")
kleindorferetal2006.fig.3$type = trimws(kleindorferetal2006.fig.3$type)
kleindorferetal2006.fig.3$land = trimws(kleindorferetal2006.fig.3$land)
kleindorferetal2006.fig.3[kleindorferetal2006.fig.3$type == "u.hl","type"] <- "upp"
kleindorferetal2006.fig.3[kleindorferetal2006.fig.3$type == "u.sm","type"] <- "upp"
kleindorferetal2006.fig.3[kleindorferetal2006.fig.3$type == "l.hl","type"] <- "low"
kleindorferetal2006.fig.3[kleindorferetal2006.fig.3$type == "l.sm","type"] <- "low"

k2006.f.3.wide = kleindorferetal2006.fig.3 %>% 
  dplyr::select(-gr) %>% 
  pivot_wider(names_from = type, values_from = beak.length) %>% 
  as.data.frame()

k2006.f.3.wide[k2006.f.3.wide$land == "Highlands", "land"] <- "Highlands K06"
k2006.f.3.wide[k2006.f.3.wide$land == "Lowlands", "land"] <- "Lowlands K06"
k2006.f.3.wide$mean.fit = NA

mean.fuli.col.merge = mean.fuli[,c("Year","land", "mean.l", "upp", "low", "mean.fit")]
colnames(mean.fuli.col.merge) = c("yr","land", "mean", "upp", "low", "mean.fit")
all.fuli.data = rbind(k2006.f.3.wide, mean.fuli.col.merge)

fuli.preci = merge(all.fuli.data, sum.peryr, by.x = "yr", by.y = "yr")


fuli.preci.klein.fig.3 = merge(k2006.f.3.wide, sum.peryr, by.x = "yr", by.y = "yr")
fuli.preci$land = factor(fuli.preci$land, levels = c("Highlands K06",  
                                                     "Lowlands K06", "Lowlands EG"))

## range traits, precipitation -------------------------------------------------------------------
range(fuli.preci$mean)
range(fuli.preci$sum.per.yr)

# Filter data for MORE precipitation in a year than the mean across years
heav.rain = fuli.preci %>% 
  dplyr::filter(sum.per.yr > mean(fuli.preci$sum.per.yr)) %>% 
  summarize(mean.bl.heavy.rainfall = mean(mean),
            sem.bl.heavy.rainfall = sem(mean),
            nb = n())

# Filter data for LESS precipitation in a year than the mean across years
low.rain = fuli.preci %>% 
  dplyr::filter(sum.per.yr < mean(fuli.preci$sum.per.yr)) %>% 
  summarize(mean.bl.low.rainfall = mean(mean),
            sem.bl.low.rainfall = sem(mean),
            nb = n())

all.rain.large = fuli.preci %>% 
  dplyr::filter(mean > 8.5) %>% 
  summarize(mean.bl.low.rainfall = mean(mean),
            sem.bl.low.rainfall = sem(mean),
            nb = n())

low.rain.small = fuli.preci %>% 
  dplyr::filter(sum.per.yr < mean(fuli.preci$sum.per.yr),
                mean < mean(fuli.preci$mean)) %>% 
  summarize(mean.bl.low.rainfall = mean(mean),
            sem.bl.low.rainfall = sem(mean),
            nb = n())

beak.meas.ful = matrix(c((heav.rain), 
                         (low.rain), 
                         (all.rain.large), 
                         (low.rain.small)), nrow = 4, byrow = TRUE)
row.names(beak.meas.ful) <- c("beak.heavey.rain", "beak.low.rain", 
                              "beak.large.all.rain", "beak.low.low.rain")
colnames(beak.meas.ful) <- c("l", "se.l", "nb")
round(apply(beak.meas.ful, 2, as.numeric),2)


fuli.preci[is.na(fuli.preci$mean.fit),"mean.fit"] <- 0 
fuli.preci$mean.fit.rel = fuli.preci$mean.fit/max(fuli.preci$mean.fit)

## GGplot: G. fuliginosa low- highlands, climate -------------------------------------------------------------------
mean(fuli.preci$mean)
mean(fuli.preci$sum.per.yr)
mean.bl.300precp = mean(fuli.preci[fuli.preci$land=="Lowlands EG" & fuli.preci$sum.per.yr > 300 , "mean"])
round(mean.bl.300precp,2)
mean.bl.300precp.less = mean(fuli.preci[fuli.preci$land=="Lowlands EG" & fuli.preci$sum.per.yr < 300 , "mean"])
round(mean.bl.300precp.less,2)
bl.more300.rain = fuli.preci[fuli.preci$land=="Lowlands EG" & fuli.preci$sum.per.yr > 300 , "mean"]
bl.less300.rain = fuli.preci[fuli.preci$land=="Lowlands EG" & fuli.preci$sum.per.yr < 300 , "mean"]
wttest = t.test(bl.more300.rain,
       bl.less300.rain)
wttest

all.pop = rbind(data.frame(bl = bl.more300.rain, rain = ">300"),
                data.frame(bl = bl.less300.rain, rain = "<300"))
all.pop$rain = as.factor(all.pop$rain)
plot(all.pop$bl~all.pop$rain)

beak.length.precipitation = fuli.preci %>% 
  ggplot(mapping = aes(x = sum.per.yr, 
                       y = mean #, size = mean.fit.rel
                       # , col = Year)) +
                       )) +
  scale_size(range = c(0, 1)) + 
  # Adding points 
  # geom_point(size = 4, col = "gray50",show.legend = FALSE) + 
  geom_pointrange(mapping = aes(x = sum.per.yr,
                                y = mean, 
                                ymin=low, ymax=upp, col = land)) +
  
  # Adding lines in the graph 
  geom_hline(yintercept = mean(fuli.preci$mean), 
             linetype = "dashed", col = "gray40") +
  geom_vline(xintercept = mean(fuli.preci$sum.per.yr), 
             linetype = "dashed", col = "gray40") +
  # Adding text to show the years 
  geom_text_repel(mapping = aes(x = sum.per.yr, 
                                y = mean, label = yr, col = land), 
                  # col = "black", 
                  size = 3, 
                  seed = 156) +
  
  # Theme
  theme_bw() + 
  theme(axis.ticks = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray90",
                                        linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "gray90",
                                        linetype = "dotted"), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position="bottom") +
  ggthemes::scale_colour_colorblind() +
  # Axes
  labs(x = "Sum precipitation per year (mm)",
       y = "Mean beak length (mm)", colour = "Sites");beak.length.precipitation

set.seed(1234)
beak.length.precipitation.bp = ggplot(data = all.pop, mapping = aes(x = rain, 
                       y = bl
                       )) +
    geom_boxplot() +
    geom_jitter(alpha = 0.8, width = 0.1, height = 0) +
  
  # Theme
  theme_bw() + 
  theme(axis.ticks = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray90",
                                        linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "gray90",
                                        linetype = "dotted"), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  ggthemes::scale_colour_colorblind() +
  # Axes
  labs(x = "Precipitation (mm)",
       y = "Mean beak length (mm)", colour = "Sites"); beak.length.precipitation.bp

plots <- align_plots(beak.length.precipitation, beak.length.precipitation.bp, align = 'h', axis = 'b')
bl.prep.bp = cowplot::plot_grid(plotlist = list(plots[[1]], plots[[2]]), 
                                nrow = 1, labels = c("A","B"), 
                                rel_widths = c(3,1));bl.prep.bp

ggsave(filename = "output/images/climate/fuliginosa.beak.length.precipitation.png", 
       plot = bl.prep.bp, device = 'png', 
       width = 9, height = 4)

