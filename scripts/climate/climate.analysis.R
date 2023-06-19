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
dir.create(path = "data/climate/ENSO", showWarnings = FALSE)

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
yr.range = 2000:2022
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


# Make polygon  ----------------------------------------------------------------------------------------------
up = max(oni.sumpre$sum.per.month)
d1 = "2010-01-20"
d2 = "2011-03-21"
d3 = "2017-02-26"
d4 = "2019-04-16"
polygon1 = data.frame(x = c(as.POSIXct(d1)-dyears(1), as.POSIXct(d1), as.POSIXct(d1), as.POSIXct(d1)-dyears(1)), 
                      y = c(0,0, up,up), poly = 1, type = "Year before")
polygon2 = data.frame(x = c(as.POSIXct(d1), as.POSIXct(d2), as.POSIXct(d2), as.POSIXct(d1)), 
                      y = c(0,0, up,up), poly = 2, type = "Year sampling")
polygon3 = data.frame(x = c(as.POSIXct(d3)-dyears(1), as.POSIXct(d3), as.POSIXct(d3), as.POSIXct(d3)-dyears(1)),
                      y = c(0,0, up,up), poly = 3, type = "Year before")
polygon4 = data.frame(x = c(as.POSIXct(d3), as.POSIXct(d4), as.POSIXct(d4), as.POSIXct(d3)),
                      y = c(0,0, up,up), poly = 4, type = "Year sampling")

all.polygons = rbind(polygon1, polygon2, polygon3, polygon4)

class(all.polygons$x)

# GGplot precipitation pattern with ONI ------------------------------------------------------------------
el.nino.la.nina.precipitation = ggplot(data = oni.sumpre.filt, 
                                       mapping = aes(x = date, 
                                                     y = sum.per.month))+
  # Add polygon for ggplot graph 
  geom_polygon(data = all.polygons,
               mapping = aes(x = as.Date(x = x, origin = "%Y-%m-%d  %H:%M"),
                             y = y, group = poly, fill = type),
               inherit.aes = FALSE, color = NA, alpha = .5)+
  geom_path(aes(colour = elnlan, group = 1), 
            linewidth = .5) + 
  # Revalue the tick marks on plot 
  scale_x_date(labels = scales::date_format("%Y"), 
               breaks = as.Date(seq(min(oni.sumpre$date), 
                                    max(oni.sumpre$date), by = "year"), 
                                format = "%Y-%m-%d"))+
  scale_y_continuous(breaks = round(seq(min(oni.sumpre$sum.per.month), 
                                        max(oni.sumpre$sum.per.month), by = 50),1)#, 
                     # limits = range(oni.sumpre$sum.per.month)
                     ) +
  # ylim(c(0,475))+
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
  labs(x = "Year", y = "Sum precipitation per month (mm)") +
  # Change the size of lineinside legend key only 
  guides(colour = guide_legend(override.aes = list(linewidth=1.2))) + 
  scale_color_manual(name = "ONI:",
                     values = c("none" = "black", 
                                "El Nino" = "red",
                                "La Nina" = "blue")) + 
  scale_fill_manual(name = "Type:",
                     values = c("Year before" = "gray80",
                                "Year sampling" = "gray50"))  ; el.nino.la.nina.precipitation

ggsave(filename = "output/images/climate/precipitation.Puerto.ayora.el.nino.la.nina.png", 
       plot = el.nino.la.nina.precipitation, device = 'png', 
       width = 9, height = 4)
