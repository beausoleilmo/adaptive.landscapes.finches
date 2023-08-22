# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created Tuesday, October 13, 2022 
# Why: 
  # Galápagos maps 
  # El Garrapatero sites from Garmin GPS tracks 
  # To be modified and improved to match our sampling design
  # Use the coordinates of the nets to find the sampling area. 
# Requires 
# NOTES: 
  # converter: 
    # https://epsg.io/transform#s_srs=4326&t_srs=32715&x=-90.2221389&y=0.6896939
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Preparation script ------------------------------------------------------
source("scripts/0.0_initialize.R")

# Load data ---------------------------------------------------------------
# Location of nets 
df2.1 = read.csv("data/maps/Galapagos.EG.nets.csv")

# EG 
CRS.gal = "+init=epsg:32715"
sputm2 = sf::st_as_sf(df2.1, coords = c("lon","lat"))
# Setting Coord. ref sys. 
st_crs(sputm2) <- 4326
spgeo2 = st_transform(sputm2, CRS(CRS.gal))
spa.finch = as_Spatial(spgeo2$geometry)


# Estimating sampling area ------------------------------------------------
## Minimum Convex Polygon Estimator ----------------------------------------
nets.mcp = mcp(spa.finch, 
               percent=100, 
               unin = c("m"),
               unout = c("ha"))

# Area size of MCP sampling area
mcp.area.nets = st_area(st_as_sf(nets.mcp))
units::set_units(mcp.area.nets, value = "ha")
units::set_units(mcp.area.nets, value = "km2")
round(units::set_units(mcp.area.nets, value = "km2") , 2)


## Estimation of Kernel Home-Range -----------------------------------------
kudl <- kernelUD(spa.finch, 
                 h="href", 
                 same4all = TRUE, 
                 grid = 500)#, g

vud_points <- getvolumeUD(kudl)

# Gives the Home ranges with a percentage of the distribution
vud_lines <- getverticeshr(kudl,percent = c(95),unin = c("m"), unout = "m2")
vud.sf = st_as_sf(vud_lines)

# Area size of Kernel sampling around the mist nets at El Garrapatero 
eg.area = st_area(vud.sf)
units::set_units(eg.area, value = "ha")
round(units::set_units(eg.area, value = "km2") , 3)

# Load map background  ----------------------------------------------------
map.back <- st_read("data/maps/map.back.gpkg")
label.ggplot = sort(unique(map.back$Nombre))

# Get colours 
col.geom = c("white","gray90","white","#969696")
map.back = map.back[!(map.back$Nombre %in% c("Pacific Ocean","Inland water")),]

# GGplot map sampling area --------------------------------------------------------------
sampling_area_EG.png = 
  ggplot(spgeo2, aes(geometry = geometry)) + 
  # Add map background 
  geom_sf(data = map.back, aes(geometry = geom, 
                               fill = factor(Nombre)), 
          color = "black", 
          inherit.aes = FALSE , 
          linewidth = 0.05)+
  scale_fill_manual(labels = map.back$Nombre,
                    values = c("gray90", "gray100", "gray80", "gray60"),
                    name = "Habitat type") +
  # Adding the KernelUD 
  geom_sf(data = vud.sf, aes(geometry = geometry), fill = alpha("black", alpha = .1), 
          linewidth = 0.7) +
  # Adding points
  geom_sf(data = spgeo2, aes(geometry = geometry)) +
  # Set x axis ticks, but not working 
  scale_x_continuous(breaks = seq(-90.226, -90.218, length.out = 3)) +
  # Determine Coordinate limits 
  coord_sf(xlim = c(808700, 809600),  
           ylim = c(9923400,9924300)
  ) + 
  # Make the raph prettier
  theme(axis.ticks = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.position="none", # REMOVES THE LEGEND 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12, face='bold'),
        plot.tag = element_text(size = 26, face="bold"),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black', face='bold',size = 15)) + 
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         height = unit(1.0, "cm"),
                         width = unit(0.8, "cm"),
                         pad_x = unit(0.25, "cm"),
                         pad_y = unit(0.55, "cm"),
                         style = north_arrow_orienteering) + 
  annotation_scale(aes(location = "br"),
                   height = unit(0.15, "cm"), 
                   text_cex = 1.2, line_width = .9) + 
  labs(tag = "B") +
  guides(fill = guide_legend(ncol = 1));sampling_area_EG.png

# Load GGplot maps_galapagos.R --------------------------------------------
# GPS Position locations --------------------------------------------------
## El Garrapatero ----------------------------------------------------------
EGpos  = c(809206.09,9923677.77) # UTM NAD83 Zone15N
EGpos_sfc = st_sfc(st_point(EGpos))
EGpos_sf = st_sf(a = 1, geom = EGpos_sfc)
st_crs(EGpos_sf) = 32715

## Weather station ---------------------------------------------------------
WSpos  = c( 800235.2628102701, 9917704.520241257) # UTM NAD83 Zone15N
WSpos_sfc = st_sfc(st_point(WSpos))
WSpos_sf = st_sf(a = 1, geom = WSpos_sfc)
st_crs(WSpos_sf) = 32715

# Distance between 2 points 
dist.EG.weather.station = st_distance(EGpos_sf, WSpos_sf)
units::set_units(dist.EG.weather.station, value = "km")

# Read layers from GeoPackage ---------------------------------------------
df_gal <- st_read(dsn = "data/maps/galapagos.gis.gpkg",  layer = "gal_gis")
df_roads <- st_read(dsn = "data/maps/galapagos.gis.gpkg",  layer = "Roads")

# Subset the map
sub.gal = df_gal[df_gal$nam %in% c("ISLA SANTA CRUZ","ISLA DAPHNE MAYOR",
                                   "ISLA DAPHNE MINOR","ISLA BALTRA","ISLA SEYMOUR"),]
isc.gal = df_gal[df_gal$nam %in% c("ISLA SANTA CRUZ"),]

col.geom = c("#FFFFD9", "#C6DBEF", "#C7E9B4", "#eefff4", "#6BAED6", "#969696")

land.colour = col.geom[1]
water.colour = col.geom[2]
topo.colour = alpha(rgb(177/256,162/256,131/256),.5)
land.colour = "white"
water.colour = "white"
topo.colour = NA

# GGplot /PUB\ (inset) galapagos.isla.archi ----------------------------------------------------
# All islands of the archipelago 
galapagos.isla.archi = ggplot() + 
  # Add map background 
  geom_sf(data = df_gal, aes(geometry = geom), color = "black", fill = land.colour) +
  geom_sf(data = isc.gal, aes(geometry = geom), color = "black", fill = 'black') +
  scale_x_continuous(breaks = seq(-92.500, -89.500, length.out = 3)) +
  scale_y_continuous(breaks = seq(1.5, -1.5, length.out = 3)) +
  # Make the raph prettier
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = alpha("white", 0.5)),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(colour = "grey98", linetype = "dotted"),
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12, face='bold'),
        plot.tag = element_text(size = 26, face="bold"),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black', face='bold', size = 15)) + 
  guides(fill = guide_legend(ncol = 1));galapagos.isla.archi


# GGplot /PUB\  galapagos.isla ----------------------------------------------------------

# Draw the map
galapagos.isla = ggplot() + 
  # Add map background 
  geom_sf(data = df_gal, aes(geometry = geom),  
          linewidth = 0.5,
          color = "black", fill = land.colour) +
  geom_sf(data = EGpos_sf, aes(geometry = geom), col = "black",size = 6)+
  scale_x_continuous(breaks = seq(-90.500, -90.200, length.out = 4)) +
  scale_y_continuous(breaks = seq(-0.4, -0.75, length.out = 6)) +
  # Determine Coordinate limits 
  coord_sf(xlim = c(771339,818000),  ylim = c(9912554,9956500)) +
  # Make the raph prettier
  theme(axis.ticks = element_line(colour = "black"), 
        axis.text = element_text(colour = "black",size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = water.colour), 
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(colour = "grey98", linetype = "dotted"),
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12, face='bold'),
        plot.tag = element_text(size = 26, face="bold"),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black', face='bold',size = 15)) + 
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         height = unit(1.0, "cm"),
                         width = unit(0.8, "cm"),
                         pad_x = unit(0.25, "cm"),
                         pad_y = unit(0.55, "cm"),
                         style = north_arrow_orienteering) + 
  annotation_scale(aes(location = "br"),height = unit(0.15, "cm"), 
                   text_cex = 1.2, line_width = .9) + 
  annotate(geom = "text",
           x = st_coordinates(EGpos_sf$geom)[1]-c(9*1000),
           y = st_coordinates(EGpos_sf$geom)[2],
           colour = "black", fontface =1, size = 7,
           label = "El Garrapatero")+
  annotate(geom = "text",
           x = 793714.8,
           y = 9930106.0 + c(3*1000),
           colour = alpha("black", .3), fontface =2, size = 12,
           label = "Santa Cruz\n Island")+
  labs(tag = "A") +
  guides(fill = guide_legend(ncol = 1));galapagos.isla


# GGplots definition ------------------------------------------------------
plot1 = sampling_area_EG.png # Sampling area in El Garrapatero (with nets locations)
plot2 = galapagos.isla       # Santa Cruz Island panel A 
plot3 = galapagos.isla.archi # Whole archipelago

# /PUB\ GGplot maps Galapagos, sampling sites with inset  ----------------------
# See inset map 
# https://geocompr.github.io/post/2019/ggplot2-inset-maps/
gg_inset_map2 = ggdraw() +
  draw_plot(plot2) +
  draw_plot(plot3, 
            x = 0.118, 
            y = 0.6515, width = 0.30, height = 0.30)

arranged.plots1 = grid.arrange(gg_inset_map2, plot1, nrow =1)
ggsave("output/images/maps/santa_Cruz_and_sampling_area_EG.png",
       plot = arranged.plots1,
       device = "png",
       width = 16, height = 7, units = "in", dpi = 300)
