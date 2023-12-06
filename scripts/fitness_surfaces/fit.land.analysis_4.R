# Description -------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created June 1, 2023

# Why: 
  # Plot a 3D fitness landscape for better visualization of the peaks. 

# Requires
  # Data from the fitness and adaptive landscape 

# NOTES:
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SETUP -------------------------------------------------------------------
## Preparation of variables and data  --------------------------------------
source("scripts/0.0_initialize.R")

## Load the data -----------------------------------------------------------
load('data/bird.data.RData', verbose = TRUE)

traits = c("avg.mbl","avg.mbd","avg.mbw")
view.traits = c("avg.mbl","avg.mbd")
zoom = 3
line.height = 1
linedat = data.frame(x =c(0,0),y = c(0,0),z = c(0,line.height))

linedat.all = bird.data %>%  
  group_by(sp2) %>%  
  summarise(x = mean(get(view.traits[1])),
            y = mean(get(view.traits[2]))) %>% 
  mutate(count = rep(2, nrow(.)),
         col = pal) %>%
  uncount(count) %>%  
  mutate(z = rep(c(0,1), nrow(.)/2))
# linedat.all$sp = linedat.all$sp2

linedat.for = linedat.all[linedat.all$sp2 %in% c("fortis small","fortis large"),]
linedat.ful = linedat.all[linedat.all$sp2 %in% c("fuliginosa"),]
linedat.mag = linedat.all[linedat.all$sp2 %in% c("magnirostris"),]
linedat.sca = linedat.all[linedat.all$sp2 %in% c("scandens"),]


fit.sum = readRDS("output/data.out/fitness.landscape.data/fitness.summary.RDS")

my.persp.data = readRDS("output/data.out/landscape_data/my.persp.data.RDS")
my.persp.data.plotly = my.persp.data
my.persp.data.plotly$z = t(my.persp.data.plotly$z)


# Data from the adative landscape
all.sp = readRDS(file = paste("output/data.out/adpt_land.sp/all.sp",ext.file,".RDS", sep = ""))
ada.data = cbind(expand.grid(all.sp$x, 
                             all.sp$y), 
                 as.vector(all.sp$z2)) 
colnames(ada.data) <- c("x", "y", "z")



# PLOT 3D Static ---------------------------------------------------------------------------------------------
## PLOT 3D: fitness landscape ----------------------------------------------------------------------------------

# Range of x and y for fitness and adaptive landscape 
range.x = range(c(my.persp.data$m1, unique(ada.data$x)))
range.y = range(c(my.persp.data$m2, unique(ada.data$y)))
range.z = range(c(my.persp.data$z, matrix(ada.data$z, nrow = 90, ncol = 90)), na.rm = TRUE)


png(filename = paste("output/images/landscape_plots/plot.3D.fit_land",ext.file,".png", sep = ""),
    width = 480,height = 480,units = "px", pointsize = 12)
par(mfrow = c(1,1), mar = c(1,2,1,2))

# Maximum of fitness or adaptive landscape 
max.z = ceiling(max(c(my.persp.data.plotly$z, ada.data$z), na.rm = TRUE))

rotation = 225
# max.z = 0.8
cat("Drawing fitness landscape in 3D", fill = TRUE)
plot3D::persp3D(x = my.persp.data$m1,
                y = my.persp.data$m2, 
                z = my.persp.data$z,
                xlab = "Beak length (mm)",
                ylab = "Beak depth (mm)",
                zlab = "Fitness",
                resfac = 1, # Resolution
                facets = TRUE, # True = surface, false = mesh  
                ticktype = "detailed", # Detailed = with numbers 
                xlim = range.x, 
                ylim = range.y, 
                # curtain = TRUE,
                # NAcol = 'yellow',
                # col.palette = heat.colors,
                # contour = list(side = c("zmin","z")), 
                colkey = list(side = 4, plot = TRUE, length = 1, width = 1, dist = 0, 
                              shift = 0, addlines = TRUE, col.clab = NULL, 
                              cex.clab = par("cex.lab"), side.clab = NULL, 
                              line.clab = NULL, adj.clab = NULL, font.clab = NULL),
                contour = list(side = c("z"), nlevels = 10), # where to add contour 
                # lighting = listt(ambient = 2),
                shade = 1, 
                alpha = 1, 
                box = T, # Draw the box
                bty = "u", #c("b", "b2", "f", "g", "bl", "bl2", "u", "n")
                # perspbox arguments
                  col.axis = "black", 
                  col.panel = "gray95", 
                  lwd.panel = 1,lwd.grid = 1,
                  col.grid = "grey", 
                # image = list(side = "zmin"),
                breaks = seq(0, (max(my.persp.data.plotly$z, na.rm = TRUE))+.1, by = .1), # Break the colours 
                zlim = c(0, max.z), # Maximum of z values on plot 
                # image = list(side = -1),
                theta = rotation, phi = 30) 

# Add line at mean beak traits to plot 
# for (spall in unique(linedat.all$sp2)) {
#   tmp.sub.line = linedat.all[linedat.all$sp2 %in% spall,]  
#   scatter3D(x = tmp.sub.line$x, 
#             y = tmp.sub.line$y, 
#             z = tmp.sub.line$z, 
#             lwd = 3, 
#             col = tmp.sub.line$col[1],
#             pch = 16, cex = 1.9,  add = TRUE,
#             colkey = FALSE, 
#             type = "l",
#             theta = rotation, phi = 30)
# }

# Add points at max fitness 
for (spall in unique(fit.sum$sp2)) {
  cat("Adding point for sp:", spall, fill = TRUE)
  tmp.sub.pt = fit.sum[fit.sum$sp2 %in% spall,]  
  points3D(
    x = tmp.sub.pt$beak.l, 
    y = tmp.sub.pt$beak.d, 
    z = tmp.sub.pt$z, # rep(.6, 5),
    # x = tmp.sub.pt$x, 
    # y = tmp.sub.pt$y, 
    # z = tmp.sub.pt$pd.z, # rep(.6, 5),
    # colkey = FALSE, 
    colkey = list(plot = FALSE),
           # colvar = fit.sum$pal,
           col = tmp.sub.pt$pal,
           pch = 19, cex = 2,  add = TRUE,
           type = "p",
           theta = rotation, phi = 30)
}

dev.off()

png(filename = paste("output/images/landscape_plots/plot.3D.fit_land_ggplot.cols_",ext.file,".png", sep = ""),
    width = 480,height = 480,units = "px", pointsize = 12)
par(mfrow = c(1,1), mar = c(1,2,1,2))

# Maximum of fitness or adaptive landscape 
max.z = ceiling(max(c(my.persp.data.plotly$z, ada.data$z), na.rm = TRUE))
rotation = 30
# max.z = 0.8
cat("Drawing fitness landscape in 3D", fill = TRUE)
plot3D::persp3D(x = my.persp.data$m1,
                y = my.persp.data$m2, 
                z = my.persp.data$z,
                xlab = "Beak length (mm)",
                ylab = "Beak depth (mm)",
                zlab = "Fitness",
                resfac = 1, # Resolution
                facets = TRUE, # True = surface, false = mesh  
                ticktype = "detailed", # Detailed = with numbers 
                xlim = range.x, 
                ylim = range.y, 
                # curtain = TRUE,
                # NAcol = 'yellow',
                # col.palette = heat.colors,
                col = hcl.colors(7, "YlOrRd", rev = TRUE, alpha = 1),
                # contour = list(side = c("zmin","z")), 
                colkey = list(side = 4, plot = TRUE, length = 1, width = 1, dist = 0, 
                              shift = 0, addlines = TRUE, col.clab = NULL, 
                              cex.clab = par("cex.lab"), side.clab = NULL, 
                              line.clab = NULL, adj.clab = NULL, font.clab = NULL),
                contour = list(side = c("z"), nlevels = 10), # where to add contour 
                # lighting = list(ambient = 2),
                # shade = 0.1,
                alpha = 1, 
                box = T, # Draw the box
                bty = "u", #c("b", "b2", "f", "g", "bl", "bl2", "u", "n")
                # perspbox arguments
                col.axis = "black", 
                col.panel = "gray95", 
                lwd.panel = 1,lwd.grid = 1,
                col.grid = "grey", 
                # image = list(side = "zmin"),
                breaks = seq(0, (max(my.persp.data.plotly$z, na.rm = TRUE))+.1, by = .1), # Break the colours 
                zlim = c(0, max.z), # Maximum of z values on plot 
                # image = list(side = -1),
                theta = rotation,
                phi = 30) 
dev.off()

png(filename = paste("output/images/landscape_plots/plot.3D.ada_land_ggplot.cols_",ext.file,".png", sep = ""),
    width = 480,height = 480,units = "px", pointsize = 12)
par(mfrow = c(1,1), mar = c(1,2,1,2))
# Maximum of fitness or adaptive landscape 
max.z = ceiling(max(c(my.persp.data.plotly$z, ada.data$z), na.rm = TRUE))
rotation = 30
# max.z = 0.8
cat("Drawing fitness landscape in 3D", fill = TRUE)
plot3D::persp3D(x = unique(ada.data$x),
                y = unique(ada.data$y), 
                z = matrix(ada.data$z, nrow = 90, ncol = 90),
                xlab = "Mean beak length (mm)",
                ylab = "Mean beak depth (mm)",
                zlab = "Fitness",
                resfac = 1, # Resolution
                facets = TRUE, # True = surface, false = mesh  
                ticktype = "detailed", # Detailed = with numbers 
                xlim = range.x, 
                ylim = range.y, 
                # curtain = TRUE,
                # NAcol = 'yellow',
                # col.palette = heat.colors,
                col = hcl.colors(7, "YlOrRd", rev = TRUE, alpha = 1),
                # contour = list(side = c("zmin","z")), 
                colkey = list(side = 4, plot = TRUE, length = 1, width = 1, dist = 0, 
                              shift = 0, addlines = TRUE, col.clab = NULL, 
                              cex.clab = par("cex.lab"), side.clab = NULL, 
                              line.clab = NULL, adj.clab = NULL, font.clab = NULL),
                contour = list(side = c("z"), nlevels = 10), # where to add contour 
                # lighting = list(ambient = 2),
                # shade = 0.1,
                alpha = 1, 
                box = T, # Draw the box
                bty = "u", #c("b", "b2", "f", "g", "bl", "bl2", "u", "n")
                # perspbox arguments
                col.axis = "black", 
                col.panel = "gray95", 
                lwd.panel = 1,lwd.grid = 1,
                col.grid = "grey", 
                # image = list(side = "zmin"),
                breaks = seq(0, (max(my.persp.data.plotly$z, na.rm = TRUE))+.1, by = .1), # Break the colours 
                zlim = c(0, max.z), # Maximum of z values on plot 
                # image = list(side = -1),
                theta = rotation,
                phi = 30) 
dev.off()


## PLOT 3D: only 1 legend --------------------------------------------------------------------------------------
png(filename = paste("output/images/landscape_plots/plot.3D.fit_land_oneleg",ext.file,".png", sep = ""),
    width = 480,height = 480,units = "px", pointsize = 12)
par(mar = c(1,2,1,2))
plot3D::persp3D(x = my.persp.data$m1,
                y = my.persp.data$m2,
                z = my.persp.data$z,
                xlab = "Beak length (mm)",
                ylab = "Beak depth (mm)",
                zlab = "Fitness",
                xlim = range.x, 
                ylim = range.y, 
                resfac = .1, # Resolution
                facets = TRUE, # True = surface, false = mesh  
                ticktype = "detailed", # Detailed = with numbers 
                colkey = list(side = 4, plot = TRUE, length = 1, width = 1, dist = 0, 
                              shift = 0, addlines = TRUE, col.clab = NULL, 
                              cex.clab = par("cex.lab"), side.clab = NULL, 
                              line.clab = NULL, adj.clab = NULL, font.clab = NULL),
                contour = list(side = c("z")), # where to add contour 
                shade = 1, 
                alpha = 1, 
                box = T, # Draw the box
                bty = "u", #c("b", "b2", "f", "g", "bl", "bl2", "u", "n")
                # perspbox arguments
                col.axis = "black", 
                col.panel = "gray95", 
                lwd.panel = 1,lwd.grid = 1,
                col.grid = "grey", 
                breaks = seq(0, (max(my.persp.data.plotly$z, na.rm = TRUE))+.1, by = .1), # Break the colours 
                zlim = c(0, max.z), # Maximum of z values on plot 
                theta = rotation, phi = 30) 
dev.off()


## PLOT 3D: adaptive landscape ----------------------------------------------------------------------------------
png(filename = paste("output/images/landscape_plots/plot.3D.adaptive.land",ext.file,".png", sep = ""),
    width = 480,height = 480,units = "px", pointsize = 12)
par(mfrow = c(1,1), mar = c(1,2,1,2))

plot3D::persp3D(x = unique(ada.data$x),
                y = unique(ada.data$y), 
                z = matrix(ada.data$z, nrow = 90, ncol = 90),
                xlab = "Beak length (mm)",
                ylab = "Beak depth (mm)",
                zlab = "Fitness",
                xlim = range.x, 
                ylim = range.y, 
                resfac = 3, # Resolution
                facets = TRUE, # True = surface, false = mesh  
                ticktype = "detailed", # Detailed = with numbers 
                colkey = list(side = 4, plot = TRUE, length = 1, width = 1, dist = 0, 
                              shift = 0, addlines = TRUE, col.clab = NULL, 
                              cex.clab = par("cex.lab"), side.clab = NULL, 
                              line.clab = NULL, adj.clab = NULL, font.clab = NULL),
                contour = list(side = c("z"), nlevels = 10), # where to add contour 
                shade = 1, 
                alpha = 1, 
                box = T, # Draw the box
                bty = "u", #c("b", "b2", "f", "g", "bl", "bl2", "u", "n")
                # perspbox arguments
                col.axis = "black", 
                col.panel = "gray95", 
                lwd.panel = 1,lwd.grid = 1,
                col.grid = "grey", 
                # Breaking colours based on fitness landscape (since adaptive landscape 
                # should be smoother, it's appropriate to use the fitness landscape as an upper 
                # bound for colour scale)
                breaks = seq(0, (max(my.persp.data.plotly$z, na.rm = TRUE))+.1, by = .1), # Break the colours 
                # breaks = seq(0, (max(ada.data$z, na.rm = TRUE))+.1, by = .1), # Break the colours 
                zlim = c(0, max.z), # Maximum of z values on plot 
                theta = rotation, phi = 30) 



dev.off()

# PLOT 3D interactive ---------------------------------------------------------------------------------------
## PLOTLY: Fit-Land Interactive --------------------------------------------------------
# https://plotly.com/r/reference/#scatter-legendgrouptitle

# PLOTLY  
zoom = 2
# Scale the fitness for better visuals 
bird.data$mxcpois.sca = bird.data$mxcpois/max(bird.data$mxcpois)
# Colour palette
n.col = 10
color_plate.dens2 = hcl.colors(n.col, "YlOrRd", rev = TRUE, alpha = 1)
seq.z = seq(min(my.persp.data.plotly$z, na.rm = TRUE), 
            max(my.persp.data.plotly$z, na.rm = TRUE), 
            length.out = n.col)
z.col = data.frame(as.character(seq.z), 
                   as.character(color_plate.dens2))
z.col.list <- split(z.col, 
                    seq(nrow(z.col)))
z.col.list = unname(z.col.list)

pg = 
  plotly::plot_ly(width = 800, height = 800) %>% 
  # Add the points 
  add_trace(data = bird.data,
            x = as.formula(paste("~",traits[1])),
            y = as.formula(paste("~",traits[2])),
            showlegend = TRUE,
            legendgroup = 'group2',
            legendgrouptitle = list(text = '<b>Individuals</b>', font = list(size = 14)), # showlegend=FALSE,
            opacity = .28,
            type="scatter3d",
            mode = "markers",
            inherit = FALSE,
            text = ~paste('Nb obs: ', mxcpois,'\nID: ', BANDFINAL),
            colors = pal,
            color = ~sp2,
            z = ~mxcpois.sca) %>%
  layout(title = '<b> Fitness landscape ground finches </b>', 
         margin = list(t = 80),
         # legend = list(title=list(text='<b> Species of finches </b>')),
         scene=list(
           xaxis = list(title = 'Avg beak median length',
                        autorange = TRUE),
           yaxis = list(title = 'Avg beak median depth'),
           # autorange = "reversed"),
           zaxis=list(title = 'Scaled apparent survival',
                      nticks = 10,
                      range = c(0,1)),# Number of years observed in population
           camera = list(eye = list(x = cos(1/1.1*pi)*zoom, 
                                    y = sin(1/2*pi)*zoom, 
                                    z= 2.00)))) %>% 
  # Add the actual surface 
  add_surface(data = my.persp.data.plotly,
              x = ~m1,
              y = ~m2, 
              z = ~z,
              type="surface",
              contours = list(
                x = list(show = TRUE, width = 1.5, start = -3, end = 2, size = 0.20, color = gray(0.2, 1)),
                y = list(show = TRUE, width = 1.5, start = -3, end = 2, size = 0.20, color = gray(0.2, 1)),
                z = list(show = TRUE, width = 2.0, start =  0, end = 1, size = 0.05, color = gray(0.5, 1))), 
              opacity = 1,
              colorscale = "Viridis",
              # colorscale = list(c(0, 1), c("yellow","red3")),
              # colorscale = list(c(0, 20), color_plate.dens2[c(1,10)]),
              # colorscale = list(c(min(z),"rgb(107,184,214)"),c(max(z),"rgb(0,90,124)")),
              colorbar = list(title = "<b>Apparent survival </b>")) %>%  
  # add_surface(x = my.persp.data$m2,
  #             y = my.persp.data$m1,
  #             z = my.persp.data$z.se.u,
  #             type="surface",
  #             showscale = FALSE,
  #             inherit = FALSE,
  #             opacity = .5,
  #             colorbar = list(title = "<b>2 SE Individuals</b>"),
  #             colorscale = "red") %>%
  # add_surface(x = my.persp.data$m2,
  #             y = my.persp.data$m1,
#             z = my.persp.data$z.se.l,
#             type="surface",
#             showscale = FALSE,
#             inherit = FALSE,
#             opacity = .5,
#             colorbar = list(title = "<b>2 SE</b>"),
#             colorscale = "red") %>%
# add_trace(data = w.data,
#           x = ~avg.trait1,
#           y = ~avg.trait2,
#           # showlegend=FALSE,
#           marker = list(color = "green",opacity = 1, size = 14),
#           type="scatter3d",
#           mode = "markers",
#           showlegend=FALSE,
#           inherit = FALSE,
#           z = ~w.bar) %>% 
  add_trace(data=linedat.for[3:4,], x=~x, y=~y, z=~z,
            type="scatter3d", mode="lines", 
            legendgroup = 'group1',
            showlegend=TRUE,
            legendgrouptitle = list(text = '<b>Average of survivors </b>', font = list(size = 14)),
            line = list(color = pal[2], 
                        width = 14),name = "<i>G. fortis</i> small")%>% 
  add_trace(data=linedat.for[1:2,], x=~x, y=~y, z=~z,
            showlegend=TRUE,
            type="scatter3d", mode="lines", legendgroup = 'group1',
            line = list(color = pal[1], 
                        width = 14),name = "<i>G. fortis</i> large")%>% 
  add_trace(data=linedat.ful, x=~x, y=~y, z=~z, 
            showlegend=TRUE,
            type="scatter3d", mode="lines", legendgroup = 'group1',
            line = list(color = pal[3], 
                        width = 14),name = "<i>G. fuliginosa</i>")%>% 
  add_trace(data=linedat.mag, x=~x, y=~y, z=~z, 
            showlegend=TRUE,
            type="scatter3d", mode="lines", legendgroup = 'group1',
            line = list(color = pal[4], 
                        width = 14),name = "<i>G. magnirostris</i>")%>% 
  add_trace(data=linedat.sca, x=~x, y=~y, z=~z, 
            showlegend=TRUE,
            type="scatter3d", mode="lines", legendgroup = 'group1',
            line = list(color = pal[5], 
                        width = 14),name = "<i>G. scandens</i>");pg

## save the output 
dir.create("output/fitland_interactive", showWarnings = FALSE)
htmlwidgets::saveWidget(pg, file = "output/fitland_interactive/interactive.fl3d.html")
saveRDS(pg,file = "output/fitland_interactive/plotlygraph.rds")
kaleido(pg, file = "output/fitland_interactive/fit.land.3d.png")





## PLOTLY: Adaptive-Land Interactive --------------------------------------------------------
surface.finch  = list(x = unique(ada.data$x),
                      y = unique(ada.data$y), 
                      z = matrix(ada.data$z, nrow = 90, ncol = 90))

# PLOTLY  
zoom = 2
# Colour palette
n.col = 10
color_plate.dens2 = hcl.colors(n.col, "YlOrRd", rev = TRUE, alpha = 1)
seq.z = seq(min(surface.finch$z, na.rm = TRUE), 
            max(surface.finch$z, na.rm = TRUE), 
            length.out = n.col)
z.col = data.frame(as.character(seq.z), 
                   as.character(color_plate.dens2))
z.col.list <- split(z.col, 
                    seq(nrow(z.col)))
z.col.list = unname(z.col.list)

pg.adapt = 
  plotly::plot_ly(width = 800, height = 800) %>% 
  # Add the points 
  add_trace(data = bird.data,
            x = as.formula(paste("~",traits[1])),
            y = as.formula(paste("~",traits[2])),
            showlegend = TRUE,
            legendgroup = 'group2',
            legendgrouptitle = list(text = '<b>Individuals</b>', font = list(size = 14)), # showlegend=FALSE,
            opacity = .28,
            type="scatter3d",
            mode = "markers",
            inherit = FALSE,
            text = ~paste('Nb obs: ', mxcpois,'\nID: ', BANDFINAL),
            colors = pal,
            color = ~sp2,
            z = ~mxcpois.sca) %>%
  layout(title = '<b> Adaptive landscape ground finches </b>', 
         margin = list(t = 80),
         # legend = list(title=list(text='<b> Species of finches </b>')),
         scene=list(
           xaxis = list(title = 'Avg beak median length',
                        autorange = TRUE),
           yaxis = list(title = 'Avg beak median depth'),
           # autorange = "reversed"),
           zaxis=list(title = 'Scaled apparent survival',
                      nticks = 10,
                      range = c(0,1)),# Number of years observed in population
           camera = list(eye = list(x = cos(1/1.1*pi)*zoom, 
                                    y = sin(1/2*pi)*zoom, 
                                    z= 2.00)))) %>% 
  # Add the actual surface 
  add_surface(data = surface.finch,
              x = ~x,
              y = ~y, 
              z = ~t(z),
              type="surface",
              contours = list(
                x = list(show = TRUE, width = 1.5, start = -3, end = 2, size = 0.20, color = gray(0.2, 1)),
                y = list(show = TRUE, width = 1.5, start = -3, end = 2, size = 0.20, color = gray(0.2, 1)),
                z = list(show = TRUE, width = 2.0, start =  0, end = 1, size = 0.05, color = gray(0.5, 1))), 
              opacity = 1,
              colorscale = "Viridis",
              colorbar = list(title = "<b>Apparent survival </b>")) %>% 
  
  add_trace(data=linedat.for[3:4,], x=~x, y=~y, z=~z,
            type="scatter3d", mode="lines", 
            legendgroup = 'group1',
            showlegend=TRUE,
            legendgrouptitle = list(text = '<b>Average of survivors </b>', font = list(size = 14)),
            line = list(color = pal[2], 
                        width = 14),name = "<i>G. fortis</i> small")%>% 
  add_trace(data=linedat.for[1:2,], x=~x, y=~y, z=~z,
            showlegend=TRUE,
            type="scatter3d", mode="lines", legendgroup = 'group1',
            line = list(color = pal[1], 
                        width = 14),name = "<i>G. fortis</i> large")%>% 
  add_trace(data=linedat.ful, x=~x, y=~y, z=~z, 
            showlegend=TRUE,
            type="scatter3d", mode="lines", legendgroup = 'group1',
            line = list(color = pal[3], 
                        width = 14),name = "<i>G. fuliginosa</i>")%>% 
  add_trace(data=linedat.mag, x=~x, y=~y, z=~z, 
            showlegend=TRUE,
            type="scatter3d", mode="lines", legendgroup = 'group1',
            line = list(color = pal[4], 
                        width = 14),name = "<i>G. magnirostris</i>")%>% 
  add_trace(data=linedat.sca, x=~x, y=~y, z=~z, 
            showlegend=TRUE,
            type="scatter3d", mode="lines", legendgroup = 'group1',
            line = list(color = pal[5], 
                        width = 14),name = "<i>G. scandens</i>");pg.adapt

# Save -------------------------------------------------------------------------------------------------------
## save HTML --------------------------------------------------------------------------------------------------
dir.create("output/fitland_interactive", showWarnings = FALSE)
htmlwidgets::saveWidget(pg.adapt, file = "output/fitland_interactive/interactive.al3d.html")

## Save data --------------------------------------------------------------------------------------------------
saveRDS(pg,file = "output/fitland_interactive/plotlygraph_adapt.land.rds")


## Save static  -----------------------------------------------------------------------------------------------
# library(reticulate)
# install kaleido in python from R 
# reticulate::py_run_string("import sys")
# system(command = "python3 -m pip install plotly")
# system(command = "python3 -m pip install kaleido")
# save_image(pg, "output/fitland_interactive/adapt.land.3d.png")
