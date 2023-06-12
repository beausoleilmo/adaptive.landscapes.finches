# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Analysis of the dynamics of Darwin's finches' fitness landscapes  
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created April 13, 2022
# Why: 
  # Model exploration 
# Requires 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Preparation of variables and data  --------------------------------------
source('scripts/0.0_initialize.R')
load('data/bird.data.RData', verbose = TRUE)

# Create directory where the models will be recorded 
dir.create("output/model.out", showWarnings = FALSE)

traits = c("avg.mbl","avg.mbd")

# MODELS ------------------------------------------------------------------
k = c(4,4,19) 
## GAM s (PC scores) ---------------------------------------------------
gam3.pc.s <- gam(mxcpois ~ s(pc1.new, bs = "tp", k = k[1]) + 
                   s(pc2.new, bs = "tp", k = k[2])+
                   s(pc1.new, pc2.new, bs = "tp", k = k[3]),
                 select=TRUE, 
                 data = bird.data, family = poisson(link = "log"))
vis.gam(gam3.pc.s, type = "response",theta=-80, phi=30,zlim = c(0,1.2),n.grid = 100)

## GAM s (traits) ------------------------------------------------------------
### EG ----------------------------------------------------------------------
if (site.check == "El Garrapatero") {
  k = 4
  gam3.p <- gam(mxcpois~s(avg.mbl, bs = "tp", k = k) + # k = 5
                  s(avg.mbd, bs = "tp", k = k)+ # k = 5
                  s(avg.mbl,avg.mbd, bs = "tp", k = 27), # k = 25
                data=bird.data,
                select=TRUE, 
                method="REML", # see https://osf.io/wgc4f/wiki/mgcv:%20model%20selection/
                family = poisson(link = "log")) 
  par(mfrow = c(2,2))
  # gam.check(gam3.p, rep = 500)
  par(mfrow = c(1,1))
  vis.gam(gam3.p, type = "response",theta=-80, phi=30,zlim = c(0,1.2),n.grid = 100)
  summary(gam3.p)
}


# Select on/off ----------------------------------------------------------------------------------------------
# Demonstration of the effect of removing the 'select' argument from the model. 
gam3.p.selectoff <- gam(mxcpois~s(avg.mbl, bs = "tp", k = k) + # k = 5
                s(avg.mbd, bs = "tp", k = k)+ # k = 5
                s(avg.mbl,avg.mbd, bs = "tp", k = 27), # k = 25
              data=bird.data,
              select=FALSE, 
              method="REML", # see https://osf.io/wgc4f/wiki/mgcv:%20model%20selection/
              family = poisson(link = "log")) 

# PNG effect of select on landscape --------------------------------------------------------------------------
png(filename = paste("output/images/landscape_plots/fitplot.3D_select",ext.file,".png", sep = ""),
    width = 780,height = 480,units = "px", pointsize = 12)
par(mfrow = c(1,2), mar = c(.1,1,1,1))
vis.gam(gam3.p,main = "with select = TRUE", type = "response",theta=-80, phi=30,zlim = c(0,1.2),n.grid = 70)
vis.gam(gam3.p.selectoff,main = "with select = FALSE", type = "response",theta=-80, phi=30,zlim = c(0,1.2),n.grid = 70)
dev.off()
par(mfrow = c(1,1))


# Diagnostics -------------------------------------------------------------
## GAM model comparison ----------------------------------------------------
gam3.pint <- gam(mxcpois~1,
                 data=bird.data,
                 select=TRUE, 
                 method="REML", 
                 family = poisson(link = "log"))
gam3.p1s <- gam(mxcpois~s(avg.mbl, bs = "tp", k = 4) +
                  s(avg.mbd, bs = "tp", k = 4),
                data=bird.data,
                select=TRUE, 
                method="REML", 
                family = poisson(link = "log"))

# Models with traits without smoothing
# For comparison with the GAM with smoothing
gam3.p1l <- gam(mxcpois~avg.mbl,
                data=bird.data,
                select=TRUE, 
                method="REML", 
                family = poisson(link = "log"))
gam3.p1ld <- gam(mxcpois~avg.mbl + avg.mbd,
                data=bird.data,
                select=TRUE, 
                method="REML", 
                family = poisson(link = "log"))
gam3.p1ldi <- gam(mxcpois~avg.mbl * avg.mbd,
                data=bird.data,
                select=TRUE, 
                method="REML", 
                family = poisson(link = "log"))

gam3.pqpois <- gam(mxcpois~s(avg.mbl, bs = "tp", k = 4) +
                     s(avg.mbd, bs = "tp", k = 4)+
                     s(avg.mbl,avg.mbd, bs = "tp", k = 27),
                   data=bird.data,
                   select=TRUE, 
                   method="REML", 
                   family = quasipoisson(link = "log"))

## AIC ---------------------------------------------------------------------
AIC(gam3.pint) - AIC(gam3.p)

AIC(gam3.pint)
AIC(gam3.p1s)
AIC(gam3.p1l)
AIC(gam3.p1ld)
AIC(gam3.p1ldi)
AIC(gam3.p)
AIC(gam3.pc.s)
AIC(gam3.pqpois)

# Absolute delta AIC
abs.diff.aic = abs(diff(c(AIC(gam3.pint),
           AIC(gam3.p1l), 
           AIC(gam3.p1ld), 
           AIC(gam3.p1ldi), 
           AIC(gam3.p1s), 
           AIC(gam3.p))));abs.diff.aic

## Deviance ----------------------------------------------------------------
deviance(gam3.pint)
deviance(gam3.p1l)
deviance(gam3.p1ld)
deviance(gam3.p1ldi)
deviance(gam3.p1s)
deviance(gam3.p)
deviance(gam3.pqpois)

abs(diff(c(deviance(gam3.pint),
           deviance(gam3.p1l),
           deviance(gam3.p1ld),
           deviance(gam3.p1ldi),
           deviance(gam3.p1s),
           deviance(gam3.p))))

deviance(gam3.pint) - deviance(gam3.p)
AIC(gam3.pint) - AIC(gam3.p)

# Anova table comparing different spline models 
anov.tab = anova(gam3.pint, 
      gam3.p1l, 
      gam3.p1ld, 
      gam3.p1ldi, 
      gam3.p1s, 
      gam3.p, 
      test = "Chisq");anov.tab
# Clean the table 
models = strsplit(x = attr(anov.tab, "heading")[2], split = "Model ")[[1]]
models = models[nchar(models)>0]
models = gsub(pattern = "\n", replacement = "", models)
models = gsub(pattern = "avg.", replacement = "", models)
models = gsub(pattern = "mxcpois", replacement = "y", models)

# Make simple model vector to append to the anova table 
models.simple = c("y ~ 1", 
  "y ~ bl", 
  "y ~ bl + bd", 
  "y ~ bl * bd", 
  "y ~ s(bl) + s(bd)", 
  "y ~ s(bl) + s(bd) + s(bl, bd)")
anov.tab.df = as.data.frame(anov.tab)
# Add delta AIC
anov.tab.df$abs.diff.aic = c(NA,abs.diff.aic)
# Combine all columns 
table.models = cbind(models.simple, 
                     family = "poisson(log)",  
                     round(anov.tab.df, 1))
# Change the p-value column 
table.models$`Pr(>Chi)` = ifelse(anov.tab.df$`Pr(>Chi)` < 0.00001, 
                                 yes = "<<0.01", 
                                 no = round(anov.tab.df$`Pr(>Chi)`, 4))
# Make column names for publication 
names(table.models) = c("Models†", "family(link)", "Residuals degrees of freedom", 
                        "Residual deviance", "Degrees of freedom", "ΔDeviance","P-value", "ΔAIC")
# Export table for publication 
write.table(x = table.models, 
            file = "output/model.out/GAM_model.selection.txt", 
            sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# With pca 
anova(gam3.pint, 
      gam3.p1l, 
      gam3.p1ld, 
      gam3.p1ldi, 
      gam3.p1s, 
      gam3.p, 
      gam3.pc.s, 
      test = "Chisq")

## Check GAM ---------------------------------------------------------------
k.check(gam3.p, subsample=5000, n.rep=400)
gam.check(gam3.p)

# See https://stats.stackexchange.com/questions/569108/gam-setting-k-values-before-or-after-testing-the-different-models
gam.gaus.bl = gam(residuals.gam(gam3.p, type = "deviance")~s(bird.data$avg.mbl, bs = "tp", k = 8), 
                  family = gaussian(), select = TRUE, method="REML")
summary(gam.gaus.bl)
gam.gaus.bd = gam(residuals.gam(gam3.p, type = "deviance")~s(bird.data$avg.mbd, bs = "tp", k = 8), 
                  family = gaussian(), select = TRUE, method="REML")
summary(gam.gaus.bd)

# SAVE --------------------------------------------------------------------
## Save model gam3.p -------------------------------------------------------
saveRDS(object = gam3.p, file = paste("output/model.out/poisson.spline.model",ext.file,".RDS", sep = ""))

## Save model gam3.pc -------------------------------------------------------
saveRDS(object = gam3.pc.s,    file = "output/model.out/poisson.spline.model_pca.RDS")
