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
AIC(gam3.p)
AIC(gam3.pc.s)
AIC(gam3.pqpois)

# Absolute delta AIC
abs(diff(c(AIC(gam3.pint),AIC(gam3.p1s), AIC(gam3.p))))

## Deviance ----------------------------------------------------------------
deviance(gam3.pint)
deviance(gam3.p1s)
deviance(gam3.p)
deviance(gam3.pqpois)

abs(diff(c(deviance(gam3.pint),deviance(gam3.p1s),deviance(gam3.p))))

deviance(gam3.pint) - deviance(gam3.p)
AIC(gam3.pint) - AIC(gam3.p)

# Anova table comparing different spline models 
anova(gam3.pint, gam3.p1s, gam3.p, test = "Chisq")

# With pca 
anova(gam3.pint, gam3.p1s, gam3.p, gam3.pc.s, test = "Chisq")

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
