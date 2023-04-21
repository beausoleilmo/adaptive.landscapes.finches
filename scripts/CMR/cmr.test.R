# Description  ------------------------------------------------------------
##########################################################################################
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created ??, 2022
# Why: 
# Requires 
# NOTES: 
# https://cran.r-project.org/web/packages/marked/marked.pdf
# https://stats.stackexchange.com/questions/38904/what-does-return-rate-mean-in-population-ecology-models/38930#38930
# https://cran.r-project.org/web/packages/marked/vignettes/markedVignette.html
##########################################################################################

# p is the recapture rate (encoutner probability)
# phi is the survival probability 
# return rate = p*phi

# Load libraries ----------------------------------------------------------
source("scripts/0.0_initialize.R")

# Load data ---------------------------------------------------------------
load('data/bird.data.RData', verbose = TRUE)

cmr.dat = bird.data
cmr.in = cmr.dat[,c("ch","sp2","avg.mbl", "avg.mbw", "avg.mbd", "Sex0")]
cmr.in$ch = cmr.in$ch$ch
cmr.in=as.data.frame(cmr.in)
cmr.in$sp2 = as.factor(cmr.in$sp2)

finches.proc = process.data(cmr.in,
                            model='cjs',
                            begin.time=2003)

cmr.in.proc = process.data(cmr.in, model='cjs', begin.time=2003)
cmr.in.dll = make.design.data(cmr.in.proc)
phi.list = list(formula=~sp2)
p.list = list(formula=~sp2+time)

mod.finch = crm(data = cmr.in.proc,
                ddl = cmr.in.dll, 
                model.parameters = list(Phi = phi.list, p = p.list), 
                accumulate = FALSE, hessian = TRUE)

mod.finch
mod.finch$results$beta

proportions.parameters = compute_real(mod.finch)

plot(proportions.parameters$p$estimate~proportions.parameters$p$sp2)
points(proportions.parameters$p$sp2, 
       proportions.parameters$p$estimate, 
       pch = 19, 
       col = viridis(n = length(levels(as.factor(proportions.parameters$p$time))))[as.factor(proportions.parameters$p$time)])

mean.pts = proportions.parameters$p %>% 
  group_by(sp2) %>% 
  summarise(mean.est = mean(estimate),
            median.est = median(estimate)) 
points(mean.pts$sp2, mean.pts$mean.est, pch =19, col = "red")
points(mean.pts$sp2, mean.pts$median.est, pch =19, col = "red")

# Reported values 
round(range(mean.pts$median.est)*100, 2)
round(mean(mean.pts$median.est)*100, 2)

# Create design data with static and time-varying covariates
finches.ddl=make.design.data(finches.proc)
mod.Phit.pt=crm(finches.proc,finches.ddl,
                model.parameters=list(Phi=list(formula=~time),
                                      p=list(formula=~time)))

mod.Phisex.pdot=crm(finches.proc,finches.ddl,groups="sp2",
                    model.parameters=list(Phi=list(formula=~sp2),
                                          p=list(formula=~1)))
