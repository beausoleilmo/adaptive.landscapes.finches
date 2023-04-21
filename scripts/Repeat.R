# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created April 13, 2022
# Why: 
  # Repeatability of remeasured individuals 
# Requires 
# NOTES: 
  # https://cran.r-project.org/web/packages/rptR/vignettes/rptR.html
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Prepare script ----------------------------------------------------------
source("scripts/0.0_initialize.R")

# Load the data -----------------------------------------------------------
load('data/bird.data.RData', verbose = TRUE)

# See who measured 
unique(bird.data$By.)
table(bird.data$By.)

# Get subset of data relevant for the analysis 
sub.dat = bird.data %>% 
  dplyr::filter(Sex0 %in% c("", "f", "f;j", "f?", "f?;j", "fledgling", 
                            "j", "m", "m;j", "m?", "nestling"), 
                Year %in% c(2003:2020), # subset years
                !(BANDFINAL %in% ""), # Remove empty bands 
                !(BANDFINAL %in% c("LF1441", "LF1442", "LF1443", "LF1444", 
                                   "LF1451", "LF5192", "SH765_L")), # Remove some other weird ID
                Site %in% c("El Garrapatero"), # Keep only 1 site
                Species1 %in% c("fortis", # Keep ground finches and species limits 
                                "fuliginosa", "magnirostris","scandens")) 

# Make longer format: as we have taken Multiple measurements (usually 3) for each individual 
# for each beak trait (Length, depth, width), we can use that data to measure repeatability. 
# It might not be 'independent' doing this
finch.recp = pivot_longer(sub.dat,
                          cols = starts_with("Beak."), 
                          names_to = "beak", 
                          values_to = "measure")
# Check the new data
head(finch.recp[,c("BANDFINAL", "beak","measure")])

# Get repeatability from all the data (three measurements taken)
reap.finch <- function(data, species, 
                       trait.search, 
                       nboot = 1000, nperm = 1000, plot=T) {
  
  data.beak.meas = data[grep(pattern = trait.search,data$beak), ]
  tmp.dat = data.beak.meas[data.beak.meas[,"Species1"] == species 
                           ,]
  my.form = formula(paste0("measure"," ~ (1|BANDFINAL)"))
  tmp=rpt(my.form,
          parallel = TRUE, # Won't show progress 
          grname = c("BANDFINAL"),
          data = tmp.dat,
          datatype ="Gaussian", 
          nboot = nboot, 
          npermut=nperm)
  if(plot){
    plot(tmp)
  }
  return(list(repeatability = tmp))
}

# Repeatability  ----------------------------------------------------------
# Get the repeatability for all finches and traits 
for.rptl  = reap.finch(finch.recp,species = "fortis",trait.search = "Beak.Length")
for.rptd  = reap.finch(finch.recp,species = "fortis",trait.search = "Beak.Depth")
for.rptw  = reap.finch(finch.recp,species = "fortis",trait.search = "Beak.Width")
ful.rptl  = reap.finch(finch.recp,species = "fuliginosa",trait.search = "Beak.Length")
ful.rptd  = reap.finch(finch.recp,species = "fuliginosa",trait.search = "Beak.Depth")
ful.rptw  = reap.finch(finch.recp,species = "fuliginosa",trait.search = "Beak.Width")
mag.rptl  = reap.finch(finch.recp,species = "magnirostris",trait.search = "Beak.Length")
mag.rptd  = reap.finch(finch.recp,species = "magnirostris",trait.search = "Beak.Depth")
mag.rptw  = reap.finch(finch.recp,species = "magnirostris",trait.search = "Beak.Width")
sca.rptl  = reap.finch(finch.recp,species = "scandens",trait.search = "Beak.Length")
sca.rptd  = reap.finch(finch.recp,species = "scandens",trait.search = "Beak.Depth")
sca.rptw  = reap.finch(finch.recp,species = "scandens",trait.search = "Beak.Width")


# Mean repeatability ------------------------------------------------------
mean.rptr = mean(unlist(c(for.rptl$repeatability$R, for.rptd$repeatability$R, for.rptw$repeatability$R,
                          ful.rptl$repeatability$R, ful.rptd$repeatability$R, ful.rptw$repeatability$R,
                          mag.rptl$repeatability$R, mag.rptd$repeatability$R, mag.rptw$repeatability$R,
                          sca.rptl$repeatability$R, sca.rptd$repeatability$R, sca.rptw$repeatability$R)))

round(mean.rptr,2)

plot(for.rptl$repeatability, type = "permut", cex.main = 1)
plot(for.rptd$repeatability, type = "permut", cex.main = 1)
plot(for.rptw$repeatability, type = "permut", cex.main = 1)
