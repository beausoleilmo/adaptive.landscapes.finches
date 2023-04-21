# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Created by Marc-Olivier Beausoleil
# McGill University 
# Created Tuesday, February 12, 2020 
# Why: Functions to compute various elements 
# Requires 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


# return sorted unique values ---------------------------------------------

id <- function(x) {
  unique(sort(as.vector(unlist(x))))
}

# standardize a vector ----------------------------------------------------
standardize <- function(x) {
  (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}

## quick scale function
sc.num <- function(x) {
  as.numeric(scale(x))
}

# Stats: mode -------------------------------------------------------------
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# SE ----------------------------------------------------------------------
stderr <- function(x) {
  sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
}

# expit and logit functions -----------------------------------------------
expit <- function(x) {
  1/(1+exp(-x))
}
logit <- function(x) {
  log(x/(1-x))
}

# invlogit and logit ------------------------------------------------------
# Function that will transforme the Y values from linear scale to [0,1]
invlogit <- function(x) {
  exp(x)/(exp(x) + 1)
}

logi <-     function(x) {
  exp(x)/(1+exp(x))
}


# load and return loaded object -------------------------------------------
load.local <- function(file) {
  v <- load(file)
  stopifnot(length(v) == 1)
  get(v)
}


# nice little pdf function ------------------------------------------------
pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

# Pairs plots -------------------------------------------------------------
# To create panels in pairs function  -------------------------------------
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gold", ...) 
}

# Pairs plot  -------------------------------------------------------------
autopairs<-function(x, ...) {
  pairs(x, 
        lower.panel = panel.smooth, 
        upper.panel = panel.cor, 
        diag.panel = panel.hist, ...)
}


# Correct capture history -------------------------------------------------
# Correct the capture history (APPARENT SURVIVAL). May take several seconds 
# Takes the first instance and the last and fills the 0s between first and last 
# (this "correct" the fact that some individuals might not be present but are indeed present)
# Otherwise, the models would think that the bird died, which is not true. 

known.state.cjs <- function(ch) {
  state <- ch
  for(i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- 1
  }
  state[state==0] <- 0
  return(state)
}

# Test data 
# test.ch = structure(list(y.1 = c(0, 0, 0, 0, 0, 0), y.2 = c(0, 0, 0, 0, 0, 1), y.3 = c(0, 0, 0, 0, 0, 0), 
#                y.4 = c(0, 0, 0, 0, 0, 1), y.5 = c(0, 0, 0, 0, 0, 0), y.6 = c(0, 0, 0, 0, 0, 1), 
#                y.7 = c(1, 1, 1, 0, 0, 1), y.8 = c(0, 0, 0, 1, 1, 0), y.9 = c(0, 0, 0, 0, 0, 0), 
#                y.10 = c(0, 0, 0, 0, 0, 0), y.11 = c(1, 1, 1, 0, 0, 0), y.12 = c(1, 1, 1, 0, 0, 0), 
#                y.13 = c(0, 0, 0, 1, 1, 0)), class = "data.frame", row.names = c(1L,2L,3L,4L,5L,6L))
# known.state.cjs(test.ch)
