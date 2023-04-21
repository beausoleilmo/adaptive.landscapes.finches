# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created 12 March 2020
# Why: 
  # Functions to plot the gams (modified from vis.gam {mgcv})
# Requires 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

plot.gam.cust <- function(mod,
                          data = mydata, survival = "y", 
                          n.grid = 50,sp = "sp",nCol=100,
                          type = "link", col.sch = "heat",
                          too.far=0, bss=NULL,
                          se = NULL,
                          view = c("X1","X2"), kk,
                          title ="", bdr = TRUE,
                          theta=40,phi=40, zlim = NULL, 
                          smooth.par = -4,
                          plots = c(TRUE,TRUE), # add 2d and 3D plot by default 
                          my.col = c('black',"pink",
                                     "darkgreen","blue"),
                          cex.surv = .7) {
  kk.t=ifelse(kk==-1,"default",kk)
  if (plots[1]) {
    cat("Plotting in 2D", fill = TRUE)
    vis.gam(mod,view=view, main ="", too.far = too.far, se = -1, 
            contour.col = "red", zlim = zlim,
            color=col.sch, n.grid=n.grid, type=type, plot.type="contour", nCol=50)  
    
    if(length(kk.t)>1){kk.t = round(kk.t[1],2)}
    title(main = paste(title,"Type:",bss,"in",kk.t,
                       "dim.,","sp=",
                       smooth.par,sep = " "),
          cex.main = .8, col.main= "black")
    points(x = data[,view[1]], y = data[,view[2]], pch = 21, 
           bg = scales::alpha(my.col[as.factor(data[,sp])],alpha = .4), 
           col = scales::alpha(my.col[as.factor(data[,sp])],alpha = .5))
    get.surv.max = max(data[,survival], na.rm = TRUE)
    if (get.surv.max==1) {
      points(x = data[data[,survival] %in% 1 ,view[1]], 
             y = data[data[,survival] %in% 1 ,view[2]], 
             pch = 21, bg = scales::alpha("yellow",.8), col = scales::alpha("yellow",.1), 
             cex = cex.surv) # plot only the one that survived
    }
    abline(h = 0, v = 0,  lty = 3)
  }
  
  cat("calculating fit...", fill = TRUE)
  v.names <- names(mod$var.summary)
  r1 <- range(mod$var.summary[[view[1]]])
  m1 <- seq(r1[1], r1[2], length = n.grid)
  
  r2 <- range(mod$var.summary[[view[2]]])
  m2 <- seq(r2[1], r2[2], length = n.grid)
  
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  cond = list()
  for (i in 1:length(mod$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- mod$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(mod$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(mod$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  zlab <- type
  fv <- predict.gam(mod, newdata = newd, se.fit = TRUE, type = type)
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, mod$model[, view[1]], 
                             mod$model[, view[2]], dist = too.far)
    fv$fit[ex.tf] <- NA
    # fv$se.fit[ex.tf] <- NA
  }
  z <- fv$fit
  
  max.z <- max(z, na.rm = TRUE)
  # z[is.na(z)] <- max.z * 10000
  z <- matrix(z, n.grid, n.grid)
  
  cat("calculating SE...", fill = TRUE)
  z.se.l <- fv$fit - fv$se.fit * se
  z.se.u <- fv$fit + se * fv$se.fit
  z.se.l <- matrix(z.se.l, n.grid, n.grid)
  z.se.u <- matrix(z.se.u, n.grid, n.grid)
  if (!is.null(zlim)) {
    if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
      stop("Something wrong with zlim")
    min.z <- zlim[1]
    max.z <- zlim[2]
  }
  
  if (plots[2]) {
    cat("Plotting in 3D", fill = TRUE)
    
    if (!(is.null(se))) {
      persp(x = m1,
            y = m2,
            z = z.se.l,
            col = NA,
            box = FALSE, axes = FALSE,
            border = scales::alpha("red",.1),
            zlim = zlim,xlab = "",ylab = "", zlab = "",
            theta=theta, phi =phi)
      par(new = TRUE)
      vis.gam(mod, view = view, main = "3D view", too.far = too.far, se = -1, 
              contour.col = "yellow",
              theta=theta, phi = phi, color=col.sch, n.grid = n.grid, zlim = zlim,
              ticktype="detailed",type=type, plot.type="persp", border = bdr)
      
      par(new = TRUE)
      persp(x = m1,
            y = m2,
            z = z.se.u,
            col = NA,
            box = FALSE, axes = FALSE,
            border = scales::alpha("blue",.1),
            zlim = zlim,xlab = "",ylab = "", zlab = "",
            theta=theta,phi =phi,
      )
    } else {
      se = -1
      vis.gam(mod, view = view, main = "3D view", too.far = too.far, 
              se = se, 
              contour.col = "yellow",
              theta=theta, phi = phi, color=col.sch, n.grid = n.grid, zlim = zlim,
              ticktype="detailed",type=type, plot.type="persp", border = bdr)
    }
  }
  # z <- matrix(fv$fit, n.grid, n.grid)
  # # z <- fv$fit - fv$se.fit * se
  # # z <- matrix(z, n.grid, n.grid)
  # z <- fv$fit
  # z <- matrix(z, n.grid, n.grid)
  # z <- fv$fit + se * fv$se.fit
  # z <- matrix(z, n.grid, n.grid)
  return(list(m1 = m1,m2 =m2,z=z, z.se.l= z.se.l,z.se.u = z.se.u))
}

# plot.gam.cust1d <- function(mod,bss=NULL,kk,title ="", bdr = TRUE) {
#   kk.t=ifelse(kk==-1,"default",kk)
#   vis.gam(mod,view=c("X1"), main ="",
#           color="heat",n.grid=50, type="link", plot.type="contour", nCol=50)
#   title(main = paste(title,"Type:",bss,"in",kk.t,"dim.,","sp=",smooth.par,sep = " "),
#         cex.main = .8, col.main= "black")
#   points(x = mydata$X1, y = mydata$y, pch = 21, bg = mydata$sp, col = mydata$sp)
#   points(x = mydata[mydata$y %in% 1 ,"X1"], y = mydata[mydata$y %in% 1 ,"X2"], pch = 21, bg = "yellow", col = "yellow", cex = .7) # plot only the one that survived
#   vis.gam(mod,view=c("X1"), main = "3D view",
#           theta=40,phi=40,color="heat",n.grid=50, ticktype="detailed",type="link", plot.type="persp", border = bdr)
# }

plot.gcv.score1 <- function(mod,data,bss,kk) {
  lambda <- exp( seq(-20,10, by=.8))        # fit a range of lambdas >0
  gcvscore <- sapply(lambda, function(lambda, data){
    gam(y~s(X1, bs = bss, k = kk)+s(X2, bs = bss, k = kk), 
        family = binomial(link = "logit"),
        sp = lambda,
        data=data, method="GCV.Cp")$gcv.ubre},
    mydata)
  
  plot(log(lambda), gcvscore, type = "l", 
       main = paste(yr.list[1], 
                    yr.list[2],sep = "-"),
       ylab ="GCV score",
       xlab = "ln(lambda)")
  abline(h = mod$gcv.ubre, 
         v = log(mod$full.sp[1]),
         lty = 3)
}

plot.gcv.score2 <- function(mod,data, bss,kk) {
  lambda <- exp( seq(-20,10, by=.8))        # fit a range of lambdas >0
  
  gcvscore <- sapply(lambda, function(lambda, data){
    gam(y~s(X1, bs = bss, k = kk) + s(X2, bs = bss, k = kk) + s(X1,X2, bs = bss, k = kk.i), 
        family = binomial(link = "logit"),
        sp = lambda,
        data=mydata, method="GCV.Cp")$gcv.ubre},
    mydata)
  
  plot(log(lambda), gcvscore, type = "l", 
       main = paste(yr.list[1], 
                    yr.list[2],sep = "-"),
       ylab ="GCV score",
       xlab = "ln(lambda)")
  abline(h = mod$gcv.ubre, 
         v = log(mod$full.sp[1]),
         lty = 3)
}



fac.seq <- function(fac, n.grid) {
  fn <- length(levels(fac))
  gn <- n.grid
  if (fn > gn) 
    mf <- factor(levels(fac))[1:gn]
  else {
    ln <- floor(gn/fn)
    mf <- rep(levels(fac)[fn], gn)
    mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
    mf <- factor(mf, levels = levels(fac))
  }
  mf
}


vis.gam.mb = function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
                       col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                       plot.type = "persp", zlim = NULL, nCol = 50, ...) 
{
  z.se.l = NULL
  z.se.u = NULL
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      }
      else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) 
        break
    }
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
  }
  else {
    if (sum(view %in% v.names) != 2) 
      stop(gettextf("view variables must be one of %s", 
                    paste(v.names, collapse = ", ")))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
                                 c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(gettextf("View variables must contain more than one value. view = c(%s,%s).", 
                  view[1], view[2]))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  if (type == "link") 
    zlab <- paste("linear predictor")
  else if (type == "response") 
    zlab <- type
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  z <- fv$fit
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  if (se <= 0) {
    old.warn <- options(warn = -1)
    av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
                 n.grid - 1)
    options(old.warn)
    max.z <- max(z, na.rm = TRUE)
    z[is.na(z)] <- max.z * 10000
    z <- matrix(z, n.grid, n.grid)
    surf.col <- t(av) %*% z %*% av
    surf.col[surf.col > max.z * 2] <- NA
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      min.z <- min(fv$fit, na.rm = TRUE)
      max.z <- max(fv$fit, na.rm = TRUE)
    }
    surf.col <- surf.col - min.z
    surf.col <- surf.col/(max.z - min.z)
    surf.col <- round(surf.col * nCol)
    con.col <- 1
    if (color == "heat") {
      pal <- heat.colors(nCol)
      con.col <- 4
    }
    else if (color == "topo") {
      pal <- topo.colors(nCol)
      con.col <- 2
    }
    else if (color == "cm") {
      pal <- cm.colors(nCol)
      con.col <- 1
    }
    else if (color == "terrain") {
      pal <- terrain.colors(nCol)
      con.col <- 2
    }
    else if (color == "gray" || color == "bw") {
      pal <- gray(seq(0.1, 0.9, length = nCol))
      con.col <- 1
    }
    else stop("color scheme not recognised")
    if (is.null(contour.col)) 
      contour.col <- con.col
    surf.col[surf.col < 1] <- 1
    surf.col[surf.col > nCol] <- nCol
    if (is.na(col)) 
      col <- pal[as.array(surf.col)]
    z <- matrix(fv$fit, n.grid, n.grid)
    if (plot.type == "contour") {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                    sep = "")
      if (color != "bw") {
        txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", 
                     ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                     sep = "")
        eval(parse(text = txt))
      }
      else {
        txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
    else {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("zlab" %in% dnm, "", ",zlab=zlab"), ",...)", 
                    sep = "")
      if (color == "bw") {
        op <- par(bg = "white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                     stub, sep = "")
        eval(parse(text = txt))
        par(op)
      }
      else {
        txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
  }
  else {
    if (color == "bw" || color == "gray") {
      subs <- paste("grey are +/-", se, "s.e.")
      lo.col <- "gray"
      hi.col <- "gray"
    }
    else {
      subs <- paste("red/green are +/-", se, "s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      max.z <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
      min.z <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
      zlim <- c(min.z, max.z)
    }
    # z <- fv$fit - fv$se.fit * se
    # z <- matrix(z, n.grid, n.grid)
    z.se.l <- fv$fit - fv$se.fit * se
    z.se.l <- matrix(z.se.l, n.grid, n.grid)
    if (plot.type == "contour") 
      warning("sorry no option for contouring with errors: try plot.gam")
    stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                  ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                                                                         dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                                                                                                        "", ",sub=subs"), ",...)", sep = "")
    txt <- paste("persp(m1,m2,z.se.l,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=lo.col"), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=\"black\""), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z.se.u <- fv$fit + se * fv$se.fit
    z.se.u <- matrix(z.se.u, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z.se.u,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=hi.col"), stub, sep = "")
    eval(parse(text = txt))
  }
  return(list(m1 = m1,m2 =m2,z=z, z.se.l= z.se.l,z.se.u = z.se.u, col = col))
  
}




addlables.persp = name <- function(x, y, z, 
                                   pmat, 
                                   tick.steps=c(2,5,0.5), # X, Y, Z
                                   pos.text = c(-0.25,0.25,- 0.5),# X, Y, Z
                                   pos.lab = c(-0.25,0.25,- 0.5),# X, Y, Z
                                   col = "black", lwd = 2, tex.cex = 1,
                                   xlab = "",
                                   ylab = "",
                                   zlab = "") {
  min.x  <- min(x, na.rm = TRUE)
  max.x  <- max(x, na.rm = TRUE)
  x.axis <- seq(floor(min.x),ceiling(max.x),by=tick.steps[1]) 
  min.y  <- min(y, na.rm = TRUE)
  max.y  <- max(y, na.rm = TRUE)
  y.axis <- seq(floor(min.y)+1,ceiling(max.y), by = tick.steps[2]) 
  min.z  <- round(min(z, na.rm = TRUE))
  max.z  <- max(z, na.rm = TRUE)
  z.axis <- seq(min.z, max.z, by=tick.steps[3]) 
  
  # X-axis
  tick.start <- trans3d(x.axis, max.y, min.z, pmat)
  tick.end   <- trans3d(x.axis, (max.y + 0.50), min.z, pmat = pmat)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, col = col, lwd = lwd)
  #Note the (min.y - 0.20) in the calculation of tick.end. 
  # This places the second line, parallel to the X axis, at the position -0.20 on the Y axis (i.e., into negative/unplotted space).
  # Y-axis
  tick.start <- trans3d(min.x, y.axis, min.z, pmat)
  tick.end   <- trans3d(min.x - 0.50, y.axis, min.z, pmat)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, col = col, lwd = lwd)
  # Z-axis
  tick.start <- trans3d(min.x, min.y, z.axis, pmat)
  tick.end <- trans3d(min.x, (min.y - 0.50), z.axis, pmat)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, col = col, lwd = lwd)
  
  # Adding the text 
  # The adj=c(0, NA) expression is used to left-justify the labels, 
  # the srt=270 expression is used to rotate the labels 270°, and 
  # the cex=0.5 expression is used to scale the label text to 75% of its original size.
  # X-axis
   labels <- as.character(x.axis)
   label.pos <- trans3d(x.axis, (max.y + pos.text[1]), min.z, pmat)
   text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), srt=0, cex=tex.cex)
   
   labels <- as.character(xlab)
   label.pos <- trans3d(mean(x.axis), (max.y + pos.lab[1]), min.z, pmat)
   text(label.pos$x, label.pos$y, labels=labels, adj=c(0.5, NA), srt=295, cex=tex.cex)
  
  # Y-axis
   labels <- as.character(y.axis)
   label.pos <- trans3d((min.x + pos.text[2]), y.axis, min.z, pmat)
   text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=0, cex=tex.cex)
 
   labels <- as.character(ylab)
   label.pos <- trans3d((min.x + pos.lab[2]), mean(y.axis), min.z, pmat)
   text(label.pos$x, label.pos$y, labels=labels, adj=c(0.5, NA), srt=30, cex=tex.cex)
  
  # Z-axis
   labels <- as.character(z.axis)
   label.pos <- trans3d(min.x, (min.y + pos.text[3]), z.axis, pmat)
   text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), srt=0, cex=tex.cex)
   
   labels <- as.character(zlab)
   label.pos <- trans3d(min.x, (min.y + pos.lab[3]), mean(z.axis), pmat)
   text(label.pos$x, label.pos$y, labels=labels, adj=c(0.5, NA), srt=80, cex=tex.cex)
  
  ######################1
  return(list(x.axis = x.axis, y.axis =y.axis, z.axis = z.axis))
}
