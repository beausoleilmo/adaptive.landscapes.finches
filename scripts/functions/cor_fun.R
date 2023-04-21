# cor_fun is found here https://stackoverflow.com/questions/37889222/change-colors-in-ggpairs-now-that-params-is-deprecated
cor_fun <- function(data, mapping, method="pearson", adj.size = FALSE, ndp=2, sz=5, stars=FALSE, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  if (adj.size) {
    lb.size <- sz* abs(est)   
  } else {
    lb.size <- sz
    
  }
  
  
  if(stars){
    stars <- c("***", "**", "*", "")[findInterval(corr$p.value, c(0, 0.001, 0.01, 0.05, 1))]
    lbl <- paste0(round(est, ndp), stars)
  }else{
    lbl <- round(est, ndp)
  }
  
  ggplot(data=data, mapping=mapping) + 
    annotate("text", x=mean(x, na.rm=TRUE), y=mean(y, na.rm=TRUE), col = "black",label=lbl, size=lb.size,...)+
    theme(panel.grid = element_blank())
}
