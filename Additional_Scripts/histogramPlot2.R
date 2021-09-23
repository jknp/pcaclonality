histogramPlot2 <- function (ptLRvec, refLRvec) 
{
  par(lwd = 2)
  
  LNtab_s <- ptLRvec[which(ptLRvec[,10] < 0.05 & ptLRvec[,10] > 0.01),]
  LNtab_ss <- ptLRvec[which(ptLRvec[,10] < 0.01 & (ptLRvec[,10] > 0.001)),]
  LNtab_sss <- ptLRvec[which(ptLRvec[,10] < 0.001),]
  LNtab_ns <- ptLRvec[which(ptLRvec[,10] >= 0.05),]
  
  res <- ptLRvec[,4]
  ref <- refLRvec[,4]
  
  a <- hist(log(ref), xlim = c(min(c(log(ref), log(res)), na.rm = TRUE), max(c(log(ref), log(res)), na.rm = TRUE)), 
            breaks = 30, main = paste("Reference distribution of logLR (black), tested pairs (red *, green **, purple***)"), xlab = "")
  
  b <- hist(log(res), plot = FALSE)
  mult <- max(1, round(max(a$counts)/max(b$counts)))
  
  a <- hist(rep(log(LNtab_ns[,4]), mult), border = "blue", col = "blue", 
            lwd = 10, add = TRUE, density = 5, breaks = 30)
  
  a <- hist(rep(log(LNtab_s[,4]), mult), border = "red", col = "red", 
            lwd = 10, add = TRUE, density = 5, breaks = (hist(log(LNtab_s[,4]), plot = FALSE))$breaks)
  
  a <- hist(rep(log(LNtab_ss[,4]), mult), border = "green", col = "green", 
            lwd = 10, add = TRUE, density = 5, breaks = 5)
  
  a <- hist(rep(log(LNtab_sss[,4]), mult), border = "purple", col = "purple", 
            lwd = 10, add = TRUE, density = 5, breaks = 60)
  
  axis(side = 4, at = sort(unique(a$counts)), labels = as.character(sort(unique(a$counts))/mult), 
       adj = 0, col = "red", col.axis = "red")
}

ptLRvec <- LNsort40
#ptLRvec2 <- LNtab40ns[,4]
refLRvec <- ref_df40
