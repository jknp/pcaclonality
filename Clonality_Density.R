#####Densityplot Clonality#####
library(ggplot2)
library(dplyr)

ref <- ref_df[,4]
res <- LNsort[,4]

refd <- data.frame(value = ref)
resd <- data.frame(res)

#Order LNsort for p value
res <- LNsort40
res <- res[order(LNsort40[,10]),]

LN40sig <- res[which(res[,10] <= 0.01),]
LN40ns <- res[which(res[,10] >= 0.01),] 

primd<- data.frame(value = LN40sig[,4])
smallp <- data.frame(value = LN40ns[,4])

#Random two groups
res <- LNsort40
res <- data.frame(res[,4])

primd <- res[sample(nrow(res),35),]
smallp <- res[sample(nrow(res),35),]

##Divide datasets
#multi dominant samples excl vec <- c(1,2,5,8,12,13,17,19,20,21,22,30,39,42,43,47,48,49,50,53,57,58,59,60,61,62,63,64,65,67)
#Exclude double prims
#[-c(23,27,35,40,41,54)]

vec <- c(1,2,5,8,12,13,17,19,20,21,22,30,39,42,43,47,48,49,50,53,57,58,59,60,61,62,63,64,65,67)
prim <- LNsort40[vec,4]
primd <- data.frame(value = prim)

#Full vector here
vec <- c(1,2,5,8,12,13,17,19,20,21,22,23,27,30,35,39,40,41,42,43,47,48,49,50,53,54,57,58,59,60,61,62,63,64,65,67)
small <- LNsort40[-vec,]
smallp <- data.frame(value = small[,4])

a <- ggplot() + geom_density(data = log(refd), aes(log(refd)), fill = "grey") + 
  geom_density(data = log(smallp), aes(log(smallp)), color = "navy", fill = "navy", alpha = 0.5) + 
  geom_density(data = log(primd), aes(log(primd)), color = "firebrick", fill = "firebrick", alpha = 0.5) +
  geom_vline(xintercept =  log(4.5e+01)) + labs(x = "log Likelihood Ratio", y = "Density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-10, 34)) + scale_y_continuous(expand= c(0,0), limits = c(0, 0.21)) +
  theme_bw() + theme(axis.text=element_text(size=rel(1.)), axis.title=element_text(size=rel(1.5)))
a

##### Boxplot of LR values #####
library(ggpubr)


refbox <- data.frame(group = "Reference", value = log(refd))
primbox <- data.frame(group = "P1-N", value = log(primd))
smallbox <- data.frame(group = "P2-n-N", value = log(smallp))

colnames(refbox) <- c("group", "value")
colnames(primbox) <- c("group", "value")
colnames(smallbox) <- c("group", "value")

Comp <- rbind(refbox, primbox, smallbox)

# comp_pairs <- list(c("Dominant", "nonDominant"), c("Reference", "Dominant"))
# stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Reference")

color <- c("grey", "firebrick", "navy")  

b <- ggplot(Comp, aes(x=group, y= value, fill = group)) + geom_boxplot() + 
  scale_fill_manual(values = color) +
  stat_compare_means(comparisons = list(c("P1-N", "P2-n-N"))) +
  theme_bw() 
b

c <- ggboxplot(Comp, x = "group", y= "value", fill = group)
c

compare_means(value ~ group,  data = Comp, ref.group = "Reference", method = "t.test")

#other pval intercepts (0.05 and 0.001 resp.) c(log(2.5e+00),log(1.8e+03)()

############################################################## V OLD V ######################################

##### Area under density #####
x0 <- log(4.5e+01)
bwidth <- 2.5

primd_d <- density(log(prim), bw = bwidth)
xx <- primd_d$x
dx <- xx[2] - xx[1]
yy <- primd_d$y
C <- sum(yy)*dx
p_unscaled <- sum(yy[xx >= x0]) * dx
scaled_p_prim <- p_unscaled / C

smallp_d <- density(log(small[,4]), bw = bwidth)
xx <- smallp_d$x
dx <- xx[2] - xx[1]
yy <- smallp_d$y
C <- sum(yy)*dx
p_unscaled <- sum(yy[xx >= x0]) * dx
scaled_p_small <- p_unscaled / C

# Result tally
Ratio <- scaled_p_prim / scaled_p_small
primd_count <- count(LNsort40[which(LNsort40[vec,10] < 0.01),])
smallp_count <- count(LNsort40[which(LNsort40[-vec,10] < 0.01),])

#bwidth originally is 3.178 for complete resd
density(log(LNsort40[,4]))


##### Old code #####
## LRpvalue histogram
# library(grid)
# resp <- data.frame(LNsort40[,10])
# 
# breaks <- hist(log(ref), plot = F)$breaks
# 
# ggplot() + geom_histogram(data = (resp), aes(resp), binwidth = 0.001, na.rm = TRUE) + xlim (0, 0.1) + ylim (0,6) + theme_minimal()
# 
# 
# b <- ggplot() + geom_histogram(data = log(resd), aes(log(resd), fill = "red"), breaks = seq(-10,30))
# b <- b + scale_y_continuous(sec.axis = sec_axis(~.*0.005))
# 
# grid.newpage()
# grid.draw(rbind(ggplotGrob(a), ggplotGrob(b), size = "last"))


#a <- hist(log(ref), xlim = c(min(c(log(ref), log(res)), na.rm = TRUE), max(c(log(ref), log(res)), na.rm = TRUE)), 
          #breaks = 30, main = paste("Reference distribution of logLR (black), p<0.01 tested pairs (red)"), xlab = "")

#How clonality calculates LR2pvalue:
#LR2pvalue <- c(LR2pvalue, mean(ptLR[i,4] <= refLR[, 4], na.rm = TRUE))