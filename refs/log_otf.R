rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)

# extract data from OTF log file
myfile = "/Users/xwang/GitHub/mydag/log_otf.txt"
nff.lines <- readLines(myfile)

cv.dedt <- nff.lines[grep("^T:\\s+(\\S+)\\s+Cv:\\s+(\\S+)\\s+dE\\/dT:\\s+(\\S+)", nff.lines)]
m2 <- t(sapply(strsplit(cv.dedt, "T:|Cv:|dE/dT:|Cv-"), function(this) this[c(2:4)]))
m2 <- apply(m2, 2, function(x) as.numeric(gsub("\\s", "", x)))
m2 <- data.frame(m2)

cost.stddev <- nff.lines[grep("^T:\\s+(\\S+)\\s+Cost:\\s+(\\S+)\\s+.*\\s+stddev:\\s+(\\S+)", nff.lines)]
m1 <- t(sapply(strsplit(cost.stddev, "T:|Cost:|stddev:|Ideal:"), function(this) this[c(2:4)]))
m1[, 2] <- gsub("\\+/-(\\d{+}.*)", "", m1[, 2])
m1 <- apply(m1, 2, function(x) as.numeric(gsub("\\s", "", x)))
colnames(m1) <- c("temps", "costs", "stddev")
m1 <- data.frame(m1)

cost.gelmc <- nff.lines[grep("^T:\\s+(\\S+)\\s+GelMC_cost:\\s+.*", nff.lines)]
m3 <- t(sapply(strsplit(cost.gelmc, "T:|R\\[1\\]:|R\\[1\\/2\\]"), function(this) this[c(2:3)]))
m3 <- apply(m3, 2, function(x) gsub("\\s", "", x))
m3 <- apply(m3, 2, function(x) as.numeric(gsub("GelMC_cost:", "", x)))
colnames(m3) <- c("temps", "GelMC")
m3 <- data.frame(m3)

plot(m1$costs[-1], xlab = "MCMC Progress", ylab = "Score")

plot(m2$CV[m2$temps < 100])
plot(m2$dEdT[m2$temps < 100])
plot(m3$GelMC[-1])
plot(m3$GelMC[m3$temps < 100])
plot(m3$temp[m3$temp < 4 & m3$temp > 2], m3$GelMC[m3$temp < 4 & m3$temp > 2])

mydt = m1[-1, ]
mydt$points = 1:nrow(mydt)

ggplot(mydt, aes(x=log10(temps), y=costs)) + geom_line() +
  geom_errorbar(aes(ymin=costs-stddev, ymax=costs+stddev), colour="black", width=.1) +
  scale_x_discrete(name = "Temperature", limits = c(0, 1, 2, 3, log10(5000)), labels = c("1", "10", "100", "1000", "5000")) +
  ylab("Cost")
