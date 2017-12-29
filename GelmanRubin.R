# 1. Make a OTF run with maxent on
# 2. Copy the OTF run in step 1 to a new directory. 
# cp -rf m1.maxent/ m1.maxent.rerun/

# 3. In new directory, run OTF from the checkpoint T = 1 file using the -R flag.
# use the -a flag to overwrite params from the original run
# Dump the ensemble every 20 npass to dump20.txt
# increase finalnpassfac to 100 from 50. (this means the total npass will be npass*100)
# with npass = 70, total walk number is 70 * 100, total sampe number saved is 70 * 100 / 20 = 350
# use 128 cores
# ./run_qsub_job -R ./checkpoint.txt-FINAL-T=1.0000000000000000 -a "-dump=20:dump20.txt -finalnpassfac=100 -maxent 0" 128

# 3. After this run, you will have a file called dump20.txt that has the entire ensemble every 20 npass 
# 4. Parse this file using a utility Karl wrote: /home/runge/lfurchtg/gelmc/dump2summary (perl)
# ./dump2summary 128 dump20.txt | tee summ20.txt

# format is:
# t: 239  net: 97 edges: 12670    score: 4644377.7877411

# 3a. re-format
# echo "t  net  edges  BIC" > summ20_tab.txt
# awk '{print $2 "  " $4 "  " $6 "  " $8}' summ20.txt >> summ20_tab.txt

# 4. Load the file into R and process locally

library(data.table)
library(ggplot2)
library(dplyr)

dt = fread("~/Github/mydag/summ20_tab.txt", data.table = T)
dt$net = as.factor(dt$net)

ggplot(dt, aes(x = t, y = BIC, color = net)) + geom_line(alpha=0.3) + theme(legend.position="none")
ggplot(dt, aes(x = t, y = edges, color = net)) + geom_line(alpha=0.3) + theme(legend.position="none")

# 6. Estimate linear drift of BIC and edge numbers
summary(lm(BIC~t, data = dt))$coefficients
summary(lm(edges~t, data = dt))$coefficients

# 7. Estimate the gelman-rubin criterion
gelmanrubin <- function(dt, X="BIC") { # X: "BIC" or "edges"
  m = length(unique(dt$net)) # number of networks
  n = length(unique(dt$t)) # number of samples per network 
  W = sapply(0:(m-1), function(x) var(dt$BIC[dt$net == x])) %>% mean # within network variance
  Bn = sapply(0:(m-1), function(x) mean(dt$BIC[dt$net == x])) %>% var # 
# W = sapply(0:(m-1), function(x) var(dt$edges[dt$net == x])) %>% mean # within network variance
# Bn = sapply(0:(m-1), function(x) mean(dt$edges[dt$net == x])) %>% var # 
  sig = (n-1)*W/n + Bn
  V = sig + Bn/m # between network variance
  return(V/W)
}

gelmanrubin(dt[t>50,], "BIC")
gelmanrubin(dt[t>50,], "edges")
