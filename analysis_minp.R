Sys.setenv(OMP_NUM_THREADS="1")

# Set number of cores here
no.cores = 30

# Working directory to the repository
setwd("Reproducible/")

source("functions.R")

setwd("data/")
dir = list.files(getwd())
dir

full.data = lapply(dir, read.csv)

setwd("..")


parallel.analysis = function(cores = 10, no.runs = 5, data = data, no.obs = 300, no.samples = 1, no.perm=10){
  
  library(rlist)
  library(ggplot2)
  library(mvtnorm)
  library(fastDummies)
  library(broom)
  library(dplyr)
  library(reshape2)
  library(viridis)
  library(doParallel)
  library(parallel)
  
  
  cl = makeCluster(cores)
  registerDoParallel(cl)
  
  
  x1 = foreach(j = 1:no.runs, .errorhandling="remove") %dopar%{
    
    source("functions.R")
    set.seed(j)
    sample_lst = analysis_effect(df = data, no.obs = no.obs, no.samples = no.samples, no.perm = no.perm) 

    return(sample_lst)
    
  }
  stopCluster(cl)
  return(x1)
}


set.seed(2)

sample50_lst = parallel.analysis(cores = no.cores, no.runs = 1000, data = full.data, no.obs = 50, no.samples = 1, no.perm=1000)
sample100_lst = parallel.analysis(cores = no.cores, no.runs = 1000, data = full.data, no.obs = 100, no.samples = 1, no.perm=1000)
sample150_lst = parallel.analysis(cores = no.cores, no.runs = 1000, data = full.data, no.obs = 150, no.samples = 1, no.perm=1000)
sample200_lst = parallel.analysis(cores = no.cores, no.runs = 1000, data = full.data, no.obs = 200, no.samples = 1, no.perm=1000)
sample250_lst = parallel.analysis(cores = no.cores, no.runs = 1000, data = full.data, no.obs = 250, no.samples = 1, no.perm=1000)
sample300_lst = parallel.analysis(cores = no.cores, no.runs = 1000, data = full.data, no.obs = 300, no.samples = 1, no.perm=1000)
sample500_lst = parallel.analysis(cores = no.cores, no.runs = 1000, data = full.data, no.obs = 500, no.samples = 1, no.perm=1000)



setwd("results/")

save(sample50_lst, file = "sample_50_effect.RData")
save(sample100_lst, file = "sample_100_effect.RData")
save(sample150_lst, file = "sample_150_effect.RData")
save(sample200_lst, file = "sample_200_effect.RData")
save(sample250_lst, file = "sample_250_effect.RData")
save(sample300_lst, file = "sample_300_effect.RData")
save(sample500_lst, file = "sample_500_effect.RData")

load("sample_50_effect.RData")
load("sample_100_effect.RData")
load("sample_150_effect.RData")
load("sample_200_effect.RData")
load("sample_250_effect.RData")
load("sample_300_effect.RData")
load("sample_500_effect.RData")


prop = function(data, sig = 0.05){
  unadj = mean(unlist((lapply(data, function(x) {  sum(x$unadjusted < sig)/nrow(x)}))))
  adj = mean(unlist((lapply(data, function(x) {  sum(x$adjusted < sig)/nrow(x)}))))
  bf = mean(unlist((lapply(data, function(x) {  sum(x$bonferroni < sig)/nrow(x)}))))
  res = c(unadj,adj,bf,sig)
  names(res) = c("unadjusted","adjusted","bonferroni","sig")
  return(res)
}

sig = as.data.frame(t(prop(sample50_lst, 0.05)))

sig = rbind(sig,
            prop(sample100_lst, 0.05),
            prop(sample150_lst, 0.05),
            prop(sample200_lst, 0.05),
            prop(sample250_lst, 0.05),
            prop(sample300_lst, 0.05),
            prop(sample500_lst, 0.05),
            prop(sample50_lst, 0.01),
            prop(sample100_lst, 0.01),
            prop(sample150_lst, 0.01),
            prop(sample200_lst, 0.01),
            prop(sample250_lst, 0.01),
            prop(sample300_lst, 0.01),
            prop(sample500_lst, 0.01),
            prop(sample50_lst, 0.1),
            prop(sample100_lst, 0.1),
            prop(sample150_lst, 0.1),
            prop(sample200_lst, 0.1),
            prop(sample250_lst, 0.1),
            prop(sample300_lst, 0.1),
            prop(sample500_lst, 0.1))

sig$size = rep(c(50,100,150,200,250,300,500), times = 3)

sig$size = as.factor(sig$size)
sig$sig = as.factor(sig$sig)

sig = sig %>%
  
  mutate(sig = recode(sig, "0.01" = "Sig. level: 1%", "0.05" = "Sig. level: 5%", "0.1" = "Sig. level: 10%"))



melted = melt(data          = sig,
              id.vars       = c("size","sig"),
              measure.vars  = c("unadjusted","adjusted","bonferroni"),
              variable.name = "method",
              value.name    = "value")
melted


################################################################################
### Create Figure 3
################################################################################


ggplot(data = melted, aes(x = size, y = value, color = method, group = method)) + 
      geom_point() +
      geom_line() +
      facet_wrap(~sig) +
      scale_color_discrete(name="Method",
                       labels=c("Unadjusted","minP-adjusted","Bonferroni-adjusted")) +
      theme_minimal(base_size = 30) +
      theme(axis.text.x = element_text(angle = 45)) +
      xlab("Sample size") +
      ylab("Proportion of sig. results") 


ggsave("figure_3.png", width = 20, height = 10, dpi = 300)

