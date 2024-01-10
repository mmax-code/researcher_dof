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
drop = lapply(full.data, function(x) `[<-`(x, 'complication', value = NULL))

setwd("..")

parallel.fwer = function(cores = cores, no.runs = 5, drop = drop, no.obs = 300, no.samples = 5, no.perm=50){
  
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
  
  
  x1 = foreach(j = 1:no.runs) %dopar%{
    
    source("functions.R")
    
    print(j)
    set.seed(j)
    
    df = lapply(drop, function(x) cbind(x, complication = rbinom(n=nrow(x),size=1,prob=0.5)))
    
    res = sample.analysis(data = df,size = no.obs,times = no.samples)
    
    data = res[[2]]
    
    p.unadj = do.call(c, unlist(res[[1]], recursive=FALSE))
    sum(p.unadj < 0.05)/length(p.unadj)
    result.unadj = c(any(p.unadj < 0.05))
    
    test = permutation(data = data, iter = no.perm)
    p.permuted = as.data.frame(test)
    min.permuted = apply(p.permuted,2,min)
    
    p.adjusted = c()
    for (i in 1:length(p.unadj)){
      p.adjusted = c(p.adjusted,sum(min.permuted < p.unadj[i])/length(min.permuted))
    }
    
    sum(p.unadj < 0.05) / length(p.unadj)
    sum(p.adjusted < 0.05) / length(p.adjusted)
    
    result.adj = c(any(p.adjusted < 0.05))
    
    
    bonferroni = p.adjust(p.unadj, method = "bonferroni")
    sum(bonferroni < 0.05) / length(bonferroni)
    
    result.bon = c(any(bonferroni < 0.05))
    
    return(list(unadj = result.unadj, adj = result.adj, bf = result.bon))
    
  }
  stopCluster(cl)
  return(x1)
}


set.seed(1)
test_100_1 = parallel.fwer(cores = no.cores, no.runs = 1000, drop = drop, no.obs = 100, no.samples = 1, no.perm = 1000)
test_200_1 = parallel.fwer(cores = no.cores, no.runs = 1000, drop = drop, no.obs = 200, no.samples = 1, no.perm = 1000)
test_300_1 = parallel.fwer(cores = no.cores, no.runs = 1000, drop = drop, no.obs = 300, no.samples = 1, no.perm = 1000)
test_500_1 = parallel.fwer(cores = no.cores, no.runs = 1000, drop = drop, no.obs = 500, no.samples = 1, no.perm = 1000)
test_2000_1 = parallel.fwer(cores = no.cores, no.runs = 1000, drop = drop, no.obs = 2000, no.samples = 1, no.perm = 1000)
test_3000_1 = parallel.fwer(cores = no.cores, no.runs = 1000, drop = drop, no.obs = 3000, no.samples = 1, no.perm = 1000)



setwd("results/")
save(test_100_1, file = "test_100_1.RData")
save(test_200_1, file = "test_200_1.RData")
save(test_300_1, file = "test_300_1.RData")
save(test_500_1, file = "test_500_1.RData")
save(test_2000_1, file = "test_2000_1.RData")
save(test_3000_1, file = "test_3000_1.RData")


res.df = rbind(as.data.frame(final.res(test_100_1)),
               as.data.frame(final.res(test_200_1)),
               as.data.frame(final.res(test_300_1)),
               as.data.frame(final.res(test_500_1)),
               as.data.frame(final.res(test_2000_1)),
               as.data.frame(final.res(test_3000_1)))

res.df$sample = c(100,200,300,500,2000,3000)
res.df



melted = melt(data          = res.df,
              id.vars       = c("sample"),
              measure.vars  = c("unadjusted","adjusted","bf"),
              variable.name = "method",
              value.name    = "value")
melted
melted$ci.lower = c(res.df$unadjusted.ci.lower,res.df$adjusted.ci.lower,res.df$bf.ci.lower)
melted$ci.upper = c(res.df$unadjusted.ci.upper,res.df$adjusted.ci.upper,res.df$bf.ci.upper)
melted$sample = as.factor(melted$sample)



################################################################################
### Create Figure 2
################################################################################

p = ggplot(melted,aes(x=sample,y=value))
final_plot = p + geom_point(aes(group=method,color=method)) +
  geom_errorbar(aes(ymin=ci.lower,ymax=ci.upper,color=method,width=0.2)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
  xlab("Sample size") +
  ylab("FWER") +
  scale_color_discrete(name="Method",
                       labels=c("Unadjusted","minP-adjusted","Bonferroni-adjusted")) +
  theme_minimal(base_size = 30)

ggsave("figure_2.png", width = 20, height = 10, dpi = 300)

