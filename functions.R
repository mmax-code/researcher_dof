# options(warn=-1)

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



analysis_p = function(data = data){
  data$cat_median = as.factor(data$cat_median)
  data$cat_mean = as.factor(data$cat_mean)
  # log_reg median
  model = glm(complication ~ median,family=binomial(link='logit'),data=data)
  sum = summary(model)
  effect_median = sum$coefficients[2,1]
  p_median = sum$coefficients[2,4]
  
  # log_reg mean
  model = glm(complication ~ mean,family=binomial(link='logit'),data=data)
  sum = summary(model)
  effect_mean = sum$coefficients[2,1]
  p_mean = sum$coefficients[2,4]
  
  # odds median
  odds = table(data$pao2_dummy_median, data$complication)
  
  odds_fisher = fisher.test(odds)
  odds_median =  as.numeric(odds_fisher$estimate)
  odds_median_p = odds_fisher$p.value
  
  # odds_mean
  odds = table(data$pao2_dummy_mean, data$complication)
  
  odds_fisher = fisher.test(odds)
  odds_fisher
  odds_mean = as.numeric(odds_fisher$estimate)
  odds_mean_p = odds_fisher$p.value
  
  
  odds = table(data$cat_median, data$complication)
  odds_fisher = fisher.test(odds)
  odds_fisher
  odds_cat_median = as.numeric(odds_fisher$estimate)
  odds_cat_median_p = odds_fisher$p.value
  
  odds = table(data$cat_mean, data$complication)
  odds_fisher = fisher.test(odds)
  odds_fisher
  odds_cat_mean = as.numeric(odds_fisher$estimate)
  odds_cat_mean_p = odds_fisher$p.value
  
  lst = list(p_median,p_mean,odds_median_p,
             odds_mean_p, odds_cat_median_p, odds_cat_mean_p)
  names(lst) = c("log_p_median","log_p_mean","odds_median_p",
                 "odds_mean_p","odds_cat_median_p","odds_cat_mean_p")
  return(lst)
}



sample.analysis = function(data = data, size = 300, times = 5){
  res.lst = list()
  data.lst = list()
  for (i in 1:times){
    tmp = lapply(data, function(x) x[sample(nrow(x),size), ] ) 
    res = lapply(tmp, function(x) analysis_p(x))
    res.lst[[i]] = res
    data.lst[[i]] = tmp
  }
  return(list(res.lst, data.lst))
}


permutation = function(data = data, iter = 500){
  p = c()
  for (i in 1:iter){
    set.seed(i)
    test = unlist(data,recursive=F)  
    tmp = lapply(test, function(x) {x[,"complication"] = sample(x[,"complication"]); x})
    tmp_p = lapply(tmp, function(x) analysis_p(x))
    res = as.vector(unlist(tmp_p))
    p = c(p,res)
  }
  matrix = matrix(p,nrow=length(p)/iter)
  return(matrix)
}



final.res = function(x1 = object){
  unadj = lapply(x1, function(l) l[[1]])
  adj = lapply(x1, function(l) l[[2]])
  bonferroni = lapply(x1, function(l) l[[3]])
  
  fwer.unadj = sum(unadj== TRUE)/length(x1)
  fwer.unadj.ci = prop.test(sum(unadj== TRUE),length(x1))$conf.int
  
  fwer.adj = sum(adj == TRUE)/length(x1)
  fwer.adj.ci = prop.test(sum(adj== TRUE),length(x1))$conf.int
  
  fwer.bf = sum(bonferroni == TRUE)/length(x1)
  fwer.bf.ci = prop.test(sum(bonferroni == TRUE),length(x1))$conf.int
  
  result = list(fwer.unadj,fwer.unadj.ci[1],fwer.unadj.ci[2], 
                fwer.adj,fwer.adj.ci[1], fwer.adj.ci[2], 
                fwer.bf, fwer.bf.ci[1],fwer.bf.ci[2])
  names(result) = c("unadjusted", "unadjusted.ci.lower","unadjusted.ci.upper",
                    "adjusted", "adjusted.ci.lower", "adjusted.ci.upper", 
                    "bf", "bf.ci.lower","bf.ci.upper")
  return(result)
}


analysis_effect = function(df, no.obs = 300, no.samples = 1, no.perm = 1000){
  res = sample.analysis(data = df,size = no.obs,times = no.samples)
  data = res[[2]]
  p.unadj = do.call(c, unlist(res[[1]], recursive=FALSE))
  sum(p.unadj < 0.05)/length(p.unadj)
  # result.unadj = c(any(p.unadj < 0.05))
  result.unadj = p.unadj 
  test = permutation(data = data, iter = no.perm)
  p.permuted = as.data.frame(test)
  min.permuted = apply(p.permuted,2,min)
  
  p.adjusted = c()
  for (i in 1:length(p.unadj)){
    p.adjusted = c(p.adjusted,sum(min.permuted < p.unadj[i])/(length(min.permuted)))
  }
  
  sum(p.unadj < 0.05) / length(p.unadj)
  sum(p.adjusted < 0.05) / length(p.adjusted)
  # result.adj = c(any(p.adjusted < 0.05))
  result.adj = p.adjusted
  bonferroni = p.adjust(p.unadj, method = "bonferroni")
  sum(bonferroni < 0.05) / length(bonferroni)
  
  # result.bon = c(any(bonferroni < 0.05))
  result.bon = bonferroni 
  final = list(unadj = result.unadj, adj = result.adj, bf = result.bon)
  unadj = as.vector(unlist(final[[1]]))
  adj = as.vector(unlist(final[[2]]))
  bonferroni = as.vector(unlist(final[[3]]))
  result = list(unadj,adj,bonferroni)
  names(result) = c("unadjusted", "adjusted", "bonferroni")
  df.final = as.data.frame(result)
  df.final$size = rep(no.obs, nrow(df.final))
  return(df.final)
}


proportions = function(df){
  unadj = sum(df$unadjusted < 0.05)/nrow(df)
  adj = sum(df$adjusted < 0.05)/nrow(df)
  bonferroni = sum(df$bonferroni < 0.05)/nrow(df)
  return(c(unadj,adj,bonferroni))
}