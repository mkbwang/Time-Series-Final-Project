library(pomp)
library(doParallel)
registerDoParallel()
library(doRNG)
registerDoRNG(3899882)

# setwd("~/UM/Biostatistics/Courses/STATS531/Time-Series-Final-Project")

source("eagle_model_v1.R")

run_level <- 2
switch(run_level, {
  eagle_Np=100; eagle_Nmif=10; bsflu_Neval=10;
  eagle_Nglobal=10; eagle_Nlocal=10
},{
  eagle_Np=20000; eagle_Nmif=100; eagle_Neval=10;
  eagle_Nglobal=10; eagle_Nlocal=10
},{
  eagle_Np=60000; eagle_Nmif=300; eagle_Neval=10;
  eagle_Nglobal=100; eagle_Nlocal=20}
)

# set starting parameters for local search

data = read.csv("https://raw.githubusercontent.com/skybullbobby/Time-Series-Final-Project/master/eagle_421.csv")

fit.gamma <- fitdist(data$msa, distr = "gamma", method = "mle")
startshape <- fit.gamma$estimate[[1]]
startscale <- 1/fit.gamma$estimate[[2]] 

eagle_startval <- c(p0=0.6, p1=0.4, shape0=startshape, scale0=startscale,
                    shape1=startshape, scale1=startscale)


eagle_box <- rbind(
  p0=c(-2,-0.1),
  p1=c(-2,-0.1),
  shape0=c(-1.5, -0.5),
  scale0=c(-1.5, -0.5),
  shape1=c(-1.5, -0.5),
  scale1=c(-1.5, -0.5)
)


eagle_rw.sd <- 0.02; eagle_cooling.fraction.50 <- 0.5

stew(file=sprintf("global_search-v1-%d.rda",run_level),{
  t_global <- system.time({
    mifs_global <- foreach(i=1:eagle_Nglobal,
                          .packages='pomp', .combine=c) %dopar%  {
                            mif2(eagle_0,
                                 params=apply(eagle_box,1,function(x)exp(runif(1,x[1],x[2]))),
                                 Np=eagle_Np,
                                 Nmif=eagle_Nmif,
                                 cooling.fraction.50=eagle_cooling.fraction.50,
                                 rw.sd=rw.sd(
                                   p0=eagle_rw.sd,
                                   p1=eagle_rw.sd,
                                   shape0=eagle_rw.sd,
                                   scale0=eagle_rw.sd,
                                   shape1=eagle_rw.sd,
                                   scale1=eagle_rw.sd)
                            )
                          }
  })
},seed=12345,kind="L'Ecuyer")

