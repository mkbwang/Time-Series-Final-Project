library(pomp)
library(doParallel)
set.seed(2020,kind="L'Ecuyer")
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset='8'))
cl <- makeCluster(cores)
registerDoParallel(cl)

setwd(getwd())
# setwd("~/UM/Biostatistics/Courses/STATS531/Time-Series-Final-Project")

source("eagle_model_v1.R")

run_level <- 1
switch(run_level, {
  eagle_Np=100; eagle_Nmif=10; eagle_Neval=10;
  eagle_Nglobal=10; eagle_Nlocal=10
},{
  eagle_Np=20000; eagle_Nmif=100; eagle_Neval=10;
  eagle_Nglobal=10; eagle_Nlocal=10
},{
  eagle_Np=60000; eagle_Nmif=300; eagle_Neval=10;
  eagle_Nglobal=100; eagle_Nlocal=20}
)

# set starting parameters for local search

# fit.gamma <- fitdist(eagle_data$msa, distr = "gamma", method = "mle")
# startshape <- fit.gamma$estimate[[1]]
# startscale <- 1/fit.gamma$estimate[[2]] 
# 
# eagle_startval <- c(p0=0.6, p1=0.4, shape0=startshape, scale0=startscale,
#                     shape1=startshape, scale1=startscale)

# box for global search
eagle_box <- rbind(
  p0=c(-2,-0.1),
  p1=c(-2,-0.1),
  shape0=c(-1.5, -0.1),
  scale0=c(-1.5, -0.1),
  shape1=c(-1.5, -0.1),
  scale1=c(-1.5, -0.1)
)

eagle_rw.sd <- 0.02; eagle_cooling.fraction.50 <- 0.5

# Iterative filtering
stew(file=sprintf("global_search-v1-%d.rda",run_level),{
  t_global <- system.time({
    mifs_global <- foreach(i=1:eagle_Nglobal,
                          .packages='pomp', .combine=c,
                          .inorder=FALSE,
                          .options.multicore=list(set.seed=TRUE),
                          .export=c("eagle_0", "eagle_box", "eagle_Np", "eagle_Nmif",
                                    "eagle_cooling.fraction.50", "eagle_rw.sd")) %dopar%  {
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


# log likelihood
stew(file=sprintf("lik_global_eval-v1-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:eagle_Nglobal,
                           .combine=rbind, .packages='pomp',
                           .inorder=FALSE,
                           .options.multicore=list(set.seed=TRUE),
                           .export=c("eagle_Neval", "eagle_0", "mifs_global", "eagle_Np")) %dopar% {
                             evals <- replicate(eagle_Neval,
                                                logLik(pfilter(eagle_0,
                                                               params=coef(mifs_global[[i]]),Np=eagle_Np)))
                             logmeanexp(evals, se=TRUE)
                           }
  })
},seed=56789,kind="L'Ecuyer")

