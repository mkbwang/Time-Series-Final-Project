library(pomp)
library(doParallel)
set.seed(2020,kind="L'Ecuyer")
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset='8'))
cl <- makeCluster(cores)
registerDoParallel(cl)

setwd(getwd())
# setwd("~/UM/Biostatistics/Courses/STATS531/Time-Series-Final-Project")

source("eagle_model_v2.R")

run_level <- 2
switch(run_level, {
  eagle_Np=100; eagle_Nmif=10; eagle_Neval=10;
  eagle_Nglobal=10; eagle_Nlocal=10
},{
  eagle_Np=30000; eagle_Nmif=150; eagle_Neval=10;
  eagle_Nglobal=50; eagle_Nlocal=10
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
eagle_box_p <- rbind(
  p0=c(-2,-0.1),
  p1=c(-2,-0.1),
  pmix=c(-2,-0.1)
)

eagle_box_s <- rbind(
  shape01=c(0.5, 10),
  scale01=c(0.01, 1),
  shape02=c(0.5, 10),
  scale02=c(0.01, 1),
  shape1=c(0.5, 10),
  scale1=c(0.01, 1)
)
# eagle_fixed_params <- c(shape1=8, scale1=0.002)

eagle_rw.sd <- 0.02; eagle_cooling.fraction.50 <- 0.6

# Iterative filtering
stew(file=sprintf("global_search-v2-%d.rda",run_level),{
  t_global <- system.time({
    mifs_global <- foreach(i=1:eagle_Nglobal,
                           .packages='pomp', .combine=c,
                           .inorder=FALSE,
                           .options.multicore=list(set.seed=TRUE),
                           .export=c("eagle_1", "eagle_box_p","eagle_box_s", "eagle_Np", "eagle_Nmif", "eagle_fixed_params",
                                     "eagle_cooling.fraction.50", "eagle_rw.sd")) %dopar%  {
                                       mif2(eagle_1,
                                            params=c(apply(eagle_box_p,1,function(x)exp(runif(1,x[1],x[2]))),
                                                     apply(eagle_box_s,1,function(x)runif(1,x[1],x[2]))),
                                            Np=eagle_Np,
                                            Nmif=eagle_Nmif,
                                            cooling.fraction.50=eagle_cooling.fraction.50,
                                            rw.sd=rw.sd(
                                              p0=eagle_rw.sd,
                                              p1=eagle_rw.sd,
                                              pmix=eagle_rw.sd,
                                              shape01=eagle_rw.sd,
                                              scale01=eagle_rw.sd,
                                              shape02=eagle_rw.sd,
                                              scale02=eagle_rw.sd,
                                              shape1=eagle_rw.sd,
                                              scale1=eagle_rw.sd)
                                       )
                                     }
  })
},seed=54321,kind="L'Ecuyer")


# log likelihood
stew(file=sprintf("lik_global_eval-v2-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:eagle_Nglobal,
                           .combine=rbind, .packages='pomp',
                           .inorder=FALSE,
                           .options.multicore=list(set.seed=TRUE),
                           .export=c("eagle_Neval", "eagle_1", "mifs_global", "eagle_Np")) %dopar% {
                             evals <- replicate(eagle_Neval,
                                                logLik(pfilter(eagle_1,
                                                               params=coef(mifs_global[[i]]),Np=eagle_Np)))
                             logmeanexp(evals, se=TRUE)
                           }
  })
},seed=98765,kind="L'Ecuyer")

