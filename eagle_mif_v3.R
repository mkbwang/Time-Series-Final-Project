library(pomp)
library(doParallel)
set.seed(2020,kind="L'Ecuyer")
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset='8'))
cl <- makeCluster(cores)
registerDoParallel(cl)

setwd(getwd())
# setwd("~/UM/Biostatistics/Courses/STATS531/Time-Series-Final-Project")

source("eagle_model_v3.R")

run_level <- 2
switch(run_level, {
  eagle_Np=100; eagle_Nmif=10; eagle_Neval=10;
  eagle_Nglobal=10; eagle_Nlocal=10
},{
  eagle_Np=20000; eagle_Nmif=100; eagle_Neval=10;
  eagle_Nglobal=30; eagle_Nlocal=10
},{
  eagle_Np=60000; eagle_Nmif=300; eagle_Neval=10;
  eagle_Nglobal=100; eagle_Nlocal=20}
)


# fixed parameters
eagle_fixed_params <- read.csv("fixed_params.csv")
fixed_params <- eagle_fixed_params$x
names(fixed_params) <- eagle_fixed_params$X

# box for global search
eagle_box_beta <- rbind(
  beta00 = c(-4,-1),
  beta01 = c(0,2),
  beta10 = c(0,3),
  beta11 = c(0,2)
)

# eagle_fixed_params <- c(shape1=8, scale1=0.002)

eagle_rw.sd <- 0.02; eagle_cooling.fraction.50 <- 0.6

# Iterative filtering
stew(file=sprintf("global_search_covar1_mix-%d.rda",run_level),{
  t_global <- system.time({
    mifs_global_covar1_mix <- foreach(i=1:eagle_Nglobal,
                                       .packages='pomp', .combine=c,
                                       .inorder=FALSE,
                                       .options.multicore=list(set.seed=TRUE),
                                       .export=c("eagle_2", "eagle_box_beta", "eagle_Np", "eagle_Nmif", "fixed_params",
                                                 "eagle_cooling.fraction.50", "eagle_rw.sd")) %dopar%  {
                                                   mif2(eagle_2,
                                                        params=c(apply(eagle_box_beta,1,function(x)runif(1,x[1],x[2])),
                                                                 fixed_params),
                                                        Np=eagle_Np,
                                                        Nmif=eagle_Nmif,
                                                        cooling.fraction.50=eagle_cooling.fraction.50,
                                                        rw.sd=rw.sd(
                                                          beta00=eagle_rw.sd,
                                                          beta01=eagle_rw.sd,
                                                          beta10=eagle_rw.sd,
                                                          beta11=eagle_rw.sd)
                                                   )
                                                 }
  })
},seed=54321,kind="L'Ecuyer")


# log likelihood
stew(file=sprintf("lik_global_eval_covar1_mix-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global_covar1_mix <- foreach(i=1:eagle_Nglobal,
                                       .combine=rbind, .packages='pomp',
                                       .inorder=FALSE,
                                       .options.multicore=list(set.seed=TRUE),
                                       .export=c("eagle_Neval", "eagle_2", "mifs_global_covar1_mix", "eagle_Np")) %dopar% {
                                         evals <- replicate(eagle_Neval,
                                                            logLik(pfilter(eagle_2,
                                                                           params=coef(mifs_global_covar1_mix[[i]]),Np=eagle_Np)))
                                         logmeanexp(evals, se=TRUE)
                                       }
  })
},seed=98765,kind="L'Ecuyer")

