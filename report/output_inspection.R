
library(pomp)
setwd("/home/wangmk/UM/Biostatistics/Courses/STATS531/Time-Series-Final-Project/finished_models")

# load("global_search-v1-2.rda")
# load("lik_global_eval-v1-2.rda")

mif_nocovar_nomix = mifs_global
liks_global_nocovar_nomix = liks_global

save(mif_nocovar_nomix, liks_global_nocovar_nomix, file="result_nocovar_nomix.rda")

# load("global_search-v2-2.rda")
# load("lik_global_eval-v2-2.rda")

# load("global_search-v2-2-allrandom.rda")
# load("lik_global_eval-v2-2-allrandom.rda")

mif_nocovar_mix = mifs_global
liks_global_nocovar_mix = liks_global

save(mif_nocovar_mix, liks_global_nocovar_mix, file="result_nocovar_mix.rda")

coefs = coef(mif_nocovar_mix)

selected_coef = coefs[,38]

write.csv(selected_coef[3:9], file="fixed_params.csv")

load("global_search_nocovar_mix-2.rda")
load("lik_global_eval_nocovar_mix-2.rda")

coefs = coef(mifs_global_nocovar_mix)
selected_coef = coefs[,41]

load("global_search_covar1_mix-2.rda")
load("lik_global_eval_covar1_mix-2.rda")

load("global_search_covar2_mix-2.rda")
load("lik_global_eval_covar2_mix-2.rda")

load("global_search_covar3_mix-2.rda")
load("lik_global_eval_covar3_mix-2.rda")


load("global_search_covar4_mix-2.rda")
load("lik_global_eval_covar4_mix-2.rda")


# simulation code

source("selected_model.R")
selected_params = coef(mifs_global_covar1_mix)[,which.max(liks_global_covar1_mix)]

set.seed(20)

sims <- simulate(eagle_2, params=selected_params, 
                 nsim=4, format="data.frame", include=TRUE, t0=1)
ggplot(sims,mapping=aes(x=timestamp,y=msa,color=.id=="data"))+
  geom_line()+guides(color=FALSE)+facet_grid(rows=vars(.id))

