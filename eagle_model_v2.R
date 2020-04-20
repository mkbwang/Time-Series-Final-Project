library(pomp)
library(fitdistrplus)


setwd(getwd())

# the state is binary, either active or inactive
eagle_rinit <- "S=1;"

# state transition and set the corresponding parameters for gamma distribution
eagle_rprocess <- "
  if(S==0){
    S=rbinom(1,p0);
  } else{
    S=rbinom(1,p1);
  }
"


# gamma distribution, mixed gamma distribution for slow state
eagle_rmeasure <- "
  if (S==0){
    double msa1 = rgamma(shape01, scale01);
    double msa2 = rgamma(shape02, scale02);
    double msas[] = {msa1, msa2};
    int id = rbinom(1, pmix);
    msa = msas[id];
  } else{
    msa = rgamma(shape1, scale1);
  }
"

eagle_dmeasure <- "
  if (S==0){
    lik = (1-pmix)*dgamma(msa, shape01, scale01, give_log) + pmix*dgamma(msa, shape02, scale02, give_log);
  } else{
    lik = dgamma(msa, shape1, scale1, give_log);
  }
"

eagle_statenames <- c("S")

eagle_paramnames <- c("p0", "p1", "pmix", "shape01", "scale01", 
                      "shape02", "scale02", "shape1", "scale1")

eagle_data <- read.csv("eagle_421.csv")

eagle_1 <- pomp(
  data=subset(eagle_data,select=c(timestamp, msa)),
  times="timestamp",
  t0=1,
  rprocess=discrete_time(Csnippet(eagle_rprocess),delta.t=1),
  rmeasure=Csnippet(eagle_rmeasure),
  dmeasure=Csnippet(eagle_dmeasure),
  partrans=parameter_trans(
    log=c("shape01", "scale01",
          "shape02", "scale02",
          "shape1", "scale1"),
    logit=c("p0", "p1", "pmix")),
  statenames=eagle_statenames,
  paramnames=eagle_paramnames,
  rinit=Csnippet(eagle_rinit)
)

