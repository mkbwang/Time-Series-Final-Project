library(pomp)

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


# gamma distribution
eagle_rmeasure <- "
  if (S==0){
    msa = rgamma(shape0, scale0);
  } else{
    msa = rgamma(shape1, scale1);
  }
"

eagle_dmeasure <- "
  if (S==0){
    lik = dgamma(msa, shape0, scale0, give_log);
  } else{
    lik = dgamma(msa, shape1, scale1, give_log);
  }
"

eagle_statenames <- c("S")

eagle_paramnames <- c("p0", "p1", "shape0", "scale0", "shape1", "scale1")

eagle_data <- read.csv("https://raw.githubusercontent.com/skybullbobby/Time-Series-Final-Project/master/eagle_421.csv")

eagle_0 <- pomp(
  data=subset(eagle_data,select=c(timestamp, msa)),
  times="timestamp",
  t0=1,
  rprocess=discrete_time(Csnippet(eagle_rprocess),delta.t=1),
  rmeasure=Csnippet(eagle_rmeasure),
  dmeasure=Csnippet(eagle_dmeasure),
  partrans=parameter_trans(
    log=c("shape0", "scale0", "shape1", "scale1"),
    logit=c("p0", "p1")),
  statenames=eagle_statenames,
  paramnames=eagle_paramnames,
  rinit=Csnippet(eagle_rinit)
)

