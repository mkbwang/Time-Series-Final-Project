library(ggplot2)
library(lubridate)
library(fitdistrplus)


# plot the data
data = read.csv("https://raw.githubusercontent.com/skybullbobby/Time-Series-Final-Project/master/eagle_421.csv")
ggplot(data, aes(x=timepoint, y=msa))+geom_line()+xlab("Timepoint")+ylab("MSA")

# fit a gamma distribution
fit.gamma <- fitdist(data$msa, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)

