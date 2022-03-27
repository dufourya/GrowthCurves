#load STAN library
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#read data from file, time = first column, OD = remaining colums
data <- read.table('data.txt', header = T)
#sort data by replicate and time
data <- data[order(data[,'Replicate'], data[,'Time']),]
#calculate the number of data points for each replicate
replicate_size <- array(as.vector(table(data[,1])))
#calculate minumum estimated growth rate
maxr <- (max(data[,'Time'])-min(data[,'Time']))/log(max(data[,'OD'])/min(data[data[,'OD']>0,'OD']));
#format data in a list structure for STAN
data_stan <- list(N = nrow(data), K = length(replicate_size), idx = replicate_size, time = data[,'Time'], OD = data[,'OD'],maxr = maxr)

#compile code
stan_model <- stan(data = data_stan, file = 'logistic_growth.stan',iter = 1, chain = 1)
#stan_model <- stan(data = data_stan, file = 'logistic_growth_nodead.stan',iter = 1, chain = 1)
#stan_model <- stan(data = data_stan, file ='logistic_growth_ode.stan',iter=1, chain=1)
#stan_model <- stan(data = data_stan, file ='gompertz_growth.stan',iter=1, chain=1)

#run parameter sampling
fit <- stan(data = data_stan, fit = stan_model ,iter = 5500, warmup = 500, thin = 20, chain = 4)

#extract data from fit
out <- extract(fit)

#parameters:
#A: maximum OD
#L: number of viable cells in the inoculum
#B: number of dead cells in the inoculum
#mu: growth peroid

#calculate the mean of the parameter posterior distributions
means <- get_posterior_mean(fit,pars = c('A', 'L', 'mu', 'B', 'A_rep', 'L_rep', 'mu_rep', 'B'))
g <- log(2)*means['mu',5]

#plot data and model fit
pos <- 1;
i <- 1;
data_color <- rainbow(data_stan$K);
idx <- data[,1] == data[pos,1];
plot(data[idx,2],data[idx,3],col = data_color[i], pch = 20)
points(colMeans(out$time_pred), colMeans(out$OD_pred[,,i]),'l',col = data_color[i])
pos <- pos + data_stan$idx[i];
if (data_stan$K > 1) {
  for (i in 2:data_stan$K) {
    idx <- data[,1] == data[pos,1];
    points(data[idx,2],data[idx,3],col = data_color[i],pch = 20);
    points(colMeans(out$time_pred), colMeans(out$OD_pred[,,i]),'l',col = data_color[i])  
    pos <- pos + data_stan$idx[i];
  }
}
points(colMeans(out$time_pred), colMeans(out$OD_pred_all),'l',col = 'black',lwd = 5)
title(paste('generation time = ',round(g),' minutes'))

#plot diagnostic traces to check sampling stability
traceplot(fit,pars = c('A','L','mu','B'));

#plot parameter sampling summary
plot(fit, pars = c('A', 'L', 'mu', 'B','A_rep', 'L_rep', 'mu_rep','B_rep'))

#scatter plot of parameters to check for correlations
pairs(out[c('A', 'L', 'mu', 'B')])

#histograms of parameter estimations
hist(out$mu,50)
hist(out$A,50)
hist(out$L,50)
hist(out$B,50)

