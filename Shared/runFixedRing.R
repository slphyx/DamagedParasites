library(checkpoint)
checkpoint("2016-07-01", checkpointLocation = tempdir())

library(rstan)
setwd("D:\\works\\ArtemisininResistance\\Stan\\FixedRing")

#initial parameter function
initf<-function(){
  list(N0=12.,mu=10.,sig=6.,pmr=8.,
       #gammaR=6.5,gammaT=6.5,gammaS=6.5,
       ec50R=20.5,ec50T=20.5,ec50S=20.5,
       emaxR=99.,emaxT=99.0,emaxS=99.99,SD=1.5)
}

#Pailin
source("bigcounts.R")
fit <- stan(file = 'SompobFixedRing.stan', data = bigcounts[[1]], iter = 100000, chains = 1, init = initf, algorithm = "NUTS");
saveRDS(fit, file = paste0("fit", ids[[1]], ".rds"))

## run the model
for(indx in 2:length(bigcounts)){
  tmpfit<-stan(fit = fit, data=bigcounts[[indx]],iter = 100000, chains = 1, init = initf,algorithm = "NUTS");
  saveRDS(tmpfit,file = paste0("fit",ids[[indx]],".rds"))
}        
  

#for WangPha
source("bigcountsWP.R")
fit <- stan(file = 'SompobFixedRingWP.stan', data = bigcountsWP[[1]], iter = 100000, chains = 1, init = initf, algorithm = "NUTS");
saveRDS(fit, file = paste0("fit", idsWP[[1]], ".rds"))

## run the model
for (indx in 2:length(bigcountsWP)) {
    tmpfit <- stan(fit = fit, data = bigcountsWP[[indx]], iter = 100000, chains = 1, init = initf, algorithm = "NUTS");
    saveRDS(tmpfit, file = paste0("fit", idsWP[[indx]], ".rds"))
}


