library(rstan)
library(loo)

extractPars<-function(rdsfile,parms,id){
  #parms must be input as a vector type
  parlen<-length(parms)
  fit<-readRDS(rdsfile)
  extr<-extract(fit,pars=parms)
  for(i in parms){
    cat(paste0("writing ",i,"-",id,".csv"),"\n")
    write.csv(eval(parse(text = paste0("extr$",i))), file=paste0(i,"-",id,".csv"))
  }
}
