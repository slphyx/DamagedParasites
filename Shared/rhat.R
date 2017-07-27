
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

RhatSummary<-function(stanrdsfilename,outcsvname){
  tmp<-readRDS(stanrdsfilename)
  tmp.sum<-summary(tmp)$summary
  
  write.csv(tmp.sum[,10],file = outcsvname)
}

