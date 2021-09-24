files.all <- list.files()
files.all <- files.all[grep("result",files.all)]

results.all <- NULL
for(i in files.all){
  load(i)
  results.all <- rbind(results.all,results.run)
}

apply(results.all,2,FUN=function(z) tapply(z, INDEX=results.all[,1],mean))