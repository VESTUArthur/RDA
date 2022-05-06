source("JTK_CYCLEv3.1.R")

options(stringsAsFactors=FALSE)

filename <- "cene.csv"
descriptor_filename <- "cene_descriptor.txt"
fixed_period <- 1


data <- read.delim(filename)
x = read.csv(descriptor_filename, header = FALSE, sep = "\t")

Names <- data[,1]
data <- data[,-1]


timepoints = x[1,1]
reps = x[1,2]
dt= x[1,3]


#jtkdist(ncol(data), reps=4) # number of replicates
jtkdist(timepoints=timepoints, reps=reps)

if (fixed_period == 1) {
  periods <- (24%/%dt):(24%/%dt)  # looking for rhythms of exactly 24 hours
} else {
  periods <- (18%/%dt):min(36%/%dt, timepoints-1)   # looking for rhythms between 18-36 hours
}


jtk.init(periods,dt); # 3 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
	jtkx(z)
  c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(Names,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

results_all <- data.frame(matrix(ncol = 6, nrow = 0))
results_all <- rbind(results_all, results[1:6])

if (fixed_period == 1) {
  write.table(results_all,file=paste("JTK_on_file",strsplit(filename,"\\.")[[1]][1],"fix_per","txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")
} else {
  write.table(results_all,file=paste("JTK_on_file",strsplit(filename,"\\.")[[1]][1],"arb_per","txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")
}


