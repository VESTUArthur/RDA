library(rain)

### 
# Options
###


filename <- "U2OS_lipidi.xlsx"
diff <- 1
trim_data <- 0
fixed_period <- 0
dt <- 1
independent <- FALSE

################
################
################
################
################


options(stringsAsFactors=FALSE)

library("openxlsx")

sheetNames <- getSheetNames(filename)
n_sheets <- length(sheetNames)


results_all <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(results_all) <- c("Names", "BH.Q", "ADJ.P", "PER", "LAG", "AMP")




for (i in 1:n_sheets) {
#for (i in 8:8) {
  # get current sheet
  mydf <- read.xlsx(filename, sheet = i, colNames = FALSE)
  sheetName <- sheetNames[i]
  x <- as.vector(t(mydf["X1"]))
  y <- as.vector(t(mydf["X2"]))

  if (trim_data == 1){
    x_start = x[1]
    x = x[-1]
    x_end = y[1]
    y = y[-1]
    if (x_start != -1)
    {
      valid_idx = which(TRUE == (x>x_start & x<x_end))
      x = x[valid_idx];
      y = y[valid_idx];
    }
  }
  
  if (diff == 1){
    for (j in 1:(length(x)-1)) {
      if (x[j+1] >= x[j]) {
        y[j] = y[j+1]-y[j]
      }
      else {
        y[j] = NA
      }
    }
    y[length(y)] = NA
  }
    
  idxs <- c(1, which(TRUE == (c(0, x) > c(x,Inf))),length(x)+1)
  
  reps <- length(idxs)-1
  timepoints = max(x)%/%dt+1
  
  data <- rep(NA, timepoints * reps)
  columns <- rep(NA, timepoints * reps)
  
 
  for (t in 1:timepoints){
    for (r in 1:reps){
      i1 <- idxs[r]
      i2 <- idxs[r+1]-1
      
      slicex <- x[i1:i2]
      slicey <- y[i1:i2]
      
      if (((t-1)*dt) %in% slicex) {
        idxs2 = which(((t-1)*dt) == slicex)
        data[(t-1)*reps + r] <- median(slicey[idxs2], na.rm = TRUE)
        
      }
      
      #columns[(t-1)*reps + r] <- paste("ZT_", toString((t-1)*dt), "_", toString(r), sep="")
      columns[(t-1)*reps + r] <- paste((t-1)*dt)
      
    }
  }
  
  
  
  df <- as.data.frame(t(data))
  colnames(df) <- columns
  #df <-cbind(gene = sheetName, df)
  #names(df)[1] <- ""
  
  #https://stat.ethz.ch/R-manual/R-devel/library/base/html/row.names.html
  row.names(df) <- sheetName
  
  
  #colnames(df) <- NULL
  
  ####
  # rain
  if (independent == TRUE) {
    method = 'independent'
  }
  else{
    method = 'longitudinal'
  }
  
  osc <- rain(t(df), deltat = dt, period = 24,  nr.series = reps, verbose=TRUE, method = method, na.rm=TRUE)
  
  if (fixed_period == 1) {
    #osc <- rain(t(df), deltat = dt, period = 24,  nr.series = reps, peak.border = c(0.3, 0.7), verbose=TRUE, na.rm=TRUE)
    osc <- rain(t(df), deltat = dt, period = 24,  nr.series = reps, verbose=TRUE, method = method, na.rm=TRUE)
  }
  else {
    osc <- rain(t(df), deltat = dt, period.delta = 8, period = 24,  nr.series = reps, peak.border = c(0.3, 0.7), verbose=TRUE, na.rm=TRUE)
  }

  
  
  
  
  
  
  
  project <- sheetName
  #data <- read.delim("Huh7_data.txt")
  
  Names <- df[,1]
  data <- df[,-1]

  #jtkdist(ncol(data), reps=4) # number of replicates
  jtkdist(timepoints=timepoints, reps=reps) # number of timepoints per replication, number of replications
  
  if (fixed_period == 1) {
    periods <- (24%/%dt):(24%/%dt)  # looking for rhythms of exactly 24 hours
  }
  else {
    periods <- (18%/%dt):min(50%/%dt, timepoints-1)   # looking for rhythms between 18-50 hours
  }
  
  
  jtk.init(periods, dt); # dt is the number of hours between time points

  cat("JTK analysis of", sheetName,"started on",date(),"\n")
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

  #write.table(results[1:6],file=paste("JTK",project,"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")
  
  results_all <- rbind(results_all, results[1:6])
  
  cat("JTK analysis of", sheetName,"finished on",date(),"\n")
  
  
 
}

p <- as.vector(t(results_all["ADJ.P"]))
bf <- p.adjust(p, method = "bonferroni")
bh <- p.adjust(p, method = "BH")

results_all["ADJ.P"] <- bf
results_all["BH.Q"] <- bh

write.table(results_all[1:6],file=paste("JTK",strsplit(filename,"\\.")[[1]][1],"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")




