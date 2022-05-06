library(MetaCycle)
annot <- cycMouseLiverRNA[1]
data <- cycMouseLiverRNA
print(data)
data$geneName <- NULL
head(data)
colnames(data) <- sub("..", "", colnames(data))
write.csv(data, "C:/Users/arthu/Documents/GitHub/RDA/df.csv")
meta2d("C:/Users/arthu/Documents/GitHub/RDA/df.csv",
    timepoints = "line1", filestyle = "csv")

library(tidyverse)
df <- read.csv(file = "df.csv")
head(df)
summary(df)

citation("Biobase")
library(rain)
#deltat : sampling interval
#period : to search for
data(menetRNASeqMouseLiver)
data <- menetRNASeqMouseLiver
rownames(data) <- NULL
write.csv(data, "data.csv")
data <- read.csv("data.csv")
head(data)
data$X <- NULL
res <- rain(t(data), deltat = 4, nr.series = 2, period = 24, verbose = TRUE)
write.csv(res, "RAINresult_")
plot(res$pVal)
print(df)




data <- read.csv("data.csv", row.names = 1)
head(data)
sampleRate <- 4
nbReplicate <- 2
period <- 24
res <- rain(t(data), deltat = sampleRate, 
    nr.series = nbReplicate, period = period, verbose = TRUE)
write.csv(res, "RAINresult_data.csv")

dir.create("rainout", showWarnings = FALSE)



data <- menetRNASeqMouseLiver
rownames(data) <- NULL
write.csv(data, "{filename}")



meta2d("C:/Users/arthu/Documents/GitHub/RDA/formated_spreadsheet.txt",
    timepoints = "line1", filestyle = "txt")
