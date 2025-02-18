library(psych)
library(GPArotation)
library(paran)
library(GGally)
library(plyr)

# Creating dataframe from relevant tests in the Human Connectome Project 1200 Subjects Data Release (public domain)
HCP.Data<-read.csv("file_name.csv",header=TRUE) # replace file name your download from HCP website
flanker<-HCP.Data[,"Flanker_Unadj"]
sequence<-HCP.Data[,"PicSeq_Unadj"]
wmTask<-HCP.Data[,"WM_Task_Acc"]
cardSort<-HCP.Data[,"CardSort_Unadj"]
listSort<-HCP.Data[,"ListSort_Unadj"]
wordMem<-HCP.Data[,"IWRD_TOT"]
wordMemReac<-HCP.Data[,"IWRD_RTC"]
orientation<-HCP.Data[,"VSPLOT_TC"]
HCP<-data.frame(Flanker=flanker,PictureSequence=sequence,WorkingMem=wmTask,CardSorting=cardSort,ListSorting=listSort,WordMemory=wordMem,SpatialOrient=orientation,WordMemoryReact=wordMemReac)

# Removing missing data
OmitHCP<-na.omit(HCP)

# Examining correlations
CorHCP<-cor(OmitHCP)
CorHCP
ggcorr(OmitHCP)

# Examining factorability
cortest.bartlett(OmitHCP)
KMO(OmitHCP)

# Checking factors to retain
pc<-principal(OmitHCP,nfactors=8,rotate="none")
plot(pc$values,type="b")
paran(OmitHCP, iterations = 240, centile = 0, quietly = FALSE, status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, file = "", width = 640, height = 640, grdevice = "png", seed = 0)

# Examining residuals with 3 factors
pcunrot<-principal(OmitHCP,nfactors=3,rotate="none")
residuals<-factor.residuals(CorHCP,pcunrot$loadings)
residuals<-as.matrix(residuals[upper.tri(residuals)])
count(abs(residuals)>.05)
sqrt(mean(residuals^2))
hist(residuals)

# Running the analysis with 2 factors
principal(OmitHCP,nfactors=2,rotate="oblimin",scores=TRUE)
principal(OmitHCP,nfactors=2,rotate="varimax",scores=TRUE)

# Running the analysis with 3 factors
pc1<-principal(OmitHCP,nfactors=3,rotate="oblimin",scores=TRUE)
pc1
principal(OmitHCP,nfactors=3,rotate="varimax",scores=TRUE)

# Examining factor loadings
scores<-pc1$scores
plot(scores)
boxplot(scores)

# Running factor analyses for comparison
factanal(OmitHCP,3,"oblimin")
factanal(OmitHCP,3,"varimax")