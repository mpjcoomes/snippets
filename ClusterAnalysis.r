# Cluster Analysis

###Part (a)
#Data Input
library(foreign)
dir <- "Juvenile offenders_Cluster.sav"
x <- read.spss(dir, use.value.labels = F, to.data.frame = TRUE, max.value.labels = Inf,
	trim.factor.names = FALSE, trim_values = TRUE, use.missings = TRUE)
x <- x[complete.cases(x),]
set.seed(43534)
y <- x[sample(1:nrow(x),100,replace=FALSE),]
y['sexabse'] <- NULL
y['Sample_50'] <- NULL
y['ID'] <- NULL
names(y) <- c(	'1.Neglect',
			'2.ICG Control',
			'3.Sibling Crim.',
			'4.Parental Alco.',
			'5.Drug Abuse',
			'6.Offense Type',
			'7.Days to ECD',
			'8.Violence',
			'9.Robber',
			'10.Burglar',
			'11.Other Crime',
			'12.Delinquency',
			'13.Drugs',
			'14.Offense Type 2',
			'15.Age 1st Arrest',
			'16.Violent Crime',
			'17.Total Arrests' )
summary(y) #Checking min/max values for Missing Values codes

#HCA Variances Assumption
options(scipen = 6)
data.frame(sapply(y, var))#variances differ greatly
summary(data.frame(sapply(y, var)))
library(vegan)
y <- decostand(y, method='normalize')
sapply(y, var)

#Correlations
corrs <- round(cor(y),2)
corrs[upper.tri(corrs)] <- ""
corrs <- data.frame(corrs)
names(corrs) <- 1:ncol(corrs)
setwd("C:/Users/user/Desktop")
#write.csv(corrs, "R_output.csv")

#KMO and Bartlett Stats
library(rela)
base <- paf(as.matrix(y), eigcrit=1, convcrit=.001)
KMO <- paste("Kaiser-Meyer-Olkin Measure of Sampling Adequacy =", round(as.numeric(base[5]),3))
df <- (dim(y)[2]^2-dim(y)[2])/2
Bartlett <-  paste("Bartlett’s Test of Sphericity =", round(as.numeric(base[7]),3),
	", df =", df,", Sig.p =", pchisq(as.numeric(base[7]), lower.tail=FALSE, df), sep=" ")
KMO;Bartlett

#Communalities
library(psych)
library(rela)
base <- paf(as.matrix(y), eigcrit=1, convcrit=.001)
Initial <- base$Communalities[,1] #s2 explained in single variable by all others

library(nFactors)
par(family="serif", ps=18, mar=c(5,4.2,0.2,2))
plotnScree(nScree(y), main="", legend=T, xlab="Factor Number") #Scree Plot
unrot <- fa(y, nfactors=5, rotate="none", residuals=TRUE)
Extraction <- unrot$communality #s2 explained in a single variable by extracted factors
Comms <- data.frame(round(cbind(Initial, Extraction ), 3))
Comms <- Comms[order(Comms[,2],decreasing=T),]
Comms

#Pattern Matrix
library(psych)
fit <- fa(y, nfactors=5, n.iter=1, rotate="oblimin", residuals=TRUE, SMC=TRUE, fm="ml")
print(fit$loadings, digits=3, cutoff=.2, sort=TRUE)
y2 <- y[,c(15,1,17,2,10)]
names(y2)

#Euclidean Distances:
distances <- dist(y2, method = "euclidean", diag = FALSE, upper = FALSE)

#HCA
hc <- hclust(distances, method = "ward")
par(family="serif", ps=18)
plot(hc, labels = FALSE)
groups <- cutree(hc, k=4)
rect.hclust(hc, k=4, border="red")
table(groups)
y3 <- data.frame(y2, groups)
HCAmeans <- aggregate(y3,by=list(y3$groups),FUN=mean)


###Part (b)
#Kmeans
library(foreign)
dir <- "E:/Main/University/HMS794 STATISTICAL MARKETING TOOLS/A2/Juvenile offenders_Cluster.sav"
x <- read.spss(dir,use.value.labels=F,to.data.frame=T,max.value.labels=Inf,trim.factor.names=F,
	trim_values=T,use.missings=T)
x <- x[complete.cases(x),]
set.seed(43534)
y <- x
y['sexabse'] <- NULL
y['Sample_50'] <- NULL
y['ID'] <- NULL
names(y) <- c(	'1.Neglect',
			'2.ICG Control',
			'3.Sibling Crim.',
			'4.Parental Alco.',
			'5.Drug Abuse',
			'6.Offense Type',
			'7.Days to ECD',
			'8.Violence',
			'9.Robber',
			'10.Burglar',
			'11.Other Crime',
			'12.Delinquency',
			'13.Drugs',
			'14.Offense Type 2',
			'15.Age 1st Arrest',
			'16.Violent Crime',
			'17.Total Arrests' )
options(scipen = 6)
library(vegan)
y1 <- decostand(y, method='normalize')
y2 <- y1[,c(15,1,17,2,10)]
names(y2)

fit <- kmeans(y2, as.matrix(HCAmeans[,2:6]),iter.max=100, algorithm="Hartigan-Wong")
table(fit$cluster)
ktable <- aggregate(y2,by=list(fit$cluster),FUN=mean)
setwd("C:/Users/user/Desktop")
#write.csv(ktable, "R_output.csv")

library(cluster)
library(fpc)
dev.new(width=10, height=5)
par(family="serif", ps=18, mfrow=c(1,2),adj=0,cex=.6)
clusplot(y2, fit$cluster, color=TRUE, shade=TRUE,  labels=2, lines=0)
plotcluster(y2, fit$cluster)

distances <- dist(y2, method = "euclidean", diag = FALSE, upper = FALSE)
face <- rFace(nrow(y2),dMoNo=2,dNoEy=0,p=2)
cstat <- cluster.stats(distances, fit$cluster, alt.clustering=as.integer(attr(face,"grouping")))
cstat$clus.avg.silwidths
cstat$avg.silwidth

#Cluster Means, original scale
new <- data.frame(y, fit$cluster)
new <- new[,c(15,1,17,2,10,18)]
names(new)
ktable <- aggregate(new,by=list(fit$cluster),FUN=mean)
#write.csv(ktable, "R_output.csv")
