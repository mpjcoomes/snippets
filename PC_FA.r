#____________________Principal Components:______________________________
#Use cor=FALSE to base the principal components on the covariance matrix.
# na.action=na.exclude

data <- USArrests

fit <- princomp(data, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

#KMO and Bartlett:
library(rela)
data <- Seatbelts[,1:7]
paf.belt <- paf(data)
summary(paf.belt)

data <- as.matrix(USArrests)
KMOB <- paf(data)
summary(KMOB)



library(psych)
library(GPArotation)
data <- Harman74.cor$cov
# Varimax Rotated Principal Components, retaining 5 components
#Pairwise deletion of missing data
#"none", "varimax", "quatimax", "promax", "oblimin", "simplimax", "cluster"
fit <- principal(data, nfactors=5, rotate="varimax")
fit
plot(fit,type="lines")

#nFactor can do PCA


#________________principal-axis factor analysis______________
library(psych)
library(GPArotation)
library(MASS)
# runs a principal-axis factor analysis (fm=”pa”) 
# with oblique rotation (rotate=”oblimin”)
pa<-fa(USArrests, nfactors=3, rotate="oblimin", fm="pa")

# creates the diagram with arrows for factor loadings greater than .3
fa.diagram(pa, cut=.3, digits=2, main="Oblique Factor Model")

# runs and creates diagram for a hierarchical factor model
# sl=FALSE overrides default
hier <- omega(USArrests, nfactors=3, sl=FALSE)
omega.diagram(hier, digits=2, main="Hierarchical Factor Model", sl=FALSE)

#runs and creates diagram for bifactor model
# default is Schmid-Leiman bifactor model
bifactor<-omega(USArrests, nfactors=3)
omega.diagram(bifactor, digits=2, main="Bifactor Model")

#____________________Exploratory Factor Analysis:________________

# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors, 
# with varimax rotation and 
#Scores can only be produced if a data matrix is supplied and used: scores="regression"/"Bartlett"
#If entering a covariance matrix, include the option n.obs=

#EFA: USArrests
data <- USArrests
ev <- eigen(cor(USArrests)) #nFactors can bootstrap data.frame for eigens
plot(ev$values, type="both") #scree
fit <- factanal(data, 1, rotation="promax")
loadings(fit)
print(fit, digits=2, cutoff=.3, sort=TRUE)

#EFA: HARMAN74
library(nFactors)
data <- Harman74.cor$cov
plotnScree(nScree(data)) #object is a covariance matrix (or use cor=TRUE), model="components"/"factors"
nBartlett(data, N=145, details=TRUE, alpha=0.05) #correction=T uses a correction for the degree of freedom after the ?rst eigen
str(results); results$detail$lawley.p
summary(nSeScree(eigen(Harman74.cor$cov)$values, model="components", details=TRUE, r2limen=0.75)) #SE&R2 factor diags
fit <- factanal(factors=5, covmat=data, n.obs=145, rotation="varimax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=rownames(data),cex=.7) # add variable names

#FactoMineR:
#use of both quantitative and qualitative variables
#inclusion of supplimentary variables and observations
library(FactoMineR)
fit <- PCA(data)
fit <- CA(data)
plot.PCA(fit)
plot.CA(fit)


###Chapter 2: EFA
#DATA INPUT
library(foreign)
dir <- "E:/SEMdata/Ex2.1_EFA_0.sav"
x <- read.spss(dir, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
	trim.factor.names = FALSE, trim_values = TRUE, use.missings = TRUE)
str(x); head(x); tail(x); summary(x)

#Correlation matrix, KMO & Bartlett’s Sphericity (for cov matrices)
library(rela)
base <- paf(as.matrix(x), eigcrit=1, convcrit=.001)
print(round(base$Correlation,3))
base[c(5,7)] #KMO & Bartlett
pchisq(as.numeric(base[7]), lower.tail=FALSE, 28) #Bartlett Sig

#Scree Plots & Eigen Values
library(nFactors)
nScree(x); plotnScree(nScree(x)) #for CovMat, nScree(cov(x)); plotnScree(nScree(cov(x)))

#Extraction: Maximum Likelihood, based on eigenvalues greater than 1, Direct Oblimin
#Display unrotated and rotated factor solution, loading plots, Regression scores?
#OR for cov matrix: factanal(covmat=cov(x), n.obs=300, factors=3, rotation="none")
library(GPArotation)
unrot <- factanal(x, factors=2, rotation="none")
print(loadings(unrot), digits=3, cutoff=0, sort=FALSE) #unrotated loadings

#Communalities
Initial <- paf(as.matrix(x), eigcrit=1, convcrit=.001)$Communalities[,1]
Extraction <- 1-factanal(x, factors=2, rotation="oblimin")$uniquenesses
print("Communalities"); round(cbind(Initial, Extraction ), 3) #Communalities

#Total Variance Explained
Total <- eigen(cor(x))$values #or paf(as.matrix(x), eigcrit=1, convcrit=.001)$Eigenvalues[1:8] #eigen method 2
PercentVariance <- 100*eigen(cor(x))$values/ sum(unlist(eigen(cor(x))[1]))
CumulativePercent <- 100*cumsum(eigen(cor(x))$values/ sum(unlist(eigen(cor(x))[1])))
print("Total Variance Explained"); round(cbind(Total, PercentVariance, CumulativePercent), 3)
unrot <- factanal(x, factors=2, rotation="none")
print(loadings(unrot), digits=3, cutoff=0, sort=FALSE) #unrotated loadings, factor matrix

#Reproduced Correlations
obrot <- factanal(x, factors=2, scores="regression", rotation="oblimin")
print(loadings(obrot), digits=3, cutoff=0.3, sort=TRUE)
str(obrot)
obrot[c(1,6,8,9,11:13)]






obrot$scores #Regression scores







#IF data has NAs:
y <- x
y[3,1] <- NA; head(y)
factanal(~y$Beh01+y$Beh02+y$Beh03+y$Beh10+y$Beh14+y$Beh05+y$Beh12+y$Beh13,
	factors=3,scores ="regression", rotation="oblimin", na.action(na.omit))


#____________________Confirmatory Factor Analysis:________________
library(lavaan)
data <- Harman74.cor$cov #covariance matrix: cov( na.om(dataframe) )

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
fit <- cfa(HS.model, data=HolzingerSwineford1939)
summary(fit, fit.measures=TRUE)



library(lavaan)
model <- ‘
# latent variable definitions
factor_1 =~ y1 + y2 + y3 + y4
factor_2 =~ y5 + y6 + y7 + y8
# covariance between factor_1 and factor_2
factor_1 ~~ factor_2
# residual covariances
y1 ~~ y5
‘

data <- as.matrix(USArrests)
model <- '	factor_1 =~ y1 + y2 + y3 + y4
		factor_2 =~ y5 + y6 + y7 + y8
		factor_1 ~~ factor_2
		y1 ~~ y5 '
fit <- cfa(model, data=data)
summary(fit)

#as a bonus lavaan has a dedicated function that lets you run a multiple-group
#confirmatory factor analysis to test for measurement invariance. Something that
#took me a while in AMOS.

measurement.invariance(model, data=ex_data, group =”school” )




## Example with missing data
## use package missMDA
## Not run: 
require(missMDA)
data(orange)
nb <- estim_ncpPCA(orange,ncp.min=0,ncp.max=5,method.cv="Kfold",nbsim=50)
imputed <- imputePCA(orange,ncp=nb$ncp)
res.pca <- PCA(imputed$completeObs)

