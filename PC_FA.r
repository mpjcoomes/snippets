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
fit <- PCA(x)
fit <- CA(x)
plot.PCA(fit)
plot.CA(x)

#IF data has NAs:
y <- x
y[3,1] <- NA; head(y)
factanal(~y$Beh01+y$Beh02+y$Beh03+y$Beh10+y$Beh14+y$Beh05+y$Beh12+y$Beh13,
	factors=3,scores ="regression", rotation="oblimin", na.action(na.omit))



###########################################################################################
#						Chapter 2: EFA
###########################################################################################
#Input data: Ex2.1

library(foreign)
dir <- "E:/SEMdata/Ex2.1_EFA_0.sav"
x <- read.spss(dir, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
	trim.factor.names = FALSE, trim_values = TRUE, use.missings = TRUE)

#__________________________________________________________________________________________
#Correlation Matrix, KMO & Bartlett’s Sphericity for cov matrices:

library(rela)
base <- paf(as.matrix(x), eigcrit=1, convcrit=.001)
Corr <- as.matrix(round(base$Correlation,3))
Corr[upper.tri(Corr, diag = TRUE)] <- ""
Corr <- data.frame(Corr)
KMO <- paste("Kaiser-Meyer-Olkin Measure of Sampling Adequacy =", round(as.numeric(base[5]),3))
df <- (dim(x)[2]^2-dim(x)[2])/2
Bartlett <-  paste("Bartlett’s Test of Sphericity =", round(as.numeric(base[7]),3),
	", df =", df,", Sig.p =", pchisq(as.numeric(base[7]), lower.tail=FALSE, df), sep=" ")
Corr
KMO;Bartlett

#__________________________________________________________________________________________
#Communalities

unrot <- factanal(x, factors=2, rotation="none")

Initial <- base$Communalities[,1] 
Extraction <- 1-unrot$uniquenesses #same as unrotated uniquenesses
#Initial: the variance explained in a single variable by all other variables
#Extraction: variance explained in a single variable by the extracted factors
print("Communalities (aka. R^2): Factor Extraction"); round(cbind(Initial, Extraction ), 3)

#__________________________________________________________________________________________
#Factor Number Selection: Total Variance Explained Table, Scree Plot

#Total Variance Explained: Eigen Values
Total <- eigen(cor(x))$values
PercentVariance <- 100*Total/sum(Total)
CumulativePercent <- 100*cumsum(Total/sum(Total))
TotalSSLoadings <- c(sum(unrot$loadings[,1]^2), sum(unrot$loadings[,2]^2),NA,NA,NA,NA,NA,NA)
SSPercentVariance <- 100*c(sum(unrot$loadings[,1]^2)/length(Total),
				sum(unrot$loadings[,2]^2)/length(Total),NA,NA,NA,NA,NA,NA)
CumSSPercentV <- c(SSPercentVariance[1], SSPercentVariance[1]+SSPercentVariance[2],NA,NA,NA,NA,NA,NA)
TVAtable <- cbind(Total, PercentVariance, CumulativePercent,TotalSSLoadings,SSPercentVariance,CumSSPercentV )
rownames(TVAtable) <- paste("Factor", 1:length(Total))
TVAtable <- round(TVAtable,3)
print("Total Variance Explained"); TVAtable
#SPSS also shows SS Loadings for extraction, see p#Structure, and "colSums((p$Structure)^2)"

#Scree Plot
library(nFactors)
nScree(x); plotnScree(nScree(x)) #nScree(cov(x)) for Covariance Matrices

#__________________________________________________________________________________________
#Print Unrotated Factor Loadings & (MLE extraction, Direct Oblimin, eigenvalues>1)

library(GPArotation)
unrot <- factanal(x, factors=2, rotation="none")
obrot <- factanal(x, factors=2, scores="regression", rotation="oblimin")
#for cov matrix: factanal(covmat=cov(x), n.obs=300, factors=3, rotation="none")
print("Unrotated (orthogonal) Factor Loadings"); print(loadings(unrot), digits=3, cutoff=0, sort=FALSE)
#Note. Unrotated factanal() loading will be identical to an unrotated fa() structure matrix.

#__________________________________________________________________________________________
#Factor Matrix Table plus communalities (unrotated solution)

Factor1 <- loadings(unrot)[,1]
Factor2 <- loadings(unrot)[,2]
Ini2plusExt2 <- loadings(unrot)[,1]^2 +loadings(unrot)[,2]^2
FactorMatrix <- round(cbind(Factor1, Factor2, Ini2plusExt2),3)
colnames(FactorMatrix) <- c("Factor 1", "Factor 2", "(Ini^2)+(Ext^2)")
print("Factor Matrix and Communalities"); FactorMatrix 

#__________________________________________________________________________________________
#Reproduced Correlations, communalities and residuals from unrotated solution

ReproducedCorrs <-(unrot$loadings[,1] %o% unrot$loadings[,1]) + (unrot$loadings[,2] %o% unrot$loadings[,2])
Residuals <- unrot$correlation - ReproducedCorrs #Also: residuals from fa()
diag(ReproducedCorrs) <- (unrot$loadings[,1]^2) + (unrot$loadings[,2]^2)
ReproducedCorrs[upper.tri(ReproducedCorrs)] <- Residuals[upper.tri(Residuals)]
print("TR=Residuals; BL=Reproduced Correlations; DG=Reproduced Communalities"); print(round(ReproducedCorrs,3))

#__________________________________________________________________________________________
#Goodness-of-fit test: MLE extraction only

paste("H0:Residuals = 0 (good model), H1: Residuals > 0 (poor model)")
paste("Goodness-of-fit test:","Chi-Square =",round(unrot$STATISTIC,3),
	", df =", unrot$dof, ", Sig.",round(unrot$PVAL,3))

#__________________________________________________________________________________________
#Loading Plots

#Unrotated
loadings <- unrot$loadings[,1:2]
plot(loadings,type="n",xlim=c(-1, 1), ylim=c(-1, 1), main="Loadings Plot: Unrotated")
text(loadings,labels=rownames(loadings),cex=.7)
abline(0,0,v=0)
#rotmat <- varimax(unrot$loadings) #sample orthogonal rotation
#abline(0,rotmat$rotmat[2,1]/rotmat$rotmat[1,1],lty=2)
#abline(0,rotmat$rotmat[2,2]/rotmat$rotmat[1,2],lty=2)

#Rotated
loadings <- obrot$loadings[,1:2]
plot(loadings,type="n",xlim=c(-1, 1), ylim=c(-1, 1), main="Loadings Plot: Oblimin Rotation")
text(loadings,labels=rownames(loadings),cex=.7)

#__________________________________________________________________________________________
#Pattern & Structure Matrix, Factor Correlation Matrix, Score Regression Weights
library(psych)
p <- fa(x, nfactors=2, n.iter=1, rotate="oblimin", scores="regression", residuals=TRUE,
SMC=TRUE, min.err=0.001, max.iter=50, warnings=TRUE, fm="mle",alpha=.1,p=.05,oblique.scores=FALSE)
print("Pattern Matrix");obrot$loadings;p$loadings
print("Structure Matrix");p$Structure #always check that p$loadings is same as for factanal
print("Factor Correlation Matrix");round(p$Phi,3) #Only for Oblique rotations where factors may correlate
print("Factor score coefficient matrix");print(p$weights,cutoff=0) #can be used to compute an estimated factor score on each factor for each student

#Manual prediction based on Factor score coefficient matrix coefficients:
F1 <- x; i <- 1;
for(i in i:length(p$weights[,1])) { F1[,i] <- F1[,i]*p$weights[i,1] }
F1[,c(length(p$weights[,1])+1)] <- rowSums(F1)
F2 <- x; i <- 1;
for(i in i:length(p$weights[,2])) { F2[,i] <- F2[,i]*p$weights[i,2] }
F2[,c(length(p$weights[,2])+1)] <- rowSums(F2)
#print("Estimated Factor 1 Scores"); round(F1[,9],3)
#print("Estimated Factor 2 Scores"); round(F2[,9],3)
#plot(F1[,9],F2[,9]) #Both plots are the same, as P$scores are normalised F1&F2 scores (e.g. [F1-mean]/SD)
#plot(p$scores[,1],p$scores[,2])

#Factor score covariance matrix
var(F1[,9]);var(p$scores[,1]) #Different from SPSS: How is it calculating this?
var(F2[,9]);var(p$scores[,2])
cov(F1[,9],F2[,9]); cov(p$scores[,1], p$scores[,2])

#Final Model
p
fa.diagram(p)
fa.graph(p)


#_________________________________________________________________________________________

###########################################################################################
#						OpenMx								#
###########################################################################################
#Install and Load OpenMx
#repos <- c('http://openmx.psyc.virginia.edu/packages/')
#install.packages(pkgs=c('OpenMx'), repos=repos)
library(OpenMx)

names(x) <- c("id","v1","v2","v3","v4",
"nv1","nv2","nv3","nv4","momed",
"wisc1","wisc2","wisc4","wisc6","cte")
myCov <- cov(x[,c("wisc1", "wisc6")])
myMeans <- mean(x[,c("wisc1", "wisc6")])
n <- length(x$wisc1)

##	Build Model
model_1 <-	mxModel(		#create new model
		"My Title",		#title
		mxData(		#4 data args: 1. observed= , 2. type= , 3. numObs= , 4. means= 
			observed=myCov	,			#a data frame or matrix
			type=c("cov" "raw" "cor" "SSCP") ,	#type of data input
			means=myMeans	,			#for non-raw data, may specify a means vector
			numObs=n	)				#for non-raw data, must specify observations
			)

##	Prepare OpenMx for our paths
model_1 <-	mxModel("My Title", mxData(myCov, "cov", myMeans, n)  , #build model
			type="RAM"	,	#tells model paths used, and RAM(ML) objecitve function, not matrices
			manifestVars=c("wisc1", "wisc6")	,	#lists (vectors) of variables in our model
			latentVars="diff"	)				#lists (vectors) of variables in our model
##	Add paths
model_1 <-	mxModel("My Title", mxData(myCov, "cov", myMeans, n),	type="RAM",	manifestVars=c("wisc1", "wisc6"),
			latentVars="diff"	,					#build model and prepare for paths
			mxPath(from=c("wisc1", "diff"), to="wisc6",	#'from' the regressor 'to' the predicted
			connect=c("single", "all.pairs", "unique.pairs", "all.bivariate", "unique.bivariate") #?mxPath
			arrows=1,		#1 single-headed indicates regression, 2 double-head for cov or variance
#regression of the ‘to’ variable on the ‘from’ variable (arrow points at ‘to’ variable)
#multiple paths in one mxPath function, arrows may take a vector of 1s and 2s (app. to set of created paths)
			free=F,		#boolean vectore c("T", "T", "F"). Indicates whether paths are free or fixed. 
			values=1, 		#numeric vector c(1, 2, 1, 1). The starting values of the parameters.
#‘values’ gives a starting value for estimation
			labels=NA)		#character vector. The names of the paths.
			lbound = NA, ubound = NA # specify lower and upper bounds for the created paths
			)


model_1 <-	mxModel("My Title", mxData(x, "raw"), type="RAM",
			manifestVars=c("Number", "Meas"), latentVars="diff",
		mxPath(from=c("Number", "diff"), to="Meas",
			connect="all.pairs", arrows=1, free=T, values=1, labels=NA, lbound = NA, ubound = NA)
		mxPath(from="one", to=c("Number", "diff"),
			connect="all.pairs", arrows=1, free=F, values=1, labels=NA, lbound = NA, ubound = NA)
			)
results <- mxRun(model_1)
summary(results)

manifests <- names(x) #vector of manifest variables
latents <- c("Please") #vector of latent variables

factorModel <- mxModel(name="One Factor",	type="RAM",
	manifestVars = manifests,
	latentVars = latents,
	mxPath(from=latents, to=manifests),
	mxPath(from=manifests, arrows=2),
	mxPath(from=latents, arrows=2, free=F, values=1.0),
	mxData(cov(x), type="cov", numObs=300)
	)

#is "free=T" the equivalent of AMOS path=1?
#mxData(x, "raw"),

#Run Model
factorFit <- mxRun(factorModel)
summary(factorFit)

#Diagram with Parameter Estimates
setwd("E:/OpenMxDotfiles")
a <- omxGraphviz(factorFit, dotFilename="")
p <- summary(factorFit)$parameters
i <- 1
for(i in i:length(p[,1]))
	{
	a <- sub(
		paste(p[i,4], " -> ", p[i,3], "[",sep=""),
		paste(p[i,4], " -> ", p[i,3], "[label=", round(p[i,5],3), ",", sep="") ,
		a, fixed=TRUE)
	}
a <- sub("Arial","Times-Roman",a, fixed=TRUE)
write(a,file="SEMdotcode.dot")
#shell("notepad E:/OpenMxDotfiles/SEMdotcode.dot")
shell("dot -Tjpeg -o SEMplot.jpeg SEMdotcode.dot")
shell("E:/OpenMxDotfiles/SEMplot.jpeg")
p[,1:8]

#headlabel = ""
#taillabel = "1"
#label = "gamma? 2.698"
#labeldistance=5, taillabel="",
#arrowhead=integer, arrowtail=integer
#orientation=landscape, 


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

