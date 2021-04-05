# Modelling anonymised power consumption data while accounting for weather induced variance
setwd('~/')
set.seed(2785621)
options(scipen=2)
library(ggplot2)
library(scales)
library(grid)
library(knitr)
apa_style <- theme_bw(base_size=10,base_family='sans') %+replace%
  theme(axis.text=element_text(size=10),
        legend.text=element_text(size=10),
        axis.ticks=element_line(colour='black'),
        axis.line=element_line(),
        axis.title.x=element_text(vjust=-.1),
        axis.title.y=element_text(vjust=1.3, angle=90),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.margin=unit(4,'mm'),
        plot.background=element_blank(),
        legend.background=element_blank(),
        legend.key=element_rect(fill=NA,colour=NA),
        legend.key.width=unit(12.5,'mm'),
        axis.ticks.margin=unit(1.5,'mm'),
        strip.background=element_rect(fill='white'),
        legend.title=element_text(face=NULL))

# Import & Format Data Set
x <- read.csv('RawConsumptionData.csv')
x$row <- 1:nrow(x)
x$ConsumptionDate <- strptime(x$Consumption.Date,'%d/%m/%Y')
x$InterventionDate <- strptime(x$Intervention.date,'%d/%m/%Y')
x$InterventionType <- x$Intervention
z <- x[x$InterventionType=='C',]
z$InterventionDate <- ave(z$ConsumptionDate,z$NMI,FUN=median)
x[x$InterventionType=='C',]$InterventionDate <- z$InterventionDate
x$Intervention <- as.numeric(x$InterventionDate<x$ConsumptionDate)
x$LkWh <- log(x$kWh+10)
x$B1LkWh <- ave(x$LkWh,x$NMI,FUN=function(x)c(NA,x[1:(length(x)-1)]))
x$B2LkWh <- ave(x$LkWh,x$NMI,FUN=function(x)c(NA,NA,x[1:(length(x)-2)]))
x$d <- substr(weekdays(x$ConsumptionDate),1,2)
x <- cbind(x,model.matrix(~d-1,data=x))

# Import BOM Data Sets & Merge http://www.bom.gov.au/climate/data/index.shtml
y <- data.frame(ConsumptionDate=seq(as.Date('2012/1/1'),length.out=900,by='days'),
                MaxC=subset(read.csv('IDCJAC0010_086282_1800_Data.csv'),Year>=2012)[1:900,6],
                MinC=subset(read.csv('IDCJAC0011_086282_1800_Data.csv'),Year>=2012)[1:900,6],
                Rain=subset(read.csv('IDCJAC0009_086282_1800_Data.csv'),Year>=2012)[1:900,6],
                SunE=subset(read.csv('IDCJAC0016_086282_1800_Data.csv'),Year>=2012)[1:900,6])
y$ConsumptionDate <- strptime(y$ConsumptionDate, '%Y-%m-%d')
x <- merge(x,y,sort=F,all.x=T)

# Create & Export Modelling Data
x <- x[order(x$row),]
row.names(x) <- x$row
data1 <- x[,c('NMI','ConsumptionDate','Intervention','InterventionType',
              'kWh','LkWh','B1LkWh','B2LkWh','MaxC','MinC','Rain','SunE',
              'dMo','dTu','dWe','dTh','dFr','dSa','InterventionDate')]
write.csv(data1,'ANON.Experimental.Data.csv',row.names=F)

ca <- 'Sample of Modelling Data Set'
a <- rep('c',10)
kable(data.frame(NMI=data1[,1],Date=data1[,2],Group=data1[,4],round(data1[,c(3,5:6,9:12)],2))[c(1:4,5988:5991,35502:35505),],format='markdown',row.names=T,caption=ca,align=a)

ca <- 'Intervention Counts'
a <- rep('c',4)
z <- aggregate(InterventionType~NMI,data1,function(x) substr(paste0(x,collapse=''),1,1))
p <- t(table(z$InterventionType))
rownames(p) <- 'Count'
kable(cbind(p,Total=sum(p)),format='markdown',caption=ca,align=a)

# kWh across Time by Intervention
library(mgcv)
ds <- as.POSIXct(c('2012-07-01','2012-11-01','2013-03-01','2013-07-01','2013-11-01','2014-03-01'))
ggplot(data=data1,aes(x=ConsumptionDate,y=LkWh)) + apa_style + geom_line(col='gray', lwd=.01) + scale_x_datetime(expand=c(0,0),breaks=ds,limits=as.POSIXct(c('2012-06-25','2014-04-05')),labels=date_format('%b-%y')) + scale_y_continuous(expand=c(0,0),limits=c(2,5.07),breaks=c(2:5)) + facet_grid(InterventionType~.,scales='free_y') + geom_smooth(method='gam',formula=y~s(x,bs='cs'),col='black',lwd=1,lty=2) + geom_vline(aes(xintercept=as.numeric(InterventionDate)),lwd=.1)

# RegModels Formula
model <- formula(LkWh~Intervention+MaxC+MinC+Rain+SunE+dMo+dTu+dWe+dTh+dFr+dSa+B1LkWh+B2LkWh)

# Training
library(plyr)
y <- dlply(data1, 'NMI', function(DataFrame) {lm(model,data=DataFrame, na.action=na.omit)})

# Coefficients, Intervention CIs and R2
Coefficients <- function(fit) {exp(coef(fit))-1}
Coefficients <- ldply(y, Coefficients)
Coefficients['(Intercept)'] <- Coefficients['(Intercept)']+1
CIs <- function(fit) {exp(confint(fit,'Intervention'))-1}
CIs <- ldply(y, CIs)
Rsq <- function(fit) c(rsq=summary(fit)$adj.r.squared)
Rsq <- ldply(y, Rsq)
library(lmtest)
DWp <- function(fit) c(dwp=dwtest(fit,alt='two.sided')$p.value)
DWp <- ldply(y, DWp)
Results <- cbind(Coefficients$Intervention,CIs['2.5 %'],CIs['97.5 %'],Rsq$rsq,DWp$dwp,
                 aggregate(InterventionType~NMI,data1,function(x) substr(paste0(x,collapse=''),1,1))[,2])
names(Results) <- c('Intervention','IntLowerCI','IntUpperCI',
                    'ExplainedVar','DurbinWp','InterventionType')

# Model Quality and Back-Transformed Intervention Coefficients (95%CI)
l <- c('Percent Change in kWh (Intervention Coefficients)',expression(paste('Explained Variance ',italic((Adj-R^2)))))
ggplot(data=Results, aes(x=ExplainedVar,y=Intervention,shape=InterventionType)) + apa_style +
  geom_hline(aes(yintercept=0),lty=2) + facet_grid(~InterventionType) + labs(y=l[1],x=l[2]) +
  geom_vline(aes(xintercept=0)) + geom_vline(aes(xintercept=1),colour='grey') +
  geom_errorbar(aes(ymax=IntUpperCI,ymin=IntLowerCI),width=.05) + theme(legend.position='none') +
  geom_point(size=3,fill='grey') + scale_shape_manual(values=rep(21,4)) +
  scale_x_continuous(expand=c(0,0),limits=c(0,1),breaks=c(0,.25,.5,.75),labels=percent) +
  scale_y_continuous(expand=c(0,0),limits=c(-.31,.2),breaks=seq(-.2,.2,by=.1),labels=percent)

ca <- 'Intervention A & C Summary'
cb <- 'Intervention H & W Summary'
z <- rep('c',6)
library(psych)
a <- data.frame(row.names=names(Results)[-6])
for(i in colnames(p)){ b <- describeBy(Results[,-6],Results[,6])[[i]][,c(3,8:9)]
names(b) <- paste0(names(b),i); a <- cbind(a,b)}
kable(round(a[,1:6],2),format='markdown',caption=ca,align=z)
kable(round(a[,7:12],2),format='markdown',caption=cb,align=z)

ca <- 'Intervention A & C Summary'
cb <- 'Intervention H & W Summary'
z <- rep('c',6)
library(psych)
a <- data.frame(row.names=names(Results)[-6])
for(i in colnames(p)){ b <- describeBy(Results[,-6],Results[,6])[[i]][,c(3,8:9)]
names(b) <- paste0(names(b),i); a <- cbind(a,b)}
kable(round(a[,1:6],2),format='markdown',caption=ca,align=z)
kable(round(a[,7:12],2),format='markdown',caption=cb,align=z)

ca <- 'One-Way ANOVA of Step Variable Coefficients'
z <- rep('c',5)
SlopeANOVA <- aov(I(log(Intervention+1))~InterventionType,data=Results)
kable(data.frame(anova(SlopeANOVA)),format='markdown',caption=ca,align=z)

ca <- 'Tukey HSD Post-Hoc Tests'
z <- rep('c',5)
kable(data.frame(TukeyHSD(SlopeANOVA)$InterventionType),format='markdown',caption=ca,align=z)

# ANOVA Hypothesis Testing
ca <- 'Residual Variances by NMI'
a <- rep('c',7)
resvar <- ldply(y,function(fit) {data.frame(ResidualVar=round(var(resid(fit)),3))})
z <- cbind(resvar[1:28,],resvar[29:56,],resvar[57:84,])
kable(z,format='markdown',caption=ca,align=a)


# Full Results Table
ca <- 'Results Table (back-transformed)'
a <- rep('c',8)
Results[names(Results[,-6])] <- lapply(Results[,-6],round,3)
kable(Results,format='markdown',caption=ca,align=a,row.names=T)

# Remaining Coefficients Table
ca <- 'Remaining Coefficients Table'
a <- rep('c',11)
z <- round(Coefficients[,-c(1,3)],2)
names(z) <-  substr(c('B0',names(z)[-1]),1,3)
kable(z,format='markdown',caption=ca,align=a,row.names=T)

# Variance Inflation Factors
ca <- 'Variance Inflation Factors (Standardised)'
a <- rep('c',13)
library(car)
varif <- function(fit) c(varif= sqrt(vif(fit)));varif <- round(ldply(y, varif)[,-1],1)
names(varif) <- names(z)
kable(varif,format='markdown',caption=ca,align=a,row.names=T)

# Studentized Residual & Squared Mahalanobis Distance
hvs <- ldply(y,function(fit) {data.frame(hvs=hatvalues(fit))})
SMD <- (hvs$hvs-1/32)*(32-1)
smd <- sum(SMD>qchisq(.999,11))
fvs <- ldply(y,function(fit) {data.frame(StudentResid=c(0,0,rstudent(fit)),StandardPred=c(0,0,scale(predict(fit))))})
fvs <- cbind(fvs,data1)
ggplot(data=fvs, aes(x=StandardPred,y=StudentResid)) + geom_hline(aes(yintercept=c(-3,3)),lty=2) +
  geom_point(pch=1,colour='#999999') + apa_style + geom_density2d(colour='black') +
  facet_wrap(~InterventionType,scales='free')

# Supplemental Time-Series Control by NMI
ggplot(data=data1[data1$InterventionType=='C',],aes(x=ConsumptionDate,y=kWh)) + apa_style +
  geom_line(lwd=.01) + scale_x_datetime(expand=c(0,0),breaks=ds,labels=date_format('%b-%y')) +
  scale_y_continuous(expand=c(0,0)) + facet_grid(NMI~.,scales='free') +
  geom_vline(aes(xintercept=as.numeric(InterventionDate)),lwd=.7)

# Supplemental Time-Series Assessor by NMI
ggplot(data=data1[data1$InterventionType=='A',],aes(x=ConsumptionDate,y=kWh)) + apa_style +
  geom_line(lwd=.01) + scale_x_datetime(expand=c(0,0),breaks=ds,labels=date_format('%b-%y')) +
  scale_y_continuous(expand=c(0,0)) + facet_grid(NMI~.,scales='free') +
  geom_vline(aes(xintercept=as.numeric(InterventionDate)),lwd=.7)

# Supplemental Time-Series Workshop by NMI
ggplot(data=data1[data1$InterventionType=='W',],aes(x=ConsumptionDate,y=kWh)) + apa_style +
  geom_line(lwd=.01) + scale_x_datetime(expand=c(0,0),breaks=ds,labels=date_format('%b-%y')) +
  scale_y_continuous(expand=c(0,0)) + facet_grid(NMI~.,scales='free') +
  geom_vline(aes(xintercept=as.numeric(InterventionDate)),lwd=.7)

# Supplemental Time-Series Home Visits by NMI
y <- data1[data1$InterventionType=='H',]
y1 <- y[y$NMI %in% unique(y$NMI)[1:20],]
y2 <- y[y$NMI %in% unique(y$NMI)[21:40],]
y3 <- y[y$NMI %in% unique(y$NMI)[41:60],]

ggplot(data=y1,aes(x=ConsumptionDate,y=kWh)) + apa_style + geom_line(lwd=.01) +
  scale_x_datetime(expand=c(0,0),breaks='6 months',labels=date_format('%m')) +
  scale_y_continuous(expand=c(0,0)) + facet_wrap(~NMI,scales='free') +
  geom_vline(aes(xintercept=as.numeric(InterventionDate)),lwd=.1,lty=5)

ggplot(data=y2,aes(x=ConsumptionDate,y=kWh)) + apa_style + geom_line(lwd=.01) +
  scale_x_datetime(expand=c(0,0),breaks='6 months',labels=date_format('%m')) +
  scale_y_continuous(expand=c(0,0)) + facet_wrap(~NMI,scales='free') +
  geom_vline(aes(xintercept=as.numeric(InterventionDate)),lwd=.1,lty=5)

ggplot(data=y3,aes(x=ConsumptionDate,y=kWh)) + apa_style + geom_line(lwd=.01) +
  scale_x_datetime(expand=c(0,0),breaks='6 months',labels=date_format('%m')) +
  scale_y_continuous(expand=c(0,0)) + facet_wrap(~NMI,scales='free') +
  geom_vline(aes(xintercept=as.numeric(InterventionDate)),lwd=.1,lty=5)

# Bibliography
library(bibtex)
write.bib(rev(.packages())[-c(2:7)],verbose=F)
read.bib('Rpackages.bib')
