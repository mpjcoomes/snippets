#Correspondence Analysis 2:

#Chi-Square Tests of Independence
x <- data.frame(cbind(Male=c(60, 40), Female=c(20,80)), row.names=c("Attend college", "Do not attend college"))
x[,"Total"] <- rowSums(x)
x["Total",] <- colSums(x)
x <- data.frame(cbind(Male=c(60, 40), Female=c(20,80)), row.names=c("Attend college", "Do not attend college"))
x[,"Proportion"] <- rowSums(x)/sum(colSums(x)); x


#Raw data table:
library(foreign)
library(gmodels)
dir <- "College.sav"
x <- read.spss(dir, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
	trim.factor.names = FALSE, trim_values = TRUE, use.missings = TRUE)
y <- x[rep(1:nrow(x),x$freq),]
CrossTable(y[,1], y[,2], digits=2, expected=T, prop.c=T, chisq=T, prop.chisq=T, prop.r=F, prop.t=F,
	 asresid=T, format="SPSS")

z <- data.frame(y[y[,1]=="Attended college",2])
ggplot(z, aes(factor(z[,1]))) +  ylim(0,100) +
	geom_bar(width=.5,fill="white", colour="black") +
	labs(title='Gender and College', x = "Group", y = "Proportion Attending College") +
	theme_bw(base_size=12,base_family='serif') +
	theme(legend.justification=c(0,0),legend.position=c(0,0),legend.title=element_blank(),
	panel.grid.major =  element_line(colour = NA)) +
	stat_bin(aes(label=paste(..count..,"%",sep='')), vjust=-1,geom="text")

library(vcd)
assoc(y[,1:2], shade=TRUE)

#Correspondence Analysis
library(ca)
library(foreign)
dir <- "smoking groups.sav"
x <- read.spss(dir, use.value.labels = TRUE, to.data.frame = TRUE, max.value.labels = Inf,
	trim.factor.names = FALSE, trim_values = TRUE, use.missings = TRUE)
y <- x[rep(1:nrow(x),x$freq),]
mytable <- with(y, table(y[,1],y[,2])) # create a 2 way table
mytable_Totals <- data.frame(rbind(mytable, colSums(mytable)))
mytable_Totals[,'Total'] <- rowSums(mytable_Totals)
rownames(mytable_Totals)[length(rownames(mytable_Totals))] <- 'Total'
mytable_Totals
CrossTable(y[,1], y[,2], digits=2, expected=T, prop.c=T, chisq=T, prop.chisq=T, prop.r=F, prop.t=F,
	 asresid=T, format="SPSS")

library(vcd)
#mosaic(mytable, split_vertical=T, shade=TRUE, legend=TRUE)
assoc(y[,1:2], shade=TRUE)

prop.table(mytable, 1) # row percentages
prop.table(mytable, 2) # column percentages
chisq.test(mytable)
chisq.test(mytable)$statistic/sum(x$freq) #Inertia = chi^2/n
fit <- ca(mytable)
fit$sv
summary(fit)


#Contribution Tables
attach(summary(fit))
RowPnts <- data.frame(StaffGroup=rows$name, Mass=rows$mass, Inertia=round(fit$rowinertia,3),
		PointIn=rows[c(7,10)], DimIn=rows[c(6,9)], Total=rowSums(rows[c(6,9)]))
detach(summary(fit))
RowPnts
a <- RowPnts[order(-RowPnts$PointIn.ctr),"PointIn.ctr"]/10
paste('Senior & Junior Employees account for ', sum(a[1:2]), '% of Dim1 inertia.', sep='')
a <- RowPnts[order(-RowPnts$PointIn.ctr.1),"PointIn.ctr.1"]/10
paste('Senior & Junior Managers account for ', sum(a[1:2]), '% of Dim2 inertia.', sep='')

attach(summary(fit))
ColPnts <- data.frame(StaffGroup=columns$name, Mass=columns$mass, Inertia=round(fit$colinertia,3),
		PointIn=columns[c(7,10)], DimIn=columns[c(6,9)], Total=rowSums(columns[c(6,9)]))
detach(summary(fit))
ColPnts
a <- ColPnts[order(-ColPnts$PointIn.ctr),"PointIn.ctr"]/10
paste('Non-smokers account for ', a[1], '% of Dim1 inertia.', sep='')
a <- ColPnts[order(-ColPnts$PointIn.ctr.1),"PointIn.ctr.1"]/10
paste('Light and Heavy smokers account ', sum(a[1:2]), '% of Dim2 inertia.', sep='')


#Standard Package Plots:
dev.new(width=10, height=5)
par(family="serif",mfrow=c(1,2),adj=0,cex=.8,mar=c(4,2,.5,.5))
plot(fit)
plot(fit,mass=T,contrib="absolute",map="symbiplot",arrows=c(T,T),xlim=c(-.7, .7),ylim=c(-1, 1))


#Basic Plotting: Linear Regression
a <- data.frame(fit$rowcoord[,1:2])
a["Group"] <- "Employee Level"
a["ID"] <- fit$rownames
b <- data.frame(fit$colcoord[,1:2])
b["Group"] <- "Smoking Level"
b["ID"] <- fit$colnames
caplot <- rbind(a,b)

dev.new(width=5, height=5)
par(family="serif",cex=.8)
plot(fit$colcoord, col="red", xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5),
	xlab='Dimension 1', ylab='Dimension 2', main='Symmetrical Normalisation')
abline(lm(fit$colcoord[,2] ~ fit$colcoord[,1]), col="red")
points(fit$rowcoord, pch=3)
abline(lm(fit$rowcoord[,2] ~ fit$rowcoord[,1]), lty="longdash")
text(fit$colcoord[,1], fit$colcoord[,2], labels=fit$colnames, pos=4)
text(fit$rowcoord[,1], fit$rowcoord[,2], labels=fit$rownames, pos=4)

lm.mod <- lm(fit$rowcoord[,2] ~ fit$rowcoord[,1])
a <- coef(lm.mod)[1] 
b <- coef(lm.mod)[2]
x <- caplot[caplot[,3]=="Employee Level",1]
y <- caplot[caplot[,3]=="Employee Level",2]
x1 <- (x+b*y-a*b)/(1+b^2)
y1 <- a + b*x1
ss <- data.frame(x0=x, y0=y, x1=x1, y1=y1)
segments(x0=ss$x0, x1=ss$x1, y0=ss$y0, y1=ss$y1, lty=2)

lm.mod <- lm(fit$colcoord[,2] ~ fit$colcoord[,1])
a <- coef(lm.mod)[1] 
b <- coef(lm.mod)[2]
x <- caplot[caplot[,3]=="Smoking Level",1]
y <- caplot[caplot[,3]=="Smoking Level",2]
x1 <- (x+b*y-a*b)/(1+b^2)
y1 <- a + b*x1
ss <- data.frame(x0=x, y0=y, x1=x1, y1=y1)
segments(x0=ss$x0, x1=ss$x1, y0=ss$y0, y1=ss$y1, col='red')


#Advanced Plotting: Linear Regression
a <- data.frame(fit$rowcoord[,1:2])
a["Group"] <- "Employee Level"
a["ID"] <- fit$rownames
b <- data.frame(fit$colcoord[,1:2])
b["Group"] <- "Smoking Level"
b["ID"] <- fit$colnames
caplot <- rbind(a,b)

lm.mod <- lm(fit$rowcoord[,2] ~ fit$rowcoord[,1])
a <- coef(lm.mod)[1] 
b <- coef(lm.mod)[2]
x <- caplot[caplot[,3]=="Employee Level",1]
y <- caplot[caplot[,3]=="Employee Level",2]
x1 <- (x+b*y-a*b)/(1+b^2)
y1 <- a + b*x1
ss <- data.frame(x1=x1, y1=y1)

lm.mod <- lm(fit$colcoord[,2] ~ fit$colcoord[,1])
a <- coef(lm.mod)[1] 
b <- coef(lm.mod)[2]
x <- caplot[caplot[,3]=="Smoking Level",1]
y <- caplot[caplot[,3]=="Smoking Level",2]
x1 <- (x+b*y-a*b)/(1+b^2)
y1 <- a + b*x1
dd <- data.frame(x1=x1, y1=y1)

caplot ['X1_lm'] <- rbind(ss,dd)[,'x1']
caplot ['X2_lm'] <- rbind(ss,dd)[,'y1']

library(ggplot2)
p <- ggplot(caplot, aes(x=X1, y=X2, colour=factor(Group), label=ID))
p + geom_point(size=4,aes(shape=factor(Group))) + xlim(-3,3) + ylim(-3,3) +
	geom_text(size=3,hjust=1.2,vjust=0,family='serif') +
	labs(title='Row/Column Points Symmetrical Normalisation',x='Dimension 1',y='Dimension 2') +
 	theme_bw(base_size=12,base_family='serif') +
	theme(legend.justification=c(0,0),legend.position=c(0,0),legend.title=element_blank(),
	panel.grid.major =  element_line(colour = NA)) +
	geom_hline(yintercept = 0, linetype=3) + geom_vline(xintercept = 0, linetype=3) +
	stat_smooth(se=F,method='lm',fullrange = T) +
	geom_segment(aes(x=X1, y=X2, xend=X1_lm, yend=X2_lm)) +
	scale_colour_manual(values=c("#0072B2","#D55E00"))	#+ coord_flip()


#Advanced Plotting: Total Least Squares
a <- data.frame(fit$rowcoord[,1:2])
a["Group"] <- "Employee Level"
a["ID"] <- fit$rownames
b <- data.frame(fit$colcoord[,1:2])
b["Group"] <- "Smoking Level"
b["ID"] <- fit$colnames
caplot <- rbind(a,b)

library(MethComp)
TLS.modE <- Deming(fit$rowcoord[,1],fit$rowcoord[,2],boot=T)
a <- TLS.modE[1] 
b <- TLS.modE[2]
x <- caplot[caplot[,3]=="Employee Level",1]
y <- caplot[caplot[,3]=="Employee Level",2]
x1 <- (x+b*y-a*b)/(1+b^2)
y1 <- a + b*x1
ss <- data.frame(x1=x1, y1=y1)

TLS.modS <- Deming(fit$colcoord[,1],fit$colcoord[,2],boot=T)
a <- TLS.modS[1] 
b <- TLS.modS[2]
x <- caplot[caplot[,3]=="Smoking Level",1]
y <- caplot[caplot[,3]=="Smoking Level",2]
x1 <- (x+b*y-a*b)/(1+b^2)
y1 <- a + b*x1
dd <- data.frame(x1=x1, y1=y1)

caplot ['X1_Dem'] <- rbind(ss,dd)[,'x1']
caplot ['X2_Dem'] <- rbind(ss,dd)[,'y1']
df <- data.frame(a=c(TLS.modE[1],TLS.modS[1]), b=c(TLS.modE[2],TLS.modS[2]))

library(ggplot2)
p <- ggplot(caplot, aes(x=X1, y=X2, colour=factor(Group), label=ID))
p + geom_point(size=4,aes(shape=factor(Group))) + xlim(-3,3) + ylim(-3,3) +
	geom_text(size=3,hjust=1.2,vjust=0,family='serif') +
	labs(title='Row/Column Points Symmetrical Normalisation (TLS Regression)',
	x='Dimension 1',y='Dimension 2') +
 	theme_bw(base_size=12,base_family='serif') +
	theme(legend.justification=c(0,0),legend.position=c(0,0),legend.title=element_blank(),
	panel.grid.major =  element_line(colour = NA)) +
	geom_hline(yintercept = 0, linetype=3) + geom_vline(xintercept = 0, linetype=3) +
	#stat_smooth(se=F,method='lm',fullrange = T) +
	geom_abline(aes(intercept=a, slope=b),colour='#0072B2', data=df[1,]) +
	geom_abline(aes(intercept=a, slope=b),colour='#D55E00', data=df[2,]) +
	geom_segment(aes(x=X1, y=X2, xend=X1_Dem, yend=X2_Dem)) +
	scale_colour_manual(values=c("#0072B2","#D55E00"))

#Distances:
z <- sqrt((caplot$X1_Dem - caplot$X1)^2 + (caplot$X2_Dem - caplot$X2)^2)
data.frame(caplot["ID"], z)[order(data.frame(caplot["ID"],z)[,2]),]





