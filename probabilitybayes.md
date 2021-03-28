Probability
=

 Sample space is the collection of possible outcomes.

    o = {1,2,3,4,5,6}, 6-sided die roll
    Ø or {} = empty set null event
    e & ~e element & not element of ...

Event is a subset of sample space, w of O.

    E = {2,4,6}, die roll is even

Elementary or simple event is one specific result.

    w = {4}, die roll is four

Set Operations

    w e E , E occurs when w occurs, w is an element of E
    w ~e E , E does not occurs when w occurs, is not an element
    E s F , occurence of E implies occurence of F
    E n F , intersection, both E and F occur
    E n F = Ø , E and F are mutually exclusive
    E U F , union, at least one of E and F occur
    ¬E ~E ?E E` Ec E_bar, event E does not occur

 All events must have a probability greater than or equal to zero.

    P(Ei)=0

For an event E c O , 0 LE P(E) LE 1 & P(o) = 1.

    P(E1UE2) = P(E1) + P(E2) , while E1nE2 = Ø
    P(Ø) = 0
    P(E) = 1-P(~E)

 If A c B then P(A) LE P(B)
 
    P(An~B) = P(A)-P(AnB)

 Common Probability Relationships
 
	P(AnB) = P(A|B)P(B)
	P(AUB) = P(A)+P(B)-P(AnB) = 1-P(~An~B)
	P(A|B) = P(AnB)/P(B)
	P(~DnT) = P(T|~D)P(~D)
	P(~A|B) = P(~AnB)/P(B)
	P(A|~B) = P(An~B)/P(~B) p74Daniel
	P(An~B) = P(A)-P(AnB) s5w1
	P(A|~B) = [P(A)-P(AnB)] / P(~B)

Two events are mutually exclusive if they cannot occur simultaneously.

Probabilities calculated with a marginal total at the numerator and the total group as the denominator result in a marginal probability Given some variable that can be broken down into m categories designated by A1, A2, ..., Ai, ..., Am, and another jointly occurring variable that is broken down into n categories, B1, B2, ..., Bj, ..., Bn, the marginal probability of Ai, mean P(Ai), equals the sum of the joint probabilities of Ai with all the categories of B. P(Ai) = Sigma(P(AinBj)), for all values of j.

Complementary events: probability of an event A is equal to 1 minus the probability of its complement.

    P(~A) = 1-P(A)

Exhaustiveness in a probabilistic process allows for all possible events taken together their total probability is 1 sum of probabilities of mutually exclusive outcomes equals 1

    P(E1)+P(E2)+...+P(En)=1

 The probability that a subject picked at random possesses two characteristics at the same time is the intersection or joint probability .

    P(AnB) = joint occurance of A and B

While EinEj = Ø any two mutually exclusive events, the probability of union is their sum note. EinEj = Ø is mutually exclusive and EinEj = 0 is independent.

	P(EiUEj) = P(Ei)+P(Ej) , 
	P(AnB) = P(B)P(A), if P(A)!=0, P(B)!=0,
	P(A|B) = P(A)
	P(B|A) = P(B)

When EinEj!=Ø, addition rule.

    P(AUB) = P(A)+P(B)-P(AnB) = 1-P(~An~B)

Probabilities calculated with a subset of the total group as the denominator result in a conditional probability.

    P(A|B) easier with Bayes

Multiplication rule: a joint probability may be computed as the product of a marginal probability and a conditional probability.

	P(AnB) = P(B)P(A|B) = P(A)P(B|A), if P(A):P(B)!=0
	P(A|B) = P(AnB)/P(B), if P(B)!=0
	P(B|A) = P(AnB)/P(A), if P(A)!=0

Bayesian
=
Bayes' Theorem assuming binary outcomes.

            P(B|A)P(A)
    P(A|B) =  ----------
                P(B)

Or alternatively.

                    P(B|A)P(A)
    P(A|B) =  --------------------------
                P(B|A)P(A)+P(B|~A)P(~A)

Bayes' Theorem Applied: Disease Testing

    Prevalence = P(D)

Sensitivity (i.e. T Pos Rate, 'hits') = P(T|D).

	proportion of real positives correctly classified as positive
	probability of testing positive, given that has disease
	P(+|D)=P(T|D)=a/(a+c), o = {all real positives}
	P(T|D)= P(DnT) / P(D)

Specificity (i.e. T Neg Rate) P(~T|~D).

	proportion of real negatives correctly classified as negative
	probability of testing negative, given that does not have disease
	P(-|~D)=P(~T|~D)=d/(b+d), o = {all real negatives}
	P(~T|~D) = P(~Tn~D)/P(~D)

False Negative rate = 1-Sensitivity, a test indicates a negative status when status is positive (e.g. 1-.8 = 20% of the time a person they have the disease they test negative).

    P(~T|D) = P(~TnD)/P(D)

False Positive rate = 1-Specificity, a test indicates a positive status when status is negative (e.g. 1-.74 = 26% of the time a person tests positive they don't have the disease).

    P(T|~D) = P(Tn~D)/P(~D)

Predictive Value Positive, probability of diseased D given + test.

    P(D|+)=P(D|T)= P(T|D)P(D) / P(T|D)P(D)+P(T|~D)P(~D)

Predictive Value Negative, probability of not diseased ~D given - test.

    P(~D|-)=P(~D|~T)=  P(~T|~D)P(~D) / P(~T|~D)P(~D)+P(~T|D)P(D)

Disease and Test Relationships.

	P(D|T) = sensitivity * prevalence / P(T)
	P(T) = P(DnT) + P(~DnT)
	P(D) = P(DnT) + P(Dn~T)
	P(~DnT) = P(T) - P(DnT)
	P(~DnT) = P(T|~D)P(~D) = FPR* 1-P(D) = 1-Specificity * 1-Prevalence

	P(DnT) = P(T|D)P(D)
	P(A|B) = P(AnB)/P(B)
	P(D|T) = P(DnT)/P(T)
	P(A|B) = P(B|A)P(A) / P(B|A)P(A) +P(B|~A)P(~A)
	P(D|T) = P(T|D)P(D) / P(T|D)P(D) +P(T|~D)P(~D)
	P(D|T) = P(DnT) / P(DnT) + P(Tn~D)
	P(D|+) = P(+|D)P(D) / P(+|D)P(D)+P(+|~D)P(~D)
	P(~D|+)= P(+|~D)P(~D)/P(+|D)P(D)+P(+|~D)P(~D)
	P(~D|~T) = P(~T|~D)P(~D) / P(~T|~D)P(~D)+P(~T|D)P(D)
	P(T|~D) = FPR = 1-Specificity = = P(Tn~D)/P(~D)

Therefore.

    P(D|T) = P(DnT)/P(T)

Numerator comes from P(DnT) = P(T|D)P(D)
Denominator () comes from P(T): How?

	P(T) = P(+ test) with or without the disease
	if P(T) = P(DnT) + P(~DnT) or P(TnD)+P(Tn~D) with mutually exclusive addition rule
	P(DnT) = P(T|D)P(D) by multiplication rule
	P(~DnT) = P(T|~D)P(~D) rewrites
	P(DnT) + P(~DnT) =  P(T|D)P(D) + P(T|~D)P(~D)
	P(DnT) = P(T|D)P(D) = Sensitivity*Prevalence
	P(T) = P(T|D)P(D) + P(T|~D)P(~D)

### Example
```r
Sens <- .84        #hits
Spec <- .77        #1-false hits, 1-P(T|~D)
Prev <- .2        #disease rate P(D)
Infreq <- 1-Prev    #P(~D)
FHits <- 1-Spec    #1-P(T|D)
FNR <- 1-Sens    # P(~T|D)

# P(D|T) =    P(DnT)      /   P(DnT)  +  P(Tn~D)
#           P(T|D)P(D)    / P(T|D)P(D) + P(T|~D)P(~D)

Sens*Prev/(Sens*Prev+FHits*Infreq)  #0.4772727

# P(~D|T) = P(~Tn~D)    /  P(~Tn~D) +  P(DnT)
#           P(~T|~D)P(~D)/( P(~T|~D)P(~D) +  P(~T|D)P(D) )

Spec*Infreq/(Spec*Infreq+FNR*Prev)  #0.9506173
```

Bayesian Updating
=
Iterated Dice problem: five die with 4,6,8,12 & 20 faces.

	If you roll a 6, what is the likelihood of each die?
	pA = p(of each die type)
	pB = p(of getting a '6' through any pathway)
	pBgA = p(of getting a '6', given one die type)
	pAnB = p(of getting '6' with each die type)
```r
pA <- rep(1/5,5) # Priors
pBgA <- c(0/4,1/6,1/8,1/12,1/20)

# P(AnB) = P(B|A)P(A)

pAnB <- pBgA*pA

# P(A|B) = P(AnB)/P(B)

Posterior <- pAnB/sum(pAnB) #0.000 0.392 0.294 0.196 0.118
```
Iterated, roll an 8.
```r
pA <- Posterior

# pB2 = p(of getting an '8' through any pathway)

pBgA <- c(0/4,0/6,1/8,1/12,1/20)
pAnB <- pBgA*Posterior
Posterior <- pAnB/sum(pAnB)
```
Roll 7, 7, 5, 4.
```r
pBgA <- c(0/4,0/6,1/8,1/12,1/20)
pAnB <- pBgA*Posterior
Posterior  <-  pAnB/sum(pAnB)
pBgA <- c(0/4,0/6,1/8,1/12,1/20)
pAnB <- pBgA*Posterior
Posterior  <-  pAnB/sum(pAnB)
pBgA <- c(0/4,1/6,1/8,1/12,1/20)
pAnB <- pBgA*Posterior
Posterior  <-  pAnB/sum(pAnB)
pBgA <- c(1/4,1/6,1/8,1/12,1/20)
pAnB <- pBgA*Posterior
Posterior <- pAnB/sum(pAnB)
```

Trainspotting Problem
=
Freight carrier operates 100 to 1000 locomotives. You see locomotive '#321'. How many locomotives do they have?
```r
pA <- rep(1/900,900) # Priors
plot(100:999,pA,type='l')
pBgA <- c(rep(0,220),1/321:1000)
pAnB <- pBgA*pA
Posterior <- data.frame(101:1000,pAnB/sum(pAnB))
plot(Posterior,type='l')
Posterior[order(Posterior[,2]),] #MLE
pMean <- sum((Posterior[,1]*Posterior[,2])) #Posterior Mean
for (i in 1:900) {
    Posterior[i,'Total']<-sum(Posterior[1:i,2])
    }
subset(Posterior, Total<.05) #Posterior lower CI
subset(Posterior, Total>.95) #Posterior upper CI
```

Coin Toss Problem
=
pA = coin probability of one side (e.g. 'H')
the hypothetical probability of 'H'
pB = coin probability of one side given 0 to 100 flips
```r
pA <- rep(1/101,101)
pBgA <- 0:100/100
for (i in 1:140) {
    pAnB <- pBgA*pA
    pA <- pAnB/sum(pAnB)
    }
for (i in 1:110) {
    pAnB <- (1-pBgA)*pA
    pA <- pAnB/sum(pAnB)
    }
Posterior <- data.frame(0:100,pA)
plot(Posterior,type='l')
round(Posterior,3)
```
.021 prob of 50, or fair, but this isn't answering the question we want because it's dependent on the 100-unit granularity of the cutup.
```r
tail(Posterior[order(Posterior[,2]),]) #MLE
sum((Posterior[,1]*Posterior[,2])) #Posterior Mean
for (i in 1:100){
    Posterior[i,'Total']<-sum(Posterior[1:i,2])
    }
tail(subset(Posterior, Total<.05),1) #Posterior lower CI
head(subset(Posterior, Total>.95),1) #Posterior upper CI
```

Alternate calculation method.
```r
pA <- rep(1/101,101)
h <- 1
pBgA <- c(rep(0,h),h:100/100)
pAnB <- pBgA*pA
Posterior1 <- data.frame(0:100,pAnB/sum(pAnB))
h <- 2
pBgA <- c(rep(0,h),h:100/100)
pAnB <- pBgA*Posterior1[,2]
Posterior2 <- data.frame(0:100,pAnB/sum(pAnB))
h <- 3
pBgA <- c(rep(0,h),h:100/100)
pAnB <- pBgA*Posterior2[,2]
Posterior3 <- data.frame(0:100,pAnB/sum(pAnB))

#Then flipped a 'T'
pBgA <-  c(100:h/100,rep(0,h))
pAnB <- pBgA*Posterior3[,2]
Posterior4 <- data.frame(0:100,pAnB/sum(pAnB))

par(mfrow=c(2,2))
plot(Posterior1,type='l',main='H, Hyp:1% chance of H')
plot(Posterior2,type='l',main='HH, Hyp:2% chance of H')
plot(Posterior3,type='l',main='HHH, Hyp:3% chance of H')
plot(Posterior4,type='l',main='HHHT, Hyp:3% chance of H, 1% chance tail')
```

Maximum Likelihood Estimation
=
After 1 'H', the most likely conclusion is that it's a double-sided coin.
```r
tail(Posterior1[order(Posterior1[,2]),])
sum((Posterior1[,1]*Posterior1[,2])) #Posterior Mean
tail(Posterior4[order(Posterior4[,2]),])
sum((Posterior4[,1]*Posterior4[,2])) #Posterior Mean
```

Bayes Factor
=
Likelihood Ratio

    p(D|H)/p(D|~H)>1

That the prob of the data under the hypothesis is larger than the prob of the data under the not hypothesis.
Coin:

    Easy: p(D|F) that the coin is fair (e.g. 50%)
    Hard: p(D|B) that the coin is bias, hard because 'bias' is unspecified which uses the same data twice.

```r
heads <- 140
tails <- 110
x <- 0:100/100
like <- x^heads*(1-x)^tails

#bias set
pA <- c(rep(1,50),0,rep(1,50)) #two sided uniform
pA <- c(0:49,0,49:0) #alternate pyramid prior
pA <- c(rep(0,51),49:0) #one sided hypothesis
pA <- pA*1/sum(pA) #normalise
H1 <- sum(like*pA)

#fair set
pA <- c(rep(0,50),1,rep(0,50))
H0 <- sum(like*pA)  

H1/H0 #Very important how you define the priors
```

Note how binomial distribution method equivalent to the previous loop method when priors are as they were (i.e. uniform).
```r
pA <- rep(1/101,101)
pAnB <- like*pA #Bin prob of each outcome * priors
PosteriorBin <- data.frame(0:100,pAnB/sum(pAnB))
data.frame(Posterior,PosteriorBin)
plot(Posterior,type='l')
lines(PosteriorBin,lwd=2)
```

SAT Problem
=
```r
priors <- read.csv('~/sat_ranks.csv',header=F,skip=3)[-(62:64),1:2]
pA <- rev(priors[,2]*1/sum(priors[,2])) #normalise and reverse to sort ascending
x <- seq(0,1,length=61) #prob of getting a question right

raws <- read.csv('~/sat_scale.csv')
#raws[5:17,3:4]
#data.frame(seq(730,800,10),seq(50,54,length=8))
```
Alice 780
```r
yes <- 52.85714
no <- 54-yes
like <- x^yes*rev(x)^no
pAnB <- like*pA
Alice <- data.frame(x,A=pAnB/sum(pAnB))
```
Bob 760
```r
yes <- 51.71429
no <- 54-yes
like <- x^yes*rev(x)^no
pAnB <- like*pA
Bob <- data.frame(x,pAnB/sum(pAnB))
```
Posteriors: Probability of getting any random question correct (x-axis).
```r
plot(Alice,type='l')
lines(Bob,col=2)
total <- data.frame(Alice,B=Bob[,2])
```
Assume independence between Alice and Bob's 'x' (i.e. prob correct).
    Every outcome can be calculated eg. ~3600 possible 'x' combinations
    Sum prob of all combinations where Alice's 'x' is higher
    Sum prob of all combinations where Bob's 'x' is higher
    Result: naturally normalised probabilities of 'smarter'
```r
x <- matrix(0)
for (i in total$A) {
    for (b in total$B) {
        if (i>b) {
            x <- rbind(x,i*b)
            }
        }
    }
y <- matrix(0)
for (i in total$B) {
    for (b in total$A) {
        if (i>b) {
            y <- rbind(y,i*b)
            }
        }
    }
sum(x) #Prob Alice is 'smarter'
sum(y) #Prob Bob is 'smarter'
```
Cumulative Distribution Function
=
Area under the curve less than or equal 'domain number' x.
```r
x <- seq(-4,4,length=1000)
n <- 3
plot(x, pnorm(x, mean=0, sd=1),type='l')
lines(x,pnorm(x, mean=0, sd=1, #area equal or above x
              lower.tail=F),col='blue')
lines(x,dnorm(x, mean=0, sd=1),col='green')
```

Bayesian
```r
x <- seq(0,1,length=1000) #Probability that x or fewer % occurs
a <- 2 # Heads?
B <- 2 # Tails?
plot(x, pbeta(x,1,.1),type='l')
lines(x,pbeta(x,.1,1),col='blue')
lines(x,pbeta(x,10,10),col='red')
lines(x,pbeta(x,100,100),col='green')
lines(x,pbeta(x,1,100),col='purple')
lines(x,pbeta(x,2,1),col='grey')
```
Find Beta(a,ß) shape parameters matching knowledge of two quantiles p=value of first and second probability, x=value of first and second qauntile.
```r
library(LearnBayes)
quantile1 <- list(p=.7, x=0.75) 
quantile2 <- list(p=.8, x=0.85)
params <- beta.select(quantile1,quantile2)
alpha = params[1]
beta = params[2]
dbeta(.5, alpha, beta)
pbeta(.5, alpha, beta)

x <- seq(-4,4,length=50)
cex <- .7
plot(x,pt(x, df=1),type='l')
points(x,pt(x, df=5),pch='5',cex=cex)
points(x,pt(x, df=100),pch='?',cex=cex)
lines(x,pnorm(x, mean=0, sd=1),col='blue')

x <- c(seq(0,1,length=20),seq(1,4,length=10))
cex <- .7
plot(x,pf(x, df1=100, df2=100),type='l')
points(x,pf(x, df1=1, df2=1),pch='?',cex=cex)
points(x,pf(x, df1=5, df2=1),pch='?',cex=cex)
points(x,pf(x, df1=100, df2=1),pch='?',cex=cex)
points(x,pf(x, df1=100, df2=4),pch='4',cex=cex)

x <- c(seq(0,1,length=50),seq(1,8,length=100))
cex <- .7
n <- 1:9
plot(x,pchisq(x, df=n),pch=as.character(n),cex=cex)
abline(h=c(0,1))

pnorm(1.96);pnorm(1.96,lower.tail=F)
pnorm(-1.96);pnorm(-1.96,lower.tail=F)
pnorm(1.64);pnorm(1.64,lower.tail=F)
pnorm(-1.64);pnorm(-1.64,lower.tail=F)
pt(1.96,df=n-1);pt(1.96,df=n-1,lower.tail=F)
pf(1.96, df1=n-1, df2=n-1)

x <- seq(-1,2,length=1000)
plot(x, punif(x, min=0, max=1),type='l')
lines(x,punif(x, min=0, max=1, #area equal or above x
              lower.tail=F),col='blue')
lines(x,dunif(x, min=0, max=1),col='red')
punif(.75,0,1);punif(.75,0,1, lower.tail=F)
quantile(runif(x, min=0, max=1),.75)

flips <- 10
kHeads <- 0:flips
pHeads <- .5
pbinom(5, size=flips, prob=pHeads)
plot(kHeads ,pbinom(kHeads, size=flips, prob=pHeads),type='s')
lines(kHeads ,pbinom(kHeads, size=flips, prob=pHeads,
                     lower.tail=F),col='blue',type='s')
lines(kHeads ,dbinom(kHeads, size=flips, prob=pHeads),col='red',type='s')

Hits <- 0:20 #p of 0,1,2... when ? = X or when ExpectedHits = X
ExpectedHits <- 10 #positive real number ? = equal expected value & var of X
plot(Hits ,ppois(Hits , lambda=ExpectedHits),type='s')
points(Hits , ppois(Hits , lambda=ExpectedHits,
                    lower.tail=F),col='red',type='s')
points(Hits , dpois(Hits , lambda=ExpectedHits),col='green',type='s')
```

Multivariate Normal Anomaly Detection
=
```r
library(mnormt)
dmnorm(0, mean = rep(0, 2),  matrix(c(1,0,0,1),2,2)) #p(x1 on two variables, no corr)
dnorm(0, mean = 0, sd=sqrt(1))^2

dmnorm(0, mean = rep(0, 3),  matrix(c(1,0,0,0,1,0,0,0,1),3,3)) #p(x1 on three variables, no corr)
dnorm(0, mean = 0, sd=sqrt(1))^3

dmnorm(0, mean = rep(0, 2),  matrix(c(1,.5,.5,1),2,2)) #p(x1 on two corr variables)

data <- mvrnorm(n=1000, rep(0, 2),  matrix(c(1,.5,.5,1),2,2))
dmnorm(data[23,], mean = apply(data,2,mean),  cov(data))
```

Quantile Function
=
Inverse CDF, percent point function, value at which the p(random variable) will be less than or equal to q.
```r
x <- seq(0,1,length=1000)
n <- 3
plot(x, qnorm(x, mean=0, sd=1),type='l')
lines(x, qt(x, df=n-1),col='green')
abline(h=c(-1.96,1.96),v=c(0.025,.975),col='blue')
abline(h=c(-1.64,1.64),v=c(0.05,.95),col='red')

qnorm(.95);qnorm(.95,lower.tail=F)
qnorm(.05);qnorm(.05,lower.tail=F)
qnorm(.975);qnorm(.975,lower.tail=F)
qnorm(.025);qnorm(.025,lower.tail=F)

plot(x,qnorm(x),type='l',xlab='p Y-val or less occurred by chance')
plot(rev(x),qnorm(x),type='l',xlab='p Y-val or greater occurred by chance')

samp1 <- seq(0,20,length=1000)
plot(qnorm(x,mean=10,sd=2),x,type='l',ylab='p x-val or less occurred by chance')
lines(samp1,dnorm(samp1,mean=10,sd=2))

plot(qnorm(x,mean=10,sd=2,lower.tail=F),x,type='l'
     ,ylab='p x-val or greater occurred by chance')
lines(samp1,dnorm(samp1,mean=10,sd=2))

n <- 3 #t-critical values
qt(.95,df=n-1);qt(.95,df=n-1,lower.tail=F)
qt(.05,df=n-1);qt(.05,df=n-1,lower.tail=F)
qt(.975,df=n-1);qt(.975,df=n-1,lower.tail=F)
qt(.025,df=n-1);qt(.025,df=n-1,lower.tail=F)

#non-invertable cdf so quantile=smallest x such that F(x) >= p 
quants <- seq(0,1,.1) 
flips <- 10 
pHeads <- .5
qbinom(quants, size=flips, prob=pHeads)
```

Density
=
```r
x <- seq(-6,6,length=1000)
n <- 3
plot(x,dunif(x),type='l',col='red')
lines(x,dnorm(x, mean=0, sd=1))
lines(x,dt(x, df=n-1),col='green')
lines(x,df(x, df1=n-1, df2=n-2),col='purple')
flips <- 5
kHeads <- 0:flips
pHeads <- .5
lines(kHeads , pbinom(kHeads , size=flips, prob=pHeads),type='s',col='blue')
Hits <- 0:6 #p of 0,1,2... when ? = X or when ExpectedHits = X
ExpectedHits <- 4 #positive real number ? = equal expected value & var of X
lines(Hits ,ppois(Hits , lambda=ExpectedHits),type='s',col='skyblue')

plot(seq(-4,8,length=1000),dnorm(seq(-4,8,length=1000)),type='l')
lines(seq(-4,8,length=1000),dnorm(seq(-4,8,length=1000),mean=3),col='red')

x <- c(seq(0,1,length=30),seq(1,2,length=10))
cex <- .7
plot(x,df(x, df1=100, df2=100),type='l')
points(x,df(x, df1=1, df2=1),pch='?',cex=cex)
points(x,df(x, df1=5, df2=1),pch='?',cex=cex)
points(x,df(x, df1=100, df2=1),pch='?',cex=cex)
points(x,df(x, df1=100, df2=4),pch='4',cex=cex)

x <- c(seq(0,1,length=50),seq(1,8,length=100))
cex <- .7
n <- 1:9
plot(x,dchisq(x, df=n),pch=as.character(n),cex=cex)
abline(h=c(0,1))
```

Random Number Generation
=
```r
n <- 100
mX <- 14.3
stdev <- 3.6
x <- rnorm(n, mean=mX, sd=stdev)
h <- hist(x,xlab='Random V');rug(x)
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mX,sd=stdev)
yfit<-yfit*diff(h$mids[1:2])*length(x)
lines(xfit,yfit,col='red')

x <- runif(x, min=0, max=1)
h <- hist(x);rug(x)
xfit<-seq(min(x),max(x),length=40)
yfit<-dunif(xfit,min=0,max=1)
yfit<-yfit*diff(h$mids[1:2])*length(x)
lines(xfit,yfit,col='red')

par(mfrow=c(2,2))
plot.ecdf(rnorm(n),main='Normal')
plot.ecdf(rt(n,df=n-1),main='t')
plot.ecdf(runif(n, min=0, max=1),main='Uniform')
plot.ecdf(rbinom(n, 10, .5),main='Binomial')
```

Significance Testing
=
```r
mu <- 30
mX <- 32
n <- 16
s <- 10
qt(.95,df=n-1) #need GT 1.75305
TS <- (mX-mu)/(s/sqrt(n)) #achieved 0.8
pt(TS,df=n-1) #p(TS below achieved by chance) = 0.781901
1-pt(TS,df=n-1)  #p(mX>mu)=0.218099

S2 <- 10
s <- sqrt(10)
n <- 100
SE <- s/sqrt(n)
mu <- 30
CritVal <- qnorm(.95, mean=mu,sd=SE)
pnorm(CritVal, mean=mu,sd=SE)
```

Confidence Intervals
=
```r
mu <- 6          #population mean
mX <- 5        #sample mean
n <- 20        #sample size
S2 <- 4        #sample variance
s <- sqrt(S2)    #standard deviation
```
One-sided.
```r
SE <- qnorm(.95)*s/sqrt(n)
-Inf;mX+SE    #H0: mX is less than than mu
sv <- (mX-mu)/(s/sqrt(n))
pnorm(sv)
pnorm(mX,mean=mu,sd=s/sqrt(n)) #p(mX<mu) by chance=0.01267366
1-pnorm(mX,mean=mu,sd=s/sqrt(n)) #p(mX>mu) by chance =0.987
```
Two-sided.
```r
SE <- s/sqrt(n)
zMult <- qnorm(.975) #95%
mX-zMult*SE;mX+zMult*SE #H0: mX is not equal to mu
sv <- (mX-mu)/(s/sqrt(n))
2*pnorm(sv);2*pnorm(mX,mean=mu,sd=s/sqrt(n))
2*pt(sv,df=n-1)

qnorm(c(.025,.975),mean=mX,sd=s/sqrt(n)) #what to expect
sv <- (mX-mu)/(s/sqrt(n))
2*pt(sv,df=2)
SE <- s/sqrt(n)
tMult <- qt(.975,df=n-1)
mX-tMult*SE;mX+tMult*SE

#back-calc SE for two-sided
upper <- mX+zMult*SE
lower <- mX-zMult*SE
totalCIspace <- upper-lower
oneCIspace <- totalCIspace/2
#zval <- qnorm(.95) #90%
zval <- qnorm(.975) #95%
#zval <- qnorm(.995) #99%
SE2 <- oneCIspace/zval
SE;SE2

#back-calc SE for two-sided Z-dist
upper <- 1123
lower <- 1077
totalCIspace <- upper-lower
oneCIspace <- totalCIspace/2
zval <- qnorm(.95) #90%
SE2 <- oneCIspace/zval
newMult <- qnorm(.975) #95%
mean(upper,lower)-SE2*newMult 
mean(upper,lower)+SE2*newMult 

#back-calc SE for two-sided t-dist
n <- 9
upper <- 1123
lower <- 1077
totalCIspace <- upper-lower
oneCIspace <- totalCIspace/2
tval <- qt(.95,df=n-1) #90% t
SE2 <- oneCIspace/tval
newMult <- qt(.975,df=n-1) #95% t
mean(upper,lower)-SE2*newMult 
mean(upper,lower)+SE2*newMult 
```

Multiple Groups
=
```r
mX1 <- 10
mX2 <- 10.5
s1 <- 3
s2 <- 2.5
n1 <- 300
n2 <- 230

#common error variance
SE <- mean(s1,s2) * sqrt(1/n1 + 1/n2)
tv <- (mX1-mX2) / SE
tMult <- qt(0.975,df=n1+n2-2)
(mX1-mX2)-tMult*SE;(mX1-mX2)+tMult*SE
tv;pt(tv,df=n1+n2-2)*2

library(MASS)
test1 <- mvrnorm(n1,mu=mX1,Sigma=s1^2,empirical=T)
test2 <- mvrnorm(n2,mu=mX2,Sigma=s2^2,empirical=T)
t.test(test1,test2,var.equal=T)

#different error variance
SE <- sqrt(s1^2/n1 + s2^2/n2)
tv <- (mX1-mX2)/SE
approxdf <- (s1^2/n1+s2^2/n2)^2 / {(s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1)}
tMult <- qt(0.975,df=approxdf)
(mX1-mX2)-tMult*SE;(mX1-mX2)+tMult*SE
tv;pt(tv,df=approxdf)*2
t.test(test1,test2,var.equal=F)
tMult <- qt(0.995,df=approxdf) #expand to 99% CI
(mX1-mX2)-tMult*SE;(mX1-mX2)+tMult*SE
t.test(test1,test2,var.equal=F,conf.level=.99)

#paired (i.e. difference)
n <- 10
test1 <- mvrnorm(n,mu=10,Sigma=1,empirical=T)
test2 <- mvrnorm(n,mu=11,Sigma=1,empirical=T)
t.test(test2,test1,alternative='two.sided',paired=T)
dX <- test2-test1
SE <- sd(dX) / sqrt(n)
tv <- mean(dX)/SE  #or: sqrt(n) * mean(dX)/sd(dX)
tMult <- qt(0.975,df=n-1)
mean(dX)-tMult*SE;mean(dX)+tMult*SE
t.test(G1,alternative='two.sided',var.equal=T)
2*pt(abs(tv),df=n-1,lower.tail=F)
tv;qt(0.975,df=n-1)

# alternative='greater' that x has a larger mean than y.
t.test(test2,test1,alternative='greater',paired=T)
tMult <- qt(0.95,df=n-1) #one-sided quantile
mean(dX)-tMult*SE;Inf
pt(abs(tv),df=n-1,lower.tail=F)
tv;qt(0.95,df=n-1)

# alternative='less' that x has a smaller mean than y.
t.test(test2,test1,alternative='less',paired=T)
-Inf;mean(dX)+tMult*SE
pt(abs(tv),df=n-1,lower.tail=T)
tv;qt(0.95,df=n-1,lower.tail=F)
```

Combinations
=
Urn with 3 balls. Pick 2 balls from the urn with replacement, get all permutations.
```r
library(gtools)
x <- c('red', 'blue', 'black')
permutations(n=3,r=2,v=x,repeats.allowed=T)
```
| |V1|V2|
|---|---|---|
|1|black|black
|2|black|blue|
|3|black|red|
|4|blue|black
|5|blue|blue|
|6|blue|red|
|7|red|black
|8|red|blue|
|9|red|red|

Number of permutations without replacement.
```r
nrow(permutations(n=3,r=2,v=x,repeats.allowed=F)) #6
factorial(3)/factorial(3-2) #6
```

Pick 2 balls from the urn without replacement, get all permutations.
```r
permutations(n=3,r=2,v=x,repeats.allowed=F)
```
| |V1|V2|
|---|---|---|
|1|black|blue|
|2|black|red|
|3|blue|black|
|4|blue|red|
|5|red|black|
|6|red|blue|

All possible combinations.
```r
combinations(n=3,r=2,v=x,repeats.allowed=F)
choose(3,2)
factorial(3)/factorial(2)*factorial(3-2)
```
| |V1|V2|
|---|---|---|
|1|black|blue|
|2|black|red|
|3|blue|red|

Binomial
=

Two heads out of three flips
```r
permutations(n=2,r=3,v=c('T','H'),repeats.allowed=T)
```
| |V1|V2|V3|test|
|---|---|---|---|---|
|1|H|H|H||
|2|H|H|T|*|
|3|H|T|H|*|
|4|H|T|T||
|5|T|H|H|*|
|6|T|H|T||
|7|T|T|H||
|8|T|T|T||
```r
3/8 #0.375
choose(3,2)*(.5^2*(1-.5)^(3-2)) #0.375
dbinom(2,3,.5) # 0.375
pbirthday(2,classes=3,coincident=2) #0.3333333
```

Bernoulli
=
```r
1 * (.5^1*(1-.5)^(1-1)) #[1] 0.5
dbinom(1,1,.5) #[1] 0.5
permutations(n=2,r=3,v=c('T','H'),repeats.allowed=T)
combinations(n=2,r=2,v=c('T','H'),repeats.allowed=T)

library(ggplot2)
x = 0:19  # range of values to display in plot
# dbinom(x,n,theta) computes the pmf of
# Bin(n,theta) at x
y1 = dbinom(x, 5, 0.3)
y2 = dbinom(x, 5, 0.5)
y3 = dbinom(x, 5, 0.7)
y4 = dbinom(x, 20, 0.3)
y5 = dbinom(x, 20, 0.5)
y6 = dbinom(x, 20, 0.7)
D = data.frame(mass = c(y1, y2, y3, y4, y5, y6), x = x,
    n = 0, theta = 0)
D$n[1:60] = "$n=5$"
D$n[61:120] = "$n=20$"
D$theta[c(1:20, 61:80)] = "$\\theta=0.3$"
D$theta[c(21:40, 81:100)] = "$\\theta=0.5$"
D$theta[c(41:60, 101:120)] = "$\\theta=0.7$"
qplot(x, mass, data = D, main = "Binomial pmf", geom = "point",
    stat = "identity", facets = theta ~ n, xlab = "$x$",
    ylab = "$p_X(x)$") + geom_linerange(aes(x = x,
    ymin = 0, ymax = mass))

pbirthday(10,classes=365,coincident=2)

permutations(n=3,r=3,v=c('T','H','6cyl'),repeats.allowed=F)

#likelihood plot
x <- seq(0,1,by=.001)
y <- x^3*(1-x)^1
z <- y/max(y)
plot(x,z,type='l')
abline(h=1/8) #parameter values above 1/8 ref line are such that no other point
# is more than 8 times better supported given the data

1-pbinom(6,8,.5)
choose(8,7)*.5^7*.5^1+choose(8,8)*.5^8*.5^0
```

Power
=
```r
mu <- 6.5 #guess of true mean
mu <- seq(5,8,length=100) #guess range
mX <- 5
s <- 2
n <- 20
SE <- s/sqrt(n)
zMult <- qnorm(0.975)
lower <- mX-zMult*SE
upper <- mX+zMult*SE
zvl <- (lower-mu)/SE
zvu <- (upper-mu)/SE
pBTII <- pnorm(zvu)-pnorm(zvl) #from minimum in mu1 distribution to m0 crit region 
1-pBTII #power
plot(seq(4,9,length=100),1-pBTII,type='l')

#t-dist
tMult <- qt(0.975,df=n-1)
lower <- mX-tMult*SE
upper <- mX+tMult*SE
tvl <- (lower-mu)/SE
tvu <- (upper-mu)/SE
pBTII <- pt(tvu,df=n-1)-pt(tvl,df=n-1)
1-pBTII #power
lines(seq(4,9,length=100),1-pBTII,type='l')

#Non-centrality function
ncp <- (mu-mX)/SE # diff / s/sqrt(n)
tCrit <- qt(.975,df=n-1)
pBTII <- pt(tCrit,df=n-1,ncp=ncp)-pt(-tCrit,df=n-1,ncp=ncp)
1-pBTII #power

#Power without intervals: NCP, power.t.test, MC sims
mu <- .01 #mean brain volume loss, delta if mX=0
s <- .04
n <- c(20,50,100,140)
mX <- 0
SE <- s/sqrt(n)

#one-sided .95
ncp <- (mu-mX)/SE # diff / s/sqrt(n)
tCrit <- qt(.95,df=n-1) # .95
pBTII <- pt(tCrit,df=n-1,ncp=ncp)
1-pBTII #power
power.t.test(n=n,delta=mu-mX,sd=s,strict=T,type="one.sample",alt='one.sided')

#two-sided .975
ncp <- (mu-mX)/SE # diff / s/sqrt(n)
tCrit <- qt(.975,df=n-1) # .975
pBTII <- pt(tCrit,df=n-1,ncp=ncp)-pt(-tCrit,df=n-1,ncp=ncp)
1-pBTII #power
power.t.test(n=n,delta=mu-mX,sd=s,strict=T,type="one.sample",alt='two.sided')
power.t.test(delta=mu-mX,sd=s,power=.8359,strict=T,type="one.sample",alt='two.sided')
```

Monte Carlo Simulation
=
```r
mu <- .01 #mean brain volume loss
s <- .04
n <- 15
mX <- 0
SE <- s/sqrt(n)
MCsims <- 100000
z <- rnorm(MCsims)
xsq <- rchisq(MCsims,df=n-1)
tCrit <- qt(.95,df=n-1) #one-sided
mean(z + sqrt(n) * (mu-mX) / s > tCrit / sqrt(n-1) * sqrt(xsq))
tCrit <- qt(.95,df=n-1) #one-sided
ncp <- (mu-mX)/SE # diff / s/sqrt(n)
pBTII <- pt(tCrit,df=n-1,ncp=ncp)-pt(-tCrit,df=n-1,ncp=ncp)
1-pBTII #power
power.t.test(n=n,delta=mu-mX,sd=s,strict=T,type="one.sample",alt='one.sided')

tCrit <- qt(.975,df=n-1) #two-sided
mean(z + sqrt(n) * (mu-mX) / s > tCrit / sqrt(n-1) * sqrt(xsq))
ncp <- (mu-mX)/SE # diff / s/sqrt(n)
tCrit <- qt(.975,df=n-1) # .975  #two-sided
pBTII <- pt(tCrit,df=n-1,ncp=ncp)-pt(-tCrit,df=n-1,ncp=ncp)
1-pBTII #power
power.t.test(n=n,delta=mu-mX,sd=s,strict=T,type="one.sample",alt='two.sided')

(X-mu)/(s/sqrt(n))
# Z value > z-alpha - scaled population value
muN <- 30
qnorm(.95,mean=mu,sd=(s/sqrt(n))) - (muN-mu)/(s/sqrt(n))
```

Sig value alpha .05
=
```r
qnorm(.95) # 1.645
#(muN-mu)/(s/sqrt(n)) two events above muN
(32-30) /(4/sqrt(16)) # = 2
# so P(Z>1.645-2) or  P(Z> -0.355
pnorm(-0.355, lower.tail=F) #Beta 64% above, away from sample mean direction
#so
plot(25:33,pnorm(qnorm(.95) - ((32-25:33)/(4/sqrt(16))),
                 lower.tail=F),type='l')
plot(25:33,pnorm(qnorm(.95) - ((25:33-30)/(4/sqrt(16))),
                 lower.tail=F),type='l')
#value of 64% is bound for all values above 32

#Sample size calculation: What power? say, 80% power (or higher)
# .80 = Z > z1-alpha - (muN-mu)/(s/sqrt(n)) two events above muN
# set z1-alpha - (muN-mu)/(s/sqrt(n)) to z0.2 and solve for n
qnorm(.2) # -0.8416212
```
