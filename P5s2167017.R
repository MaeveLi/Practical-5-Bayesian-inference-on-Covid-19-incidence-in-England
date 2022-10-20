######## Practical 5 - Maeve Li (Minqing Li) s2167017 ########
#### Overview: This R code along with the model.jags file implements the MCMC Covid model stated in the practical pdf file
#### and samples the expected number of daily deaths and the inferred number of daily infections from the MCMC model.

#### The stated model models the log of daily infection x_i (i.e., n_i=exp(x_i)) and x_i~N(2*x_(i-2)-x_(i-2),tau).
#### We assume that x_1~N(0,tau_0=0.01) and x_2~N(x_1,tau). We also specify that precision tau~Gamma(4,0.04).
#### n_i is related to the available number of deaths. We have that the fatal disease duration follows a 
#### log normal distribution (3.235,0.4147^2). For the vector of the expected number of deaths m, there's the relationship 
#### that m = B %*% n, where B is a matrix that has the value B_ij = PI_D(i-j) when i>=j and B_ij = 0 when i<j. 
#### Finally, we have that the observed number of deaths on day i, y_i ~ Poi(m_i).

#### Using rjags and coda and the provide data of y, samples for n_i and m_i are generated using 10000 iterations.  
#### Several diagnostic plots are plotted with brief comments on their meanings, and a recommended iteration number for this model 
#### with a desired burn-in was provided. 

#### A final plot is produced where the x axis represent the number of days in 2020, 
#### which has the posterior mean for the n_i trajectory (discarding the last 20 days of data), 
#### the 95% credible interval of n_i overlaid on the plot, the observed number of deaths, the posterior mean for 
#### the m_i trajectory, and the first day of UK lockdown.

## Reference: https://r-charts.com/evolution/area-between-lines/
##            https://stats.stackexchange.com/questions/203281/number-of-markov-chain-monte-carlo-samples
library(rjags)
library(coda)

#### Import data
y <- c(1,2,0,2,2,0,4,4,1,9,14,20,22,27,40,46,66,63,105,103,149,159,204,263,326,353,360,437,498,576,645,647,
       700,778,743,727,813,900,792,740,779,718,699,646,686,639,610,571,524,566,486,502,451,438,386,380,345,341,
       325,313,306,268,251,259,252,266,256,215,203,196,166,183,162,180,172,167,138,160,144,153,150,122,128,116,
       133,139,122,124,117,92,83,94,109,111,83,86,83,80,73,69)
append0 <- rep(0,20)
newy <- append(append0,y)  ## append 20 zeros to the start of the data

#### Create B matrix
nlen = length(newy)
B = matrix(0,nlen,nlen) ## initialize B matrix
for (i in 1:nlen){
  ## B's values are calculated using these two loops
  for (j in 1:nlen){
    if (j <= i)
      B[i,j] = dlnorm((i-j),meanlog=3.235,sdlog=0.4147)
    else
      B[i,j] = 0
  }
}

#### Jags sampling
mod <- jags.model("model.jags",data=list(y=newy,N=nlen,B=B))
sam.coda <- coda.samples(mod,c("n","m"),n.iter=10000)  ## sample for the n_i and m_i values

#### Calculate the posterior means for n and m and 95% credible interval for n
mean_n <- apply(sam.coda[[1]][,121:240],2,mean)
mean_m <- apply(sam.coda[[1]][,1:120],2,mean)
CI_n <- apply(sam.coda[[1]][,121:240],2,quantile,prob=(c(0.025,.975)))

#### Diagnostics
## To create a sensible small subset, 4 nodes are randomly picked to monitor the mixing status (m[50],m[100],n[30],n[80])
plot(sam.coda[[1]][,c(50,100,150,200)],type="l") 
## There's no long initializing period of one or more parameters drifting systematically in one direction 
## for n[30] and m[50] and they seem to mix well and we speculated it's because they're at the first half of the data, 
## however m[100] and n[80] do seem to have regular fluctuations around the mean and don't seem to have converged
effectiveSize(sam.coda[[1]][,c(50,100,150,200)])
## Effective sizes are all very low relative to the number of iterations - all under 500, 
## especially low for m[100] and n[80] as they are at the end of the data
crosscorr(sam.coda[[1]][,c(50,100,150,200)])
## We can see that n[80] and m[100] are moderately correlated (which accounts for the low effective sample size)
acfplot(sam.coda[[1]][,c(50,100,150,200)],aspect=1)
## From the autocorrelation plot, we can see that mixing is quite slow, especially for n[80] and m[100]

## There is clear correlation between the nodes and effective samples sizes are quite low, so I would recommend to
## first reduce the correlation between certain n & m nodes and then run a second iteration.
## I would recommend an iteration number of 20000 to provide more effective samples and more space for nodes at the
## end of the data to converge. Low effective sample sizes may result from a large burn-in period, so 
## I would recommend discarding around 8000 burn-in since the largest effective sample size is around 600
## and we can see from the trace plot that data do seem to become more stable from around 8000 (without sudden jumps and decreases)

#### Final plot preparation
## y, m
ystart <- julian(as.Date("2020-3-02"),origin=as.Date("2019-12-31"))[1] ## calculate the julian day from the start of the year for the start of y 
newystart = ystart - 20; newyend = newystart + 119 
ydate <- seq(newystart,newyend,1) ## form julian dates for y
yd <- data.frame(ydate,newy); md <- data.frame(ydate,mean_m)
## n, credible intervals for n
newn <- mean_n[1:100]; newCI <- CI_n[,1:100] ## discard the last 20 days of n as they're unstable
CI_upper <- newCI[2,]; CI_lower <- newCI[1,]
nstart = newystart; nend = nstart + 99
ndate <- seq(nstart,nend,1)
nd <- data.frame(ndate,newn)
nCIupper <- data.frame(ndate,CI_upper); nCIlower <- data.frame(ndate,CI_lower)
## the first day of lock down
lockdown <- julian(as.Date("2020-3-24"),origin=as.Date("2019-12-31"))[1] ## calculate the julian day of lockdown

#### The Final Plot
## Data starts from 12/2/2020
plot(yd,xlim=c(0,newyend),ylim=c(0,max(nCIupper)),
     xlab="Day in 2020",ylab="Number of Incidence",pch=20,
     main="Daily Deaths and Bayesian Inferrence on Daily Infections in England") ## daily deaths
grid(5,5,col="gray66")
lines(md,col="red",lwd=2) ## posterior means for m_i
lines(nd,col="orange",lwd=2) ## posterior means for n_i
lines(nCIupper,col="cornflowerblue")
lines(nCIlower,col="cornflowerblue") ## 95% credible interval for n_i
polygon(c(nCIupper$ndate,rev(nCIupper$ndate)),c(nCIlower$CI_lower,rev(nCIupper$CI_upper)),
        density=20,col ="cornflowerblue",border=NA)  ## illustrate the 95% credible interval
abline(v=lockdown,lwd=1,lty=5) ## mark the lockdown date
text(x=lockdown,y=1800,label='Lockdown Day',pos=4,cex=0.8)
legend("topleft",bg="transparent",
       legend=c("daily hospital deaths","posterior mean for daily infection","posterior mean for daily death"),
       fill=c("black","orange","red"),cex=0.6,border=NA)
legend("topright",bg="transparent",legend="95% CI for daily infection",
       fill="cornflowerblue",density=50,cex=0.6,border=NA)