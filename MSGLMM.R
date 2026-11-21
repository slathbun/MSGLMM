#
#
# List of regression coefficients for years education
#
betaEdyrpar<-list("betas[1,2]","betas[2,2]","betas[3,2]","betas[4,2]","betas[5,2]","betas[6,2]","betas[7,2]","betas[8,2]")
#
# List of regression coefficients for male gender
#
betaMalepar<-list("betas[1,3]","betas[2,3]","betas[3,3]","betas[4,3]","betas[5,3]","betas[6,3]","betas[7,3]","betas[8,3]")
#
# List of regression coefficients for age
#
betaAgepar<-list("betas[1,4]","betas[2,4]","betas[3,4]","betas[4,4]","betas[5,4]","betas[6,4]","betas[7,4]","betas[8,4]")
#
# List of parameters for black race
#
betaBlackpar<-list("betas[1,5]","betas[2,5]","betas[3,5]","betas[4,5]","betas[5,5]","betas[6,5]","betas[7,5]","betas[8,5]")
#
# List of parameters for APOE
#
betaAPOEpar<-list("betas[1,6]","betas[2,6]","betas[3,6]","betas[4,6]","betas[5,6]","betas[6,6]","betas[7,6]","betas[8,6]")

#
# Create parameter lists for coregionalization matrix A with p=8 dependent variables
#
# q=2 latent random fields
A82par<-list("A[1,1]","A[1,2]","A[2,1]","A[2,2]","A[3,1]","A[3,2]","A[4,1]","A[4,2]","A[5,1]","A[5,2]","A[6,1]","A[6,2]","A[7,1]","A[7,2]","A[8,1]","A[8,2]")
# q=3 latent random fileds
A83par<-list("A[1,1]","A[1,2]","A[1,3]","A[2,1]","A[2,2]","A[2,3]","A[3,1]","A[3,2]","A[3,3]","A[4,1]","A[4,2]","A[4,3]","A[5,1]","A[5,2]","A[5,3]","A[6,1]","A[6,2]","A[6,3]","A[7,1]","A[7,2]","A[7,3]","A[8,1]","A[8,2]","A[8,3]")
# q=4 latent random fields
A84par<-list("A[1,1]","A[1,2]","A[1,3]","A[1,4]","A[2,1]","A[2,2]","A[2,3]","A[2,4]","A[3,1]","A[3,2]","A[3,3]","A[3,4]","A[4,1]","A[4,2]","A[4,3]","A[4,4]","A[5,1]","A[5,2]","A[5,3]","A[5,4]","A[6,1]","A[6,2]","A[6,3]","A[6,4]","A[7,1]","A[7,2]","A[7,3]","A[7,4]","A[8,1]","A[8,2]","A[8,3]","A[8,4]")

#
# Create parameter lists  for correlations at lag zero
#
corrpar<-list("corr[1,2]","corr[1,3]","corr[1,4]","corr[1,5]","corr[1,6]","corr[1,7]","corr[1,8]","corr[2,3]","corr[2,4]","corr[2,5]","corr[2,6]","corr[2,7]","corr[2,8]","corr[3,4]","corr[3,5]","corr[3,6]","corr[3,7]","corr[3,8]","corr[4,5]","corr[4,6]","corr[4,7]","corr[4,8]","corr[5,6]","corr[5,7]","corr[5,8]","corr[6,7]","corr[6,8]","corr[7,8]")

#
# Create range parameter lists
#
# q=2 latent random fields
range2par<-list("range[1]","range[2]")
# q=3 latent random fields
range3par<-list("range[1]","range[2]","range[3]")
# q=4 latent random fields
range4par<-list("range[1]","range[2]","range[3]","range[4]")

#
# Carry out Wald's test of the null hypothesis that all regression coefficients associated
# with a predictor specified by "parm" are equal to zero against the alternative that at
# least one coefficient is not equal to zero. In addition to the results of Wald's test
# given at the end of the output, the output includes estimates of regression coefficients
# together with simultaneous confidence intervals and estimates of exponentiated
# regression coefficients and their confidence intervals
#
# stanobj = Stan object containing Markov chains generated using the Hamiltonian Markov
# Chain algorithm.
#
# parm = list object containing list of regression coefficients associated with a given
# predictor; e.g., betaEdyrpar for the regression coefficients for years education defined
# in the above.
#
# Example: To test the null hypothesis that black race has zero effect on all cognitive
# function variables, and to compute simultaneous 95% confidence intervals for exponentiated
# regression coefficients, call
#
# waldconf(MatMod4thin,betaBlackar)
#

waldconf<-function (stanobj,parm)
{
  data<-as.data.frame(stanobj,par=parm)
  # Compute sample covariance matrix of regression coefficients
  cov<-as.matrix(var(data))
  # Compute estimated regression coefficients
  beta<-as.vector(colMeans(data))
  # Compute Wald's statistic, its d.f., and p-value
  w<-t(beta)%*%solve(cov)%*%beta
  df<-length(beta)
  p<-1-pchisq(w,df)
  wald<-cbind(w,df,p)
  colnames(wald)<-c("Wald's Chi-Square","d.f.","p-value")
  # Compute SE of regression coefficient and lower and upper bounds of 95%
  # confidence intervals
  crit<-sqrt(qchisq(p=0.95,df=df))
  se=sqrt(as.vector(diag(cov)))
  lower<-beta-crit*se
  upper<-beta+crit*se
  conf<-cbind(beta,se,lower,upper)
  colnames(conf)<-c("Beta","SE","Lower","Upper")
  # Exponentiated values
  expbeta<-exp(beta)
  explower=exp(lower)
  expupper=exp(upper)
  exparm<-cbind(expbeta,explower,expupper)
  colnames(exparm)<-c("Exp(Beta)","Lower","Upper")
  result<-list(wald=list(),conf=list())
  result$wald<-wald
  result$conf<-conf
  result$exparm<-exparm
  return(result)
}

#
# The matrix of covariance functions is computed using matcovplot which depends on the
# functions corrfun and crosscovplot
#

#
# corrfun computes the correlation at lag distances specified by the vector "d" for an
# exponential correlation function with a specified range. For example,
d<-0.1*c(0:100)

corrfun<-function (d,range)
{
  n<-length(d)
  cov<-rep(0,n)
  for(i in 1:n)
  {cov[i]<-exp(-3*d[i]/range)}
  return(cov)
}

#
# crosscovplot plots the cross covariance function together with 95% confidence bands
# between dependent variables i and j.
# mcmc = Stan object containing Markov chains generated using the Hamiltonian Markov
# Chain algorithm.
# apar = Parameter list for coregionalization matrix
# rangepar = Parameter list for range parameters
# d = vector of distances at which the cross-covariance function is to be evaluated
# plotbands=TRUE if you wish to plot the 95% confidence bands; otherwise set plotbands=FALSE
# color = color of plotted cross-covariance function
#

crosscovplot <- function(mcmc,apar,rangepar,d,i,j,plotbands=TRUE, color="darkgreen", add=FALSE)
{
  mcmc_a<-as.data.frame(mcmc,par=apar)
  mcmc_range<-as.data.frame(mcmc,par=rangepar)
  n<-nrow(mcmc_a)
  m=length(d)
  q=ncol(mcmc_range)
  covs<-matrix(0.0,n,m)
  mincov <- rep(NA,m)
  maxcov <- rep(NA,m)
  for(k in 1:n)
  {
    for(l in 1:q)
    {
      temp<-as.vector(corrfun(d,mcmc_range[k,l]))
      # line below picking out the wrong element  <-  fix
      covs[k,]<-covs[k,]+mcmc_a[k, q*(i-1) + l]*mcmc_a[k, q*(j-1)+l]*t(temp)
    }
  }
  meancov<-apply(covs,2,mean)
  for(k in 1:m)
  {
    mincov[k]<-quantile(covs[,k],probs=0.025)
    maxcov[k]<-quantile(covs[,k],probs=0.975)
  }
  result<-list(covs=list(),meancov=list(),mincov=list(),maxcov=list())
  result$covs<-covs
  result$meancov<-meancov
  result$mincov<-mincov
  result$maxcov<-maxcov
  if (plotbands) {
    lines(d,maxcov,lty=2,lwd=1)
    lines(d,mincov,lty=2,lwd=1)
    polygon(c(d,rev(d)),c(mincov,rev(maxcov)),col="#6BD7AF",border=NA)

  }
  if (add) {
    lines(d,meancov,lty=1, type="l", col=color,lwd=1.3)
  } else {
    plot(d,meancov,lty=1, type="l", col=color, ylim=c(-1,1),lwd=1.3)
  }

  return(result)
}

#
# matcovplot plots the matrix of direct- and cross-covariance functions for three MSGLMM
# models. The 95% confidence band is only plotted for the last model.
# mod2, mod3, mod4 are the STAN outputs from the three models
# A2, A3, A4 = Parameter lists for the corregionalization matrices of the three
# models.
# range2, range3, range4 = Range lists for the three models
# d = vector of distances at which the cross-covariance function is to be evaluated
# p = number of dependent variables
#
# Example:
#
# matcovplot(MatMod2thin,MatMod3thin,MatMod4thin,A82par,A83par,A84par,range2par,range3par,range4par,d,8,d8=T)
#
#Here, MatMod2thin,MatMod3thin,MatMod4thin are STAN outputs for q=2,3, and 4 models, respectively
#

matcovplot <- function(mod2,mod3,mod4,A2,A3,A4,range2,range3,range4,d,p,d8)
{
  if (d8) {
    pdf("msglim_ccs_d8.pdf")
    d2 <- d3 <- d4 <- 1:200/20
  } else {
    #pdf("msglim_ccs.pdf")
    d2 <- d3 <- d4 <- 1:200/25
  }
  par(mfrow=c(p,p))
  par(mar=c(0,0,0,0))
  par(oma=c(4,4,2,2))
  par(mgp=c(3,1,0))

  md <- max(d2)

  # manually set ranges for y axes
  yrans <- list(c(-0.1, 0.4), c(-0.1, 1.0),
                c(-0.1, 0.6), c(-0.2, 1.0),
                c(-0.3, 0.6), c(-0.4, 1.0),
                c(-0.2, 0.6),c(-0.2, 0.5))
  # Variable names for dependent variables
  varnames <- c("MOCA", "WordRec", "ReyFig", "BNT", "Fluency", "Executive", "Clock", "Logical")

  i <- j <- 1
  for (i in 1:p) {
    for (j in 1:p) {
      if (j > i) {
        # no plots on upper triangle, except panel 1,2 has the legend
        plot.new()
        plot.window(c(0,1),c(0,1))

        #  plot legend
        if (i == 3 & j == 5) {
          plot.window(c(0,1),c(0,1))
          # manually enter contents of legend for three models
          legend(x=0.0,y=0.6, c("q=4", "q=3", "q=2"),
                 col=c("darkgreen", "purple", "orange"), cex=0.8,
                 lty=c(1,2,2), lwd=c(2,2,2), xpd=0)
        }

        # add variable names to edges of upper triangle
        if (j == i+1) {
          corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          #par(xpd = TRUE) #Draw outside plot area
          text(x = 0.1, y = mean(corners[3:4]),
               varnames[j-1], srt = 270, font=par("font.main"), cex=0.85)

          corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          #par(xpd = TRUE) #Draw outside plot area
          text(x = mean(corners[1:2]), y = 0.1,
               varnames[j], srt = 0, font=par("font.main"), cex=0.85)
        }
      } else {

        # we're on lower triangle, set up plots with or without
        # axis ticks, depending if we're on the left col or bottom row
        if (i == p & j == 1) {
          plot(1:10, type="n",xlim=c(0,md), ylim=yrans[[i]], las=1)
        } else if (i == p) {
          # supress only y axis
          plot(1:10, type="n",xlim=c(0,md), ylim=yrans[[i]], yaxt="n")
        } else if (j == 1) {
          # supress only x axis
          plot(1:10, type="n",xlim=c(0,md), ylim=yrans[[i]], xaxt="n", las=1)
        } else {
          # supress axes
          plot(1:10, type="n",xlim=c(0,md), ylim=yrans[[i]], xaxt="n", yaxt="n")
        }

        # add varnames to first column, pth row
        if (i==1 & j == 1)
          text(varnames[1], x=md/2, y=yrans[[1]][2]*1.15, xpd=NA, cex=0.85, font=par("font.main"))
        if (i==p & j == p)
          text(varnames[p], x=md*1.15, y=mean(yrans[[p]]), xpd=NA, cex=0.85, font=par("font.main"),srt=270)

        # add lines from cross covs
        crosscovplot(mod4,A4,range4,d,i,j, add=TRUE, color="darkgreen")
        crosscovplot(mod3,A3,range3,d,i,j, add=TRUE, color="purple", plotbands=FALSE)
        crosscovplot(mod2,A2,range2,d,i,j, add=TRUE, color="orange", plotbands=FALSE)
        crosscovplot(mod4,A4,range4,d,i,j, add=TRUE, color="darkgreen",plotbands=FALSE)
        abline(0,0, col="gray20", xpd=0, lty=3, lwd=1.2)

      }
    }
  }
  mtext(text="Distance (km)", side=1, outer=TRUE, line=3)
  dev.off()
}

#
# corrsum obtains the estimate of the correlations between the dependent variables at lag
# zero.
# mcmc = Stan object containing Markov chains generated using the Hamiltonian Markov
# Chain algorithm.
# corrpar = parameter list for correlations
# p = number of dependent variables
#
# Example:
# corrsum(MatMod4thin,corrpar,8)
#

corrsum <- function (mcmc,corrpar,p)
{
  data<-as.data.frame(mcmc,par=corrpar)
  corrmeanvec<-as.vector(colMeans(data))
  corrmat<-matrix(0.0,p,p)
  k=0
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if(i==j) {corrmat[i,j]=1.0}
      if(j>i)
      {
        k=k+1
        corrmat[i,j]=corrmeanvec[k]
        corrmat[j,i]=corrmeanvec[k]
      }
    }
  }
  return(corrmat)
}

#
# practrange estimates the practical range of spatial correlation and variance of the latent
# random field (varepsilon) for dependent variable i. 95% confidence intervals are provided
# for both estimates.
# mcmc = Stan object containing Markov chains generated using the Hamiltonian Markov
# Chain algorithm.
# apar = list of coregionalization matrix parameters
# rangepar = list of range parameters
#
# Example
# practrange(MatMod4thin,A84par,range4par,1)
#

practrange <- function(mcmc,apar,rangepar,i)
{
  mcmc_a<-as.data.frame(mcmc,par=apar)
  mcmc_range<-as.data.frame(mcmc,par=rangepar)
  n<-nrow(mcmc_a)
  q=ncol(mcmc_range)
  pract_range=rep(0.0,n)
  var_epsilon=rep(0.0,n)
  for(j in 1:n)
  {
    dmin<-mcmc_range[j,1]
    dmax<-mcmc_range[j,q]
    dmid=0.5*(dmin+dmax)
    pract_range[j]=dmid
    diff<-dmax-dmin
    cov0<-0
    for(k in 1:q)
    {
      cov0<-cov0+mcmc_a[j, q*(i-1) + k]*mcmc_a[j, q*(i-1)+k]
    }
    var_epsilon[j]<-cov0
    while(diff>0.0001)
    {
      covs=0.0
      for(k in 1:q)
      {
        covs<-covs+mcmc_a[j, q*(i-1) + k]*mcmc_a[j, q*(i-1)+k]*exp(-3.0*dmid/mcmc_range[j,k])
      }
      corr<-covs/cov0;
      if(corr>0.05) {dmin<-dmid}
      else {dmax<-dmid}
      dmid=0.5*(dmin+dmax)
      pract_range[j]=dmid
      diff<-dmax-dmin
    }
  }
  mean_range<-mean(pract_range)
  lower_range<-quantile(pract_range,probs=0.025)
  upper_range<-quantile(pract_range,probs=0.975)
  mean_var_epsilon<-mean(var_epsilon)
  lower_var_epsilon<-quantile(var_epsilon,probs=0.025)
  upper_var_epsilon<-quantile(var_epsilon,probs=0.975)
  stats<-cbind(mean_range,lower_range,upper_range,mean_var_epsilon,lower_var_epsilon,upper_var_epsilon)
  colnames(stats)<-c("Practical Range","Lower","Upper","Var_Epsilon","Lower","Upper")
  result<-list(stats=list())
  result$stats<-stats
  return(result)
}
