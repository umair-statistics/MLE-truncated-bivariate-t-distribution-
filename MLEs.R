#load the required packages
library(reshape2)
library(ggplot2)
library(ggpubr)
library(latex2exp)
#devtools::install_github('umair-statistics/RTM')
library(RTM)
theme_set(theme_bw(base_size = 11))
theme_update(plot.title = element_text(hjust = 0.5))

#### MAXIMUM LIKELIHOOD ESTIMATION ####
  #### EQUAL DEGREE OF FREEDOMS ####


#Generate data from bivariate t distribution.
set.seed(123)
n<-10000
df<- c(9,9)
theta<-pi/6
cutoff<-0.5
T.Pop<-rbvt(n=n,df=df,theta=theta)
#Select those observations whose baseline measurement is greater than 0.05
T.trun<-T.Pop[(T.Pop[,1]>=0.5),]

NLL.Trunc = function(pars, data){
  #Extract parameters from the vector
  df = 9
  theta = pars[1]
  cutoff=0.5
  t1<-data[,1]
  t2<-data[,2]
  # Calculate Negative Log-Likelihood
  -sum(log(dbvt(t1=t1,t2=t2,df=c(df,df),theta=theta,cutoff=cutoff,truncation="right")))
}

mle.trun = optim(par = c(theta = 0.5),
                 method = "Brent",lower = 0,upper = 1,
                 fn = NLL.Trunc, data = T.trun,
                 control = list(parscale = c(theta = 0.5)))

mle.trun$par

theta.grid<-seq(0,pi/2,by=0.01)
ll.theta<-NA

for (i in 1:length(theta.grid)) {
  t1<-T.Pop[,1]
  t2<-T.Pop[,2]
  ll.theta[i] <-  sum(log(dbvt(t1=t1,t2=t2,df=df,theta=theta.grid[i],cutoff = 0.5, truncation = "right")))
}

### Plot log likelihood

theta_ll_df<-data.frame("theta"=theta.grid,"loglik"=ll.theta)

ll.plot <-ggplot(theta_ll_df,aes(theta,loglik))+geom_line()+geom_vline(xintercept=theta.grid[which.max(ll.theta)],lty=2)+
  xlab(TeX(r"($\theta$ (In radians))"))+ylab("Log likelihood")

ll.plot


#### MAXIMUM LIKELIHOOD ESTIMATION ####
  #### UNEQUAL DEGREE OF FREEDOMS ####


### Negative log-likelihood for truncated data

#Generate data from bi-variate t distribution.
set.seed(123)
n<-10000
df<- c(15,9)
theta<-pi/6
cutoff<-0.5
T.Pop<-rbvt(n=n,df=df,theta=theta)
#Select those observations whose baseline measurement is greater than 0.05
T.trun<-T.Pop[(T.Pop[,1]>=0.5),]

NLL.Trunc.ue = function(pars, data) {
  # Extract parameters from the vector
  n1<-15
  n2<-9
  theta = pars[1]
  t1<-data[,1]
  t2<-data[,2]
  # Calculate Negative Log-Likelihood
  suppressWarnings(-sum(log(dbvt(t1,t2,df=c(n1,n2),theta=theta,cutoff = 0.5,truncation = "right"))))
}


mle = optim(par = c(theta=0.5),
            method = "Brent",lower = 0,upper = 1,
            fn = NLL.Trunc.ue, data =T.trun,
            control = list(parscale = c(theta=0.5)))

theta.est<-mle$par

theta.grid<-seq(0,pi/2,by=0.01)
ll.theta<-NA

for (i in 1:length(theta.grid)) {
  df<-c(15,9)
  ll.theta[i] <- sum(log(dbvt(t1=T.Pop[,1],t2=T.Pop[,2],df=df,theta =  theta.grid[i],cutoff = 0.5,truncation = "right")))
  print(i)
}

# Plot log likelihood (Truncated)

thetall<-data.frame("theta"=theta.grid,"loglik"=ll.theta)
ll.plot2<-ggplot(thetall,aes(theta,loglik))+geom_line()+geom_vline(xintercept=theta.grid[which.max(ll.theta)],lty=2)+
  xlab(TeX(r"($\theta$ (Radians))"))+ylab("Log likelihood")

ggarrange(ll.plot,ll.plot2,ncol = 2, nrow = 1,labels = c("(a)","(b)"),
          font.label = list(size = 11, color = "black", face = "plain", family = NULL))
