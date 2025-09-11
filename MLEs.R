#load the required packages
library(reshape2)
library(ggplot2)
library(ggpubr)
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

### Negative log-likelihood for truncated data
NLL.Trunc = function(pars, data){
  #Extract parameters from the vector
  df = 9
  theta = pars[1]
  cutoff=0.5
  t1<-data[,1]
  t2<-data[,2]
  # Calculate Negative Log-Likelihood
  nllik <- -sum(dbvt(t1=t1,t2=t2,df=c(df,df),theta=theta,
                   cutoff=cutoff,truncation="right", log = TRUE))
  return(nllik)
}
mle.trun = optim(par = c(theta = 0.5),
                 method = "Brent",lower = -pi/2,upper = pi/2,
                 fn = NLL.Trunc, data = T.trun,
                 control = list(parscale = c(theta = 0.5)))
mle.trun$par

theta.grid<-seq(-1.47,1.47,by=0.01)
t1 <- T.trun[,1]
t2 <- T.trun[,2]

ll.theta <- sapply(theta.grid, function(th) {
  sum(dbvt(t1 = t1, t2 = t2, df = df, theta = th,
           cutoff = 0.5, truncation = "right",log = TRUE))
})

df <- data.frame("theta"=theta.grid, "loglik"=ll.theta)

### Plot log likelihood
ll.plot <- ggplot(data = df, aes(theta, loglik)) +
  geom_line() +
  geom_vline(xintercept = theta.grid[which.max(ll.theta)], lty = 2) +
  labs(
    title = expression("Case 1:" ~ n[1] == ~ n[2] == 9),
    x     = expression(theta ~ "(In radians)"),
    y     = "Log likelihood"
  )
ll.plot

#### MAXIMUM LIKELIHOOD ESTIMATION ####
#### UNEQUAL DEGREE OF FREEDOMS ####

#Generate data from bi-variate t distribution.
set.seed(123)
n<-10000
df<- c(15,9)
theta<-pi/6
cutoff<-0.5
T.Pop<-rbvt(n=n,df=df,theta=theta)
#Select those observations whose baseline measurement is greater than 0.05
T.trun<-T.Pop[(T.Pop[,1]>=0.5),]

### Negative log-likelihood for truncated data
NLL.Trunc.ue = function(pars, data) {
  # Extract parameters from the vector
  n1<-15
  n2<-9
  theta = pars[1]
  t1<-data[,1]
  t2<-data[,2]
  # Calculate Negative Log-Likelihood
  nllik <- -sum(dbvt(t1,t2,df=c(n1,n2),theta=theta,
                   cutoff = 0.5,truncation = "right",log = TRUE))
  return(nllik)
}
mle = optim(par = c(theta = 0.5),
            method = "Brent",lower = -pi/2,upper = pi/2,
            fn = NLL.Trunc.ue, data = T.trun,
            control = list(parscale = c(theta = 0.5)))
theta.est<-mle$par
theta.est

theta.grid<-seq(-1.47,1.47,by=0.01)
t1 <- T.trun[,1]
t2 <- T.trun[,2]

ll.theta <- sapply(theta.grid, function(th) {
  sum(dbvt(t1 = t1, t2 = t2, df = df, theta = th,
           cutoff = 0.5, truncation = "right", log = TRUE))
})

# Plot log likelihood (Truncated)
df2<-data.frame("theta"=theta.grid,"loglik"=ll.theta)
ll.plot2<-ggplot(data =  df2,aes(theta,loglik))+geom_line()+
  geom_vline(xintercept=theta.grid[which.max(ll.theta)],lty=2)+
  labs(
    title = expression("Case 2:" ~ n[1] == 15 ~ "," ~ n[2] == 9),
    x = expression(theta ~ "(In radians)"),
    y = "Log likelihood"
  )
ll.plot2

ggarrange(ll.plot,ll.plot2,ncol = 2, nrow = 1,labels = c("(a)","(b)"),
          font.label = list(size = 11, color = "black", face = "plain", family = NULL))

