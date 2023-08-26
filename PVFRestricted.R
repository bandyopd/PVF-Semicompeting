## R code for estimation in semi-competing endpoints with PVF frailty: Restricted Model 
## Written by: Elizbeth Bedia
## Edited by: Dipankar Bandyopadhyay, on 08/26/2023



############################ CALL DATA INTO R #####################################################
library(survival)
data(colon)
dados1=colon[colon$etype==1,]
dados2=colon[colon$etype==2,]
dados=colon
Y1=dados1$time/365 # years
d1=dados1$status
Y2=dados2$time/365 # years
d2=dados2$status
x1=as.factor(dados1$rx)
#x9 =as.factor(dados1$extent)	
#x10 =dados1$surg
#x11=dados2$node4

X=model.matrix(~1+x1)
nc=ncol(X)
nrow(X)


############## LIKELIHOOD FUNCTION FOR PVF FRAILTY (RESTRICTED MODEL) ############
LogLike.PVFRestC = function(par){
  gama = par[1]
  theta = par[2]
  beta1 = par[3]
  beta2 = par[4]
  nc=ncol(X)
  b1 = par[5:(nc+4)]
  b2 = par[(nc+5):(2*nc+4)]
  alpha1=exp(X%*%b1)
  alpha2=exp(X%*%b2)
  
  S1=alpha1*Y1^(beta1)+alpha2*Y2^(beta2)
  S2=alpha1*Y1^(beta1)+alpha2*Y1^(beta2)
  a1=(1+(theta*S1)/(1-gama))
  a2=(1+(theta*S2)/(1-gama))
  K1=exp((1-gama)/(theta*gama)*(1-(a1)^gama))
  K2=exp((1-gama)/(theta*gama)*(1-(a2)^gama))
  
  temp = (alpha1*beta1*Y1^(beta1-1))^d1*
    (alpha2*beta2*Y2^(beta2-1))^d2*
    (theta+a1^gama)^(d1*d2)*
    a1^(d1*(gama-1-d2))*K1^d1*
    a2^((gama-1)*d2*(1-d1))*K2^(1-d1)
  lglik = sum(log(temp))
  return(lglik)
}
################################################################################
par=c(0.9,1,.01,.01,rep(-0.01,2*nc))
LogLike.PVFRestC(par)
fit=optim(par,LogLike.PVFRestC,control=list(fnscale=-1),method="L-BFGS-B",hessian = T)
#lower = c(-1, rep(0.0001,3),rep(-Inf,3)),upper = c(1,rep(Inf, length(par)-1)))
EMV=fit$par
EP=sqrt(diag(solve(-fit$hessian)))
L=EMV-qnorm(0.975)*EP
U=EMV+qnorm(0.975)*EP
Stat=abs(EMV/EP)
Saida=cbind(EMV, EP,Stat,L, U)
rownames(Saida)=c("gama","log_theta",paste("log_beta_", 1:2, sep = ""),paste("ph1_", 1:nc, sep = ""), paste("ph2_", 1:nc, sep = ""))
print(Saida, digits = 2)

############## For profile Likelihood function ############
LogLike.PVFRest = function(par){
  theta = par[1]
  beta1 = par[2]
  beta2 = par[3]
  nc=ncol(X)
  b1 = par[4:(nc+3)]
  b2 = par[(nc+4):(2*nc+3)]
  alpha1=exp(X%*%b1)
  alpha2=exp(X%*%b2)
  
  S1=alpha1*Y1^(beta1)+alpha2*Y2^(beta2)
  S2=alpha1*Y1^(beta1)+alpha2*Y1^(beta2)
  a1=(1+(theta*S1)/(1-gama))
  a2=(1+(theta*S2)/(1-gama))
  K1=exp((1-gama)/(theta*gama)*(1-(a1)^gama))
  K2=exp((1-gama)/(theta*gama)*(1-(a2)^gama))
  
  temp = (alpha1*beta1*Y1^(beta1-1))^d1*
    (alpha2*beta2*Y2^(beta2-1))^d2*
    (theta+a1^gama)^(d1*d2)*
    a1^(d1*(gama-1-d2))*K1^d1*
    a2^((gama-1)*d2*(1-d1))*K2^(1-d1)
  lglik = sum(log(temp))
  return(lglik)
}
gama=0.5
par=c(1,.01,.01,.01,rep(-.01,2*nc))
LogLike.PVFRest(par)

################################################################################
rgamma=seq(-0.2,0.2,0.003) # Parameters of interest
rgamma
parPVFR=c(1,.01,.01,.01,rep(-.01,2*nc))
ng=length(rgamma)
EMVRpvf=NULL
llR=numeric(ng)
for(k in 1:ng){
  gama=rgamma[k]  
  LogLike.PVFRest(parPVFR)
  fitPVFR=optim(parPVFR,LogLike.PVFRest,control=list(fnscale=-1),method="L-BFGS-B",hessian = T,
                lower = c(rep(0.0001,4),rep(-Inf,3)),upper = rep(Inf, length(parPVFR)))
  EMVRpvf=rbind(EMVRpvf,fitPVFR$par)
  llR[k]=fitPVFR$value
}
llR
max(llR)
gamaestR=rgamma[65]
gamaestR
plot(rgamma[1:107],llR[1:107], type ='l',lwd=2, xlab="Gamma", ylab="Profile likelihood function")
abline(v=-0.008,col='black',lty=2, lwd= 2)
EMVRper=EMVRpvf[65,]
EPRpvf=sqrt(diag(solve(-fitPVFR$hessian)))
Lpvfr=EMVRper-qnorm(0.975)*EPRpvf
Upvfr=EMVRper+qnorm(0.975)*EPRpvf
Statpvfr=abs(EMVper/EPpvf)
Saidapvfr=cbind(EMVRper, EPRpvf,Statpvfr,Lpvfr, Upvfr)
rownames(Saidapvfr)=c("theta",paste("beta_", 1:3, sep = ""),paste("ph1_", 1:nc, sep = ""), paste("ph2_", 1:nc, sep = ""), paste("ph3_", 1:nc, sep = ""))
print(Saidapvfr, digits = 2)

###############################################################################
############## LIKELIHOOD FUNCTION FOR GAMMA FRAILTY ##########################
LogLike.GamaRest = function(par){
  theta = par[1]
  beta1 = par[2]
  beta2 = par[3]
  nc=ncol(X)
  b1 = par[4:(nc+3)]
  b2 = par[(nc+4):(2*nc+3)]
  alpha1=exp(X%*%b1)
  alpha2=exp(X%*%b2)
  
  temp = (alpha1*beta1*Y1^(beta1-1))^d1*(alpha2*beta2*Y2^(beta2-1))^d2*(theta + 1)^(d1*d2)*
    (1 + theta*(alpha1*Y1^(beta1)+alpha2*Y2^(beta2)))^(-d1-d2- 1/theta)
  lglik = sum(log(temp))
  return(lglik)
}
##########################################################################################
pargr=c(1,.01,.01,rep(-0.01,2*nc))
LogLike.GamaRest(pargr)
fitgr=optim(pargr,LogLike.GamaRest,control=list(fnscale=-1),method="L-BFGS-B",hessian = T)
      # lower = c(rep(0.0001,3),rep(-Inf,3)),upper = rep(Inf, length(par)))
EMVgr=fitgr$par
EPgr=sqrt(diag(solve(-fitgr$hessian)))
Lgr=EMVgr-qnorm(0.975)*EPgr
Ugr=EMVgr+qnorm(0.975)*EPgr
Statgr=abs(EMVgr/EPgr)
Saidagr=cbind(EMVgr, EPgr,Statgr,Lgr, Ugr)
rownames(Saidagr)=c("log_theta",paste("log_beta_", 1:2, sep = ""),paste("ph1_", 1:nc, sep = ""), paste("ph2_", 1:nc, sep = ""))
print(Saidagr, digits = 2)

############################## LIKELIHOOD FUNCTION FOR IG FRAILTY #####################
LogLike.IGRest = function(par){
  theta = exp(par[1])
  beta1 = (par[2])
  beta2 = (par[3])
  nc=ncol(X)
  b1 = par[4:(nc+3)]
  b2 = par[(nc+4):(2*nc+3)]
  alpha1=exp(X%*%b1)
  alpha2=exp(X%*%b2)
  
  K1=sqrt(1+2*theta*(alpha1*Y1^(beta1)+alpha2*Y2^(beta2)))
  K2=sqrt(1+2*theta*(alpha1*Y1^(beta1)+alpha2*Y1^(beta2)))
  
  temp = (alpha1*beta1*Y1^(beta1-1))^d1*(alpha2*beta2*Y2^(beta2-1))^d2*(K1+theta)^(d1*d2)*
    exp(1/theta*(1-K1))^d1*exp(1/theta*(1-K2))^(1-d1)*K1^(-d1*(1+2*d2))*K2^(-d2*(1-d1))
  lglik = sum(log(temp))
  return(lglik)
  
}
#################################################################################
parIGr=c(1,.01,.01,rep(-0.01,2*nc))
LogLike.IGRest(parIGr)
fitIGr=optim(parIGr,LogLike.IGRest,control=list(fnscale=-1),method="L-BFGS-B",hessian = T)
           #lower = c(rep(0.0001,3),rep(-Inf,3)),upper = rep(Inf, 9))
EMVigr=fitIGr$par
EPigr=sqrt(diag(solve(-fitIGr$hessian)))
Ligr=EMVigr-qnorm(0.975)*EPigr
Uigr=EMVigr+qnorm(0.975)*EPigr
Statigr=abs(EMVigr/EPigr)
Saidaigr=cbind(EMVigr, EPigr,Statigr,Ligr, Uigr)
rownames(Saidaigr)=c("theta",paste("beta_", 1:2, sep = ""),
                     paste("ph1_", 1:nc, sep = ""), paste("ph2_", 1:nc, sep = ""))
print(Saidaigr, digits = 2)
########################### AIC and BIC method  #################################
AIC= function(Logf, K){
  AIC=2*K-2*Logf
  AIC
}
BIC = function(Logf, K, n){
  BIC=log(n)*K-2*Logf
  BIC
}
n=nrow(X)
gama=gamaestR
Kpvfr=length(fit$par)
Kgr=length(fitgr$par)
Kigr=length(fitIGr$par)

Logfpvfr=LogLike.PVFRestC(EMV) 
Logfgr=LogLike.GamaRest(EMVgr)
Logfigr=LogLike.IGRest(EMVigr)

AICpvfr=AIC(Logfpvfr,Kpvfr) 
AICgr=AIC(Logfgr,Kgr)
AICigr=AIC(Logfigr,Kigr)

BICpvfr=BIC(Logfpvfr,Kpvfr,n)
BICgr=BIC(Logfgr,Kgr,n)
BICigr=BIC(Logfigr,Kigr,n)

Saidapvfr=cbind(Logfpvfr, AICpvfr,BICpvfr)
rownames(Saidapvfr)=c("PVF Model")
print(Saidapvfr, digits = 7)

Saidagr=cbind(Logfgr, AICgr,BICgr)
rownames(Saidagr)=c("Gamma Model")
print(Saidagr, digits = 7)

Saidaigr=cbind(Logfigr, AICigr,BICigr)
rownames(Saidaigr)=c("IG Model")
print(Saidaigr, digits = 7)

