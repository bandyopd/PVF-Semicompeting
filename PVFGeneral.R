## R code for profile likelihood based estimation in semi-competing endpoints with PVF frailty: General Model 
## Written by: Elizbeth Bedia
## Edited by: Dipankar Bandyopadhyay

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
##################### PROFILE LIKELIHOOD APPROACH  #################################

############## LIKELIHOOD FUNCTION FOR PVF FRAILTY #################################
LogLikelihood.PVFCPer = function(par){
  theta = par[1]
  beta1 = par[2]
  beta2 = par[3]
  beta3 = par[4]
  nc=ncol(X)
  b1 = par[5:(nc+4)]
  b2 = par[(nc+5):(2*nc+4)]
  b3 = par[-c(1:(2*nc+4))]
  alpha1=exp(X%*%b1)
  alpha2=exp(X%*%b2)
  alpha3=exp(X%*%b3)
  
  S1=alpha1*Y1^(beta1)+alpha2*Y1^(beta2)+(alpha3*Y2^(beta3)-alpha3*Y1^(beta3))
  S2=alpha1*Y1^(beta1)+alpha2*Y1^(beta2)
  a1=(1+(theta*S1)/(1-gama))
  a2=(1+(theta*S2)/(1-gama))
  K1=exp((1-gama)/(theta*gama)*(1-(a1)^gama))
  K2=exp((1-gama)/(theta*gama)*(1-(a2)^gama))
  
  temp = (alpha1*beta1*Y1^(beta1-1))^d1*
    (alpha2*beta2*Y2^(beta2-1))^(d2*(1-d1))*
    (alpha3*beta3*Y2^(beta3-1))^(d1*d2)*
    (theta+a1^gama)^(d1*d2)*
    a1^(d1*(gama-1-d2))*K1^d1*
    a2^((gama-1)*d2*(1-d1))*K2^(1-d1)
  lglik = sum(log(temp))
  return(lglik)
}

################################################################################
mgamma=seq(-0.99,0.05,0.004) # Parameters of interest
parPVF=c(1,.01,.01,.01,rep(-.01,3*nc))
ng=length(mgamma)
EMVpvf=NULL
ll=numeric(ng)
for(k in 1:ng){
  gama=mgamma[k]  
  LogLikelihood.PVFCPer(parPVF)
  fitPVF=optim(parPVF,LogLikelihood.PVFCPer,control=list(fnscale=-1),
               method="L-BFGS-B",hessian = t)
  EMVpvf=rbind(EMVpvf,fitPVF$par)
  ll[k]=fitPVF$value
}
ll
max(ll)
gamaest=mgamma[195]
gamaest
plot(mgamma,ll, type ='l',lwd=2, xlab="Gamma", ylab="Profile likelihood function")
abline(v=-0.214, col='blue',lty=2, lwd= 2)
EMVper=EMVpvf[195,]
EPpvf=sqrt(diag(solve(-fitPVF$hessian)))
Lpvf=EMVper-qnorm(0.975)*EPpvf
Upvf=EMVper+qnorm(0.975)*EPpvf
Statpvf=abs(EMVper/EPpvf)
Saidapvf=cbind(EMVper, EPpvf,Statpvf,Lpvf, Upvf)
rownames(Saidapvf)=c("theta",paste("beta_", 1:3, sep = ""),paste("ph1_", 1:nc, sep = ""), paste("ph2_", 1:nc, sep = ""), paste("ph3_", 1:nc, sep = ""))
print(Saidapvf, digits = 2)


############## LIKELIHOOD FUNCTION FOR GAMMA FRAILTY #################################
LogLikelihood.powerGamaC = function(par){
  theta = par[1]
  beta1 = par[2]
  beta2 = par[3]
  beta3 = par[4]
  nc=ncol(X)
  b1 = par[5:(nc+4)]
  b2 = par[(nc+5):(2*nc+4)]
  b3 = par[-c(1:(2*nc+4))]
  alpha1=exp(X%*%b1)
  alpha2=exp(X%*%b2)
  alpha3=exp(X%*%b3)
  
  temp = (alpha1*beta1*Y1^(beta1-1))^d1*(alpha2*beta2*Y2^(beta2-1))^(d2*(1-d1))*(alpha3*beta3*Y2^(beta3-1))^(d1*d2)*(theta + 1)^(d1*d2)*
    (1 + theta*(alpha1*Y1^(beta1)+alpha2*Y1^(beta2) + alpha3*Y2^(beta3)-alpha3*Y1^(beta3)))^(-d1-d2- 1/theta)
  lglik = sum(log(temp))
  return(lglik)
}
##########################################################################################
parg=c(1,.01,.01,.01,rep(-0.01,3*nc))
LogLikelihood.powerGamaC(parg)
fitg=optim(parg,LogLikelihood.powerGamaC,control=list(fnscale=-1),method="L-BFGS-B",hessian = t,
           lower = c(rep(0.0001,3),rep(-Inf,9)),upper = rep(Inf, length(parg)))
EMVg=fitg$par
EPg=sqrt(diag(solve(-fitg$hessian)))
Lg=EMVg-qnorm(0.975)*EPg
Ug=EMVg+qnorm(0.975)*EPg
Statg=abs(EMVg/EPg)
Saida=cbind(EMVg, EPg,Statg,Lg, Ug)
rownames(Saida)=c("log_theta",paste("log_beta_", 1:3, sep = ""),paste("ph1_", 1:nc, sep = ""), paste("ph2_", 1:nc, sep = ""), paste("ph3_", 1:nc, sep = ""))
print(Saida, digits = 2)

############################## LIKELIHOOD FUNCtION FOR IG FRAILtY #####################
LogLikelihood.powerIGC = function(par){
  theta = par[1]
  beta1 = par[2]
  beta2 = par[3]
  beta3 = par[4]
  nc=ncol(X)
  b1 = par[5:(nc+4)]
  b2 = par[(nc+5):(2*nc+4)]
  b3 = par[-c(1:(2*nc+4))]
  alpha1=exp(X%*%b1)
  alpha2=exp(X%*%b2)
  alpha3=exp(X%*%b3)
  
  K1=sqrt(1+2*theta*(alpha1*Y1^(beta1)+alpha2*Y1^(beta2)
                     + (alpha3*Y2^(beta3)-alpha3*Y1^(beta3))))
  K2=sqrt(1+2*theta*(alpha1*Y1^(beta1)+alpha2*Y1^(beta2)))
  
  temp = (alpha1*beta1*Y1^(beta1-1))^d1*(alpha2*beta2*Y2^(beta2-1))^(d2*(1-d1))*
    (alpha3*beta3*Y2^(beta3-1))^(d1*d2)*(K1+theta)^(d1*d2)*
    exp(1/theta*(1-K1))^d1*exp(1/theta*(1-K2))^(1-d1)*K1^(-d1*(1+2*d2))*K2^(-d2*(1-d1))
  lglik = sum(log(temp))
  return(lglik)
  
}
#################################################################################
parIG=c(1,.01,.01,.01,rep(-0.01,3*nc))
LogLikelihood.powerIGC(parIG)
fitIG=optim(parIG,LogLikelihood.powerIGC,control=list(fnscale=-1),method="L-BFGS-B",hessian = t,
            lower = c(rep(0.0001,3),rep(-Inf,9)),upper = rep(Inf, length(parIG)))
EMVig=fitIG$par
EPig=sqrt(diag(solve(-fitIG$hessian)))
Lig=EMVig-qnorm(0.975)*EPig
Uig=EMVig+qnorm(0.975)*EPig
Statig=abs(EMVig/EPig)
Saidaig=cbind(EMVig, EPig,Statig,Lig, Uig)
rownames(Saidaig)=c("log_theta",paste("log_beta_", 1:3, sep = ""),paste("ph1_", 1:nc, sep = ""), paste("ph2_", 1:nc, sep = ""), paste("ph3_", 1:nc, sep = ""))
print(Saidaig, digits = 2)
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
gama=gamaest
Kpvf=length(fitPVF$par)
Kg=length(fitg$par)
Kig=length(fitIG$par)

Logfpvf=LogLikelihood.PVFCPer(EMVper) 
Logfg=LogLikelihood.powerGamaC(EMVg)
Logfig=LogLikelihood.powerIGC(EMVig)

AICpvf=AIC(Logfpvf,Kpvf) 
AICGam=AIC(Logfg,Kg)
AICig=AIC(Logfig,Kig)

BICpvf=BIC(Logfpvf,Kpvf,n)
BICGam=BIC(Logfg,Kg,n)
BICig=BIC(Logfig,Kig,n)

Saidapvf=cbind(Logfpvf, AICpvf,BICpvf)
rownames(Saidapvf)=c("PVF Model")
print(Saidapvf, digits = 7)

Saidag=cbind(Logfg, AICGam,BICGam)
rownames(Saidag)=c("Gamma Model")
print(Saidag, digits = 7)

Saidaig=cbind(Logfig, AICig,BICig)
rownames(Saidaig)=c("IG Model")
print(Saidaig, digits = 7)
