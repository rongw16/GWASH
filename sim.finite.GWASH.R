#-----------------------------------------------------------------------
# Finite sample simulations
#-----------------------------------------------------------------------
rm(list=ls())
library(MASS)
library(psych)
library(copula)

m1.m2_est=function(n,m,S.tilde,mu23.version,my.mask,I2){

	m1=1;
	if (mu23.version==1){ #Eq(23)
		m2=sum(S.tilde^2)/m - (m-1)/(n-1)  
		return(list(m1=m1,m2=m2))
	}else{ #Eq(24)
		S.tilde.mask = S.tilde * my.mask
		m2 = sum(S.tilde.mask^2)/m - (m1^2*I2/m)/(n-1)
		
		return(list(m1=m1,m2=m2,S.tilde.mask=S.tilde.mask))

	}
	
}

m3_est=function(n,m,S.tilde,m1,m2,mu23.version,I2,I3){
	if(mu23.version==1){ #Eq (27)
		S2=S.tilde%*%S.tilde;
		m3=sum(S2*S.tilde)/m-3*(m-1)*m2/(n-1)-(m-1)*(m-2)/(n-1)^2;
	}else{ #Eq (28)	: S.tilde is actually S.tilde.mask;
		
		m3=sum(diag(S.tilde%*%S.tilde%*%S.tilde))/m - 3*m2*(I2/m)/(n-1) - (I3/m)/(n-1)^2
		
	}
	
	return(list(m3=m3))
}

h2.MM.est.indep=function(z,n,m){ 
	s2=mean(z^2);
    h2.MM =(m/n)*(s2 - 1)
	return(h2.MM)
}

h2.MM.variance.CI.indep=function(n,m,h2){
	psi2=2*(m/n + 2*h2 -h2^2);	
	upper_h2=h2+1.96*sqrt(psi2/n)
	lower_h2=h2-1.96*sqrt(psi2/n)
	return(list(variance=psi2,upper_h2=upper_h2,lower_h2=lower_h2))
	
}

h2.MM.est.dep=function(u,m1,m2,n,m,t.score){ 
	if(t.score == 1){ #use t-score, adjust
		u2=((n-1)/(n-2))*u^2/(1+u^2/(n-2)) #Eq.30
		s2=mean(u2)
	}else{ #use z-score
		s2=mean(u^2);	
	}

	h2.MM=(m*m1/(n*m2))*(s2 - m1) #m=1 anyways

	return(h2.MM)
}

h2.MM.variance.CI.dep=function(n,m,h2,m1,m2,m3){
	psi2=2*(m*m1^2/(n*m2) + 2*m1*m3*h2/m2^2 -h2^2);
	upper_h2=h2+1.96*sqrt(psi2/n)
	lower_h2=h2-1.96*sqrt(psi2/n)
	return(list(variance=psi2,upper_h2=upper_h2,lower_h2=lower_h2))
	
}

# LD regression
LD.reg=function(n,m,z,S){ 
	
	ell.hat = colSums(S^2)	
	my.lm=coef(summary(lm(z^2 - 1 ~ I(ell.hat*n/m - 1) - 1)))
	return(list(coef=my.lm[,1],se=my.lm[, 2]))
	
	
}


#Dicker's estimate 
h2.II.est=function(n,m,y,x){ #eq.11
	S_II=t(x)%*%x/(n-1);

	m1=(1/m)*tr(S_II)

	m2=(1/m)*sum(S_II^2)-(m*m1^2)/(n-1);

	h2.II.estimate=(m*m1^2/(n*m2))*(sum((t(x)%*%y)^2)/(m*m1*t(y)%*%y) -1)
	
	return(h2.II.estimate)		
}


#######################
# Parameters
h2 = 0.5
m = 200 #number of SNP
n = 100  #number of subject
nsim = 10 #number of dataset simulated


#######################
#Don't change 
RANDOM.b = F  #fixed b
RANDOM.X = T  #always
SCALE.X = T   #always
SCALE.y = T   #always 
#######################
b.dist=2; #1:normal beta; 2:beta is the mixture of normal and binominal 

X.indep = F 
X.norm = F; 

AR=T; #Autoregressive correlation matrix in finite sample simulation

n.subdiag=3; #when rho=0.2, n.subdiag=3 (see paper Table 1 )
diag.1 = NA;
p=NA;
rho=NA;
null.b.prop=NA;
mu2.true=1;mu3.true=1;

if(!X.norm){ #binom
	p=0.1;q=1-p;
	rep.num=10; #for covariance matrix;
	rep.num.2=1000; #for correlation matrix, approximate true mu2 and mu3;
}else{
	diag.1=T; #identity matrix, diagnal =1
}

if(AR){
	rho=0.2; #rho can take 0, 0.2, 0.4, 0.6, 0.8,etc
}

if(b.dist==2){ #if beta is mixed
	null.b.prop=0.9;
}


###########################
#varince-covariance matrix
if (X.norm){
	if(X.indep){
		SIGMA=diag(1,m,m)
	}else{ 
		if(AR){
			Sigma = rho^abs(outer(1:m, 1:m, "-"))
			if (diag.1){
				SIGMA=Sigma;
			}else{
				SIGMA = sqrt(diag(1:m)) %*% Sigma %*% sqrt(diag(1:m))
			}
		}else{
			#Dicker 2014, only generate once SIGMA
			W= array(rnorm(2*m^2), dim=c(2*m, m))
			SIGMA=(t(W)%*%W)/(2*m) #only generate one SIGMA
		}
		
	}
}else{ #binom
	if(X.indep){
		SIGMA=diag(2*p*q,m,m)
	}else{		
		cov.sum=array(0,dim=c(m,m))
		for ( i in 1:rep.num){
			normal.cop=normalCopula(rho,dim=m,dispstr="ar1")
			u <- rCopula(n, normal.cop)
			correlated.binom=array(qbinom(u,size=2,prob=p), dim=c(n, m))
			cov.sum=cov.sum+cov(correlated.binom)					
		}
		SIGMA=cov.sum/rep.num #SIGMA is an average of rep.num
	}

}
SIGMA[1:10,1:10]

#################################################
#This section may be time consuming if m is large
if(X.norm){
	mu2.true=sum(SIGMA^2)/m;
	mu3.true=sum(SIGMA%*%SIGMA*SIGMA)/m
}else{ #binom: obtain "true" correlation matrix;
	if (!X.indep){
		cor.sum=array(0,dim=c(m,m))
		for ( i in 1:rep.num.2){
			normal.cop=normalCopula(rho,dim=m,dispstr="ar1")
			u <- rCopula(n, normal.cop)
			correlated.binom=array(qbinom(u,size=2,prob=p), dim=c(n, m))
			cor.sum=cor.sum+cor(correlated.binom)			
		}
		SIGMA.cor=cor.sum/rep.num.2
		print(SIGMA.cor[1:10,1:10])
		
		mu2.true=sum(SIGMA.cor^2)/m;
		mu3.true=sum(SIGMA.cor%*%SIGMA.cor*SIGMA.cor)/m
	}
	
}

mu2.true
mu3.true
###################################################
#sigma2.b = h2/tr(SIGMA); #only need this if RANDOM.b;
sigma2.eps =1-h2


if(!RANDOM.b){ #fixed b
		if(b.dist==1){
			b_star=rnorm(m);		 	
						
		}else{
			b_star=rep(0,m)
			null.b.indx=sort(sample(seq(1:m),size=m*null.b.prop,replace=F))			
			b_star[-null.b.indx]=rnorm(round((1-null.b.prop)*m))
			
		}
		b=b_star*sqrt(h2)/sqrt(t(b_star)%*%SIGMA%*%b_star)
}

b[1:100]


#--------------------------------------
z.mat=matrix(0,nrow=m,ncol=nsim)
t.mat=z.mat;
mu2.v1=rep(0, nsim);
mu2.v2=mu2.v1;

mu3.v1=mu2.v1;
mu3.v2=mu3.v1;

h2.LD = rep(0, nsim) 
h2.LD.se=h2.LD;

h2.II = rep(0, nsim) 

h2.MM.1 = rep(0, nsim)
h2.MM.2.v1 = h2.MM.1
h2.MM.2.v2 = h2.MM.1

#h2.GWAS.1=h2.MM.1;
#h2.GWAS.2.v1=h2.MM.1;
#h2.GWAS.2.v2=h2.MM.1;

h2.var.1= rep(0, nsim)
h2.upper.1=h2.var.1
h2.lower.1=h2.var.1

h2.var.2.v1= h2.var.1
h2.upper.2.v1=h2.var.1
h2.lower.2.v1=h2.var.1

h2.var.2.v2= h2.var.1
h2.upper.2.v2=h2.var.1
h2.lower.2.v2=h2.var.1

time_used.mm.1<-NULL;
time_used.mm.2.v1<-NULL;
time_used.mm.2.v2<-NULL;
time_used.LD<-NULL;
time_used.II<-NULL;
time_used.S<-NULL;
time_used.m3.v1<-NULL;
time_used.m3.v2<-NULL;
time_used.I2<-NULL;
time_used.I3<-NULL;

####################
ptm.I2<-proc.time()
mask = abs(outer(1:m, 1:m, "-")) <= n.subdiag
I2 = sum(mask) - m
array.mask=array(mask, dim=c(m, m));
ptm2.I2<-proc.time()
time_used.I2=ptm2.I2-ptm.I2;
####################
ptm.I3<-proc.time()
ind = expand.grid(1:m, 1:m, 1:m)
valid.ind = (ind[,1] != ind[,2]) & (ind[,1] != ind[,3]) & (ind[,2] != ind[,3]) & (abs(ind[,1] - ind[,2]) <= n.subdiag) & (abs(ind[,1] - ind[,3]) <= n.subdiag) & (abs(ind[,2] - ind[,3]) <= n.subdiag)
I3 = sum(valid.ind)
ptm2.I3<-proc.time()
time_used.I3=ptm2.I3-ptm.I3;		
####################		
		

ptm<-proc.time()
for (i in 1:nsim){
	# if(RANDOM.b & b.dist==1 & X.indep){ #won't use;
  		# b = rnorm(m, sd = sqrt(sigma2.b))   # random normal b for independ X 
	# }
	
	if(RANDOM.X){
		if(X.indep){ #identity matrix;
			if(X.norm){ #X normal distribution, identity matrix
  				X = array(rnorm(n*m), dim=c(n, m))   # random X, each of nsim dataset has different values;
  			}else{ #X binom 				
  				X =array(rbinom(n*m,size=2,prob=p), dim=c(n, m))
  			}
  		}else{ #SIGMA is NOT Identitiy matrix;
  			if(X.norm){ #generate correlated normal X;
  				X = mvrnorm(n, rep(0, m), SIGMA)
  			}else{ #generate correlated Binominal X;
  				normal.cop=normalCopula(rho,dim=m,dispstr="ar1")
				u <- rCopula(n, normal.cop)
				X=array(qbinom(u,size=2,prob=p), dim=c(n, m))
  			}
  		}
  	}
  	
	X=scale(X,center = TRUE, scale = F) #center only
	
	eps = rnorm(n, sd = sqrt(sigma2.eps))
	y = rep(0,n)
	y = X %*% b + eps
	y=scale(y,center = TRUE, scale = F)
	
	# GWAS
	if(SCALE.X) X.s = scale(X) else X.s = X  #center & scale
	if(SCALE.y) y.s = scale(y) else y.s = y
	
  	b.GWAS = (1/colSums(X^2))*t(X) %*% y   			
  	s2.GWAS = colSums((matrix(y,nrow=n,ncol=m) - X*matrix(b.GWAS,nrow=n,ncol=m,byrow=T))^2)/(n-2)

	z.mat[,i] = t(X.s) %*% y.s/sqrt(n-1); 
	t.mat[,i] = sqrt(colSums(X^2)) * b.GWAS / sqrt(s2.GWAS)
   
    ptm.S<-proc.time();
    S=(t(X.s)%*%(X.s))/(n-1); #sample covariance matrix,time consuming; only do it once;
    ptm2.S<-proc.time();
    time_used.S= rbind(time_used.S,(ptm2.S - ptm.S))
   
    ###################
	#h2.MM: Independ 
	###################
	ptm.mm.1 <- proc.time() 
	h2.MM.1[i]=h2.MM.est.indep(z.mat[,i],n,m)
	ptm2.mm.1<-proc.time()
	time_used.mm.1=rbind(time_used.mm.1,(ptm2.mm.1 - ptm.mm.1))
	
		
	est.obj=h2.MM.variance.CI.indep(n,m,h2.MM.1[i])
	h2.var.1[i]=est.obj$variance;
	h2.upper.1[i]=est.obj$upper_h2;
	h2.lower.1[i]=est.obj$lower_h2;
		
	#h2.GWAS.1[i]=h2.MM.est.indep(t.mat[,i],n,m)
		
		
	if(!X.indep){ #non-Identity	
		###################
		#version 1: Full 
		###################
		ptm.mm.v1 <- proc.time()
		m.obj=m1.m2_est(n,m,S,1,0,0);
		m1=m.obj$m1
		mu2.v1[i]=m.obj$m2

		h2.MM.2.v1[i]=h2.MM.est.dep(z.mat[,i],m1,mu2.v1[i],n,m,0)
		ptm2.mm.v1<-proc.time()
		time_used.mm.2.v1=rbind(time_used.mm.2.v1,(ptm2.mm.v1 - ptm.mm.v1))
		
		ptm.m3.v1<-proc.time();
		mu3.v1[i]=m3_est(n,m,S,m1,mu2.v1[i],1,0,0)$m3
		ptm2.m3.v1<-proc.time();
		time_used.m3.v1=rbind(time_used.m3.v1,(ptm2.m3.v1-ptm.m3.v1))
		
		est.obj=h2.MM.variance.CI.dep(n,m,h2.MM.2.v1[i],m1,mu2.v1[i],mu3.v1[i])
		h2.var.2.v1[i]=est.obj$variance;
		h2.upper.2.v1[i]=est.obj$upper_h2;
		h2.lower.2.v1[i]=est.obj$lower_h2;
		
		#h2.GWAS.2.v1[i]=h2.MM.est.dep(t.mat[,i],m1,mu2.v1[i],n,m,1)
		
		
		###################
		#version 2: partial
		###################
		ptm.mm.v2 <- proc.time()
		m.obj=m1.m2_est(n,m,S,2,array.mask,I2);
		mu2.v2[i]=m.obj$m2
		S.mask=m.obj$S.tilde.mask		


		h2.MM.2.v2[i]=h2.MM.est.dep(z.mat[,i],m1,mu2.v2[i],n,m,0)
		ptm2.mm.v2<-proc.time()
		time_used.mm.2.v2=rbind(time_used.mm.2.v2,(ptm2.mm.v2 - ptm.mm.v2))
		
		ptm.m3.v2<-proc.time();
		mu3.v2[i]=m3_est(n,m,S.mask,m1,mu2.v2[i],2,I2,I3)$m3
		ptm2.m3.v2<-proc.time();
		time_used.m3.v2=rbind(time_used.m3.v2,(ptm2.m3.v2-ptm.m3.v2))
		
		est.obj=h2.MM.variance.CI.dep(n,m,h2.MM.2.v2[i],m1,mu2.v2[i],mu3.v2[i])
		h2.var.2.v2[i]=est.obj$variance;
		h2.upper.2.v2[i]=est.obj$upper_h2;
		h2.lower.2.v2[i]=est.obj$lower_h2;
		
		
		#h2.GWAS.2.v2[i]=h2.MM.est.dep(t.mat[,i],m1,mu2.v2[i],n,m,1) #not use

		
	}
	

	
	#h2.LD 
	ptm.LD<-proc.time()
	my.LD.reg=LD.reg(n,m,z.mat[,i],S)
	h2.LD[i]=my.LD.reg$coef;
	ptm2.LD<-proc.time()
	time_used.LD=rbind(time_used.LD,(ptm2.LD - ptm.LD))
	h2.LD.se[i]=my.LD.reg$se;
	
	
	#h2.Dickers
	ptm.II<-proc.time()
	h2.II[i]=h2.II.est(n,m,y,X) #y and X are centered only, not scaled;
	ptm2.II=proc.time()
	time_used.II=rbind(time_used.II,(ptm2.II - ptm.II))
}

#process ends;
ptm2<-proc.time()

#overall process time needed 
time_used=ptm2 - ptm
time_used

#-----------------------------------------------------------------------
#This section for independent X only
mean(h2.MM.1) #X.indep
sd(h2.MM.1)
if(X.indep){
	sqrt(2/n*(1+m/n-(1-h2)^2)) 
}
mean(h2.lower.1[!is.na(h2.lower.1)])
mean(h2.upper.1[!is.na(h2.upper.1)])

#bias
mean(h2.MM.1) - h2

# mean(h2.GWAS.1)
# sd(h2.GWAS.1)

mean(h2.var.1)
mean(sqrt(h2.var.1/n))
#-----------------------------------------------------------------------


if(!X.indep){
	
	print("###############################")
	print("v1")
	print("###############################")
	print(paste("h2.MM.2.v1:",mean(h2.MM.2.v1)))
	print(paste("S.E:",sd(h2.MM.2.v1)))
	
	#mu2
	print(paste("mu2.v1:",mean(mu2.v1)))
	print(paste("S.E:",sd(mu2.v1)))
	
	#mu3
	print(paste("mu3.v1:",mean(mu3.v1)))
	print(paste("S.E:",sd(mu3.v1)))
	

	#variance
	print(paste("asymptotic variance:",mean(h2.var.2.v1)))
	print(paste("asymptotic S.E:",mean(sqrt(h2.var.2.v1/n))))
	
	print(paste("lower:",mean(h2.lower.2.v1[!is.na(h2.lower.2.v1)])))
	print(paste("upper:",mean(h2.upper.2.v1[!is.na(h2.upper.2.v1)])))

	#bias
	print(paste("bias:",mean(h2.MM.2.v1) - h2))

	# print(paste("h2.GAWS.2.v1:",mean(h2.GWAS.2.v1)))
	# print(paste("S.E:",sd(h2.GWAS.2.v1)))
	
	print("###############################")
	print("v2")
	print("###############################")
	print(paste("h2.MM.2.v2:",mean(h2.MM.2.v2)))
	print(paste("S.E:",sd(h2.MM.2.v2)))
	
	#mu2
	print(paste("mu2.v2:",mean(mu2.v2)))
	print(paste("S.E:",sd(mu2.v2)))
	
	#mu3
	print(paste("mu3.v2:",mean(mu3.v2)))
	print(paste("S.E:",sd(mu3.v2)))
	

	#variance
	print(paste("asymptotic variance:",mean(h2.var.2.v2)))
	print(paste("asymptotic S.E:",mean(sqrt(h2.var.2.v2/n))))
	
	print(paste("lower:",mean(h2.lower.2.v2[!is.na(h2.lower.2.v2)])))
	print(paste("upper:",mean(h2.upper.2.v2[!is.na(h2.upper.2.v2)])))

	#bias
	print(paste("bias:",mean(h2.MM.2.v2) - h2))

	# print(paste("h2.GAWS.2.v2:",mean(h2.GWAS.2.v2)))
	# print(paste("S.E:",sd(h2.GWAS.2.v2)))
	
}

#LD.regression
mean(h2.LD)
sd(h2.LD)
mean(h2.LD.se)

#DICKER'S 
mean(h2.II)
sd(h2.II)

mu2.true
mu3.true;

