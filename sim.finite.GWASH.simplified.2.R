#-----------------------------------------------------------------------
# Finite sample simulations
#-----------------------------------------------------------------------
rm(list=ls())
library(MASS)
library(psych)
library(copula)

##############################################
# function sections
##############################################
#mu2 estimation Eq (22) or Eq (27)
mu2_est=function(n,m,S.tilde,mu23.version,my.mask,I2){ 
	mu1=1; #mu1 is always 1
	if (mu23.version==1){ 
		#Eq(22), full sample
		mu2=sum(S.tilde^2)/m - (m-1)/(n-1)  
		return(list(mu1=mu1,mu2=mu2))

	}else{ 
		
		#Eq(27), mu.2,I2, subset
		S.tilde.mask = S.tilde * my.mask
		mu2 = sum(S.tilde.mask^2)/m - (mu1^2*I2/m)/(n-1)
		
		return(list(mu1=mu1,mu2=mu2,S.tilde.mask=S.tilde.mask))

	}
	
}
##############################################
#mu3 estimation Eq (30) or Eq (31)
mu3_est=function(n,m,S.tilde,mu1,mu2,mu23.version,I2,I3){ 
	if(mu23.version==1){ #Eq (30), full sample
		S2=S.tilde%*%S.tilde;
		mu3=sum(S2*S.tilde)/m-3*(m-1)*mu2/(n-1)-(m-1)*(m-2)/(n-1)^2;
	}else{ 
		
		#Eq (31): mu.3,I3, subset		
		mu3=sum(diag(S.tilde%*%S.tilde%*%S.tilde))/m - 3*mu2*(I2/m)/(n-1) - (I3/m)/(n-1)^2
		
	}
	
	return(list(mu3=mu3))
}

##############################################	
#GWASH h^2 estimate,Eq (20)
h2.GWASH.est.dep=function(u,mu1,mu2,n,m){  
	s2=mean(u^2); #Eq. (21)
	h2.GWASH=(m*mu1/(n*mu2))*(s2 - mu1) #mu=1 always #Eq.20, GWASH

	return(h2.GWASH)
}

##############################################
# Asymptotic 95% CI and variance of GWASH h^2 estimate;
h2.GWASH.variance.CI.dep=function(n,m,h2,mu1,mu2,mu3){ #Based on Eq (25), mu1=1 always;
	psi2=2*(m*mu1^2/(n*mu2) + 2*mu1*mu3*h2/mu2^2 -h2^2);
	upper_h2=h2+1.96*sqrt(psi2/n)
	lower_h2=h2-1.96*sqrt(psi2/n)
	return(list(variance=psi2,upper_h2=upper_h2,lower_h2=lower_h2))
	
}


##############################################
# Information needed for simulating dataset
##############################################
h2 = 0.2
m = 250 #number of SNP
n = 100  #number of subject
nsim = 50 #number of dataset simulated, if nsim=1, only 1 dataset is generated


b.dist=2; #1: beta is from N(0,1); 2:beta is the mixture of N(0,1) and 0's 
X.norm = F; #X.norm=T if X is normal distribution; X.norm = F if X is binomial distribution.
n.subdiag=3; #equivalent to q in Table 1; for example, when rho=0.2, n.subdiag=3 (see paper Table 1 )
rho=0.2; #rho>0 is the AR correlation coefficient

if(!X.norm){ #if X is binomial, 
	p=0.1;q=1-p; #binomial success and failure rates, could be other values;
	rep.num=10; #for obtaining covariance matrix;
}

if(b.dist==2){ #if beta is mixed
	null.b.prop=0.9; #90% of beta are 0, the other 10% beta are from N(0,1)  
}

##############################################
#Generate varince-covariance matrix for X
#assuming AR
##############################################
if (X.norm){ #X is normal distribution
	Sigma = rho^abs(outer(1:m, 1:m, "-")) #correlation matrix
	SIGMA = sqrt(diag(1:m)) %*% Sigma %*% sqrt(diag(1:m)) #fake a covariance matrix

}else{ #X is binomial	
	cov.sum=array(0,dim=c(m,m))
	for ( i in 1:rep.num){
		normal.cop=normalCopula(rho,dim=m,dispstr="ar1")
		u <- rCopula(n, normal.cop)
		correlated.binom=array(qbinom(u,size=2,prob=p), dim=c(n, m))
		cov.sum=cov.sum+cov(correlated.binom)					
	}
	SIGMA=cov.sum/rep.num #SIGMA is an average of rep.num

}


##############################################
#Generate beta
##############################################
if(b.dist==1){ #beta is from N(0,1)
	b_star=rnorm(m);		 							
}else{ #beta is mixture
	b_star=rep(0,m)
	null.b.indx=sort(sample(seq(1:m),size=m*null.b.prop,replace=F))			
	b_star[-null.b.indx]=rnorm(round((1-null.b.prop)*m))
			
}

#normalize beta
b=b_star*sqrt(h2)/sqrt(t(b_star)%*%SIGMA%*%b_star)
sigma2.eps =1-h2

##############################################
#Prepare for subset for computational efficiency, 
#when using full dateset, this isn't needed;
#n.subdiag is equivalent to q in Table 1
##############################################
mask = abs(outer(1:m, 1:m, "-")) <= n.subdiag
I2 = sum(mask) - m
array.mask=array(mask, dim=c(m, m));

ind = expand.grid(1:m, 1:m, 1:m)
valid.ind = (ind[,1] != ind[,2]) & (ind[,1] != ind[,3]) & (ind[,2] != ind[,3]) & (abs(ind[,1] - ind[,2]) <= n.subdiag) & (abs(ind[,1] - ind[,3]) <= n.subdiag) & (abs(ind[,2] - ind[,3]) <= n.subdiag)
I3 = sum(valid.ind)	

##############################################
#Create variables to hold the results;
##############################################
z.mat=matrix(0,nrow=m,ncol=nsim)

#For all the following:
#v1: full sample; v2: subset

mu2.v1=rep(0, nsim); 
mu2.v2=rep(0, nsim);  
mu3.v1=rep(0, nsim); 
mu3.v2=rep(0, nsim); 

#GAWSH h^2;
h2.GWASH.v1 = rep(0, nsim)
h2.GWASH.v2 = rep(0, nsim)

#Asymptotic 95%CI and variance for GAWASH h^2
h2.var.v1= rep(0, nsim)
h2.upper.v1=rep(0, nsim)
h2.lower.v1=rep(0, nsim)


h2.var.v2= rep(0, nsim)
h2.upper.v2=rep(0, nsim)
h2.lower.v2=rep(0, nsim)
	
##############################################	
#Generate totally nsim dataset and obtain
#parameter estimates from each dataset
##############################################	
for (i in 1:nsim){ 
	
	######################
	#Generate data;
	######################
	if(X.norm){ #generate correlated normal X;
  		X = mvrnorm(n, rep(0, m), SIGMA)
  	}else{ #generate correlated Binomial X;
  		normal.cop=normalCopula(rho,dim=m,dispstr="ar1")
		u <- rCopula(n, normal.cop)
		X=array(qbinom(u,size=2,prob=p), dim=c(n, m))
  	}

	eps = rnorm(n, sd = sqrt(sigma2.eps))
	y = rep(0,n)
	y = X %*% b + eps

	X.s = scale(X) #centered & scaled
	y.s = scale(y) #centered & scaled
	
	z.mat[,i] = t(X.s) %*% y.s/sqrt(n-1); #Eq (16) in vector form
	#End of data generation;
	
	
	  
   	######################
   	#Parameter estimates
   	######################
    S=(t(X.s)%*%(X.s))/(n-1); #Eq (17), sample covariance matrix	
	######################
	#version 1: Full data
	######################		
	m.obj=mu2_est(n,m,S,1,0,0);
	mu1=m.obj$mu1 #always 1
	mu2.v1[i]=m.obj$mu2

	#GAWSH estimate h^2 under full data;
	h2.GWASH.v1[i]=h2.GWASH.est.dep(z.mat[,i],mu1,mu2.v1[i],n,m)

	mu3.v1[i]=mu3_est(n,m,S,mu1,mu2.v1[i],1,0,0)$mu3
				
	#SE estimate for GAWSH estimate h^2 under entire data
	est.obj=h2.GWASH.variance.CI.dep(n,m,h2.GWASH.v1[i],mu1,mu2.v1[i],mu3.v1[i])
	h2.var.v1[i]=est.obj$variance;
	h2.upper.v1[i]=est.obj$upper_h2;
	h2.lower.v1[i]=est.obj$lower_h2;
		

	###################
	#version 2: subset
	###################
	
	m.obj=mu2_est(n,m,S,2,array.mask,I2);
	mu2.v2[i]=m.obj$mu2
	S.mask=m.obj$S.tilde.mask		

	#GAWSH estimate h^2  under subset;
	h2.GWASH.v2[i]=h2.GWASH.est.dep(z.mat[,i],mu1,mu2.v2[i],n,m)
		
	mu3.v2[i]=mu3_est(n,m,S.mask,mu1,mu2.v2[i],2,I2,I3)$mu3

	#SE estimate for GAWSH estimate h^2 under subset
	est.obj=h2.GWASH.variance.CI.dep(n,m,h2.GWASH.v2[i],mu1,mu2.v2[i],mu3.v2[i])
	h2.var.v2[i]=est.obj$variance;
	h2.upper.v2[i]=est.obj$upper_h2;
	h2.lower.v2[i]=est.obj$lower_h2;

}

#################
#Report results;
#################

	################################
	#Based on full sample, v.1
	################################
	print(paste("h2.GWASH.v1:",mean(h2.GWASH.v1)))
	print(paste("Empirical S.E of GWASH h^2:",sd(h2.GWASH.v1)))
	
	#mu2
	print(paste("mu2.v1:",mean(mu2.v1)))
	print(paste("Empirical S.E of mu2:",sd(mu2.v1)))
	
	#mu3
	print(paste("mu3.v1:",mean(mu3.v1)))
	print(paste("Empirical S.E of mu3:",sd(mu3.v1)))
	

	#Asymptotic variance and 95% CI
	print(paste("asymptotic variance of GWASH h^2:",mean(h2.var.v1)))
	print(paste("asymptotic S.E of of GWASH h^2:",mean(sqrt(h2.var.v1/n))))
	
	print(paste("lower:",mean(h2.lower.v1[!is.na(h2.lower.v1)])))
	print(paste("upper:",mean(h2.upper.v1[!is.na(h2.upper.v1)])))


	################################	
	#Based on subset, v2
	################################
	print(paste("h2.GWASH.v2:",mean(h2.GWASH.v2)))
	print(paste("Empirical S.E of GWASH h^2:",sd(h2.GWASH.v2)))
	
	#mu2
	print(paste("mu2.v2:",mean(mu2.v2)))
	print(paste("Empirical S.E of mu2:",sd(mu2.v2)))
	
	#mu3
	print(paste("mu3.v2:",mean(mu3.v2)))
	print(paste("Empirical S.E of mu3:",sd(mu3.v2)))
	

	#Asymptotic variance and 95% CI
	print(paste("asymptotic variance of of GWASH h^2:",mean(h2.var.v2)))
	print(paste("asymptotic S.E of GWASH h^2:",mean(sqrt(h2.var.v2/n))))
	
	print(paste("lower:",mean(h2.lower.v2[!is.na(h2.lower.v2)])))
	print(paste("upper:",mean(h2.upper.v2[!is.na(h2.upper.v2)])))
