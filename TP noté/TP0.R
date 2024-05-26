# TP0
# ENSAI
# Library utiles  
library(mvtnorm) # pour simuler des vecteurs gaussiens
library(mnormt)
library(invgamma) # pour simuler une inverse gamma

####################################################### 
# Simulations des donn�es que l'on �tudiera 
nobs=200
simuvarexpl <- matrix(rep(0,300*nobs),ncol=300,nrow=nobs)
for (i in 1 :300) simuvarexpl[,i] <- runif(nobs,-5,5)
#simuvarexpl = as.data.frame(simuvarexpl)
simuvarexpl=as.matrix(simuvarexpl) 
trueind <- c(10,20,30,40,50)
beta <- c(1,-1,2,-2,3)
ysimu <- simuvarexpl[,trueind]%*% beta + rnorm(nobs,0,2)

###############################################################
# Fonction de pr�diction 
predictions <- function(matableTest,muchap,betachap)
{
ychap <- muchap * rep(1,dim(matableTest)[1])+ as.matrix(matableTest[,]) %*% betachap
return(ychap)
}

########################################################################
########################################################################
# Question 2
##############Simulation d'un coefficient aléatoire
# on a déjà calculé la loi de beta | y (c'est une gaussienne)
# on voit dans la variance a posteriori que le bay�sien fonctionne comme du ridge 
# la matrice devient inversible 
# on va la simuler (ici le Gibbs n'a qu'une seule variable multivariée ! 
# sinon il faudrait simuler composante par composante conditionnellement ce vecteur gaussien)
# Initialisaton des paramètres 
s2eps=4  # supposé connu ici 
s2beta= 1 # a priori ici ça influence le résultat 
# un s2beta trop petit va trop contraindre les beta 
# et ils auront du mal à s'éloigner de zéro
# pour aller vite on prend des tailles de burn in et d'it�rations petites
niter=1000  # nbre total itérations 
nburn=200 # nbre itérations burn in 

###############################################################
# ci dessous la fonction qui simule le beta a posteriori 
# ici on suppose connu le sigma2eps 
# le sigma2beta est consid�r� comme un hyperparam�tre
# on consid�re que mu = 0 connue 

BayesBasic  <- function(y,X,sigma2beta,sigma2eps,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
identite=diag(1,p,p)
### resultats a garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,
ncol=(nbiter-nburn))
### initialisation
beta <- rep(0,p)
### iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- solve(t(X)%*%X/sigma2eps + identite/sigma2beta)
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/sigma2eps,
Sigmabeta))
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
}
}
return(resbeta)
}
#### Illustrations 
# si on prend s2beta trop grand on a des beta trop variables
# et inversement 
s2beta=100
nburn=10
niter=210
res=BayesBasic(ysimu,simuvarexpl,s2beta,4,niter,nburn)
plot((rowMeans(res)))
plot(sort(abs(rowMeans(res))))
res=BayesBasic(ysimu,simuvarexpl,200,4,300,10)
# => trop dispers� 
res=BayesBasic(ysimu,simuvarexpl,0.1,4,300,10)
plot(sort(abs(rowMeans(res))))
which(abs(rowMeans(res))>0.8)
plot(density(res[50,]))
plot(res[50,])
mean(res[50,])
plot(density(res[51,]))
# si on prend s2beta = 0.1 alors tous les coef vont être très resserrés 
# et les vrais seront sous estim�s en valeur absolue 
# Inversement si on prend s2beta= 1000 alors les coef vont �tre plus dispers�s 
# et ils seront surestim�s en valeur absolue  
niter=500
nburn=100
# attention on devrait prendre seulement les 100 premi�res observations (training)
resbayes=BayesBasic(ysimu[1:100],simuvarexpl[1:100,],s2beta,s2eps,niter,nburn)
plot(density(resbayes[1,])) # bruit 
plot(density(resbayes[5,])) # bruit 
mean(resbayes[5,])
plot(density(resbayes[40,])) # vrai coef
mean(resbayes[40,]) # vrai coef 
plot(density(resbayes[50,])) # vrai coef
mean(resbayes[50,]) # vrai coef
#rowMeans(resBayes)
#colMeans(t(resbayes))
plot(sort(abs(colMeans(t(resbayes))))) # on voit les coef qui se d�gagent 
# on voit les 5 coef sortir du lot 

########################################################################
########################################################################
# Question 3 
# On peut utiliser le programe pr�c�dent simplement en changeant la variance
# s2beta est remplac�e par c s2eps solve(t(X)%*%X)
s2eps=4 # suppos�e connue 
c=20  # entre 10 et 100 
szellner= c*s2eps*solve(t(simuvarexpl)%*%simuvarexpl)  # variance de Zellner
# �a ne marche pas, 
# la matrice (X'X) n'est pas inversible
ridge=0.1 # param�tre ridge
p=300 
identite=diag(1,p,p)
sridge= c*s2eps*solve(t(simuvarexpl)%*%simuvarexpl + ridge*identite)
nburn=10
niter=210
res=BayesBasic(ysimu,simuvarexpl,sridge,4,niter,nburn)
plot(sort(abs(rowMeans(res))))
# prenons un param�tre ridge plus important 
ridge=10 # param�tre ridge
p=300 
identite=diag(1,p,p)
sridge= c*s2eps*solve(t(simuvarexpl)%*%simuvarexpl + ridge*identite)
nburn=10
niter=210
res=BayesBasic(ysimu,simuvarexpl,sridge,4,niter,nburn)
plot(sort(abs(rowMeans(res))))




#######################################################
#######################################################
# Question 4
# Algorithme avec  s2 al�atoire ? 
# prenons sigma2eps al�atoire de loi InvGamma(c,d) 
# On montre (voir poly) que sa loi a posteriori est une InvGamma de param�tre 
# (c+n/2,d+(Y-mu-X beta)'(Y-mu-X beta)/2)
# on programme le Gibbs ainsi : 
library(invgamma)
BayesBasic2  <- function(y,X,sigma2beta,c,d,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
identite=diag(1,p,p)
### resultats a garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn)) # param�tre vectoriel
reseps <- rep(0,(nbiter-nburn)) # param�tre univari� 
### initialisation
beta <- rep(0,p) # on aurait pu simuler un vecteur gaussien N(0,s2beta I)
# beta <- rmnorm(1,rep(0,p), sigma2beta * identite)
sigma2eps <- var(y)   # on intialise avec l'empirique
# on aurait pu aussi simuler une inv gamma (c,d)
# sigma2eps <- rinvgamma(1,c,d)
### iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- solve(t(X)%*%X/as.numeric(sigma2eps) + identite/sigma2beta)
# premi�re marge conditionnelle sachant s2eps
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/as.numeric(sigma2eps),Sigmabeta))
s=t(y-X%*%beta)%*%(y-X%*%beta)
# deuxi�me loi conditionnelle sachant beta
sigma2eps <- rinvgamma(1,c+n/2,d+s/2)
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
reseps[iter-nburn] <- sigma2eps
}
}
resu=list(resbeta,reseps)
return(resu)
}
############ Illustrations 
# on fixe c=1=d pour l'inverse Gamma => vague car esp�rance et variance = infini
c=d=1  # => la variance du prior vaut l'infini !
# moyenne inverse gamma(c,d) = c/(d-1)
# variance inverse gamma(c,d) = c^2/((d-1)^2(d-2)) 
sigma2beta=10 # a priori assez vague sur la variance de beta
res2=BayesBasic2(ysimu,simuvarexpl,sigma2beta,c,d,600,300)
resbeta=res2[[1]]
reseps=res2[[2]] 
plot(resbeta[1,])
plot(resbeta[50,])
plot(reseps)
plot(density(resbeta[1,]))
plot(density(resbeta[50,]))
plot(density(resbeta[40,]))
plot(density(reseps))
mean(resbeta[50,]) # moyenne a posteriori qui est l'estimateur bay�sien
mean(resbeta[40,])
mean(reseps) # moyenne a posteriori 
plot(sort(abs(colMeans(t(resbeta)))))
plot(reseps)
plot(resbeta[50,])
plot(resbeta[110,])
# on peut regarder la distribution par exemple de beta50 
# il faudrait extraire une valeur sur 10 (par ex) pour rendre les donn�es 
# moins d�pendantes 
plot(density(resbeta[50,]))
########################################################################
########################################################################
# Question 5
# Ici on est tr�s proche du Bayes A qui sera �tudier au TP1 
# ce qui change ici c'est que mu est suppos�e �gale � 0
# Mais on va le faire 
BayesBasic3  <- function(y,X,a,b,c,d,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
identite=diag(1,p,p)
### resultats a garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn)) # param�tre vectoriel
reseps <- rep(0,(nbiter-nburn)) # param�tre univari� 
resvarbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn))
### initialisation
beta <- rep(0,p) # car les beta centr�es 
sigma2beta = rinvgamma(1,a,b)
sigma2eps <- var(y)   # on intialise avec l'empirique
### iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- solve(t(X)%*%X/as.numeric(sigma2eps) + identite/sigma2beta)
# premi�re marge conditionnelle sachant s2eps (et Y)
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/as.numeric(sigma2eps),Sigmabeta))
s=t(y-X%*%beta)%*%(y-X%*%beta)
# deuxi�me loi conditionnelle sachant beta (et Y)
sigma2eps <- rinvgamma(1,c+n/2,d+s/2)
# troisi�me loi conditionnelle sachant beta (et Y)
sigma2beta  = rinvgamma(p,a+1/2,b+1/2*beta^2)
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
reseps[iter-nburn] <- sigma2eps
resvarbeta[,iter-nburn] <- sigma2beta
}
}
resu=list(resbeta,reseps,resvarbeta)
return(resu)
}
######## on va essayer avec bayesAA
a=c=1
b=d=3
# a=b=c=d=1  # tr�s vague car exp�rance et variance infinies 
res3=BayesBasic3(ysimu,simuvarexpl,a,b,c,d,800,300)
res4=BayesAA(ysimu,simuvarexpl,a,b,c,d,800,300)
res3=res4
resbeta=res3[[1]]
reseps=res3[[2]]
resvarbeta=res3[[3]]
resmu=res3[[4]]
par(mfrow=c(2,2))
plot(resbeta[1,])
plot(resbeta[10,])
plot(resbeta[40,])
plot(resbeta[50,])
plot(reseps)
plot(resmu)
plot(resvarbeta[2,])
plot(resvarbeta[40,])
plot(resvarbeta[200,])
plot(density(resbeta[1,]))
plot(density(resbeta[50,]))
plot(density(resbeta[40,]))
plot(density(reseps))
plot(density(resvarbeta[1,]))
plot(density(resvarbeta[50,]))
plot(density(resvarbeta[40,]))
mean(resbeta[50,]) # moyenne a posteriori qui est l'estimateur bay�sien
mean(resbeta[40,])
mean(reseps) # moyenne a posteriori 
resbeta=resbeta[,300:799]
plot(sort(abs(colMeans(t(resbeta)))))
# on rep�re les 5 coefficients les plus significatifs du mod�le 
# on peut regarder leurs trajectoires, leurs densit�s estim�es 
par(mfrow=c(1,1))
plot(resbeta[10,])
plot(resbeta[20,])
plot(resbeta[30,])
plot(resbeta[40,])
plot(resbeta[50,])

plot(density(resbeta[10,]))
plot(density(resbeta[20,]))
plot(density(resbeta[30,]))
plot(density(resbeta[40,]))
plot(density(resbeta[50,]))
# Critiques : il faudrait prendre 1) un burn-in plus grand (10.000) 
# 2) prendre une simu sur 10 (pour diminuer la d�pdendance) et avoir 
# ainsi une �chantillon "presque" "i"id pour construire par exemepl des intervalles
# de confiance a posteriori => c'est un moyen de s�lectionner les betas 
# significatifs (ceux dont 0 n'est pas dans l'intervalle). 
# Le mod�le repose sur la normalit� des Y|beta,sigma 
# mais c'est assez robuste � la non-normalit�  
 # une analyse de sensibilit� serait de faire varier les hyperparam�tres : 
# a,b,c,d (les prendre plus grands...) 

########################################################################
# Question 5 bis
# Ici on fait tout le Bayes A  
BayesAA  <- function(y,X,a,b,c,d,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
identite=diag(1,p,p)
vec1=rep(1,n)
### resultats a garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn)) # param�tre vectoriel
reseps <- rep(0,(nbiter-nburn)) # param�tre univari� 
resmu <- rep(0,(nbiter-nburn))  ##############
resvarbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn))
### initialisation
beta <- rep(0,p) # car les beta centr�es 
sigma2beta = rinvgamma(1,a,b)
sigma2eps <- var(y)   # on intialise avec l'empirique
mu <- 0 
### iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- solve(t(X)%*%X/as.numeric(sigma2eps) + identite/sigma2beta)
# premi�re marge conditionnelle sachant s2eps (et Y)
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%(y-mu*vec1)/as.numeric(sigma2eps),Sigmabeta)) ######
s=t(y-mu*vec1-X%*%beta)%*%(y-mu*vec1-X%*%beta)
# deuxi�me loi conditionnelle sachant beta (et Y)
sigma2eps <- rinvgamma(1,c+n/2,d+s/2)
# troisi�me loi conditionnelle sachant beta (et Y)
sigma2beta  = rinvgamma(p,a+1/2,b+1/2*beta^2)
# loi de mu 
mu = rnorm(1,vec1%*%(y-X%*%beta)/n,sigma2eps/n)
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
reseps[iter-nburn] <- sigma2eps
resvarbeta[,iter-nburn] <- sigma2beta
resmu[iter-nburn] <- mu
}
}
resu=list(resbeta,reseps,resvarbeta,resmu)
return(resu)
}


########################################################################
# Question 6
########################### # R�gression bay�sienne et ABC 
# Illustration avec un mod�le de scale mixture
# Simulation du mod�le
# on pose 
# Y_i = mu_i + X_i beta_i + e_i 
# avec e ~ N(0,s)
# mu ~uniforme (-a,a)
# beta_i ~ N(m,v) 
#  
#####################################################
# Les param�tres 
n=200 # nbre observations (simulations)
d=1 # taille beta 
m=6
v=1
s=1
m0=0
a=5
b=5
x <- matrix(rep(0,d*n),ncol=d,nrow=n)
for (i in 1:d){
x[,i] <- runif(n,-5,5)}
x <- as.data.frame(x)
X=as.matrix(x) 
beta <- rnorm(n,m,v) # chaque beta est al�atoire 
m0 <- runif(n,-a,a)
s <- rchisq(n,b)
y <- m0 + X%*%beta + rnorm(n,0,s)  # nos y simul�s 



############################################################################################
# ABC basique 
############################################################################################
ABCbasic  = function(y, X, a, b, V, seuil, K)
{ 
d = dim(X)[2]
identite=diag(rep(1,d)) # matrice identit� =1 (dimension =1)
identn=diag(rep(1,n))    # matrice identit� n x n 
initbeta=rep(0,(d*K))
resbeta=matrix(initbeta,ncol=K)
resmu=rep(0,K)
ress2=rep(0,K)
n=length(y)
compteur=0
for (j in 1:(K-1))
{
mu=runif(1,-a,a) # a priori uniforme 
s2=rchisq(1,b) # a priori sur la variance du bruit 
s=sqrt(s2)
beta=rmnorm(1,0,V*identite ) # a priori sur beta 
m=mu+X%*%beta  # moyenne des y sachant les param�tres
res = rmnorm(1,m,s*identn) # simulation du y sachant les param�tres 
if (sqrt(sum((y-res)^2)/n) < seuil)  # acceptation rejet 
{
compteur = compteur+1  
resbeta[,compteur] =beta 
ress2[compteur] =s2 
resmu[compteur] = mu 
}
}
accept=sum(resmu!=0) # nbre d'acceptation 
betafinal=resbeta[,c(1:accept)]
resmufinal=resmu[c(1:accept)]
ress2final=ress2[c(1:accept)]
return(list(mu=resmufinal,beta=betafinal,s2=ress2final, taux=100*accept/(K-1)))
}
##################################################
a=5
b=5
V=40
seuil=8
K=20000
res = ABCbasic(y, X, a, b, V, seuil, K)
res$taux
mean(res$mu)
mean(res$s2)
mean(res$beta)
plot(density(res$beta))
# le d�s�quilibre peut venir du fait que l'on pr�sente plus souvent des beta 
# � gauche qu'� droite de 5 
# il s'estompe si on diminue le seuil ou si on augmente V
# il y a un pb d'identifiabilit� avec mu !!! 
# il faut enlever mu ou alors reporter la moyenne de beta sur mu !!! 




# Question 3 (suite / compl�ment)
#### Algorithme avec variance de Zellner
# Pour extraire une  matrice inversible on peut prendre les premi�res variables
# on peut aussi se r�f�rer � l'indice de conditionnement
# Par exemple calculer les valeurs propres
X=as.matrix(simuvarexpl[1:100,]) 
valpX=eigen(t(X)%*%X, symmetric =TRUE, only.values = TRUE)
valpX
valmax=max(valpX$values)
valmin=min(valpX$values)
indcond=valmax/valmin
indcond
# on trouve ici quelque chose de n�gatif car ce sont les approximations 
# (normalement �a vaut z�ro) 
# Prenons une sous matrice de taille 50*50 par exemple 
indalea=sample(c(1:300),50,replace=FALSE)
X2=X[,indalea] 


# c = coefficiet d'�chelle 
# L = coefficient ridge
BayesZellner  <- function(y,X,c ,L, sigma2beta,sigma2eps,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
identite=diag(1,p,p)
### resultats a garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,
ncol=(nbiter-nburn))
### initialisation
beta <- rep(0,p)
### iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- c*solve(t(X)%*%X/sigma2eps + L*identite/sigma2beta)
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/sigma2eps,
Sigmabeta))
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
}
}
return(resbeta)
}
# Algorithme avec Zelnner (+ Ridge ? => question suppl�mentaire) 
s2eps=4
s2beta=200
# si on prend s2beta = 0.1 alors tous les coef vont �tre tr�s resserr�s 
# et les vrais seront sous estim�s en valeur absolue 
# Inversement si on prend s2beta= 1000 alors les coef vont �tre plus dispers�s 
# et ils seront surestim�s en valeur absolue  
niter=300
nburn=100
c=1
L=100
# attention ici on prend seulement les 100 premi�res observations (training)
resbayes=BayesZellner(ysimu[1:100],as.matrix(simuvarexpl[1:100,]),c,L, s2beta,s2eps,niter,nburn)
######################################################################
plot(density(resbayes[1,])) # bruit 
plot(density(resbayes[5,])) # bruit 
mean(resbayes[5,])
plot(density(resbayes[40,])) # vrai coef
mean(resbayes[40,]) # vrai coef 
plot(density(resbayes[50,])) # vrai coef
mean(resbayes[50,]) # vrai coef
#colMeans(t(resbayes))
plot(sort(abs(colMeans(t(resbayes))))) # on voit les coef qui se d�gagent 
# Rermarque : on pourrait ajouter un coefficient ridge pour stabiliser l'inverse.
