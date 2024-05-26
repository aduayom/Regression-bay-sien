# TP1
# ENSAI
library(mvtnorm)
library(mnormt)
# Remarque : il faut normaliser le design X en faisant X <- scale(X)
# Ici les X sont toutes de m�me loi 
# Simulations des donn�es 
simuvarexpl <- matrix(rep(0,300*200),ncol=300,nrow=200)
for (i in 1 :300) simuvarexpl[,i] <- runif(200,-5,5)
# on les normalise
simuvarexpl=scale(simuvarexpl)
simuvarexpl = as.data.frame(simuvarexpl)
# vrais indices de la r�gression 
trueind <- c(10,20,30,40,50)
# les valeurs non nulles des beta sont 
beta <- c(1,-1,2,-2,3)
# simulation des Y
ysimu <- as.matrix(simuvarexpl)[,trueind]%*% beta + rnorm(200,0,2)
###############################################################
# Fonction de pr�diction 
# Elle servra � pr�dire des valeurs de Y 
# � partir des beta estim�s et des x du jeu de donn�es test 
predictions <- function(matableTest,muchap,betachap)
{
ychap <- muchap * rep(1,dim(matableTest)[1])+ as.matrix(matableTest[,]) %*% betachap
return(ychap)
}
############################
# Programme pour s�lectionner les coefficients 
# les plus significatifs, c�d ceux qui sortent d'un intervalle (mini, maxi) 
# et ce programme renvoie les noms des variables associ�es 
subsetSelected <- function(resAlgo, varexpl, mini, maxi) 
{
numselected <- c(which(resAlgo < mini), which(resAlgo > maxi))
selected <- character()
valeurs <- numeric()
for (i in 1 :length(numselected)) 
{
selected[i] <- names(varexpl)[numselected[i]]
valeurs[i] <- resAlgo[numselected[i]] 
}
subset <- cbind(selected, valeurs)
subset <- as.data.frame(subset)
subset$valeurs <- as.numeric(as.vector(subset$valeurs))
return(subset) 
}

# Remarque : il faudrait faire de la validation crois�e en 
# changeant le jeu de donn�es TEST 

#####################Mod�le RR 
# Y = 1*mu + X beta + eps # notation du cours 
# avec les notations du package mu = beta 
# beta = u 
# X = Z 
# le 1 = X 
# Y = X*beta + Z u + e  
library(rrBLUP)

######## le mod�le sur le training (les 100 premi�res observations)
resBLUP <- mixed.solve(ysimu[1 :100], Z = as.matrix(simuvarexpl[1 :100, ]), X = as.matrix(rep(1,
100)), method = "REML")
# les notations de mixed.solve sont : 
# notre mu = beta 
# notre beta = u 
# notre epsilon = e
######## les estimations 
muchap <- resBLUP$beta  # moyenne  mu
muchap
betachap <- resBLUP$u   # les coefficients beta
# c'est un vecteur de longueur p=300
par(mfrow=c(1,1))
plot(sort(abs(betachap)))   # on voit un fort shrinkage 
betachap[50]
resBLUP$Ve  # variance des epsilon
resBLUP$Vu # variance des betas
# on a des variances tr�s instables  
# des beta sous-etim�s !!! 

################### Les pr�dictions sur les donn�es du jeu test (les 100 derni�res valeurs)
predBLUP <- predictions(simuvarexpl[101 :200, ], muchap, betachap)
# on r�cup�re 100 pr�dictions de Y 

################ corr�lation en y pr�dits et y observ�s (sur jeu test)
cor(ysimu[101 :200], predBLUP)
# 0.4461686
plot(predBLUP ~ ysimu[101 :200])
abline(0,1)
# on s'aper�oit de l'effet shrinkage 

############ S�lection des coefficients 
plot(betachap,xlab="Indice variable",ylab="Param�tre de r�gression associ�",main="rrBLUP")
plot(sort(abs(betachap)))
boxplot(betachap)
bb <- boxplot(betachap)
varselectedBLUP <- subsetSelected(betachap,simuvarexpl,bb$stats[1,], bb$stats[5,])
varselectedBLUP
#On retient les valeurs suivantes 


1      V40 -0.4682636
2      V10  0.3532256
3      V30  0.4055062
4      V50  0.6666112
5     V112  0.3073959


################ Bayes A
library(MCMCpack) ; library(LearnBayes) ; 
library(invgamma)
################## le mod�le
# y|beta, sigmaSeps, mu ~N(mu + X beta , sigma2eps)
# u uniforme 
# beta ~N(0, sigma2(1), ..., sigma2(p)) 
# sigma2(i) ~ InvGamma(a,b) iid
# sigma2eps ~InvGamma(c,d)
BayesA <- function(y,X,a,b,c,d,muinit,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
# resultats a garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn))
ressigma2beta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn))
resmu <- rep(0,nbiter-nburn)
ressigma2eps <- rep(0,nbiter-nburn)
# initialisation
beta <- rep(0,p)
mu <- muinit
sigma2beta <- rinvgamma(p,a,b) # initialisation des valeurs 
sigma2eps <- rinvgamma(1,c,d)
#iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- solve(t(X)%*%X/sigma2eps + diag(1/sigma2beta))
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%(y-mu*rep(1,n))/sigma2eps, Sigmabeta))
mu <- rnorm(1,t(rep(1,n))%*%(y-X%*%beta)/n,sqrt(sigma2eps/n))
for (j in 1 :p)
{
sigma2beta[j] <- rinvgamma(1,a+1/2,b+1/2*beta[j]^2)
}
sigma2eps <- rinvgamma(1,c+n/2,d+1/2*t(y-mu*rep(1,n)-X%*%beta)%*%(y-mu*rep(1,n)-X%*%beta))
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
ressigma2beta[,iter-nburn] <- sigma2beta
resmu[iter-nburn] <- mu
ressigma2eps[iter-nburn] <- sigma2eps 
}
}
return(list(resbeta,ressigma2beta,resmu,ressigma2eps)) 
}

####################### Application sur nos donn�es
# les a priori sont des inverses gamma 
# de param�tres (a,b) et (c,d) 
# leurs moyennes et variances sont 
a=c=1 ; b=d=1
b/(a-1) ; d/(c-1)   # esp�rance d'une inv gamma
b^2/((a-1)^2*(a-2)) ; d^2/((c-1)^2*(c-2)) # variance inv gamma
# avec ces valeurs ont met des a priori tr�s vague, de variance infinie 
priora <- 2
priorb <- 1
priorc <- 2
priord <- 1
resBAYESA <- BayesA(ysimu[1 :100], as.matrix(simuvarexpl[1 :100, ]), priora, priorb, priorc,
priord, mean(ysimu[1 :100]), 2000, 1500)
moybeta <- apply(resBAYESA[[1]], 1, mean) # estimateurs des beta
moysigma2beta <- apply(resBAYESA[[2]], 1, mean)
moymu <- mean(resBAYESA[[3]])
moysigma2eps <- mean(resBAYESA[[4]])
moybeta[50]
plot(moybeta)
plot(sort(abs(moybeta)))
moymu
plot(moysigma2beta)
moysigma2eps
# regardons les trajectoires : 
par(mfrow=c(2,2))
plot(resBAYESA[[1]][50, ])  # pour beta50
plot(resBAYESA[[1]][10, ])  # pour beta10
plot(resBAYESA[[1]][20, ])  # pour beta20
plot(resBAYESA[[1]][1, ])  # pour beta1
plot(resBAYESA[[2]][1, ])   # pour s2beta1
plot(resBAYESA[[3]])        # pour mu
plot(resBAYESA[[4]])        # pour s2eps 
# en conclusion il faut un plus grand burn-in. 50.000 par exemple. 
# ou alors lancer plusieurs cha�nes et comparer les trajectoires (empirique)
# dans tous les cas il peut y avoir avoir une d�pendance des variables g�n�r�es 
# par l'algo. 

# on peut relancer avec des a priori plus ou moins vagues 
# par exemple en prenant 
priora <- 10   # variance finie    
priorb <- 10   # encore plus variable
priorc <- 10   # variance finie 
priord <- 10   # plus variable 
# alors les variances sont petites ou grandes 
############# Comparaison des lois a priori et a posteriori
tracesigmabeta1 <- resBAYESA[[2]][30,]
plot(density(tracesigmabeta1),col="blue",lwd=2, main="sigma_beta_1 : Bayes A")
# attention il faudrait prendre "une variable sur 10" car ce n'est pas iid (c'est juste id)
curve(dinvgamma(x,priora,priorb), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))

tracesigmabeta50 <- resBAYESA[[2]][50,]
plot(density(tracesigmabeta50),col="blue",lwd=2, main="sigma_beta_50 : Bayes A", ylim=c(0,2))
curve(dinvgamma(x,2,1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))

tracesigmaeps <- resBAYESA[[4]]
plot(density(tracesigmaeps),col="blue",lwd=2, main="sigma_eps : Bayes A")
curve(dinvgamma(x,priorc,priord), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))
# WARNING !! il faudrait utiliser un lag pour estimer les distributions car 
# les données sont dépendantes ce qui peut fausser l'estimation 

# Int�r�t en bay�sien s�quentiel : si on fait le m�me mod�le avec des donn�es diff�rentes 
# on peut utiliser les lois a posteriori comme a priori ensuite... 
# => �a rend l'algorithme compl�tement automatique 
# => pas d'hyperparam�tre � r�gler 

# On pourrait faire la m�me chose si l'�chantillon est assez grand 
# on en prend 10% par exemple pour r�cup�rer des lois a posteriori 
# Et ensuite on travaille avec ces lois comme a posteriori sur les 90% des donn�es restantes 

##################### Pr�dictions 
predBAYESA <- predictions(simuvarexpl[101 :200, ], moymu, moybeta)  # sur le jeu test
cor(ysimu[101 :200], predBAYESA)
# 0.4772729
plot(predBAYESA ~ ysimu[101 :200])  # pr�dictions versus observations sur jeu test 
abline(0,1)

#################### S�lection des betas
plot(moybeta,xlab="Indice variable",ylab="Param�tre de r�gression associ�",main="Bayes
A")
plot(sort(abs(moybeta)))
bb <- boxplot(moybeta)
varselectedBAYESA <- subsetSelected(moybeta,simuvarexpl, bb$stats[1,],bb$stats[5,])
varselectedBAYESA

selected    valeurs
1      V40 -0.8353826
2     V228 -0.5354466
3      V10  0.4885054
4      V50  1.1065191
5     V197  0.4886490

# on a autant de var mais avec une diff�rence 
# et moins de "shrinkage" 

####################### LASSO Bay�sien 
# On utilise le package BLR(Y,XF,XR,XL,Z,prior,nIter, burnIn
# car Y = mu*1 + XF betaF + XR betaR + XL betaL + Z u + epsilon 
# avec betaF = fixe (= mod�le classique, ici non) 
# mu = constante 
# betaR = ridge (= bayes A, avec plusieurs variance, ici non)  
# betaL = Lasso (ici oui) 
# u = effet al�atoire du RR (avec une seule variance, ici non)

############# le modèle 
 library(BLR)
LASSO_BLR <- BLR(y=ysimu[1 :100],XL=as.matrix(simuvarexpl[1 :100,]), prior=list(varE=list(df=2,S=1),
lambda=list(shape=10, rate=0.1, type='random',value=2)), nIter=3000,burnIn=2000, saveAt="BLR_",
thin=1)#,thin2=1)
# thin = 1 signifie qu'il prend toutes les simus (1/10 par d�faut)
# thin = 5 => il garde une simu sur 5 
# attention il garde même le burn in 
# 
# les résultats sont sauvés avec saveAt="BLR_" 
# Attention : si on relance le programme les résultats sont concaténés (soit changer de nom, 
# soit le supprimer)

# remarque : la variance de lambda^2 ici vaut shape/rate^2 
# plus on prend rate petit et plus la variance va �tre �lev�e et lambda pourra �voluer librement
# VarE est la variance de epsilon qui suit une inverse khi deux 
# lambda de type 'random', value=2, signifie que lambda^2 suit une gamma 
# on peut m�me fixer la valeur de lambda, par exemple issue d'un LASSO classique
LASSO_BLR$mu # muchap = moyenne des trajectoires pour mu = estimation de la moyenne a posteriori
LASSO_BLR$varE # = moyenne des trajectoires pour sigma2epsilon= estimation de la moyenne a posteriori
LASSO_BLR$bL   # betachap = moyenne des trajectoires pour les betas = estimation de la moyenne a posteriori
LASSO_BLR$tau2  # = moyenne des trajectoires pour les tau = estimation de la moyenne a posteriori
LASSO_BLR$lambda # = moyenne des trajectoires pour lambda^2 = estimation de la moyenne a posteriori

# remarque : il serait int�ressant de faire varier le shape et le rate de la loi (a priori) gamma du 
# param�tre LASSO lambda. On regarderait alors si la corr�lation est meilleure

################## Trajectoires simul�es 
# Attention : 
# si on relance alors les trajecgoires se cumulent (concat�nation des r�sultats)
# il faut soit effac� BLR_, soit changer le nom...
# Par ailleurs, on ne garde qu'une valeur sur 10 (par d�faut, pour garantir une "ind�pendance") 

tracevarE <- scan('BLR_varE.dat')
plot(tracevarE,type='o',xlab="iteration",ylab="varE",main="trace varE")
# on observe une rupture => � explorer sur une cha�ne plus longue
# ou bien sur plusieurs cha�nes 
plot(density(tracevarE),col="blue",lwd=2, main="varE : LASSO bay�sien")
curve(dinvgamma(x,1,1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))


tracelambda <- scan('BLR_lambda.dat')
plot(tracelambda,type='o',xlab="iteration",ylab="lambda",main="trace lambda")
plot(density(tracelambda),xlim=c(0,15),col="blue",lwd=2, main="lambda : LASSO bay�sien")
curve(dgamma(x,10,1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))
# il faudrait peut être choisir une autre loi pour lambda 
# la loi a posteriori est très différente de l'a priori


#####################Pr�diction 
predLASSOBLR <- predictions(simuvarexpl[101 :200, ], LASSO_BLR$mu, LASSO_BLR$bL)
cor(ysimu[101 :200], predLASSOBLR)
#  0.7771987
plot(predLASSOBLR ~ ysimu[101 :200])
abline(0,1)

# La corr�lation est bien meilleure que pour les deux pr�c�dentes m�thodes. 
# Pour l'instant LASSO > BAYES A > Random Regression (en termes de corr�lation sur le TEST)

###################### S�lection
plot(LASSO_BLR$bL,xlab="Indice variable",ylab="Param�tre de r�gression associ�", main="LASSO
bay�sien")
plot(sort(abs(LASSO_BLR$bL)))
bb <- boxplot(LASSO_BLR$bL)
bb
varselectedLASSO <- subsetSelected(LASSO_BLR$bL,simuvarexpl, bb$stats[1,],bb$stats[5,])
varselectedLASSO
# on va faire du seuillage � la main pour �liminer quelques variables 
varselectedLASSO <- subsetSelected(LASSO_BLR$bL,simuvarexpl,-0.4,0.4)
varselectedLASSO

 selected    valeurs
1      V20 -0.5266782
2      V40 -1.2790562
3     V251 -0.4326029
4      V10  0.4750113
5      V30  0.8649100
6      V50  2.3414589
7     V112  0.4152102

# Si on sélectionnait les 5 meilleures ce serait les bonnes !! 

which(abs(LASSO_BLR$bL)>0.26)
selec5=order(abs(LASSO_BLR$bL),decreasing=TRUE)[1:5]
selec5
# Les bonnes variables sont retenues en premier... 

# Pour s�lectionner moins de variable il faudrait un srinkage plus fort, 
# ce qui correspond � un param�tre lambda du LASSO plus �lev�
# Ici lambda^2 suit une Gamma(e,f) avdc (e,f) = (1,0.1)
# E(lambda^2)=e/f  et V(lamdba^2) = e^2/f
# Prenons e + grand (=10) 

################# SSVS 
####### Mod�lisation 
 selection_SSVS <- function (vardep, varexpl, nbiter, nburn, lec, nbSelecInit,nbToChange,Pi)
 {
     X <- as.matrix(varexpl)
     y <- as.numeric(as.vector(vardep))
     y <- y - mean(y)
     nind <- dim(X)[1]
     nvar <- dim(X)[2]
     sumgamma <- rep(0, nvar)
     nbactu <- 0
     nbselecactu <- numeric()
     indgamma10 <- sample(c(1:nvar), nbSelecInit, replace = FALSE)
     gamma0 <- rep(0, nvar)
     for (i in 1:nbSelecInit) {
         gamma0[indgamma10[i]] <- 1
     }
     indgamma1 <- indgamma10
     gamma <- gamma0
     nbSelec <- nbSelecInit
     
     for (iter in 1:nbiter) {
         print(iter)
         gammaprop <- gamma
         indgamma1prop <- indgamma1
         indToChange <- sample(c(1:nvar), nbToChange, replace = FALSE)
         for (i in 1:nbToChange){
           if (gamma[indToChange[i]]==0){
             gammaprop[indToChange[i]] <- 1
             indgamma1prop <- c(indgamma1prop,indToChange[i])
           }
           else {
             gammaprop[indToChange[i]] <- 0
             indremove <- which(indgamma1prop==indToChange[i])
             indgamma1prop <- indgamma1prop[-indremove]
           }
         }
         nbSelecprop <- length(indgamma1prop)
         if (nbSelecprop==0){ # condition pour empecher gamma avec que des 0
           cond <- 0
           while(cond==0){
             gammaprop <- gamma
             indgamma1prop <- indgamma1
             indToChange <- sample(c(1:nvar), nbToChange, replace = FALSE)
             for (i in 1:nbToChange){
               if (gamma[indToChange[i]]==0){
                 gammaprop[indToChange[i]] <- 1
                 indgamma1prop <- c(indgamma1prop,indToChange[i])
               }
               else {
                 gammaprop[indToChange[i]] <- 0
                 indremove <- which(indgamma1prop==indToChange[i])
                 indgamma1prop <- indgamma1prop[-indremove]
               }
             }
             nbSelecprop <- length(indgamma1prop)
             if (nbSelecprop>0){cond <- 1}
           }
         }
         indgamma1 <- which(gamma == 1)
         nbSelec <- length(indgamma1)
         Xgamma <- X[, indgamma1]
         Xgammaprop <- X[, indgamma1prop]
         temp <- (t(y)%*%(diag(rep(1,nind))-lec/(1+lec)*Xgammaprop%*%solve(t(Xgammaprop)%*%Xgammaprop)%*%t(Xgammaprop))%*%y)/(t(y)%*%(diag(rep(1,nind))-lec/(1+lec)*Xgamma%*%solve(t(Xgamma)%*%Xgamma)%*%t(Xgamma))%*%y)
         A <- (1+lec)^((nbSelec-nbSelecprop)/2)*(Pi/(1-Pi))^(nbSelecprop-nbSelec)*temp^(-(nind-1)/2)
         probaccept1 <- min(1,A)
         seuil <- runif(1)
         if (seuil < probaccept1){
           gamma <- gammaprop
           indgamma1 <- indgamma1prop
           nbSelec <- nbSelecprop
           nbactu <- nbactu+1
           nbselecactu <- c(nbselecactu,nbSelec)
         }
         if (iter > nburn) {
             sumgamma <- sumgamma + gamma
           }
     }
     return(list(sumgamma,nbactu,nbselecactu))
 }
#### Initialisation 
nbinit <- 10 # nombre de beta non nuls au d�part 
nbToChange <- 2 # nombre de beta que l'on propose de changer � chaque fois 
Pi <- 10/300   # la probabilit� de choisir un beta (en moyenne 10 car il y en a 300)
lec <- 50  # pr�conisation : entre 10 et 100 
resSSVS <- selection_SSVS(ysimu[1 :100], as.matrix(simuvarexpl[1 :100, ]), 8000, 7000, lec,
nbinit, nbToChange, Pi)
resSSVS[[2]]  # on voit combien il y a eu d'�changes dans les 1000 derniers run
resSSVS[[3]]  # on voit le nbre de variables s�lectionn�es lors de ces �changes 
resSSVS[[1]]
####### S�lection 
plot(sort(abs(resSSVS[[1]])))
boxplot(resSSVS[[1]])
plot(resSSVS[[1]],xlab="Indice variable",ylab= "Nombre de s�lections post-burn-in", main="SSVS")
varselectedSSVS <- subsetSelected(resSSVS[[1]],simuvarexpl,0,999)
varselectedSSVS
which(resSSVS[[1]]>600)
# 10  20  30  40  50  54 163 185 191 253
# On a s�lectionn� les 5 bonnes variables. 
# Les variables suppl�mentairtes ne sont pas celles du LASSO, ni Random Regression, ni BAYES A. 
# Or on sait que le LASSO a tendance � ne pas s�lectionner des variables 
# corr�l�es. C'est peut �tre une explication... 

# on peut changer pi 
# en fait le burn-in permet de chercher les bons betas. Une fois qu'ils 
# sont s�l�ctionn�s ils restent car les autres beta ont une densit�s conditionnelle a 
# posteriori moins forte.  
# CONSEIL : ne pas hésiter à prendre un long burn in !! 

### CONCLUSION 


##################### 
# Mod�le lin�aire simple avec les 10 ou 14 variables retenues LASSO bay�sien et par SSVS 
indLASSO=c() 
indSSVS=c() 
resLASSO=lm(ysimu~as.matrix(simuvarexpl[,indLASSO])) 
resSSVS=lm(ysimu~as.matrix(simuvarexpl[,indSSVS])) 
summary(resLASSO)
summary(resSSVS)
# on peut comparer les pr�dictions sur le TEST et dons les r�sidus des deux mod�les 
# et choisir le mod�le qui a le plus petit r�sidu.
# Remarque : il faudrait faire de la validation crois�e en changeant le jeu de donn�es TEST 

# Une probable explication de la diff�rence entre la s�lection LASSO et SSVS est que le LASSO 
# a tendance � moins retenir les groupes de covariables tr�s corr�l�es





