###############################################################
# Fonction de prédiction 
# Elle servra a prédire des valeurs de Y 
# à partir des beta estimés et des x du jeu de données test 
predictions <- function(matableTest,muchap,betachap)
{
  ychap <- muchap * rep(1,dim(matableTest)[1])+ as.matrix(matableTest[,]) %*% betachap
  return(ychap)
}

############################
# Programme pour sélectionner les coefficients 
# les plus significatifs, càd ceux qui sortent d'un intervalle (mini, maxi) 
# et ce programme renvoie les noms des variables associ?es 
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

############################
# Programme pour la mise en place du modèle Bayes A
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


##### SSVS
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

##### Amélioration ssvss

selection_SSVS_2 <- function (vardep, varexpl, nbiter, nburn, lec, nbSelecInit, Pi) {
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
  nbToChange <- round(sqrt(nvar))  # Modification : Utilisation d'un nombre adaptatif de variables à changer
  
  for (iter in 1:nbiter) {
    print(iter)
    gammaprop <- gamma
    indgamma1prop <- indgamma1
    indToChange <- sample(c(1:nvar), nbToChange, replace = FALSE)
    for (i in 1:nbToChange){
      if (gamma[indToChange[i]]==0){
        gammaprop[indToChange[i]] <- 1
        indgamma1prop <- c(indgamma1prop,indToChange[i])
      } else {
        gammaprop[indToChange[i]] <- 0
        indremove <- which(indgamma1prop==indToChange[i])
        indgamma1prop <- indgamma1prop[-indremove]
      }
    }
    nbSelecprop <- length(indgamma1prop)
    if (nbSelecprop==0){
      cond <- 0
      while(cond==0){
        gammaprop <- gamma
        indgamma1prop <- indgamma1
        indToChange <- sample(c(1:nvar), nbToChange, replace = FALSE)
        for (i in 1:nbToChange){
          if (gamma[indToChange[i]]==0){
            gammaprop[indToChange[i]] <- 1
            indgamma1prop <- c(indgamma1prop,indToChange[i])
          } else {
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

