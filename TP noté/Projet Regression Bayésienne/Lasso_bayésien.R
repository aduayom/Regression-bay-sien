library(BLR)

# Lasso Bayésien vs Approche Bayes A
# ----------------------------------

## Lasso Bayésien :
# - Utilise une régularisation de type Lasso pour la sélection et la régularisation des variables.
# - Chaque coefficient βj suit une distribution a priori Laplace doublement exponentielle.
# - La régularisation favorise des solutions éparses où certains coefficients peuvent être exactement égaux à zéro.
# - Les paramètres de la distribution a priori peuvent être ajustés pour contrôler la force de la régularisation.
# - Estimation des paramètres via des algorithmes MCMC ou des méthodes variationnelles.

## Approche Bayes A :
# - Chaque coefficient βj suit une loi normale multivariée avec sa propre variance.
# - Les variances des coefficients et de l'erreur suivent des lois inverse Gamma.
# - L'a priori sur µ est vague.
# - Estimation des paramètres via un Gibbs sampler.

## Avantages du Lasso Bayésien :
# - Sélection automatique des variables grâce à la régularisation Lasso.
# - Facilite l'interprétation des résultats en produisant des modèles plus simples et plus épars.
# - Permet un contrôle fin de la régularisation et de la sélection des variables par l'ajustement des hyperparamètres.

## Avantages de l'Approche Bayes A :
# - Plus de flexibilité en permettant à chaque coefficient d'avoir sa propre variance.
# - Meilleure adaptation aux données avec des variations hétérogènes entre les coefficients.
# - Approche bayésienne permettant d'incorporer des informations a priori et de traiter l'incertitude de manière cohérente.

## En résumé :
#Le Lasso Bayésien offre une régularisation et une sélection automatique des variables,
#produisant des modèles épars et interprétables, mais il impose une structure stricte de
#régularisation identique pour tous les coefficients. L'Approche Bayes A, quant à elle,
#offre une flexibilité accrue en permettant des variations individuelles des coefficients 
#avec des variances spécifiques, au prix d'une estimation plus complexe via des méthodes MCMC

# Role du paramètre lambda et son choix
# -----------------------------------------------------

#Le paramètre lambda dans le Lasso Bayésien contrôle la force de la régularisation
#appliquée aux coefficients du modèle. Un lambda élevé augmente l'effet de shrinkage,
#réduisant davantage les coefficients vers zéro, ce qui peut conduire à des modèles plus épars.
#À l'inverse, un lambda faible diminue l'effet de shrinkage, permettant aux coefficients de
# rester plus proches de leurs valeurs non régularisées.

### Choix de lambda :
# - Utiliser la validation croisée pour sélectionner le lambda qui minimise l'erreur de prédiction sur un ensemble de validation.

# Estimation du modèle avec avec des hyperparamètre autour de 1 et 2
# -----------------------------------------------------
# Sélection des colonnes pour Y et X
Y <- training$Y  # Remplacez "Y" par le nom de la colonne contenant votre variable cible
X <- training[, -(1:2)]  # Sélectionnez toutes les colonnes sauf les deux premières (Y et "X")

LASSO_BLR <- BLR(y = Y,
                 XL = as.matrix(X),
                 prior = list(
                   varE = list(df = 2, S = 1),
                   lambda = list(shape = 2, rate = 1, type = 'random', value = 1.5)
                 ),
                 nIter = 3000,
                 burnIn = 2000,
                 saveAt = "BLR_",
                 thin = 1)


LASSO_BLR$mu # muchap = moyenne des trajectoires = estimation de la moyenne a posteriori
LASSO_BLR$varE # = moyenne des trajectoires = estimation de la moyenne a posteriori
LASSO_BLR$bL   # betachap = moyenne des trajectoires = estimation de la moyenne a posteriori
LASSO_BLR$tau2  # = moyenne des trajectoires = estimation de la moyenne a posteriori
LASSO_BLR$lambda # = moyenne des trajectoires = estimation de la moyenne a posteriori


# Distribution a posteriori VS priori
# -----------------------------------------------------

################## Trajectoires simul?es 
# Attention : 
# si on relance alors les trajectoires se cumulent (concat?nation des r?sultats)
# il faut soit effac? BLR_, soit changer le nom...
# Par ailleurs, on ne garde qu'une valeur sur 10 (par d?faut, pour garantir une "ind?pendance") 

tracevarE <- scan('BLR_varE.dat')
plot(tracevarE,type='o',xlab="iteration",ylab="varE",main="trace varE")
# on observe une rupture => ? explorer sur une cha?ne plus longue
# ou bien sur plusieurs cha?nes 
plot(density(tracevarE),col="blue",lwd=2, main="varE : LASSO bay?sien")
curve(dinvgamma(x,1,1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))


tracelambda <- scan('BLR_lambda.dat')
plot(tracelambda,type='o',xlab="iteration",ylab="lambda",main="trace lambda")
plot(density(tracelambda),xlim=c(0,5),col="blue",lwd=2, main="lambda : LASSO bay?sien")
curve(dgamma(x,1,0.1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))

# Prédiction et corrélation
# -----------------------------------------------------
predLASSOBLR <- predictions(test[, -(1:2)], LASSO_BLR$mu, LASSO_BLR$bL)
# Tracer la série prédite
plot(predLASSOBLR, type = "l",
     col = "blue", lwd = 2,
     ylim = range(c(predLASSOBLR, test$Y)), xlab = "Observations", ylab = "Valeurs")
# Ajouter la série réelle
lines(test$Y, col = "red", lwd = 2)  # Ajoute la série réelle sur le même graphique

cor(test$Y, predLASSOBLR)
# 0.9025047
plot(predLASSOBLR ~ test$Y)  # pr?dictions versus observations sur jeu test 
abline(0,1)

# Selection des variables
# -----------------------------------------------------
plot(LASSO_BLR$bL,xlab="Indice variable",ylab="Param?tre de r?gression associ?", main="LASSO
bay?sien")
plot(sort(abs(LASSO_BLR$bL)))
bb <- boxplot(LASSO_BLR$bL)
bb
varselectedLASSO <- subsetSelected(LASSO_BLR$bL,data[, -(1:2)], bb$stats[1,],bb$stats[5,])
varselectedLASSO

#   selected    valeurs
#1  Science.9 -0.7279170
#2    Jeux.13 -1.0970720
#3     Film.3  1.2799031
#4     Film.8  2.0021597
#5    Film.13  1.1466155
#6    Serie.8  2.1890472
#7   Sport.10  4.0140078
#8   Sport.11  2.4337661
#9   Sport.15  2.1019194
#10  Music.13  3.2808854
#11      sexe  0.8304813

# on va faire du seuillage ? la main pour eliminer quelques variables 
varselectedLASSO <- subsetSelected(LASSO_BLR$bL,data[, -(1:2)],-0.3,0.3)
varselectedLASSO


# Différences entre les variables retenues par le LASSO Bayésien et Bayes A (Ridge Bayésien)

# LASSO Bayésien :
# - Régularisation L1 favorisant la sparsité, certains coefficients sont réduits à zéro.
# - Sélectionne automatiquement un sous-ensemble de variables pertinentes.

# Bayes A (Ridge Bayésien) :
# - Régularisation L2 réduisant les coefficients sans les mettre à zéro.
# - Tient compte de toutes les variables avec des contributions plus faibles pour les moins importantes.

# Explications des différences :
# - LASSO élimine les variables non pertinentes, Ridge les garde mais réduit leur impact.
# - LASSO plus efficace pour modèles épars, Ridge pour modèles où toutes les variables sont pertinentes.

