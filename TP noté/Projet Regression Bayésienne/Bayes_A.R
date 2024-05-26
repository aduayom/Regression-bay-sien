library(MCMCpack) 
library(LearnBayes) 
library(invgamma)

# Modèle Bayes A vs Modèle Random Regression
# ------------------------------------------
# Modèle Bayes A :
# - Chaque coefficient βj suit une loi normale multivariée avec sa propre variance.
# - Les variances des coefficients et de l'erreur suivent des lois inverse Gamma.
# - L'a priori sur µ est vague.
# - Estimation des paramètres via un Gibbs sampler.

# Modèle Random Regression :
# - Tous les coefficients partagent la même variance.
# - Estimation des paramètres via d'autres méthodes, comme la méthode des moindres carrés.

# Avantages du modèle Bayes A :
# - Plus de flexibilité en permettant à chaque coefficient d'avoir sa propre variance.
# - Meilleure adaptation aux données avec des variations hétérogènes entre les coefficients.
# - Approche bayésienne permettant d'incorporer des informations a priori et de traiter l'incertitude de
#   manière cohérente.

# Avantages du modèle Random Regression :
# - Simplicité et facilité d'interprétation avec une seule variance pour tous les coefficients.
# - Estimation des paramètres peut être plus rapide et moins coûteuse en calcul.

# En résumé, le modèle Bayes A offre une flexibilité accrue en permettant des variations individuelles
# des coefficients, mais il nécessite une méthode d'estimation plus complexe. Le modèle Random Regression,
# quant à lui, est plus simple mais peut manquer de précision dans la modélisation des variations individuelles
#des coefficients.



# Mise en place du modèle Bayes
# -----------------------------------------------------
# les a priori sont des inverses gamma 
# de paramètres (a,b) et (c,d) 
# leurs moyennes et variances sont 
a=c=1 ; b=d=1
b/(a-1) ; d/(c-1)   # espérance d'une inv gamma
b^2/((a-1)^2*(a-2)) ; d^2/((c-1)^2*(c-2)) # variance inv gamma
# avec ces valeurs ont met des a priori très vague, de variance infinie 
priora <- 2
priorb <- 1
priorc <- 2
priord <- 1
resBAYESA <- BayesA(Y, as.matrix(X), priora, priorb, priorc,
                    priord, mean(Y), 2000, 1500)

# Estimation des paramètres
# -----------------------------------------------------

moybeta <- apply(resBAYESA[[1]], 1, mean) # Estimations moyennes des coefficients β.
moysigma2beta <- apply(resBAYESA[[2]], 1, mean) # Estimations moyennes des variances des coefficients β.
moymu <- mean(resBAYESA[[3]]) # Moyenne des valeurs échantillonnées pour µ.
moysigma2eps <- mean(resBAYESA[[4]]) # Moyenne des valeurs échantillonnées pour la variance de l'erreur ε.
moybeta[50] # Valeur estimée du 50ème coefficient β.

plot(moybeta) # Visualisation des estimations moyennes des coefficients β.
plot(sort(abs(moybeta))) # Visualisation des valeurs absolues des estimations moyennes des coefficients β.

plot(moysigma2beta) # Visualisation des estimations moyennes des variances des coefficients β.

moymu
moysigma2eps

# Modification des paramètres pour influencer le modèle
# -----------------------------------------------------


# Prédiction VS Valeur réelle coorélation
# -----------------------------------------------------
predBAYESA <- predictions(test[, -(1:2)], moymu, moybeta)  # sur le jeu test

# Tracer la série prédite
plot(predBAYESA, type = "l",
     col = "blue", lwd = 2,
     ylim = range(c(predBAYESA, test$Y)), xlab = "Observations", ylab = "Valeurs")
# Ajouter la série réelle
lines(test$Y, col = "red", lwd = 2)  # Ajoute la série réelle sur le même graphique

cor(test$Y, predBAYESA)
# 0.8878437
plot(predBAYESA ~ test$Y)  # pr?dictions versus observations sur jeu test 
abline(0,1)

# selection des varibales
# -----------------------------------------------------
plot(moybeta,xlab="Indice variable",ylab="Parametre de regression associe",main="Bayes
A")
plot(sort(abs(moybeta)))
bb <- boxplot(moybeta)
varselectedBAYESA <- subsetSelected(moybeta,data[, -(1:2)], bb$stats[1,],bb$stats[5,])
varselectedBAYESA

# selected   valeurs
#1  Jeux.13 -1.226395
#2   Film.3  1.105470
#3   Film.8  1.613347
#4  Serie.8  1.710540
#5 Sport.10  3.793498
#6 Sport.11  2.125721
#7 Sport.15  2.056701
#8 Music.13  3.215987
