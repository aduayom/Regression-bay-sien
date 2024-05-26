library(rrBLUP)

# Random Regression (RR) : 
# -------------------------

# La Random Regression (RR) est un modèle statistique utilisé pour explorer
# les relations entre une variable dépendante (Y) et un ensemble de variables
# indépendantes (X) normalisées. Dans ce modèle, les variables explicatives sont
# supposées être normalisées, avec des moyennes nulles et des variances égales à 1. 

# Le modèle RR, introduit par Henderson en 1975, est le modèle le plus simple de la RR.
# Il modélise Y en fonction de Xβ et ε, où β suit souvent une loi normale multivariée centrée
# en zéro avec une matrice de covariance G égale à σβ²I, et ε suit souvent une loi normale multivariée
# centrée en zéro avec une matrice de covariance R égale à σε²I.

# La variance de Y est donnée par V(Y) = XGX' + R. La RR est utile pour modéliser
# des relations complexes entre des variables dépendantes et indépendantes, 
# en tenant compte de leur covariance et de leur dépendance, 
# notamment dans les cas où les données présentent une certaine structure hiérarchique ou de corrélation.

# Mise en place du modèle RR
# ---------------------------------------------------------------
# Sélection des colonnes pour Y et X
Y <- training$Y  # Remplacez "Y" par le nom de la colonne contenant votre variable cible
X <- training[, -(1:2)]  # Sélectionnez toutes les colonnes sauf les deux premières (Y et "X")

# Construction du modèle avec mixed.solve()
resBLUP <- mixed.solve(Y, Z = as.matrix(X), X = as.matrix(rep(1, 100)), method = "REML")



# Représenter les estimations des paramètres aléatoires
# -----------------------------------------------------
muchap <- resBLUP$beta  # moyenne  mu
muchap
betachap <- resBLUP$u   # les coefficients beta

# c'est un vecteur de longueur p=300
par(mfrow=c(1,1))
plot(sort(abs(betachap)))   
betachap[50]
resBLUP$Ve  # variance des epsilon
resBLUP$Vu # variance des betas


# Mise en place des prédictions sur les données test
# -----------------------------------------------------
predBLUP <- predictions(test[, -(1:2)], muchap, betachap)
# Supposons que predBLUP est votre série prédite et test$Y est votre série réelle

# Tracer la série prédite
plot(predBLUP, type = "l",
     col = "blue", lwd = 2,
     ylim = range(c(predBLUP, test$Y)), xlab = "Observations", ylab = "Valeurs")

# Ajouter la série réelle
lines(test$Y, col = "red", lwd = 2)  # Ajoute la série réelle sur le même graphique

###### corrélation en y prédits et y observés (sur jeu test)
cor(test$Y, predBLUP)
# 0.5032711
plot(predBLUP ~ test$Y)
abline(0,1)
# on s'aperçoit de l'effet shrinkage 


# Sélection des variables les plus pertinentes
# --------------------------------------------
### Sélection des coefficients 
plot(betachap,xlab="Indice variable",ylab="Paramètre de régression associé",main="rrBLUP")
plot(sort(abs(betachap)))
boxplot(betachap)
bb <- boxplot(betachap)
varselectedBLUP <- subsetSelected(betachap,data[, -(1:2)],bb$stats[1,], bb$stats[5,])
varselectedBLUP
#On retient les valeurs suivantes 

# selected   valeurs
#  Jeux.13 -1.072094
#  Serie.8  1.201440
# Sport.10  1.781806
# Sport.15  1.717928
# Sport.17  1.363854
# Music.13  2.025700