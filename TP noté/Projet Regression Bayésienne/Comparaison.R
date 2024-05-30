library(glmnet)
library(stargazer)

# Classez les différentes approches
# -----------------------------------------------------


# Comparaison correlations entre valeurs predites et observees
# -----------------------------------------------------


# Comparer les variables obtenues par les 4 methodes
# -----------------------------------------------------

# Approches non bayésien: Lasso et Elastinet
# -----------------------------------------------------
# Sélection des colonnes pour Y et X
Y_data <- data$Y  
X_data <- data[, -(1:2)]
# 1. Regression LASSO 
res=cv.glmnet(as.matrix(X_data),Y_data,family='gaussian',alpha=1) 
plot(res)

# on va selectionner le meilleur lambda
# on prend le lambda.1se car c'est du LASSO (=> moins de s?lection)
seuil=res$lambda.1se
reslasso=glmnet(as.matrix(X_data),Y_data,family='gaussian',alpha=1,lambda=seuil)
# on extrait les coefficients moins la constante 
coeflasso=coefficients(reslasso)[-1]
plot(sort(abs(coeflasso)))
# regardons le boxplot des coef non nuls 
coefnonnul=coeflasso[coeflasso!=0]
boxplot(coefnonnul)

# selection des variables 
selec1=order(coeflasso,decreasing=TRUE)[1:8]
colnames(X_data[,selec1])
cor(X_data[,selec1])

# 2. Regression ELASTINET
res=cv.glmnet(as.matrix(X_data),Y_data,family='gaussian',alpha=0.5) 
plot(res)

# on va sElectionner le meilleur lambda
# on prend le lambda.1se car c'est du LASSO (=> moins de sElection)
seuil=res$lambda.1se
resnet=glmnet(as.matrix(X_data),Y_data,family='gaussian',alpha=0.5,lambda=seuil)
# on extrait les coefficients moins la constante 
coefnet=coefficients(resnet)[-1]
plot(sort(abs(coefnet)))
# on a arbitrairement envie de conserver 1 + 4 + 3  coefficients
# regardons le boxplot des coef non nuls 
coefnonnul=coefnet[coefnet!=0]
boxplot(coefnonnul)
selecnet=order(coefnet,decreasing=TRUE)[1:10]
colnames(X_data[,selecnet])

# 3. Regression Standard
# Noms des variables sélectionnées
selected_vars <- colnames(X_data[, selecnet])

# Créez la formule de régression
# Y ~ var1 + var2 + var3 + ...
formula <- as.formula(paste("Y_data ~", paste(selected_vars, collapse = " + ")))

# Ajustez le modèle de régression OLS
model <- lm(formula, data = data)

# Résumé du modèle
summary(model)
