set.seed(123)

# Chargement du fichier CSV
data <- read.csv("C:/Users/Daniel/Desktop/Mastère ENSAI/Regression Bayésienne/TP noté/telecat.csv")

# Affichage des premières lignes pour vérifier le chargement
head(data)

# Résumé statistique des variables numériques
summary(data)

# Description de la base de données "telecat.csv":
# Le jeu de données « telecat.csv » comporte les résultats d'une enquête de satisfaction menée auprès
# 150 clients d'une chaîne câblée.
# Pour chaque client un score de satisfaction a été établi, noté Y, sur une échelle qui va de moins
# satisfait (valeurs négatives) à plus satisfait (valeurs positives).
# Les variables explicatives étudiées sont les temps passés sur différentes chaînes, combinés aux
# nombres de visites sur ces chaînes. Ces covariables ont été normalisées. Il y en a p=160.

View(data)


# Division de la base de données en jeu d'apprentissage et jeu de test
# ---------------------------------------------------------------

# Étape 1 : Définition de la taille des échantillons
# On spécifie le nombre d'observations à inclure dans le jeu d'apprentissage et le jeu de test.
n_training <- 100  # Taille de l'échantillon d'apprentissage
n_test <- 50       # Taille de l'échantillon de test

# Étape 2 : Sélection aléatoire des indices pour l'échantillon d'apprentissage
# On utilise la fonction sample() pour sélectionner aléatoirement les indices des observations pour l'échantillon d'apprentissage.
indices_training <- sample(1:nrow(data), n_training)

# Étape 3 : Création du jeu de données training
# On extrait les lignes correspondant aux indices sélectionnés pour l'échantillon d'apprentissage et on les stocke dans un nouveau data frame.
training <- data[indices_training, ]

# Étape 4 : Sélection des indices pour l'échantillon de test
# On utilise la fonction setdiff() pour obtenir les indices qui ne sont pas dans l'échantillon d'apprentissage.
indices_test <- setdiff(1:nrow(data), indices_training)

# Étape 5 : Création du jeu de données test
# On extrait les lignes correspondant aux indices sélectionnés pour l'échantillon de test et on les stocke dans un nouveau data frame.
test <- data[indices_test, ]

