# Stochastic Search Variable Selection (SSVS)

# SSVS, introduit par George et McCulloch (1993), utilise un paramètre γ (spike and slab) pour indiquer si un coefficient β est sélectionné.
# γ est un vecteur de dimension p avec des 1 pour les coefficients βj non nuls et des 0 pour les coefficients βj nuls.
# βγ représente les coefficients non nuls et Xγ la matrice design associée.
# Le modèle SSVS est utilisé pour la sélection de variables parmi de nombreux candidats.

# Modèle :
# Y = μI + Xβ + ϵ
# ϵ ∼ N(0, σϵ²I)
# γk = 1 si βk ≠ 0 (sélection)
# γk = 0 si βk = 0 (non sélection)
# βγ | γ, σϵ² ∼ N(0, σϵ²c (Xγ'Xγ)⁻¹)
# P(γj = 1) = π
# f(σϵ²) = 1/σϵ²
# μ ∼ uniform

# Le paramètre c (>0) est un facteur d'échelle, généralement entre 10 et 100.
# Le paramètre π (0<π<1) contrôle le nombre de βj sélectionnés.
# Un π petit limite les sélections, un π grand en sélectionne trop.
# β a un a priori de Zellner avec une structure de covariance gérée par σϵ²c (Xγ'Xγ)⁻¹.
# En cas de matrice Xγ'Xγ mal conditionnée, on utilise Ridge-Zellner : βγ | γ, σϵ² ∼ N(0, σϵ²c(Xγ'Xγ + λI)⁻¹).

# SSVS est utile pour sélectionner les meilleures variables dans des modèles complexes, comme en actuariat ou modélisation de la mortalité.


# Selection des variables
# -----------------------------------------------------
Y <- training$Y  # Remplacez "Y" par le nom de la colonne contenant votre variable cible
X <- training[, -(1:2)]  # Sélectionnez toutes les colonnes sauf les deux premières (Y et "X")

nbinit <- 10 # nombre de beta non nuls au d?part 
nbToChange <- 2 # nombre de beta que l'on propose de changer a chaque fois 
Pi <- 10/160   # la probabilite de choisir un beta (en moyenne 10 car il y en a 300)
lec <- 50  # preconisation : entre 10 et 100 
resSSVS <- selection_SSVS(Y, as.matrix(X), 8000, 7000, lec,
                          nbinit, nbToChange, Pi)
resSSVS[[2]]  # on voit combien il y a eu d'echanges dans les 1000 derniers run
resSSVS[[3]]  # on voit le nbre de variables selectionnees lors de ces echanges 
resSSVS[[1]]

####### S?lection 
plot(sort(abs(resSSVS[[1]])))
boxplot(resSSVS[[1]])
plot(resSSVS[[1]],xlab="Indice variable",ylab= "Nombre de s?lections post-burn-in", main="SSVS")
varselectedSSVS <- subsetSelected(resSSVS[[1]],X,0,999)
varselectedSSVS
