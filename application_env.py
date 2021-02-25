import environnement_exp as env

# Test pour voir de beaux graphes de base et cool 
"""
NN = (10**5) 
beta,pi,mu = 2.5,1,0.01
VIRUS = env.maladie(beta,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN-40,40,0,0,0)
proportion_test= 0.1 # proportionde la population testée
duree = 20 # en jours 
Lorraine.crache_un_graphe(proportion_test,duree)
"""


# Les test pour savoir si notre modèle est cohérent 
"""
NN = (10**5) 
R,pi,mu = 20,1,1/14
VIRUS = env.maladie(R,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
proportion_test = 0.5 # proportionde la population testée
duree = 50 # en jours 
Lorraine.crache_un_graphe(proportion_test,duree)

"""
#Enquêtons sur les valeurs négatives ! 
"""
NN = (10**5) 
R,pi,mu = 8,1,1/14
VIRUS = env.maladie(R,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
proportion_test = 0 
duree = 50 # en jours 
for i in range (100) :
    Lorraine.evol_local_seule (proportion_test)
"""

"""
# Définition du virus identique aux deux environnement
NN = (10**5) 
R,pi,mu = 10,1,1
VIRUS = env.maladie(R,pi,mu)
proportion_test = 0.1 
duree = 100 # en jours 

# Graphique avec raisonnement discret 

Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
Lorraine.crache_un_graphe(proportion_test,duree)

# Graphique  avec raisonnement continu
VIRUS = env.maladie(R,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
Lorraine.crache_un_graphe_continu(proportion_test,duree)
"""

# Test de la fonction qui détermine si l'épidémie est contrôllable 

#Définition de l'env 
NN = (10**5) 
R,pi,mu = 10,1,1/14
VIRUS = env.maladie(R,pi,mu)
proportion_test = 0.1 
duree_max = 100 # en jours 


Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
Lorraine.evol_local_seule_continu(proportion_test)
print (Lorraine.history)
# Phase de test de l'épidémie 

Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
vecteur_controle = Lorraine.determiner_controllabilite(proportion_test,duree_max)
print(vecteur_controle)