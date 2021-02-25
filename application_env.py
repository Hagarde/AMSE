import environnement_exp as env

# Test pour voir de beaux graphes de base et cool 
"""
NN = (10**5) 
beta,pi,mu = 2.5,1,0.01
VIRUS = env.maladie(beta,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN-40,40,0,0,0)
epsilon = 0.1 # proportionde la population testée
duree = 20 # en jours 
Lorraine.crache_un_graphe(epsilon,duree)
"""


# Les test pour savoir si notre modèle est cohérent 

NN = (10**5) 
R,pi,mu = 20,1,1/14
VIRUS = env.maladie(R,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
proportion_de_pop_testee = 0.5 # proportionde la population testée
duree = 50 # en jours 
Lorraine.crache_un_graphe(proportion_de_pop_testee,duree)


#Enquêtons sur les valeurs négatives ! 
"""
NN = (10**5) 
R,pi,mu = 8,1,1/14
VIRUS = env.maladie(R,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
proportion_de_la_pop_testee = 0 
duree = 50 # en jours 
for i in range (100) :
    Lorraine.evol_local_seule (proportion_de_la_pop_testee)
"""