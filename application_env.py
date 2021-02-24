import environnement_exp as env

# Test pour voir de beaux graphes de base et cool 
"""

NN = (10**6) * 2
beta,pi,mu = 2.5,1,0.01
VIRUS = env.epidemie(beta,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN-40,40,0,0,0)
epsilon = 0.1 # proportionde la population testée
duree = 20 # en jours 
Lorraine.crache_un_graphe(epsilon,duree)

"""


# Les test pour savoir si notre modèle est cohérent 

NN = (10**5) 
R,pi,mu = 8,1,1/14
VIRUS = env.maladie(R,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
epsilon = 0.1 # proportionde la population testée
duree = 50 # en jours 
Lorraine.crache_un_graphe(epsilon,duree)






