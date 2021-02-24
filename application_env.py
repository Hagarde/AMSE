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


# Un appel de la fonction qui renvoie le nombre d'infecté au bout d'une durée pour un env 


NN = (10**6) * 2
beta,pi,mu = 5,1,0.5
VIRUS = env.epidemie(beta,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN-40,40,0,0,0)
epsilon = 0.5 # proportionde la population testée
duree = 15 # en jours 
Lorraine.crache_un_graphe(epsilon,duree)
print (Lorraine.resultat_politique_epidemie(epsilon, duree))


