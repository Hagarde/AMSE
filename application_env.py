import environnement_exp as env
import  numpy as np 
import pandas as pd
import time 
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
R,pi,mu = 10,1,1/15
VIRUS = env.maladie(R,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
proportion_test = 0.1 # proportionde la population testée
duree = 100 # en jours 
Lorraine.crache_un_graphe_continu(proportion_test,duree)
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
R,pi,mu = 5,0,1/14
VIRUS = env.maladie(R,pi,mu)
proportion_test = 0
duree = 100 # en jours 

# Graphique avec raisonnement discret 

Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
Lorraine.crache_un_graphe(proportion_test,duree)

# Graphique  avec raisonnement continu
VIRUS = env.maladie(R,pi,mu)
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
Lorraine.crache_un_graphe_continu(proportion_test,duree)
"""
#  La fonction de controllabilité est plutot ok don con va faire un petit grille search 

# Définissons les maladies que nous allons étudiés
"""
NN = 10 ** 4
grille = []
for proportion_test in [0,0.01,0.05] :
    for R in range (3,21) :
        print(R)
        for PI in np.arange(0.1,1.1,0.1) :
            for MU in np.arange (5,30,5) : 
                MU = 1/MU
                coronavirus = env.maladie(R,PI,MU)
                region = env.env_minimal("Lorraine",NN,coronavirus,0,NN*0.95,NN*0.05,0,0,0)
                grille.append(region.determiner_controllabilite (proportion_test,100))
grille = np.array(grille)
np.savetxt("grille_search.txt",grille,delimiter=";",header='R,mu,pi,proportion testée,infecté max, date, Infecté par maladie , infectée au bout de 100 jours ')
"""


#Je teste le monstre de calcul avec 4 régions et je vais corriger petit à petit les milliard ed bug 
# Il faut multiplier par 4 beta car la population 4 fois plus importante il y a qqch que je ne comprends pas de ouf 
R,pi,mu = 10,1,1/15
peste = env.maladie(R,pi,mu)
matrice_influence = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
NN1,NN2,NN3,NN4 = 10**5,10**5,10**5,10**5
S01,U01,P01,R0_U1,R0_P1 = NN1 * 0.95 ,NN1* 0.05 ,0,0,0
S02,U02,P02,R0_U2,R0_P2 = NN2 * 0.95 ,NN2* 0.05 ,0,0,0
S03,U03,P03,R0_U3,R0_P3 = NN3 * 0.95 ,NN3* 0.05 ,0,0,0
S04,U04,P04,R0_U4,R0_P4 = NN4 * 0.95 ,NN4* 0.05 ,0,0,0
# question sur le beta qui doit-être divisé par pop  totale ou locale 
monde = env.env_total(NN1,NN2,NN3,NN4,peste,matrice_influence,S01,U01,P01,R0_U1,R0_P1,S02,U02,P02,R0_U2,R0_P2,S03,U03,P03,R0_U3,R0_P3,S04,U04,P04,R0_U4,R0_P4)
monde.graphe_sur_charge(0.1,0.25,0.4,0.5)


"""
NN = (10**5) 
R,pi,mu = 15,1,1/14
VIRUS = env.maladie(R,pi,mu)
proportion_test = 0.25
VIRUS = env.maladie(R,pi,mu)
duree = 100
Lorraine = env.env_minimal("Lorraine",NN,VIRUS,0,NN*0.95,NN*0.05,0,0,0)
Lorraine.crache_un_graphe_continu(proportion_test,duree)
print(Lorraine.determiner_controllabilite (proportion_test,100))
"""

