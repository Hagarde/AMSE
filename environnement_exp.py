from random import uniform
from random import gauss
import numpy as np 
import matplotlib.pyplot as plt
 
class maladie :
    def __init__(self,R,pi,mu) :
        self.R = R
        self.beta = R*mu
        self.pi = pi 
        self.mu = mu 

    def get_beta(self) :
            return self.beta

    def get_pi(self) :
        return self.pi

    def get_mu(self) :
        return self.mu

    def get_R(self) :
        return self.R
   
class env_total: 
    def __init__(self,name1,name2,NN1,NN2,beta,pi,mu,influence_inter_regionale,S01,U01,R0_U1,R0_P1,P01,S02,U02,R0_U2,R0_P2,P02) :
        self.region1 = env_minimal(name1,NN1,beta,pi,mu,influence_inter_regionale,S01,U01,R0_U1,R0_P1,P01)
        self.region2 = env_minimal(name2,NN2,beta,pi,mu,influence_inter_regionale,S02,U02,R0_U2,R0_P2,P02)
        self.population = NN1+NN2
        self.beta = beta / self.population
        self.pi = pi 
        self.mu = mu 
        self.influence_inter_regionale = influence_inter_regionale
    def action_global(self,test1,test2) : #j'ai mis les régions en caractéristique de env_global puis j'y accède par Bou..region1.method
        # problème de calcul l'un après l'autre je crois pq les caractéristiques de region1 vont s'actualisé avant que region2 ait pu calculé 
        Bouches_Rhones_et_Lorraine.region1.evol_local(test1,Bouches_Rhones_et_Lorraine.region2)
        Bouches_Rhones_et_Lorraine.region2.evol_local(test2,Bouches_Rhones_et_Lorraine.region1)

class env_minimal :
    def __init__(self,name,NN,virus,influence_inter_regionale,S0,U0,P0,R0_U,R0_P) :
        self.history = [[S0],[U0],[P0],[R0_U],[R0_P]]
        self.name = name 
        self.population = NN
        self.virus = virus
        virus.beta = virus.beta / self.population
        self.influence_inter_regionale = influence_inter_regionale
        self.S0 = S0
        self.P0 =P0
        self.U0 = U0
        self.R0_P = R0_P
        self.R0_U = R0_U 

# Les fonctions un peu chiantes des get qui sont utiles mais c'est surtout pour m'entrainer 
    def get_population(self) :
        return self.population
    def get_susceptible(self) :
        return self.S0
    def get_name(self) :
        return self.name
    def get_positif(self) :
        return self.P0
    def get_recovered_positif(self) :
        return self.RP_0
    def get_recovered_undetected(self) :
        return self.RU_0
    def get_undetected(self) :
        return self.U0
    def get_history(self) :
        return self.history

#Les fonctions de modifications des valeurs de la classe ,je vous jure après les fonctions sont mieux ! 
    def evol_local(self,nmbr_test,region_autre) :
        dS = - self.virus.beta* self.S0 * (self.U0 + (1-self.virus.pi)*self.P0) #  on enlève les personnes infectées 
        dU = -dS - self.virus.mu*self.U0 - (nmbr_test)*(self.U0/(self.U0+self.R0_U+self.S0)) + self.virus.beta*(self.influence_inter_regionale * region_autre.U0)# (epsilon * NN )*(U0/(U0+R0+S0) le nombr de test * la proportion d'infectée dans la population qu'il reste à tester donc c'est le nombre de personne détectées positives à la fin des test
        dP1 = (nmbr_test)*(self.U0/(self.U0+self.R0_U+self.S0)) # on enlève les guéries et nmbr_test*Prbl de tomber sur un +
        dP = (nmbr_test )*(self.U0/(self.U0+self.R0_U+self.S0)) - self.virus.mu * self.P0 
        dR_U = self.virus.mu * self.virus.U0  # on ajoute les personnes guéries et qui vont se refaire testées
        dR_P = self.virus.mu * self.virus.P0 # on ajoute les personnes guéries sans le savoir
        dS = (dP+dU+dR_U+dR_P) * (-1)
        self.S0 = self.S0 + dS 
        self.U0 = self.U0 + dU 
        self.P0 = self.P0 + dP 
        self.R0_P = self.R0_P + dR_P
        self.R0_U = self.R0_U + dR_U
        self.history[0].append(self.S0)
        self.history[1].append(self.U0)
        self.history[2].append(self.P0)
        self.history[3].append(self.R0_U)
        self.history[4].append(self.R0_P)
        if  (self.S0 < 0): # ça va poser prbl pour les simulation donc peut-être renvoyer des données nulles 
            return "il y a un problème soit la pandémie est finie soit tout le monde est contaminé"
        else : 
            print("Bonjour ! Dans " + self.name+ ", il y a actuellement d'après les tests effectués "+ str(self.get_name())+" individus positifs ! ")
            print ("U0 = " + str(self.get_U0()))
            print ("S0 = " + str(self.get_S0()))
            print ("P0 = " + str(self.get_P0()))
            print ("R0_P = " + str(self.get_R0_P()))
            print ("R0_U = " + str(self.get_R0_U()))
            print ("beta = " + str(self.virus.get_beta()))
            print ("pi = " + str(self.virus.get_pi()))
            print ("mu = " + str(self.virus.get_mu()))

    def evol_local_seule (self,proportion_de_la_pop_testee):
        # Condition pour pas faire des calculs inutiles ou incohérent
        if  (self.U0 < (self.population)*(10**-10)) or(self.U0 > self.population):
            nmbr_test = proportion_de_la_pop_testee * self.population
            dS = - self.virus.beta* self.S0 * (self.U0 + (1-self.virus.pi)*self.P0) # *(self.S0/self.population) 
            dU = -dS - self.virus.mu*self.U0 - (nmbr_test)*(self.U0/(self.U0+self.R0_U+self.S0)) # (epsilon * NN )*(U0/(U0+R0+S0) le nombr de test * la proportion d'infectée dans la population qu'il reste à tester donc c'est le nombre de personne détectées positives à la fin des test
            dP1 = (nmbr_test)*(self.U0/(self.U0+self.R0_U+self.S0)) # on enlève les guéries et nmbr_test*Prbl de tomber sur un +
            dP = (nmbr_test )*(self.U0/(self.U0+self.R0_U+self.S0)) - self.virus.mu * self.P0
            dR_U = self.virus.mu * self.U0  # on ajoute les personnes guéries et qui vont se refaire testées
            dR_P = self.virus.mu * self.P0 # on ajoute les personnes guéries sans le savoir
            dS = (dP+dU+dR_U+dR_P) * (-1)

                #Actualistaion des facteurs 
            self.S0 = self.S0 + dS 
            self.U0 = self.U0 + dU 
            self.P0 = self.P0 + dP 
            self.R0_P = self.R0_P + dR_P
            self.R0_U = self.R0_U + dR_U
            self.history[0].append(self.S0)
            self.history[1].append(self.U0)
            self.history[2].append(self.P0)
            self.history[3].append(self.R0_U)
            self.history[4].append(self.R0_P)
            # AVERTISSEMENT DE NEGATIVITE IMPORTNANT 
            print ([row[-2] for row in self.history])
            print ([row[-1] for row in self.history])
            print(len(self.history[0]))


        else : 
                # Calcul des variations 
            nmbr_test = proportion_de_la_pop_testee * self.population
            dS = - self.virus.beta* self.S0 * (self.U0 + (1-self.virus.pi)*self.P0) # *(self.S0/self.population) 
            dU = -dS - self.virus.mu*self.U0 - (nmbr_test)*(self.U0/(self.U0+self.R0_U+self.S0)) # (epsilon * NN )*(U0/(U0+R0+S0) le nombr de test * la proportion d'infectée dans la population qu'il reste à tester donc c'est le nombre de personne détectées positives à la fin des test
            dP1 = (nmbr_test)*(self.U0/(self.U0+self.R0_U+self.S0)) # on enlève les guéries et nmbr_test*Prbl de tomber sur un +
            dP = (nmbr_test )*(self.U0/(self.U0+self.R0_U+self.S0)) - self.virus.mu * self.P0
            dR_U = self.virus.mu * self.U0  # on ajoute les personnes guéries et qui vont se refaire testées
            dR_P = self.virus.mu * self.P0 # on ajoute les personnes guéries sans le savoir
            dS = (dP+dU+dR_U+dR_P) * (-1)

                #Actualistaion des facteurs 
            self.S0 = self.S0 + dS 
            self.U0 = self.U0 + dU 
            self.P0 = self.P0 + dP 
            self.R0_P = self.R0_P + dR_P
            self.R0_U = self.R0_U + dR_U
            self.history[0].append(self.S0)
            self.history[1].append(self.U0)
            self.history[2].append(self.P0)
            self.history[3].append(self.R0_U)
            self.history[4].append(self.R0_P)

    def crache_un_graphe (self,proportion_de_la_pop_testee,duree_de_experience) :
            # duree doit être un nombre de jour entier 
        TEMPS = np.arange(0,duree_de_experience,1)
        for i in range (duree_de_experience-1) :
            self.evol_local_seule(proportion_de_la_pop_testee) 
        data = np.array(self.get_history())
        data = data / self.population
        fig, ax = plt.subplots(nrows=1, ncols=1)
        S, = ax.plot(TEMPS, data[0], marker='+', color='blue', label='S0')
        U ,= ax.plot(TEMPS, data[1], marker='+', color='red', label='U0')
        P ,= ax.plot(TEMPS, data[2], marker='+', color='yellow', label='P0')
        RU, = ax.plot(TEMPS, data[3], marker='+', color='green', label='RO_U')
        RP ,= ax.plot(TEMPS, data[4], marker='+', color='violet', label='R0_P')
        R, = ax.plot(TEMPS, data[4] + data [3], marker='+', color='orange', label='RP+RU')
        ax.set(xlabel='Temps (en jours)', ylabel='Proportion de la population', title="Evolution pour beta = "+ str(self.virus.beta)+ " / pi = " + str(self.virus.pi) + " / mu = " + str(self.virus.mu))
        plt.legend([S,U,P,RU,RP,R], ['S', 'U', 'P','R_U','R_P','R total'], loc='best')
        plt.show()

    def resultat_politique_epidemie (self,proportion_de_la_pop_testee,duree_de_experience) :
        for i in range(duree_de_experience) :
            self.evol_local_seule(proportion_de_la_pop_testee)

        if self.P0 > 0 :
            return [self.U0 + self.P0*(1-self.virus.pi),self.virus.beta,self.virus.pi,self.virus.mu]
        else :
            return [self.U0 ,self.virus.beta,self.virus.pi,self.virus.mu]
    
    def determiner_controllabilite (self,proportion_de_la_pop_testee,time_limit) :
        # renvoie un vecteur comprenant la controlabilité la durée pour controlable et les paramètres du virus 
        i = 0
        ok = 0

        while self.P0+self.U0 > self.population*0.001 and i < time_limit and self.P0+ self.U0Ò< self.population: 
            self.evol_local_seule (proportion_de_la_pop_testee)
            print(self.P0)
            print(self.U0)
            i = i + 1
        if (self.P0+ self.U0<w <self.population* 0.00001 ) : 
            ok = 1
        duree = i
        print (i)
        return [i,ok,self.virus.R,self.virus.mu,proportion_de_la_pop_testee,self.virus.pi]