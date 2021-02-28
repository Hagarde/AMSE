import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint 

class maladie :
    def __init__(self,R,pi,mu) :
        self.R = R
        self.beta = R*mu
        self.pi = pi 
        self.mu = mu 
#definition des fonctions usuelles inutiles
    def get_beta(self) :
            return self.beta

    def get_pi(self) :
        return self.pi

    def get_mu(self) :
        return self.mu

    def get_R(self) :
        return self.R
   
class env_total: 
    def __init__(self,NN1,NN2,NN3,NN4,virus,matrice_influence,S01,U01,P01,R0_U1,R0_P1,S02,U02,P02,R0_U2,R0_P2,S03,U03,P03,R0_U3,R0_P3,S04,U04,P04,R0_U4,R0_P4) :
        self.population = NN1+NN2+NN3+NN4
        self.pi = virus.pi 
        self.mu = virus.mu 
        self.virus = virus 
        self.beta1 = virus.beta/NN1
        self.beta2 = virus.beta/NN2
        self.beta3 = virus.beta/NN3
        self.beta4 = virus.beta/NN4
        self.influence_1_2 = matrice_influence[0][1]
        self.influence_1_3 = matrice_influence[0][2]
        self.influence_1_4 = matrice_influence[0][3]
        self.influence_2_1 = matrice_influence[1][0]
        self.influence_2_3 = matrice_influence[1][2]
        self.influence_2_4 = matrice_influence[1][3]
        self.influence_3_1 = matrice_influence[2][0]
        self.influence_3_2 = matrice_influence[2][1]
        self.influence_3_4 = matrice_influence[2][3]
        self.influence_4_1 = matrice_influence[3][0]
        self.influence_4_2 = matrice_influence[3][1]
        self.influence_4_3 = matrice_influence[3][2]
        self.S1 = S01
        self.S2 = S02
        self.S3 = S03
        self.S4 = S04
        self.U1 = U01
        self.U2 = U02
        self.U3 = U03
        self.U4 = U04
        self.P1 = P01
        self.P2 = P02
        self.P3 = P03
        self.P4 = P04
        self.R_U1 = R0_U1
        self.R_U2 = R0_U2
        self.R_U3 = R0_U3
        self.R_U4 = R0_U4
        self.R_P1 = R0_P1
        self.R_P2 = R0_P2
        self.R_P3 = R0_P3
        self.R_P4 = R0_P4
        self.NN1= NN1 
        self.NN2=NN2
        self.NN3=NN3
        self.NN4=NN4
        self.history = [[S01,U01,P01,R0_U1,R0_P1],[S02,U02,P02,R0_U2,R0_P2],[S03,U03,P03,R0_U3,R0_P3],[S04,U04,R0_U4,R0_P4,P04]]

    def je_veux_juste_le_lendemain (self,prop_test1,prop_test2,prop_test3,prop_test4) :
        test1 = prop_test1 * self.NN1
        test2 = prop_test2 * self.NN2
        test3 = prop_test3 * self.NN3
        test4 = prop_test4 * self.NN4
        def systeme_diff (vecteur_condition_initiale,t) :
            S01,U01,P01,R0_U1,R0_P1,S02,U02,P02,R0_U2,R0_P2,S03,U03,P03,R0_U3,R0_P3,S04,U04,P04,R0_U4,R0_P4 = vecteur_condition_initiale
            dS1 = - self.beta1 * S01 * (U01 + (1-self.virus.pi)*P01)
            dS2 = - self.beta2 * S02 * (U02 + (1-self.virus.pi)*P02)
            dS3 = - self.beta3 * S03 * (U03 + (1-self.virus.pi)*P03)
            dS4 = - self.beta4 * S04 * (U04 + (1-self.virus.pi)*P04) 
            dR_U1 = self.virus.mu * U01 
            dR_U2 = self.virus.mu * U02 
            dR_U3 = self.virus.mu * U03 
            dR_U4 = self.virus.mu * U04  
            dR_P1 = self.virus.mu * P01
            dR_P2 = self.virus.mu * P02
            dR_P3 = self.virus.mu * P03
            dR_P4 = self.virus.mu * P04
            dP1 = (test1)*((U01/(U01+R0_U1+S01))) - self.virus.mu * P01
            dP2 = (test2)*((U02/(U02+R0_U2+S02))) - self.virus.mu * P02
            dP3 = (test3)*((U03/(U03+R0_U3+S03))) - self.virus.mu * P03
            dP4 = (test4)*((U04/(U04+R0_U4+S04))) - self.virus.mu * P04
            dU1 = -dS1 - self.virus.mu*U01 - (test1)*(U01/(U01+R0_U1+S01)) + U02 * self.influence_1_2* self.beta2 + U03 * self.influence_1_3 * self.beta3 + U04 * self.influence_1_4 *self.beta4
            dU2 = -dS2 - self.virus.mu*U02 - (test2)*(U02/(U02+R0_U2+S02)) + U01 * self.influence_2_1 * self.beta1 + U03 * self.influence_2_3 *self.beta3+ U04 * self.influence_2_4*self.beta4 
            dU3 = -dS3 - self.virus.mu*U03 - (test3)*(U03/(U03+R0_U3+S03)) + U01 * self.influence_3_1 * self.beta1+ U02 * self.influence_3_2 * self.beta2 + U04 * self.influence_3_4*self.beta4
            dU4 = -dS4 - self.virus.mu*U04 - (test4)*(U04/(U04+R0_U4+S04)) + U01 * self.influence_4_1 * self.beta1+ U02 * self.influence_4_2 * self.beta2 + U03 * self.influence_4_3*self.beta3 
            return [dS1,dU1,dP1,dR_U1,dR_P1,dS2,dU2,dP2,dR_U2,dR_P2,dS3,dU3,dP3,dR_U3,dR_P3,dS4,dU4,dP4,dR_U4,dR_P4]
        T  = np.arange(0,100,1)
        data = odeint (systeme_diff,[self.S1,self.U1,self.P1,self.R_U1,self.R_P1,self.S2,self.U2,self.P2,self.R_U2,self.R_P2,self.S3,self.U3,self.P3,self.R_U3,self.R_P3,self.S4,self.U4,self.P4,self.R_U4,self.R_P4],T)
    # Actualisation du milion et demi de valeur !!
        self.S1 = data[1,0] 
        self.U1 = data[1,1]
        self.P1 = data[1,2]
        self.R1_U = data[1,3]
        self.R1_P = data[1,4]
        self.history[0][0].append(self.S1)
        self.history[0][1].append(self.U1)
        self.history[0][2].append(self.P1)
        self.history[0][3].append(self.R1_U)
        self.history[0][4].append(self.R1_P)

        self.S2 = data[1,5] 
        self.U2 = data[1,6]
        self.P2 = data[1,7]
        self.R2_U = data[1,8]
        self.R2_P = data[1,9]
        self.history[1][0].append(self.S2)
        self.history[1][1].append(self.U2)
        self.history[1][2].append(self.P2)
        self.history[1][3].append(self.R2_U)
        self.history[1][4].append(self.R2_P)

        self.S3 = data[1,10] 
        self.U3 = data[1,11]
        self.P3 = data[1,12]
        self.R3_U = data[1,13]
        self.R3_P = data[1,14]
        self.history[2][0].append(self.S3)
        self.history[2][1].append(self.U3)
        self.history[2][2].append(self.P3)
        self.history[2][3].append(self.R3_U)
        self.history[2][4].append(self.R3_P)

        self.S4 = data[1,15] 
        self.U4 = data[1,16]
        self.P4 = data[1,17]
        self.R4_U = data[1,18]
        self.R4_P = data[1,19]
        self.history[3][0].append(self.S4)
        self.history[3][1].append(self.U4)
        self.history[3][2].append(self.P4)
        self.history[3][3].append(self.R4_U)
        self.history[3][4].append(self.R4_P)

    def graphe_sur_charge(self,prop_test1,prop_test2,prop_test3,prop_test4) :
        test1 = prop_test1 * self.NN1
        test2 = prop_test2 * self.NN2
        test3 = prop_test3 * self.NN3
        test4 = prop_test4 * self.NN4
        def systeme_diff (vecteur_condition_initiale,t) :
            S01,U01,P01,R0_U1,R0_P1,S02,U02,P02,R0_U2,R0_P2,S03,U03,P03,R0_U3,R0_P3,S04,U04,P04,R0_U4,R0_P4 = vecteur_condition_initiale
            dS1 = - self.beta1 * S01 * (U01 + (1-self.virus.pi)*P01)
            dS2 = - self.beta2 * S02 * (U02 + (1-self.virus.pi)*P02)
            dS3 = - self.beta3 * S03 * (U03 + (1-self.virus.pi)*P03)
            dS4 = - self.beta4 * S04 * (U04 + (1-self.virus.pi)*P04) 
            dR_U1 = self.virus.mu * U01 
            dR_U2 = self.virus.mu * U02 
            dR_U3 = self.virus.mu * U03 
            dR_U4 = self.virus.mu * U04  
            dR_P1 = self.virus.mu * P01
            dR_P2 = self.virus.mu * P02
            dR_P3 = self.virus.mu * P03
            dR_P4 = self.virus.mu * P04
            dP1 = (test1)*((U01/(U01+R0_U1+S01))) - self.virus.mu * P01
            dP2 = (test2)*((U02/(U02+R0_U2+S02))) - self.virus.mu * P02
            dP3 = (test3)*((U03/(U03+R0_U3+S03))) - self.virus.mu * P03
            dP4 = (test4)*((U04/(U04+R0_U4+S04))) - self.virus.mu * P04
            dU1 = -dS1 - self.virus.mu*U01 - (test1)*(U01/(U01+R0_U1+S01)) + U02 * self.influence_1_2* self.beta2 + U03 * self.influence_1_3 * self.beta3 + U04 * self.influence_1_4 *self.beta4
            dU2 = -dS2 - self.virus.mu*U02 - (test2)*(U02/(U02+R0_U2+S02)) + U01 * self.influence_2_1 * self.beta1 + U03 * self.influence_2_3 *self.beta3+ U04 * self.influence_2_4*self.beta4 
            dU3 = -dS3 - self.virus.mu*U03 - (test3)*(U03/(U03+R0_U3+S03)) + U01 * self.influence_3_1 * self.beta1+ U02 * self.influence_3_2 * self.beta2 + U04 * self.influence_3_4*self.beta4
            dU4 = -dS4 - self.virus.mu*U04 - (test4)*(U04/(U04+R0_U4+S04)) + U01 * self.influence_4_1 * self.beta1+ U02 * self.influence_4_2 * self.beta2 + U03 * self.influence_4_3*self.beta3 
            return [dS1,dU1,dP1,dR_U1,dR_P1,dS2,dU2,dP2,dR_U2,dR_P2,dS3,dU3,dP3,dR_U3,dR_P3,dS4,dU4,dP4,dR_U4,dR_P4]
        T  = np.arange(0,100,1)
        data = odeint (systeme_diff,[self.S1,self.U1,self.P1,self.R_U1,self.R_P1,self.S2,self.U2,self.P2,self.R_U2,self.R_P2,self.S3,self.U3,self.P3,self.R_U3,self.R_P3,self.S4,self.U4,self.P4,self.R_U4,self.R_P4],T)
    # On extrait les données pour préparer la représentation graphique 
        s1 = data[:,0] 
        u1 = data[:,1]
        p1 = data[:,2]
        r_u1 = data[:,3] 
        r_p1 =  data[:,4] 

        s2 = data[:,5] 
        u2 = data[:,6]
        p2 = data[:,7]
        r_u2 = data[:,8] 
        r_p2 =  data[:,9] 

        s3 = data[:,10] 
        u3 = data[:,11]
        p3 = data[:,12]
        r_u3 = data[:,13] 
        r_p3 =  data[:,14] 

        s4 = data[:,15] 
        u4 = data[:,16]
        p4 = data[:,17]
        r_u4 = data[:,18] 
        r_p4 =  data[:,19] 
        
    # Création des graphiques
        plt.figure(1)
        plt.suptitle("Simulation d'une épidémie dans un environnement de 4 régions connectée")

        plt.subplot(2,2,1)
        plt.plot(T,s1, color = 'blue')
        plt.plot(T,u1, color = 'red')
        plt.plot(T,p1, color = 'yellow')
        plt.plot(T,r_u1, color = 'green')
        plt.plot(T,r_p1, color = 'violet')
        plt.title("Région 1 Jeune : beta = "+ str(self.virus.beta)+"/ pi = "+ str(self.virus.pi) + "/ mu = "+ str(self.virus.mu))
        
        plt.subplot(2,2,2)
        plt.plot(T,s2, color = 'blue')
        plt.plot(T,u2, color = 'red')
        plt.plot(T,p2, color = 'yellow')
        plt.plot(T,r_u2, color = 'green')
        plt.plot(T,r_p2, color = 'violet')
        plt.title("Région 1 Vieux : beta = "+ str(self.virus.beta)+"/ pi = "+ str(self.virus.pi) + "/ mu = "+ str(self.virus.mu))

        plt.subplot(2,2,3)
        plt.plot(T,s3, color = 'blue')
        plt.plot(T,u3, color = 'red')
        plt.plot(T,p3, color = 'yellow')
        plt.plot(T,r_u3, color = 'green')
        plt.plot(T,r_p3, color = 'violet')
        plt.title("Région 2 Jeune : beta = "+ str(self.virus.beta)+"/ pi = "+ str(self.virus.pi) + "/ mu = "+ str(self.virus.mu))
        plt.subplot(2,2,4)
        plt.plot(T,s4, color = 'blue')
        plt.plot(T,u4, color = 'red')
        plt.plot(T,p4, color = 'yellow')
        plt.plot(T,r_u4, color = 'green')
        plt.plot(T,r_p4, color = 'violet')
        plt.title("Région 2 Vieux : beta = "+ str(self.virus.beta)+"/ pi = "+ str(self.virus.pi) + "/ mu = "+ str(self.virus.mu))

        plt.show()

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
    def evol_local_seule (self,proportion_de_la_pop_testee):
        # Condition pour pas faire des calculs inutiles ou incohérent
        if  (self.U0 < (self.population)*(10**-10)) or(self.U0 > self.population):

            self.history[0].append(self.S0)
            self.history[1].append(self.U0)
            self.history[2].append(self.P0)
            self.history[3].append(self.R0_U)
            self.history[4].append(self.R0_P)
        else : 
                # Calcul des variations 
            nmbr_test = proportion_de_la_pop_testee * self.population
            dS = - self.virus.beta* self.S0 * (self.U0 + (1-self.virus.pi)*self.P0) # *(self.S0/self.population) 
            dU = -dS - self.virus.mu*self.U0 - (nmbr_test)*(self.U0/(self.U0+self.R0_U+self.S0)) # (epsilon * NN )*(U0/(U0+R0+S0) le nombr de test * la proportion d'infectée dans la population qu'il reste à tester donc c'est le nombre de personne détectées positives à la fin des test
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
        ax.set(xlabel='Temps (en jours)', ylabel='Proportion de la population', title="Evolution discrète pour beta = "+ str(self.virus.beta)+ " / pi = " + str(self.virus.pi) + " / mu = " + str(self.virus.mu))
        plt.legend([S,U,P,RU,RP,R], ['S', 'U', 'P','R_U','R_P','R total'], loc='best')
        plt.show()
    
    def  crache_un_graphe_continu(self,proportion_test,time_limit) : 
            # définition du système d'équation diff 
        def evolution_continue(vecteur_condition_initiale,t) :# avec t = 1 dans notre cas
            population = self.population
            nmbr_test = proportion_test * population
            S,U,P,R_U,R_P = vecteur_condition_initiale
            dS = - self.virus.beta* S * (U + (1-self.virus.pi)*P) 
            dU = -dS - self.virus.mu*U - (nmbr_test)*(U/(U+R_U+S)) 
            dP = (nmbr_test )*((U/(U+R_U+S))) - self.virus.mu * P
            dR_U = self.virus.mu * U  # on ajoute les personnes guéries et qui vont se refaire testées
            dR_P = self.virus.mu * P # on ajoute les personnes guéries sans le savoir
            dS = (dP+dU+dR_U+dR_P) * (-1)
            return [dS,dU,dP,dR_U,dR_P]
        
        T = np.arange(0,time_limit,1) # liste de 0 à 99 
            # On extrait les données de la résolution 
        data = odeint (evolution_continue,[self.S0,self.U0,self.P0,self.R0_U,self.R0_P],T)
        data = data /self.population
        Susceptible = data[:,0] 
        Undetected = data[:,1]
        Positive = data[:,2]
        R_Undectected = data[:,3]
        R_Positive = data[:,4]
                # Création du Graphique 
        fig, ax = plt.subplots(nrows=1, ncols=1)
        S, = ax.plot(T, Susceptible, marker='+', color='blue', label='S0')
        U ,= ax.plot(T, Undetected , marker='+', color='red', label='U0')
        P ,= ax.plot(T, Positive, marker='+', color='yellow', label='P0')
        RU, = ax.plot(T, R_Undectected, marker='+', color='green', label='RO_U')
        RP ,= ax.plot(T, R_Positive, marker='+', color='violet', label='R0_P')
        R, = ax.plot(T, R_Undectected + R_Positive, marker='+', color='orange', label='RP+RU')
        ax.set(xlabel='Temps (en jours)', ylabel='Proportion de la population', title="Evolution continu pour beta = "+ str(self.virus.beta)+ " / pi = " + str(self.virus.pi) + " / mu = " + str(self.virus.mu))
        plt.legend([S,U,P,RU,RP,R], ['S', 'U', 'P','R_U','R_P','R total'], loc='best')
        plt.show()

    def evol_local_seule_continu (self,proportion_test) :
            # définition du système d'équation diff 
        def evolution_continue(vecteur_condition_initiale,t) :
            population = self.population
            nmbr_test = proportion_test * population
            S,U,P,R_U,R_P = vecteur_condition_initiale
            dS = - self.virus.beta* S * (U + (1-self.virus.pi)*P) 
            dU = -dS - self.virus.mu*U - (nmbr_test)*(U/(U+R_U+S)) 
            dP = (nmbr_test )*((U/(U+R_U+S))) - self.virus.mu * P
            dR_U = self.virus.mu * U  # on ajoute les personnes guéries et qui vont se refaire testées
            dR_P = self.virus.mu * P # on ajoute les personnes guéries sans le savoir
            dS = (dP+dU+dR_U+dR_P) * (-1)
            return [dS,dU,dP,dR_U,dR_P]
            #Extraction et actualisation des données 
        T = np.arange(0,100,1)
        data = odeint (evolution_continue,[self.S0,self.U0,self.P0,self.R0_U,self.R0_P],T)
        self.S0 = data[1,0] 
        self.U0 = data[1,1]
        self.P0 = data[1,2]
        self.R0_U = data[1,3]
        self.R0_P = data[1,4]
        self.history[0].append(self.S0)
        self.history[1].append(self.U0)
        self.history[2].append(self.P0)
        self.history[3].append(self.R0_U)
        self.history[4].append(self.R0_P)

    
    def evol_local_seule_continu_controle (self,proportion_test) :
            # définition du système d'équation diff 
        def evolution_continue(vecteur_condition_initiale,t) :
            population = self.population
            nmbr_test = proportion_test * population
            S,U,P,R_U,R_P = vecteur_condition_initiale
            dS = - self.virus.beta* S * (U + (1-self.virus.pi)*P) 
            dU = -dS - self.virus.mu*U - (nmbr_test)*(U/(U+R_U+S)) 
            dP = (nmbr_test )*((U/(U+R_U+S))) - self.virus.mu * P
            dR_U = self.virus.mu * U  # on ajoute les personnes guéries et qui vont se refaire testées
            dR_P = self.virus.mu * P # on ajoute les personnes guéries sans le savoir
            dS = (dP+dU+dR_U+dR_P) * (-1)
            return [dS,dU,dP,dR_U,dR_P]
            #Extraction et actualisation des données 
        T = np.arange(0,100,1)
        data = odeint (evolution_continue,[self.S0,self.U0,self.P0,self.R0_U,self.R0_P],T)
        self.S0 = data[1,0] 
        self.U0 = data[1,1]
        self.P0 = data[1,2]
        self.R0_U = data[1,3]
        self.R0_P = data[1,4]
        self.history[0].append(self.S0)
        self.history[1].append(self.U0)
        self.history[2].append(self.P0)
        self.history[3].append(self.R0_U)
        self.history[4].append(self.R0_P)
        return []

    def determiner_controllabilite (self,proportion_de_la_pop_testee,time_limit) :
            # renvoie un vecteur des info cool et les parmètres 
            #Ici on imule l'épidémie sur une durée donnée 
            for i in range (time_limit) :
                self.evol_local_seule_continu(proportion_de_la_pop_testee)   
            # Extraction des données de pic 
            data = self.history
            data_U = np.array(self.history[1])
            data_P = np.array(self.history [2])
            data_I = data_P + data_U
            data_I = list(data_I)
            pic_value = max(data_I)
            pic_time = data_I.index(pic_value)
            I_last_time = data[1][-1] + data [2][-1] 
            I_all_time = data[1][-1] + data [2][-1] + data [3][-1] + data [4][-1 ] # je somme au dernier temps les recovered les infecte à ce temps pour avoir le noimbre d'infécté jusque là par l'épidémie 
            return [self.virus.R,self.virus.mu,self.virus.pi,proportion_de_la_pop_testee,pic_value,pic_time,I_all_time,I_last_time]

