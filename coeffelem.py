#PROJET : SIMULATION DESCENTE D'UNE FUSEE
#Travail fait par :
#ELBAHRAOUI Imade
#KHADRAOUI Mohamed El Bachir
#CHPS

import numpy as np
from fonctions import *
from fichier import *





class Solveur:
    def __init__(self, NomduFichier, V, e):
        self.NomduFichier = NomduFichier
        self.V = V
        self.e = e


    def coeffelem_P1_rigid(self,l): 
        nbn, nbe, nba, coord, tri, ar, refn, reft, refa = lit_fichier_msh(self.NomduFichier)
        
        M1, M2, M3 = [int(j) for j in tri[l]] #Les indices de chaque triangle

        M1 = [coord[M1,0],coord[M1,1]]
        M2 = [coord[M2,0],coord[M2,1]]
        M3 = [coord[M3,0],coord[M3,1]]
        
        K = np.zeros((3,3))
        mesTl = 1/2*abs((M2[0]-M1[0])*(M3[1]-M2[1])-(M3[0]-M2[0])*(M2[1]-M1[1])) 

        K[0,0] = (fct_k1()/(4*mesTl)) * ((M2[0]-M3[0])**2 + (M2[1]-M3[1])**2)
        K[1,1] = (fct_k1()/(4*mesTl)) * ((M3[0]-M1[0])**2 + (M3[1]-M1[1])**2)
        K[2,2] = (fct_k1()/(4*mesTl)) * ((M1[0]-M2[0])**2 + (M1[1]-M2[1])**2)

        K[0,1] = (fct_k1()/(4*mesTl)) * (-(M1[0]-M3[0])*(M2[0]-M3[0])) - ((M1[1]-M3[1])*(M2[1]-M3[1]))
        K[1,0] = K[0,1]

        K[0,2] = (fct_k1()/(4*mesTl)) * (-(M3[0]-M2[0])*(M1[0]-M2[0])) - ((M3[1]-M2[1])*(M1[1]-M2[1]))
        K[2,0] = K[0,2]

        K[1,2] = (fct_k1()/(4*mesTl)) * (-(M2[0]-M1[0])*(M3[0]-M1[0])) - ((M2[1]-M1[1])*(M3[1]-M1[1]))
        K[2,1] = K[1,2]

        #print("k",l,"=\n",K)
        return K

    def coeffelem_P1_source(self,l): 
        nbn, nbe, nba, coord, tri, ar, refn, reft, refa = lit_fichier_msh(self.NomduFichier)
        M1, M2, M3 = [int(j) for j in tri[l]] #Les indices de chaque triangle

        M1 = [coord[M1,0],coord[M1,1]]
        M2 = [coord[M2,0],coord[M2,1]]
        M3 = [coord[M3,0],coord[M3,1]]

        mesTl = 1/2*abs(((M2[0]-M1[0])*(M3[1]-M2[1]))-((M3[0]-M2[0])*(M2[1]-M1[1]))) 
        
        val = mesTl / 3

        ff = fct_f( ((M1[0] + M2[0] + M3[0]) / 3), ((M1[1] + M2[1] + M3[1]) / 3) )

        mat = [1,1,1]
        f = []

        for i in range(0,3):
            f.append(val * ff * mat[i])

        #print("f",l,"=\n",f)
        return f 
            
    ##
    # Coefficients des termes de bord
    ## 

    def coeffelem_P1_poids_fourier(self,alpha): 

        nbn, nbe, nba, coord, tri, ar, refn, reft, refa = lit_fichier_msh(self.NomduFichier)

        ai1, ai2 = [int(j) for j in ar[alpha]] #Les indices de chaque sommets de l'arete

        a1 = [coord[ai1,0], coord[ai1,1]]

        a2 = [coord[ai2,0], coord[ai2,1]]

        mes = sqrt(abs((a2[0] - a1[0])**2 - (a2[1] - a1[1])**2))

        val = mes/6


        mat = [[2,1],
                [1,2]]

        p = [[0,0],[0,0]]


        for i in range(0,2):
            for j in range(0,2):

                p[i][j] = (val * fct_alpha(self.e) * mat[i][j])

        
        return p


    def coeffelem_P1_transf_fourier(self,alpha): 
        nbn, nbe, nba, coord, tri, ar, refn, reft, refa = lit_fichier_msh(self.NomduFichier)

        ai1, ai2 = [int(j) for j in ar[alpha]] #Les indices de chaque sommets de l'arete

        a1 = [coord[ai1,0], coord[ai1,1]]

        a2 = [coord[ai2,0], coord[ai2,1]]

        mes = sqrt(abs((a2[0] - a1[0])**2 - (a2[1] - a1[1])**2))

        val = mes/2

        mat = [1,1]
        p =[]

        

        for i in range(0,2):

            p.append(val * fct_alpha(self.e) * fct_uE(self.V) * mat[i])

        #print("ea=\n",p,"\n")
        return p


    def coeffelem_P1_poids_dirichlet(self,alpha): 

        nbn, nbe, nba, coord, tri, ar, refn, reft, refa = lit_fichier_msh(self.NomduFichier)

        ai1, ai2 = [int(j) for j in ar[alpha]] #Les indices de chaque sommets de l'arete

        a1 = [coord[ai1,0], coord[ai1,1]]

        a2 = [coord[ai2,0], coord[ai2,1]]

        mes = sqrt(abs((a2[0] - a1[0])**2 - (a2[1] - a1[1])**2))

        val = mes/6


        mat = [[2,1],
                [1,2]]

        p = [[0,0],[0,0]]


        for i in range(0,2):
            for j in range(0,2):

                p[i][j] = (val * fct_beta() * mat[i][j])

        
        return p

    def coeffelem_P1_transf_dirichlet(self,alpha): 
        nbn, nbe, nba, coord, tri, ar, refn, reft, refa = lit_fichier_msh(self.NomduFichier)

        ai1, ai2 = [int(j) for j in ar[alpha]] #Les indices de chaque sommets de l'arete

        a1 = [coord[ai1,0], coord[ai1,1]]

        a2 = [coord[ai2,0], coord[ai2,1]]

        mes = sqrt(abs((a2[0] - a1[0])**2 - (a2[1] - a1[1])**2))

        val = mes/2

        mat = [1,1]
        p =[]

        

        for i in range(0,2):

            p.append(val * fct_beta() * fct_uD() * mat[i])

        
        return p


    def assemblage_EF_P1(self):
        
        nbn, nbe, nba, coord, tri, ar, refn, reft, refa = lit_fichier_msh(self.NomduFichier)

        A = np.zeros((nbn,nbn)) #Matrice de rigidité

        f = np.zeros(nbn)

        print("nbn, nba, nbe : ",nbn,nba,nbe)
        

        for l in range(0,nbe):
            
            k = self.coeffelem_P1_rigid(l)
            ff = self.coeffelem_P1_source(l)

            i1 = tri[l][0]
            i2 = tri[l][1] 
            i3 = tri[l][2]

            
            A[i1][i1] += k[0][0]
            A[i1][i2] += k[0][1]
            A[i1][i3] += k[0][2]
            f[i1] += ff[0]

            A[i2][i1] += k[1][0]
            A[i2][i2] += k[1][1]
            A[i2][i3] += k[1][2]
            f[i2] += ff[1]
            
            A[i3][i1] += k[2][0]
            A[i3][i2] += k[2][1]
            A[i3][i3] += k[2][2]
            f[i3] += ff[2]

        K = np.copy(A)
        
        for i in range(0,len(ar)):
            
            if(int(refa[i]) == 2):
                
                p = self.coeffelem_P1_poids_fourier(i)
                e = self.coeffelem_P1_transf_fourier(i)

                i1 = ar[i][0]
                i2 = ar[i][1]
                
                A[i1][i1] += p[0][0]
                A[i1][i2] += p[0][1]
                f[i1] += e[0]

                A[i2][i1] += p[1][0]  
                A[i2][i2] += p[1][1]
                f[i2] += e[1]

            if(int(refa[i]) == 1):
                
                p = self.coeffelem_P1_poids_dirichlet(i)
                e = self.coeffelem_P1_transf_dirichlet(i)

                i1 = ar[i][0]
                i2 = ar[i][1]
                
                A[i1][i1] += p[0][0]
                A[i1][i2] += p[0][1]
                f[i1] += e[0]

                A[i2][i1] += p[1][0]  
                A[i2][i2] += p[1][1]
                f[i2] += e[1]

        print("================= Validation =================== * Resultats sur",self.NomduFichier)
        print("A=\n",A,"\n")
        print("F=\n",f,"\n")

        Uh = np.linalg.solve(A,f)
       
        print("Vitesse en m/s : \n", self.V)
        print("Epaisseur en cm :\n", self.e*100)
        print("min Uh =\n",min(Uh))
        print("max Uh =\n",max(Uh))
        
        





    