#PROJET : SIMULATION DESCENTE D'UNE FUSEE
#Travail fait par :
#ELBAHRAOUI Imade
#KHADRAOUI Mohamed El Bachir
#CHPS

from math import *
from pylab import *

def lit_fichier_msh(fichier_msh):

    file  = open(fichier_msh,'r')
    
    nbs = [int(i) for i in file.readline().split()]
    #nbn, nbe, nba


    m_coord = zeros((nbs[0], 2),double)
    m_tri = zeros((nbs[1], 3),int)
    m_ar = zeros((nbs[2], 2),int)
    
    refn = []
    reft = []
    refa = []
    

    #Remplissage matrice de coordonées
    for i in range(0,nbs[0]):
        temp = file.readline().split()
        
        m_coord [i][0] = temp[0]
        m_coord [i][1] = temp[1]
        refn.append(temp[2])

    

    #Remplissage matrice de triangles
    for i in range(0,nbs[1]):
        temp = file.readline().split()
        
        m_tri [i][0] = double(temp[0]) -1
        m_tri [i][1] = double(temp[1]) -1
        m_tri [i][2] = double(temp[2]) -1
        reft.append(temp[3])

    

    #Remplissage de la matrice arêtes
    for i in range(0,nbs[2]):
        temp = file.readline().split()
        
        m_ar [i][0] = double(temp[0]) -1
        m_ar [i][1] = double(temp[1]) -1
        
        refa.append(temp[2])


    #Attribution des nbn, nbe, nba afin de pouvoir les retrouner
    nbn, nbe, nba = [i for i in nbs]
    
    return nbn, nbe, nba, m_coord, m_tri, m_ar, refn, reft, refa


def trace_maillage_ind(nbn, nbe, nba, m_coord, m_tri, m_ar) :

    for i in range(nbn):
        text(m_coord[i][0] - 0.03, m_coord[i][1] + 0.01, str(i), color='black')

    for i in range(nbe):
        x, y, z = [int(j) for j in m_tri[i]]
        text((m_coord[x, 0] + m_coord[y, 0] + m_coord[z, 0]) / 3, (m_coord[x, 1] + m_coord[y, 1] + m_coord[z, 1]) / 3, str(i),
             color='purple')

    for i in range(nba):
        x, y = [int(j) for j in m_ar[i]]
        text((m_coord[x, 0] + m_coord[y, 0]) / 2, (m_coord[x, 1] + m_coord[y, 1]) / 2, str(i), color='cyan')


def trace_maillage_ref(nbn, nbe, nba, coord, tri, ar, refn, reft, refa):
    for i in range(nbn):
        text(coord[i][0] - 0.03, coord[i][1] - 0.04, str(refn[i]), color='red')

    for i in range(nbe):
        x, y, z = [int(j) for j in tri[i]]
        text((coord[x, 0] + coord[y, 0] + coord[z, 0]) / 3, (coord[x, 1] + coord[y, 1] + coord[z, 1]) / 3 - 0.04,
             str(reft[i]), color='green')

    for i in range(nba):
        x, y = [int(j) for j in ar[i]]
        text((coord[x, 0] + coord[y, 0]) / 2, (coord[x, 1] + coord[y, 1]) / 2 - 0.04, str(refa[i]), color='orange')


def charge_et_affiche_maillage(FichierMaillage) : 
    fichier_msh = FichierMaillage
    nbn, nbe, nba, m_coord, m_tri, m_ar, refn, reft, refa = lit_fichier_msh(fichier_msh)
    #trace_maillage_ind(nbn, nbe, nba, m_coord, m_tri, m_ar)
    #trace_maillage_ref(nbn, nbe, nba, m_coord, m_tri, m_ar, refn, reft, refa)
    triplot(m_coord[:,0], m_coord[:,1], m_tri)
    pas_qualite_maillage(m_coord, m_tri, nbe,0)
    text(0.0,0,"Qualité : " + pas_qualite_maillage(m_coord, m_tri, nbe,0),bbox=dict(alpha=0.8),fontsize=8)
    text(0.0,2.0,"Pas : " + pas_qualite_maillage(m_coord, m_tri, nbe,1),bbox=dict(alpha=0.8),fontsize=8)
    show()



def pas_qualite_maillage(m_coord, m_tri, nbe,case): 
    
    
    pas_max = 0
    q_max = 0

    #Pour chaque triangle
    for i in range(nbe):
        
        s1 = m_tri[i][0]
        s2 = m_tri[i][1]
        s3 = m_tri[i][2]

        
        l1 = sqrt((m_coord[s2][0] - m_coord[s1][0])**2 +(m_coord[s2][1] - m_coord[s1][1])**2)
        l2 = sqrt((m_coord[s3][0] - m_coord[s2][0])**2 +(m_coord[s3][1] - m_coord[s2][1])**2)
        l3 = sqrt((m_coord[s3][0] - m_coord[s1][0])**2 +(m_coord[s3][1] - m_coord[s1][1])**2)
        
        current_p = max(l1,l2,l3)
        if(current_p > pas_max) :
            pas_max = current_p  #pas_max est le pas du maillage
        
        base = 0.5*(l1+l2+l3)
        mes = abs(0.5*((m_coord[s2][0] - m_coord[s1][0]) * (m_coord[s3][1]- m_coord[s2][1]) - (m_coord[s3][0] - m_coord[s2][0]) * (m_coord[s2][0] - m_coord[s1][0]) ))

        rayon = mes / base

        current_q = sqrt(3)/6 * (current_p/rayon)
        if(current_q > q_max):
            q_max = current_q #q_max est la qualité du maillage


    print(pas_max)
    print(q_max)

    if(case == 0):
        return str(q_max)
    else :
        return str(pas_max)

if __name__ == "__main__" :
    FichierMaillage = "./Maillages/m1.msh"
    charge_et_affiche_maillage(FichierMaillage)