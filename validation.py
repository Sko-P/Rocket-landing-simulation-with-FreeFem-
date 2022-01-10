#PROJET : SIMULATION DESCENTE D'UNE FUSEE
#Travail fait par :
#ELBAHRAOUI Imade
#KHADRAOUI Mohamed El Bachir
#CHPS


#Nous vous invitons à faire varier la vitesse et l'épaisseur afin d'observer les résultats

from coeffelem import Solveur

if __name__ == "__main__":
    Vitesse = 3200
    Epaisseur = 1
    S = Solveur("fusee.msh", Vitesse, Epaisseur)
    S.assemblage_EF_P1()
