#PROJET : SIMULATION DESCENTE D'UNE FUSEE
#Travail fait par :
#ELBAHRAOUI Imade
#KHADRAOUI Mohamed El Bachir
#CHPS

from math import exp

def fct_uE(V):
    return 0.001 * (V**2)

def fct_f(x : float, y : float):
    result = -0.5*exp(-(x**2 + y**2))
    return result

def fct_k1():
    return 0.168

def fct_alpha(e):
    k1 = 0.168
    return k1/e

def fct_beta():
    return 10**8

def fct_uD():
    return 20

def fct_k0():
    return 0.026

