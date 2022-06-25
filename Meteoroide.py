#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import sys
import os

epsilon = sys.float_info.epsilon

par1 = -6.65125
par2 = 1.39813e-06
grav = 980.0
km2cm = 100000

# ID | x_0 | y_0 | u_0 | v_0 | M_0| K_1 | K_2 | tau

Dados = [['M1', 0.0, 160.0, 20.0, -60.0, 1.00, 1.00, 1.0E-11, 0.02], 
	 ['M2', 0.0, 99.0, 20.0, -25.0, 0.010, 1.00, 1.0E-11, 0.02], 
         ['M3', 0.0, 90.0, 20.0, -20.0, 0.050, 1.00, 1.0E-11, 0.02], 
         ['M4', 0.0, 100.0 , 20.0, -25.0, 0.025 , 1.00, 1.0E-11, 0.02], 
         ['D6028', 0.0, 160.0, 20.0, -59.80, 0.48, 1.00, 1.0E-11, 0.02], 
	 ['H1139', 0.0, 99.0, 20.0, -26.69, 0.00810, 1.00, 1.0E-11, 0.02], 
         ['H1073', 0.0, 92.2, 20.0, -24.20, 0.01510, 1.00, 1.0E-11, 0.02], 
         ['H1134', 0.0, 96.5 , 20.0, -20.85, 0.01870 , 1.00, 1.0E-11, 0.02]] 
      

# compara float
def compara(a,b):    
    print ("epsilon {0} ".format(epsilon))
    return (abs(a-b) <= epsilon)

# Calcula a densidade atmosferica a uma altitude y(cm)
def dens_atmosf(y):
    return np.exp(par1 - (par2 * y))

# Calcula fx, fy, fu, fv, fm e s (velocidade) para a posicao (x,y)
# a velocidades (u,v) e a massa (m), com os parametros k1, k2 e tau
def dif_eq(x, y, u, v, m, k1, k2, tau):#, fx, fy, fu, fv, fm, s):
    fx = u
    fy = v
    s = np.sqrt(u*u + v*v)
    rho = dens_atmosf(y)
    
    fu = -k1 * rho * s * u * np.exp((-1/3) * np.log(m))
    fv = -k1 * rho * s * v * np.exp((-1/3) * np.log(m)) - grav
    fm = -k2 * rho * (s*s*s) * np.exp((2/3) * np.log(m))
    return (fx, fy, fu, fv, fm, s)


"""
# Main
"""
def main(argv = None):
    for i in range(0,8):
        fp = None

        fx = 0.0
        fy = 0.0
        fu = 0.0
        fv = 0.0
        fm = 0.0
        s = 0.0
        rho = 0.0

        fx1 = 0.0
        fy1 = 0.0
        fu1 = 0.0
        fv1 = 0.0
        fm1 = 0.0
        s1 = 0.0
        rho1 = 0.0

        t = 0.0
        dt = 0.0
        k = 1

        print("METEOROIDE")
        print("Condicoes iniciais")
        print("                          Meteoroide : {0}".format(Dados[i][0])) # ID
        print("                Posicao inicial (km) : {0}".format(Dados[i][1])) # x_0
        print("               Altitude inicial (km) : {0}".format(Dados[i][2])) # y_0
        print("Velocidade horizontal inicial (km/s) : {0}".format(Dados[i][3])) # u_0
        print("  Velocidade vertical inicial (km/s) : {0}".format(Dados[i][4])) # v_0
        print("               Massa inicial (grama) : {0}".format(Dados[i][5])) # M_0
        print("                        Parametro K1 : {0}".format(Dados[i][6]))
        print("                        Parametro K2 : {0}".format(Dados[i][7]))
        print("                       Parametro tau : {0}".format(Dados[i][8]))

        # cria nome do arquivo de valores calculados
        arq = Dados[i][0]+".dat"

        # Transforma dados de entrada x, y, u e v de km para cm
        # certifica que a velocidade e negativa.
        x = Dados[i][1] * km2cm
        y = Dados[i][2] * km2cm
        u = Dados[i][3] * km2cm

        v = -np.fabs(Dados[i][4] * km2cm)
        m = Dados[i][5]
        massa_inicial = m
        k1 = Dados[i][6]
        k2 = Dados[i][7]
        tau = Dados[i][8]

        # Cria arquivo de dados
        fp = open(arq,"w")

        print("  #    t          x          y          u          v          m      mag")
        fp.write("# id,t,x,y,u,v,massa,mag\n")
    
        #  This is the main loop of the program
        while (m >= (massa_inicial * 0.1)):
            # Os ifs a seguir selecionam o incremento adequado para o tempo dt
            if (m > (0.8 * massa_inicial)):
                dt = 0.1
            else:
                if(m > (0.5 * massa_inicial)):
                    dt = 0.05
                else:
                    if (m > (0.35 * massa_inicial)):
                        dt = 0.02
                    else:
                        dt = 0.01
            t = t + dt

            # Calcula 5 equacoes diferenciais (estado i)
            fx, fy, fu, fv, fm, s = dif_eq(x, y, u, v, m, k1, k2, tau) #, fx, fy, fu, fv, fm, s)

            # Calcula a predicao (estado i+1)
            x1 = x + dt * fx
            y1 = y + dt * fy
            u1 = u + dt * fu
            v1 = v + dt * fv
            m1 = m + dt * fm
            rho1 =  dens_atmosf(y1)
            s1 = np.sqrt((u1*u1) + (v1*v1))

            # Calcula 5 equacoes diferenciais (estado i+1)
            fx1, fy1, fu1, fv1, fm1, s1 = dif_eq(x1, y1, u1, v1, m1, k1, k2, tau)#, fx1, fy1, fu1, fv1, fm1, s1)

            # Calcula a predicao do estado i + 1
            x = x + 0.5 * dt * (fx + fx1)
            y = y + 0.5 * dt * (fy + fy1)
            u = u + 0.5 * dt * (fu + fu1)
            v = v + 0.5 * dt * (fv + fv1)
            m = m + 0.5 * dt * (fm + fm1)            
            
            # Calcula a magnitude aparente do meteoroide
            esp = -0.5 * tau * fm1 * ((u1*u1) + (v1*v1))
            mag = 5.0 * (np.log(y) / np.log(10)) - 2.5 * (np.log(esp) / np.log(10)) - 5.465 

            # Mostra resultados na tela
            # if (compara((y/km2cm),0.0) == False):
            if ((y) > 0.0):
                print("{0:3d} {1:4.2f} {2:10.4f} {3:10.4f} {4:10.5f} {5:10.5f} {6:10.7f} {7:8.2f}".format(k, t, x/km2cm, y/km2cm, u/km2cm, v/km2cm, m, mag))
                fp.write("{0:4.2f},{1:10.4f},{2:10.4f},{3:10.5f},{4:10.5f},{5:10.7f},{6:8.2f}\n".format(t, x/km2cm, y/km2cm, u/km2cm, v/km2cm, m, mag))
            else:
                print("Meteoroide atingiu o solo.")
                m = 0

            k+=1

        fp.close()

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
