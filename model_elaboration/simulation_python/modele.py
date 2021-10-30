import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from control.matlab import *
from sympy import *
from sympy.abc import x, y, z, v, c


''' MODÉLISATION '''

g = 9.81
m = 0.028
k_1 = 2.130295*10**(-11)
k_2 = 1.032633*10**(-6)
k_3 = 5.484560*10**(-4)
scale = 1000
base = 36000
K_pz = 2
K_iz = 0.5
K_pv = 25
K_iv = 15
K_p = K_pv*scale
K_i = K_iv*scale


delta = (4*k_2)**2 - 4*(4*k_1)*(4*k_3-m*g)
pwm_e = (-4*k_2 + np.sqrt(delta))/(2*(4*k_1))
u_1_e = (pwm_e-base)/K_i
print("Le PWM d'équilibre est : " + str(pwm_e) + "\n")
print("Le u_1 d'équilibre est : " + str(u_1_e) + "\n")

alpha = 8*pwm_e*k_1/m + 4*k_2/m
A = np.array([[           0,          1,         0,             0],
			  [-2*K_p*alpha, -K_p*alpha, K_i*alpha, 0.5*K_p*alpha],
			  [          -2,         -1,         0,           0.5],
			  [          -1,          0,         0,             0]])


print("Nous avons une équation de la forme dx/dt = Ax avec A :")
print(str(A) + "\n")
lbd, B = LA.eig(A)
print("\nLes valeurs propres de A sont " + str(lbd))

Q = np.eye(4)
# X = lyap(A, Q) résout AX+XA^t=-Q mais je veux A^tP+PA=-Q
A_T = np.array([[0, 	 -2*K_p*alpha,	 -2,	-1],
				[1,       -K_p*alpha,	 -1,	 0],
				[0,        K_i*alpha,	  0,	 0],
				[0,		0.5*K_p*alpha,	0.5,	 0]])
P = lyap(A_T,Q)
# epsilon = 0.1
# perturbation = np.array([epsilon*np.random.uniform(-0.5, 0.5, 4),\
# 						 epsilon*np.random.uniform(-0.5, 0.5, 4),\
# 						 epsilon*np.random.uniform(-0.5, 0.5, 4),\
# 						 epsilon*np.random.uniform(-0.5, 0.5, 4)])
# P = lyap(A_T + perturbation,Q)
print("Résolution de l'équation de Lyapunov pour Q :")
print(Q)
print("\nSolution P de l'équation de Lyapunov :")
print(P)


''' DE CANDIDAT À FONCTION DE LYAPUNOV '''

lbd, B = LA.eig(P)
print("\nLes valeurs propres de P sont " + str(lbd)+"\n")




''' TESTE DE LA RÉSISTANCE AUX PERTURBATIONS '''

# epsilon = 0.1
# perturbation = np.array([epsilon*np.random.uniform(-0.5, 0.5, 4),\
# 						 epsilon*np.random.uniform(-0.5, 0.5, 4),\
# 						 epsilon*np.random.uniform(-0.5, 0.5, 4),\
# 						 epsilon*np.random.uniform(-0.5, 0.5, 4)])
# lbd, B = LA.eig(P+perturbation)
# print("\nAvec perturbations :")
# print(P + perturbation)
# print("\nLes valeurs propres de P avec perturbations sont " + str(lbd) + " et la base des vecteurs propres :")
# print(B)




''' TESTE DE LA SOLUTION DE L'ÉQUATION DE LYAPUNOV '''

# print("\nTeste de l'équation de Lyapunov (A^tP+PA = -Q) :")
# print((A_T).dot(P)+P.dot((A_T).T))
# lbd, B = LA.eig((A_T).dot(P)+P.dot((A_T).T))
# print("\nLes valeurs propres de A^tP+PA sont " + str(lbd))




''' CORRESPONDANCE ENTRE LES VARIABLES
z   -> z
vz  -> v
u1  -> x
u2  -> y
'''

''' LINÉARISÉ '''

dz  = Poly(v, z, v, x, y, domain='RR')
dvz = Poly(2*(-K_p*alpha)*z+(-K_p*alpha)*v+(K_i*alpha)*x-1/2*y, z, v, x, y, domain='RR')
du1 = Poly(-2*z-v+1/2*y, z, v, x, y, domain='RR')
du2 = Poly(-z, z, v, x, y, domain='RR')

lyap = Poly((2.01656529*z+0.08342414*v+0.11486412*x-1.29171207*y)*z+(0.08342414*z+0.07182875*v-0.08901516*x-0.04613902*y)*v+(0.11486412*z-0.08901516*v+1.35166103*x-0.56806025*y)*x+(-1.29171207*z-0.04613902*v-0.56806025*x+2.44792802*y)*y, z, v, x, y, domain='RR')
print("Polynome candidat de Lyapunov pour la dynamique linéarisée :")
print(lyap)
print()
# lyap_dz  = lyap.diff((0,1))
# lyap_dvz = lyap.diff((1,1))
# lyap_du1 = lyap.diff((2,1))
# lyap_du2 = lyap.diff((3,1))
# lyap_dt = (((dz.mul(lyap_dz)).add(dvz.mul(lyap_dvz))).add(du1.mul(lyap_du1))).add(du2.mul(lyap_du2))
# print(lyap_dt)



''' NON LINÉARISÉ '''

dz  = Poly(v, z, v, x, y, domain='RR')
dvz = Poly(4/m*(k_1*(K_p*(-2*z-v+1/2*y) + K_i*x + base)**2 + k_2*(K_p*(-2*z-v+1/2*y) + K_i*x + base) + k_3) - g, z, v, x, y, domain='RR')
du1 = Poly(-2*z-v+1/2*y, z, v, x, y, domain='RR')
du2 = Poly(-z, z, v, x, y, domain='RR')

lyap = Poly((2.01656529*z+0.08342414*v+0.11486412*(x-u_1_e)-1.29171207*y)*z+(0.08342414*z+0.07182875*v-0.08901516*(x-u_1_e)-0.04613902*y)*v+(0.11486412*z-0.08901516*v+1.35166103*(x-u_1_e)-0.56806025*y)*(x-u_1_e)+(-1.29171207*z-0.04613902*v-0.56806025*(x-u_1_e)+2.44792802*y)*y, z, v, x, y, domain='RR')
print("Polynome candidat de Lyapunov pour la dynamique non linéarisée :")
print(lyap)
# lyap_dz  = lyap.diff((0,1))
# lyap_dvz = lyap.diff((1,1))
# lyap_du1 = lyap.diff((2,1))
# lyap_du2 = lyap.diff((3,1))
# lyap_dt = (((dz.mul(lyap_dz)).add(dvz.mul(lyap_dvz))).add(du1.mul(lyap_du1))).add(du2.mul(lyap_du2))
# print(lyap_dt)