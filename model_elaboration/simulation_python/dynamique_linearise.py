import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

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

delta = (4*k_2)**2 - 16*k_1*(4*k_3-m*g)
pwm_e = (-4*k_2 + np.sqrt(delta))/(8*k_1)
u_1_e = (pwm_e-base)/K_i

alpha = 8*pwm_e*k_1/m + 4*k_2/m
A = np.array([[           0,          1,         0,             0],
			  [-2*K_p*alpha, -K_p*alpha, K_i*alpha, 0.5*K_p*alpha],
			  [          -2,         -1,         0,           0.5],
			  [          -1,          0,         0,             0]])

def f(Dx):
	"""
	Avec Dx un vecteur contenant :
	- Dz   : écart à l'équilibre de l'altitude du drône
	- Dv_z : écart à l'équilibre de la vitesse verticale du drône
	- Du1  : écart à l'équilibre de \int(-2z+0.5*\int(-z)-v_z)
	- Du2  : écart à l'équilibre de \int(-z)
	"""
	return A.dot(Dx)

def EulerMethod(f, x0, t0, tend, h):
    T = np.arange(t0, tend, h)
    X = [0 for i in T]
    x = x0
    i = 0
    for t in T :
    	X[i] = x
    	x = x + h*f(x)
    	i += 1
    X = np.array(X)

    for j in range(i) :
    	X[j,2] = X[j,2] + u_1_e

    return X,T

[X, T] = EulerMethod(f, np.array([1,0,0,0]), 0, 10, 0.0002)
plt.plot(T, X[:,0], label='z')
plt.plot(T, X[:,1], label='vz')
plt.plot(T, X[:,2], label='u1')
plt.plot(T, X[:,3], label='u2')
plt.xlabel('Temps [s]')
plt.legend()
plt.title('Contrôleur d\'altitude, système linéarisé')
plt.grid()
plt.show()