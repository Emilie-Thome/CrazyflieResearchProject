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


def f(x, pwm):
	"""
	Avec x un vecteur contenant :
	- z   : l'altitude du drône
	- v_z : la vitesse verticale du drône
	- u1  : correspondant à \int(-2z+0.5*\int(-z)-v_z)
	- u2  : correspondant à \int(-z)
	"""
	dz = x[1]
	dv_z = 4/m*(k_1*pwm**2 + k_2*pwm + k_3) - g
	du_1 = -2*x[0] + 0.5*x[3] -x[1]
	du_2 = -x[0]
	
	dx = np.array([dz, dv_z, du_1, du_2])

	return dx

def EulerMethod(f, x0, t0, tend, h):
    T = np.arange(t0, tend, h)
    X = [0 for i in T]
    x = x0
    pwm = 0

    i = 0
    for t in T :
        X[i] = x
        dx = f(x, pwm)
        x = x + h*dx
        pwm = K_p*dx[2] + K_i*x[2] + base
        i += 1

    X = np.array(X)
    return X,T

[X, T] = EulerMethod(f, np.array([1,0,u_1_e,0]), 0, 10, 0.002)
plt.plot(T, X[:,0], label='z')
plt.plot(T, X[:,1], label='vz')
plt.plot(T, X[:,2], label='u1')
plt.plot(T, X[:,3], label='u2')
plt.xlabel('Temps [s]')
plt.legend()
plt.title('Contrôleur d\'altitude, système non linéarisé')
plt.grid()
plt.show()