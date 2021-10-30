import numpy as np
import matplotlib.pyplot as plt
from control.matlab import *

Kp = 5
Kd = 1.25
g = 9.81
m = 0.33

A = np.array([[0, 1], [-Kp/m, -Kd/m]])
Q = np.array([[1, 0], [0, 1]])
P = lyap(A, Q)

print("Soit x = [z v_z]t.")
print("dx/dt = A*x, avec A = [[0, 1], [-Kp/m, -Kd/m]].\n")
print("On cherche à résoudre l'équation de Lyapunov : AtP + PA + Q = 0.")
print("Si P > 0 et Q > 0 alors le système est asymptoticalement stable.\n")

print("Pour Q = I2 > 0, on a :")
print("P[0,0] = " + str(P[0,0]) + " > 0")
print("P[0,0]*P[1,1]-P[0,1]² = " + str(P[0,0]*P[1,1]-P[0,1]*P[0,1]) + " > 0\n")

print("Nous avons donc P > 0 et donc le système est asymptoticalement stable.\n")

# Transfer Function
num = [Kd, Kp]
den = [m, Kd, Kp]
sys = TransferFunction(num,den)
print("La fonction de transfert du système est :")
print(sys)

plt.figure(1)
bode(tf(sys))

plt.figure(2)
rlocus(tf(sys))

y,t = step(sys, T=10)
plt.figure(3)
plt.plot(t,y)
plt.xlabel('Temps')
plt.ylabel('Réponse (z)')
plt.show()