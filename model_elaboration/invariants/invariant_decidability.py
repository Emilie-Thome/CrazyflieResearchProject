import numpy as np
from sympy import *
from sympy.abc import x, y, z, v, a, b



''' CONSTANTES DE LA DYNAMIQUE '''
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

delta = (4*k_2)**2-16*k_1*(k_3-m*g)
pwm_e = (-4*k_2+np.sqrt(delta))/(8*k_1)


''' RAPPELS DE LA DYNAMIQUE LINÉAIRE
alpha = 8*pwm_e*k_1/m+4*k_2/m
A = np.array([[           0,          1,         0,             0],
			  [-2*K_p*alpha, -K_p*alpha, K_i*alpha, 0.5*K_p*alpha],
			  [          -2,         -1,         0,           0.5],
			  [          -1,          0,         0,             0]])
'''


''' CORRESPONDANCE ENTRE LES VARIABLES
z   -> z
vz  -> v
u1  -> x
u2  -> y
dz   -> k
dvz  -> l
du1  -> m
du2  -> n
'''

c = (8*k_1*pwm_e+4*k_2)/m  # = alpha
a = -K_p*c
b = K_i*c

dz  = Poly(v, z, v, x, y, domain='QQ[a, b]')
dvz = Poly(2*a*z+a*v+b*x-1/2*y, z, v, x, y, domain='QQ[a, b]')
du1 = Poly(-2*z-v+1/2*y, z, v, x, y, domain='QQ[a, b]')
du2 = Poly(-z, z, v, x, y, domain='QQ[a, b]')


inv = Poly(	b**2*(8*b**3*x**4-12*a*b**2*x**3*y-16*b**3*x**3*y+\
			6*a**2*b*x**2*y**2+16*a*b**2*x**2*y**2+4*b**3*x**2*y**2-\
			a**3*x*y**3-4*a**2*b*x*y**3-2*a*b**2*x*y**3+b**2*y**4+\
			8*a*b**2*x**3*v-8*a**2*b*x**2*y*v-16*b**2*x**2*y*v-\
			16*a*b**2*x**2*y*v+2*a**3*x*y**2*v+10*a*b*x*y**2*v+\
			8*a**2*b*x*y**2*v+24*b**2*x*y**2*v+4*a*b**2*x*y**2*v-\
			a**2*y**3*v-4*a*b*y**3*v-4*b**2*y**3*v-16*a*b*x**2*v**2+\
			8*b**2*x**2*v**2+8*a**2*x*y*v**2+16*a*b*x*y*v**2-\
			16*b**2*x*y*v**2+2*a**2*y**2*v**2+4*b*y**2*v**2+\
			8*a*b*y**2*v**2+4*b**2*y**2*v**2+4*a*x*v**3-16*b*x*v**3+\
			8*a*y*v**3+24*b*y*v**3+4*v**4+32*a*b**2*x**3*z+\
			16*b**3*x**3*z-32*a**2*b*x**2*y*z-72*a*b**2*x**2*y*z-\
			16*b**3*x**2*y*z+8*a**3*x*y**2*z+32*a**2*b*x*y**2*z+\
			8*b**2*x*y**2*z+16*a*b**2*x*y**2*z-a*b*y**3*z-\
			12*b**2*y**3*z-12*a*b*x**2*v*z+16*a**2*b*x**2*v*z+\
			48*b**2*x**2*v*z+16*a*b**2*x**2*v*z+6*a**2*x*y*v*z-\
			8*a**3*x*y*v*z-32*a*b*x*y*v*z-32*a**2*b*x*y*v*z-\
			112*b**2*x*y*v*z-16*a*b**2*x*y*v*z+8*a**2*y**2*v*z+\
			34*a*b*y**2*v*z+32*b**2*y**2*v*z-8*a**2*x*v**2*z-\
			16*b*x*v**2*z+16*b**2*x*v**2*z+6*a*y*v**2*z-\
			8*a**2*y*v**2*z+8*b*y*v**2*z-32*a*b*y*v**2*z-\
			16*b**2*y*v**2*z-8*a*v**3*z-16*b*v**3*z+36*a**2*b*x**2*z**2+\
			8*b**2*x**2*z**2+48*a*b**2*x**2*z**2+8*b**3*x**2*z**2-\
			18*a**3*x*y*z**2+2*a*b*x*y*z**2-72*a**2*b*x*y*z**2-\
			40*b**2*x*y*z**2-36*a*b**2*x*y*z**2+8*a*b*y**2*z**2+\
			50*b**2*y**2*z**2-8*a**2*x*v*z**2+4*a**3*x*v*z**2+\
			32*a*b*x*v*z**2+16*a**2*b*x*v*z**2+80*b**2*x*v*z**2+\
			8*a*b**2*x*v*z**2-18*a**2*y*v*z**2+8*b*y*v*z**2-\
			80*a*b*y*v*z**2-72*b**2*y*v*z**2-8*a*v**2*z**2+\
			4*a**2*v**2*z**2-12*b*v**2*z**2+16*a*b*v**2*z**2+\
			8*b**2*v**2*z**2+2*a**2*x*z**3+8*a**3*x*z**3+32*a**2*b*x*z**3+\
			40*b**2*x*z**3+16*a*b**2*x*z**3-18*a*b*y*z**3-80*b**2*y*z**3+\
			2*a*v*z**3+8*a**2*v*z**3-8*b*v*z**3+36*a*b*v*z**3+\
			32*b**2*v*z**3+2*b*z**4+8*a*b*z**4+32*b**2*z**4), z, v, x, y, domain='QQ[a, b]')


epsilon = 10**(-10)

def not_in_differential_order(p, diff_order):
	if diff_order == [] :
		return True
	(q, r) = pdiv(p, diff_order[0])
	print(r)
	if r == 0 :
		return False
	else :
		return not_in_differential_order(r, diff_order[1:])


def differential_order(p):
	lie_derivatives = []
	p_to_derivate = p

	Np = 0
	while not_in_differential_order(p_to_derivate, lie_derivatives) :
		print(Np)
		lie_derivatives.append(p_to_derivate)

		p_dz  = p_to_derivate.diff((0,1))
		p_dvz = p_to_derivate.diff((1,1))
		p_du1 = p_to_derivate.diff((2,1))
		p_du2 = p_to_derivate.diff((3,1))
		Np += 1

		p_to_derivate = (((dz.mul(p_dz)).add(dvz.mul(p_dvz))).add(du1.mul(p_du1))).add(du2.mul(p_du2))
	
	return lie_derivatives, Np

def invariance_decidability(inv):
	lie_derivatives, Np = differential_order(inv)
	for i in range(Np-1):
		(q, r) = pdiv(lie_derivatives[i+1], lie_derivatives[i])
		for j in range(1, i):
			(q, r) = pdiv(r, lie_derivatives[j])
		if r == nul :
			roots = nroots(lie_derivatives[0], n=30)
			ans = True
			for j in range(1, i):
				ans = ans and (np.abs(lie_derivatives[j].eval(roots))<epsilon)
			if ans :
				return True
	return False

x=0.1
y=0.1
v=0.1
print(nroots(Poly(	b**2*(8*b**3*x**4-12*a*b**2*x**3*y-16*b**3*x**3*y+\
			6*a**2*b*x**2*y**2+16*a*b**2*x**2*y**2+4*b**3*x**2*y**2-\
			a**3*x*y**3-4*a**2*b*x*y**3-2*a*b**2*x*y**3+b**2*y**4+\
			8*a*b**2*x**3*v-8*a**2*b*x**2*y*v-16*b**2*x**2*y*v-\
			16*a*b**2*x**2*y*v+2*a**3*x*y**2*v+10*a*b*x*y**2*v+\
			8*a**2*b*x*y**2*v+24*b**2*x*y**2*v+4*a*b**2*x*y**2*v-\
			a**2*y**3*v-4*a*b*y**3*v-4*b**2*y**3*v-16*a*b*x**2*v**2+\
			8*b**2*x**2*v**2+8*a**2*x*y*v**2+16*a*b*x*y*v**2-\
			16*b**2*x*y*v**2+2*a**2*y**2*v**2+4*b*y**2*v**2+\
			8*a*b*y**2*v**2+4*b**2*y**2*v**2+4*a*x*v**3-16*b*x*v**3+\
			8*a*y*v**3+24*b*y*v**3+4*v**4+32*a*b**2*x**3*z+\
			16*b**3*x**3*z-32*a**2*b*x**2*y*z-72*a*b**2*x**2*y*z-\
			16*b**3*x**2*y*z+8*a**3*x*y**2*z+32*a**2*b*x*y**2*z+\
			8*b**2*x*y**2*z+16*a*b**2*x*y**2*z-a*b*y**3*z-\
			12*b**2*y**3*z-12*a*b*x**2*v*z+16*a**2*b*x**2*v*z+\
			48*b**2*x**2*v*z+16*a*b**2*x**2*v*z+6*a**2*x*y*v*z-\
			8*a**3*x*y*v*z-32*a*b*x*y*v*z-32*a**2*b*x*y*v*z-\
			112*b**2*x*y*v*z-16*a*b**2*x*y*v*z+8*a**2*y**2*v*z+\
			34*a*b*y**2*v*z+32*b**2*y**2*v*z-8*a**2*x*v**2*z-\
			16*b*x*v**2*z+16*b**2*x*v**2*z+6*a*y*v**2*z-\
			8*a**2*y*v**2*z+8*b*y*v**2*z-32*a*b*y*v**2*z-\
			16*b**2*y*v**2*z-8*a*v**3*z-16*b*v**3*z+36*a**2*b*x**2*z**2+\
			8*b**2*x**2*z**2+48*a*b**2*x**2*z**2+8*b**3*x**2*z**2-\
			18*a**3*x*y*z**2+2*a*b*x*y*z**2-72*a**2*b*x*y*z**2-\
			40*b**2*x*y*z**2-36*a*b**2*x*y*z**2+8*a*b*y**2*z**2+\
			50*b**2*y**2*z**2-8*a**2*x*v*z**2+4*a**3*x*v*z**2+\
			32*a*b*x*v*z**2+16*a**2*b*x*v*z**2+80*b**2*x*v*z**2+\
			8*a*b**2*x*v*z**2-18*a**2*y*v*z**2+8*b*y*v*z**2-\
			80*a*b*y*v*z**2-72*b**2*y*v*z**2-8*a*v**2*z**2+\
			4*a**2*v**2*z**2-12*b*v**2*z**2+16*a*b*v**2*z**2+\
			8*b**2*v**2*z**2+2*a**2*x*z**3+8*a**3*x*z**3+32*a**2*b*x*z**3+\
			40*b**2*x*z**3+16*a*b**2*x*z**3-18*a*b*y*z**3-80*b**2*y*z**3+\
			2*a*v*z**3+8*a**2*v*z**3-8*b*v*z**3+36*a*b*v*z**3+\
			32*b**2*v*z**3+2*b*z**4+8*a*b*z**4+32*b**2*z**4), z, domain = 'QQ[a,b]'), n=30))
# print(invariance_decidability(inv))