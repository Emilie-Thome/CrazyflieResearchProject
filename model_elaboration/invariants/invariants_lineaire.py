import numpy as np
from numpy import linalg as LA

# Traitement des coefficients des invariants générés.

# [u1^2, u2^2, vz^2, z^2, u1u2, u1vz, u1z, u2vz, u2z, vzz]
linear1 = np.zeros((10,10))
# [u1, u2, vz, z]
linear2 = np.zeros((4,4))

linear1[0] = [-5.21132, 0.214624, -0.534143, 1.72909, \
-3.80898, 5.83535, 15.0162, 0.131419, -1.27244, -0.979337]
linear1[1] = [1.31988, -0.770285, 0.041812, -1.72909, \
-0.139783, -0.831863, -2.80116, 0.569732, 3.30339, -0.188717]
linear1[2] = [-5.94753, 0.222855, -0.344033, 1.72909, \
-3.92687, 4.01903, 14.9345, 0.246536, -1.2898, -1.18894]
linear1[3] = [1.1565, -0.741834, 0.0649169, -1.72909, \
-0.040769, -1.24237, -2.88279, 0.915108, 3.28602, -0.398323]
linear1[4] = [31.2101, 1.3432, 0.172869, 1.72909, \
-27.1931, -4.85767, 19.0546, 2.66012, -3.65302, -1.75886]
linear1[5] = [-0.220388, -0.123081, -0.129193, -1.72909, \
-0.329764, 0.351377, 1.23731, 0.259424, 0.922806, -0.968237]
linear1[6] = [141.33, 0.390524, 1.42941, 1.73762, \
-14.8584, -28.4266, 31.3419, 1.49428, -1.64752, -3.15199]
linear1[7] = [6.96035, 4.66558, 0.0211133, 1.73762, \
-11.3972, -0.766697, 6.95542, 0.627713, -5.69457, -0.383077]
linear1[8] = [0.252764, 0.128432, 0.0836221, 1.73762, \
0.360349, -0.290769, -1.32546, -0.207265, -0.944809, 0.762374]
linear1[9] = [0.19406, 0.119119, 0.201574, 1.73762, \
0.304081, -0.395563, -1.16138, -0.309912, -0.909911, 1.18366]

linear2[0] = [9.01861, -0.474074, -0.906985, 1]
linear2[1] = [2.00142, -1.63861, -0.11023, 1]
linear2[2] = [-0.381399, -0.271868, 0.219373, 1]
linear2[3] = [-0.334187, -0.261826, 0.340596, 1]

def print_invariants():
	def print_invariant1(l) :
		print(str(l[0]) + " u1^2 + " \
			+ str(l[1]) + " u2^2 + " \
			+ str(l[7]) + " u2 vz + " \
			+ str(l[2]) + " vz^2 + u1 (" \
			+ str(l[4]) + " u2 + " \
			+ str(l[5]) + " vz + " \
			+ str(l[6]) + " z) + " \
			+ str(l[8]) + " u2 z + " \
			+ str(l[9]) + " vz z + " \
			+ str(l[3]) + " z^2")


	def print_invariant2(l) :
		print(str(l[0]) + " u1^2 + " \
			+ str(l[1]) + " u2^2 + " \
			+ str(l[7]) + " u2 vz + " \
			+ str(l[2]) + " vz^2 + " \
			+ str(l[8]) + " u2 z + " \
			+ str(l[9]) + " vz z + " \
			+ str(l[3]) + " z^2 + u1 (" \
			+ str(l[4]) + " u2 + " \
			+ str(l[5]) + " vz + " \
			+ str(l[6]) + " z)")

	def print_invariant3(l) :
		print(str(l[0]) + " u1 + " \
			+ str(l[1]) + " u2 + " \
			+ str(l[2]) + " vz + " \
			+ str(l[3]) + " z)")

	print_invariant2(linear1[0])
	print_invariant1(linear1[1])
	print_invariant2(linear1[2])
	print_invariant1(linear1[3])
	print_invariant2(linear1[4])
	print_invariant2(linear1[5])
	print_invariant2(linear1[6])
	print_invariant2(linear1[7])
	print_invariant1(linear1[8])
	print_invariant1(linear1[9])
	print_invariant3(linear2[0])
	print_invariant3(linear2[1])
	print_invariant3(linear2[2])
	print_invariant3(linear2[3])



eigvals1, P1 = LA.eig(linear1)
eigvals2, P2 = LA.eig(linear2)

for eig in eigvals2 :
	if LA.norm(eig) < 0.1 :
		print(eig)