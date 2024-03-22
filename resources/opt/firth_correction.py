from sympy import *
import numpy
import time

alpha, beta, maf, g = symbols('alpha beta maf g')
p = exp(alpha + beta * g) / (1 + exp(alpha + beta * g))
e_0 = Sum(p * (1 - p) * g**0 * (1 - maf)**2, (g, 0, 3))
e_1 = Sum(p * (1 - p) * g**1 * (2 * maf * (1 - maf)), (g, 0, 3))
e_2 = Sum(p * (1 - p) * g**2 * maf**2, (g, 0, 3))
db_0 = Sum(p * (1 - p) * (1 - 2 * p) * (g**(0+1)) * (1 - maf)**2, (g, 0, 3))
db_1 = Sum(p * (1 - p) * (1 - 2 * p) * (g**(1+1)) * (2 * maf * (1 - maf)), (g, 0, 3))
db_2 = Sum(p * (1 - p) * (1 - 2 * p) * (g**(2+1)) * maf**2, (g, 0, 3))
da_0 = Sum(p * (1 - p) * (1 - 2 * p) * (g**(0)) * (1 - maf)**2, (g, 0, 3))
da_1 = Sum(p * (1 - p) * (1 - 2 * p) * (g**(1)) * (2 * maf * (1 - maf)), (g, 0, 3))
da_2 = Sum(p * (1 - p) * (1 - 2 * p) * (g**(2)) * maf**2, (g, 0, 3))

K_alpha = (da_0 * e_2 + da_2 * e_0 - 2 * da_1 * e_1) / (2 * (e_0 * e_2 - (e_1)**2))
K_beta = (db_0 * e_2 + db_2 * e_0 - 2 * db_1 * e_1) / (2 * (e_0 * e_2 - (e_1)**2))
K_alpha_wrt_alpha = diff(K_alpha, alpha)
K_alpha_wrt_beta = diff(K_alpha, beta)
K_beta_wrt_beta = diff(K_beta, beta)

Ka = lambdify([alpha, beta, maf], K_alpha, "numpy")
Kb = lambdify([alpha, beta, maf], K_beta, "numpy")
Kaa = lambdify([alpha, beta, maf], K_alpha_wrt_alpha, "numpy")
Kab = lambdify([alpha, beta, maf], K_alpha_wrt_beta, "numpy")
Kbb = lambdify([alpha, beta, maf], K_beta_wrt_beta, "numpy")


