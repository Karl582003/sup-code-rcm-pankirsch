#En este codigo, dado un punto, evalua y determina valores propios y vectores propios
import numpy as np
import math

#Parámetros del modelo
c = 0.025
mu_2 = 0.03
p_1 = 0
a = 1
r_2 = 0.18
mu_3 = 10
p_2 = 5
g_1 = 2*10**7
g_2 = 10**5
g_3 = 10**3
b = 1*10**(-9)
s_1 = 0
s_2 = 0
k_1 = c/mu_2
k_2 = math.sqrt(a**2 * k_1**2 + 2*a*b *g_2 * k_1 * r_2 - 2*a*k_1*r_2 +b**2 * g_2**2 *r_2**2 + 2*b*g_2*r_2**2 + r_2**2)

#Punto de equilibrio 2
x_2 = -( k_1 * (k_2 - r_2 + a*k_1 + b*g_2 *r_2) )/(2*b*r_2)
y_2 =  -(k_2 - r_2 + a*k_1 + b*g_2 * r_2 )/(2*b*r_2)
z_2 = ( k_1 * p_2 * (k_2 - r_2 + a*k_1 + b*g_2 *r_2)**2   ) / (4*b**2 * mu_3 * r_2**2 * (g_3 - (k_2 - r_2 + a*k_1 +b*g_2 * r_2)/(2*b*r_2)  ) )
#print(x_2, y_2, z_2)

#Punto de equilibrio 3
x3 = (k_1 * (k_2 + r_2 - a * k_1 - b * g_2 * r_2)) / (2 * b * r_2)
y3 = (k_2 + r_2 - a * k_1 - b * g_2 * r_2) / (2 * b * r_2)
z3 = (k_1 * p_2 * (k_2 + r_2 - a * k_1 - b * g_2 * r_2)**2) / (4 * b**2 * mu_3 * r_2**2 * (g_3 + (k_2 + r_2 - a * k_1 - b * g_2 * r_2) / (2 * b * r_2)))


def ecdif(x,y,z):
    dx = c*y - mu_2 * x
    dy = r_2*y*(1-b*y) - (a*x*y)/(g_2 + y)
    dz = (p_2 * x * y) / (g_3 + y ) - mu_3 * z
    return (dx,dy,dz)

#print(x3, y3, z3)
print("Punto de equilibrio analizado:")
print('x =',x3)
print('y =',y3)
print('z =',z3)

def matriz(x,y,z):
    a_21 = - (a*y)/(g_2 + y)
    a_22 = r_2 * (1-2*b*y) - (a*g_2*x)/((g_2 +y)**2 )
    a_31 = (p_2 * y) / (g_3 + y)
    a_32 = (g_3 * p_2 * x)/(  (g_3 + y)**2 )
    A = np.array([[-mu_2,c, 0],
                  [ a_21 , a_22, 0],
                  [a_31, a_32, -mu_3]])
    return A

A=matriz(x3,y3,z3)

# Calcula los valores propios
valores_propios, vectores_propios = np.linalg.eig(A)

# Imprime los resultados
print("Valores propios:")
print(valores_propios)
print("\nVectores propios:")
print(vectores_propios)
