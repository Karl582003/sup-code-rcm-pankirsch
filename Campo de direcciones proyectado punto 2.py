# Con este código graficamos la representación en los planos coordenados del campo de direcciones
# del punto de equilibrio.

import numpy as np
import matplotlib.pyplot as plt

# Parámetros del sistema
c = 0.025 #Vale la pena decir que se bifurcó con este
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

xo,yo,zo = 22958,27550,11077
eps = 100
# Definir el campo vectorial
def campo_vectorial(x, y, z):
    dx = c * y - mu_2 * x
    dy = y * (r_2 * (1 - b * y) - (a * x / (g_2 + y)))
    dz = (p_2 * x * y) / (g_3 + y) - mu_3 * z
    return dx, dy, dz

# Generar una malla de puntos
x = np.linspace(xo-eps, xo+eps, 10)
y = np.linspace(yo-eps, yo+eps, 10)
z = np.linspace(zo-eps, zo+eps, 10)

X, Y = np.meshgrid(x, y)
X, Z = np.meshgrid(x, z)
Y, Z = np.meshgrid(y, z)


# Calcular el campo de direcciones para cada plano
DX_XY, DY_XY, _ = campo_vectorial(X, Y, zo)
_, DY_YZ, DZ_YZ = campo_vectorial(xo, Y, Z)
DX_XZ, _, DZ_XZ = campo_vectorial(X, yo, Z)

# Graficar el campo de direcciones en el plano XY
plt.figure(figsize=(15, 6))

plt.subplot(1, 3, 1)
plt.plot(xo, yo, 'ko')
plt.quiver(X, Y, DX_XY, DY_XY, color='blue')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Campo de Direcciones en el Plano XY')

# Graficar el campo de direcciones en el plano YZ
plt.subplot(1, 3, 2)
plt.plot(yo, zo, 'ko')
plt.quiver(Y, Z, DY_YZ, DZ_YZ, color='green')
plt.xlabel('Y')
plt.ylabel('Z')
plt.title('Campo de Direcciones en el Plano YZ')

# Graficar el campo de direcciones en el plano XZ
plt.subplot(1, 3, 3)
plt.plot(xo, zo, 'ko')
plt.quiver(X, Z, DX_XZ, DZ_XZ, color='red')
plt.xlabel('X')
plt.ylabel('Z')
plt.title('Campo de Direcciones en el Plano XZ')

plt.suptitle("Proyecciones del Campo de direcciones del sistema")
plt.show()
