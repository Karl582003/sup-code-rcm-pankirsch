import numpy as np
import matplotlib.pyplot as plt

#Definimos los valores de los parámetros
c = 0.025
mu_2 = 0.03
r_2 = 0.18

#Definimos el rango de tiempo
t = np.linspace(0, 10, 400)

#Definimos las funciones
x = (1 - c / (mu_2 + r_2)) * np.exp(-mu_2 * t) + (c / (mu_2 + r_2)) * np.exp(r_2 * t)
y = np.exp(r_2 * t)
z = np.exp(-10 * t)


fig = plt.figure(figsize=(10,6))
ax_1 = fig.add_subplot(1,3,1)
ax_2 = fig.add_subplot(1,3,2)
ax_3 = fig.add_subplot(1,3,3)


ax_1.set_title("u(t)")
ax_2.set_title("v(t)")
ax_3.set_title("w(t)")

# Graficar x(t)
fig.suptitle("Solución analítica del sistema linealizado.")
ax_1.plot(t, x, label='u(t)', color='b')

# Graficar y(t)
ax_2.plot(t, y, label='v(t)', color='g')

# Graficar z(t)
ax_3.plot(t, z, label='w(t)', color='r')

plt.show()