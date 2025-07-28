import numpy as np
import matplotlib.pyplot as plt


#Definimos los parámetros

#Primero, los que nos servirán para el escalado
t_s = 1
E_0 = 1
I_0 = 1
T_0 = 1

#Ahora, los parámetros como tal del modelo
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


def sist(w): #Definimos el sistema de EDOs
    x = w[0]
    y = w[1]
    z = w[2]
    dx = -mu_2 * x + c*y
    dy = r_2 * y
    dz = -mu_3*z
    output = np.array((dx,dy,dz))
    return output


#Empezamos el algoritmo
#Parámetros del algoritmo
n = 50000 #Pasos del algoritmo
t_0 = 0 #Tiempo inicial
t_f = 10 #Tiempo final
h= (t_f-t_0)/n #Calculamos el paso del algoritmo
pto_inicial = np.array((1,1,1)) #Los valores iniciales de x,y,z

#Metodo RK4
t = np.arange(t_0,t_f+h,h) #Intervalo de tiempo
w_output = [] #En este array guardaremos los valores del sistema a traves del tiempo

i=0

while i<len(t):
    if i==0:
        w_output.append(pto_inicial)
    else:
        k_1 = sist(w_output[-1])
        k_2 = sist(w_output[-1]+h*k_1/2)
        k_3 = sist(w_output[-1]+h*k_2/2)
        k_4 = sist(w_output[-1]+h*k_3)
        w = w_output[-1] + (h/6)*(k_1+2*k_2+2*k_3+k_4)
        w_output.append(w)
    i+=1



print("Paso:",h)    
print("Resultado:",w_output[-1])

#Vamos a plottear esta funcion
x = []
y = []
z = []
for i in range(len(t)):
    x.append(w_output[i][0])
    y.append(w_output[i][1])
    z.append(w_output[i][2])

fig = plt.figure(figsize=(10,5))

ax_1 = fig.add_subplot(1,3,1)
plt.grid()
ax_2 = fig.add_subplot(1,3,2)
plt.grid()
ax_3 = fig.add_subplot(1,3,3)
plt.grid()

ax_1.set_title("u(t)")
ax_2.set_title("v(t)")
ax_3.set_title("w(t)")

ax_1.plot(t,x,color='b')
ax_2.plot(t,y,color='g')
ax_3.plot(t,z,color='r')
plt.show()