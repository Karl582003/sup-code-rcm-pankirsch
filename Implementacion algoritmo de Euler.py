#En este código está la implementación del método de Euler para el sistema
import numpy as np
import matplotlib.pyplot as plt


#Definimos los parámetros

#Primero los que nos serviran para el escalado
#Esto realmente no lo llegamos a utilizar
t_s = 1
E_0 = 1
I_0 = 1
T_0 = 1

#Ahora los parametros como tal del modelo
c = 0.025 #Vale la pena decir que este fue el principal que modificamos para bifurcación
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


def sist(x,y,z): #Definimos el sistema de EDOs
    dx = (c*y - mu_2*x + s_1)/(t_s*E_0)
    dy = (r_2*y*(1-b*y) - ((a*x*y)/(g_2+y)))/(t_s*T_0)
    dz = (((p_2*x*y)/(g_3+y)) - mu_3*z +s_2)/(t_s*I_0)
    return dx,dy,dz


#Empezamos el algoritmo
#Parametros del algoritmo
n = 500000 #Pasos del algoritmo
t_0 = 0 #Tiempo inicial
t_f = 10 #Tiempo final
h= (t_f-t_0)/n #Calculamos el paso del algoritmo
pto_inicial = np.array((1,1,1)) #Los valores iniciales de x,y,z

#Metodo de Euler
t = np.arange(t_0,t_f+h,h) #Intervalo de tiempo
w_output = [] #En este array guardaremos los valores del sistema a traves del tiempo

i=0

while i<len(t):
    if i==0:
        w_output.append(pto_inicial)
    else:
        w_ant = w_output[-1]
        dw = np.array(sist(w_ant[0],w_ant[1],w_ant[2]))
        w = w_ant+ h*dw
        w_output.append(w)
    i+=1



print("Paso:",h)    
print("Resultado:",w_output[-1])

#Vamos a plottear el resultado
x = []
y = []
z = []
for i in range(len(t)):
    x.append(w_output[i][0])
    y.append(w_output[i][1])
    z.append(w_output[i][2])

"""
fig = plt.figure(figsize=(10,5))
ax_1 = fig.add_subplot(1,3,1)
ax_2 = fig.add_subplot(1,3,2)
ax_3 = fig.add_subplot(1,3,3)


ax_1.set_title("EULER X:Celulas Inmunes")
ax_2.set_title("Y:Celulas Tumorosas")
ax_3.set_title("Z: IL")

ax_1.plot(t,x)
ax_2.plot(t,y)
ax_3.plot(t,z)

plt.show()
"""


#Los tres juntos

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot()
ax.set_title("Euler ")
ax.plot(t,x,"r", label="X:Celulas Inmunes")
ax.plot(t,y,"b", label="Y:Celulas Tumorosas")
ax.plot(t,z,'black', label="Z:IL")
ax.legend()

plt.show()