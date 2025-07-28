# Aqui esta implementado el algoritmo descrito en el proyecto para la aproximacion de los puntos criticos en la 
# bifurcacion
import numpy as np
import sympy as sp


mu_2 = 0.03
p_1 = 0
a = 1
r_2 = 0.18
mu_3 = 10
p_2 = 5
g_1 = 2*10**7
g_2 = 10**5
g_3 = 10**3
b = 10**(-9)
s_1 = 0
s_2 = 0

t = sp.symbols('t')
x,y,z = sp.symbols("x,y,z")
c = sp.symbols('c')

dx = c*y - mu_2*x
dy = r_2*y*(1-b*y) - ((a*x*y)/(g_2+y))
dz = ((p_2*x*y)/(g_3+y)) - mu_3*z

M = sp.Matrix([dx,dy,dz])
J = M.jacobian([x,y,z])

#Estos los determine en otro script usando Matlab
x_1=-(2500000000000*c*((100*c)/3 + ((10000*c**2)/9 - (29997*c)/2500 + 8101620081/250000000000)**(1/2) - 89991/500000))/27
y_1= 499950000 - (25000000000*((10000*c**2)/9 - (29997*c)/2500 + 8101620081/250000000000)**(1/2))/9 - (2500000000000*c)/27
z_1=-(125000000*c*(50000000*c + (2500000000000000*c**2 - 26997300000000*c + 72914580729)**(1/2) - 269973)**2)/(81*(2500000000*c + 50*(2500000000000000*c**2 - 26997300000000*c + 72914580729)**(1/2) - 13498677))
x_2=(2500000000000*c*(((10000*c**2)/9 - (29997*c)/2500 + 8101620081/250000000000)**(1/2) - (100*c)/3 + 89991/500000))/27
y_2=(25000000000*((10000*c**2)/9 - (29997*c)/2500 + 8101620081/250000000000)**(1/2))/9 - (2500000000000*c)/27 + 499950000
z_2=(125000000*c*((2500000000000000*c**2 - 26997300000000*c + 72914580729)**(1/2) - 50000000*c + 269973)**2)/(81*(50*(2500000000000000*c**2 - 26997300000000*c + 72914580729)**(1/2) - 2500000000*c + 13498677))

Jacob_1 = J.subs([(x,x_1),(y,y_1),(z,z_1)])
Jacob_2 = J.subs([(x,x_2),(y,y_2),(z,z_2)])

#Voy a crear un algoritmo, que coja varios valores de c en un intervalo, sustituya y halle los valores propios, 
#y analice si ocurre un cambio de signos respecto al c anterior

#Con este par de funciones me permite poder sacar la estabilidad de los puntos
def all_positive(list):
    t = [complex(x) for x in list]
    return all(x.real>0 for x in t)
def all_negative(list):
    t = [complex(x) for x in list]
    return all(x.real<0 for x in t)
def estabilidad(list):
    if all_positive(list):
        return 'P'
    elif all_negative(list):
        return 'N'
    else:
        return 'M'

c0 = 0
cf = 0.05
n = 1000 #Para el resultado publicado utilizamos n=100000
c_vals = np.linspace(c0,cf,n)
c_camb_1 = []
c_camb_2 = []

#En resumidas cuentas te pongo los valores de c, valores propios, P si todos son positivos, N con negativos y M mixto

print("Procediendo a analizar para el punto 1.")
i=0
est = ''
est_ant = ''
while i < len(c_vals):
    Jacob_i = Jacob_1.subs(c,c_vals[i])
    Jacob_i_eigenvals = list(Jacob_i.eigenvals().keys()) #Saco los valores propios, cojo los valores propios y los meto en una lista
    est = estabilidad(Jacob_i_eigenvals)
    print("c=",c_vals[i],Jacob_i_eigenvals,est)
    if est_ant == '':
        est_ant = est
    elif est != est_ant:
        print("AQUI HA HABIDO UN CAMBIO para c=",c_vals[i])
        c_camb_1.append((c_vals[i],est_ant,est))
        est_ant = est
    else:
        est_ant = est
    i+=1


print("Procediendo a analizar para el punto 2.")
i=0
est = ''
est_ant = ''
while i < len(c_vals):
    Jacob_i = Jacob_2.subs(c,c_vals[i])
    Jacob_i_eigenvals = list(Jacob_i.eigenvals().keys()) #Saco los valores propios, cojo los valores propios y los meto en una lista
    est = estabilidad(Jacob_i_eigenvals) #Usando el criterio del signo de la parte real determinamos esto
    print("c=",c_vals[i],Jacob_i_eigenvals,est)
    if est_ant == '':
        est_ant = est
    elif est != est_ant:
        print("AQUI HA HABIDO UN CAMBIO para c=",c_vals[i])
        c_camb_2.append((c_vals[i],est_ant,est))
        est_ant = est
    else:
        est_ant = est
    i+=1

#Aqui va a sacar los puntos donde hay cambio, que ocurria antes y que ocurre en el siguiente punto
print("Para el punto 1",c_camb_1)
print("Para el punto 2",c_camb_2)
print("FIN")


