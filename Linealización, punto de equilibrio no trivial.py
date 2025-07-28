import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

#PUNTO POSITIVO
x_eq = 2.2958*10**4
y_eq = 2.7550*10**4
z_eq = 1.1077*10**4

t = sp.symbols('t')
x = sp.Function('x')(t)
y = sp.Function('y')(t)
z = sp.Function('z')(t)

#El jacobiano lo determinamos en otro código utilizando Matlab
jacob = np.array([[-0.0300,0.0250,0],[-0.2160,0.0389,0],[4.8249,0.1408,-10]])



dx = jacob[0][0]*(x-x_eq) + jacob[0][1]*(y-y_eq) + jacob[0][2]*(z-z_eq)
dy = jacob[1][0]*(x-x_eq) + jacob[1][1]*(y-y_eq) + jacob[1][2]*(z-z_eq)
dz = jacob[2][0]*(x-x_eq) + jacob[2][1]*(y-y_eq) + jacob[2][2]*(z-z_eq)

eq1 = sp.Eq(sp.Derivative(x,t),dx)
eq2 = sp.Eq(sp.Derivative(y,t),dy)
eq3 = sp.Eq(sp.Derivative(z,t),dz)

#Soluciones generales
system = [eq1,eq2,eq3]
sol = sp.dsolve(system)
print("Soluciones generales")
for s in sol:
        print(s)


#Solucion particular para un punto cercano al equilibrio
p_inic = (x_eq+10,y_eq+10,z_eq+10)

C1, C2, C3 = sp.symbols('C1,C2,C3')
constantes = sp.solve([sol[0].rhs.subs(t,0) - p_inic[0], sol[1].rhs.subs(t,0) - p_inic[1], sol[2].rhs.subs(t,0) - p_inic[2]], [C1,C2,C3])

sol_part = [s.subs(constantes) for s in sol]
print("Soluciones generales")
for s in sol_part:
    print(s)

#Vamos a graficar

t_vals = np.linspace(0,10,40000)
x_sol = sp.lambdify(t, sol_part[0].rhs, 'numpy')
y_sol = sp.lambdify(t, sol_part[1].rhs, 'numpy')
z_sol = sp.lambdify(t, sol_part[2].rhs, 'numpy')

x_vals = x_sol(t_vals)
y_vals = y_sol(t_vals)
z_vals = z_sol(t_vals)

fig = plt.figure(figsize=(10,5))
ax_1 = fig.add_subplot(1,3,1)
ax_2 = fig.add_subplot(1,3,2)
ax_3 = fig.add_subplot(1,3,3)


ax_1.set_title("(Lineal p)X:Celulas Inmunes")
ax_2.set_title("Y:Celulas Tumorosas")
ax_3.set_title("Z: IL")

ax_1.plot(t_vals,x_vals)
ax_2.plot(t_vals,y_vals)
ax_3.plot(t_vals,z_vals)

plt.show()