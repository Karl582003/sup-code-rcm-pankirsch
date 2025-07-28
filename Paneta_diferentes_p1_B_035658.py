import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Fixed parameters (according to the table)
mu2 = 0.03
mu3 = 10
p2 = 5
a = 1
r2 = 0.18
b = 1e-9
c = 0.05  # Using the maximum value of c
g1 = 2e7
g2 = 1e5
g3 = 1e3
s1 = 0    # Assuming no external source of effector cells
s2 = 0    # Assuming no external source of IL-2

# Values of p1 to test
p1_values = 1e1 * np.array([12.45, 1.245, 0.1245, 0.01245, 0.001245, 0.0001245, 0])

# Formatear las etiquetas de p1 para la leyenda
p1_labels = [f'p₁ = {p:.5f}'.rstrip('0').rstrip('.') if p != 0 else 'p₁ = 0' for p in p1_values]

# Initial conditions
x0 = 1
y0 = 1
z0 = 1
initial_conditions = [x0, y0, z0]

# Time range (adjustable)
t_span = (0, 100)
t_eval = np.linspace(0, 100, 1000)

# Configuration for solve_ivp
rtol = 1e-6
atol = 1e-9

# Curve style configuration
line_styles = ['-', '-', '-', '-.', '-', '-', '-']
markers = ['s', 'x', 'd', '^', '<', '+', 'o']
colors = ['m', 'b', 'y', 'r', 'g', 'c', 'k']
line_width = 1.0
marker_spacing = 50  # Marker spacing (points every 'marker_spacing' steps)

# Initialize solution storage
Y1 = np.zeros((len(t_eval), 3, len(p1_values)))
Y_markers = np.zeros((len(t_eval[::marker_spacing]), 3, len(p1_values)))

# Create figure with larger size for better visualization
plt.figure(figsize=(18, 6))

# Solve the system for each p1 value
for i, p1 in enumerate(p1_values):
    # Define the ODE system
    def ode_system(t, Y):
        x, y, z = Y
        dxdt = c * y - mu2 * x + (p1 * x * z) / (g1 + z) + s1
        dydt = r2 * y * (1 - b * y) - (a * x * y) / (g2 + y)
        dzdt = (p2 * x * y) / (g3 + y) - mu3 * z + s2
        return [dxdt, dydt, dzdt]
    
    # Numerical solution
    sol = solve_ivp(ode_system, t_span, initial_conditions, 
                    t_eval=t_eval, rtol=rtol, atol=atol, method='LSODA')
    Y1[:, :, i] = sol.y.T

# Get markers data
t_markers = t_eval[::marker_spacing]
Y_markers = Y1[::marker_spacing, :, :]

# Plot with distinct styles for each p1
# Subplot 1: Effector cells (x(t))
plt.subplot(1, 3, 1)
for i in range(len(p1_values)):
    # Plot solid line
    plt.plot(t_eval, Y1[:, 0, i], color=colors[i], linestyle=line_styles[i], 
             linewidth=line_width, label=p1_labels[i])
    # Plot markers
    plt.plot(t_markers, Y_markers[:, 0, i], markers[i], color=colors[i], 
             markersize=8, markeredgewidth=1, linestyle='None')
plt.xlabel('Tiempo (días)', fontsize=14)
plt.ylabel('Células efectoras (x(t))', fontsize=14)
plt.title('Dinámica de las células efectoras', fontsize=14)
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.box(True)

# Subplot 2: Tumor cells (y(t))
plt.subplot(1, 3, 2)
for i in range(len(p1_values)):
    plt.plot(t_eval, Y1[:, 1, i], color=colors[i], linestyle=line_styles[i], 
             linewidth=line_width, label=p1_labels[i])
    plt.plot(t_markers, Y_markers[:, 1, i], markers[i], color=colors[i], 
             markersize=8, markeredgewidth=1, linestyle='None')
plt.xlabel('Tiempo (días)', fontsize=14)
plt.ylabel('Células tumorales (y(t))', fontsize=14)
plt.title('Dinámica de las células tumorales', fontsize=14)
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.box(True)

# Subplot 3: IL-2 concentration (z(t))
plt.subplot(1, 3, 3)
for i in range(len(p1_values)):
    plt.plot(t_eval, Y1[:, 2, i], color=colors[i], linestyle=line_styles[i], 
             linewidth=line_width, label=p1_labels[i])
    plt.plot(t_markers, Y_markers[:, 2, i], markers[i], color=colors[i], 
             markersize=8, markeredgewidth=1, linestyle='None')
plt.xlabel('Tiempo (días)', fontsize=14)
plt.ylabel('Concentración de IL-2 (z(t))', fontsize=14)
plt.title('Dinámica de la concentración de IL-2', fontsize=14)
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.box(True)

# Create a single legend for all subplots
handles, labels = plt.gca().get_legend_handles_labels()
plt.figlegend(handles, labels, loc='lower center', ncol=len(p1_values), 
              fontsize=10, bbox_to_anchor=(0.5, -0.1), title='Valores del parámetro p₁')

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)  # Adjust bottom margin for the legendplt.legend(title='Área Principal')
plt.legend()
plt.show()