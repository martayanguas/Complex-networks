import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import os
import numpy as np
from scipy.integrate import solve_ivp

fig = plt.figure(figsize=(8, 4))  
gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05])  
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharey=ax1)
fig4, ax4 = plt.subplots()

lambda_values = [0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 25.0, 50.0]
cmap = cm.get_cmap('jet')
norm = colors.Normalize(vmin=0, vmax=50)
mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
colores = list(mcolors.TABLEAU_COLORS.values())[:5]

average_final_values = []

i=0
previous_steady_state=0
lambdas_ploted=[]
for lambda_val in lambda_values:
    filename = f"sim_{int(lambda_val*100)}.txt"
    if os.path.isfile(filename):
        time_values = []
        infection_density_values = []
        with open(filename, "r") as file:
            for line in file:
                time, infection_density = map(float, line.split())
                time_values.append(time)
                infection_density_values.append(infection_density)
                    
                if infection_density_values[-1] == 0:
                    time_values.append(25)
                    infection_density_values.append(0)

        img1 = ax1.plot(time_values, infection_density_values, label=f"\u03BB={lambda_val:.2f}", color = mapper.to_rgba(i))
        
        if infection_density_values[-1] == 0:
            average_final_value = 0
        else:
            average_final_value = np.mean(infection_density_values[-100:])  
        average_final_values.append(average_final_value)
        
        if average_final_value == 0:
            ax4.cla()
            ax4.plot(time_values, infection_density_values, label=f"\u03BB={lambda_val:.2f}", color=colores[0])
            lambdas_ploted = []
            lambdas_ploted.append(lambda_val)   
                 
        elif average_final_value - previous_steady_state > 0.1975:
            ax4.plot(time_values, infection_density_values, label=f"\u03BB={lambda_val:.2f}", color=colores[len(lambdas_ploted)])
            previous_steady_state = average_final_value
            lambdas_ploted.append(lambda_val)
                
        i+=50/len(lambda_values)

img1 = ax1.set_title("Evolution of the density of infected agents \n SIS simulation")
img1 = ax1.set_xlabel("Time")
img1 = ax1.set_ylabel("Density of infected agents, $\u03C1_i$")
#plt.legend(ncol=2, loc='center left', bbox_to_anchor=(1, 0.5))
#color_bar = plt.colorbar(mapper)
#color_bar.set_label('\u03BB')
#color_bar.set_ticks([])

def read_edge_list(filename):
    with open(filename, 'r') as file:
        next(file)
        next(file)
        edge_list = [tuple(map(int, line.strip().split())) for line in file]
    return edge_list

def initialize_vectors(edge_list):
    N = max(max(edge_list, key=max)) + 1  #number of nodes
    D = [0] * N
    E = len(edge_list) #number of links
    V = [None] * 2 * E
    
    for edge in edge_list:
        D[edge[0]] += 1
        D[edge[1]] += 1

    return N, V, D

def initialize_pointers(N, D):
    first_pointers = [0] * N
    second_pointers = [0] * N

    pos = 0
    for i in range(N):
        first_pointers[i] = pos
        second_pointers[i] = pos
        pos += D[i]

    return first_pointers, second_pointers

def store_neighbors(edge_list, V, first_pointers, second_pointers):
    for edge in edge_list:
        i, j = edge

        V[second_pointers[i]] = j
        second_pointers[i] += 1

        V[second_pointers[j]] = i
        second_pointers[j] += 1
    
    return V

def degree_distribution(D):
    max_k = max(D)
    count = [0] * (max_k + 1)

    for k in D:
        count[k] += 1

    total_nodes = len(D)
    P_k = [k_count / total_nodes for k_count in count]
    
    cc_P_k = [sum(P_k[i:]) for i in range(len(P_k))]
    
    return P_k, cc_P_k

delta = 1.0
filename = 'out.opsahl-powergrid'
edge_list = read_edge_list(filename)
N, V, D = initialize_vectors(edge_list)
first_pointers, second_pointers = initialize_pointers(N, D)
V = store_neighbors(edge_list, V, first_pointers, second_pointers)

def differential_equations(t, y, lambda_val, N, V, delta):
    drho_dt = np.zeros(N)
    
    for i in range(N):
        sum_neighbors = 0
        for neighbor in V[first_pointers[i]:second_pointers[i]]:
            sum_neighbors += y[neighbor]
            
        drho_dt[i] = -delta * y[i] + lambda_val * (1 - y[i]) * sum_neighbors
    
    return drho_dt

time = (0, 25)

infected_nodes = np.random.choice(N, size=int(N*0.01), replace=False)  
rho_init = np.zeros(N)
for node in infected_nodes:
    rho_init[node] = 1
    
steady_states = []
i=0
n=0
for lam in lambda_values:
    sol = solve_ivp(lambda t, y: differential_equations(t, y, lam, N, V, delta), time, rho_init, dense_output=True)
    
    rho = np.mean(sol.y, axis=0)
    
    img2 = ax2.plot(sol.t, rho, label=f'\u03BB={lam:.2f}', color = mapper.to_rgba(i))

    steady_states.append(rho[-1])
    i+=50/len(lambda_values)
    
    if lam in lambdas_ploted:
        ax4.plot(sol.t, rho, '--', color=colores[n])
        n+=1

img2 = ax2.set_title("Evolution of the density of infected agents \n Numerical solution")
img2 = ax2.set_xlabel("Time")
#img2 = ax2.set_ylabel("Density of infected agents, $\u03C1_{i}$")
#plt.legend(ncol=2, loc='center left', bbox_to_anchor=(1, 0.5))
cax = fig.add_subplot(gs[2])
color_bar = fig.colorbar(mapper, cax=cax)
color_bar.set_label('\u03BB')
color_bar.set_ticks([])

ax4.set_title("Evolution of the density of infected agents")
ax4.set_xlabel("Time")
ax4.set_ylabel("Density of infected agents, $\u03C1_{i}$")
ax4.legend(title = 'Infection rate', loc='center left', bbox_to_anchor=(1, 0.5))

fig3, ax3=plt.subplots()
ax3.plot(lambda_values, steady_states, 'o-', label='Numerical solution')
ax3.plot(lambda_values, average_final_values, 'o-',label='SIS simulation')
ax3.set_xlabel('Infection rate, \u03BB')
ax3.set_ylabel('Steady State, $\u03C1_{st}$')
ax3.legend(loc='center left', bbox_to_anchor=(0.55, 0.75))
ax3.axhline(y=1, color='gray', linestyle='dashed')

ax3_zoom = plt.axes([0.45, 0.2, 0.4, 0.4])  
ax3_zoom.plot(lambda_values[0:15], steady_states[0:15], 'o-', label='Zoom')
ax3_zoom.plot(lambda_values[0:15], average_final_values[0:15], 'o-', label='Zoom')

ax3_zoom.set_xlim(0, 0.5)
ax3_zoom.set_ylim(-0.01, 0.1)
    
plt.show()