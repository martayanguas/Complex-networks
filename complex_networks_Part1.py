# -*- coding: utf-8 -*-
"""
Created on Fri May 12 11:21:15 2023

@author: Marta
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random

def read_edge_list(filename):
    with open(filename, 'r') as file:
        next(file)
        next(file)
        edge_list = [tuple(map(int, line.strip().split())) for line in file]
    return edge_list

def initialize_vectors(edge_list):
    N = max(max(edge_list, key=max))   #number of nodes
    D = [0] * (N+1)
    E = len(edge_list) #number of links
    V = [None] * 2 * E
    
    for edge in edge_list:
        D[edge[0]] += 1
        D[edge[1]] += 1

    return N, V, D

def initialize_pointers(N, D):
    first_pointers = [0] * (N+1)
    second_pointers = [0] * (N+1)

    pos = 0
    for i in range(N+1):
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

def average_nearest_neighbor_degree(first_pointers, D, V):
    max_degree = max(D)
    k_nn_sum = [0] * (max_degree + 1)
    k_nn_count = [0] * (max_degree + 1)

    for node in range(len(D)):
        node_degree = D[node]
        neigh = V[first_pointers[node]:first_pointers[node]+node_degree]
        k_nn_sum[node_degree] += sum(D[neighbor] for neighbor in neigh)
        k_nn_count[node_degree] += len(neigh)

    k_nn = [k_nn_sum[i] / k_nn_count[i] if k_nn_count[i] != 0 else 0 for i in range(len(k_nn_sum))]
    
    return k_nn

def clustering_coefficient(first_pointers, D, V):
    clustering = [0] * len(D)
    triangles = [0] * len(D)
    
    max_degree = max(D)
    C_sum = [0] * (max_degree + 1)
    C_count = [0] * (max_degree + 1)

    for node in range(len(D)):
        node_degree = D[node]
        if node_degree < 2:
            continue

        neighbors = V[first_pointers[node]:first_pointers[node]+node_degree]
        for i in range(len(neighbors)):
            for j in range(i+1, len(neighbors)):
                node_i_neighbors = V[first_pointers[neighbors[i]]:first_pointers[neighbors[i]]+D[neighbors[i]]]
                if neighbors[j] in node_i_neighbors:
                    triangles [node] += 1
                
        clustering[node] = 2 * triangles[node] / (node_degree * (node_degree-1))
        C_sum[node_degree] += clustering [node]
        C_count[node_degree] += 1
        
    C_avg = [C_sum[i] / C_count[i] if C_count[i] != 0 else 0 for i in range(len(C_sum))]
    
    return C_avg
    
def try_configuration_model(D):
    nodes_with_stubs = [(node, degree) for node, degree in enumerate(D) if degree > 0]
    edge_list = []

    while nodes_with_stubs:
        random.shuffle(nodes_with_stubs)
        node, degree = nodes_with_stubs.pop(0)

        if degree > len(nodes_with_stubs):
            return []

        num_edges_added = 0
        for other_node, other_degree in nodes_with_stubs[:]:
            if node != other_node and (node, other_node) not in edge_list and (other_node, node) not in edge_list:
                edge_list.append((node, other_node))
                nodes_with_stubs.remove((other_node, other_degree))
                if other_degree > 1:
                    nodes_with_stubs.append((other_node, other_degree - 1))

                num_edges_added += 1
                if num_edges_added == degree:
                    break

        if num_edges_added < degree:
            return []

    return edge_list


def configuration_model(D):
    result = try_configuration_model(D)

    while not result:
        result = try_configuration_model(D)

    return result


def write_output(N, num_links, V, first_pointers, second_pointers):
    with open('data_network.txt', 'w') as f:
        f.write(f'{N}\n')
        f.write(f'{num_links}\n')
        for v in V:
            f.write(f'{v}\n')
        for fp in first_pointers:
            f.write(f'{fp}\n')
        for sp in second_pointers:
            f.write(f'{sp}\n')

def main():
    ################ Assigment 1
    filename = 'out.opsahl-powergrid'
    edge_list = read_edge_list(filename)
    N, V, D = initialize_vectors(edge_list)
    first_pointers, second_pointers = initialize_pointers(N, D)
    V = store_neighbors(edge_list, V, first_pointers, second_pointers)
    num_links = len(edge_list)
    avg_degree = np.mean(D)

    print("Number of nodes:", N)
    print("Number of links:", num_links)
    print("List of degrees:")
    for i in range(1,N):
        print("Node", i+1, "Degree", D[i])
    print("Average degree:", avg_degree)
    
    write_output(N, num_links, V, first_pointers, second_pointers)
    
    for i in range(N):
        if D[i] != V.count(i):
            print('ERROR', D[i], V.count(i))
            
    G = nx.Graph()
    G.add_edges_from(edge_list)
    
    for i in range(1, N):
        nei = list(G.neighbors(i))
        if not np.all(nei == V[first_pointers[i]:first_pointers[i]+D[i]]):
            print('ERROR', "Los vecinos de", i, " en G son {}".format(list(G.neighbors(i))), 'not', V[first_pointers[i]:first_pointers[i]+D[i]])
            
    
    ################ Assigment 2
    P_k, cc_P_k = degree_distribution(D)
    fig1, ax1 = plt.subplots()
    k_vec = np.arange(1, len(P_k)+1)
    ax1.plot(k_vec, P_k, 'o-', label='Original Network')
    ax1.set_xlabel('Degree (k)')
    ax1.set_ylabel('P(k)')
    ax1.set_title('Degree Distribution')
    ax1.set_xscale('log')
    ax1.set_yscale('log')   
    fig1.savefig('degree_distribution.png')
    
    fig2, ax2 = plt.subplots()
    ax2.plot(k_vec, cc_P_k, 'o-', label='Original Network')
    ax2.set_xlabel('Degree (k)')
    ax2.set_ylabel('cc_P(k)')
    ax2.set_title('Complementary Cumulative Degree Distribution')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    fig2.savefig('complementary_cumulative_degree_distribution.png')
    
    k_nn = average_nearest_neighbor_degree(first_pointers, D, V)
    fig3, ax3 = plt.subplots()
    ax3.plot(k_vec, k_nn, 'o-', label='Original Network')
    ax3.set_xlabel('Degree (k)')
    ax3.set_ylabel('k_nn(k)')
    ax3.set_title('Average Nearest Neighbors Degree')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    fig3.savefig('average_nearest_neighbours_degree.png')
    
    C = clustering_coefficient(first_pointers, D, V)
    fig4, ax4 = plt.subplots()
    ax4.plot(k_vec, C, 'o-', label='Original Network')
    ax4.set_xlabel('Degree (k)')
    ax4.set_ylabel('C(k)')
    ax4.set_title('Clustering Coefficient')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    fig4.savefig('clustering_coefficient.png')
    
    clustering = nx.clustering(G)
    degrees = dict(G.degree())
    
    # Inicializa un diccionario para guardar los coeficientes de clustering de cada grado
    clustering_per_degree = {}

    for node, clust_coeff in clustering.items():
        degree = degrees[node]
        
        # Añade el coeficiente al grado correspondiente en el diccionario
        if degree in clustering_per_degree:
            clustering_per_degree[degree].append(clust_coeff)
        else:
            clustering_per_degree[degree] = [clust_coeff]
            
    # Calcula el promedio de los coeficientes de clustering para cada grado
    avg_clustering_per_degree = {degree: sum(values)/len(values) for degree, values in clustering_per_degree.items()}

    for node, clust_coeff in avg_clustering_per_degree.items():
        if not np.isclose(clust_coeff, C[node], atol=1e-4):
            print(node, clust_coeff, C[node])
    
    ############# Assignment 3
    if max(D) < np.sqrt(avg_degree*N):
        # Generate a configuration model graph
        cm_edge_list = configuration_model(D)
    
        # Initialize vectors and pointers for the new graph
        N_cm, V_cm, D_cm = initialize_vectors(cm_edge_list)
        num_links = len(cm_edge_list)
        avg_degree = np.mean(D_cm)
        
        print("Number of nodes:", N_cm)
        print("Number of links:", num_links)
        print("Average degree:", avg_degree)
        
        first_pointers_cm, second_pointers_cm = initialize_pointers(N_cm, D_cm)
        V_cm = store_neighbors(cm_edge_list, V_cm, first_pointers_cm, second_pointers_cm)
        
        # Compute the topological properties
        P_k_cm, cc_P_k_cm =degree_distribution(D_cm)
        k_vec = np.arange(1, len(P_k)+1)
        ax1.plot(k_vec, P_k_cm, 'o-', label='Configuration Model')
        ax1.legend()
        fig1.savefig('degree_distribution2.png')
        ax2.plot(k_vec, cc_P_k_cm, 'o-', label='Configuration Model')
        ax2.legend()
        fig2.savefig('complementary_cumulative_degree_distribution2.png')
        
        k_nn_cm = average_nearest_neighbor_degree(first_pointers_cm, D_cm, V_cm)
        ax3.plot(k_vec, k_nn_cm, 'o-', label='Configuration Model')
        ax3.legend()
        fig3.savefig('average_nearest_neighbor2.png')
        
        C_cm = clustering_coefficient(first_pointers_cm, D_cm, V_cm)
        ax4.plot(k_vec, C_cm, 'o-', label='Configuration Model')
        ax4.legend()
        fig4.savefig('clustering_coefficient2.png')
        
        CM_realizations=100 #This gives the number of Configurational Models we will create
        P_k_total = np.zeros(len(k_vec))
        cc_P_k_total = np.zeros(len(k_vec))
        k_nn_total = np.zeros(len(k_vec))
        C_total = np.zeros(len(k_vec))
        for i in range(CM_realizations):
            print('CM realization number', i)
            cm_edge_list = configuration_model(D)
    
            N_cm, V_cm, D_cm = initialize_vectors(cm_edge_list)
            first_pointers_cm, second_pointers_cm = initialize_pointers(N_cm, D_cm)
            V_cm = store_neighbors(cm_edge_list, V_cm, first_pointers_cm, second_pointers_cm)

            #Topological properties
            P_k_cm, cc_P_k_cm = degree_distribution(D_cm)
            k_nn_cm = average_nearest_neighbor_degree(first_pointers_cm, D_cm, V_cm)
            C_cm = clustering_coefficient(first_pointers_cm, D_cm, V_cm)
            
            for j in range(len(P_k_total)):
                P_k_total[j] += P_k_cm[j]
                cc_P_k_total[j] += cc_P_k_cm[j]
                k_nn_total[j] += k_nn_cm[j]
                C_total[j] += C_cm[j]
                
        P_k_avg = P_k_total / CM_realizations
        k_vec = np.arange(1, len(P_k)+1)
        ax1.plot(k_vec, P_k_avg, 'o-', label='Average 100 CM')
        ax1.legend()
        fig1.savefig('degree_distribution3.png')
        cc_P_k_avg = cc_P_k_total / CM_realizations
        ax2.plot(k_vec, cc_P_k_avg, 'o-', label='Average 100 CM')
        ax2.legend()
        fig2.savefig('complementary_cumulative_degree_distribution3.png')
        
        k_nn_avg = k_nn_total / CM_realizations
        ax3.plot(k_vec, k_nn_avg, 'o-', label='Average 100 CM')
        ax3.legend()
        fig3.savefig('average_nearest_neighbor3.png')
        
        C_avg = C_total / CM_realizations
        ax4.plot(k_vec, C_avg, 'o-', label='Average 100 CM')
        ax4.legend()
        fig4.savefig('clustering_coefficient3.png')
                    
if __name__ == "__main__":
    main()
