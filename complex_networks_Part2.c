#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int N;
int num_links;
int *V;
int *first_pointers;
int *second_pointers;

int *infected_nodes;
int num_infected_nodes;

int **active_links;
int num_active_links;

int *active_link_positions;

void read_data(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%d", &N);
    N+=1;

    fscanf(file, "%d", &num_links);

    V = (int*) malloc(num_links * 2 * sizeof(int));
    first_pointers = (int*) malloc(N * sizeof(int));
    second_pointers = (int*) malloc(N * sizeof(int));

    for (int i = 0; i < num_links * 2; i++) {
        fscanf(file, "%d", &V[i]);
    }

    for (int i = 0; i < N; i++) {
        fscanf(file, "%d", &first_pointers[i]);
    }

    for (int i = 0; i < N; i++) {
        fscanf(file, "%d", &second_pointers[i]);
    }

    fclose(file);
}

void infect_random_nodes() {
    infected_nodes = (int*) malloc(N * sizeof(int));

    num_infected_nodes = 0;

    srand(time(NULL));

    int num_nodes_to_infect = (N-1) / 100;
    for (int i = 0; i < num_nodes_to_infect; i++) {
        int node_to_infect;
        do {
            node_to_infect = rand() % (N-1) + 1;

            int already_infected = 0;
            for (int j = 0; j < num_infected_nodes; j++) {
                if (infected_nodes[j] == node_to_infect) {
                    already_infected = 1;
                    break;
                }
            }

            if (!already_infected) {
                infected_nodes[num_infected_nodes] = node_to_infect;
                num_infected_nodes++;
                break;
            }
        } while (1);
    }
}

void find_active_links() {
    active_links = (int**) malloc(num_links * 2 * sizeof(int*));  // Over-allocate for simplicity

    num_active_links = 0;

    for (int i = 0; i < num_infected_nodes; i++) {
        int infected_node = infected_nodes[i];

        int first_pointer = first_pointers[infected_node];
        int second_pointer = second_pointers[infected_node];
        
        for (int j = first_pointer; j < second_pointer; j++) {
            int connected_node = V[j];
            
            int is_infected = 0;
            for (int k = 0; k < num_infected_nodes; k++) {
                if (infected_nodes[k] == connected_node) {
                    is_infected = 1;
                    break;
                }
            }

            if (!is_infected) {
                active_links[num_active_links] = (int*) malloc(2 * sizeof(int));
                active_links[num_active_links][0] = infected_node;
                active_links[num_active_links][1] = connected_node;
                num_active_links++;
            }
        }
    }
}

void create_active_link_positions() {
    active_link_positions = (int*) malloc(num_links * 2 * sizeof(int));

    for (int i = 0; i < num_links * 2; i++) {
        active_link_positions[i] = -1;
    }

    for (int i = 0; i < num_active_links; i++) {
        int node1 = active_links[i][0];
        int node2 = active_links[i][1];

        for (int j = first_pointers[node1]; j < second_pointers[node1]; j++) {
            if (V[j] == node2) {
                active_link_positions[j] = i;
                break;
            }
        }
        for (int j = first_pointers[node2]; j < second_pointers[node2]; j++) {
            if (V[j] == node1) {
                active_link_positions[j] = i;
                break;
            }
        }
    }
}

void gillespie_simulation(double lambda) {
    srand(time(NULL));

    double total_time = 0.0;

    char filename[50];
    sprintf(filename, "sim_%d.txt", (int)(lambda * 100));
    
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    
    while (total_time < 25 && num_infected_nodes > 0 && num_infected_nodes < N) {
        double total_rate = num_infected_nodes + lambda * num_active_links;
        
        double recovery_prob = num_infected_nodes / total_rate;

        double rand_num = (double)rand() / RAND_MAX;
        
        if (rand_num <= recovery_prob) { 
            int recover_index = rand() % num_infected_nodes;
            int recover_node = infected_nodes[recover_index];
            
            infected_nodes[recover_index] = infected_nodes[num_infected_nodes - 1];
            num_infected_nodes--;

            for (int i = first_pointers[recover_node]; i < second_pointers[recover_node]; i++) {
                int connected_node = V[i];
                if (active_link_positions[i] >= 0) { 
                    int active_link_index = active_link_positions[i];
                    if (active_link_index == num_active_links -1){
                        free(active_links[active_link_index]);
                    } else {
                        active_links[active_link_index][0] = active_links[num_active_links - 1][0];
                        active_links[active_link_index][1] = active_links[num_active_links - 1][1];
                        
                        for (int j = first_pointers[active_links[num_active_links - 1][0]]; j < second_pointers[active_links[num_active_links - 1][0]]; j++){
                            if (V[j] == active_links[num_active_links - 1][1]) {
                                active_link_positions[j] = active_link_index;
                            }
                        }
                        for (int j = first_pointers[active_links[num_active_links - 1][1]]; j < second_pointers[active_links[num_active_links - 1][1]]; j++){
                            if (V[j] == active_links[num_active_links - 1][0]) {
                                active_link_positions[j] = active_link_index;
                            }
                        }
                    }
                    num_active_links--;

                    active_link_positions[i] = -1;

                    for (int j = first_pointers[V[i]]; j < second_pointers[V[i]]; j++) {
                        if (V[j] == recover_node) {
                            active_link_positions[j] = -1;
                            break;
                        }
                    }
                } else {
                    active_links[num_active_links] = (int*) malloc(2 * sizeof(int));
                    active_links[num_active_links][0] = connected_node;
                    active_links[num_active_links][1] = recover_node;
                    
                    for (int m = first_pointers[recover_node]; m < second_pointers[recover_node]; m++) {
                        if (V[m] == connected_node) {
                            active_link_positions[m] = num_active_links;
                            break;
                        }
                    }
                    for (int m = first_pointers[connected_node]; m < second_pointers[connected_node]; m++) {
                        if (V[m] == recover_node) {
                            active_link_positions[m] = num_active_links;
                            break;
                        }
                    }
                    num_active_links++;
                }
            }

        } else { 
            int infect_index = rand() % num_active_links;
            int infect_node = active_links[infect_index][1];
            
            infected_nodes[num_infected_nodes] = infect_node;
            num_infected_nodes++;

            for (int i = first_pointers[infect_node]; i < second_pointers[infect_node]; i++) {
                int connected_node = V[i];
                
                if (active_link_positions[i] >= 0) { 
                    int active_link_index = active_link_positions[i];
                    if (active_link_index == num_active_links -1){
                        free(active_links[active_link_index]);
                    } else {
                        active_links[active_link_index][0] = active_links[num_active_links - 1][0];
                        active_links[active_link_index][1] = active_links[num_active_links - 1][1];

                        for (int j = first_pointers[active_links[num_active_links - 1][0]]; j < second_pointers[active_links[num_active_links - 1][0]]; j++){
                            if (V[j] == active_links[num_active_links - 1][1]) {
                                active_link_positions[j] = active_link_index;
                            }
                        }
                        for (int j = first_pointers[active_links[num_active_links - 1][1]]; j < second_pointers[active_links[num_active_links - 1][1]]; j++){
                            if (V[j] == active_links[num_active_links - 1][0]) {
                                active_link_positions[j] = active_link_index;
                            }
                        }
                    }    
                    num_active_links--;
                    active_link_positions[i] = -1;

                    for (int j = first_pointers[V[i]]; j < second_pointers[V[i]]; j++) {
                        if (V[j] == infect_node) {
                            active_link_positions[j] = -1;
                            break;
                        }
                    }   
                } else {
                    active_links[num_active_links] = (int*) malloc(2 * sizeof(int));
                    active_links[num_active_links][0] = infect_node;
                    active_links[num_active_links][1] = connected_node;
                    
                    for (int m = first_pointers[infect_node]; m < second_pointers[infect_node]; m++) {
                        if (V[m] == connected_node) {
                            active_link_positions[m] = num_active_links;
                            break;
                        }
                    }
                    for (int m = first_pointers[connected_node]; m < second_pointers[connected_node]; m++) {
                        if (V[m] == infect_node) {
                            active_link_positions[m] = num_active_links;
                            break;
                        }
                    }

                    num_active_links++;
                }
            }
        }
        double u = (double)rand() / RAND_MAX + 1e-10; 
        double time_step = -log(u) / total_rate;
        total_time += time_step;
        double infection_density = (double) num_infected_nodes / N;
        fprintf(file, "%f %f\n", total_time, infection_density);
        }
    fclose(file);
}

int main() {
    // int result = system("python complex_networks.py"); this line call the python code you need previously

    read_data("data_network.txt");
    
    double lambdas[] = {0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 25.0, 50.0};
    int num = sizeof(lambdas) / sizeof(lambdas[0]);
    for (int i = 0; i < num; i++) {

        infect_random_nodes();
        find_active_links();
        create_active_link_positions();
        
        double lambda = lambdas[i];
        printf("lambda: %f\n", lambda);
        gillespie_simulation(lambda);
    }
    
    // int result = system("results.py"); this line call the python code for the results

    // free the memory 
    free(V);
    free(first_pointers);
    free(second_pointers);
    free(infected_nodes);
    for (int i = 0; i < num_active_links; i++) {
        free(active_links[i]);
    }
    free(active_links);
    free(active_link_positions);

    return 0;
}
