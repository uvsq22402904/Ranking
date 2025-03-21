#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DAMPING 0.85
#define EPSILON 1e-6
#define MAX_ITER 1000

typedef double proba;
typedef int indice;

indice N, M;  // Number of nodes and edges
struct elem *P;
proba *x, *y, *dangling;  // Changed 'f' to 'dangling'

struct elem {
    indice i, j;
    proba val;
};

// Allocate memory for the structures
void allocateMemory(indice N, indice M) {
    P = malloc(M * sizeof(struct elem));
    x = malloc(N * sizeof(proba));
    y = malloc(N * sizeof(proba));
    dangling = calloc(N, sizeof(proba)); // Initialize to 0
    if (!P || !x || !y || !dangling) {
        perror("Memory allocation failed.");
        exit(EXIT_FAILURE);
    }
}

// Read the Matrix Market file (.mtx format)
void lectureFichier(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        exit(EXIT_FAILURE);
    }
    
    char line[256];
    do {
        fgets(line, sizeof(line), file);
    } while (line[0] == '%');  // Skip comments in .mtx format

    sscanf(line, "%d %d %d", &N, &N, &M); // Square matrix
    allocateMemory(N, M);

    for (indice k = 0; k < M; k++) {
        fscanf(file, "%d %d %lf", &P[k].i, &P[k].j, &P[k].val);
        P[k].i;  // Convert to zero-based indexing
        P[k].j;

        // Mark nodes that have outbound links
        dangling[P[k].i] = 1;  
    }
    fclose(file);
}

// Set y to zero
void mettre_a_zero(proba *y) {
    for (indice i = 0; i < N; i++) {
        y[i] = 0.0;
    }
}

// Multiply x * P (sparse matrix multiplication)
void mult(proba *x, struct elem *P, proba *y) {
    mettre_a_zero(y);
    for (indice k = 0; k < M; k++) {
        y[P[k].j] += x[P[k].i] * P[k].val;
    }
}

// Compute norm for convergence check
double norme_diff(proba *x, proba *y) {
    double sum = 0.0;
    for (indice i = 0; i < N; i++) {
        sum += fabs(x[i] - y[i]);
    }
    return sum;
}

// Initialize x(0) = 1/N
void init_x() {
    for (indice i = 0; i < N; i++) {
        x[i] = 1.0 / N;
    }
}

// Apply damping factor and handle dangling nodes
void apply_damping() {
    double S = 0.0;  // Sum of contributions from dangling nodes
    for (indice i = 0; i < N; i++) {
        if (!dangling[i]) S += x[i];  
    }
    S = (DAMPING / N) * S;

    for (indice i = 0; i < N; i++) {
        y[i] = DAMPING * y[i] + (1.0 - DAMPING) / N + S;
    }
}

void recopie(proba *x, proba *y){
    for (indice i = 0; i < N; i++){
        x[i] = y[i];
    }
}

// Compute PageRank
void pageRank() {
    int iter;
    init_x();
    for (iter = 0; iter < MAX_ITER; iter++) {
        mult(x, P, y);
        apply_damping();
        if (norme_diff(x, y) < EPSILON) break;
        recopie(x, y);
    }
}

// Display results
void afficher_resultats() {
    printf("PageRank values:\n");
    float sum = 0;
    for (indice i = 0; i < N; i++) {
        printf("Page %d: %lf\n", i + 1, x[i]);
        sum += x[i];
    }
    printf("Sum: %f\n", sum);
}



// Main function
int main() {
    lectureFichier("matrices/yanis.txt");
    pageRank();
    afficher_resultats();
    
    // Free memory
    free(P);
    free(x);
    free(y);
    free(dangling);

}