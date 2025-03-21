#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define C 6  // Nombre de pages (colonnes et lignes de la matrice M)
#define alpha 0.85  // Facteur de téléportation
#define EPSILON 1e-6  // Seuil de convergence
#define MAX_ITER 100  // Nombre maximal d'itérations

typedef double proba;  
typedef int indice;

void initx(proba *x, indice taille) {
    for (indice i = 0; i < taille; i++) {
        x[i] = 1.0 / taille;
    }
}

void multiplication(proba M[C][C], proba *x, proba *y) {
    proba somme = 0.0;

    // Initialisation de y à zéro
    for (int i = 0; i < C; i++) {
        y[i] = 0.0;
    }

    // Multiplication de la matrice M par le vecteur x
    for (int i = 0; i < C; i++) {
        for (int j = 0; j < C; j++) {
            y[i] += M[i][j] * x[j];
        }
        somme += y[i];
    }

    // Ajout du facteur de téléportation
    proba teleporte = (1.0 - alpha) / C;
    for (int i = 0; i < C; i++) {
        y[i] = alpha * y[i] + teleporte;
    }

    // Normalisation
    somme = 0.0;
    for (int i = 0; i < C; i++) {
        somme += y[i];
    }

    if (fabs(somme - 1.0) > EPSILON) {
        for (int i = 0; i < C; i++) {
            y[i] /= somme;
        }
    }
}

double erreur(proba *x, proba *y) {
    double maxDiff = 0.0;
    for (int i = 0; i < C; i++) {
        double diff = fabs(x[i] - y[i]);
        if (diff > maxDiff) {
            maxDiff = diff;
        }
    }
    return maxDiff;
}

void affiche_vecteur(proba *v) {
    for (int i = 0; i < C; i++) {
        printf("%.6f ", v[i]);
    }
    printf("\n");
}

void pagerank(proba M[C][C]) {
    proba x[C], y[C];
    initx(x, C);

    printf("Iteration 0 : ");
    affiche_vecteur(x);

    int iter = 0;
    while (iter < MAX_ITER) {
        multiplication(M, x, y);
        double diff = erreur(x, y);

        printf("Iteration %d : ", iter + 1);
        affiche_vecteur(y);

        if (diff < EPSILON) {
            break;
        }

        // Mise à jour de x pour la prochaine itération
        for (int i = 0; i < C; i++) {
            x[i] = y[i];
        }

        iter++;
    }
}

int main() {
    proba M[C][C] = {
        {0.25, 0.25, 0, 0.25, 0, 0.25},
        {0, 0, 0.25, 0.25, 0.25, 0.25},
        {0.25, 0.25, 0, 0.25, 0.25, 0},
        {0, 0, 0, 0, 1, 0},
        {0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0} // Page morte
    };

    pagerank(M);

    return 0;
}
