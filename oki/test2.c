#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>  // Pour fabs()

typedef float proba;
typedef int indice;

proba alpha = 0.85;  // Facteur de téléportation
proba sigma = 1e-6;  // Critère de convergence (epsilon)

struct elem {
    indice i, j;
    proba val;
};

struct matrice {
    indice N;
    proba **V;
};

// Fonction de norme ||x - y||
float norme(proba *x, proba *y, indice taille) {
    float somme = 0.0;
    for (indice i = 0; i < taille; i++) {
        somme += fabs(x[i] - y[i]);  // Utilisation de fabs() pour la sûreté
    }
    return somme;
}

// Initialiser le vecteur x avec 1/N
void initx(proba *x, indice taille) {
    for (indice i = 0; i < taille; i++) {
        x[i] = 1.0 / taille;
    }
}

void mult(proba *x, struct elem *p, proba *y, indice M, indice C) {
    indice i, j, k;
    
    // Réinitialisation du tableau y à zéro
    for (i = 0; i < C; i++) {
        y[i] = 0.0;
    }

    // Détection des dangling nodes (pages sans liens sortants)
    int *has_outlinks = calloc(C, sizeof(int));
    for (k = 0; k < M; k++) {
        has_outlinks[p[k].i] = 1;
    }

    // Calcul de la masse des dangling nodes
    proba dangling_mass = 0.0;
    for (i = 0; i < C; i++) {
        if (!has_outlinks[i]) {
            dangling_mass += x[i];
        }
    }
    free(has_outlinks);

    // Multiplication matrice-vecteur (mise à jour de PageRank)
    for (k = 0; k < M; k++) {
        i = p[k].i;
        j = p[k].j;
        y[j] += alpha * x[i] * p[k].val;
    }

    // Ajout de la masse des dangling nodes (redistribuée uniformément)
    proba dangling_contribution = alpha * dangling_mass / C;
    for (i = 0; i < C; i++) {
        y[i] += dangling_contribution;
    }

    // Facteur de téléportation
    proba somme_x = 0.0;
    for (i = 0; i < C; i++) {
        somme_x += x[i];
    }

    proba teleporte = (1.0 - alpha) * somme_x / C;
    
    // Ajout du facteur de téléportation
    for (i = 0; i < C; i++) {
        y[i] += teleporte;
    }

    // 🔥 **Normalisation pour s'assurer que sum(y) == 1**
    proba somme_y = 0.0;
    for (i = 0; i < C; i++) {
        somme_y += y[i];
    }
    if (somme_y > 0) {  // Éviter une division par zéro
        for (i = 0; i < C; i++) {
            y[i] /= somme_y;
        }
    }
}


void recopie(proba *x, proba *y, indice taille) {
    for (indice i = 0; i < taille; i++) {
        x[i] = y[i];
    }
}

float une_iteration(proba *x, struct elem *p, proba *y, indice M, indice C) {
    mult(x, p, y, M, C);
    float s = norme(x, y, C);
    recopie(x, y, C);
    return s;
}

void iterer(proba *x, struct elem *p, proba *y, indice M, indice C) {
    indice k = 0;
    float s = 1.0;
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    while (s > sigma) {
        k++;
        s = une_iteration(x, p, y, M, C);
    }

    gettimeofday(&t2, NULL);
    printf("Convergence en %d itérations, temps total = %f s\n", k, 
           (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6);
}

// Lecture du fichier MatrixMarket
void lire_matrice_mtx(char *nom_fichier, struct elem **p, indice *C, indice *M) {
    FILE *f = fopen(nom_fichier, "r");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    char ligne[256];

    // Ignorer les lignes de commentaire
    while (fgets(ligne, sizeof(ligne), f)) {
        if (ligne[0] != '%') break;
    }

    // Lire les dimensions
    indice taille, nb_elements;
    sscanf(ligne, "%d %d %d", &taille, &taille, &nb_elements);
    
    *C = taille;
    *M = nb_elements;

    *p = malloc(nb_elements * sizeof(struct elem));
    if (*p == NULL) {
        printf("Erreur d'allocation mémoire pour la matrice creuse\n");
        exit(1);
    }

    // Lecture des éléments
    for (indice i = 0; i < nb_elements; i++) {
        indice ligne_elem, col_elem;
        float val_elem;

        fscanf(f, "%d %d %f", &ligne_elem, &col_elem, &val_elem);

        (*p)[i].i = ligne_elem - 1;  // 1-based to 0-based
        (*p)[i].j = col_elem - 1;
        (*p)[i].val = val_elem;
    }

    fclose(f);
    // Normalisation des liens sortants
    proba *sum_out = calloc(*C, sizeof(proba));
    for (indice k = 0; k < *M; k++) {
        sum_out[(*p)[k].i] += (*p)[k].val;
    }

    for (indice k = 0; k < *M; k++) {
        if (sum_out[(*p)[k].i] > 0) {
            (*p)[k].val /= sum_out[(*p)[k].i];
        }
    }

    free(sum_out);

}

//int main() {
//    char *nom_fichier = "small_test_matrix.mtx";
//    indice C, M;
//    struct elem *p = NULL;
//
//    lire_matrice_mtx(nom_fichier, &p, &C, &M);
//
//    printf("Matrice : %d x %d, Éléments non nuls : %d\n", C, C, M);
//
//    // Allocation mémoire pour PageRank
//    proba *x = malloc(C * sizeof(proba));
//    proba *y = malloc(C * sizeof(proba));
//
//    if (!x || !y) {
//        printf("Erreur d'allocation mémoire\n");
//        free(p);
//        return 1;
//    }
//
//    // Initialisation
//    initx(x, C);
//
//    // Exécution de l'algorithme PageRank
//    iterer(x, p, y, M, C);
//
//    // Affichage des scores
//    printf("\nScores de PageRank :\n");
//    for (indice i = 0; i < C; i++) {
//        if (x[i] > 0.001) {
//            printf("Nœud %d : %.6f\n", i + 1, x[i]);
//        }
//    }
//
//    // Libération de la mémoire
//    free(x);
//    free(y);
//    free(p);
//
//    return 0;
//}