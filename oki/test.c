#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Inclure directement le fichier d'implémentation
#include "test2.c"

void test_mult() {
    struct elem test_matrix[] = {
        {0, 1, 0.5},  
        {1, 2, 0.5},  
        {2, 0, 1.0}   
    };
    
    indice M = 3;  
    indice C = 3;  
    
    proba x[3] = {1.0/3, 1.0/3, 1.0/3};  
    proba y[3] = {0};  
    
    mult(x, test_matrix, y, M, C);
    
    printf("Test mult() - Résultats :\n");
    for (int i = 0; i < C; i++) {
        printf("y[%d] = %f\n", i, y[i]);
    }
}

void test_pagerank_sum() {
    char *nom_fichier = "webbase-1M.mtx";
    indice C, M;
    struct elem *p = NULL;

    lire_matrice_mtx(nom_fichier, &p, &C, &M);

    printf("Matrice : %d x %d, Éléments non nuls : %d\n", C, C, M);

    proba *x = malloc(C * sizeof(proba));
    proba *y = malloc(C * sizeof(proba));

    if (!x || !y) {
        printf("Erreur d'allocation mémoire\n");
        free(p);
        return;
    }

    initx(x, C);
    iterer(x, p, y, M, C);

    proba somme_totale = 0.0;
    for (indice i = 0; i < C; i++) {
        somme_totale += x[i];
    }

    printf("\nSomme totale des PageRanks : %f\n", somme_totale);
    
    if (fabs(somme_totale - 1.0) < 0.001) {
        printf("✅ Somme CORRECTE (proche de 1)\n");
    } else {
        printf("❌ Somme INCORRECTE (doit être proche de 1)\n");
    }

    //printf("\nScores de PageRank :\n");
    //for (indice i = 0; i < C; i++) {
    //    printf("Nœud %d : %.6f\n", i + 1, x[i]);
    //}

    free(x);
    free(y);
    free(p);
}

int main() {
    printf("Test de la fonction mult() :\n");
    test_mult();

    printf("\n--- Test PageRank Sum ---\n");
    test_pagerank_sum();

    return 0;
}