#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

typedef float proba;
typedef int indice;

#define abs(x) ((x) > 0 ? (x) : -(x))

indice L=0,C=0,M=0;//L cest le nbr de lignes C cest le nbr de colonnes M cest le nbr de valeur non nulles
proba alpha = 0.85;  // Facteur de téléportation
proba sigma = 1e-6;  // Critère de convergence (epsilon)
indice *est_dangling = NULL; // Tableau booléen : 1 si le nœud est un dangling node


struct elem{
    indice i,j;
    proba val;
};

struct elem *p;
proba *x;
proba *y;




struct matrice {
    indice N;
    proba **V;
};



void recopie(proba *x, proba *y, indice taille){
    for (indice i = 0; i < taille; i++) {
        x[i] = y[i];
    }
}


void mult(proba *x, struct elem *p, proba *y, indice taille, indice M) {
    // Initialize y to zeros
    for (indice i = 0; i < taille; i++) {
        y[i] = 0.0;
    }

    // 1. Multiplier les contributions normales (liens sortants existants)
    for (indice k = 0; k < M; k++) {
        y[p[k].j] += alpha * x[p[k].i] * p[k].val;
    }

    // 2. Ajouter l'effet des nœuds dangling (distribuent leur poids à tous les nœuds)
    proba somme_dangling = 0.0;
    for (indice i = 0; i < taille; i++) {
        if (est_dangling[i]) {
            somme_dangling += x[i];
        }
    }
    
    // Distribute dangling node weight evenly
    if (taille > 0 && somme_dangling > 0) {
        proba contribution_dangling = alpha * somme_dangling / taille;
        for (indice i = 0; i < taille; i++) {
            // Remove the debug print that's causing issues
            y[i] += contribution_dangling;
        }
    }

    // 3. Ajouter la partie téléportation (facteur 1-alpha)
    proba teleporte = (1.0 - alpha) / taille;
    for (indice i = 0; i < taille; i++) {
        y[i] += teleporte;  // Each node gets equal teleportation probability
    }
    
    // Ensure the sum is 1.0 (optional normalization step)
    proba sum = 0.0;
    for (indice i = 0; i < taille; i++) {
        sum += y[i];
    }
    if (sum > 0) {
        for (indice i = 0; i < taille; i++) {
            y[i] /= sum;
        }
    }
}


// Fonction de norme ||x - y||
float norme(proba *x, proba *y, indice taille){
    float somme = 0.0;
    for (indice i = 0; i < taille; i++) {
        somme += abs(x[i] - y[i]);
    }
    return somme;
}


void lire_matrice_mtx(char *nom_fichier, struct elem *p, struct matrice *t) {
    FILE *f;
    char ligne[256];
    int lignes_entete = 0;
    
    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    // Ignorer les lignes de commentaire et d'en-tête
    while (fgets(ligne, sizeof(ligne), f)) {
        if (ligne[0] != '%') {
            break;
        }
        lignes_entete++;
    }

    // Lecture des dimensions et du nombre d'éléments
    sscanf(ligne, "%d %d %d", &t->N, &t->N, &M);
    
    // Initialisation de la matrice
    // Modification : allocation manuelle plutôt que copie
    t->V = malloc(t->N * sizeof(proba*));
    for (int i = 0; i < t->N; i++) {
        t->V[i] = calloc(t->N, sizeof(proba));
    }

    // Variables pour la lecture des triplets
    int ligne_elem, col_elem;
    float val_elem;
    int index_triplet = 0;

    // Lecture des éléments
    while (fscanf(f, "%d %d %f", &ligne_elem, &col_elem, &val_elem) == 3) {
        // Ajustement des indices (1-based to 0-based)
        ligne_elem--;
        col_elem--;

        // Remplissage de la matrice dense
        t->V[ligne_elem][col_elem] = val_elem;

        // Remplissage des triplets
        p[index_triplet].i = ligne_elem;
        p[index_triplet].j = col_elem;
        p[index_triplet].val = val_elem;
        
        index_triplet++;
    }

    fclose(f);
}


// Initialiser le vecteur x avec 1/N
void initx(proba *x, indice taille) {
    for (indice i = 0; i < taille; i++) {
        x[i] = 1.0 / taille;
    }
}

float une_iteration(proba *x, struct elem *p, proba *y, indice taille, indice M){
    mult(x, p, y, taille, M);
    float s = norme(x, y, taille);
    recopie(x, y, taille);
    return s;
}

void iterer(proba *x, struct elem *p, proba *y, indice taille, indice M){
    indice k = 0;
    float s = 1.0;
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    while (s > sigma) {
        k++;
        s = une_iteration(x, p, y, taille, M);
    }
    gettimeofday(&t2, NULL);
    printf("Convergence en %d iterations, temps total = %f\n", k,
           (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6);
}



void detecter_dangling_noeuds(struct matrice *matrice, int *est_dangling, indice taille) {
    // Reset all nodes to not dangling
    for (indice i = 0; i < taille; i++) {
        est_dangling[i] = 1;  // Assume all are dangling initially
    }
    
    // For each non-zero element, mark its source node as non-dangling
    for (indice i = 0; i < taille; i++) {
        for (indice j = 0; j < taille; j++) {
            if (matrice->V[i][j] != 0.0) {
                est_dangling[i] = 0;  // Node i has an outgoing link, not dangling
                break;
            }
        }
    }
}



int main() {
    // Nom du fichier en dur dans le code
    char *nom_fichier = "webbase-1M.mtx";

    // Lire le fichier MatrixMarket
    struct matrice t;
    t.N = 0;  // Initialiser à 0 pour que la fonction puisse définir la taille

    // Calculer le nombre d'éléments non nuls
    FILE *f = fopen(nom_fichier, "r");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        return 1;
    }

    char ligne[256];
    int lignes_entete = 0;
    while (fgets(ligne, sizeof(ligne), f)) {
        if (ligne[0] != '%') {
            break;
        }
        lignes_entete++;
    }

    int taille, taille2, nb_elements;
    sscanf(ligne, "%d %d %d", &taille, &taille2, &nb_elements);
    fclose(f);

    // Définir les variables globales
    L = C = taille;
    M = nb_elements;

    printf("Matrice : %d x %d, Éléments non nuls : %d\n", L, C, M);

    // Allouer la mémoire
    x = malloc(C * sizeof(proba));
    y = malloc(C * sizeof(proba));
    p = malloc(M * sizeof(struct elem));

    if (!x || !y || !p) {
        printf("Erreur d'allocation mémoire\n");
        return 1;
    }

    // Lire la matrice
    lire_matrice_mtx(nom_fichier, p, &t);

    // Détecter les nœuds sans lien sortant
    est_dangling = malloc(C * sizeof(indice));
    if (!est_dangling) {    
        printf("Erreur d'allocation mémoire pour les dangling nodes\n");
        return 1;
    }
    detecter_dangling_noeuds(&t, est_dangling, C);

    


    // Initialiser le vecteur de PageRank
    initx(x, C);

    // Effectuer l'algorithme PageRank
    iterer(x, p, y, L, M);

    // Afficher les résultats (optionnel)
    printf("\nScores de PageRank :\n");
    for (int i = 0; i < C; i++) {
        printf("Nœud %d : %.6f\n", i+1, x[i]);
    }

    double sum = 0;
    for (int i = 0; i < C; i++) sum += x[i];
    printf("Somme du PageRank = %.10f\n", sum);


    // Libérer la mémoire
    free(x); 
    free(y); 
    free(p);
    for (indice i = 0; i < t.N; i++) {
        free(t.V[i]);
    }
    free(t.V);
    free(est_dangling);
    return 0;
}