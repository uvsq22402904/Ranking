#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> // Pour fabs

// Utiliser double pour une meilleure précision 

typedef double proba;
typedef int indice;

indice N = 0; // Nombre de nœuds (lignes/colonnes)
indice M = 0; // Nombre d'éléments non nuls (liens)
proba alpha[] = {0.85, 0.90, 0.95, 0.99};  // Facteur de téléportation
proba pourcentage_suppression[] = {0.0, 0.1, 0.2}; // Pourcentage de suppression des liens 
proba sigma = 1e-9; // Critère de convergence (epsilon) - plus strict pour double
indice *est_dangling = NULL; // Tableau booléen : 1 si le nœud est un dangling node

// Structure pour les éléments de matrice creuse (triplets)
struct elem {
    indice i, j; // Lien de i vers j
    proba val;  // Probabilité (sera normalisée à 1/degre_sortant(i))
};

struct elem *p = NULL; // Tableau des éléments non nuls
proba *x = NULL;       // Vecteur PageRank courant
proba *y = NULL;       // Vecteur PageRank suivant (vecteur de travail)
indice *degre_sortant = NULL; // Stocke le degré sortant de chaque nœud




// --- Déclarations de Fonctions ---
void lire_et_normaliser_mtx(char *nom_fichier);

// ajout de la fonction pour suppression aléatoire
void lire_et_normaliser_mtx_avec_suppression(char *nom_fichier, proba pourcentage_suppression);


void initialiser_vecteurs_pagerank();
void identifier_noeuds_dangling();

// Modification de la fonction pour retourner le nombre d'itérations et le temps de convergence via des pointeurs
void iteration_puissance(proba alpha, indice *iterations, proba *temps);


proba norme_L1(proba *v1, proba *v2, indice taille);
void recopier(proba *dest, proba *src, indice taille);
void multiplication_pagerank(proba *z, proba *w, proba alpha);

// Fonction pour enregistrer les résultats dans un fichier pour le rapport
void enregistrer_resultats(char *nom_fichier, char *nom_matrice, proba pourcentage_suppression, proba alpha, indice iterations, proba temps, proba somme_pagerank);


void liberer_memoire();

// --- Fonction Principale ---
int main() {
    char *nom_fichier = "Matrix\\wikipedia-20051105.mtx";
    char *fichier_resultats = "wikipedia-20051105_resultats.csv";
    indice nb_suppressions = sizeof(pourcentage_suppression) / sizeof(pourcentage_suppression[0]);
    indice nb_alphas = sizeof(alpha) / sizeof(alpha[0]);

    // Écrire l'en-tête du fichier CSV
    FILE *f = fopen(fichier_resultats, "w");
    if (f) {
        fprintf(f, "Matrice,Pourcentage_Suppression,Alpha,Iterations,Temps,Somme_PageRank\n");
        fclose(f);
    }

    for (indice i = 0; i < nb_suppressions; i++) {
        for (indice j = 0; j < nb_alphas; j++) {
            printf("\n==> Suppression: %.2f, Alpha: %.2f\n", pourcentage_suppression[i], alpha[j]);

            lire_et_normaliser_mtx_avec_suppression(nom_fichier, pourcentage_suppression[i]);
            if (N == 0 || M == 0 || p == NULL) {
                fprintf(stderr, "Erreur lors de la lecture et normalisation de la matrice.\n");
                continue;
            }

            printf("Matrice : %d x %d, Éléments non nuls : %d\n", N, N, M);

            initialiser_vecteurs_pagerank();
            if (x == NULL || y == NULL) {
                fprintf(stderr, "Erreur d'allocation mémoire pour les vecteurs PageRank.\n");
                liberer_memoire();
                continue;
            }

            identifier_noeuds_dangling();
            if (est_dangling == NULL) {
                fprintf(stderr, "Erreur d'allocation mémoire pour est_dangling.\n");
                liberer_memoire();
                continue;
            }

            indice iterations;
            proba temps;
            iteration_puissance(alpha[j], &iterations, &temps);

            printf("\nScores de PageRank (premiers 20 ou moins):\n");
            indice nb_affichage = (N < 20) ? N : 20;
            for (indice k = 0; k < nb_affichage; k++) {
                printf("Noeud %d : %.8e\n", k, x[k]);
            }
            if (N > nb_affichage) {
                printf("...\n");
            }

            proba somme = 0.0;
            for (indice k = 0; k < N; k++) {
                somme += x[k];
            }
            printf("Somme du PageRank = %.15f\n", somme);

            // Enregistrer les résultats
            enregistrer_resultats(fichier_resultats, nom_fichier, pourcentage_suppression[i], alpha[j], iterations, temps, somme);

            liberer_memoire();
        }
    }

    return 0;
}




void lire_et_normaliser_mtx(char *nom_fichier) {
    FILE *f;
    char ligne[256];
    indice L; 

    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        fprintf(stderr, "Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }


    while (fgets(ligne, sizeof(ligne), f)) {
        if (ligne[0] != '%') {
            break;
        }
    }


    if (sscanf(ligne, "%d %*d %d", &N, &L) != 2) {
         rewind(f);
         while (fgets(ligne, sizeof(ligne), f)) {
            if (ligne[0] != '%') {
                break;
            }
         }
         if (sscanf(ligne, "%d %d %d", &N, &N, &L) != 3) {
            fprintf(stderr, "Erreur: Impossible de lire les dimensions et le nombre d'éléments depuis l'en-tête.\n");
            fclose(f);
            exit(1);
         }
    }
    M = L; 

    if (N <= 0 || M <= 0) {
        fprintf(stderr, "Erreur: Dimensions invalides N=%d, M=%d.\n", N, M);
        fclose(f);
        exit(1);
    }


    p = malloc(M * sizeof(struct elem));
    degre_sortant = calloc(N, sizeof(indice)); 
    if (!p || !degre_sortant) {
        fprintf(stderr, "Erreur d'allocation mémoire pour p ou degre_sortant.\n");
        fclose(f);
        free(p); free(degre_sortant); p = NULL; degre_sortant = NULL; N=0; M=0;
        return;
    }

    indice m_courant = 0;
    indice ligne_elem, col_elem;
    proba val_elem; 
    char buffer[256];

    while (m_courant < M && fgets(buffer, sizeof(buffer), f)) {
        // Tenter de lire avec trois valeurs (i, j, val)
        int items_read = sscanf(buffer, "%d %d %lf", &ligne_elem, &col_elem, &val_elem);
        if (items_read == 2) {
            // Format (i, j) sans valeur, assigner 1.0
            val_elem = 1.0;
        } else if (items_read != 3) {
            fprintf(stderr, "Erreur de lecture à l'élément %d: format invalide.\n", m_courant + 1);
            fclose(f);
            free(p); free(degre_sortant); p = NULL; degre_sortant = NULL; N=0; M=0;
            return;
        }

        p[m_courant].i = ligne_elem - 1;
        p[m_courant].j = col_elem - 1;
        p[m_courant].val = 1.0; 

        if (p[m_courant].i < 0 || p[m_courant].i >= N || p[m_courant].j < 0 || p[m_courant].j >= N) {
             fprintf(stderr, "Attention: Index (%d, %d) hors limites [0, %d) lu à l'élément %d.\n",
                     ligne_elem, col_elem, N, m_courant + 1);

        } else {
            degre_sortant[p[m_courant].i]++;
        }
        m_courant++;
    }
    fclose(f);

    if (m_courant != M) {
        fprintf(stderr, "Attention: Nombre d'éléments lus (%d) différent de M annoncé (%d).\n", m_courant, M);
        M = m_courant; 
        struct elem * p_realloc = realloc(p, M * sizeof(struct elem));
        if (M > 0 && p_realloc == NULL){
            fprintf(stderr, "Erreur realloc p.\n");
            free(p); free(degre_sortant); p = NULL; degre_sortant = NULL; N=0; M=0;
            return;
        }
        p = p_realloc;
    }


    for (indice k = 0; k < M; k++) {
        indice i = p[k].i;
        if (degre_sortant[i] > 0) {
            p[k].val = 1.0 / (proba)degre_sortant[i];
        } else {

            p[k].val = 0.0;
        }
    }

}

void lire_et_normaliser_mtx_avec_suppression(char *nom_fichier, proba pourcentage_suppression) {
    FILE *f;
    char ligne[256];
    indice L;
    
    srand(time(NULL));
    
    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        fprintf(stderr, "Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    // Passer les commentaires
    while (fgets(ligne, sizeof(ligne), f)) {
        if (ligne[0] != '%') {
            break;
        }
    }

    // Lire les dimensions
    if (sscanf(ligne, "%d %*d %d", &N, &L) != 2) {
         rewind(f);
         while (fgets(ligne, sizeof(ligne), f)) {
            if (ligne[0] != '%') {
                break;
            }
         }
         if (sscanf(ligne, "%d %d %d", &N, &N, &L) != 3) {
            fprintf(stderr, "Erreur: Impossible de lire les dimensions et le nombre d'éléments depuis l'en-tête.\n");
            fclose(f);
            exit(1);
         }
    }
    
    // Allouer de la mémoire pour stocker tous les éléments lus
    struct elem *temp_p = malloc(L * sizeof(struct elem));
    degre_sortant = calloc(N, sizeof(indice));
    if (!temp_p || !degre_sortant) {
        fprintf(stderr, "Erreur d'allocation mémoire.\n");
        fclose(f);
        free(temp_p); free(degre_sortant);
        temp_p = NULL; degre_sortant = NULL; N=0; M=0;
        return;
    }

    // Lire tous les éléments du fichier
    indice elements_lus = 0;
    indice ligne_elem, col_elem;
    proba val_elem;
    char buffer[256];

    while (elements_lus < L && fgets(buffer, sizeof(buffer), f)) {
        // Tenter de lire avec trois valeurs (i, j, val)
        int items_read = sscanf(buffer, "%d %d %lf", &ligne_elem, &col_elem, &val_elem);
        if (items_read == 2) {
            // Format (i, j) sans valeur, assigner 1.0
            val_elem = 1.0;
        } else if (items_read != 3) {
            fprintf(stderr, "Erreur de lecture à l'élément %d: format invalide.\n", elements_lus + 1);
            fclose(f);
            free(temp_p); free(degre_sortant);
            temp_p = NULL; degre_sortant = NULL; N=0; M=0;
            return;
        }
        indice i = ligne_elem - 1;
        indice j = col_elem - 1;
        
        if (i < 0 || i >= N || j < 0 || j >= N) {
            fprintf(stderr, "Attention: Index (%d, %d) hors limites [0, %d) lu à l'élément %d.\n",
                    ligne_elem, col_elem, N, elements_lus + 1);
        } else {
            // Décider aléatoirement si on garde cet arc
            proba u = (proba)rand() / RAND_MAX;
            
            if (u >= pourcentage_suppression) {
                // On garde l'arc
                temp_p[elements_lus].i = i;
                temp_p[elements_lus].j = j;
                temp_p[elements_lus].val = 1.0;
                degre_sortant[i]++;
                elements_lus++;
            }
            // Sinon, on saute l'arc (suppression)
        }
    }
    fclose(f);
    
    // Mettre à jour M avec le nombre réel d'éléments conservés
    M = elements_lus;
    printf("Arcs après suppression aléatoire de %.1f%% : %d (sur %d initialement)\n", 
           pourcentage_suppression * 100, M, L);
    
    // Réallouer p à la taille exacte nécessaire
    p = malloc(M * sizeof(struct elem));
    if (!p) {
        fprintf(stderr, "Erreur d'allocation mémoire pour p.\n");
        free(temp_p); free(degre_sortant);
        p = NULL; degre_sortant = NULL; N=0; M=0;
        return;
    }
    
    // Copier les éléments de temp_p vers p
    for (indice k = 0; k < M; k++) {
        p[k] = temp_p[k];
    }
    free(temp_p);
    
    // Normaliser les valeurs
    for (indice k = 0; k < M; k++) {
        indice i = p[k].i;
        if (degre_sortant[i] > 0) {
            p[k].val = 1.0 / (proba)degre_sortant[i];
        } else {
            p[k].val = 0.0;
        }
    }
}

void initialiser_vecteurs_pagerank() {
    x = malloc(N * sizeof(proba));
    y = malloc(N * sizeof(proba));
    if (!x || !y) {
        free(x); free(y); x = NULL; y = NULL;
        return;
    }
    proba rang_initial = 1.0 / (proba)N;
    for (indice i = 0; i < N; i++) {
        x[i] = rang_initial;
    }
}

void identifier_noeuds_dangling() {
    est_dangling = malloc(N * sizeof(indice)); 
    if (!est_dangling) {
        free(degre_sortant); degre_sortant = NULL;
        return;
    }
    for (indice i = 0; i < N; i++) {
        est_dangling[i] = (degre_sortant[i] == 0);
    }

    free(degre_sortant);
    degre_sortant = NULL;
}

// Norme L1: somme(|v1[i] - v2[i]|)
proba norme_L1(proba *v1, proba *v2, indice taille) {
    proba somme = 0.0;
    for (indice i = 0; i < taille; i++) {
        somme += fabs(v1[i] - v2[i]); 
    }
    return somme;
}

void recopier(proba *dest, proba *src, indice taille) {
    for (indice i = 0; i < taille; i++) {
        dest[i] = src[i];
    }
}

// Effectue une étape: w = alpha * P * z + contrib_dangling + contrib_teleport
void multiplication_pagerank(proba *z, proba *w, proba alpha) {

    for (indice i = 0; i < N; i++) {
        w[i] = 0.0;
    }

    for (indice k = 0; k < M; k++) {

        w[p[k].j] += z[p[k].i] * p[k].val;
    }

    proba somme_dangling = 0.0;
    for (indice i = 0; i < N; i++) {
        if (est_dangling[i]) {
            somme_dangling += z[i];
        }
    }


    proba dist_dangling_par_noeud = 0.0;
    if (N > 0) {
         dist_dangling_par_noeud = somme_dangling / (proba)N;
    }
    proba teleport_par_noeud = (1.0 - alpha) / (proba)N;

    for (indice i = 0; i < N; i++) {
        w[i] = alpha * (w[i] + dist_dangling_par_noeud) + teleport_par_noeud;
    }


}

void iteration_puissance(proba alpha, indice *iterations, proba *temps) {
    indice k = 0;
    proba diff = 1.0;
    struct timeval t1, t2;

    gettimeofday(&t1, NULL);

    while (diff > sigma) {
        k++;
        multiplication_pagerank(x, y, alpha);
        diff = norme_L1(x, y, N);
        recopier(x, y, N);

        if (k % 10 == 0) {
            printf("Itération %d, Diff = %e\n", k, diff);
        }
        if (k > 10000) {
            printf("Attention: Dépassement de 10000 itérations.\n");
            break;
        }
        if (isnan(diff) || isinf(diff)) {
            fprintf(stderr, "Erreur: Convergence échouée (NaN ou Inf détecté) à l'itération %d.\n", k);
            break;
        }
    }

    gettimeofday(&t2, NULL);
    *temps = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
    *iterations = k;

    printf("Convergence en %d itérations (Diff = %e), temps total = %.4f sec\n", k, diff, *temps);
}

void enregistrer_resultats(char *nom_fichier, char *nom_matrice, proba pourcentage_suppression, proba alpha, indice iterations, proba temps, proba somme_pagerank) {
    FILE *f = fopen(nom_fichier, "a");
    if (f == NULL) {
        fprintf(stderr, "Erreur d'ouverture du fichier de résultats '%s'\n", nom_fichier);
        return;
    }
    fprintf(f, "%s,%.2f,%.2f,%d,%.4f,%.15f\n", nom_matrice, pourcentage_suppression, alpha, iterations, temps, somme_pagerank);
    fclose(f);
}

void liberer_memoire() {
    free(x);
    free(y);
    free(p);
    free(est_dangling); 

    x = NULL;
    y = NULL;
    p = NULL;
    est_dangling = NULL;
}