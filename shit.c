#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> // Pour fabs

// Utiliser double pour une meilleure précision
typedef double proba;
typedef int indice;

typedef struct {
    int *groupes;       // Affectation des nœuds aux groupes
    int nb_groupes;     // Nombre total de groupes
    double connectivite;// Mesure de décomposabilité
} structure_NCD;

structure_NCD sncd;

indice N = 0; // Nombre de nœuds (lignes/colonnes)
indice M = 0; // Nombre d'éléments non nuls (liens)
proba alpha = 0.85;  // Facteur de téléportation
proba sigma = 1e-9; // Critère de convergence (epsilon) - plus strict pour double
int *est_dangling = NULL; // Tableau booléen : 1 si le nœud est un dangling node

// Structure pour les éléments de matrice creuse (triplets)
struct elem {
    indice i, j; // Lien de i vers j
    proba val;  // Probabilité (sera normalisée à 1/degre_sortant(i))
};

struct elem *p = NULL; // Tableau des éléments non nuls
proba *x = NULL;       // Vecteur PageRank courant
proba *y = NULL;       // Vecteur PageRank suivant (vecteur de travail)
int *degre_sortant = NULL; // Stocke le degré sortant de chaque nœud

proba pourcentage_suppression = 0.1; // Pourcentage de suppression des liens 

// --- Déclarations de Fonctions ---
void lire_et_normaliser_mtx(char *nom_fichier);
void lire_et_normaliser_mtx_avec_suppression_NCD(char *nom_fichier);
void initialiser_vecteurs_pagerank();
void identifier_noeuds_dangling();
void iteration_puissance();
proba norme_L1(proba *v1, proba *v2, indice taille);
void recopier(proba *dest, proba *src, indice taille);
void multiplication_pagerank(proba *z, proba *w);
void liberer_memoire();
void partitionner_graphe_NCD();
void analyser_connectivite_NCD();
void tester_convergence_alpha_proche_1();
void sauvegarder_resultats(const char *fichier);
void lire_et_normaliser_mtx(char *nom_fichier) {
    FILE *f;
    char ligne[256];
    int L; 

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
    degre_sortant = calloc(N, sizeof(int)); 
    if (!p || !degre_sortant) {
        fprintf(stderr, "Erreur d'allocation mémoire pour p ou degre_sortant.\n");
        fclose(f);
        free(p); free(degre_sortant); p = NULL; degre_sortant = NULL; N=0; M=0;
        return;
    }

    indice m_courant = 0;
    int ligne_elem, col_elem;
    double val_elem; 
    while (m_courant < M && fscanf(f, "%d %d %lf", &ligne_elem, &col_elem, &val_elem) == 3) {
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

// --- Fonction Principale ---
int main() {
    // Initialisation de la structure NCD
        char *nom_fichier = "webbase-1M copy.mtx";
    
        // First pass: read matrix to get degrees
        lire_et_normaliser_mtx(nom_fichier);
        if (N == 0 || M == 0 || p == NULL) {
            fprintf(stderr, "Erreur lecture initiale de la matrice.\n");
            return 1;
        }
    
        // Cluster based on initial degrees
        partitionner_graphe_NCD();
        
        // Free initial read data
        free(p);
        free(degre_sortant);
        p = NULL;
        degre_sortant = NULL;
        N = 0;
        M = 0;
    
        // Second pass: NCD-aware reading
        lire_et_normaliser_mtx_avec_suppression_NCD(nom_fichier);
        if (N == 0 || M == 0 || p == NULL) {
            fprintf(stderr, "Erreur lecture NCD de la matrice.\n");
            liberer_memoire();
            return 1;
        }
    
        // Only analyze connectivity now (don't cluster again)
        analyser_connectivite_NCD();
    
        printf("Matrice : %d x %d, Éléments non nuls : %d\n", N, N, M);
    
        initialiser_vecteurs_pagerank();
        if (x == NULL || y == NULL) {
            fprintf(stderr, "Erreur d'allocation mémoire pour les vecteurs PageRank.\n");
            liberer_memoire();
            return 1;
        }
    
        identifier_noeuds_dangling();
        if (est_dangling == NULL) {
            fprintf(stderr, "Erreur d'allocation mémoire pour est_dangling.\n");
            liberer_memoire();
            return 1;
        }
    
        // Tests de convergence
        tester_convergence_alpha_proche_1();
    
        printf("\nScores de PageRank (premiers 20 ou moins):\n");
        int nb_affichage = (N < 20) ? N : 20;
        for (indice i = 0; i < nb_affichage; i++) {
            printf("Noeud %d : %.8e\n", i, x[i]); 
        }
        if (N > nb_affichage) {
            printf("...\n");
        }
    
        proba somme = 0.0;
        for (indice i = 0; i < N; i++) {
            somme += x[i];
        }
        printf("Somme du PageRank = %.15f\n", somme);
    
        sauvegarder_resultats("resultats_ncd.csv");
        liberer_memoire();
    
        return 0;
    }

// Partitionnement du graphe pour NCD
void partitionner_graphe_NCD() {
    // Déterminer le nombre de groupes (entre 3 et 10)
    sncd.nb_groupes = (int)fmax(3, fmin(10, N*0.05));
    sncd.groupes = malloc(N * sizeof(int));
    
    // Initialisation aléatoire
    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        sncd.groupes[i] = rand() % sncd.nb_groupes;
    }

    // Simple clustering basé sur les degrés (k-means simplifié)
    double *degrees = calloc(N, sizeof(double));
    for (int k = 0; k < M; k++) {
        degrees[p[k].i]++;
    }

    // Quelques itérations de raffinement
    for (int iter = 0; iter < 5; iter++) {
        // Calcul des centroïdes
        double *centroids = calloc(sncd.nb_groupes, sizeof(double));
        int *counts = calloc(sncd.nb_groupes, sizeof(int));
        
        for (int i = 0; i < N; i++) {
            centroids[sncd.groupes[i]] += degrees[i];
            counts[sncd.groupes[i]]++;
        }
        
        for (int g = 0; g < sncd.nb_groupes; g++) {
            if (counts[g] > 0) centroids[g] /= counts[g];
        }
        
        // Réaffectation des nœuds
        int changes = 0;
        for (int i = 0; i < N; i++) {
            int best_group = 0;
            double min_dist = INFINITY;
            
            for (int g = 0; g < sncd.nb_groupes; g++) {
                double dist = fabs(degrees[i] - centroids[g]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_group = g;
                }
            }
            
            if (sncd.groupes[i] != best_group) {
                sncd.groupes[i] = best_group;
                changes++;
            }
        }
        
        free(centroids);
        free(counts);
        
        if (changes < N*0.01) break; // Arrêt si peu de changements
    }
    
    free(degrees);
}

// Analyse de la connectivité NCD
void analyser_connectivite_NCD() {
    int liens_internes = 0, liens_externes = 0;
    
    for (int k = 0; k < M; k++) {
        int i = p[k].i;
        int j = p[k].j;
        if (sncd.groupes[i] == sncd.groupes[j]) {
            liens_internes++;
        } else {
            liens_externes++;
        }
    }
    
    sncd.connectivite = (liens_externes == 0) ? 0 : (double)liens_internes/liens_externes;
    printf("Connectivité NCD: %.2f (ratio liens internes/externes)\n", sncd.connectivite);
}

// Version NCD de la lecture de matrice
void lire_et_normaliser_mtx_avec_suppression_NCD(char *nom_fichier) {
    FILE *f;
    char ligne[256];
    int L;
    
    srand(time(NULL));
    
    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        fprintf(stderr, "Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    // Passer les commentaires
    while (fgets(ligne, sizeof(ligne), f)) {
        if (ligne[0] != '%') break;
    }

    // Lire les dimensions
    if (sscanf(ligne, "%d %*d %d", &N, &L) != 2) {
        rewind(f);
        while (fgets(ligne, sizeof(ligne), f)) {
            if (ligne[0] != '%') break;
        }
        if (sscanf(ligne, "%d %d %d", &N, &N, &L) != 3) {
            fprintf(stderr, "Erreur: Impossible de lire les dimensions\n");
            fclose(f);
            exit(1);
        }
    }
    
    // Allocation mémoire
    struct elem *temp_p = malloc(L * sizeof(struct elem));
    degre_sortant = calloc(N, sizeof(int));
    if (!temp_p || !degre_sortant) {
        fprintf(stderr, "Erreur d'allocation mémoire.\n");
        fclose(f);
        free(temp_p); free(degre_sortant);
        exit(1);
    }

    // Lire les éléments avec suppression NCD-aware
    indice elements_lus = 0;
    int ligne_elem, col_elem;
    double val_elem;
    
    while (elements_lus < L && fscanf(f, "%d %d %lf", &ligne_elem, &col_elem, &val_elem) == 3) {
        int i = ligne_elem - 1;
        int j = col_elem - 1;
        
        if (i < 0 || i >= N || j < 0 || j >= N) {
            fprintf(stderr, "Index hors limites (%d, %d)\n", ligne_elem, col_elem);
            continue;
        }
        
        double u = (double)rand() / RAND_MAX;
        int meme_groupe = (sncd.groupes[i] == sncd.groupes[j]);
        
        // Supprimer préférentiellement les liens entre groupes
        if (meme_groupe || u >= pourcentage_suppression*2) {
            temp_p[elements_lus].i = i;
            temp_p[elements_lus].j = j;
            temp_p[elements_lus].val = 1.0;
            degre_sortant[i]++;
            elements_lus++;
        }
    }
    fclose(f);
    
    // Mise à jour du nombre d'éléments
    M = elements_lus;
    printf("Arcs après suppression NCD: %d (initial: %d)\n", M, L);
    
    // Réallocation à la taille exacte
    p = realloc(temp_p, M * sizeof(struct elem));
    if (!p && M > 0) {
        fprintf(stderr, "Erreur de réallocation\n");
        free(temp_p); free(degre_sortant);
        exit(1);
    }
    
    // Normalisation
    for (indice k = 0; k < M; k++) {
        indice i = p[k].i;
        if (degre_sortant[i] > 0) {
            p[k].val = 1.0 / (proba)degre_sortant[i];
        } else {
            p[k].val = 0.0;
        }
    }
}

// Test de convergence pour alpha proche de 1
void tester_convergence_alpha_proche_1() {
    proba alphas[] = {0.85, 0.90, 0.95, 0.99, 0.999};
    int nb_tests = sizeof(alphas)/sizeof(proba);
    
    printf("\n=== Test de convergence pour alpha proche de 1 ===\n");
    
    for (int i = 0; i < nb_tests; i++) {
        alpha = alphas[i];
        printf("\nAlpha = %.3f\n", alpha);
        
        // Réinitialiser les vecteurs
        proba init_val = 1.0 / (proba)N;
        for (indice j = 0; j < N; j++) {
            x[j] = init_val;
        }
        
        iteration_puissance();
    }
}

// Sauvegarde des résultats
void sauvegarder_resultats(const char *fichier) {
    FILE *f = fopen(fichier, "w");
    if (!f) {
        perror("Erreur ouverture fichier résultats");
        return;
    }
    
    fprintf(f, "Noeud,PageRank,Groupe\n");
    for (indice i = 0; i < N; i++) {
        fprintf(f, "%d,%.10f,%d\n", i, x[i], sncd.groupes[i]);
    }
    fclose(f);
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
    est_dangling = malloc(N * sizeof(int)); 
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
void multiplication_pagerank(proba *z, proba *w) {

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

void iteration_puissance() {
    indice k = 0;
    proba diff = 1.0; 
    clock_t t1, t2;

    t1 = clock();

    while (diff > sigma) {
        k++;
        multiplication_pagerank(x, y);      
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

    t2 = clock();
    double temps_total = ((double)(t2 - t1)) / CLOCKS_PER_SEC;
    printf("Convergence en %d itérations (Diff = %e), temps total = %.4f sec\n", k, diff, temps_total);
}

void liberer_memoire() {
    
    free(x);
    free(y);
    free(p);
    free(est_dangling); 
    free(degre_sortant); 
    free(sncd.groupes);

    x = NULL;
    y = NULL;
    p = NULL;
    est_dangling = NULL;
    degre_sortant = NULL;
    sncd.groupes = NULL;
}

// [Les autres fonctions existantes (initialiser_vecteurs_pagerank, identifier_noeuds_dangling, 
//  iteration_puissance, norme_L1, recopier, multiplication_pagerank, liberer_memoire) 
//  restent inchangées par rapport à votre code original]