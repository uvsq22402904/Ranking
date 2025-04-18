#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> // Pour fabs

// Utiliser double pour une meilleure précision
typedef double proba;
typedef int indice;

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

// --- Déclarations de Fonctions ---
void lire_et_normaliser_mtx(char *nom_fichier);
void initialiser_vecteurs_pagerank();
void identifier_noeuds_dangling();
void iteration_puissance();
proba norme_L1(proba *v1, proba *v2, indice taille);
void recopier(proba *dest, proba *src, indice taille);
void multiplication_pagerank(proba *x_in, proba *y_out);
void liberer_memoire();

// --- Fonction Principale ---
int main() {
    // Nom du fichier en dur dans le code
    char *nom_fichier = "webbase-1M copy.mtx"; // Ou utilisez votre petit fichier d'exemple

    // 1. Lire la matrice, déterminer les dimensions, allouer la structure creuse et normaliser
    lire_et_normaliser_mtx(nom_fichier);
    if (N == 0 || M == 0 || p == NULL) {
        fprintf(stderr, "Erreur lors de la lecture et normalisation de la matrice.\n");
        return 1;
    }
    printf("Matrice : %d x %d, Éléments non nuls : %d\n", N, N, M);

    // 2. Allouer les vecteurs PageRank
    initialiser_vecteurs_pagerank();
    if (x == NULL || y == NULL) {
        fprintf(stderr, "Erreur d'allocation mémoire pour les vecteurs PageRank.\n");
        liberer_memoire();
        return 1;
    }

    // 3. Identifier les nœuds dangling (en utilisant degre_sortant calculé pendant la normalisation)
    identifier_noeuds_dangling();
    if (est_dangling == NULL) {
        fprintf(stderr, "Erreur d'allocation mémoire pour est_dangling.\n");
        liberer_memoire();
        return 1;
    }

    // 4. Effectuer l'itération de puissance
    iteration_puissance();

    // 5. Afficher les résultats (optionnel - peut être très long pour grand N)
    printf("\nScores de PageRank (premiers 20 ou moins):\n");
    int nb_affichage = (N < 20) ? N : 20;
    for (indice i = 0; i < nb_affichage; i++) {
        printf("Nœud %d : %.8e\n", i, x[i]); // Utilise l'index 0-based cohérent avec le code
    }
    if (N > nb_affichage) {
        printf("...\n");
    }

    // 6. Vérifier la somme (devrait être proche de 1.0)
    proba somme = 0.0;
    for (indice i = 0; i < N; i++) {
        somme += x[i];
    }
    printf("Somme du PageRank = %.15f\n", somme);

    // 7. Nettoyer la mémoire allouée
    liberer_memoire();

    return 0;
}

// --- Implémentations des Fonctions ---

void lire_et_normaliser_mtx(char *nom_fichier) {
    FILE *f;
    char ligne[256];
    int M_lu; // M lu depuis l'en-tête

    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        fprintf(stderr, "Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    // Ignorer les lignes de commentaire
    while (fgets(ligne, sizeof(ligne), f)) {
        if (ligne[0] != '%') {
            break;
        }
    }

    // Lire les dimensions et le nombre d'éléments non nuls de la première ligne non-commentaire
    if (sscanf(ligne, "%d %*d %d", &N, &M_lu) != 2) {
         // Gérer les cas où la deuxième dimension peut être manquante ou différente
         rewind(f); // Revenir au début pour réessayer la lecture après les commentaires
         while (fgets(ligne, sizeof(ligne), f)) {
            if (ligne[0] != '%') {
                break;
            }
         }
         if (sscanf(ligne, "%d %d %d", &N, &N, &M_lu) != 3) {
            fprintf(stderr, "Erreur: Impossible de lire les dimensions et le nombre d'éléments depuis l'en-tête.\n");
            fclose(f);
            exit(1);
         }
    }
    M = M_lu; // Utiliser le nombre de non-zéros lu depuis le fichier

    if (N <= 0 || M <= 0) {
        fprintf(stderr, "Erreur: Dimensions invalides N=%d, M=%d.\n", N, M);
        fclose(f);
        exit(1);
    }

    // Allouer la structure creuse et le tableau des degrés sortants
    p = malloc(M * sizeof(struct elem));
    degre_sortant = calloc(N, sizeof(int)); // Initialiser les degrés sortants à 0
    if (!p || !degre_sortant) {
        fprintf(stderr, "Erreur d'allocation mémoire pour p ou degre_sortant.\n");
        fclose(f);
        free(p); free(degre_sortant); p = NULL; degre_sortant = NULL; N=0; M=0;
        return;
    }

    // Lire les triplets et calculer les degrés sortants
    indice m_courant = 0;
    int ligne_elem, col_elem;
    double val_elem; // Lire comme double
    while (m_courant < M && fscanf(f, "%d %d %lf", &ligne_elem, &col_elem, &val_elem) == 3) {
        // Ajuster de l'index 1-based MTX à l'index 0-based C
        p[m_courant].i = ligne_elem - 1;
        p[m_courant].j = col_elem - 1;
        p[m_courant].val = 1.0; // Stocker 1.0 temporairement, sera normalisé plus tard

        if (p[m_courant].i < 0 || p[m_courant].i >= N || p[m_courant].j < 0 || p[m_courant].j >= N) {
             fprintf(stderr, "Attention: Index (%d, %d) hors limites [0, %d) lu à l'élément %d.\n",
                     ligne_elem, col_elem, N, m_courant + 1);
             // Optionnellement ignorer cet élément ou quitter
        } else {
            degre_sortant[p[m_courant].i]++; // Incrémenter le degré sortant pour le nœud source 'i'
        }
        m_courant++;
    }
    fclose(f);

    if (m_courant != M) {
        fprintf(stderr, "Attention: Nombre d'éléments lus (%d) différent de M annoncé (%d).\n", m_courant, M);
        M = m_courant; // Ajuster M au compte réel lu
        // Considérer réallouer p si significativement différent, bien que réduire soit OK
        struct elem * p_realloc = realloc(p, M * sizeof(struct elem));
        if (M > 0 && p_realloc == NULL){
            fprintf(stderr, "Erreur realloc p.\n");
            free(p); free(degre_sortant); p = NULL; degre_sortant = NULL; N=0; M=0;
            return;
        }
        p = p_realloc;
    }

    // Normaliser les valeurs dans p
    for (indice k = 0; k < M; k++) {
        indice i = p[k].i;
        if (degre_sortant[i] > 0) {
            p[k].val = 1.0 / (proba)degre_sortant[i];
        } else {
            // Ce cas ne devrait pas se produire pour les éléments réellement dans p,
            // mais la programmation défensive ne fait pas de mal.
            p[k].val = 0.0;
        }
    }
    // Le tableau degre_sortant est conservé pour identifier les nœuds dangling plus tard
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
    est_dangling = malloc(N * sizeof(int)); // Utiliser int (0 ou 1)
    if (!est_dangling) {
        free(degre_sortant); degre_sortant = NULL; // Libérer degre_sortant si dangling échoue
        return;
    }
    for (indice i = 0; i < N; i++) {
        est_dangling[i] = (degre_sortant[i] == 0);
    }
    // Nous n'avons plus besoin de degre_sortant après cela
    free(degre_sortant);
    degre_sortant = NULL;
}

// Norme L1: somme(abs(v1[i] - v2[i]))
proba norme_L1(proba *v1, proba *v2, indice taille) {
    proba somme = 0.0;
    for (indice i = 0; i < taille; i++) {
        somme += fabs(v1[i] - v2[i]); // Utiliser fabs pour double
    }
    return somme;
}

void recopier(proba *dest, proba *src, indice taille) {
    for (indice i = 0; i < taille; i++) {
        dest[i] = src[i];
    }
}

// Effectue une étape: y_out = alpha * P * x_in + contrib_dangling + contrib_teleport
void multiplication_pagerank(proba *x_in, proba *y_out) {
    // Initialiser y_out à zéro
    for (indice i = 0; i < N; i++) {
        y_out[i] = 0.0;
    }

    // 1. Calculer les contributions de base des nœuds non-dangling (partie P*x)
    //    y_j = somme_{i liens vers j} (x_i / degre_sortant(i))
    for (indice k = 0; k < M; k++) {
        // p[k].val contient déjà 1.0 / degre_sortant(p[k].i)
        y_out[p[k].j] += x_in[p[k].i] * p[k].val;
    }

    // 2. Calculer la contribution des nœuds dangling
    proba somme_dangling = 0.0;
    for (indice i = 0; i < N; i++) {
        if (est_dangling[i]) {
            somme_dangling += x_in[i];
        }
    }

    // 3. Combiner alpha, dangling, et téléportation
    //    La formulation Google: y = alpha * (P*x + dist_dangling) + (1-alpha)/N
    //    Soit P*x la somme calculée à l'étape 1.
    //    Soit dist_dangling = (somme(x_dangling)/N) ajoutée à chaque élément.
    //    Donc, y_out[i] = alpha * (y_out[i] + somme_dangling / N) + (1 - alpha) / N

    proba dist_dangling_par_noeud = 0.0;
    if (N > 0) {
         dist_dangling_par_noeud = somme_dangling / (proba)N;
    }
    proba teleport_par_noeud = (1.0 - alpha) / (proba)N;

    for (indice i = 0; i < N; i++) {
        y_out[i] = alpha * (y_out[i] + dist_dangling_par_noeud) + teleport_par_noeud;
    }

    // Optionnel: Rescaler y_out pour qu'il somme à 1.0 (aide à contrer la dérive en virgule flottante)
    // proba somme_courante = 0.0;
    // for (indice i = 0; i < N; i++) somme_courante += y_out[i];
    // if (somme_courante != 0.0) { // Éviter division par zéro si quelque chose a très mal tourné
    //     proba echelle = 1.0 / somme_courante;
    //     for (indice i = 0; i < N; i++) y_out[i] *= echelle;
    // }
    // Note: Avec la formulation Google ci-dessus, la somme devrait théoriquement rester 1.
}

void iteration_puissance() {
    indice k = 0;
    proba diff = 1.0; // Différence initiale > sigma
    struct timeval t1, t2;

    gettimeofday(&t1, NULL);

    while (diff > sigma) {
        k++;
        multiplication_pagerank(x, y);      // Calculer y basé sur x
        diff = norme_L1(x, y, N);         // Calculer la différence ||x - y||_1
        recopier(x, y, N);                // Copier y vers x pour la prochaine itération

        if (k % 10 == 0) { // Afficher la progression occasionnellement
             printf("Itération %d, Diff = %e\n", k, diff);
        }
         if (k > 10000) { // Ajouter une sécurité pour les exécutions très longues
             printf("Attention: Dépassement de 10000 itérations.\n");
             break;
         }
         // Vérifier NaN ou Inf
         if (isnan(diff) || isinf(diff)) {
             fprintf(stderr, "Erreur: Convergence échouée (NaN ou Inf détecté) à l'itération %d.\n", k);
             break;
         }
    }

    gettimeofday(&t2, NULL);
    double temps_total = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
    printf("Convergence en %d itérations (Diff = %e), temps total = %.4f sec\n", k, diff, temps_total);
}

void liberer_memoire() {
    free(x);
    free(y);
    free(p);
    free(est_dangling); // Note: degre_sortant est libéré plus tôt

    // Assigner NULL séparément
    x = NULL;
    y = NULL;
    p = NULL;
    est_dangling = NULL;
}