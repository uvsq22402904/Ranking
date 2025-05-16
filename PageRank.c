/*
 * - printf, fprintf : affichage formaté dans la console ou un fichier
 * - fopen, fclose : manipulation de fichiers
 * - fgets : lecture de lignes depuis un fichier
*/
#include <stdio.h>
/*
 * - malloc, calloc, realloc, free : allocation et libération dynamique de mémoire
 * - rand, srand : génération de nombres aléatoires (pour la suppression d'arcs)
*/
#include <stdlib.h>
/*
 * - struct timeval : structure pour stocker un horodatage précis
 * - gettimeofday : fonction pour obtenir le temps actuel avec précision microseconde
*/
#include <sys/time.h>
/*
 * - fabs : calcul de la valeur absolue pour nombres flottants
 * - isnan, isinf : détection des valeurs spéciales (NaN, Inf)
*/
#include <math.h>

// Utiliser double pour une meilleure précision numérique
typedef double proba;
// Utiliser int pour les indices, suffisant pour indexer les nœuds du graphe
typedef int indice;

/* Nombre de nœuds dans le graphe (dimension de la matrice d'adjacence) */
indice N = 0; 

/* Nombre d'éléments non nuls (arcs/liens) dans le graphe */
indice M = 0; 

/* Facteur d'amortissement (damping factor)
 * Représente la probabilité qu'un utilisateur suive un lien plutôt que de "téléporter"
 * α proche de 1 donne plus d'importance à la structure du graphe
 * α plus petit augmente l'importance du facteur de téléportation */
proba alpha[] = {0.85, 0.90, 0.95, 0.99};  

/* Pourcentages de suppression aléatoire des liens
 * Permet d'évaluer la robustesse de l'algorithme face à des modifications de la topologie */
proba pourcentage_suppression[] = {0.0, 0.1, 0.2}; 

/* Critère de convergence (epsilon)
 * L'algorithme s'arrête lorsque ||x_{k+1} - x_k||_1 < sigma
 * Valeur stricte pour garantir une précision élevée avec le type double */
proba sigma = 1e-9; 

/* Tableau booléen indiquant les noeuds sans liens sortants (dangling nodes)
 * Ces nœuds nécessitent un traitement spécial dans l'algorithme PageRank */
indice *est_dangling = NULL; 

// Structure pour les éléments de matrice creuse stocke les triplets(i,j,val)
struct elem {
    indice i, j; // Lien de i vers j
    proba val;  // Probabilité (sera normalisée à 1/degre_sortant(i))
};

/* Tableau des éléments non nuls (arcs du graphe) */
struct elem *p = NULL; 

/* Vecteur PageRank courant x_k */
proba *x = NULL;       

/* Vecteur PageRank suivant x_{k+1} */
proba *y = NULL;       

/* Stocke le degré sortant de chaque nœud
 * Utilisé pour normaliser les probabilités de transition */
indice *degre_sortant = NULL; 







//ecrire la matrice du td en format creux - lire matrice - ecrire matrice - stocker en creux - maintenant faut stocker en memoire de la meuilleur maniere on veut faire xP->y (x vecteur p matrice)  - 


// notre matrice doit etre moins de 4 milliard de valeurs 2^32=4(2^9)^3 ~ 4 milliard









// --- Déclarations de Fonctions ---


/**
 * Lit une matrice au format MTX et normalise ses colonnes
 * Pour transformer la matrice d'adjacence en matrice stochastique
 * @param nom_fichier Chemin vers le fichier MTX
 */
void lire_et_normaliser_mtx(char *nom_fichier);

/**
 * Lit une matrice MTX et supprime aléatoirement un pourcentage d'arcs
 * Utilisé pour analyser la robustesse de PageRank face aux perturbations du graphe
 * @param nom_fichier Chemin vers le fichier MTX
 * @param pourcentage_suppression Pourcentage d'arcs à supprimer aléatoirement [0,1]
 */
void lire_et_normaliser_mtx_avec_suppression(char *nom_fichier, proba pourcentage_suppression);

/**
 * Initialise les vecteurs PageRank avec une distribution uniforme
 * x = [1/N, 1/N, ..., 1/N]
 */
void initialiser_vecteurs_pagerank();

/**
 * Identifie les nœuds sans liens sortants (dangling nodes)
 * Ces nœuds créent des colonnes nulles dans la matrice stochastique
 * et nécessitent un traitement spécial dans l'algorithme
 */
void identifier_noeuds_dangling();

/**
 * Execute l'algorithme de la puissance itérée pour calculer le vecteur PageRank
 * Itère jusqu'à convergence: ||x_{k+1} - x_k||_1 < sigma
 * @param alpha Facteur d'amortissement
 * @param iterations Pointeur pour stocker le nombre total d'itérations
 * @param temps Pointeur pour stocker le temps de calcul
 */
void iteration_puissance(proba alpha, indice *iterations, proba *temps);

/**
 * Calcule la norme L1 entre deux vecteurs: ||v1 - v2||_1 = Σ|v1[i] - v2[i]|
 * Utilisée comme critère de convergence de l'algorithme
 * @param v1 Premier vecteur
 * @param v2 Second vecteur
 * @param taille Dimension des vecteurs
 * @return La norme L1 de la différence
 */
proba norme_L1(proba *v1, proba *v2, indice taille);

/**
 * Copie un vecteur source dans un vecteur destination
 * @param dest Vecteur destination
 * @param src Vecteur source
 * @param taille Dimension des vecteurs
 */
void recopier(proba *dest, proba *src, indice taille);

/**
 * Effectue une multiplication PageRank: w = α(Pz + d(z)) + (1-α)v
 * Où P est la matrice stochastique normalisée
 * d(z) est la contribution des dangling nodes
 * v est le vecteur de téléportation uniforme
 * @param z Vecteur d'entrée
 * @param w Vecteur résultat
 * @param alpha Facteur d'amortissement
 */
void multiplication_pagerank(proba *z, proba *w, proba alpha);

/**
 * Enregistre les résultats de l'expérience dans un fichier CSV
 * Pour analyse ultérieure et génération de graphiques
 * @param nom_fichier Chemin du fichier de sortie
 * @param nom_matrice Nom de la matrice traitée
 * @param pourcentage_suppression Pourcentage de suppression appliqué
 * @param alpha Facteur d'amortissement utilisé
 * @param iterations Nombre d'itérations jusqu'à convergence
 * @param temps Temps de calcul en secondes
 * @param somme_pagerank Somme des scores PageRank (devrait être proche de 1)
 */
void enregistrer_resultats(char *nom_fichier, char *nom_matrice, proba pourcentage_suppression, proba alpha, indice iterations, proba temps, proba somme_pagerank);

/**
 * Libère la mémoire allouée dynamiquement
 * Prévient les fuites mémoire lors de l'exécution des multiples expériences
 */
void liberer_memoire();

// --- Fonction Principale ---
int main() {
    char *nom_fichier = "Matrix\\wikipedia-20051105.mtx";
    char *fichier_resultats = "wikipedia-20051105_resultats.csv";
    indice nb_suppressions = sizeof(pourcentage_suppression) / sizeof(pourcentage_suppression[0]);
    indice nb_alphas = sizeof(alpha) / sizeof(alpha[0]);

    /* Création de l'en-tête du fichier CSV pour les résultats */
    FILE *f = fopen(fichier_resultats, "w");
    if (f) {
        fprintf(f, "Matrice,Pourcentage_Suppression,Alpha,Iterations,Temps,Somme_PageRank\n");
        fclose(f);
    }

    /* Double boucle pour tester toutes les combinaisons de paramètres */
    for (indice i = 0; i < nb_suppressions; i++) {
        for (indice j = 0; j < nb_alphas; j++) {
            printf("\n==> Suppression: %.2f, Alpha: %.2f\n", pourcentage_suppression[i], alpha[j]);

            /* Lecture du fichier MTX avec suppression aléatoire de certains arcs */
            lire_et_normaliser_mtx_avec_suppression(nom_fichier, pourcentage_suppression[i]);
            if (N == 0 || M == 0 || p == NULL) {
                fprintf(stderr, "Erreur lors de la lecture et normalisation de la matrice.\n");
                continue;
            }

            printf("Matrice : %d x %d, Éléments non nuls : %d\n", N, N, M);

            /* Initialisation du vecteur PageRank avec distribution uniforme */
            initialiser_vecteurs_pagerank();
            if (x == NULL || y == NULL) {
                fprintf(stderr, "Erreur d'allocation mémoire pour les vecteurs PageRank.\n");
                liberer_memoire();
                continue;
            }

            /* Identification des nœuds dangling (sans liens sortants) */
            identifier_noeuds_dangling();
            if (est_dangling == NULL) {
                fprintf(stderr, "Erreur d'allocation mémoire pour est_dangling.\n");
                liberer_memoire();
                continue;
            }

            /* Exécution de l'algorithme de la puissance itérée */
            indice iterations;
            proba temps;
            iteration_puissance(alpha[j], &iterations, &temps);

            /* Affichage des scores PageRank (top 20) */
            printf("\nScores de PageRank (premiers 20 ou moins):\n");
            indice nb_affichage = (N < 20) ? N : 20;
            for (indice k = 0; k < nb_affichage; k++) {
                printf("Noeud %d : %.8e\n", k, x[k]);
            }
            if (N > nb_affichage) {
                printf("...\n");
            }

            /* Vérification que la somme des scores est bien égale à 1 */
            proba somme = 0.0;
            for (indice k = 0; k < N; k++) {
                somme += x[k];
            }
            printf("Somme du PageRank = %.15f\n", somme);

            /* Enregistrement des résultats pour analyse comparative */
            enregistrer_resultats(fichier_resultats, nom_fichier, pourcentage_suppression[i], alpha[j], iterations, temps, somme);

            /* Libération mémoire avant de passer à l'expérience suivante */
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
    
    /* Initialisation du générateur de nombres aléatoires */
    srand(time(NULL));
    
    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        fprintf(stderr, "Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    /* Passer les commentaires */
    while (fgets(ligne, sizeof(ligne), f)) {
        if (ligne[0] != '%') {
            break;
        }
    }

    /* Lire les dimensions */
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
    
    /* Allouer de la mémoire temporaire pour tous les éléments potentiels 
     * avant de connaître le nombre d'arcs qui seront conservés */
    struct elem *temp_p = malloc(L * sizeof(struct elem));
    degre_sortant = calloc(N, sizeof(indice));
    if (!temp_p || !degre_sortant) {
        fprintf(stderr, "Erreur d'allocation mémoire.\n");
        fclose(f);
        free(temp_p); free(degre_sortant);
        temp_p = NULL; degre_sortant = NULL; N=0; M=0;
        return;
    }

    /* Lecture des éléments avec suppression aléatoire */
    indice elements_lus = 0;
    indice ligne_elem, col_elem;
    proba val_elem;
    char buffer[256];

    while (elements_lus < L && fgets(buffer, sizeof(buffer), f)) {
        /* Gestion des deux formats MTX possibles */
        int items_read = sscanf(buffer, "%d %d %lf", &ligne_elem, &col_elem, &val_elem);
        if (items_read == 2) {
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
            /* Décision aléatoire: conserver l'arc avec probabilité (1-pourcentage_suppression) */
            proba u = (proba)rand() / RAND_MAX;
            
            if (u >= pourcentage_suppression) {
                /* L'arc est conservé */
                temp_p[elements_lus].i = i;
                temp_p[elements_lus].j = j;
                temp_p[elements_lus].val = 1.0;
                degre_sortant[i]++;
                elements_lus++;
            }
            /* Sinon, l'arc est supprimé (simulant la disparition d'un lien web) */
        }
    }
    fclose(f);
    
    /* Mise à jour du nombre réel d'éléments après suppression */
    M = elements_lus;
    printf("Arcs après suppression aléatoire de %.1f%% : %d (sur %d initialement)\n", 
           pourcentage_suppression * 100, M, L);
    
    /* Allocation définitive pour le tableau p avec la taille exacte requise */
    p = malloc(M * sizeof(struct elem));
    if (!p) {
        fprintf(stderr, "Erreur d'allocation mémoire pour p.\n");
        free(temp_p); free(degre_sortant);
        p = NULL; degre_sortant = NULL; N=0; M=0;
        return;
    }
    
    /* Copie des éléments conservés de temp_p vers p */
    for (indice k = 0; k < M; k++) {
        p[k] = temp_p[k];
    }
    free(temp_p);
    
    /* Normalisation par degré sortant pour obtenir la matrice stochastique */
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
    // Calcul de la probabilité initiale uniforme pour chaque nœud (1/N)
    proba rang_initial = 1.0 / (proba)N;
    for (indice i = 0; i < N; i++) {
        x[i] = rang_initial;
    }
    // Note : y n'est pas initialisé ici car il sera rempli dans multiplication_pagerank
}

void identifier_noeuds_dangling() {
    // indique si chaque nœud est dangling (1 si oui, 0 sinon)
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
    // w[j] += z[i] * P[i,j]
    for (indice k = 0; k < M; k++) {

        w[p[k].j] += z[p[k].i] * p[k].val;
    }

    proba somme_dangling = 0.0;
    for (indice i = 0; i < N; i++) {
        if (est_dangling[i]) {
            somme_dangling += z[i];
        }
    }

    // Calcul de la contribution dangling par nœud (somme_dangling / N)
    proba dist_dangling_par_noeud = 0.0;
    if (N > 0) {
         dist_dangling_par_noeud = somme_dangling / (proba)N;
    }
    //Modelise le saut aléatoire du surfeur vers n'importe quel nœud
    proba teleport_par_noeud = (1.0 - alpha) / (proba)N;
    // w[i] = alpha * (P*z + dangling) + (1-alpha)/N
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
        // Vérifie si le nombre d'itérations dépasse 10000  Évite une boucle infinie en cas de non-convergence
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