#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

typedef float proba;
typedef int indice;

#define abs(x) ((x) > 0 ? (x) : -(x))

indice L=0,C=0,M=0;//L cest le nbr de lignes C cest le nbr de colonnes M cest le nbr de valeur non nulles
proba alpha = 0.85;  // Facteur de téléportation
proba sigma = 1e-6;  // Critère de convergence (epsilon)

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

// Fonction de norme ||x - y||
float norme(proba *x,proba *y){
    indice i;
    float somme=0.0;
    for (i=0;i<C;i++){
        somme+=abs(x[i]-y[i]);
    }
    return somme;
}

struct matrice init (indice ma_taille){
    indice i,j;
    struct matrice t;
    t.N=ma_taille;
    t.V=(proba**)malloc(ma_taille*sizeof(proba*));
    if (t.V==NULL) {
        printf("Erreur d'allocation\n");
        exit(1);
    }
    for (i=0;i<ma_taille;i++){
        t.V[i]=(proba*)malloc(ma_taille*sizeof(proba));
        if (t.V[i]==NULL) {
            printf("Erreur d'allocation\n");
            exit(1);
        }
    }
    return t;

}

struct matrice lire_matrice(char *nom_fichier){
    FILE *f;
    struct matrice t;
    indice i,j;
    f=fopen(nom_fichier,"r");
    if (f==NULL) {
        printf("Erreur d'ouverture du fichier\n");
        exit(1);
    }
    fscanf(f,"%d",&t.N);
    t=init(t.N);
    for(i=0;i<t.N;i++){
        for(j=0;j<t.N;j++){
            fscanf(f, "%f", &t.V[i][j]);
        }
    }
    fclose(f);
    return t;
}
void ecrire_matrice(char *nom_fichier, struct matrice t){
    FILE *f;
    indice i,j;
    f=fopen(nom_fichier,"w");
    if (f==NULL) {
        printf("Erreur d'ouverture du fichier\n");
        exit(1);
    }
    fprintf(f,"%d\n",t.N);
    for(i=0;i<t.N;i++){
        for(j=0;j<t.N;j++){
            fprintf(f, "%.6f ", t.V[i][j]);
        }
        fprintf(f,"\n");
    }
    fclose(f);
}

void afficher_matrice (struct matrice t){
    indice i,j;
    printf("%d\n",t.N);
    
    for(i=0;i<t.N;i++){
        for(j=0;j<t.N;j++){
            printf(" %.6f", t.V[i][j]);
        }
        printf("\n");
    }
}



struct matrice matrice_creuse_a_pleine(char *nom_fichier) {
    FILE *f;
    struct matrice t;
    indice i, j, nbValeurs, col;
    proba valeur;

    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    // Lire la taille de la matrice
    fscanf(f, "%d", &t.N);
    t = init(t.N);

    // Lecture et remplissage de la matrice
    for (i = 0; i < t.N; i++) {
        fscanf(f, "%d", &nbValeurs);  // Nombre de valeurs non nulles dans la ligne

        for (j = 0; j < nbValeurs; j++) {
            fscanf(f, "%d %f", &col, &valeur);  // Lire colonne et valeur
            t.V[i][col - 1] = valeur;  // Stocker la vraie valeur dans la matrice
        }
    }

    fclose(f);
    return t;
}


struct matrice lire_matrice_creuse(char *nom_fichier) {
    FILE *f;
    struct matrice t;
    indice i, j, nbValeurs, col;

    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    // Lire la taille de la matrice
    fscanf(f, "%d", &t.N);
    t = init(t.N);

    // Lecture et remplissage de la matrice
    for (i = 0; i < t.N; i++) {
        fscanf(f, "%d", &nbValeurs);  // Nombre de 1 sur la ligne

        for (j = 0; j < nbValeurs; j++) {
            fscanf(f, "%d", &col);  // Lire les colonnes où il y a des 1
            t.V[i][col - 1] = 1.0;  // Ajuster pour l'index 0-based en C
        }
    }

    fclose(f);
    return t;
}


void ecrire_matrice_creuse(char *nom_fichier, struct matrice t) {
    FILE *f;
    indice i, j, count;

    f = fopen(nom_fichier, "w");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    // Écrire la taille de la matrice
    fprintf(f, "%d\n", t.N);

    // Écrire la matrice au format creux
    for (i = 0; i < t.N; i++) {
        count = 0;

        // Compter les 1 dans la ligne
        for (j = 0; j < t.N; j++) {
            if (t.V[i][j] != 0.0) {
                count++;
            }
        }

        fprintf(f, "%d", count);  // Nombre de valeurs dans la ligne

        // Écrire les indices des colonnes avec des 1
        for (j = 0; j < t.N; j++) {
            if (t.V[i][j] != 0.0) {
                fprintf(f, " %d %.6f", j + 1, t.V[i][j]);  // On remet en notation 1-based
            }
        }

        fprintf(f, "\n");
    }

    fclose(f);
}

// Initialiser le vecteur x avec 1/N
void initx(proba *x, indice taille) {
    for (indice i = 0; i < taille; i++) {
        x[i] = 1.0 / taille;
    }
}

// Multiplication matrice-vecteur
void mult (proba *x,struct elem *p,proba *y){
    indice i, j, k;
    struct elem e;

    // Réinitialisation du tableau y à zéro avant le calcul
    for (i = 0; i < C; i++) {
        y[i] = 0.0;
    }

    // Multiplication matrice-vecteur (mise à jour de PageRank)
    for (k = 0; k < M; k++) {
        e = p[k];
        i = e.i;  // Indice de départ du lien
        j = e.j;  // Indice d'arrivée du lien
        y[j] += alpha * x[i] * e.val;  // Application du facteur alpha
    }

    // Calcul du facteur de téléportation
    proba somme = 0.0;
    for (i = 0; i < C; i++) {
        somme += x[i];
    }

    proba teleporte = (1.0 - alpha) * somme / C;
    
    // Ajout du facteur de téléportation au vecteur y
    for (i = 0; i < C; i++) {
        y[i] += teleporte;
    }
}

void mettreazero(proba *y){
    indice i;
    for (i=0;i<C;i++){
        y[i]=0.0;
    }
}

void recopie(proba *x,proba *y){
    indice i;
    for (i=0;i<C;i++){
        x[i]=y[i];
    }
}

void remplir_triplets(struct elem *p, struct matrice t) {
    indice k = 0;
    for (indice i = 0; i < t.N; i++)
        for (indice j = 0; j < t.N; j++)
            if (t.V[i][j] != 0.0)
                p[k++] = (struct elem){i, j, t.V[i][j]};
}

indice calculer_nombre_elements_non_nuls_creuse(char *nom_fichier) {
    FILE *f = fopen(nom_fichier, "r");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier '%s'\n", nom_fichier);
        exit(1);
    }

    indice M = 0, nbValeurs, i, N;
    fscanf(f, "%d", &N);
    for (i = 0; i < N; i++) {
        fscanf(f, "%d", &nbValeurs);
        M += nbValeurs;
        while (nbValeurs--) {  // Lire et ignorer les valeurs
            int col;
            float val;
            fscanf(f, "%d %f", &col, &val);
        }
    }

    fclose(f);
    return M;
}

float une_iteration(proba *x, struct elem *p, proba *y) {
    mult(x, p, y);
    float s = norme(x, y);
    recopie(x, y);
    return s;
}

void iterer(proba *x, struct elem *p, proba *y) {
    indice k= 0;
    float s=1.0;
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    while (s > sigma) {
        k++;
        s=une_iteration(x, p, y);
        gettimeofday(&t2, NULL);
        //printf("Iteration %d : norme = %f, temps = %f\n", k, s, (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6);
    }
    gettimeofday(&t2, NULL);
    printf("Convergence en %d iterations, temps total = %f\n", k, (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6);
}



int main() {
    struct matrice t = matrice_creuse_a_pleine("matricecreuse.txt");

    if (t.N <= 0) {
        printf("Erreur : Taille de la matrice invalide.\n");
        return 1;
    }

    L = C = t.N;
    M = calculer_nombre_elements_non_nuls_creuse("matricecreuse.txt");

    x = malloc(C * sizeof(proba));
    y = malloc(C * sizeof(proba));
    p = malloc(M * sizeof(struct elem));

    if (!x || !y || !p) {
        printf("Erreur d'allocation mémoire\n");
        return 1;
    }

    initx(x, t.N);
    remplir_triplets(p, t); // ca bug
    mettreazero(y);

    printf("Début de l'itération...\n");
    iterer(x, p, y);

    free(x);
    free(y);
    free(p);

    for (indice i = 0; i < t.N; i++) {
        free(t.V[i]);
    }
    free(t.V);

    return 0;
}
