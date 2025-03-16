#include <stdio.h>
#include <stdlib.h>

typedef float proba;
typedef int indice;

struct elem{
    indice i,j;
    proba val;
};

struct elem *p;
proba *x;
proba *y;

indice L=0,C=0,M=0;//L cest le nbr de lignes C cest le nbr de colonnes M cest le nbr de valeur non nulles

struct matrice {
    indice N;
    proba **V;
};
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


void initx(proba *x, indice taille) {
    for (indice i = 0; i < taille; i++) {
        x[i] = 1.0 / taille;
    }
}


void mult (proba *x,struct elem *p,proba *y){
    indice i,j,k;
    struct elem e;
    for (k=0;k<M;k++){
        e=p[k];
        i=e.i;
        j=e.j;
        y[j]+=x[i]*e.val;
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






int main() {
    struct matrice t;
    
    // 1. Lire la matrice pleine depuis un fichier
    printf("Lecture de la matrice pleine...\n");
    t = lire_matrice("matricepleine.txt");

    // Vérification de la taille de la matrice
    if (t.N <= 0) {
        printf("Erreur : La taille de la matrice est invalide (%d)\n", t.N);
        return 1;
    }

    // Initialisation de L et C (nombre de lignes et colonnes)
    L = C = t.N;

    // 2. Afficher la matrice lue
    printf("Matrice lue :\n");
    afficher_matrice(t);

    // 3. Sauvegarder la matrice lue dans un autre fichier
    printf("Écriture de la matrice pleine dans 'matriceecriture.txt'...\n");
    ecrire_matrice("matriceecriture.txt", t);

    // 4. Conversion en matrice creuse
    printf("Conversion en matrice creuse...\n");
    ecrire_matrice_creuse("matricecreuse.txt", t);

    // 5. Charger la matrice creuse et la reconvertir en pleine
    printf("Lecture de la matrice creuse et conversion en pleine...\n");
    struct matrice t_creuse = matrice_creuse_a_pleine("matricecreuse.txt");

    printf("Matrice convertie depuis le format creux :\n");
    afficher_matrice(t_creuse);

    // 6. Calculer M (nombre de valeurs non nulles) à partir du fichier creux
    printf("Calcul du nombre d'éléments non nuls dans la matrice creuse...\n");
    M = calculer_nombre_elements_non_nuls_creuse("matricecreuse.txt");
    printf("M (nombre d'éléments non nuls) = %d\n", M);

    if (M <= 0) {
        printf("Erreur : M est invalide (%d), impossible de continuer.\n", M);
        return 1;
    }

    // 7. Allocation dynamique du vecteur x, y et tableau p
    x = malloc(C * sizeof(proba));
    y = malloc(C * sizeof(proba));
    p = malloc(M * sizeof(struct elem));

    if (x == NULL || y == NULL || p == NULL) {
        printf("Erreur d'allocation mémoire pour x, y ou p\n");
        return 1;
    }

    // 8. Initialisation du vecteur x avec 1/N
    initx(x, t.N);

    // 9. Remplir `p` avec les éléments non nuls de la matrice
    remplir_triplets(p, t);

    // Vérification du contenu de `p`
    printf("=== Contenu de p (matrice creuse sous forme de triplets) ===\n");
    for (indice k = 0; k < M; k++) {
        printf("p[%d] : i=%d, j=%d, val=%.6f\n", k, p[k].i, p[k].j, p[k].val);
    }
    printf("============================================================\n");

    // 10. Mise à zéro de y
    mettreazero(y);

    // 11. Appliquer la multiplication
    printf("Multiplication matrice x vecteur...\n");
    mult(x, p, y);

    // 12. Afficher le résultat de la multiplication
    printf("Résultat de la multiplication :\n");
    for (int i = 0; i < C; i++) {
        printf("%f ", y[i]);
    }
    printf("\n");

    // 13. Libération de la mémoire
    free(x);
    free(y);
    free(p);

    for (int i = 0; i < t.N; i++) {
        free(t.V[i]);
        free(t_creuse.V[i]);
    }
    free(t.V);
    free(t_creuse.V);

    return 0;
}
