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
    int N;
    float **V;
};


struct matrice init (int ma_taille){
    int i, j;
    struct matrice t;
    t.N = ma_taille;
    t.V = (float**)malloc(ma_taille * sizeof(float*));
    if (t.V == NULL) {
        printf("Erreur d'allocation\n");
        exit(1);
    }
    for (i = 0; i < ma_taille; i++){
        t.V[i] = (float*)malloc(ma_taille * sizeof(float));
        if (t.V[i] == NULL) {
            printf("Erreur d'allocation\n");
            exit(1);
        }
        // Initialiser la ligne à 0.0
        for (j = 0; j < ma_taille; j++){
            t.V[i][j] = 0.0;
        }
    }
    return t;
}


struct matrice lire_matrice(char *nom_fichier){
    FILE *f;
    struct matrice t;
    int i,j;
    f=fopen(nom_fichier,"r");
    if (f==NULL) {
        printf("Erreur d'ouverture du fichier\n");
        exit(1);
    }
    fscanf(f,"%d",&t.N);
    t=init(t.N);
    for(i=0;i<t.N;i++){
        for(j=0;j<t.N;j++){
            fscanf(f,"%f",&t.V[i][j]);
        }
    }
    fclose(f);
    return t;
}
void ecrire_matrice(char *nom_fichier, struct matrice t){
    FILE *f;
    int i,j;
    f=fopen(nom_fichier,"w");
    if (f==NULL) {
        printf("Erreur d'ouverture du fichier\n");
        exit(1);
    }
    fprintf(f,"%d\n",t.N);
    for(i=0;i<t.N;i++){
        for(j=0;j<t.N;j++){
            fprintf(f,"%f ",t.V[i][j]);
        }
        fprintf(f,"\n");
    }
    fclose(f);
}

void afficher_matrice (struct matrice t){
    int i,j;
    printf("%d\n",t.N);
    
    for(i=0;i<t.N;i++){
        for(j=0;j<t.N;j++){
            printf(" %f",t.V[i][j]);
        }
        printf("\n");
    }
}

void matrice_pleine_a_creuse(char *nom_fichier, struct matrice t) {
    FILE *f;
    int i, j, M = 0;

    // Comptage du nombre de valeurs non nulles
    for (i = 0; i < t.N; i++) {
        for (j = 0; j < t.N; j++) {
            if (t.V[i][j] != 0.0) {
                M++;
            }
        }
    }

    // Ouverture du fichier pour écriture
    f = fopen(nom_fichier, "w");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier\n");
        exit(1);
    }

    // Écriture des dimensions et du nombre de valeurs non nulles
    fprintf(f, "%d %d %d\n", t.N, t.N, M);

    // Écriture des éléments non nuls
    for (i = 0; i < t.N; i++) {
        for (j = 0; j < t.N; j++) {
            if (t.V[i][j] != 0.0) {
                fprintf(f, "%d %d %.6f\n", i, j, t.V[i][j]); // On garde la précision
            }
        }
    }

    fclose(f);
}

struct matrice creuse_a_pleine(char *nom_fichier) {
    FILE *f;
    struct matrice t;
    int i, j, M;
    
    printf("Ouverture du fichier %s\n", nom_fichier);

    f = fopen(nom_fichier, "r");
    if (f == NULL) {
        printf("Erreur d'ouverture du fichier\n");
        exit(1);
    }


    fscanf(f, "%d %d %d", &t.N, &t.N, &M); 

    printf("%d %d %d\n", t.N, t.N, M);
    t = init(t.N);


    for (int k = 0; k < M; k++) {
        fscanf(f, "%d %d %f", &i, &j, &t.V[i][j]);
    }

    fclose(f);
    return t;
}



int main() {
    struct matrice t;
    //t=lire_matrice("matricepleine.txt");
    printf("Matrice pleine\n");
    t= creuse_a_pleine("matricecreuse.txt");
    //afficher_matrice(t);
    //matrice_pleine_a_creuse("matricecreuse.txt",t);
    //ecrire_matrice("matriceecriture.txt",t);


    return 0;
    
}
