
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

struct matrice {
    int N;
    float **V;
};
void liberer(struct matrice t) {
    for (int i = 0; i < t.N; i++) {
        free(t.V[i]);
    }
    free(t.V);
}

void sizer(){
    printf("float: %lu\n", sizeof(float));
    printf("double: %lu\n", sizeof(double));
    printf("long double: %lu\n", sizeof(long double));
    printf("int: %lu\n", sizeof(int));
    printf("long int: %lu\n", sizeof(long int));
    printf("char: %lu\n", sizeof(char));
    printf("short: %lu\n", sizeof(short));
    printf("long: %lu\n", sizeof(long));
    printf("long long: %lu\n", sizeof(long long));
    printf("unsigned int: %lu\n", sizeof(unsigned int));
    printf("unsigned char: %lu\n", sizeof(unsigned char));
    printf("unsigned short: %lu\n", sizeof(unsigned short));
    printf("unsigned long: %lu\n", sizeof(unsigned long));
    printf("unsigned long long: %lu\n", sizeof(unsigned long long));
}

// Fonction pour écrire une matrice en format creux de la maniere suivante: dabord premiere ligne ya la taille de la matrice, puis pour chaque ligne on a le nombre de valeurs non nulles dans cette ligne et les indices de ces valeurs
// les float sont stocker en 32 bits ya des bits de mentisse et des bits d'exposant il yen a 23 pour la mantisse et 8 pour l'exposant
// les exposants sont en complement a 2 donc au max 8 bits la plus petite valeur cest 
// en terme de temps de calcule il faudrait mettre des float que des double vrai ou faux ? : vrai
// en terme de memoire il faudrait mettre des float que des double vrai ou faux ? : faux


//ecrire la matrice du td en format creux - lire matrice - ecrire matrice - stocker en creux - maintenant faut stocker en memoire de la meuilleur maniere on veut faire xP->y (x vecteur p matrice)  - 


// notre matrice doit etre moins de 4 milliard de valeurs 2^32=4(2^9)^3 ~ 4 milliard


int main() {
    //struct matrice t;
    //t=lire_matrice("matricepleine.txt");
    //afficher_matrice(t);
    //ecrire_matrice("matriceecriture.txt",t);

    //p=malloc(M*sizeof(struct elem)); 
    //if (p==NULL) exit(21);
    //x=malloc(C*sizeof(proba)); 
    //if (x==NULL) exit(22);
    //y=malloc(L*sizeof(proba));
    //if (y==NULL) exit(23);

    //initx(x);
    //tantque(delta>epsilon)
    //mettre y a zero
    //mult(x,p,y)
    //recopie(x,y)


    return 0;
    
}