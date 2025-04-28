#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

int max = -32;
int min = 1000;
int suma_paralela = 0;
int *vector;
int n;

void buscar_max();
void sumar();

int main(int argc, char *argv[]){
    n = atoi(argv[1]);
    vector = (int *) malloc(sizeof(int) * n);
    srand(time(NULL));
    omp_set_num_threads(5);

    for(int i = 0; i < n; i++){
        vector[i] = rand() % 100;
    }

    #pragma parallel
    {
        buscar_max();
        sumar();
    }

    int aux = 0;
    for(int i = 0; i < n; i++){
        printf("%d ",vector[i]);
        aux = aux + vector[i];
    }
    printf("\n");
    printf("la suma de los elementos es: %d\n",aux);

    printf("max quedo en : %d\n", max);
    printf("Suma paralela quedo en : %d\n",suma_paralela);

    return 0;
}

void buscar_max(){
    for(int i = 0; i< n; i++){
        if(vector[i] > max){
            max = vector[i];
        }
    }
}

void sumar(){
    for(int i = 0; i < n; i++){
        suma_paralela += vector[i];
    }
}
