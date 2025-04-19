/*
Desarrolle un algoritmo paralelo que compute la multiplicación de matrices cuadradas de NxN. Primero,
considere a la versión optimizada del ejercicio 6 de la práctica anterior como algoritmo base. Luego,
paralelice la versión que computa por bloques. Mida el tiempo de ejecución para N={512, 1024, 2048, 4096}
y T={2,4,8}. Analice el rendimiento.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

double dwalltime();
void initmatrix(double *matrix, int n, int transpose, int value);
void * multmatrix(void *ptr);
void blkmul(double *ablk, double *bblk, double *cblk);

double *a;
double *b;
double *c;
int n,t, bs;

typedef struct {
    int inicio;
    int fin;
} Rango;

int main(int argc, char *argv[]){
    int i;
    double timetick;

    /*conversion de parametros a integer y asignacion a variables */
    n = atoi(argv[1]);
    t = atoi(argv[2]);
    bs = atoi(argv[3]);

    pthread_attr_t attr;
    pthread_t hilos[t];
    pthread_attr_init(&attr);

    Rango rangos[t];

    int bloques_por_hilo = bs / t;
    int bloques = n*n / bs;

    /*Reserva de memoria para las matrices */
    a = (double *) malloc(n * n * sizeof(double));
    b = (double *) malloc(n * n * sizeof(double));
    c = (double *) malloc(n * n * sizeof(double));

    /*Inicializacion de matrices a, b y c */
    initmatrix(c,n,0,0);
    initmatrix(a, n, 0,1);
    initmatrix(b, n, 1, 1);

    timetick = dwalltime();

    int aux = 0;

    for(i = 0; i < t; i++){
        rangos[i].inicio = i * bloques_por_hilo;
        aux += bloques_por_hilo;
        rangos[i].fin = aux;
        pthread_create(&hilos[i], &attr,multmatrix, &rangos[i]);
    }

    for(i = 0; i < t; i++){
        pthread_join(hilos[i], NULL);
    }

    double time = dwalltime() - timetick;

    printf("Se tardo : %f segundos \n",
        time);

    printf("Impresion de la matriz a: \n");

    for(i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%f ",
                a[i + j * n]);
        }
        printf("\n");
    }

    printf("Impresion de la matriz b: \n");

    for(i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%f ",
                b[i + j * n]);
        }
        printf("\n");
    }

    printf("Impresion de la matriz c: \n");

    for(i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%f ",
                c[i * j + n]);
        }
        printf("\n");
    }

    free(a);
    free(b);
    free(c);

    return 0;
}






double dwalltime(){
	double seconds;
	struct timeval time_value;

	gettimeofday(&time_value, NULL);
	seconds = time_value.tv_sec + time_value.tv_usec/1000000.0;
	return seconds;
}

void initmatrix(double *matrix, int n, int transpose, int value){
    int i, j;

    if(transpose == 0){
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                matrix[i*n + j] = value;
            }
        }
    }
    else{
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                matrix[j * n + i] = value;
            }
        }
    }
}

void * multmatrix(void *ptr){
    Rango *rango = (Rango *)ptr;
    int first = rango ->inicio, last = rango ->fin;
    int i, j, k;
    for(i = first; i < last; i += bs){
        for(j = first; j < last; j += bs){
            for(k = 0; k < n; k +=bs){
                blkmul(&a[i * n + k], &b[j * n + k], &c[i * n + j]);
            }
        }
    }
}

void blkmul(double *ablk, double *bblk, double *cblk){
    int i, j, k;

    for (i = 0; i < bs; i++){
        for(j = 0; j < bs; j++){
            for(k = 0; k < bs; k++){
                cblk[i * n + j] += ablk[i * n + k] * bblk[j * n + k];
            }
        }
    }
}
