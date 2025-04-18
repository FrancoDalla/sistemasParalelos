/*
Codigo que compute la suma de los vectores que sea Ai = Bi + Ci
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>


/*Variables compartidas */
double *a, *b, *c;

typedef struct {
    int inicio;
    int fin;
} Rango;

double dwalltime(){
	double seconds;
	struct timeval time_value;

	gettimeofday(&time_value, NULL);
	seconds = time_value.tv_sec + time_value.tv_usec/1000000.0;
	return seconds;
}

void * calculo(void * ptr){
    Rango *rango = (Rango*)ptr;
    int inicio = rango->inicio, fin = rango->fin;
    for(int i = inicio; i < fin;i++){
        a[i] = b[i] + c[i];
    }
}

int main(int argc, char *argv[]){
	int i, j, k, n,t;
	pthread_attr_t attr;
	double timetick;
	pthread_attr_init(&attr);

	/*Agregar despues validacion de argumentos*/
	n = atoi(argv[1]);
	t = atoi(argv[2]);

	pthread_t hilos[t];

	/*alocacion de matrices */
	a = (double*) malloc(sizeof(double) * n);
	b = (double*) malloc(sizeof(double) * n);
	c = (double*) malloc(sizeof(double) * n);

	/*Inicializacion de vectores */
	for(i = 0; i<n; i++){
		   b[i] = 1;
		   c[i] = 1;
		}

	Rango rangos[t];

	int tamaño_base = n / t;
	int extras = n % t;
	int inicio = 0;

	/*Inicializacion de los struct */
	for(i = 0; i < t; i++){
        rangos[i].inicio = inicio;
        rangos[i].fin = inicio + tamaño_base -1 + (i < extras ? 1: 0);
        inicio = rangos[i].fin + 1;
	}

	timetick = dwalltime();

	/*Creacion de hilos */
	for(i = 0; i< t; i++){
	   pthread_create(&hilos[i], &attr, calculo, &rangos[i]);
	}

	/*Espera a finalizacion de los hilos */
	for(i = 0; i < t; i++){
	   pthread_join(hilos[i], NULL);
	}

	double time = dwalltime() - timetick;

	printf("Tiempo en calcular: %f\n",
	time);

	printf("Impresion del vector A: \n");
	for(i = 0; i<n ;i++){
	   printf("%f ",
				a[i]);
	}
	printf("\n");


	return 0;
}
