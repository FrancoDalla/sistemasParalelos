#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>


int n, t, buscar, ocurrencias;
int *vector;
pthread_mutex_t search_value_lock;
typedef struct{
    int inicio;
    int fin;
    int id;
} Rango;

double dwalltime();
void * sumar(void * ptr);

void make_vector();

int main(int argc, char *argv[]){
    n = atoi(argv[1]);
    t = atoi(argv[2]);
    buscar = atoi(argv[3]);
    ocurrencias = 0;
    Rango rangos[t];
    vector = (int *) malloc(sizeof(int) * n);
    //make_vector();
    int cantidad_trabajo = n / t;

    pthread_attr_t attr;
    pthread_t threads[t];

    pthread_attr_init(&attr);

    pthread_mutex_init(&search_value_lock, NULL);


    int trabajo_actual = n;
    int aux = 0;
    int hilos_actual = t;

    double timetick = dwalltime();

    for(int i = 0; i < t; i++){

        printf("El hilo con id %d quedara con una carga de trabajo de: %d\n",i,(trabajo_actual / hilos_actual));
        rangos[i].id = i;
        rangos[i].inicio = aux;
        aux += trabajo_actual / hilos_actual;
        rangos[i].fin = aux;
        trabajo_actual = trabajo_actual - (trabajo_actual / hilos_actual);
        hilos_actual --;
        pthread_create(&threads[i], &attr, sumar, &rangos[i]);
    }

    for(int i = 0; i < t; i++){
        pthread_join(threads[i], NULL);
    }

    double work_time = dwalltime() - timetick;

    printf("Tardo %f\n",
        work_time);

    printf("Impresion del vector: \n");

    for(int i = 0; i < n; i++){
        printf("%d , ",
            vector[i]);
    }

    printf("\n");

    printf("Cantidad de ocurrencias: %d\n",
        ocurrencias);
    return 0;
}

double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

void * sumar(void * ptr){
    Rango *rango = (Rango *) ptr;
    int inicio = rango ->inicio;
    int fin = rango ->fin;

    for(int i = inicio; i<fin; i++){
        if(vector[i] == buscar){
            pthread_mutex_lock(&search_value_lock);
            ocurrencias += 1;
            pthread_mutex_unlock(&search_value_lock);
        }
    }
}

void make_vector(){
    srand(time(NULL));
    for(int i = 0; i < n; i++){
        vector[i] = rand() % 10;
    }
}
