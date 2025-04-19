#include <stdio.h>
#include <pthread.h>

#define THREADS_NUMBER 10

void * hello_world(void * ptr);

int main(){
	int i, ids[THREADS_NUMBER];
	pthread_attr_t attr;
	pthread_t threads[THREADS_NUMBER];
	
	pthread_attr_init(&attr);

	/* Creacion de los hilos */
	for (i = 0; i< THREADS_NUMBER; i++){
		ids[i] = i;
		pthread_create(&threads[i], &attr, hello_world, &ids[i]);
	}
	
	/* Espera a que los hilos terminen */
	for (i = 0; i < THREADS_NUMBER; i++){
		pthread_join(threads[i], NULL);
	}
	
	return 0;
}

void * hello_world (void * ptr){
	int *p, id;
	p = (int *) ptr;
	id = *p;
	
	printf("\n Hola mundo soy el hilo %d\n", id);
	
	pthread_exit(0);
}
