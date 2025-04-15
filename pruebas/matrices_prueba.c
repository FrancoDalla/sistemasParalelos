#include <stdio.h>
#include <stdlib.h>

#define N 4

int main(){
	float * matriz = malloc(N*N * sizeof(float));
	
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			matriz[j * N + i] = j;
		}
	}

	
	for(int i = 0; i< N*N; i++){
		printf("%f\n",matriz[i]);
	}

	
	return 0;
}

