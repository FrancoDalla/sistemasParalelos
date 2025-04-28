#include <stdio.h>
#include <omp.h>

#define N 1000

int main(){
    float v[N], sum = 0, avg;
    int i;

    for(i = 0; i < N; i++){
        v[i] = 1;
    }

    #pragma omp parallel
    {
        #pragma omp for reduction(+:sum)
        for(i = 0; i < N; i++)
            sum += v[i];
        #pragma omp single
        {
            avg = sum / N;
        }
    }

    printf("%f %f ->Resultados\n",
        sum, avg);

}
