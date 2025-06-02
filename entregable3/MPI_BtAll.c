/*Este codigo realiza el calculo haciendo que la transposición de B sea realizada por cada proceso individualmente */
/*
    Ecuacion: R = ((maxA * maxB - MinA * MinB)/PromA * PromB) * [A * B] + C * Bt
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define DBL_MAX 2.2250738585072014e-308
#define DBL_MIN -2.2250738585072014e-308
#define BS 128

double dwalltime();
void transpose(double *m, double *mt);
void blkmul(double *ablk, double *bblk, double *cblk);
void alocar_memoria_maestro();
void alocar_memoria_trabajador();
void inicializar_matriz(double *matriz, double valor);

double *a;
double *b, *bt;
double *c;
double *ab;
double *cbt;
double max_a,min_a,promedio_a;
double max_b,min_b,promedio_b;
double escalar;
double *r;
int n;
int procesos;
int identificador;
int matriz_tamaño;
int carga_trabajo;


int main(int argc, char *argv[]){
    n = atoi(argv[1]);
    matriz_tamaño = n * n;

    max_a = DBL_MIN; max_b = DBL_MIN;
    min_a = DBL_MAX, min_b = DBL_MIN;
    promedio_a = 0; promedio_b = 0;

    b = (double *) malloc(sizeof(double) * matriz_tamaño);
    bt = (double *) malloc(sizeof(double) * matriz_tamaño);

    MPI_Init(&argc, &argv);

        MPI_Comm_size(MPI_COMM_WORLD, &procesos);
        MPI_Comm_rank(MPI_COMM_WORLD, &identificador);
        carga_trabajo = matriz_tamaño / procesos;


        if(identificador == MASTER){
            alocar_memoria_maestro();
        }
        else{
            alocar_memoria_trabajador()
        }

    MPI_Finalize();

    printf("Soy el proceso %d, adios\n",identificador);

    return 0;
}



double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

void transpose(double *m, double *m_trans) {
    int j, i;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            m_trans[i + j * n] = m[i * n + j];
        }
    }
}

void blkmul(double *ablk, double *bblk, double *cblk) {
    int i, j, k;
    for (i = 0; i < BS; i++) {
        for (j = 0; j < BS; j++){
            for(k = 0; k < BS; k++){
                cblk[i*n+j] += ablk[i*n+k] * bblk[k+j*n];
            }
        }
    }
}

void alocar_memoria_maestro(){
    a = (double *) malloc(sizeof(double) * matriz_tamaño);
    r = (double * ) malloc(sizeof(double) * matriz_tamaño);
    ab = (double *) malloc(sizeof(double) * matriz_tamaño);
    cbt = (double *) malloc(sizeof(double) * matriz_tamaño);
}

void alocar_memoria_trabajador(){
    a = (double *) malloc(sizeof(double) * carga_trabajo);
    r = (double * ) malloc(sizeof(double) * carga_trabajo);
    ab = (double *) malloc(sizeof(double) * carga_trabajo);
    cbt = (double *) malloc(sizeof(double) * carga_trabajo);
}

void inicializar_matriz(double *matriz, double valor){
    for(int i = 0; i < n; i++){
        for(int j=0;j<n;j++){
            matriz[i*n+j] = valor;
        }
    }
}
