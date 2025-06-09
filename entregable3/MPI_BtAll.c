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
#define MASTER 0

double dwalltime();
void transpose(double *m, double *mt);
void blkmul(double *ablk, double *bblk, double *cblk);
void alocar_memoria_maestro();
void alocar_memoria_trabajador();
void inicializar_matrices_maestro();
void inicializar_matrices_trabajador();

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
int strip_size;
int actual_a, actual_b;


int main(int argc, char *argv[]){
    int i, j;
    n = atoi(argv[1]);
    matriz_tamaño = n * n;

    max_a = DBL_MIN; max_b = DBL_MIN;
    min_a = DBL_MAX, min_b = DBL_MIN;
    promedio_a = 0; promedio_b = 0;

    MPI_Init(&argc, &argv);

        MPI_Comm_size(MPI_COMM_WORLD, &procesos);
        MPI_Comm_rank(MPI_COMM_WORLD, &identificador);
        strip_size = n / procesos;

        if(identificador == MASTER){
            alocar_memoria_maestro();
            inicializar_matrices_maestro();
        }
        else{
            alocar_memoria_trabajador();
            inicializar_matrices_trabajador();
        }
        MPI_Barrier(MPI_COMM_WORLD);


        MPI_Bcast(b, sizeof(double) * matriz_tamaño, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        MPI_Scatter(a, sizeof(double) * (n * strip_size) , MPI_DOUBLE, a, sizeof(double) * (n * strip_size), MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        MPI_Scatter(c, sizeof(double) * (n * strip_size), MPI_DOUBLE, c, sizeof(double) * (n * strip_size), MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

        /*Cada proceso transpone la matriz para el mismo */
        transpose(b, bt);

        /*Calculo de minimo, maximo y promedio para A */
        for(i = 0; i < strip_size; i++){
            for(j = 0; j < n + j++){
                actual_a = a[i*n+j];

                if(actual_a > max_a){
                    max_a = actual_a;
                }

                if(actual_a < min_a){
                    min_a = actual_a;
                }

                promedio_a += actual_a;

            }
        }

        /*Calculo de minimo, maximo y promedio para A */
        int inicio = strip_size * identificador;
        int fin = inicio + strip_size;
        for(i = inicio; i < fin; i++){
            for(j = 0; j < n; j++){
                actual_b = b[i*n + j];
                if(actual_b > max_b){
                    max_b = actual_b;
                }

                if(actual_b < min_b){
                    min_b = actual_b;
                }

                promedio_b += actual_b;
            }
        }

        MPI_Reduce(&min_a, &min_a, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(&min_b, &min_b,1,MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(&promedio_a, &promedio_a, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(&promedio_b, &promedio_b, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

        if(identificador == MASTER){
            promedio_a = promedio_a / matriz_tamaño;
            promedio_b = promedio_b / matriz_tamaño;

            escalar = (max_a * max_b - min_a * min_b)/(promedio_a * promedio_b);
        }


    MPI_Finalize();

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
    for (i = 0; i < strip_size; i++) {
        for (j = 0; j < n; j++){
            for(k = 0; k < n; k++){
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
    bt = (double *) malloc(sizeof(double) * matriz_tamaño);
    b = (double *) malloc(sizeof(double) * matriz_tamaño);
}

void alocar_memoria_trabajador(){
    a = (double *) malloc(sizeof(double) * carga_trabajo);
    r = (double * ) malloc(sizeof(double) * carga_trabajo);
    ab = (double *) malloc(sizeof(double) * carga_trabajo);
    cbt = (double *) malloc(sizeof(double) * carga_trabajo);
    bt = (double *) malloc(sizeof(double) * matriz_tamaño);
    b = (double *) malloc(sizeof(double) * matriz_tamaño);
}

void inicializar_matrices_maestro(){
   	int i,j;
    for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i * n + j] = j;
			c[i * n + j] = j;
			b[i + j * n] = j;
			ab[i * n + j] = 0.0;
			cbt[i * n + j] = 0.0;
			r[i * n + j] = 0.0;
		}
	}
}

void inicializar_matrices_trabajador(){
    int i,j;
    for(i = 0; i < strip_size; i++){
        for(j = 0; j < strip_size; j++){
            ab[i*n+j] = 0.0;
            cbt[i*n+j] = 0.0;
            r[i*n+j] = 0.0;
        }
    }
}


/*
 * BARRIER
 * Si soy master, calculo tiempo.
 * SCATTER A
 * SCATTER C
 * BCAST B
 * GATHER
 * Si soy master, calculo tiempo.
 *
 * Cal
 */
