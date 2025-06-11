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
void recorrido_escalar();
void matmulblks(double *a, double *b, double *c);
void sumar_vector_escalar();

double *a;
double *b, *bt;
double *c;
double *ab;
double *cbt;
double suma_local_a = 0;
double suma_local_b = 0;
double promedio_a;
double promedio_b;
double escalar;
double *r;
int n;
int procesos;
int identificador;
int matriz_tamaño;
int strip_size;
int actual_a, actual_b;
double maximos[2] = { DBL_MIN, DBL_MIN };
double minimos[2] = { DBL_MAX, DBL_MAX };
double maximo[2];
double minimo[2];
MPI_Request requests[15];
MPI_Status status;

int main(int argc, char *argv[]){
    int i, j, verificacion;
    n = atoi(argv[1]);

    if(argc > 2){
        verificacion = atoi(argv[2]);
    }

    matriz_tamaño = n * n;

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


        MPI_Ibcast(b, matriz_tamaño, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[0]);
        MPI_Iscatter(a, n * strip_size, MPI_DOUBLE, a, n * strip_size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[1]);
        MPI_Iscatter(c, n * strip_size, MPI_DOUBLE, c, n * strip_size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[2]);


        /*Cada proceso transpone la matriz para el mismo */
        MPI_Wait(&requests[0], &status);
        transpose(b, bt);

        recorrido_escalar();


        MPI_Ireduce(maximos, maximo, 2, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD, &requests[3]);
        MPI_Ireduce(minimos, minimo, 2, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD, &requests[4]);
        MPI_Ireduce(&suma_local_a, &promedio_a, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD, &requests[5]);
        MPI_Ireduce(&suma_local_b, &promedio_b, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD,&requests[6]);


        MPI_Waitall(4, &requests[3], MPI_STATUSES_IGNORE);

        if(identificador == MASTER){
            promedio_a = promedio_a / matriz_tamaño;
            promedio_b = promedio_b / matriz_tamaño;
            escalar = (maximo[0] * maximo[1] - minimo[0] * minimo[1])/(promedio_a * promedio_b);
        }

        /*Calculo de A X B */
        matmulblks(a, b, ab);

        /*CALCULO DE C * Bt */
        matmulblks(c, bt, cbt);

        MPI_Ibcast(&escalar, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[7]);

        MPI_Wait(&requests[7], &status);

        /*Calculo de  escalar a (AXB + CXBt)*/
        sumar_vector_escalar();

        MPI_Igather(r, strip_size * n, MPI_DOUBLE, r, strip_size * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[8]);
        MPI_Wait(&requests[8], &status);

    MPI_Finalize();

    if(identificador == MASTER){
        printf("EL ESCALAR QUEDO: %f\n", escalar);
    }

    free(a);
    free(b);
    free(c);
    free(ab);
    free(cbt);
    free(bt);

    /* VERIFICACION */
    if(identificador == MASTER)
    {
        if(verificacion >= 1){
            if( escalar != 4.0){
                fprintf(stderr, "El valor escalar no dio el resultado esperado de 4, el valor fue: %f\n", escalar);
                exit(1);
            }

      		for (i = 0; i < n; i++) {
    			for (j = 0; j < n; j++) {
    				double ab = j * (n - 1.0) * n / 2.0;
    				double cbt =
    					(n - 1.0) * n * (2.0 * n - 1.0) / 6.0;
    				double expected = 4.0 * ab + cbt;

    				if (r[i * n + j] != expected) {
    					fprintf(stderr,
    						"Error en posición (%d,%d). Esperado valor: (%g), se obtuvo (%g).\n",
    						i, j, expected, r[i * n + j]);
    					exit(1);
    				}
    			}
    		}
        }
        printf("Los resultados en la matriz R son correctos\n");
    }

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

void matmulblks(double *a, double *b, double *c){
    int i, j, k, posicion_i, posicion_j;
    for(i = 0; i < strip_size; i+=BS){
        posicion_i = i * n;
        for(j = 0; j < n; j += BS){
            posicion_j = j * n;
            for(k = 0; k<n; k += BS){
                blkmul(&a[posicion_i + k], &b[posicion_j + k], &c[posicion_i + j]);
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
    c = (double *) malloc(sizeof(double) * matriz_tamaño);
}

void alocar_memoria_trabajador(){
    a = (double *) malloc(sizeof(double) * strip_size * n);
    r = (double * ) malloc(sizeof(double) * strip_size * n);
    ab = (double *) malloc(sizeof(double) * strip_size * n);
    cbt = (double *) malloc(sizeof(double) * strip_size * n);
    c = (double *) malloc(sizeof(double) * strip_size * n);

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
        for(j = 0; j < n; j++){
            ab[i*n+j] = 0.0;
            cbt[i*n+j] = 0.0;
            r[i*n+j] = 0.0;
        }
    }
}

void recorrido_escalar(){
    int i,j;

    /*Calculo de minimo, maximo y promedio para A */
    for(i = 0; i < strip_size; i++){
        for(j = 0; j < n; j++){
            actual_a = a[i*n+j];

            if(actual_a > maximos[0]){
                maximos[0] = actual_a;
            }

            if(actual_a < minimos[0]){
                minimos[0] = actual_a;
            }

            suma_local_a += actual_a;

        }
    }

    /*Calculo de minimo, maximo y promedio para B */
    int inicio = strip_size * identificador;
    int fin = inicio + strip_size;
    for(i = inicio; i < fin; i++){
        for(j = 0; j < n; j++){
            actual_b = b[i*n + j];
            if(actual_b > maximos[1]){
                maximos[1] = actual_b;
            }

            if(actual_b < minimos[1]){
                minimos[1] = actual_b;
            }

            suma_local_b += actual_b;
        }
    }
}

void sumar_vector_escalar(){
    int i, j, k;

    for(i = 0; i < strip_size; i++){
        for(j = 0; j < n; j++){
            r[i*n+j] = escalar * ab[i*n+j] + cbt[i*n+j];
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
