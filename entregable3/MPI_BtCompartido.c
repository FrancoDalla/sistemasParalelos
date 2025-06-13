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
#define DEFAULT_BS 32
#define MASTER 0

double dwalltime();
void transpose(double *m, double *mt, int n, int strip_size, int identificador);
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);
void recorrido_escalar(double *a, double *b,int n, int strip_size, int identificador,double *maximos,double *minimos,double *suma_local_a,double *suma_local_b);
void matmulblks(double *a, double *b, double *c, int n, int strip_size, int bs);
void sumar_vector_escalar(int strip_size, int n, double escalar, double *r, double *ab, double *cbt);


int main(int argc, char *argv[]){
    int i,j,k;
    /*Declaracion de variables para los argumentos */
    int n,verificacion;
    int bs = DEFAULT_BS;
    verificacion = 0;
    n = atoi(argv[1]);

    /*Declaración de matrices */
    double *a;
    double *b, *bt;
    double *c;
    double *ab;
    double *cbt;
    double *r;

    /*Declaración de valores asociados a la matrices */
    int matriz_tamano;
    int strip_size;

    /*Declaración de variables relacionadas al calculo del escalar */
    double suma_local_a = 0;
    double suma_local_b = 0;
    double promedio_a;
    double promedio_b;
    double escalar;
    int actual_a, actual_b;
    double maximos[2] = { DBL_MIN, DBL_MIN };
    double minimos[2] = { DBL_MAX, DBL_MAX };
    double maximo[2];
    double minimo[2];

    /*Declaración de variables para datos de MPI*/
    int procesos;
    int identificador;

    /*Declaración de variables por uso de comunicación asincrona*/
    MPI_Request requests[15];
    MPI_Status status;

    /*Declaración de variables para medición del tiempo de ejecución */
    double comm_times[20] = {0};
    double total_time;
    double comm_time = 0;
    double local_time = 0;

    if(argc > 2){
        verificacion = atoi(argv[2]);
    }

    matriz_tamano = n * n;

    promedio_a = 0; promedio_b = 0;

    MPI_Init(&argc, &argv);

        MPI_Comm_size(MPI_COMM_WORLD, &procesos);
        MPI_Comm_rank(MPI_COMM_WORLD, &identificador);
        strip_size = n / procesos;

        if((strip_size) < bs){
            bs = strip_size;
        }

        if(identificador == MASTER){
            /*Se reserva memoria para los datos del maestro */
            a = (double *) malloc(sizeof(double) * matriz_tamano);
            b = (double *) malloc(sizeof(double) * matriz_tamano);
            ab = (double *) malloc(sizeof(double) * matriz_tamano);
            bt = (double *) malloc(sizeof(double) * matriz_tamano);
            c = (double *) malloc(sizeof(double) * matriz_tamano);
            cbt = (double *) malloc(sizeof(double) * matriz_tamano);
            r = (double * ) malloc(sizeof(double) * matriz_tamano);

            /*Asignación de valores a las variables de maestro*/
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    a[i * n + j] = j;
                    c[i * n + j] = j;
                    b[i + j * n] = j;
                }
            }
        }
        else{
            /*Se reserva memoria para los datos de los procesos worker */
            a = (double *) malloc(sizeof(double) * strip_size * n);
            r = (double * ) malloc(sizeof(double) * strip_size * n);
            ab = (double *) malloc(sizeof(double) * strip_size * n);
            cbt = (double *) malloc(sizeof(double) * strip_size * n);
            c = (double *) malloc(sizeof(double) * strip_size * n);
            bt = (double *) malloc(sizeof(double) * matriz_tamano);
            b = (double *) malloc(sizeof(double) * matriz_tamano);
        }
        /*Se espera a que todos los procesos tengan los datos en condiciones */
        MPI_Barrier(MPI_COMM_WORLD);
        total_time = dwalltime();
        /*Tomamos el tiempo para el total */

        comm_times[0] = MPI_Wtime();
        MPI_Ibcast(b, matriz_tamano, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[0]);
        MPI_Iscatter(a, n * strip_size, MPI_DOUBLE, a, n * strip_size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[1]);
        MPI_Iscatter(c, n * strip_size, MPI_DOUBLE, c, n * strip_size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[2]);
        comm_times[1] = MPI_Wtime();

        comm_times[2] = MPI_Wtime();
        MPI_Wait(&requests[0], &status);
        comm_times[3] = MPI_Wtime();

        /*Cada proceso transpone la matriz en conjunto */
        transpose(b, bt, n, strip_size, identificador);

        comm_times[4] = MPI_Wtime();
        MPI_Iallgather(bt, strip_size*n, MPI_DOUBLE, bt, strip_size * n, MPI_DOUBLE, MPI_COMM_WORLD, &requests[9]);
        comm_times[5] = MPI_Wtime();

        recorrido_escalar(a, b, n, strip_size, identificador, maximos, minimos, &suma_local_a, &suma_local_b);

        comm_times[6] = MPI_Wtime();
        MPI_Ireduce(maximos, maximo, 2, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD, &requests[3]);
        MPI_Ireduce(minimos, minimo, 2, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD, &requests[4]);
        MPI_Ireduce(&suma_local_a, &promedio_a, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD, &requests[5]);
        MPI_Ireduce(&suma_local_b, &promedio_b, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD,&requests[6]);
        comm_times[7] = MPI_Wtime();

        /*Calculo de A X B */
        matmulblks(a, b, ab, n, strip_size, bs);

        /*CALCULO DE C * Bt */
        comm_times[8] = MPI_Wtime();
        MPI_Wait(&requests[9], &status);
        comm_times[9] = MPI_Wtime();
        matmulblks(c, bt, cbt, n, strip_size, bs);

        /*Esperar datos pendientes para el escalar */
        if(identificador == MASTER){
            comm_times[10] = MPI_Wtime();
            MPI_Waitall(4,&requests[3], MPI_STATUSES_IGNORE);
            comm_times[11] = MPI_Wtime();

            promedio_a = promedio_a / matriz_tamano;
            promedio_b = promedio_b / matriz_tamano;
            escalar = (maximo[0] * maximo[1] - minimo[0] * minimo[1])/(promedio_a * promedio_b);
        }

        comm_times[12] = MPI_Wtime();
        MPI_Ibcast(&escalar, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[7]);
        comm_times[13] = MPI_Wtime();

        /*Se espera a que los procesos tengan el escalar */
        comm_times[14] = MPI_Wtime();
        MPI_Wait(&requests[7], &status);
        comm_times[15] = MPI_Wtime();

        /*Calculo de  escalar a (AXB + CXBt)*/
        sumar_vector_escalar(strip_size, n, escalar, r, ab, cbt);

        comm_times[16] = MPI_Wtime();
        MPI_Igather(r, strip_size * n, MPI_DOUBLE, r, strip_size * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[8]);
        comm_times[17] = MPI_Wtime();

        if(identificador == MASTER){
            comm_times[18] = MPI_Wtime();
            MPI_Wait(&requests[8], &status);
            comm_times[19] = MPI_Wtime();
        }

        for(i = 0; i < 18; i++){
            local_time += comm_times[i+1] - comm_times[i];
        }

        MPI_Reduce(&local_time, &comm_time, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

    MPI_Finalize();

    if(identificador == MASTER){
        total_time = dwalltime() - total_time;
        comm_time = comm_time / procesos;
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
            printf("PARA N %d\n", n);
            printf("Los resultados en la matriz R son correctos\n");
        }

        printf("Tiempo de ejecucion: %f @@ Tiempo de comunicación: %f\n", total_time, comm_time);
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

void transpose(double *m, double *m_trans, int n, int strip_size, int identificador) {
    int i,j;
    int inicio = identificador * strip_size;
    int fin = inicio + strip_size;
    for (i = 0; i < n; i++){
        for(j = inicio; j < fin; j++){
            m_trans[i + j * n] = m[i * n + j];
        }
    }

}


void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs) {
    int i, j, k;
    for (i = 0; i < bs; i++) {
        for (j = 0; j < bs; j++){
            for(k = 0; k < bs; k++){
                cblk[i*n+j] += ablk[i*n+k] * bblk[k+j*n];
            }
        }
    }
}

void matmulblks(double *a, double *b, double *c, int n, int strip_size, int bs){
    int i, j, k, posicion_i, posicion_j;
    for(i = 0; i < strip_size; i+=bs){
        posicion_i = i * n;
        for(j = 0; j < n; j += bs){
            posicion_j = j * n;
            for(k = 0; k<n; k += bs){
                blkmul(&a[posicion_i + k], &b[posicion_j + k], &c[posicion_i + j], n, bs);
            }
        }
    }
}

void recorrido_escalar(
    double *a, double *b,
    int n, int strip_size, int identificador,
    double *maximos,
    double *minimos,
    double *suma_local_a,
    double *suma_local_b
    )
{
    int i,j;
    double actual_a, actual_b;

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

            *suma_local_a += actual_a;

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

            *suma_local_b += actual_b;
        }
    }
}

void sumar_vector_escalar(int strip_size, int n, double escalar, double *r, double *ab, double *cbt){
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
