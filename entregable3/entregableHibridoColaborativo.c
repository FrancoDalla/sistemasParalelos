#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>

#define DBL_MAX 2.2250738585072014e-308
#define DBL_MIN -2.2250738585072014e-308
#define BS 8
#define COORDINATOR 0

double dwalltime();
void blkmul(double *ablk, double *bblk, double *cblk);

int n = 0;

int main(int argc, char *argv[])
{
	/*
        Definición de variables.
    */
	int i, j, k, ii, jj, kk, numProcs, numThreads, rank, stripSize,
		provided;

	double *a, *c, *r, *matriz_ab, *matriz_cbt, *b_trans, *b;

	double val_a, val_b;
	double prom_a, prom_b;
	double localCommTime = 0.0, commTime, totalTime;
	double max[2];
	double min[2];

	double sum_a = 0.0, sum_b = 0.0;
	double local_max[2] = { DBL_MIN, DBL_MIN };
	double local_min[2] = { DBL_MAX, DBL_MAX };

	double escalar;

	int verificacion = 0;
	MPI_Status status;
	MPI_Request requests[9];
	// EN TEORÍA c debería hacer que todo sea 0 declarandolo de esta forma.
	double commTimes[26] = { 0 };

	if (argc < 3 || argc > 4) {
		fprintf(stderr, "Uso %s T N [verificacion]\n", argv[0]);
		fprintf(stderr, "   T: Cantidad de hilos OMP.\n");
		fprintf(stderr, "   N: Tamaño de las matrices (NxN).\n");
		fprintf(stderr,
			"   verificacion: Con >=1 activa la verificación.\n");
		exit(1);
	} else {
		numThreads = atoi(argv[1]);
		n = atoi(argv[2]);
		if (argc == 4)
			verificacion = atoi(argv[3]);
	}
	// Verificación de input

	if (n <= 0) {
		fprintf(stderr, "N debe ser mayor o igual a 0.\n");
		exit(1);
	}
	if (n % BS != 0) {
		fprintf(stderr, "N debe ser divisible por BS.\n");
		exit(1);
	}

	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

	if (provided < MPI_THREAD_FUNNELED) {
		fprintf(stderr,
			"Threading provisto por MPI no es suficiente\n");
		exit(1);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if ((n % numProcs != 0) && ((n / numProcs) % numThreads != 0)) {
		printf("El tamaño de la matriz debe ser multiplo del numero de procesos, y esta de los hilos.\n");
		exit(1);
	}

	omp_set_num_threads(numThreads);

	// calcular porcion de cada worker
	stripSize = n / numProcs;

	// Reservar memoria
	if (rank == COORDINATOR) {
		a = malloc(n * n * sizeof(double));
		c = malloc(n * n * sizeof(double));
		r = malloc(n * n * sizeof(double));
		matriz_ab = malloc(n * n * sizeof(double));
		matriz_cbt = malloc(n * n * sizeof(double));
	} else {
		a = malloc(n * stripSize * sizeof(double));
		c = malloc(n * stripSize * sizeof(double));
		r = malloc(n * stripSize * sizeof(double));
		matriz_ab = malloc(n * stripSize * sizeof(double));
		matriz_cbt = malloc(n * stripSize * sizeof(double));
	}
	b_trans = malloc(n * n * sizeof(double));
	b = malloc(n * n * sizeof(double));

	/*
        Inicialización de matrices.
    */
	if (rank == COORDINATOR) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				a[i * n + j] = j;
				c[i * n + j] = j;
				b[i + j * n] = j;
				b_trans[i + j * n] = 0.0;
				matriz_ab[i * n + j] = 0.0;
				matriz_cbt[i * n + j] = 0.0;
				r[i * n + j] = 0.0;
			}
		}
	} else {
		for (i = 0; i < stripSize; i++) {
			for (j = 0; j < n; j++) {
				a[i * n + j] = 0.0;
				c[i * n + j] = 0.0;
				matriz_ab[i * n + j] = 0.0;
				matriz_cbt[i * n + j] = 0.0;
				r[i * n + j] = 0.0;
			}
		}
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				b[i + j * n] = 0.0;
				b_trans[i + j * n] = 0.0;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	commTimes[0] = MPI_Wtime();

	/* distribuir datos*/
	// El orden de la comunicación es importante para esperar lo menos posible
	MPI_Ibcast(b, n * n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD,
		   &requests[0]);
	MPI_Iscatter(a, n * stripSize, MPI_DOUBLE, a, n * stripSize, MPI_DOUBLE,
		     COORDINATOR, MPI_COMM_WORLD, &requests[1]);
	MPI_Iscatter(c, n * stripSize, MPI_DOUBLE, c, n * stripSize, MPI_DOUBLE,
		     COORDINATOR, MPI_COMM_WORLD, &requests[2]);

	commTimes[1] = MPI_Wtime();

#pragma omp parallel default(shared) private(i, j, k)
	{ // Inicio de la región paralela (OMP)
// Transpuesta
#pragma omp barrier
#pragma omp master
		if (rank != COORDINATOR) {
			commTimes[2] = MPI_Wtime();
			MPI_Wait(&requests[0], &status);
			commTimes[3] = MPI_Wtime();
		}
#pragma omp barrier
#pragma omp for schedule(static)
		for (i = 0; i < n; i++) {
			for (j = (rank * stripSize);
			     j < (rank * stripSize + stripSize); j++) {
				b_trans[i + j * n] = b[i * n + j];
			}
		}
		// Comunicación de transpuesta
#pragma omp barrier
#pragma omp master
		{
			commTimes[4] = MPI_Wtime();
			MPI_Iallgather(b_trans, stripSize * n, MPI_DOUBLE,
				       b_trans, stripSize * n, MPI_DOUBLE,
				       MPI_COMM_WORLD, &requests[3]);
			commTimes[5] = MPI_Wtime();
		}
#pragma omp barrier

/*Obtención de MinA, MinB, MaxA, MaxB, PromA, PromB */
#pragma omp barrier
#pragma omp master
		if (rank != COORDINATOR) {
			commTimes[6] = MPI_Wtime();
			MPI_Wait(&requests[1], &status);
			commTimes[7] = MPI_Wtime();
		}
#pragma omp barrier
#pragma omp for schedule(static) reduction(+ : sum_a) \
	reduction(max : local_max[0]) reduction(min : local_min[0])
		for (i = 0; i < stripSize; i++) {
			for (j = 0; j < n; j++) {
				val_a = a[i * n + j];
				sum_a += val_a;
				local_max[0] = val_a > local_max[0] ?
						       val_a :
						       local_max[0];
				local_min[0] = val_a < local_min[0] ?
						       val_a :
						       local_min[0];
			}
		}
#pragma omp barrier
#pragma omp master
		if (rank != COORDINATOR) {
			commTimes[8] = MPI_Wtime();
			MPI_Wait(&requests[0], &status);
			commTimes[9] = MPI_Wtime();
		}
#pragma omp barrier
#pragma omp for schedule(static) reduction(+ : sum_b) \
	reduction(max : local_max[1]) reduction(min : local_min[1])
		for (i = (stripSize * rank); i < (stripSize * rank) + stripSize;
		     i++) {
			for (j = 0; j < n; j++) {
				val_b = b[i * n + j];
				sum_b += val_b;
				local_max[1] = val_b > local_max[1] ?
						       val_b :
						       local_max[1];
				local_min[1] = val_b < local_min[1] ?
						       val_b :
						       local_min[1];
			}
		}
		// Comunicación de valores para el estadístico
#pragma omp barrier
#pragma omp master
		{
			commTimes[10] = MPI_Wtime();
			// Aparentemente MPI_SUM no funciona para vectores así que lo descompongo en dos reducciones.
			MPI_Ireduce(&sum_a, &prom_a, 1, MPI_DOUBLE, MPI_SUM,
				    COORDINATOR, MPI_COMM_WORLD, &requests[4]);
			MPI_Ireduce(&sum_b, &prom_b, 1, MPI_DOUBLE, MPI_SUM,
				    COORDINATOR, MPI_COMM_WORLD, &requests[5]);
			MPI_Ireduce(local_max, max, 2, MPI_DOUBLE, MPI_MAX,
				    COORDINATOR, MPI_COMM_WORLD, &requests[6]);
			MPI_Ireduce(local_min, min, 2, MPI_DOUBLE, MPI_MIN,
				    COORDINATOR, MPI_COMM_WORLD, &requests[7]);
			commTimes[11] = MPI_Wtime();
		}
#pragma omp barrier

/*Multiplicacion a * b */
#pragma omp for schedule(static)
		for (i = 0; i < stripSize; i += BS) {
			for (j = 0; j < n; j += BS) {
				for (k = 0; k < n; k += BS) {
					blkmul(&a[i * n + k], &b[k + j * n],
					       &matriz_ab[i * n + j]);
				}
			}
		}

		// Esperar a las reducciones y hacer broadcast.

#pragma omp barrier
#pragma omp master
		{
			if (rank == COORDINATOR) {
				commTimes[12] = MPI_Wtime();
				MPI_Wait(&requests[4], &status);
				MPI_Wait(&requests[5], &status);
				commTimes[13] = MPI_Wtime();
				prom_a /= (n * n);
				prom_b /= (n * n);
				commTimes[14] = MPI_Wtime();
				MPI_Wait(&requests[6], &status);
				MPI_Wait(&requests[7], &status);
				commTimes[15] = MPI_Wtime();
				escalar = (max[0] * max[1] - min[0] * min[1]) /
					  (prom_a * prom_b);
			}
			commTimes[16] = MPI_Wtime();
			MPI_Ibcast(&escalar, 1, MPI_DOUBLE, COORDINATOR,
				   MPI_COMM_WORLD, &requests[8]);
			commTimes[17] = MPI_Wtime();
		}
#pragma omp barrier

		// Asegurarse que se tiene B Transpuesta y C.
#pragma omp barrier
#pragma omp master
		{
			commTimes[18] = MPI_Wtime();
			MPI_Wait(&requests[3], &status);
			commTimes[19] = MPI_Wtime();
			if (rank != COORDINATOR) {
				commTimes[20] = MPI_Wtime();
				MPI_Wait(&requests[2], &status);
				commTimes[21] = MPI_Wtime();
			}
		}
#pragma omp barrier
		/*Calculo de c * b_transpuesta */
#pragma omp for schedule(static)
		for (i = 0; i < stripSize; i += BS) {
			for (j = 0; j < n; j += BS) {
				for (k = 0; k < n; k += BS) {
					blkmul(&c[i * n + k],
					       &b_trans[k + j * n],
					       &matriz_cbt[i * n + j]);
				}
			}
		}

//Asegurar que se tiene el escalar
#pragma omp barrier
#pragma omp master
		if (rank != COORDINATOR) {
			commTimes[22] = MPI_Wtime();
			MPI_Wait(&requests[8], &status);
			commTimes[23] = MPI_Wtime();
		}
#pragma omp barrier
/*Calculo de R */
#pragma omp for schedule(static)
		for (i = 0; i < stripSize; i++) {
			for (j = 0; j < n; j++) {
				r[i * n + j] = escalar * matriz_ab[i * n + j] +
					       matriz_cbt[i * n + j];
			}
		}
	} // Fin de la región paralela (OMP)

	// Se guarda R en el COORDINADOR (no se puede usar no-bloqueante, porque si no puede terminarse antes de tiempo el programa).
	commTimes[24] = MPI_Wtime();
	MPI_Gather(r, stripSize * n, MPI_DOUBLE, r, stripSize * n, MPI_DOUBLE,
		   COORDINATOR, MPI_COMM_WORLD);
	commTimes[25] = MPI_Wtime();

	// Calcular el tiempo total y de las comunicaciones
	for (i = 0; i <= 24; i += 2) {
		localCommTime += commTimes[i + 1] - commTimes[i];
	}

	MPI_Reduce(&localCommTime, &commTime, 1, MPI_DOUBLE, MPI_SUM,
		   COORDINATOR, MPI_COMM_WORLD);

	MPI_Finalize();

	if (rank == COORDINATOR) {
		totalTime = commTimes[25] - commTimes[0];
		commTime /= numProcs; // Promedio

		/*
        Verificación de resultados.
    */
		if (verificacion >= 1) {
			if (escalar != 4.0) {
				fprintf(stderr,
					"Error en el valor escalar, esperado 4 se calculó %lf\n",
					escalar);
				exit(1);
			}

			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (b_trans[i + j * n] != (double)i) {
						fprintf(stderr,
							"B_t Error en posición (%d,%d). Esperado valor: (%g), se obtuvo (%g).\n",
							i, j, (double)i,
							b_trans[i + j * n]);
						exit(1);
					}
				}
			}

			// Mejora de verificación en base al código original:
			// Debido a la forma en que se hace la inicialización, se sabe de antemano los valores correctos
			// ab = j * (n-1)*n/2   (suma de 0 a n-1 de k * j)
			// cbt = (n-1)*n*(2*n-1)/6   (suma de 0 a n-1 de k * k)
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					double ab = j * (n - 1.0) * n / 2.0;
					double cbt = (n - 1.0) * n *
						     (2.0 * n - 1.0) / 6.0;
					double expected = 4.0 * ab + cbt;

					if (r[i * n + j] != expected) {
						fprintf(stderr,
							"R = Error en posición (%d,%d). Esperado valor: (%g), se obtuvo (%g).\n",
							i, j, expected,
							r[i * n + j]);
						exit(1);
					}
				}
			}
		}

		printf("Tiempo de ejecución del calculo: %lf | Tiempo de comunicaciones: %lf \n",
		       totalTime, commTime);
	}

	// Liberando memoria por si acaso.
	free(a);
	free(b);
	free(c);
	free(r);
	free(b_trans);
	free(matriz_ab);
	free(matriz_cbt);

	return 0;
}

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

/* Multiplicación de los bloques y obtener valores para el valor escalar */
void blkmul(double *ablk, double *bblk, double *cblk)
{
	int i, j, k;
	for (i = 0; i < BS; i++) {
		for (j = 0; j < BS; j++) {
			for (k = 0; k < BS; k++) {
				cblk[i * n + j] +=
					ablk[i * n + k] * bblk[k + j * n];
			}
		}
	}
}
