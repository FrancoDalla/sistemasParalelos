#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <semaphore.h>

#define MAX_THREADS 100
#define DBL_MAX 1.7976931348623157e+308
#define BS 32

double dwalltime();
void *laburar(void *ptr);
void blkmul(double *ablk, double *bblk, double *cblk);

// Como se nos indicó, declaramos variables globales para evitar usar parámetros.
int n = 0, t, cont = 0;
double prom_a = 0.0, prom_b = 0.0, escalar = 0.0;
double max_a = -DBL_MAX, max_b = -DBL_MAX;
double min_a = DBL_MAX, min_b = DBL_MAX;
pthread_mutex_t prom_a_mutex, prom_b_mutex, max_a_mutex, max_b_mutex,
	min_a_mutex, min_b_mutex, cont_escalar_mutex;
double *a, *b, *c, *b_trans, *r, *matriz_ab, *matriz_cbt;
pthread_barrier_t barrera_transpuesta_escalar, barrera_matriz_cbt;

typedef struct {
	int start;
	int job_size;
} job;

int main(int argc, char *argv[])
{
	/*
        Definición de variables.
    */
	int i = 0, j = 0, k = 0;
	double timetick;
	int verificacion = 0;

	if (argc < 2 || argc > 4) {
		fprintf(stderr, "Uso %s N [T] [verificacion]\n", argv[0]);
		fprintf(stderr, "N: Tamaño matriz.\n");
		fprintf(stderr, "T: Cantidad de threads.");
		fprintf(stderr,
			"verificacion: Si es 1 o mayor activa la verificación de los calculos.\n");
		exit(1);
	} else {
		n = atoi(argv[1]);
		t = atoi(argv[2]);
		if (argc == 4) {
			verificacion = atoi(argv[3]);
		}
	}
	// Verificación de input

	if (n <= 0) {
		fprintf(stderr, "N debe ser mayor o igual a 0.");
		exit(1);
	}
	if (t <= 0) {
		fprintf(stderr, "T debe ser mayor o igual a 0.");
		exit(1);
	}
	if (n % BS != 0) {
		fprintf(stderr, "N debe ser divisible por %d.", BS);
		exit(1);
	}

	job thread_args[t];
	pthread_attr_t attr;
	pthread_t threads[t];
	int *status;
	int jobs = (n * n) / (BS * BS);

	a = malloc(n * n * sizeof(double));
	b = malloc(n * n * sizeof(double));
	c = malloc(n * n * sizeof(double));
	r = malloc(n * n * sizeof(double));
	b_trans = malloc(n * n * sizeof(double));
	matriz_ab = malloc(n * n * sizeof(double));
	matriz_cbt = malloc(n * n * sizeof(double));

	/*
        Inicialización de matrices.
    */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i * n + j] = j;
			c[i * n + j] = j;
			b[i + j * n] = j;
			matriz_ab[i * n + j] = 0.0;
			matriz_cbt[i * n + j] = 0.0;
			r[i * n + j] = 0.0;
		}
	}

	timetick = dwalltime();

	pthread_attr_init(&attr);
	pthread_mutex_init(&prom_a_mutex, NULL);
	pthread_mutex_init(&prom_b_mutex, NULL);
	pthread_mutex_init(&max_a_mutex, NULL);
	pthread_mutex_init(&max_b_mutex, NULL);
	pthread_mutex_init(&min_a_mutex, NULL);
	pthread_mutex_init(&min_b_mutex, NULL);
	pthread_mutex_init(&cont_escalar_mutex, NULL);
	pthread_barrier_init(&barrera_transpuesta_escalar, NULL, t);
	pthread_barrier_init(&barrera_matriz_cbt, NULL, t);

	// Cálculo de A*B + Obtención de MinA, MinB, MaxA, MaxB, PromA, PromB.
	/* Crea los hilos */
	for (i = 0; i < t; i++) {
		thread_args[i].start =
			i * (jobs / t) + ((i < jobs % t) ? i : jobs % t);
		thread_args[i].job_size = i < (jobs % t) ? jobs / t + 1 :
							   jobs / t;
		pthread_create(&threads[i], &attr, laburar, &thread_args[i]);
	}

	for (i = 0; i < t; i++) {
		pthread_join(threads[i], NULL);
	}

	timetick = dwalltime() - timetick;

	/*
        Verificación de resultados por método lento.
    */
	if (verificacion >= 1) {
		// Por se incializan los valores, el valor escalar siempre es 4.
		// escalar = (MaxA * MaxB - MinA * MinB) / (PromA * PromB).
		// escalar = ([n-1] * [n-1] - 0 * 0) / ([n-1]/2) * ([n-1]/2)).
		if (escalar != 4.0) {
			fprintf(stderr,
				"Error en el valor escalar, esperado 4 se calculó %g",
				escalar);
			exit(1);
		}
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				double ab = 0.0, cbt = 0.0;
				for (k = 0; k < n; k++) {
					ab += a[i * n + k] * b[k + j * n];
					cbt += c[i * n + k] * b[i + k * n];
				}

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

	printf("Tiempo de ejecución del calculo: %f\n", timetick);

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

void *laburar(void *ptr)
{
	int i, j, k, ii, jj, kk;
	job args = *(job *)ptr;
	int start_idx = args.start;
	int job_size = args.job_size;
	int end_idx = start_idx + job_size;
	int blocks_per_row = n / BS;

	double local_prom_a = 0.0, local_prom_b = 0.0;
	double local_max_a = -DBL_MAX, local_max_b = -DBL_MAX;
	double local_min_a = DBL_MAX, local_min_b = DBL_MAX;

	for (int idx = start_idx; idx < end_idx; idx++) {
		// Calcular indices de la matriz a calcular
		i = (idx / blocks_per_row) * BS;
		j = (idx % blocks_per_row) * BS;

		for (int k = 0; k < n; k += BS) {
			for (ii = 0; ii < BS; ii++) {
				for (jj = 0; jj < BS; jj++) {
					// Obtención de MinA, MaxA, PromA
					if (j == 0) {
						double val_a = a[(i + ii) * n +
								 (k + jj)];
						if (val_a > local_max_a)
							local_max_a = val_a;
						if (val_a < local_min_a)
							local_min_a = val_a;
						local_prom_a += val_a;
					}
					// Obtención de MinB, MaxB, PromB
					if (i == 0) {
						double val_b = b[(j + jj) * n +
								 (k + ii)];
						if (val_b > local_max_b)
							local_max_b = val_b;
						if (val_b < local_min_b)
							local_min_b = val_b;
						local_prom_b += val_b;
					}

					for (kk = 0; kk < BS; kk++) {
						// Calculo A*B
						matriz_ab[(i + ii) * n +
							  (j + jj)] +=
							a[(i + ii) * n +
							  (k + kk)] *
							b[(j + jj) * n +
							  (k + kk)];
					}
				}
			}
		}
	}

	// Aportar a los resultados
	pthread_mutex_lock(&max_a_mutex);
	if (local_max_a > max_a)
		max_a = local_max_a;
	pthread_mutex_unlock(&max_a_mutex);
	pthread_mutex_lock(&min_a_mutex);
	if (local_min_a < min_a)
		min_a = local_min_a;
	pthread_mutex_unlock(&min_a_mutex);
	pthread_mutex_lock(&prom_a_mutex);
	prom_a += local_prom_a;
	pthread_mutex_unlock(&prom_a_mutex);

	pthread_mutex_lock(&max_b_mutex);
	if (local_max_b > max_b)
		max_b = local_max_b;
	pthread_mutex_unlock(&max_b_mutex);
	pthread_mutex_lock(&min_b_mutex);
	if (local_min_b < min_b)
		min_b = local_min_b;
	pthread_mutex_unlock(&min_b_mutex);
	pthread_mutex_lock(&prom_b_mutex);
	prom_b += local_prom_b;
	pthread_mutex_unlock(&prom_b_mutex);

	pthread_mutex_lock(&cont_escalar_mutex);
	cont++;
	// Si todos terminaron de sumar los promedios, calcular escalar
	if (cont == t) {
		prom_a = prom_a / (n * n);
		prom_b = prom_b / (n * n);
		escalar = (max_a * max_b - min_a * min_b) / (prom_a * prom_b);
	}
	pthread_mutex_unlock(&cont_escalar_mutex);

	// Transponer b en b_trans
	for (int idx = start_idx; idx < end_idx; idx++) {
		int ii, jj, kk;
		// Calcular indices de la matriz a calcular
		int i_blk = (idx / blocks_per_row) * BS;
		int j_blk = (idx % blocks_per_row) * BS;

		for (ii = 0; ii < BS; ii++) {
			for (jj = 0; jj < BS; jj++) {
				i = i_blk + ii;
				j = j_blk + jj;
				b_trans[j * n + i] = b[i * n + j];
			}
		}
	}

	// Barrera para esperar la transpuesta de b y el escalar
	pthread_barrier_wait(&barrera_transpuesta_escalar);

	for (int idx = start_idx; idx < end_idx; idx++) {
		// Calcular indices de la matriz a calcular
		int ii = (idx / blocks_per_row) * BS;
		int jj = (idx % blocks_per_row) * BS;

		for (int k = 0; k < n; k += BS) {
			blkmul(&c[ii * n + k], &b_trans[jj * n + k],
			       &matriz_cbt[ii * n + jj]);
		}
	}

	pthread_barrier_wait(&barrera_matriz_cbt);

	for (int idx = start_idx; idx < end_idx; idx++) {
		int ii, jj, kk;
		// Calcular indices de la matriz a calcular
		int i_blk = (idx / blocks_per_row) * BS;
		int j_blk = (idx % blocks_per_row) * BS;

		for (ii = 0; ii < BS; ii++) {
			for (jj = 0; jj < BS; jj++) {
				i = i_blk + ii;
				j = j_blk + jj;
				r[i * n + j] = escalar * matriz_ab[i * n + j] +
					       matriz_cbt[i * n + j];
			}
		}
	}

	pthread_exit(0);
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

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}
