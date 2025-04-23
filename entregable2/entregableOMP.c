#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <float.h>
#include <omp.h>

double dwalltime();
void matmulblks_and_calculate_scalar(double *a, double *b, double *c);
void blkmul_getscalar(double *ablk, double *bblk, double *cblk, int first_row, int first_column);
void transpose(double *m, double *mt);
void matmulblks(double *a, double *b, double *c);
void blkmul(double *ablk, double *bblk, double *cblk);

int n = 0, bs = 32;
double prom_a = 0.0, prom_b = 0.0, escalar = 0.0;
double max_a = -1, max_b = -1;
double min_a = 32000, min_b = 32000;

int main(int argc, char *argv[]){
    double timetick;
    int verificacion = 0;
    int i,j,t,k;

    /*Verificacion del input */
    if(argc < 3 || argc > 5){
        fprintf(stderr, "Uso %s N [BS] [verificacion]\n", argv[0]);
        fprintf(stderr, "N: Tamaño matriz.\n");
        fprintf(stderr, "T: Cantidad de hilos a utilizar\n");
        fprintf(stderr, "BS: Tamaño de los bloques, por defecto 32(x32).\n");
        fprintf(stderr, "verificacion: Si es 1 o mayor activa la verificación de los calculos.\n");
        exit(1);
    } else {
        n = atoi(argv[1]);
        t = atoi(argv[2]);
        if (argc >= 4)
            bs = atoi(argv[3]);
        if (argc == 5)
            verificacion = atoi(argv[4]);
    }

    if (n<=0) {
        fprintf(stderr, "N debe ser mayor o igual a 0.");
        exit(1);
    }
    if (bs<=0) {
        fprintf(stderr, "BS debe ser mayor o igual a 0.");
        exit(1);
    }
    if (n%bs!=0) {
        fprintf(stderr,"N debe ser divisible por BS.");
        exit(1);
    }

    /*Reserva de memoria para matrices */
    double * a = malloc(n*n * sizeof(double));
    double * b = malloc(n*n * sizeof(double));
    double * c = malloc(n*n * sizeof(double));
    double * r = malloc(n*n * sizeof(double));
    double * b_trans = malloc(n*n * sizeof(double));
    double * matriz_ab = malloc(n*n * sizeof(double));
    double * matriz_cbt = malloc(n*n * sizeof(double));

    omp_set_num_threads(t);
    /*Inicializacion de matrices */
    for(i = 0; i < n; i++){
        for(j = 0;j < n;j++){
            a[i*n+j] = j;
            c[i*n+j] = j;
            b[i+j*n] = j;
            matriz_ab[i*n+j] = 0.0;
            matriz_cbt[i*n+j] = 0.0;
            r[i*n+j] = 0.0;
        }
    }

    timetick = dwalltime();

    /*
    escalar = [(maxA * maxB - minA * minB) / (promA * promB)]
    R = (escalar * [A*B] + [C * bt])
    */
    #pragma omp parallel
    {
        matmulblks_and_calculate_scalar(a, b, c);
        transpose(b, b_trans);
        matmulblks(c,b_trans, matriz_cbt);

        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                r[i*n+j] = escalar * matriz_ab[i*n+j] + matriz_cbt[i*n+j];
            }
        }
    }


    timetick = dwalltime() - timetick;

    for(int i = 0; i < n*n; i++){
        printf("%g ",r[i]);
    }
    printf("\n");

    /*
        Verificación de resultados.
    */
    if (verificacion >= 1) {
        // Por se incializan los valores, el valor escalar siempre es 4.
        // escalar = (MaxA * MaxB - MinA * MinB) / (PromA * PromB).
        // escalar = ([n-1] * [n-1] - 0 * 0) / ([n-1]/2) * ([n-1]/2)).
        if (escalar != 4.0) {
            fprintf(stderr, "Error en el valor escalar, esperado 4 se calculó %g", escalar);
            exit(1);
        }
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                double ab = 0.0, cbt = 0.0;
                for (k = 0; k < n; k++){
                    ab += a[i*n+k] * b[k+j*n];
                    cbt += c[i*n+k] * b[i+k*n];
                }

                double expected = 4.0 * ab + cbt;
                if (r[i*n+j] != expected) {
                    fprintf(
                        stderr,
                        "Error en posición (%d,%d). Esperado valor: (%g), se obtuvo (%g).\n",
                        i,j,expected, r[i*n+j]
                    );
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


/* Multiplicación de matrices utilizando bloques. Para A x B y conseguir el valor escalar. */
void matmulblks_and_calculate_scalar(double *a, double *b, double *c){
    int i, j, k;
    for(i = 0; i < n; i += bs){
        for(j = 0; j < n;j += bs){
            for(k = 0; k < n; k += bs) {
                blkmul_getscalar(&a[i*n+k], &b[k+j*n], &c[i*n+j], j==0, i==0);
            }
        }
    }
    prom_a = prom_a / (n*n);
    prom_b = prom_b / (n*n);
    escalar = (max_a * max_b - min_a * min_b) / (prom_a * prom_b);
}

/* Multiplicación de los bloques y obtener valores para el valor escalar */
void blkmul_getscalar(double *ablk, double *bblk, double *cblk, int first_row, int first_column) {
    int i, j, k;
    for (i = 0; i < bs; i++) {
        for (j = 0; j < bs; j++){
            // Obtención de MinA, MaxA, PromA
            if (first_row) {
                double val_a = ablk[i*n+j];
                if (val_a > max_a) max_a = val_a;
                if (val_a < min_a) min_a = val_a;
                prom_a += val_a;
            }
            // Obtención de MinB, MaxB, PromB
            if (first_column) {
                double val_b = bblk[i+j*n];
                if (val_b > max_b) max_b = val_b;
                if (val_b < min_b) min_b = val_b;
                prom_b += val_b;
            }

            for(k = 0; k < bs; k++){
                // Calculo A*B
                cblk[i*n+j] += ablk[i*n+k] * bblk[k+j*n];
            }
        }
    }
}

void transpose(double *m, double *m_trans) {
    int j, i;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            m_trans[i + j * n] = m[i * n + j];
        }
    }
}

/* Multiplicación de matrices utilizando bloques. Para B_t x C */
void matmulblks(double *a, double *b, double *c) {
    int i, j, k;
    for(i = 0; i < n; i += bs){
        for(j = 0; j < n; j += bs){
            for(k = 0; k < n; k += bs) {
                blkmul(&a[i*n+k],&b[k+j*n],&c[i*n+j]);
            }
        }
    }
}

/* Multiplicación de los bloques y obtener valores para el valor escalar */
void blkmul(double *ablk, double *bblk, double *cblk) {
    int i, j, k;
    for (i = 0; i < bs; i++) {
        for (j = 0; j < bs; j++){
            for(k = 0; k < bs; k++){
                cblk[i*n+j] += ablk[i*n+k] * bblk[k+j*n];
            }
        }
    }
}

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}
