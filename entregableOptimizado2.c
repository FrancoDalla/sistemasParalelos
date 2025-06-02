#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define DBL_MAX 2.2250738585072014e-308
#define DBL_MIN -2.2250738585072014e-308
#define BS 128

double dwalltime();
void transpose(double *m, double *mt);
void blkmul(double *ablk, double *bblk, double *cblk);

int n = 0;

int main(int argc, char *argv[]) {

    /*
        Definición de variables.
    */
    int i = 0, j = 0, k = 0, ii = 0, jj = 0, kk = 0;
    double timetick;
    int verificacion = 0;
    double prom_a = 0.0, prom_b = 0.0, escalar = 0.0;
    double max_a = DBL_MIN, max_b = DBL_MIN;
    double min_a = DBL_MAX, min_b = DBL_MAX;
    double *a, *b, *c, *r, *b_trans, *matriz_ab, *matriz_cbt;

    if(argc < 2 || argc > 3){
        fprintf(stderr, "Uso %s N [verificacion]\n", argv[0]);
        fprintf(stderr, "   N: Tamaño de la matriz NxN.\n");
        fprintf(stderr, "   verificacion: Con >=1 activa la verificación.\n");
        exit(1);
    } else {
        n = atoi(argv[1]);
        if (argc == 3)
            verificacion = atoi(argv[2]);
    }

    // Verificación de input
    if (n<=0) {
        fprintf(stderr, "N debe ser mayor o igual a 0.");
        exit(1);
    }
    if (n%BS!=0) {
        fprintf(stderr,"N debe ser divisible por BS.");
        exit(1);
    }

    // Alojamiento de memoria
    a = malloc(n*n * sizeof(double));
    b = malloc(n*n * sizeof(double));
    c = malloc(n*n * sizeof(double));
    r = malloc(n*n * sizeof(double));
    b_trans = malloc(n*n * sizeof(double));
    matriz_ab = malloc(n*n * sizeof(double));
    matriz_cbt = malloc(n*n * sizeof(double));

    /*
        Inicialización de matrices.
    */
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
        escalar = [(MaxA * MaxB - MinA * MinB) / (PromA * PromB)].
        R = ( escalar * [A*B] + [C*B_t].
    */

    {
        double val_a;
        double val_b;
        /*Obtención de MinA, MinB, MaxA, MaxB, PromA, PromB */
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                val_a = a[i*n+j];

                prom_a += val_a;
                if(val_a > max_a) max_a = val_a;
                if(val_a < min_a) min_a = val_a;
            }
        }
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                val_b = b[i*n+j];

                prom_b += val_b;
                if(val_b > max_b) max_b = val_b;
                if(val_b < min_b) min_b = val_b;
            }
        }

        prom_a = prom_a / (n*n);
        prom_b = prom_b / (n*n);
        escalar = (max_a * max_b - min_a * min_b) / (prom_a * prom_b);

        /*Multiplicacion a * b */
        for(i = 0; i < n; i += BS){
            for(j = 0; j < n; j += BS){
                for(k = 0; k < n; k += BS){
                    blkmul(&a[i*n+k], &b[k+j*n],&matriz_ab[i*n+j]);
                }
            }
        }
        // Calcular la transposición de B
        transpose(b, b_trans);


        /*Calculo de c * b_transpuesta */
        for(i = 0; i < n; i += BS){
            for(j = 0; j < n; j += BS){
                for(k = 0; k < n; k += BS){
                    blkmul(&c[i*n+k], &b_trans[k+j*n], &matriz_cbt[i*n+j]);
                }
            }
        }

        /*Calculo de R */
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                r[i*n+j] = escalar * matriz_ab[i*n+j] + matriz_cbt[i*n+j];
            }
        }
    }

    timetick = dwalltime() - timetick;

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

        // Mejora de verificación en base al código original:
        // Debido a la forma en que se hace la inicialización, se sabe de antemano los valores correctos
        // ab = j * (n-1)*n/2   (suma de 0 a n-1 de k * j)
        // cbt = (n-1)*n*(2*n-1)/6   (suma de 0 a n-1 de k * k)
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                double ab = j * (n-1.0)*n/2.0;
                double cbt = (n-1.0)*n*(2.0*n-1.0)/6.0;
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

/* Multiplicación de los bloques y obtener valores para el valor escalar */
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
