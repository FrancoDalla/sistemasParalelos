#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define DBL_MAX 1.7976931348623157e+308
#define DBL_MIN 2.2250738585072014e-308

double dwalltime();

int main(int argc, char *argv[]) {

    /*
        Definición de variables.
    */
    int i = 0, n = 0, j = 0, k = 0, bs = 32, ii = 0, jj = 0, kk = 0;
    double prom_a = 0.0, prom_b = 0.0, escalar = 0.0;
    double max_a = DBL_MIN, max_b = DBL_MIN;
    double min_a = DBL_MAX, min_b = DBL_MAX;
    double timetick;
    int verificacion = 0;

    if(argc < 2 || argc > 4){
        fprintf(stderr, "Uso %s N [BS] [verificacion]\n", argv[0]);
        fprintf(stderr, "N: Tamaño matriz.\n");
        fprintf(stderr, "BS: Tamaño de los bloques, por defecto 32(x32).\n");
        fprintf(stderr, "verificacion: Si es 1 o mayor activa la verificación de los calculos.\n");
        exit(1);
    } else {
        n = atoi(argv[1]);
        if (argc >= 3)
            bs = atoi(argv[2]);
        if (argc == 4)
            verificacion = atoi(argv[3]);
    }
    // Verificación de input

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

    double * a = malloc(n*n * sizeof(double));
    double * b = malloc(n*n * sizeof(double));
    double * c = malloc(n*n * sizeof(double));
    double * r = malloc(n*n * sizeof(double));
    double * b_trans = malloc(n*n * sizeof(double));
    double * matriz_ab = malloc(n*n * sizeof(double));
    double * matriz_cbt = malloc(n*n * sizeof(double));

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


    /*
        escalar = [(MaxA * MaxB - MinA * MinB) / (PromA * PromB)].
        R = ( escalar * [A*B] + [C*B_t].
    */
    // Multiplicación de matrices e inicialización de matriz B_t.
    // Cálculo de A*B + Obtención de MinA, MinB, MaxA, MaxB, PromA, PromB + Transposición de B
    timetick = dwalltime();

    for(i = 0; i < n; i += bs){
        for(j = 0; j < n;j += bs){
            for(k = 0; k < n; k += bs) {
                double *ablk = &a[i*n+k], *bblk = &b[k+j*n], *b_transblk = &b_trans[j+k*n], *matriz_abblk = &matriz_ab[i*n+j];
                for (ii = 0; ii < bs; ii++) {
                    for (jj = 0; jj < bs; jj++){
                        // Obtención de MinA, MinB, MaxA, MaxB, PromA, PromB
                        if (j==0) {
                            double val_a = ablk[ii*n+jj];
                            if (val_a > max_a) max_a = val_a;
                            if (val_a < min_a) min_a = val_a;
                            prom_a += val_a;
                        }
                        if (i==0) {
                            double val_b = bblk[ii+jj*n];
                            if (val_b > max_b) max_b = val_b;
                            if (val_b < min_b) min_b = val_b;
                            prom_b += val_b;
                        }

                        // Calcular transposición de B
                        b_transblk[jj + ii*n] = bblk[ii + jj*n]; // Inicializar de manera ineficiente para usar de manera eficiente.

                        for(kk = 0; kk < bs; kk++){
                            // Calculo A*B
                            matriz_abblk[ii*n+jj] += ablk[ii*n+kk] * bblk[kk+jj*n];
                        }
                    }
                }            }
        }
    }
    prom_a = prom_a / (n*n);
    prom_b = prom_b / (n*n);
    escalar = (max_a * max_b - min_a * min_b) / (prom_a * prom_b);

    // Cálculo de C*B_t.
    for(i = 0; i < n; i += bs){
        for(j = 0; j < n; j += bs){
            for(k = 0; k < n; k += bs) {
                double *cblk = &c[i*n+k], *b_transblk = &b_trans[k+j*n], *matriz_cbtblk = &matriz_cbt[i*n+j];
                for (ii = 0; ii < bs; ii++) {
                    for (jj = 0; jj < bs; jj++){
                        for(kk = 0; kk < bs; kk++){
                            matriz_cbtblk[ii*n+jj] += cblk[ii*n+kk] * b_transblk[kk+jj*n];
                        }
                    }
                }
            }
        }
    }

    // Cálculo de R.
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            r[i*n+j] = escalar * matriz_ab[i*n+j] + matriz_cbt[i*n+j];
        }
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

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}
