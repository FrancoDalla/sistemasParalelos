struct job {
    int start;
    int job_size;
};

int main(int argc, char *argv[])
{
    // (...)
    timetick = dwalltime();

    int jobs = (n * n) / (bs * bs);
    pthread_attr_init(&attr);

    /* Crea los hilos */
    for (i = 0; i < t; i++) {
        thread_args[i].start =
            i * (jobs / t) + ((i < jobs % t) ? i : jobs % t);
        thread_args[i].job_size = i < (jobs % t) ? jobs / t + 1 :
                               jobs / t;
        pthread_create(&threads[i], &attr, matmulblks, &thread_args[i]);
    }

    /* Espera a que los hilos terminen */
    for (i = 0; i < t; i++) {
        pthread_join(threads[i], (void *)&status);
    }

    double workTime = dwalltime() - timetick;
    // (...)
}

/* Multiply square matrices, blocked version */
void *matmulblks(void *ptr)
{
    struct job args = *(struct job *)ptr;
    int start_idx = args.start;
    int job_size = args.job_size;
    int end_idx = start_idx + job_size;
    int blocks_per_row = n / bs;

    for (int idx = start_idx; idx < end_idx; idx++) {
        // Calcular indices de la matriz a calcular
        int i = (idx / blocks_per_row) * bs;
        int j = (idx % blocks_per_row) * bs;

        for (int k = 0; k < n; k += bs) {
            blkmul(&a[i * n + k], &b[j * n + k], &c[i * n + j]);
        }
    }

    pthread_exit(0);
}

/*****************************************************************/

/* Multiply (block)submatrices */
void blkmul(double *ablk, double *bblk, double *cblk)
{
    int i, j, k; /* Guess what... again... */

    for (i = 0; i < bs; i++) {
        for (j = 0; j < bs; j++) {
            for (k = 0; k < bs; k++) {
                cblk[i * n + j] +=
                    ablk[i * n + k] * bblk[j * n + k];
            }
        }
    }
}
