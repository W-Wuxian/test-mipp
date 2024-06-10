#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <mpi.h>
// A is m x n
// B is n x p
// C is m x p

#define xstr(s) str(s)
#define str(s) #s

#ifndef REAL
#define REAL double
#endif

#ifndef LB
#define LB -1000
#endif

#ifndef UB
#define UB 1000
#endif

#if defined(LB) && defined(UB)
#define GENRAND() \
( rand()% (UB - LB +1) \
+ LB \
)
#endif

void sequential_matmul(int m, int n, int p, double* A, double* B, double* C ) {
    int i = 0, j = 0, k = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            for (k = 0; k < n; k++) {
                /*
                if ((i * p + j) >= (m*p))printf("CSUP\n");
                if ((i * n + k) >= (m*n))printf("ASUP\n");
                if ((k * p + j) >= (n*p))printf("BSUP\n");
                */
                C[i * p + j] += A[i * n + k] * B[k * p +j];
            }
        }
    }
}

void axis_workload(int rank, int Np, int dim1, int* s1, int* s2) {
    double CO_RE = 0.;
    int m = 0;
    if (Np == 1) {
        *s1 = 0;
        *s2 = dim1 - 1;
    } else {
        CO_RE = dim1 / Np;
        m = dim1 % Np;
        if (rank < m) {
            //printf("rank %d m %d\n", rank, m);
            *s1 = rank * (CO_RE + 1);
            *s2 = (rank + 1) * (CO_RE + 1) - 1;
        } else {
            *s1 = m + rank * CO_RE;
            *s2 = *s1 + CO_RE - 1;
        }
    }
}



int main(int argc, char** argv) {

    
    MPI_Init(&argc, &argv);
    int rank = 0; int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 5) {
        printf("Usage: %s <m> <n> <p> <s>\n", argv[0]);
        return 1;
    }

    srand(1234567);
    char* optim = "mpi";
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int p = atoi(argv[3]);
    int sample = atoi(argv[4]);
    int device = size;
    clock_t  start, end;
    double cpu_time_used = 0;
    double Mem = ( ( (double)(m*n) + (double)(n*p) + (double)(m*p) )  * sizeof(REAL) ) / (1024.0 * 1024.0); // MB

    int m_lc = 0, n_lc = 0, p_lc = 0;
    int lenA_lc = 0, lenB_lc = 0;
    int Am1 = 0, AmN = 0;
    int Bp1 = 0, BpN = 0;
 
    axis_workload(rank, size, m, &Am1, &AmN);
    axis_workload(rank, size, p, &Bp1, &BpN);
    lenA_lc = (AmN-Am1+1) * n;
    lenB_lc = n * (BpN-Bp1+1);
    m_lc = (AmN-Am1+1);
    p_lc = (BpN-Bp1+1);
    n_lc = n + (0 *n_lc);

    REAL* A = calloc((AmN-Am1+1) * n, sizeof(REAL));
    REAL* B = calloc(n * (BpN-Bp1+1), sizeof(REAL));
    REAL* C = calloc((AmN-Am1+1) * (BpN-Bp1+1), sizeof(REAL));

    for ( int i = 0; i < lenA_lc; i++)A[i]=((REAL)i * GENRAND()) / (REAL)(m*n);
    for ( int i = 0; i < lenB_lc; i++)B[i]=((REAL)i * GENRAND()) / (REAL)(n*p);

    for (int i = 0; i < sample; i++){
        start = MPI_Wtime();
        sequential_matmul(m_lc, n, p_lc, A, B, C);
        end = MPI_Wtime();
        cpu_time_used += (end-start);
        memset(C, 0, sizeof(double) * m_lc * p_lc);
    }

    if (rank == 0){
        FILE* fptr =  fopen("MATMUL_BENCH_MPI.csv", "a+");
        if (fptr == NULL) {
            free(A);
            free(B);
            free(C);
            printf("The file is not opened. Program will exit.\n");
            exit(0);
        }
        //header = "type,Total,Time,NSample,Mem,optim,device";
        fprintf(fptr, "%s,%lf,%lf,%d,%lf,%s,%d\n", xstr(REAL),cpu_time_used,cpu_time_used/sample,sample,Mem,optim,device);
        fclose(fptr);     
    }
    
    

    if (m<10 && p<10 && rank==0) {
        printf("Result:\n");
        for (int i = 0; i < m_lc; i++) {
            for (int j = 0; j < p_lc; j++) {
                printf("%f ", C[i * p_lc + j]);
            }
            printf("\n");
        }
    }
    

    if (A != NULL) {
        free(A);
        A = NULL;
    }
    if (B != NULL) {
        free(B);
        B = NULL;
    }
    if (C != NULL) {
        free(C);
        C = NULL;
    }

    MPI_Finalize();
    return 0;
}