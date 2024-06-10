#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#define RUN dynamic
// dynamic static guided
#endif // _OPENMP
// A is ni x nk
// B is nk x nj
// C is ni x nj
// D is ni x nj

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

#define MIN(a,b) ( ((a)<(b))?(a):(b) )

#if defined(LB) && defined(UB)
#define GENRAND() \
( rand()% (UB - LB +1) \
+ LB \
)
#endif

void sequential_matmul(int ni, int nj, int nk, REAL* a, REAL* b, REAL* c, REAL* d){
    #pragma omp parallel
    {
    #pragma omp for schedule(RUN)
    for (int i=0; i<ni; i++){
        for (int j=0; j<nj; j++){
            for (int k=0; k<nk; k++){
                d[i*nj +j] = d[i*nj +j] + a[k+i*nk] * b[j+k*nj];
            }
            d[j+i*nj]= d[j+i*nj] + c[j+i*nj];
        }
    }
    }
    
}

REAL checksum(int ni, int nj, REAL* d) {
    REAL dsum = 0.;
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            dsum = dsum + d[j + i * nj];
        }
    }
    return dsum;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        printf("Usage: %s <ni> <nj> <nk> <s>\n", argv[0]);
        return 1;
    }

    //srand(1234567);
    char* optim = "omp";
    int ni = atoi(argv[1]);
    int nj = atoi(argv[2]);
    int nk = atoi(argv[3]);
    int tile = 1;
    int sample = atoi(argv[4]);
    int device = 1;

    double start = 0., end = 0.;
    struct timeval t_elapsed_0, t_elapsed_1;
    double t_elapsed = 0.;
    clock_t t_cpu_0, t_cpu_1;
    double  t_cpu = 0.;
    double cpu_time_used = 0;

    double Mem = ( ( (double)(ni*nk) + (double)(nk*nj) + (double)(ni*nj) + (double)(ni*nj) )  * sizeof(REAL) ) / (1024.0 * 1024.0); // MB
    REAL test = 0.;
    #ifdef _OPENMP
    int nb_taches;
    #pragma omp parallel
    { nb_taches = omp_get_num_threads(); }
    fprintf(stdout, "\n\n   Execution prod_mat en parallele avec %d threads\n", nb_taches);
    device = nb_taches;
    #endif // _OPENMP


    REAL* A = calloc(ni * nk, sizeof(REAL));
    REAL* B = calloc(nk * nj, sizeof(REAL));
    REAL* C = calloc(ni * nj, sizeof(REAL));
    REAL* D = calloc(ni * nj, sizeof(REAL));

    #pragma omp parallel
    {
    #pragma omp for schedule(RUN) nowait
    for ( int i = 0; i < ni*nk; i++){
        A[i]=(REAL)i/(ni*nk);//((REAL)i * GENRAND()) / (REAL)(ni*nj);
    }
    #pragma omp for schedule(RUN) nowait
    for ( int i = 0; i < nk*nj; i++){
        B[i]=(REAL)i/(nk*nj);//((REAL)i * GENRAND()) / (REAL)(nj*nk);
    }
    #pragma omp for schedule(RUN) nowait
    for ( int i = 0; i < ni*nj; i++){
        C[i]=(REAL)2;
    }
    #pragma omp for schedule(RUN)
    for ( int i = 0; i < ni*nj; i++){
        D[i]=(REAL)0;
    }
    }

    

    for (int i = 0; i < sample; i++){
        // Temps CPU de calcul initial.
        //t_cpu_0 = clock();
        // Temps elapsed de reference.
        //gettimeofday(&t_elapsed_0, NULL);
        start = omp_get_wtime();
        sequential_matmul(ni, nj, nk, A, B, C, D);
        end = omp_get_wtime();
        // Temps elapsed final
        //gettimeofday(&t_elapsed_1, NULL);
        //t_elapsed += (t_elapsed_1.tv_sec - t_elapsed_0.tv_sec) + (t_elapsed_1.tv_usec - t_elapsed_0.tv_usec) / (double)1000000;
        // Temps CPU de calcul final
        //t_cpu_1 = clock();
        //t_cpu += (t_cpu_1 - t_cpu_0) / (double)CLOCKS_PER_SEC;
        cpu_time_used += (end - start);
        test = checksum(ni, nj, D);
        //printf("MAX_INT %d iter %d checksum %f Time %lf TotalTime %lf elapse %lf TotalElapse %lf t_cpu %lf Totalt_cpu %lf\n", INT_MAX, i, test, cpu_time_used/(i+1), cpu_time_used, t_elapsed/(i+1), t_elapsed, t_cpu/(i+1), t_cpu);
        printf("MAX_INT %d iter %d checksum %f Time %lf TotalTime %lf\n", INT_MAX, i, test, cpu_time_used/(i+1), cpu_time_used);
        memset(D, 0, sizeof(REAL) * ni * nk);
    }

    FILE* fptr =  fopen("MATMUL_BENCH_OMP.csv", "a+");
    if (fptr == NULL) {
        free(A);
        A = NULL;
        free(B);
        B = NULL;
        free(C);
        C = NULL;
        free(D);
        D = NULL;
        printf("The file is not opened. Program will exit.\n");
        exit(0);
    }
    //header = "type,Total,Time,NSample,Mem,optim,tile,device,test";
    fprintf(fptr, "%s,%lf,%lf,%d,%lf,%s,%d,%d,%f\n", xstr(REAL),cpu_time_used,cpu_time_used/sample,sample,Mem,optim,tile,device,test);
    fclose(fptr);
    

    if (ni<10 && nk<10) {
        printf("Result:\n");
        for (int i = 0; i < ni; i++) {
            for (int j = 0; j < nk; j++) {
                printf("%f ", D[i * nk + j]);
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
    if (D != NULL) {
        free(D);
        D = NULL;
    }

    return 0;
}