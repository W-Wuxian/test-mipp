#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <string.h>
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
    #pragma acc parallel loop default(present)
    for (int i=0; i<ni; i++){
        #pragma acc loop
        for (int j=0; j<nj; j++){
            for (int k=0; k<nk; k++){
                d[i*nj +j] = d[i*nj +j] + a[k+i*nk] * b[j+k*nj];
            }
            d[j+i*nj]= d[j+i*nj] + c[j+i*nj];
        }
    }
}

void naive_matmul_acc_tiled(int ni, int nj, int nk, double* a, double* b, double* c, double* d){
    #pragma acc parallel loop tile(32,32) default(present)
    for (int i=0; i<ni; i++){
        for (int j=0; j<nj; j++){
            for (int k=0; k<nk; k++){
                d[i*nj +j] = d[i*nj +j] + a[k+i*nk] * b[j+k*nj];
            }
        d[j+i*nj]= d[j+i*nj] + c[j+i*nj];
        }
    }
}

void seq_tiled_matmul(int tile, int ni, int nj, int nk, REAL* a, REAL* b, REAL* c, REAL* d){
    #pragma acc parallel loop default(present) num_workers(8) vector_length(128)
    #pragma acc loop gang collapse(2)
    for (int i=0; i<ni; i+=tile){
        for (int j=0; j<nj; j+=tile){
            #pragma acc loop worker
            for (int ii=i; ii< MIN(i+tile,ni); ii++){
                #pragma acc loop vector
                for (int jj=j; jj<MIN(j+tile,nj); jj++){
                    #pragma acc loop seq
                    for (int k=0; k<nk; k++){
                        d[ii*nj +jj] = d[ii*nj +jj] + a[k+ii*nk] * b[jj+k*nj];
                    }
                }
            }
        }
    }
    #pragma acc parallel loop default(present)
    for (int i=0; i<ni; i++){
        #pragma acc loop
        for (int j=0; j<nj; j++){
            d[j+i*nj]= d[j+i*nj] + c[j+i*nj];
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
    char* optim = "oacc";
    int ni = atoi(argv[1]);
    int nj = atoi(argv[2]);
    int nk = atoi(argv[3]);
    int tile = 1;
    int sample = atoi(argv[4]);
    int device = 1;
    clock_t  start, end;
    double cpu_time_used = 0;
    double Mem = ( ( (double)(ni*nk) + (double)(nk*nj) + (double)(ni*nj) + (double)(ni*nj) )  * sizeof(REAL) ) / (1024.0 * 1024.0); // MB
    REAL test = 0.;

    REAL* A = calloc(ni * nk, sizeof(REAL));
    REAL* B = calloc(nk * nj, sizeof(REAL));
    REAL* C = calloc(ni * nj, sizeof(REAL));
    REAL* D = calloc(ni * nj, sizeof(REAL));

    for ( int i = 0; i < ni*nk; i++)A[i]=(REAL)i/(ni*nk);//((REAL)i * GENRAND()) / (REAL)(ni*nj);
    for ( int i = 0; i < nk*nj; i++)B[i]=(REAL)i/(nk*nj);//((REAL)i * GENRAND()) / (REAL)(nj*nk);
    for ( int i = 0; i < ni*nj; i++)C[i]=(REAL)2;

    #pragma acc data copyin(A[0:ni*nk], B[0:nk*nj], C[0:ni*nj]) create(D[0:ni*nj])
    {
    for (int i = 0; i < sample; i++){
        start = clock();
        sequential_matmul(ni, nj, nk, A, B, C, D);
        end = clock();
        #pragma acc update self(D[0:ni*nj])
        cpu_time_used += ((double) (end-start)) / CLOCKS_PER_SEC;
        test = checksum(ni, nj, D);
        printf("MAX_INT %d iter %d checksum %f Time %lf TotalTime %lf\n", INT_MAX, i, test, cpu_time_used/(i+1), cpu_time_used);
        memset(D, 0, sizeof(REAL) * ni * nk);
        #pragma acc update device(D[0:ni*nj])
    }
    }
    

    FILE* fptr =  fopen("MATMUL_BENCH_OACC.csv", "a+");
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