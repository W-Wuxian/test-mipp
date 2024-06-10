#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <openacc.h>
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

void openacc_matmul(int m, int n, int p, double* A, double* B, double* C) {
    #pragma acc parallel loop
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                C[i * p + j] += A[i * n + k] * B[k * p + j];
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argc != 5) {
        printf("Usage: %s <m> <n> <p> <s>\n", argv[0]);
        return 1;
    }

    srand(1234567);
    char* optim = "openacc";
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int p = atoi(argv[3]);
    int sample = atoi(argv[4]);
    int device = 1;
    
    double gpu_time_used = 0;
    double Mem = ( ( (double)(m*n) + (double)(n*p) + (double)(m*p) )  * sizeof(REAL) ) / (1024.0 * 1024.0); // MB

    REAL* A = calloc(m * n, sizeof(REAL));
    REAL* B = calloc(n * p, sizeof(REAL));
    REAL* C = calloc(m * p, sizeof(REAL));

    for ( int i = 0; i < m*n; i++)A[i]=((REAL)i * GENRAND()) / (REAL)(m*n);
    for ( int i = 0; i < n*p; i++)B[i]=((REAL)i * GENRAND()) / (REAL)(n*p);

    acc_device_t device_type = acc_get_device_type();
    int num_devices = acc_get_num_devices(device_type);
    const char* device_name = acc_get_device_name(device_type, 0);

    acc_event_t start, end;

    acc_init(acc_device_nvidia);
    acc_set_device_num(num_devices, acc_device_nvidia);

    for (int i = 0; i < sample; i++){
        acc_get_start_event(&start);
        openacc_matrix_multiply(m, n, p, A, B, C);
        acc_get_end_event(&end);
        acc_wait_all_async();
        gpu_time_used += acc_get_elapsed_time(start, end) / 1000.0;
        memset(C, 0, sizeof(double) * m * p);
    }

    FILE* fptr =  fopen("MATMUL_BENCH_OPENACC.csv", "a+");
    if (fptr == NULL) {
        free(A);
        free(B);
        free(C);
        printf("The file is not opened. Program will exit.\n");
        exit(0);
    }
    //header = "type,Total,Time,NSample,Mem,optim,device";
    fprintf(fptr, "%s,%lf,%lf,%d,%lf,%s,%d\n", xstr(REAL),gpu_time_used,gpu_time_used/sample,sample,Mem,optim,device);
    fclose(fptr);
    

    if (m<10 && p<10) {
        printf("Result:\n");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < p; j++) {
                printf("%f ", C[i * p + j]);
            }
            printf("\n");
        }
    }

    printf("Device: %s\n", device_name);
    printf("Time taken: %f seconds\n", gpu_time_used);
    

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

    acc_shutdown(acc_device_nvidia);

    return 0;
}