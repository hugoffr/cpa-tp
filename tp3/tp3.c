#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <papi.h>

void lu_factorization(double** A, int n) {
    int i, j, k;
    for (k = 0; k < n-1; k++) {
        for (i = k+1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (j = k+1; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            A[i][k] = factor;
        }
    }
}

void solve_lu(double** A, double* b, double* x, int n) {
    int i, j;
    double sum;
    double* y = (double*) malloc(n * sizeof(double));
    
    // Forward substitution
    y[0] = b[0];
    for (i = 1; i < n; i++) {
        sum = 0.0;
        for (j = 0; j < i; j++) {
            sum += A[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }
    
    // Backward substitution
    x[n-1] = y[n-1] / A[n-1][n-1];
    for (i = n-2; i >= 0; i--) {
        sum = 0.0;
        for (j = i+1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / A[i][i];
    }
    
    free(y);
}

void init_matrix_and_vector(double** A, double* b, int n) {
    int i, j;
    srand(time(NULL));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = (double) rand() / RAND_MAX;
        }
        b[i] = (double) rand() / RAND_MAX;
    }
}

// This function should only be used to test the validity of the code
void print_results(double** A, double* b, double* x, int n) {

    int i, j;

    // Print results
    printf("Matrix A:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }

    printf("Vector b:\n");
    for (i = 0; i < n; i++) {
        printf("%f ", b[i]);
    }

    printf("\nSolution x:\n");
    for (i = 0; i < n; i++) {
        printf("%f ", x[i]);
    }
    printf("\n");
}

void sequential_solution(int n) {

    //---------------------------------------------------------
    // Initialize Counters

    int ret, EventSet;

    long long values[4];

    clock_t tic, toc;

    ret = PAPI_library_init(PAPI_VER_CURRENT);
    if (ret != PAPI_VER_CURRENT) printf("FAIL\n");

    ret = PAPI_create_eventset(&EventSet);
    if (ret != PAPI_OK) printf("ERRO: create eventset\n");

    ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK) printf("ERRO: PAPI_L1_DCM\n");

    ret = PAPI_add_event(EventSet, PAPI_L2_DCM);
    if (ret != PAPI_OK) printf("ERRO: PAPI_L2_DCM\n");

    ret = PAPI_add_event(EventSet, PAPI_FP_INS);
    if (ret != PAPI_OK) printf("ERRO: PAPI_FP_INS\n");

    ret = PAPI_add_event(EventSet, PAPI_TOT_INS);
    if (ret != PAPI_OK) printf("ERRO: PAPI_TOT_INS\n");

    //
    //---------------------------------------------------------

    double** A = (double**) malloc(n * sizeof(double*));
    double* b = (double*) malloc(n * sizeof(double));
    double* x = (double*) malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        A[i] = (double*) malloc(n * sizeof(double));
    }

    init_matrix_and_vector(A, b, n);

    //---------------------------------------------------------
    // Start Counting
    tic = clock();

    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK) printf("ERRO: Start PAPI\n");
    //
    //---------------------------------------------------------
    
    // Perform LU factorization
    lu_factorization(A, n);
    
    // Solve for x
    solve_lu(A, b, x, n);

    //---------------------------------------------------------
    // Stop Counting
    toc = clock();

    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK) printf("ERRO: Stop PAPI\n");
    printf("L1 DCM: %lld \n",values[0]);
    printf("L2 DCM: %lld \n",values[1]);
    printf("Flops: %lld \n",values[2]);
    printf("Total Instructions: %lld \n",values[3]);

    ret = PAPI_reset( EventSet );
    if ( ret != PAPI_OK ) printf("FAIL reset\n"); 

    printf("Sequential time for %d: %3.6f\n", n,(double)(toc - tic) / CLOCKS_PER_SEC);

    ret = PAPI_remove_event( EventSet, PAPI_L1_DCM );
    if ( ret != PAPI_OK ) printf("FAIL remove event\n"); 

    ret = PAPI_remove_event( EventSet, PAPI_L2_DCM );
    if ( ret != PAPI_OK ) printf( "FAIL remove event\n");

    ret = PAPI_remove_event( EventSet, PAPI_FP_INS );
    if ( ret != PAPI_OK ) printf( "FAIL remove event\n"); 

    ret = PAPI_remove_event( EventSet, PAPI_TOT_INS );
    if ( ret != PAPI_OK ) printf( "FAIL remove event\n"); 

    ret = PAPI_destroy_eventset( &EventSet );
    if ( ret != PAPI_OK ) printf("FAIL destroy\n");

    //
    //---------------------------------------------------------

}

int main() {

    for (int n = 1024; n <= 8192; n = n + 1024) {
        sequential_solution(n);
    }
}
