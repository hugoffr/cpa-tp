#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

    clock_t tic, toc;

    double** A = (double**) malloc(n * sizeof(double*));
    double* b = (double*) malloc(n * sizeof(double));
    double* x = (double*) malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        A[i] = (double*) malloc(n * sizeof(double));
    }

    init_matrix_and_vector(A, b, n);

    tic = clock();
    
    // Perform LU factorization
    lu_factorization(A, n);
    
    // Solve for x
    solve_lu(A, b, x, n);

    toc = clock();

    printf("Sequential time for %d: %3.6f\n", n,(double)(toc - tic) / CLOCKS_PER_SEC);
}

int main() {

    for (int n = 1024; n <= 8192; n = n + 1024) {
        sequential_solution(n);
    }
}
=======

void LU_factorization(double** A, double** L, double** U, int n) {
    int i, j, k;
    double sum;

    // Initialize L and U matrices
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            U[i][j] = A[i][j];
            if (i == j) {
                L[i][j] = 1.0;
            } else {
                L[i][j] = 0.0;
            }
        }
    }

    // Perform Gaussian elimination to obtain L and U matrices
    for (k = 0; k < n-1; k++) {
        for (i = k+1; i < n; i++) {
            if (U[k][k] == 0) {
                printf("Error: Division by zero\n");
                exit(1);
            }
            L[i][k] = U[i][k] / U[k][k];
            for (j = k; j < n; j++) {
                U[i][j] = U[i][j] - L[i][k] * U[k][j];
            }
        }
    }
}

int main() {
    int n = 3;
    double** A = (double**) malloc(n * sizeof(double*));
    double** L = (double**) malloc(n * sizeof(double*));
    double** U = (double**) malloc(n * sizeof(double*));

    for (int i = 0; i < n; i++) {
        A[i] = (double*) malloc(n * sizeof(double));
        L[i] = (double*) malloc(n * sizeof(double));
        U[i] = (double*) malloc(n * sizeof(double));
    }

    // Initialize matrix A
    A[0][0] = 4; A[0][1] = 3; A[0][2] = 2;
    A[1][0] = 2; A[1][1] = 1; A[1][2] = 3;
    A[2][0] = 1; A[2][1] = 1; A[2][2] = 1;

    LU_factorization(A, L, U, n);

    // Print the matrices L and U
    printf("L:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", L[i][j]);
        }
        printf("\n");
    }

    printf("U:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", U[i][j]);
        }
        printf("\n");
    }

    // Free memory
    for (int i = 0; i < n; i++) {
        free(A[i]);
        free(L[i]);
        free(U[i]);
    }
    free(A);
    free(L);
    free(U);

    return 0;
}
