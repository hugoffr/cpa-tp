#include <stdio.h>
#include <stdlib.h>

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
