#include <stdio.h>
#include <stdlib.h>
#include <time.h>



void OnMultLine(double** A, double** B)
{
	char st[100];
	double temp;
	int i, j, k;
	for(i=0; i<m_ar; i++)
	{	for( j=0; j<m_br; j++)
		{
			for( k=0; k<m_ar; k++)
			{	
				phc[i*m_ar+k] += pha[i*m_ar+k] * phb[j*m_br+k];
			}
		}
	}
}


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

void block_lu_factorization(double** A, int n, int block_size){

    //TODO Check for n non multiple of block_size

    if (block_size > n){
        lu_factorization(A, n);
        return;
    }
    for (int pivot_block = 0; pivot_block < n/block_size+1; pivot_block++){
        //We have to do this for each pivot block
        for (int pivot_i = 0; pivot_i < block_size; pivot_i++){
            double pivot = A[pivot_block*block_size+pivot_i][pivot_block*block_size+pivot_i];

            //Solve pivot block
            //Divide left column by pivot
            for (int left_column_i = pivot_block*block_size+pivot_i; left_column_i < block_size; left_column_i++){
                A[pivot_block * block_size + pivot_i][left_column_i] /= pivot;
            }
            //Update all remainig values
            for (int inner_mat_x = pivot_block*block_size+pivot_i+1; inner_mat_x < (pivot_block+1) * block_size; inner_mat_x++){
                for (int inner_mat_y = pivot_block*block_size+pivot_i+1; inner_mat_y < (v_block_i+1) * block_size; inner_mat_y++){
                    A[inner_mat_x][inner_mat_y] -= A[inner_mat_x][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_y];
                }
            }


            //Check if it is the last block

            //For each vertical block under the current pivot block
            for(int v_block_i = pivot_block + 1 ; v_block_i < n/block_size+1;v_block_i ++){
                //Divide left column by pivot
                for (int left_column_i = 0; left_column_i < block_size; left_column_i++){
                    A[pivot_block*block_size+pivot_i][pivot_block*block_size+left_column_i] /= pivot;
                }
                //Update all remainig values
                for (int inner_mat_x = pivot_block*block_size+pivot_i+1; inner_mat_x < (pivot_block+1) * block_size; inner_mat_x++){
                    for (int inner_mat_y = v_block_i*block_size; inner_mat_y < (v_block_i+1) * block_size; inner_mat_y++){
                        A[inner_mat_x][inner_mat_y] -= A[inner_mat_x][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_y];
                    }
                }
            }

            //For each horizontal block to the right of the pivot block

            for(int h_block_i = pivot_block + 1 ; h_block_i < n/block_size+1;h_block_i ++){
                //Update all remainig values
                for (int inner_mat_x = h_block_i*block_size; inner_mat_x < (h_block_i+1) * block_size; inner_mat_x++){
                    for (int inner_mat_y = pivot_block*block_size+pivot_i+1; inner_mat_y < (pivot_block+1) * block_size; inner_mat_y++){
                        A[inner_mat_x][inner_mat_y] -= A[inner_mat_x][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_y];
                    }
                }
            }

            //Set each inner block as the multiplication of its upper most and left most blocks
            for(int h_block_i = pivot_block+1;h_block_i < n/block_size+1;h_block_i ++){
                for(int v_block_i = pivot_block+1;v_block_i < n/block_size+1;v_block_i ++){

                    //Start by zeroing the inner block
                    for(int block_x = h_block_i * block_size; block_x < (h_block_i+1)*block_size; block_x++){
                        for(int block_y = h_block_i * block_size; block_y < (h_block_i+1)*block_size; block_y++){
                            A[block_x][block_y] = 0;
                        }
                    }

                    //Multiply the matrices by line

                    for(int line_i = 0; line_i < block_size; line_i++){
                        for (int left_mat_i = 0; left_mat_i < block_size; left_mat_i++){
                            for (int top_mat_i = 0; top_mat_i < block_size; top_mat_i++){
                                A[pivot_block*block_size+left_mat_i][pivot_block*block_size+line_i + left_mat_i] += A[pivot_block*block_size+left_mat_i][block_y+line_i] * A[block_x+top_mat_i][pivot_block*block_size+line_i+left_mat_i];
                            }
                        }
                    }
                }
            }

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
