#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/*
void OnMultLine(double** A, double** B, int n)
{
	char st[100];
	double temp;
	int i, j, k;
	for(i=0; i<n; i++)
	{	for( j=0; j<n; j++)
		{
			for( k=0; k<n; k++)
			{	
				phc[i*m_ar+k] += pha[i*m_ar+k] * phb[j*m_br+k];
			}
		}
	}
}*/


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


void printMatrix(double** m, int n){
    printf("Matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", m[i][j]);
        }
        printf("\n");
    }
}

void block_lu_factorization(double** A, int n, int block_size){

    if (block_size > n){
        lu_factorization(A, n);
        return;
    }
    for (int pivot_block = 0; pivot_block < n/block_size; pivot_block++){
        //printf("Starting pivot block %d\n",pivot_block);
        //printMatrix(A,n);
        //We have to do this for each pivot block
        for (int pivot_i = 0; pivot_i < block_size; pivot_i++){
            double pivot = A[pivot_block*block_size+pivot_i][pivot_block*block_size+pivot_i];
            //printf("Working on pivot %d: %lf\n",pivot_i,pivot);

            //Solve pivot block
            //Divide left column by pivot
            for (int left_column_i = pivot_block*block_size+pivot_i; left_column_i < (pivot_block+1)*block_size; left_column_i++){
                A[left_column_i][pivot_block * block_size + pivot_i] /= pivot;
                //printf("Dividing element at (%d,%d) by pivot %f\n",left_column_i,pivot_block*block_size+pivot_i,pivot);
            }
            //printMatrix(A,n);
            //printf("Divided left line by pivot\n");
            //Update all remainig values
            for (int inner_mat_x = pivot_block*block_size+pivot_i+1; inner_mat_x < (pivot_block+1) * block_size; inner_mat_x++){
                for (int inner_mat_y = pivot_block*block_size+pivot_i+1; inner_mat_y < (pivot_block+1) * block_size; inner_mat_y++){
                    A[inner_mat_x][inner_mat_y] -= A[inner_mat_x][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_y];
                }
            }
            //printMatrix(A,n);
            //printf("Updated pivot block inner matrix\n");


            //Check if it is the last block

            //For each vertical block under the current pivot block
            //printf("Updating vertical blocks\n");
            for(int v_block_i = pivot_block + 1 ; v_block_i < n/block_size;v_block_i ++){
                //printf("Solving vertical block %d\n",v_block_i);
                //Divide left column by pivot
                for (int left_column_i = 0; left_column_i < block_size; left_column_i++){
                    A[v_block_i*block_size+left_column_i][pivot_block*block_size+pivot_i] /= pivot;
                    //printf("Divinding element at %d,%d (%f) by %f:  pivot_block=%d  v_block_i=%d  left_column_i=%d\n",pivot_block*block_size+pivot_i,v_block_i*block_size+left_column_i,A[pivot_block*block_size+pivot_i][pivot_block*block_size+left_column_i],pivot,pivot_block,v_block_i,left_column_i);
                }
                //Update all remainig values
                for (int inner_mat_x = pivot_block*block_size+pivot_i+1; inner_mat_x < (pivot_block+1) * block_size; inner_mat_x++){
                    for (int inner_mat_y = v_block_i*block_size; inner_mat_y < (v_block_i+1) * block_size; inner_mat_y++){
                        //printf("Solving element at %d,%d\n",inner_mat_x,inner_mat_y);
                        A[inner_mat_y][inner_mat_x] -= A[inner_mat_y][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_x];
                    }
                }
                //printMatrix(A,n);
            }

            //For each horizontal block to the right of the pivot block

            //printf("Updating horizontal blocks\n");
            for(int h_block_i = pivot_block + 1 ; h_block_i < n/block_size;h_block_i ++){
                //printf("Solving horizontal block %d\n",h_block_i);
                //Update all remainig values
                for (int inner_mat_x = (h_block_i+pivot_block)*block_size; inner_mat_x < (pivot_block + h_block_i+1) * block_size; inner_mat_x++){
                    for (int inner_mat_y = pivot_block*block_size+pivot_i+1; inner_mat_y < (pivot_block+1) * block_size; inner_mat_y++){
                        //printf("Analysing %d,%d\n",inner_mat_x,inner_mat_y);
                        A[inner_mat_y][inner_mat_x] -= A[inner_mat_y][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_x];
                    }
                }
                //printMatrix(A,n);
            }

            //printf("Updating the rest of the matrix\n");
            //printMatrix(A,n);
            //Set each inner block as the multiplication of its upper most and left most blocks
            for(int h_block_i = pivot_block+1;h_block_i < n/block_size;h_block_i ++){
                for(int v_block_i = pivot_block+1;v_block_i < n/block_size;v_block_i ++){
                    //printf("Solving block %d,%d\n",h_block_i,v_block_i);
                    for (int inner_mat_x = (h_block_i)*block_size; inner_mat_x < (h_block_i+1) * block_size; inner_mat_x++){
                        for (int inner_mat_y = v_block_i*block_size; inner_mat_y < (v_block_i+1) * block_size; inner_mat_y++){
                            A[inner_mat_y][inner_mat_x] -= A[inner_mat_y][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_x];
                        }
                    }

                    /*
                    //Start by zeroing the inner block
                    for(int block_x = h_block_i * block_size; block_x < (h_block_i+1)*block_size; block_x++){
                        for(int block_y = h_block_i * block_size; block_y < (h_block_i+1)*block_size; block_y++){
                            A[block_y][block_x] = 0;
                        }
                    }

                    //Multiply the matrices by line
                    for(int line_i = 0; line_i < block_size; line_i++){
                        for (int left_mat_i = 0; left_mat_i < block_size; left_mat_i++){
                            for (int top_mat_i = 0; top_mat_i < block_size; top_mat_i++){
                                A[pivot_block*block_size+line_i + left_mat_i][pivot_block*block_size+left_mat_i] += A[pivot_block*block_size+left_mat_i][v_block_i*block_size+line_i] * A[h_block_i*block_size+top_mat_i][pivot_block*block_size+line_i+left_mat_i];
                            }
                        }
                    }

                    for(int left_matrix_line = 0; left_matrix_line < block_size; left_matrix_line++){
                        for(int top_matrix_line = 0; top_matrix_line < block_size; top_matrix_line++){
                            for(int column_i = 0; column_i < block_size; column_i++){
                                A[v_block_i*block_size + left_matrix_line][h_block_i*block_size + column_i] += A[v_block_i*block_size+left_matrix_line][pivot_block* block_size + top_matrix_line]
                            }
                        }
                    }
                    */
                }
            }
            //printf("Updated Matrix:\n");
            //printMatrix(A,n);
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

    /*for (int n = 1024; n <= 8192; n = n + 1024) {
        sequential_solution(n);
    }*/

    int n=9;
    double** A = (double**) malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++){
        A[i] = (double*)malloc(n*sizeof(double));
        for(int j = 0; j < n; j++){
            A[i][j] = (double)(rand()%100-50);
            //A[i][j] = (double)(1+2*i+j);
        }
    }

    block_lu_factorization(A, n, 3);

    printMatrix(A,n);
}
