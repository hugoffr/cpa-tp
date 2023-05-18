#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>    
#include <mpi.h>
#include <math.h>


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

void distributed_block_lu_factorization(double** A, int n, int block_size){



    double t1,t2,t3,t4,t5,t6,t7,t8,t9,st,et;
    int c1,c2,c3,c4,c5,c6,c7,c8,c9;

    int verbose = 0;
    // if (n == 512){
    //     verbose = 1;
    // }

    if (block_size > n){
        lu_factorization(A, n);
        return;
    }
    //For each pivot block

    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        printf("[%d]Process %d: Now starting solution of matrix with n=%d (block_size=%d)\n",n,rank,n,block_size);
    }

    for (int pivot_block = 0; pivot_block < n/block_size; pivot_block++){
        //For each pivot inside a pivot block
        printf("\n[%d]Process %d: Solving pivot block %d\n",n,rank,pivot_block);
        for (int pivot_i = 0; pivot_i < block_size; pivot_i++){
            st = omp_get_wtime();
            double pivot = A[pivot_block*block_size+pivot_i][pivot_block*block_size+pivot_i];
            //Divide left column by pivot
            for (int left_column_i = pivot_block*block_size+pivot_i; left_column_i < (pivot_block+1)*block_size; left_column_i++){
                A[left_column_i][pivot_block * block_size + pivot_i] /= pivot;
            }
            //Update all remainig values
            for (int inner_mat_x = pivot_block*block_size+pivot_i+1; inner_mat_x < (pivot_block+1) * block_size; inner_mat_x++){
                for (int inner_mat_y = pivot_block*block_size+pivot_i+1; inner_mat_y < (pivot_block+1) * block_size; inner_mat_y++){
                    A[inner_mat_x][inner_mat_y] -= A[inner_mat_x][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_y];
                }
            }
            et = omp_get_wtime();
            t1 += (et-st);
            c1++;
            st = omp_get_wtime();
            //For each vertical block under the current pivot block
            int process_block_number = ceil(((n/block_size) - (pivot_block + 1)) / (double)numProcesses);
            //printf("DEBUG: n/block_size=%d pivot_block+1=%d numProcesses=%d Process_block_number=%d\n",n/block_size,pivot_block + 1,numProcesses,process_block_number);
            // Send data from process 0 to other processes
            if (rank == 0) {
                //printf("[%d]Process %d: Now sending vertical blocks\n",n,rank);
                for (int process = 1; process < numProcesses; process++) { //Send to process 'process'
                    //printf("[%d]Process %d: Sending to %d rows %d-%d\n",n,rank,process,(pivot_block+1+process*process_block_number)*block_size,(pivot_block+1+(process+1)*process_block_number)*block_size);
                    for (int row_i = (pivot_block+1+process*process_block_number)*block_size;row_i < (pivot_block+1+(process+1)*process_block_number)*block_size && row_i < n;row_i++){
                        //printf("[%d]Process %d: Sending row %d to process %d (tag %d)\n",n,rank,row_i,process,process);
                        MPI_Send(&A[row_i][pivot_block*block_size], block_size, MPI_DOUBLE, process, process, MPI_COMM_WORLD);
                    }
                }
            } else {
                // Receive data from process 0
                //printf("[%d]Process %d: Now receiving vertical blocks\n",n,rank);
                for (int row_i = (pivot_block+1+rank*process_block_number)*block_size;row_i < (pivot_block+1+(rank+1)*process_block_number)*block_size && row_i < n;row_i++) {
                    //printf("[%d]Process %d: Receiving row %d (tag%d)\n",n,rank,row_i,rank);
                    
                    MPI_Recv(&A[row_i][pivot_block*block_size], block_size, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            et = omp_get_wtime();
            t2 += (et-st);
            c2++;
            st = omp_get_wtime();


            //MPI_Barrier(MPI_COMM_WORLD);

            if(rank==0 && verbose){
                printf("[%d:%d]STATUS: Successfully transmitted vertical blocks\n",n,pivot_block);
            }


            int start_v = pivot_block + 1 + rank * process_block_number;
            int end_v = pivot_block + 1 + (rank + 1) * process_block_number;
            if (start_v * block_size < n){
                for(int v_block_i = start_v; v_block_i < end_v; v_block_i++){
                    if (v_block_i >= n/block_size){
                        break;
                    }
                    //Divide left column by pivot
                    if(verbose){
                        printf("[%d:%d]Process %d:Dividing Left Column\n",n,pivot_block,rank);
                    }
                    for (int left_column_i = 0; left_column_i < block_size; left_column_i++) {
                        if(verbose){
                            //printf("Analysing x=%d _y=%d\n",v_block_i*block_size+left_column_i,pivot_block*block_size+pivot_i);
                        }
                        A[v_block_i*block_size+left_column_i][pivot_block*block_size+pivot_i] /= pivot;
                    }
                    if(verbose){
                        printf("[%d:%d]Process %d: Updating matrix in h=%d-%d v=%d-%d\n\n",n,pivot_block,rank,pivot_block*block_size+pivot_i+1,(pivot_block+1) * block_size,v_block_i*block_size,(v_block_i+1) * block_size);
                    }
                    // Update all remaining values
                    for (int inner_mat_x = pivot_block*block_size+pivot_i+1; inner_mat_x < (pivot_block+1) * block_size; inner_mat_x++) {
                        for (int inner_mat_y = v_block_i*block_size; inner_mat_y < (v_block_i+1) * block_size; inner_mat_y++) {
                            if(inner_mat_x >= n || inner_mat_y >= n || pivot_block*block_size+pivot_i >=n){
                                printf("Error! inner_mat_x=%d inner_mat_y=%d const=%d\n",inner_mat_x,inner_mat_y,pivot_block*block_size+pivot_i);
                            }
                            A[inner_mat_y][inner_mat_x] -= A[inner_mat_y][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_x];
                        }
                    }

                    if(verbose){
                        printf("[%d:%d]Process %d: Updated rest of matrix\n",n,pivot_block,rank);
                    }
                }
            }

            
            et = omp_get_wtime();
            t3 += (et-st);
            c3++;
            st = omp_get_wtime();

            //MPI_Barrier(MPI_COMM_WORLD);

            if(rank==0 && verbose){
                printf("[%d:%d]STATUS: Successfully calculated vertical blocks\n",n,pivot_block);
            }


            if (rank == 0) {
                for (int process = 1; process < numProcesses; process++) { //Send to process 'process'
                    //printf("[%d]Process %d: Sending to %d rows %d-%d\n",n,rank,process,(pivot_block+1+process*process_block_number)*block_size,(pivot_block+1+(process+1)*process_block_number)*block_size);
                    for (int row_i = pivot_block*block_size;row_i < (pivot_block+1)*block_size && row_i<n;row_i++){
                        //printf("[%d]Process %d: Sending row %d to process %d (tag %d)\n",n,rank,row_i,process,process*100000+row_i);
                        MPI_Send(&A[row_i][(pivot_block+1+(process_block_number)*process)*block_size], block_size*process_block_number, MPI_DOUBLE, process, process*100000+row_i, MPI_COMM_WORLD);
                    }
                    //printf("[%d]Process %d: Done sending rows to process %d\n", n,rank,process);
                }
            } else {
                // Receive data from process 0
                //printf("[%d]Process %d: Starting to receive data(from %d to %d)\n",n,rank,pivot_block*block_size,(pivot_block+1)*block_size);
                for (int row_i = pivot_block*block_size;row_i < (pivot_block+1)*block_size && row_i<n;row_i++) {
                    //printf("[%d]Process %d: Waiting to receive row %d (tag %d)\n",n,rank,row_i,rank*100000+row_i);
                    MPI_Recv(&A[row_i][(pivot_block+1+(process_block_number)*rank)*block_size], block_size*process_block_number, MPI_DOUBLE, 0, rank*100000+row_i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //printf("[%d]Process %d: Received row %d\n",n,rank,row_i);
                }
            }


            et = omp_get_wtime();
            t4 += (et-st);
            c4++;
            st = omp_get_wtime();
            
            //MPI_Barrier(MPI_COMM_WORLD);

            if (rank == 0 && verbose){
                printf("[%d:%d]STATUS: Successfully calculated horizontal blocks\n",n,pivot_block);
            }


            //For each horizontal block to the right of the pivot block
            int start_h = pivot_block + 1 + rank * process_block_number;
            int end_h = pivot_block + 1 + (rank + 1) * process_block_number;

            // Parallelize the second loop for the horizontal blocks
            for (int h_block_i = start_h; h_block_i < end_h; h_block_i++) {
                // Update all remaining values
                for (int inner_mat_x = (h_block_i+pivot_block)*block_size; inner_mat_x < (pivot_block + h_block_i+1) * block_size; inner_mat_x++) {
                    for (int inner_mat_y = pivot_block*block_size+pivot_i+1; inner_mat_y < (pivot_block+1) * block_size; inner_mat_y++) {
                        A[inner_mat_y][inner_mat_x] -= A[inner_mat_y][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_x];
                    }
                }
            }


            //MPI_Barrier(MPI_COMM_WORLD);


            et = omp_get_wtime();
            t5 += (et-st);
            c5++;
            st = omp_get_wtime();

            if (rank == 0) {
                //Receive data from all processes
                for (int process = 1; process < numProcesses; process++) { //Send to process 'process'
                    for (int row_i = (pivot_block+1+process*process_block_number)*block_size;row_i < (pivot_block+1+(process+1)*process_block_number);row_i++){
                        MPI_Recv(&A[row_i][pivot_block*block_size], block_size, MPI_DOUBLE, process, process, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            } else {
                // Send data to process 0
                for (int row_i = (pivot_block+1+rank)*process_block_number;row_i < (pivot_block+1+rank+1)*process_block_number;row_i++) {
                        MPI_Send(&A[row_i][pivot_block*block_size], block_size, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
                }
            }
            if (rank == 0) {
                //Receive data from all processes
                for (int process = 1; process < numProcesses; process++) { //Send to process 'process'
                    for (int row_i = pivot_block*block_size;row_i < block_size;row_i++){
                        MPI_Recv(&A[row_i][(pivot_block+1+(process_block_number)*rank)*block_size], block_size*process_block_number, MPI_DOUBLE, process, process, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            } else {
                // Send data to process 0
                for (int row_i = pivot_block*block_size;row_i < block_size;row_i++) {
                    MPI_Send(&A[row_i][(pivot_block+1+(process_block_number)*rank)*block_size], block_size*process_block_number, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
                }
            }


            et = omp_get_wtime();
            t6 += (et-st);
            c6++;
            st = omp_get_wtime();
            


            if (rank == 0 && verbose){
                printf("[%d:%d]STATUS: Successfully transmited edge blocks results\n",n,pivot_block);
            }

            MPI_Barrier(MPI_COMM_WORLD);



            if (rank == 0) {
                for (int process = 1; process < numProcesses; process++) { //Send to process 'process'
                    for (int row_i = (pivot_block+1+process*process_block_number)*block_size;row_i < (pivot_block+1+(process+1)*process_block_number)*block_size && row_i<n;row_i++){
                        MPI_Send(&A[row_i][(pivot_block+1)*block_size], n/block_size - pivot_block - 1, MPI_DOUBLE, process, process, MPI_COMM_WORLD);
                    }
                }
            } else {
                // Receive data from process 0
                    for (int row_i = (pivot_block+1+rank*process_block_number)*block_size;row_i < (pivot_block+1+(rank+1)*process_block_number)*block_size && row_i<n;row_i++){
                        MPI_Recv(&A[row_i][(pivot_block+1)*block_size], n/block_size - pivot_block - 1, MPI_DOUBLE, 0, rank ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
            }


            et = omp_get_wtime();
            t7 += (et-st);
            c7++;
            st = omp_get_wtime();

            //MPI_Barrier(MPI_COMM_WORLD);

            if (rank == 0 && verbose){
                printf("[%d:%d]STATUS: Successfully transmited rest blocks\n",n,pivot_block);
            }

            //Update other blocks
            // Parallelize the nested loops for h_block_i and v_block_i

            //printf("[%d]Process %d: I'm responsible for calculating blocks vertically from %d to %d and horizontally from %d to %d (process_block_number=%d)\n",n,rank,pivot_block+1+process_block_number*rank,pivot_block+1+process_block_number*(rank+1),pivot_block + 1,n/block_size,process_block_number);
            for (int v_block_i = pivot_block+1+process_block_number*rank; v_block_i < pivot_block+1+process_block_number*(rank+1); v_block_i++) {
                #pragma omp parallel for
                for (int h_block_i = pivot_block + 1; h_block_i < n/block_size; h_block_i++) {
                    for (int inner_mat_x = h_block_i*block_size; inner_mat_x < (h_block_i+1) * block_size && inner_mat_x<n; inner_mat_x++) {
                        for (int inner_mat_y = v_block_i*block_size; inner_mat_y < (v_block_i+1) * block_size && inner_mat_y<n; inner_mat_y++) {
                            A[inner_mat_y][inner_mat_x] -= A[inner_mat_y][pivot_block*block_size+pivot_i] * A[pivot_block*block_size+pivot_i][inner_mat_x];
                        }
                    }
                }
            }


            et = omp_get_wtime();
            t8 += (et-st);
            c8++;
            st = omp_get_wtime();

            //MPI_Barrier(MPI_COMM_WORLD);

            if (rank == 0 && verbose){
                printf("[%d:%d]STATUS: Successfully calculated rest blocks\n",n,pivot_block);
            }

            if (rank == 0) {
                for (int process = 1; process < numProcesses; process++) { //Send to process 'process'
                    for (int row_i = (pivot_block+1+process*process_block_number)*block_size;row_i < (pivot_block+1+(process+1)*process_block_number)*block_size && row_i < n;row_i++){
                        //printf("[%d]Process %d: Receiving rest block row %d from %d (tag %d,size %d)\n",n,rank,row_i,process,process*10000+row_i,n - (pivot_block +1)*block_size);
                        MPI_Recv(&A[row_i][(pivot_block+1)*block_size], n - (pivot_block +1)*block_size, MPI_DOUBLE, process, process*10000 + row_i ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            } else {
                // Receive data from process 0
                for (int row_i = (pivot_block+1+rank*process_block_number)*block_size;row_i < (pivot_block+1+(rank+1)*process_block_number)*block_size && row_i<n;row_i++){
                    //printf("[%d]Process %d: Sending rest block row %d to %d (tag %d,size %d)\n",n,rank,row_i,0,rank*10000+row_i,n - (pivot_block +1)*block_size);
                    MPI_Send(&A[row_i][(pivot_block+1)*block_size], n - (pivot_block +1)*block_size, MPI_DOUBLE, 0, rank*10000+row_i, MPI_COMM_WORLD);
                }
                //printf("[%d]Process %d: finished sending rest blocks\n",n,rank);
            }

            
            //printf("STATUS: Successfully transmited rest blocks results (%d)\n",rank);



            et = omp_get_wtime();
            t9 += (et-st);
            c9++;

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    printf("Avg pivot block solving:             %fms\n",t1/c1*1000);
    printf("Avg vertical block transmission      %fms\n",t2/c2*1000);
    printf("Avg vertical block solving:          %fms\n",t3/c3*1000);
    printf("Avg horizontal block transmission:   %fms\n",t4/c4*1000);
    printf("Avg horizontal block solving:        %fms\n",t5/c5*1000);
    printf("Avg edge block transmission:         %fms\n",t6/c6*1000);
    printf("Avg rest block transmission:         %fms\n",t7/c7*1000);
    printf("Avg rest block solving:              %fms\n",t8/c8*1000);
    printf("Avg rest block results transmission: %fms\n",t9/c9*1000);
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

void distributed_block_solution(int n,int block_size) {

    double tic, toc;

    double** A = (double**) malloc(n * sizeof(double*));
    double* b = (double*) malloc(n * sizeof(double));
    double* x = (double*) malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        A[i] = (double*) malloc(n * sizeof(double));
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    init_matrix_and_vector(A, b, n);

    tic = omp_get_wtime();
    
    // Perform LU factorization
    distributed_block_lu_factorization(A, n, block_size);
    
    // Solve for x
    solve_lu(A, b, x, n);

    toc = omp_get_wtime();

    if (rank == 0){
        printf("Distributed block time for %d: %3.6f\n", n,(toc - tic));
    }
}


int main(int argc, char **argv) {

    MPI_Init(&argc,&argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Starting LU Decomposition Script (reporting from node %d)\n",rank);
    for (int n = 256; n <= 2048; n = n + 256) {
        distributed_block_solution(n, 64);
    }

    MPI_Finalize();


    /* CODE FOR TESTING WITH A 9x9 MATRIX
    int n=9;
    double** A = (double**) malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++){
        A[i] = (double*)malloc(n*sizeof(double));
        for(int j = 0; j < n; j++){
            A[i][j] = (double)(rand()%100-50);
            //A[i][j] = (double)(1+2*i+j);
        }
    }
    */
}
