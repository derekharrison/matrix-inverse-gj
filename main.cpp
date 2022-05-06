//
//  main.cpp
//  gauss-jordan
//
//  Created by mndx on 17/04/2022.
//

#include <math.h>
#include <stdio.h>
#include <time.h>

const int MAX_INT = 1215752192;
const double SMALL_NUM = 1e-10;

typedef struct order_array_elem {
    int old_row;
    double val;
} oa_elem_t;

double ** mat2D(int n) {

    double ** mat = new double * [n];

    for(int i = 0; i < n; ++i)
        mat[i] = new double[n];

    return mat;
}

void free_mat2D(double ** mat, int n) {

    for(int i = 0; i < n; ++i)
        delete [] mat[i];

    delete [] mat;
}

void get_order(double ** mat, int n, double * order_arr) {
    for(int row = 0; row < n; ++row) {
        int order = 0;
        while(fabs(mat[row][order]) <= SMALL_NUM && order < n) {
            mat[row][order] = 0.0;
            order++;
        }
        order_arr[row] = order;
    }
}

void merge(oa_elem_t A[], int p, int q, int r) {
    int size_r, size_l;
    int i, j;
    size_l = q - p + 1;
    size_r = r - q;
    oa_elem_t L[size_l + 1];
    oa_elem_t R[size_r + 1];
    i = 0;
    j = 0;
    for(int n = p; n < q + 1; ++n) {
        L[i] = A[n];
        ++i;
    }
    L[size_l].val = MAX_INT;
    for(int n = q + 1; n < r + 1; ++n) {
        R[j] = A[n];
        ++j;
    }
    R[size_r].val = MAX_INT;
    i = 0;
    j = 0;
    for(int n = p; n < r + 1; ++n) {
        if(L[i].val < R[j].val) {
            A[n] = L[i];
            ++i;
        }
        else {
            A[n] = R[j];
            ++j;
        }
    }
}

void merge_sort(oa_elem_t A[], int p, int r) {
    int q;
    if(p < r) {
        q = (p + r)/2;
        merge_sort(A, p, q);
        merge_sort(A, q + 1, r);
        merge(A, p, q, r);
    }
}

void mergesort(oa_elem_t A[], int size) {
    merge_sort(A, 0, size - 1);
}

void print_mat(double ** mat, int n) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            printf("%.3f ", mat[i][j]);
        }
        printf("\n");
    }
}

void sort_mat(double ** mat, int n, double * order_arr, double ** ordered_mat) {

    oa_elem_t * order_array = new oa_elem_t[n];

    for(int row = 0; row < n; ++row) {
        order_array[row].old_row = row;
        order_array[row].val = order_arr[row];
    }

    mergesort(order_array, n);

    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            int old_row = order_array[row].old_row;
            ordered_mat[row][c] = mat[old_row][c];
            if(fabs(ordered_mat[row][c]) <= SMALL_NUM) {
                ordered_mat[row][c] = 0.0;
            }
        }
    }
    
    delete [] order_array;
}

int count_leading_zeros(double ** mat, int n, int row) {

    int count = 0;

    while(fabs(mat[row][count]) <= SMALL_NUM && count < n) {
        count++;
    }

    return count;
}

void mat_mult_sq(double ** A, double ** A_inv, int n, double ** mat_res) {

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double sum_loc = 0;

            for(int k = 0; k < n; ++k) {
                sum_loc = sum_loc + A[i][k] * A_inv[k][j];
            }

            mat_res[i][j] = sum_loc;
        }
    }
}

void singularity_check(double ** mat_ref, int n, bool & is_singular) {
    
    for(int row = 0; row < n; ++row) {
        bool all_zeros_c = true;
        for(int col = 0; col < n; ++col) {
            if(fabs(mat_ref[row][col]) > SMALL_NUM) {
                all_zeros_c = false;
            }
        }
        if(all_zeros_c && !is_singular) {
            printf("Matrix is singular\n");
            is_singular = true;
        }
    }
    
    for(int col = 0; col < n; ++col) {
        bool all_zeros_r = true;
        for(int row = 0; row < n; ++row) {
            if(fabs(mat_ref[row][col]) > SMALL_NUM) {
                all_zeros_r = false;
            }
        }
        if(all_zeros_r && !is_singular) {
            printf("Matrix is singular\n");
            is_singular = true;
        }
    }
}

void cut_low_vals(double ** mat, int n) {
    
    for(int row = 0; row < n; ++row) {
        for(int col = 0; col < n; ++col) {
            if(fabs(mat[row][col]) <= SMALL_NUM) {
                mat[row][col] = 0.0;
            }
        }
    }
}

void gauss_jordan(double ** mat, int n, double ** mat_inv) {

    double ** mat_ref = mat2D(n);
    double ** mat_ordered = mat2D(n);
    double ** mat_inv_ordered = mat2D(n);
    double * order_arr = new double[n];

    // Initialize matrix inverse
    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            if(c == row) {
                mat_inv[row][c] = 1.0;
            }
            else {
                mat_inv[row][c] = 0.0;
            }
        }
    }
    
    // Sort the input matrix
    get_order(mat, n, order_arr);

    sort_mat(mat, n, order_arr, mat_ordered);
    
    sort_mat(mat_inv, n, order_arr, mat_inv_ordered);

    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            mat_ref[row][c] = mat_ordered[row][c];
            mat_inv[row][c] = mat_inv_ordered[row][c];
        }
    }
    
    // Initialize singularity flag
    bool is_singular = false;

    // Check if input matrix is singular
    singularity_check(mat_ref, n, is_singular);

    // Convert to row echelon form
    for(int c = 0; c < n; ++c) {
        
        // Sort if under threshold
        if(fabs(mat_ref[c][c]) <= SMALL_NUM) {
            get_order(mat_ref, n, order_arr);

            sort_mat(mat_ref, n, order_arr, mat_ordered);

            sort_mat(mat_inv, n, order_arr, mat_inv_ordered);

            for(int row = 0; row < n; ++row) {
                for(int c = 0; c < n; ++c) {
                    mat_ref[row][c] = mat_ordered[row][c];
                    mat_inv[row][c] = mat_inv_ordered[row][c];
                }
            }
        }

        // Normalize matrix row
        for(int col = c + 1; col < n; ++col) {
            mat_ref[c][col] = mat_ref[c][col] / (mat_ref[c][c] + SMALL_NUM);
        }

        // Update row matrix inverse
        for(int col = 0; col < n; ++col) {
            mat_inv[c][col] = mat_inv[c][col] / (mat_ref[c][c] + SMALL_NUM);
        }

        mat_ref[c][c] = 1.0;

        // Delete elements in rows below
        for(int row = c + 1; row < n; ++row) {
            if(mat_ref[row][c] != 0) {
                for(int col = c + 1; col < n; ++col) {
                    mat_ref[row][col] = -1.0 * mat_ref[row][c] * mat_ref[c][col] + mat_ref[row][col];
                }
                for(int col = 0; col < n; ++col) {
                    mat_inv[row][col] = -1.0 * mat_ref[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat_ref[row][c] = 0;
            }
            
            int num_lead_zeros = count_leading_zeros(mat_ref, n, row);
            
            if(num_lead_zeros >= n && !is_singular) {
                printf("Matrix is singular\n");
                is_singular = true;
            }
        }

    }

    // Backtrace to convert to reduced row echelon form
    for(int c = n - 1; c > 0; --c) {
        for(int row = c - 1; row > -1; --row) {
            if(mat_ref[row][c] != 0) {
                for(int col = 0; col < n; ++col) {
                    mat_inv[row][col] = -1.0 * mat_ref[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat_ref[row][c] = 0;
            }
        }
    }
    
    // Check if matrix is singular
    singularity_check(mat_ref, n, is_singular);
    
    // Free allocated space
    free_mat2D(mat_ref, n);
    free_mat2D(mat_ordered, n);
    free_mat2D(mat_inv_ordered, n);
    delete [] order_arr;
}

double rand_num(double min, double max) {

    double val = (double) rand() / (RAND_MAX + 1.0);

    return val * (max - min) - (max - min) / 2;
}

double rand_di(int min, int max) {
    
    return rand() % (max - min) + min;
}

void init_mat(int n, double ** mat) {
    
    srand((unsigned) time(NULL));
    
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double rand_num_loc = rand_num(-25, 25);
            if(fabs(rand_num_loc) <= SMALL_NUM) { rand_num_loc = 0.0; }
            mat[i][j] = rand_num_loc;
        }
    }
}

int main(int argc, char * argv[]) {

    // Declarations
    int n = 10;

    // Allocate space for matrices
    double ** mat = mat2D(n);
    double ** mat_inv = mat2D(n);
    double ** mat_prod = mat2D(n);

    // Populate matrix mat with some data
    init_mat(n, mat);

    // Compute inverse using Gauss-Jordan method
    gauss_jordan(mat, n, mat_inv);

    // Verify compuation
    mat_mult_sq(mat, mat_inv, n, mat_prod);

    // Print results
    print_mat(mat_prod, n);
    
    // Free allocated space
    free_mat2D(mat, n);
    free_mat2D(mat_inv, n);
    free_mat2D(mat_prod, n);

    return 0;
}



