/*
 * lib_gauss.cpp
 *
 *  Created on: May 7, 2022
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "lib_mat.hpp"
#include "lib_mem.hpp"
#include "lib_sort.hpp"
#include "user_types.hpp"

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

int count_leading_zeros(double ** mat, int n, int row) {

    int count = 0;

    while(fabs(mat[row][count]) <= SMALL_NUM && count < n) {
        count++;
    }

    return count;
}

void init_mat_inv(double ** mat_inv, int n) {
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
}

void check_leading_zeros(double ** mat, int n, bool & is_singular) {
    // Check if matrix is singular
    for(int row = 0; row < n; ++row) {
        int num_lead_zeros = count_leading_zeros(mat, n, row);

        if(num_lead_zeros >= row + 1 && !is_singular) {
            printf("Matrix is singular\n");
            is_singular = true;
        }
    }
}

void sort_mat(double * order_arr, int n, double ** mat) {

    double ** mat_ordered = mat2D(n);

    mergesort_mat(mat, n, order_arr, mat_ordered);

    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            mat[row][c] = mat_ordered[row][c];

            // Cut numerically low values
            if(fabs(mat[row][c]) <= SMALL_NUM) {
                mat[row][c] = 0.0;
            }
        }
    }

    free_mat2D(mat_ordered, n);
}

void gauss_jordan(double ** mat, int n, double ** mat_inv) {

    double * order_arr = new double[n];

    // Initialize matrix inverse
    init_mat_inv(mat_inv, n);

    // Initialize singularity flag
    bool is_singular = false;

    // Convert to row echelon form
    for(int c = 0; c < n; ++c) {

        // Sort if under threshold
        if(fabs(mat[c][c]) <= SMALL_NUM) {
            get_order(mat, n, order_arr);

            sort_mat(order_arr, n, mat);

            sort_mat(order_arr, n, mat_inv);

            check_leading_zeros(mat, n, is_singular);
        }

        // Normalize matrix row
        for(int col = c + 1; col < n; ++col) {
            mat[c][col] = fabs(mat[c][c]) <= SMALL_NUM ? 0.0 : mat[c][col] / mat[c][c];;
        }

        // Update row matrix inverse
        for(int col = 0; col < n; ++col) {
            mat_inv[c][col] = fabs(mat[c][c]) <= SMALL_NUM ? 0.0 : mat_inv[c][col] / mat[c][c];;
        }

        mat[c][c] = 1.0;

        // Delete elements in rows below
        for(int row = c + 1; row < n; ++row) {
            if(mat[row][c] != 0) {
                for(int col = c + 1; col < n; ++col) {
                    mat[row][col] = -1.0 * mat[row][c] * mat[c][col] + mat[row][col];
                }
                for(int col = 0; col < n; ++col) {
                    mat_inv[row][col] = -1.0 * mat[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat[row][c] = 0;
            }
        }
    }

    // Backtrace to convert to reduced row echelon form
    for(int c = n - 1; c > 0; --c) {
        for(int row = c - 1; row > -1; --row) {
            if(mat[row][c] != 0) {
                for(int col = 0; col < n; ++col) {
                    mat_inv[row][col] = -1.0 * mat[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat[row][c] = 0;
            }
        }
    }

    // Free allocated space
    delete [] order_arr;
}
