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
#include "user_types.hpp"

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

void swap(double ** mat, int r1, int r2, int n) {
    
    double * row = new double[n];
    
    for(int j = 0; j < n; j++) {
        row[j] = mat[r1][j];
        mat[r1][j] = mat[r2][j];
    }
    
    for(int j = 0; j < n; j++)
        mat[r2][j] = row[j];
    
}

int find(double ** mat, int c, int n) {
    
    for(int row = c + 1; row < n; row++) {
        if(fabs(mat[row][c]) > SMALL_NUM)
            return row;
    }
    
    return -1;
}

void gauss_jordan(double ** mat, int n, double ** mat_inv) {

    double * order_arr = new double[n];

    // Initialize matrix inverse
    init_mat_inv(mat_inv, n);

    // Convert to row echelon form
    for(int c = 0; c < n; ++c) {

        // Sort if under threshold
        if(fabs(mat[c][c]) <= SMALL_NUM) {
            int row = find(mat, c, n);
            swap(mat, row, c, n);
            swap(mat_inv, row, c, n);
        }

        // Normalize matrix row
        for(int col = c + 1; col < n; ++col) {
            mat[c][col] = fabs(mat[c][c]) <= SMALL_NUM ? 0.0 : mat[c][col] / mat[c][c];
        }

        // Update row matrix inverse
        for(int col = 0; col < n; ++col) {
            mat_inv[c][col] = fabs(mat[c][c]) <= SMALL_NUM ? 0.0 : mat_inv[c][col] / mat[c][c];
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

    // Backtrack to convert to reduced row echelon form
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
