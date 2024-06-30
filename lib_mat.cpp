/*
 * lib_mat.cpp
 *
 *  Created on: May 7, 2022
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "user_types.hpp"

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
            double rand_num_loc = rand_di(-25, 25);
            if(fabs(rand_num_loc) <= SMALL_NUM) { rand_num_loc = 0.0; }
            mat[i][j] = rand_num_loc;
        }
    }
}

void set_mat(double ** mat, int n, double ** mat_store) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            mat_store[i][j] = mat[i][j];
        }
    }
}

void print_mat(double ** mat, int n) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            printf("%.3f ", fabs(mat[i][j]));
        }
        printf("\n");
    }
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
