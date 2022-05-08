/*
 * lib_sort.cpp
 *
 *  Created on: May 7, 2022
 *      Author: d-w-h
 */

#include <math.h>

#include "lib_mem.hpp"
#include "user_types.hpp"

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

void mergesort_mat(double ** mat, int n, double * order_arr, double ** ordered_mat) {

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
