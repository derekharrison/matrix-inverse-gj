/*
 * lib_mem.cpp
 *
 *  Created on: May 7, 2022
 *      Author: d-w-h
 */

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

