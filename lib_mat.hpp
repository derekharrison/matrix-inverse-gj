/*
 * lib_mat.hpp
 *
 *  Created on: May 7, 2022
 *      Author: d-w-h
 */

#ifndef LIB_MAT_HPP_
#define LIB_MAT_HPP_


void init_mat(int n, double ** mat);
void set_mat(double ** mat, int n, double ** mat_store);
void print_mat(double ** mat, int n);
void mat_mult_sq(double ** A, double ** A_inv, int n, double ** mat_res);


#endif /* LIB_MAT_HPP_ */
