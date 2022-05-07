//
//  main.cpp
//  gauss-jordan
//
//  Created by mndx on 17/04/2022.
//

#include "lib_gauss.hpp"
#include "lib_mat.hpp"
#include "lib_mem.hpp"
#include "user_types.hpp"

int main(int argc, char * argv[]) {

    // Declarations
    int n = 10;

    // Allocate space for matrices
    double ** mat = mat2D(n);
    double ** mat_inv = mat2D(n);
    double ** mat_prod = mat2D(n);
    double ** mat_store = mat2D(n);

    // Populate matrix mat with some data
    init_mat(n, mat);

    // Store initial matrix mat
    set_mat(mat, n, mat_store);

    // Compute inverse using Gauss-Jordan method
    gauss_jordan(mat, n, mat_inv);

    // Verify computation
    mat_mult_sq(mat_store, mat_inv, n, mat_prod);

    // Print results
    print_mat(mat_inv, n);
    
    // Free allocated space
    free_mat2D(mat, n);
    free_mat2D(mat_inv, n);
    free_mat2D(mat_prod, n);
    free_mat2D(mat_store, n);

    return 0;
}



