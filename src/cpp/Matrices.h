/**
    Matrices header
    Matrices.h
    Purpose: Functions to matrices calculations.

    @author J.M.G. Salmerón
    @version 1.1 21/09/2018
    
    Created on: 20/09/2018

    Help: 
    - BLAS and Sparse BLAS Routines
    https://software.intel.com/en-us/mkl-developer-reference-c-blas-and-sparse-blas-routines
*/

#ifndef MATRICES_H_
#define MATRICES_H_

#include <tuple>
#include <map> // DmMatrix class
#include <sstream> // to_string
#include <algorithm> // min_element
#include <limits> // Max double value

#include "mkl.h" // blas

#include "Compare.h" // lt

// #include <iostream>

void product_matrix_matrix(double* a, int ar, int ac, double* b, int br, int bc, double* c_ar_bc);
double product_vector_vector(double* a, double* b, int n);
double* matrix_minus_matrix(double* a, double* b, int r, int c);
double* product_scalar_matrix(double s, double* a, int r, int c);

double* matrix_diagonal(double* a, int r, int c);
bool ge_matrix_zero(double* a, int r, int c);
bool matrix_semipositive_definite(double* hk, int m); // TODO Cuando uso memoria temporal sale mal
double min_real_eig(double* hk, int m, double* w);

double f(double* x, double* a, int an, double* temp_1_an); // xAx
double xAy(double* x, double* a, double* y, int n, double* temp_1_an); // xAy
void get_matrix_vav(double* v, int vr, int vc, double* a, int ar, int ac, double* vav, double* temp_vr_ac); // VAV

double** get_matrix_dm(int m);
void get_matrix_ak(double* a, int n, std::string z, int m, double* ak);

std::string matrix_to_string(double* a, int r, int c);
std::string vector_to_string(double* v, int n);

class DmMatrix {
    std::map<int, double*> dm_matrices;

public:
    void free();
    double* get_dm(int m);
};

#endif /* MATRICES_H_ */