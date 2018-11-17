/**
    Matrices body
    Matrices.h
    Purpose: Functions to matrices calculations.

    @author J.M.G. Salmerón
    @version 1.1 21/09/2018
    
    Created on: 20/09/2018
*/

#include "Matrices.h"

/**
 * Performs matrix matrix multiplication
 * 
 * @param a: First matrix called A
 * @param ar: Number of rows of matrix A
 * @param ac: Number of columns of matrix A
 * @param b: Second matrix called B
 * @param br: Number of rows of matrix B
 * @param bc: Number of columns of matrix B
 * @return c, ar, bc: Tuple with result C matrix with size ar rows and bc columns
 */
void product_matrix_matrix(double* a, int ar, int ac, double* b, int br, int bc, double* c_ar_bc) {
    // C := alpha*op(A)*op(B) + beta*C
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ar, bc, br, 1.0, a, ac, b, bc, 1.0, c_ar_bc, bc);
}

/**
 * Performs vector vector multiplication
 * 
 * @param a: Vector A
 * @param b: Vector B
 * @param n: Size of vector A and B
 * @return The result of the dot product of A and B, if n is positive. Otherwise, returns 0.
 */
double product_vector_vector(double* a, double* b, int n) {
	return cblas_ddot(n, a, 1, b, 1);
}

/**
 * Performs matrix matrix substraction
 * 
 * @param a: Matrix A
 * @param b: Matrix B
 * @param r: Rows A and B
 * @param c: Columns A and B
 * @return Matrix A - B
 */
double* matrix_minus_matrix(double* a, double* b, int r, int c) {
	double* res = new double[r*c];
    for (int i = 0; i < r*c; i++) {
        res[i] = a[i] - b[i];
    }
    return res;
}

/**
 * Performs scalar matrix multiplication
 * 
 * @param s: Scalar
 * @param a: Matrix A
 * @param r: Rows A 
 * @param c: Columns A
 * @return Matrix s * A
 */
double* product_scalar_matrix(double s, double* a, int r, int c) {
    double* res = new double[r*c];
    for (int i = 0; i < r*c; i++) {
        res[i] = s * a[i];
    }
    return res;
}

/**
 * Get matrix diagonal
 * 
 * @param a: Matrix A
 * @param r: Rows A 
 * @param c: Columns A
 * @return Array matrix diagonal
 */
double* matrix_diagonal(double* a, int r, int c) {
    double* diag = new double[c];
    for (int i = 0, j = 0; i < c; i++, j++) {
        diag[i] = a[i*r+j];
    }
    return diag;
}

/**
 * Matrix elements are positives
 * 
 * @param a: Matrix A
 * @param r: Rows A 
 * @param c: Columns A
 * @return True if all elements of the matrix are greater than accuracy
 */
bool ge_matrix_zero(double* a, int r, int c) {
    for (int i = 0; i < r*c; i++)
        if (lt(a[i], 0.0))
            return false;
    return true;
}

/**
 * Matrix is semipositive definite
 * 
 * @param hk: Matrix Hk
 * @param m: Hk size
 * @param eye_m: Temp matrix mxm to store identity
 * @param ready: Temp matrix mxm to store result
 * @return True if A si semipositive_definite
 */
bool matrix_semipositive_definite(double* hk, int m) {
    //TODO Improve eye and ready
    double* eye_m = new double[m*m];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i == j)
                eye_m[i*m+j] = 1.0 * Chol_accuracy;
            else
                eye_m[i*m+j] = 0.0;
        }
    }

    double* ready = new double[m*m];
    for (int i = 0; i < m*m; i++)
        ready[i] = 0.0;
    
    vdAdd(m*m, hk, eye_m, ready);

    int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', m, ready, m);

    delete[] eye_m;
    delete[] ready;

    return info == 0;
}

/**
 * Minimum eigenvalue
 * 
 * @param hk: Matrix Hk
 * @param m: Hk size
 * @return Min real eigenvalue of Hk
 */
double min_real_eig(double* hk, int m, double* w) {
    int info;
    int lwork;
    double wkopt;
    double* work;

    int n = m;
    int lda = m;

    lwork = -1;

    dsyev("Vectors", "Upper", &n, hk, &lda, w, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = new double[lwork];
    dsyev("Vectors", "Upper", &n, hk, &lda, w, work, &lwork, &info);

    return w[0];
}

/**
 * Evaluate function in point x with matrix A
 * 
 * @param x: Point to evaluate
 * @param a: Matrix A
 * @param an: Size of square matrix A
 * @return c: Function value on point x
 */
double f(double* x, double* a, int an, double* temp_1_an) {
	product_matrix_matrix(x, 1, an, a, an, an, temp_1_an); // temp_1_an (1 x m)
    return product_vector_vector(temp_1_an, x, an); // xax
}

/**
 * 
 */
double xAy(double* x, double* a, double* y, int n, double* temp_1_n) {
    for (int i = 0; i < n; i++)
        temp_1_n[i] = 0.0;
    product_matrix_matrix(x, 1, n, a, n, n, temp_1_n); // temp_1_an (1 x m)
    return product_vector_vector(temp_1_n, y, n); // xay
}

/**
 * Performs matrix matrix matrix multiplication
 * 
 * r = c = m
 * 
 * @param v: Matrix V
 * @param a: Matrix A
 * @param r: Rows
 * @param c: Columns
 * @param vav: Matrix VAV of same size
 */
void get_matrix_vav(double* v, int vr, int vc, double* a, int ar, int ac, double* vav, double* temp_vr_ac) {
	product_matrix_matrix(v, vr, vc, a, ar, ac, temp_vr_ac); // temp_r_c (vr x ac)
	product_matrix_matrix(temp_vr_ac, vr, ac, v, vr, vc, vav); // vav (vr x vc)
}

/**
 * Calculate D_m matrix
 * 
 * @param m: Facet dimension
 * @return Square matrix D_m of m x m
 */
double** get_matrix_dm(int m) {
    double** dm = new double*[m];
    for (int i = 0; i < m; i++)
        dm[i] = new double[m];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i == j)
                dm[i][j] = (m - 1.0) / m; //dm[i*m+j] = (m - 1.0) / m
            else
                dm[i][j] = - 1.0 / m; //dm[i*m+j] = - 1.0 / m
        }
    }

    return dm;
}

/**
 * Calculate A_m matrix
 * 
 * @param a: Matrix A
 * @param n: Size of A
 * @param z: Facet vector
 * @param m: Facet dimension
 * @param ak: Square matrix A_k of m x m
 */
void get_matrix_ak(double* a, int n, std::string z, int m, double* ak) {    
    int ik = 0, jk = 0;
    for (int i = 0; i < n; i++) {
        jk = 0;
        for (int j = 0; j < n; j++) {
            if (z[i] == '1' && z[j] == '1') {
                ak[ik*m+jk] = a[i*n+j];
                jk++;
            }
        }
        if (z[i] == '1')
            ik++;
    }
}

/**
 * Matrix to string
 * 
 * @param a: Matrix A
 * @param r: Rows A 
 * @param c: Columns A
 * @return Printed matrix
 */
std::string matrix_to_string(double* a, int r, int c) {
    std::stringstream ss;
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            if (ge(a[i*r + j], 0.0))
                ss << " " << a[i*r + j] << " ";
            else
                ss << a[i*r + j] << " ";
        }
        ss << std::endl;
    }
    return ss.str();
}

/**
 * vector to string
 * 
 * @param v: Vector V
 * @param n: Size V 
 * @return Printed vector
 */
std::string vector_to_string(double* v, int n) {
    std::stringstream ss;
    for (int i = 0; i < n; i++)
        ss << v[i] << " ";
    return ss.str();
}

/**
 * Dm matrices delete
 */
void DmMatrix::free() {
    for (auto [m, dm] : dm_matrices) { // m is a must, sorry for the warning
        delete[] dm;
    }
}

/**
 * Calculate and store D_m matrix as row major
 * 
 * @param m: Facet dimension
 * @return Square matrix D_m of m x m
 */
double* DmMatrix::get_dm(int m) {
    if (dm_matrices.find(m) != dm_matrices.end())
        return dm_matrices[m];

    double* dm = new double[m*m];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i == j)
                dm[i*m+j] = (m - 1.0) / m;
            else
                dm[i*m+j] = - 1.0 / m;
        }
    }

    dm_matrices[m] = dm;

    return dm;
}
