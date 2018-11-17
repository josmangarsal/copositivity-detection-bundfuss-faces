/**
    Simplex Facet body
    SimplexFacet.cpp
    Purpose: Simplex object that allow use Facets of any dimension.

    @author J.M.G. Salmerón
    @version 1.0 18/09/2018
        
    Created on: 18/09/2018
*/

#include "SimplexFacet.h"

/**
 * SimplexFacet constructor and initial simplex
 * 
 * @param n: Dimension
 */
SimplexFacet::SimplexFacet(int n) {
    // Set class variables
    this->n = n;

    // Initial simplex
    this->v = new double[n*n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                this->v[i*n+j] = 1.0;
            else
                this->v[i*n+j] = 0.0;
        }
    }

    this->z = std::string(n, '1');
    this->z_int = bin_to_int(this->z);
    this->m = n;
}

/**
 * SimplexFacet constructor
 * 
 * @param n: Dimension
 * @param v: Vertices matrix
 * @param m: Facet dimension
 * @param z: Facets binary vector
 */
SimplexFacet::SimplexFacet(int n, double* v, int m, std::string z, ULLI z_int) {
    // Set class variables
    this->n = n;
    this->v = v;
    this->m = m;
    this->z = z;
    this->z_int = z_int;
}

/**
 * SimplexFacet destructor
 */
SimplexFacet::~SimplexFacet() {
    // Delete vertices array (matrix)
    delete[] v;
}

/**
 * Check monotony
 * 
 * @param DmAk: Matrix a
 * @return new z
 */
std::tuple<std::string, int> SimplexFacet::monotonicity(double* dmak) {
    std::string aux_z = std::string(n, '1');
    std::string new_z = z;

    int sum_col_j;
    for (int j = 0; j < m; j++) {
        sum_col_j = 0;
        for (int i = 0; i < m; i++) {
            if (ge(dmak[i * m + j], 0.0))
                sum_col_j += 1;
        }

        if (sum_col_j == m)
            aux_z[j] = '0';
    }

    int aux_i = 0;
    for (int i = 0; i < n; i++) {
        if (z[i] == '1') {
            new_z[i] = aux_z[aux_i];
            aux_i++;
        }
    }

    return std::make_tuple(new_z, count1s(new_z));
}

/**
 * Get new facet from this facet and the z
 * 
 * @param new_z: New facet binary vector
 * @param new_m: New number of vertices in the facet
 * @return New facet from this and z
 */
SimplexFacet* SimplexFacet::get_facet(std::string new_z, int new_m) {
    double* new_v = new double[new_m*n];

    for (int i = 0; i < new_m; i++) {
        for (int j = 0; j < n; j++) {
            new_v[i*new_m + j] = 2;
        }
    }

    int new_i = 0;
    for (int i = 0; i < n; i++) {
        if (new_z[i] == '1') {
            for (int j = 0; j < m; j++) {
                new_v[new_i*new_m+j] = v[i*m+j];
            }
            new_i++;
        }
    }

    return new SimplexFacet(n, new_v, new_m, new_z, bin_to_int(new_z));
}

/**
 * SimplexFacet to string
 * 
 * @return String with SimplexFacets features
 */
std::string SimplexFacet::to_string() {
    std::stringstream ss;

    ss << "n: " << n << std::endl;
    ss << "m: " << m << std::endl;
    ss << "z: " << z << std::endl;

    ss << "v:" << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ss << v[i*m + j] << " ";
        }
        ss << std::endl;
    }

    return ss.str();
}