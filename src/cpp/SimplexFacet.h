/**
    Simplex Facet header
    SimplexFacet.h
    Purpose: Simplex object that allow use Facets of any dimension.

    @author J.M.G. Salmerón
    @version 1.0 18/09/2018
    
    Created on: 18/09/2018
*/

#ifndef SIMPLEXFACET_H_
#define SIMPLEXFACET_H_

#include <tuple>
#include <sstream> // to_string
#include <string> // z
#include <algorithm> // min_element
#include <iterator> // min_element

#include "Compare.h" // gt
#include "Matrices.h" // dm, vav
#include "Binary.h" // count1s

class SimplexFacet {
    int n; // Dimension space
    double* v; // Vertices
    int m; // Simplex dimension
    std::string z; // Binary vector of specific dimension
    ULLI z_int; // Binary vector of specific dimension

public:
    //SimplexFacet(int n, double* v, int m, bool* z) : n(n), v(v), m(m), z(z) {}
    SimplexFacet(int n);
    SimplexFacet(int n, double* v, int m, std::string z, ULLI z_int);
    ~SimplexFacet();

    std::tuple<std::string, int> monotonicity(double* dmak);
    SimplexFacet* get_facet(std::string new_z, int new_m);

    inline double* get_v() { return v; }
    inline std::string get_z() { return z; }
    inline ULLI get_z_int() { return z_int; }
    inline int get_m() { return m; }

    std::string to_string();
};

#endif /* SIMPLEXFACET_H_ */