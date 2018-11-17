/**
    Simplex Bundfuss header
    SimplexBundfuss.h
    Purpose: Simplex Bundfuss object that allow use Facets of any dimension.

    @author J.M.G. Salmerón
    @version 1.0 11/10/2018
    
    Created on: 11/10/2018
*/

#ifndef SIMPLEXBUNDFUSS_H_
#define SIMPLEXBUNDFUSS_H_

#include <tuple>
#include <sstream> // to_string
#include <string> // z

#include "Compare.h" // gt
#include "Matrices.h" // dm, vav
#include "Binary.h" // count1s

#include <vector> // is_copos
#include <map>
#include <algorithm>
#include <tuple>

#include <exception>

#include <limits> // Mono

class BundfussException: public std::exception {
    const char* e;

public:
    BundfussException(const char* e) : e(e) {}

    virtual const char* what() const throw() {
        return e;
    }
};

typedef std::pair<std::tuple<int, int>, double> MapType;
struct CompareSecond {
    bool operator() (const MapType &left, const MapType &right) const {
        return left.second < right.second;
    }
};

class SimplexBundfuss {
    int n; // Dimension space
    double** v; // Vertices

    double* Q; // Eval matrix
    
    int vertex_i, vertex_j; // Edge
    double var_alpha; // Alpha
    double var_beta; // Beta
    double var_gamma; // Gamma
    
    double copos; // Copos/LB

    int m; // Simplex dimension
    std::string z; // Binary vector of specific dimension

public:
    //SimplexBundfuss(int n, double* v, int m, bool* z) : n(n), v(v), m(m), z(z) {}
    SimplexBundfuss(int n);
    SimplexBundfuss(int n, double** v, int m, std::string z);
    ~SimplexBundfuss();

    std::tuple<SimplexBundfuss*, SimplexBundfuss*> divide();
    bool is_copos(double* matrix_a, double eps, double* temp_1_n, std::map<std::string, int> &counters);
    std::tuple<std::string, int> monotonicity(double* a, double* temp_1_m);
    SimplexBundfuss* get_hat(std::string new_z, int new_m);

    inline double** get_v() { return v; }
    inline double get_copos() { return copos; }
    inline std::string get_z() { return z; }
    inline int get_m() { return m; }

    std::string to_string();
};

#endif /* SIMPLEXBUNDFUSS_H_ */