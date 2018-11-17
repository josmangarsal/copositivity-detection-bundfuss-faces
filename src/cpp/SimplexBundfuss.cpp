/**
    Simplex Bundfuss body
    SimplexBundfuss.cpp
    Purpose: Simplex Bundfuss object that allow use Facets of any dimension.

    @author J.M.G. Salmerón
    @version 1.0 11/10/2018
    
    Created on: 11/10/2018
*/

#include "SimplexBundfuss.h"

SimplexBundfuss::SimplexBundfuss(int n) {
    this->n = n;

    // Initial simplex
    this->v = new double*[n];
    for (int i = 0; i < n; i++) {
        this->v[i] = new double[n];
        for (int j = 0; j < n; j++) {
            if (i == j)
                this->v[i][j] = 1.0;
            else
                this->v[i][j] = 0.0;
        }
    }

    this->z = std::string(n, '1');
    //this->z_int = bin_to_int(this->z);
    this->m = n;

    this->Q = new double[n*n];
}

SimplexBundfuss::SimplexBundfuss(int n, double** v, int m, std::string z) {
    this->n = n;
    this->v = v;
    this->m = m;
    this->z = z;

    this->Q = new double[n*n];
}

SimplexBundfuss::~SimplexBundfuss() {
    for (int i = 0; i < m; i++)
        delete[] v[i];
    delete[] v;

    delete[] Q;
}

std::tuple<SimplexBundfuss*, SimplexBundfuss*> SimplexBundfuss::divide() {
    double var_lambda1 = var_gamma / (var_gamma-var_alpha);
    double var_lambda3 = var_beta / (var_beta-var_gamma);
    double var_lambda2 = (var_beta-var_gamma) / (var_alpha-2*var_gamma+var_beta);

    double var_lambda_min23 = std::min(var_lambda2, var_lambda3);
    double var_lambda = std::max(var_lambda1, var_lambda_min23);

    // New point
    double* sigma = new double[n];
    for (int k = 0; k < n; k++)
        sigma[k] = var_lambda * v[vertex_i][k] + (1 - var_lambda) * v[vertex_j][k];

    // New vertices. Copy previous
    double** v1 = new double*[n];
    double** v2 = new double*[n];
    for (int i = 0; i < n; i++) {
        v1[i] = new double[n];
        v2[i] = new double[n];
        for (int j = 0; j < n; j++) {
            v1[i][j] = v[i][j];
            v2[i][j] = v[i][j];
        }
    }

    // New vertices. Insert new vertex
    for (int k = 0; k < n; k++) {
        v1[vertex_i][k] = sigma[k];
        v2[vertex_j][k] = sigma[k];
    }

    delete[] sigma;

    return std::make_tuple(new SimplexBundfuss(n, v1, m, z), new SimplexBundfuss(n, v2, m, z));
}

bool SimplexBundfuss::is_copos(double* matrix_a, double eps, double* temp_1_n, std::map<std::string, int> &counters) {
    bool is_copositive = true;

    double res = 0.0;

    std::vector<double> q_vertex;
    std::map<std::tuple<int, int>, double> q_edge;

    for (int i = 0; i < m; i++) {
        for (int j = i; j < m; j++) {
            res = xAy(v[i], matrix_a, v[j], n, temp_1_n);
            counters["xAx"] += 1;

            if (i == j) {
                if (lt(res, 0.0))
                    throw BundfussException("A is not copositive (is_copos)");

                if (i == 0) // && j == 0
                    copos = res; // Init copos

                q_vertex.push_back(res); // q_vertex[i]
                Q[i*m+j] = res;
            } else {
                q_edge[std::make_tuple(i,j)] = res;
                Q[i*m+j] = res;
                Q[j*m+i] = res;
            }
                
            if (lt(res, copos))
                copos = res; // Update copos

            if (lt(res, eps))
                is_copositive = false;
        }
    }

    MapType min_edge = *std::min_element(q_edge.begin(), q_edge.end(), CompareSecond());

    vertex_i = std::get<0>(min_edge.first);
    vertex_j = std::get<1>(min_edge.first);
    var_gamma = min_edge.second;
    var_alpha = q_vertex[vertex_i];
    var_beta = q_vertex[vertex_j];

    if (ge(var_alpha, 0.0) && ge(var_beta, 0.0) && lt(var_gamma, 0.0)) {
        if (eq(var_alpha, 0.0) and eq(var_beta, 0.0))
            throw BundfussException("A is not copositive (Bundfuss Lemma 4 a)");

        double var_lambda_final = var_alpha*var_beta - var_gamma*var_gamma;

        if (lt(var_lambda_final, 0.0))
            throw BundfussException("A is not copositive (Bundfuss Lemma 4 b + Eligius)");
    }

    return is_copositive;
}

std::tuple<std::string, int> SimplexBundfuss::monotonicity(double* a, double* temp_1_m) {
    std::string aux_z = std::string(m, '1');
    std::string new_z = z;

    double** matrix_dm = get_matrix_dm(m);
    double mindiQ = std::numeric_limits<double>::max();

    for (int i = 0; i < m; i++) {
        for (int k = 0; k < m; k++)
            temp_1_m[k] = 0.0;
        product_matrix_matrix(matrix_dm[i], 1, m, Q, m, m, temp_1_m); // temp_1_m = diQ

        mindiQ = std::numeric_limits<double>::max();
        for (int k = 0; k < m; k++) {
            if (lt(temp_1_m[k], mindiQ))
                mindiQ = temp_1_m[k];
        }

        if (ge(mindiQ, 0.0))
            aux_z[i] = '0';
    }

    for (int i = 0; i < m; i++)
        delete[] matrix_dm[i];
    delete[] matrix_dm;

    int aux_i = 0;
    for (int i = 0; i < n; i++) {
        if (z[i] == '1') {
            new_z[i] = aux_z[aux_i];
            aux_i++;
        }
    }

    return std::make_tuple(new_z, count1s(new_z));
}

SimplexBundfuss* SimplexBundfuss::get_hat(std::string new_z, int new_m) {
    double** hat_v = new double*[n];
    int hat_i = 0;
    for (int i = 0; i < n; i++) {
        hat_v[i] = new double[n];
        if (z[i] == '1') {
            for (int j = 0; j < n; j++) {
                hat_v[hat_i][j] = v[i][j];
            }
            hat_i++;
        }
    }

    return new SimplexBundfuss(n, hat_v, new_m, new_z);
}

std::string SimplexBundfuss::to_string() {
    std::stringstream ss;

    ss << "n: " << n << std::endl;
    ss << "m: " << m << std::endl;
    ss << "z: " << z << std::endl;

    ss << "v:" << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ss << v[i][j] << " ";
        }
        ss << std::endl;
    }

    return ss.str();
}