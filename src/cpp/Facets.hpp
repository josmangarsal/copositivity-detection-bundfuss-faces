/**
    Facets
    Facets.hpp
    Purpose: Facets refinement

    @author J.M.G. Salmerón
    @version 1.0 11/10/2018
        
    Created on: 11/10/2018
*/

// Libs
#include <iostream> // cout

#include <Python.h> // Python
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "termcolor.hpp" // Output color

// Headers
#include "SimplexFacet.h"
#include "Matrices.h"
#include "Binary.h"
#include "ChronoTime.hpp"

#ifndef FACETS_H_
#define FACETS_H_

/**
 * //TODO Add doc
 */
bool py_init() { // return bool is a must here, sorry for the warning
    Py_Initialize();

    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    import_array();
}

/**
 * //TODO Add doc
 */
void py_end() {
    Py_Finalize();
}

/**
 * //TODO Add doc
 */
void py_lp(int null_space_matrix_rows, int m_facet, int m_level, double* matrix_dmak, double* x_star, bool &solution) {
    PyObject *py_name, *py_module, *py_function;
    PyObject *py_args, *py_value, *py_rslt;

    py_name = PyUnicode_FromString("ScipyOptimize");

    py_module = PyImport_Import(py_name);
    Py_DECREF(py_name);

    if (py_module != NULL) {
        py_function = PyObject_GetAttrString(py_module, "scipy_optimize_linprog");

        if (py_function && PyCallable_Check(py_function)) {
            py_args = PyTuple_New(2);
            //Arg 0 m
            py_value = PyLong_FromLong(m_facet);
            if (!py_value) {
                Py_DECREF(py_args);
                Py_DECREF(py_module);
                std::cerr << "Py++: Cannot convert argument" << std::endl;
                exit(1);
            }
            PyTuple_SetItem(py_args, 0, py_value);
            //Arg 1 DmAk
            const int nd = 1;
            npy_intp vdim[1]{m_level*m_level};

            py_value = PyArray_SimpleNewFromData(nd, vdim, NPY_DOUBLE, reinterpret_cast<void*>(matrix_dmak));
            if (!py_value) {
                Py_DECREF(py_args);
                Py_DECREF(py_module);
                std::cerr << "Py++: Cannot convert argument" << std::endl;
                exit(1);
            }
            PyTuple_SetItem(py_args, 1, py_value);

            py_rslt = PyObject_CallObject(py_function, py_args);
            Py_DECREF(py_args);
            if (py_rslt != NULL) {
                PyObject *py_x, *py_g;

                PyArg_UnpackTuple(py_rslt, "ref", 1, 2, &py_x, &py_g);

                if (PyList_Check(py_x) == 0) {
                    std::cout << "LP: Not solution found" << std::endl;
                    std::cerr << "Py++: Cannot convert to list" << std::endl;
                    exit(1);
                }

                if (null_space_matrix_rows != (int) PyList_Size(py_x)) {
                    std::cerr << "Py++: Lists size does not equal" << std::endl;
                    exit(1);
                }

                for (int i = 0; i < (int) PyList_Size(py_x); i++) {
                    py_value = PyList_GetItem(py_x, i);
                    x_star[i] = PyFloat_AsDouble(py_value);
                }
                solution = ge(PyFloat_AsDouble(py_g), 0.0);
            } else {
                Py_DECREF(py_function);
                Py_DECREF(py_module);
                PyErr_Print();
                std::cerr << "Py++: Call failed" << std::endl;
                exit(1);
            }
        } else {
            if (PyErr_Occurred())
                PyErr_Print();
            std::cerr << "Py++: Cannot find function scipy_optimize_linprog" << std::endl;
            exit(1);
        }
        Py_XDECREF(py_function);
        Py_DECREF(py_module);
    } else {
        PyErr_Print();
        std::cerr << "Py++: Failed to load ScipyOptimize" << std::endl;
        exit(1);
    }
}

/**
 * //TODO Add doc
 */
void py_null(int m_facet, int m_level, double* matrix_dmak, int &nn, double* matrix_b) {
    PyObject *py_name, *py_module, *py_function;
    PyObject *py_args, *py_value, *py_rslt;

    py_name = PyUnicode_FromString("ScipyOptimize");

    py_module = PyImport_Import(py_name);
    Py_DECREF(py_name);

    if (py_module != NULL) {
        py_function = PyObject_GetAttrString(py_module, "scipy_optimize_null");

        if (py_function && PyCallable_Check(py_function)) {
            py_args = PyTuple_New(2);
            //Arg 0 m
            py_value = PyLong_FromLong(m_facet);
            if (!py_value) {
                Py_DECREF(py_args);
                Py_DECREF(py_module);
                std::cerr << "Py++: Cannot convert argument" << std::endl;
                exit(1);
            }
            PyTuple_SetItem(py_args, 0, py_value);
            //Arg 1 DmAk
            const int nd = 1;
            npy_intp vdim[1]{m_level*m_level};

            py_value = PyArray_SimpleNewFromData(nd, vdim, NPY_DOUBLE, reinterpret_cast<void*>(matrix_dmak));
            if (!py_value) {
                Py_DECREF(py_args);
                Py_DECREF(py_module);
                std::cerr << "Py++: Cannot convert argument" << std::endl;
                exit(1);
            }
            PyTuple_SetItem(py_args, 1, py_value);

            py_rslt = PyObject_CallObject(py_function, py_args);
            Py_DECREF(py_args);
            if (py_rslt != NULL) {
                PyObject *py_b, *py_nn;

                PyArg_UnpackTuple(py_rslt, "ref", 1, 2, &py_b, &py_nn);

                if (PyList_Check(py_b) == 0) {
                    std::cout << "Null: Not solution found" << std::endl;
                    std::cerr << "Py++: Cannot convert to list" << std::endl;
                    exit(1);
                }

                for (int i = 0; i < (int) PyList_Size(py_b); i++) {
                    py_value = PyList_GetItem(py_b, i);
                    matrix_b[i] = PyFloat_AsDouble(py_value);
                }
                nn = PyLong_AsLong(py_nn);
            } else {
                Py_DECREF(py_function);
                Py_DECREF(py_module);
                PyErr_Print();
                std::cerr << "Py++: Call failed" << std::endl;
                exit(1);
            }
        } else {
            if (PyErr_Occurred())
                PyErr_Print();
            std::cerr << "Py++: Cannot find function scipy_optimize_null" << std::endl;
            exit(1);
        }
        Py_XDECREF(py_function);
        Py_DECREF(py_module);
    } else {
        PyErr_Print();
        std::cerr << "Py++: Failed to load ScipyOptimize" << std::endl;
        exit(1);
    }
}

/**
 * //TODO Add doc
 */
void py_chol(int m_facet, int m_level, double* matrix_hk, int &p) {
    PyObject *py_name, *py_module, *py_function;
    PyObject *py_args, *py_value, *py_rslt;

    py_name = PyUnicode_FromString("ScipyOptimize");

    py_module = PyImport_Import(py_name);
    Py_DECREF(py_name);

    if (py_module != NULL) {
        py_function = PyObject_GetAttrString(py_module, "scipy_optimize_chol");

        if (py_function && PyCallable_Check(py_function)) {
            py_args = PyTuple_New(2);
            //Arg 0 m
            py_value = PyLong_FromLong(m_facet);
            if (!py_value) {
                Py_DECREF(py_args);
                Py_DECREF(py_module);
                std::cerr << "Py++: Cannot convert argument" << std::endl;
                exit(1);
            }
            PyTuple_SetItem(py_args, 0, py_value);
            //Arg 1 DmAk
            const int nd = 1;
            npy_intp vdim[1]{m_level*m_level};

            py_value = PyArray_SimpleNewFromData(nd, vdim, NPY_DOUBLE, reinterpret_cast<void*>(matrix_hk));
            if (!py_value) {
                Py_DECREF(py_args);
                Py_DECREF(py_module);
                std::cerr << "Py++: Cannot convert argument" << std::endl;
                exit(1);
            }
            PyTuple_SetItem(py_args, 1, py_value);

            py_rslt = PyObject_CallObject(py_function, py_args);
            Py_DECREF(py_args);
            if (py_rslt != NULL) {
                p = PyLong_AsLong(py_rslt);
            } else {
                Py_DECREF(py_function);
                Py_DECREF(py_module);
                PyErr_Print();
                std::cerr << "Py++: Call failed" << std::endl;
                exit(1);
            }
        } else {
            if (PyErr_Occurred())
                PyErr_Print();
            std::cerr << "Py++: Cannot find function scipy_optimize_null" << std::endl;
            exit(1);
        }
        Py_XDECREF(py_function);
        Py_DECREF(py_module);
    } else {
        PyErr_Print();
        std::cerr << "Py++: Failed to load ScipyOptimize" << std::endl;
        exit(1);
    }
}

bool chol = true;

/**
 * //TODO Add doc
 * //Input: See declaration
 */
int eval_facet(
    double* matrix_a, int n, int m_level, DmMatrix &dm_matrices,                     // General
    SimplexFacet* facet_k, int m_facet, std::string z_facet_k,                       // Facet
    double* matrix_ak, double* matrix_hk, double* temp_1_m, double* matrix_dmak,     // Temps matrices
    double* matrix_dm, double* x_star,                                               // Temps matrices
    std::map<std::string, int> &counters, int &m_level_chol_evals,                   // Counters
    std::string &mono_z, int &mono_m
) {
    for (int i = 0; i < m_facet*m_facet; i++) {
        matrix_ak[i] = 0.0;
        matrix_hk[i] = 0.0;
        matrix_dmak[i] = 0.0;
    }
    for (int i = 0; i < m_facet; i++) {
        temp_1_m[i] = 0.0;
        x_star[i] = 0.0;
    }

    get_matrix_ak(matrix_a, n, z_facet_k, m_facet, matrix_ak);
    if (ge_matrix_zero(matrix_ak, m_facet, m_facet)) {
        counters["matrixpositive"] += 1;
        return 2;
    }

    matrix_dm = dm_matrices.get_dm(m_facet);
    product_matrix_matrix(matrix_dm, m_facet, m_facet, matrix_ak, m_facet, m_facet, matrix_dmak);

    auto [new_z, new_m] = facet_k->monotonicity(matrix_dmak);
    if (m_facet > new_m) { // Significa que he reducido la dimension, el resto de facetas de Fk las quito
        counters["monotonicity"] += 1;
        mono_z = new_z;
        mono_m = new_m;
        //Leo
        //std::cout << "Facet monotonous at level " << m_level << std::endl;
        //exit(1);
        return 1;
    }

    if (chol) {
        product_matrix_matrix(matrix_dmak, m_facet, m_facet, matrix_dm, m_facet, m_facet, matrix_hk); // DmAkDm
        m_level_chol_evals += 1;
        if (!matrix_semipositive_definite(matrix_hk, m_facet)) {
            return 0;
        }
    }

    int nn = 0; // Cols of null space matrix
    // temp_1_m = matrix_b -> vector 1 m
    py_null(m_facet, m_level, matrix_dmak, nn, temp_1_m);

//    std::cout << "nn= " << nn << std::endl;

    if (gt(nn, 0)) {
        // We have a nullspace
        // Look at signs of onedim nullspace
        int sum_sinnp = 0; // +1 if b_i is > 0
        int sum_sinnm = 0; // +1 if b_i is < 0
        double sum_b = 0;
        for (int i = 0; i < m_facet; i++) {
            if (ge(temp_1_m[i], 0.0))
                sum_sinnp++;
            else if (le(temp_1_m[i], 0.0))
                sum_sinnm++;

            sum_b += temp_1_m[i];
        }
        bool solution = (sum_sinnp == m_facet || sum_sinnm == m_facet);

//        std::cout << "Sol1= " << solution << std::endl;

        if (solution) {
            // We have an easy to find optimum on S_m
            for (int i = 0; i < m_facet; i++)
                x_star[i] = temp_1_m[i] / sum_b;
        } else if (gt(nn, 1)) {
            // LP
            counters["lp"] += 1;
            py_lp(m_facet, m_facet, m_level, matrix_dmak, x_star, solution);

//            std::cout << "Sol2= " << solution << std::endl;
        }

        if (solution) {
            counters["xAx"] += 1;
            for (int i = 0; i < m_facet; i++) // Reset vector
                temp_1_m[i] = 0.0;
            double f_value = f(x_star, matrix_ak, m_facet, temp_1_m);
            if (lt(f_value, 0.0)) {
                return -1;
            } else {
                if (chol)
                    return 2;

                product_matrix_matrix(matrix_dmak, m_facet, m_facet, matrix_dm, m_facet, m_facet, matrix_hk); // DmAkDm

                double maxsqrdist = 0.0;
                double* eye_m = new double[m_facet*m_facet];
                double sum_col;
                for (int i = 0; i < m_facet; i++) {
                    sum_col = 0.0;
                    for (int j = 0; j < m_facet; j++) {
                        if (i == j)
                            eye_m[i*m_facet+j] = (x_star[j] - 1.0) * (x_star[j] - 1.0);
                        else
                            eye_m[i*m_facet+j] = (x_star[j] - 0.0) * (x_star[j] - 0.0);
                        
                        sum_col += eye_m[i*m_facet+j];
                    }

                    if (gt(sum_col, maxsqrdist))
                        maxsqrdist = sum_col;
                }

                for (int i = 0; i < m_facet; i++) // Reset vector
                    temp_1_m[i] = 0.0;
                double muu1 = min_real_eig(matrix_hk, m_facet, temp_1_m);

                if (ge((f_value - muu1)*maxsqrdist, 0.0))
                    return 0;
            }
        }
    } // gt(nn, 0)

    return 0;
}

/**
 * Facetas copos sorted by level [live]
 * 
 * @param matrix_a: Square matrix A to check copositivity
 * @param n: Matrix A size
 * @return True if A is copositive False if not
 */
std::tuple< bool, std::string, std::map<std::string, int>, std::map<int, int> > facetas_copos_sorted_bylevel_live(double* matrix_a, int n) {
    // Counters
    std::map<std::string, int> counters; // Counter -> Count
    counters["monotonicity"] = 0;
    counters["lp"] = 0;
    counters["matrixpositive"] = 0;
    counters["xAx"] = 0;
    std::map<int, int> chol_evals; // Level -> Evals
    for (int i = 1; i <= n; i++)
        chol_evals[i] = 0;
    int m_level_chol_evals = 0;

    double* temp_1_m;
    temp_1_m = new double[n];
    for (int i = 0; i < n; i++)
        temp_1_m[i] = 0.0;

    // Check centre point is positive
    double centre_value = 1 / n;
    double* centre_point = new double[n];
    for (int i = 0; i < n; i++)
        centre_point[i] = centre_value;
    double centre_f_value = f(centre_point, matrix_a, n, temp_1_m);
    counters["xAx"] += 1;
    delete[] centre_point;
    if (lt(centre_f_value, 0.0))
        return std::make_tuple(false, "A no es copositiva (cAc)", counters, chol_evals);

    // Check main diagonal is positive
    double* diag = matrix_diagonal(matrix_a, n, n);
    double mindiag = *std::min_element(diag, diag+n);
    delete[] diag;
    if (lt(mindiag, 0.0))
        return std::make_tuple(false, "A no es copositiva (DIagonal)", counters, chol_evals);

    // Check matrix is possitive (according accuracy)
//    if (ge_matrix_zero(matrix_a, n, n)) {
//        counters["matrixpositive"] += 1;
//        return std::make_tuple(true, "A es copositiva (1)", counters, chol_evals);
//    }
    
    // Check A is positive semidefinite
    chol_evals[n] += 1;
    if (matrix_semipositive_definite(matrix_a, n))
        return std::make_tuple(true, "A es copositiva (Cholesky)", counters, chol_evals);

    // Binary strings manager
    //Binary binaries(n);

    // Init python
    py_init();

    // D_m matrix manager
    DmMatrix dm_matrices;

    // Initial simplex
    SimplexFacet initial_simplex_n(n);

    chronotime::start_second_time();

    // Level
    std::vector<ULLI> level;
    level.push_back(bin_to_int(initial_simplex_n.get_z()));

    // All SubFacets Done
    std::vector<ULLI> all_subfacets_done;

    // Monos
    std::vector<ULLI> monos;

    // Actual facet
    SimplexFacet* facet_k;
    int m_level = n;

    int m_facet;
    std::string z_facet_k;
    double* matrix_ak;
    double* matrix_hk;
    double* temp_m_m;
    double* matrix_dm;
    double* x_star; // 1 x m

    matrix_ak = new double[m_level*m_level];
    matrix_hk = new double[m_level*m_level];
    temp_m_m = new double[m_level*m_level];
    matrix_dm = NULL;
    x_star = new double[m_level];

    // Check level
    int i = 0;
    int limite = size(level);
    int n_facets_with_min = 0;
    
    int state;
    std::string mono_z;
    int mono_m;
    std::map<int, std::vector<ULLI>> facets_with_min;

    // Stats
    std::cout   << "·" << termcolor::blue << "L: " << m_level << termcolor::reset
                << termcolor::cyan << " Eval: " << limite << termcolor::reset
                << termcolor::magenta << " ASFD: " << all_subfacets_done.size() << termcolor::reset << std::endl;

    std::cout   << termcolor::red << "\t Generation " << termcolor::reset
                << termcolor::bold << termcolor::yellow << "\t T: " << chronotime::format_time(chronotime::final_second_time()) << termcolor::reset << " | "
                << termcolor::bold << chronotime::format_time(chronotime::final_time()) << termcolor::reset << std::endl;
    chronotime::start_second_time();

    while (i < limite) {
        z_facet_k = int_to_bin(level[i], n);
        facet_k = initial_simplex_n.get_facet(z_facet_k, m_level); // level[i] == Fk->get_z()
        m_facet = facet_k->get_m();

        m_level_chol_evals = 0;

        state = eval_facet(matrix_a, n, m_level, dm_matrices,
                            facet_k, m_facet, z_facet_k,
                            matrix_ak, matrix_hk, temp_1_m, temp_m_m, matrix_dm, x_star,
                            counters, m_level_chol_evals,
                            mono_z, mono_m);

        chol_evals[m_level] += m_level_chol_evals;        

//        std::reverse(z_facet_k.begin(), z_facet_k.end());
//        std::cout << termcolor::cyan << bin_to_int(z_facet_k) << " " << state << termcolor::reset << std::endl;

        switch(state) {
            case -1: // f(x) < 0
                return std::make_tuple(false, "A no es copositiva (4)", counters, chol_evals);
                break;
            case 0:
                break;
            case 1: // Mono
                monos.push_back(level[i]);
                facets_with_min[mono_m].push_back(bin_to_int(mono_z));
                break;
            case 2: // All SubFacets Done per level
                all_subfacets_done.push_back(level[i]);
                break;
            default:
                std::cerr << "eval_facet returns " << state << std::endl;
                exit(state);
        }

        delete facet_k;

        i++;
        if (i == limite) { // Next level
            // Stats
            std::cout   << termcolor::green << "\t Complete    " << termcolor::reset
                        << termcolor::bold << termcolor::yellow << "\t T: " << chronotime::format_time(chronotime::final_second_time()) << termcolor::reset << " | "
                        << termcolor::bold << chronotime::format_time(chronotime::final_time()) << termcolor::reset << std::endl;
            chronotime::start_second_time();

            // Remove Mono and ASFD from level
            for (auto f : all_subfacets_done) {
                level.erase(std::remove(level.begin(), level.end(), f), level.end());
            }
            for (auto f : monos) {
                level.erase(std::remove(level.begin(), level.end(), f), level.end());
            }

            // New level
            std::vector<ULLI> next_level;
            std::vector<ULLI> next_all_subfacets_done;

            if (all_subfacets_done.size() == 0) {
                // Full level
                next_level = gen_next_level(level, n);
            } else {
                // Next level from level
                std::vector<std::string> temp;
                for (auto facet : level) { // UNION next_level
                    temp = get_nextzs(int_to_bin(facet, n));
                    //next_level.insert(next_level.end(), temp.begin(), temp.end());
                    for (auto z : temp)
                        next_level.push_back(bin_to_int(z));
                    temp.clear();
                }

                sort(next_level.begin(), next_level.end());
                next_level.erase(unique(next_level.begin(), next_level.end()), next_level.end());

                for (auto facet : all_subfacets_done) {
                    temp = get_nextzs(int_to_bin(facet, n));
                    for (auto z : temp)
                        next_all_subfacets_done.push_back(bin_to_int(z));
                    temp.clear();
                }

                sort(next_all_subfacets_done.begin(), next_all_subfacets_done.end());
                next_all_subfacets_done.erase(unique(next_all_subfacets_done.begin(), next_all_subfacets_done.end()), next_all_subfacets_done.end());
            }

            if (all_subfacets_done.size() > 0) {
                for (auto facet : facets_with_min[m_level-1])
                    next_level.push_back(facet);

                sort(next_level.begin(), next_level.end());
                next_level.erase(unique(next_level.begin(), next_level.end()), next_level.end());

                std::vector<ULLI> intersect;
                std::set_intersection(next_level.begin(), next_level.end(), next_all_subfacets_done.begin(), next_all_subfacets_done.end(), back_inserter(intersect));

                for (auto facet_to_removed : intersect)
                    next_level.erase(std::remove(next_level.begin(), next_level.end(), facet_to_removed), next_level.end());
            }

            level = next_level;
            all_subfacets_done = next_all_subfacets_done;

            // Update counters
            i = 0;
            limite = size(level);
            m_level--;

            facets_with_min.erase(m_level);
            n_facets_with_min = facets_with_min.size();

            if (limite == 0 && n_facets_with_min == 0)
                break;

            // Stats
            std::cout   << "·" << termcolor::blue << "L: " << m_level << termcolor::reset
                        << termcolor::cyan << " Eval: " << limite << termcolor::reset
                        << termcolor::magenta << " ASFD: " << all_subfacets_done.size() << termcolor::reset << std::endl;

            std::cout   << termcolor::red << "\t Generation " << termcolor::reset
                        << termcolor::bold << termcolor::yellow << "\t T: " << chronotime::format_time(chronotime::final_second_time()) << termcolor::reset << " | "
                        << termcolor::bold << chronotime::format_time(chronotime::final_time()) << termcolor::reset << std::endl;
            chronotime::start_second_time();
        }
    }

    // Delete dm
    dm_matrices.free();

    delete[] matrix_ak;
    delete[] matrix_hk;
    delete[] temp_m_m;
    delete[] temp_1_m;
    delete[] x_star;

    py_end();

    return std::make_tuple(true, "A es copositiva (End)", counters, chol_evals);
}

#endif /* FACETS_H_ */
