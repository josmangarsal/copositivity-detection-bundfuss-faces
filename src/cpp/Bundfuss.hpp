/**
    Bundfuss
    Bundfuss.hpp
    Purpose: Bundfuss refinement

    @author J.M.G. Salmerón
    @version 1.0 11/10/2018
        
    Created on: 11/10/2018
*/

// Libs
#include <iostream> // cout
#include "termcolor.hpp" // Output color

// Headers
#include "SimplexBundfuss.h"
#include "Matrices.h"
#include "Binary.h"
#include "ChronoTime.hpp"

#ifndef BUNDFUSS_H_
#define BUNDFUSS_H_

/**
 * Bundfuss bisection
 * 
 * @param matrix_a: Square matrix A to check copositivity
 * @param n: Matrix A size
 * @param eps: Accuracy
 * @return True if A is copositive False if not
 */
std::tuple< bool, std::string, std::map<std::string, int> > bundfuss_bisection(double* matrix_a, int n, double eps) {
    // Temporal arrays
    double* temp_1_n = new double[n];

    // Counters
    std::map<std::string, int> counters; // Counter -> Count
    counters["generated_simplices"] = 0;
    counters["evaluated_simplices"] = 0;
    counters["divided_simplices"] = 0;
    counters["final_simplices"] = 0;
    //counters["max_tree_level"] = 0;
    counters["monotonicity"] = 0;
    counters["Q>-eps"] = 0;
    counters["xAx"] = 0;

    // Initial simplex
    SimplexBundfuss* s0 = new SimplexBundfuss(n);
    counters["generated_simplices"] += 1;
    //counters["max_tree_level"] = 0;

    // Evaluate initial simplex
    counters["evaluated_simplices"] += 1;
    try {
        if (s0->is_copos(matrix_a, eps, temp_1_n, counters))
            counters["final_simplices"] += 1;
    } catch(BundfussException &e) {
        return std::make_tuple(false, e.what(), counters);
    }

    std::vector<SimplexBundfuss*> working_list;

    // Check initial simplex
    if (lt(s0->get_copos(), eps))
        working_list.push_back(s0);
    else
        counters["Q>-eps"] += 1;

    SimplexBundfuss* working_simplex;
    while (!working_list.empty()) {
        // Get simplex
        working_simplex = working_list.back(); // Depth search
        working_list.pop_back();

        // Divide
        auto [s1, s2] = working_simplex->divide();
        counters["divided_simplices"] += 1;
        counters["generated_simplices"] += 2;
        //counters["max_tree_level"] = s1->get_level();

        delete working_simplex;

        // Evaluate child 1
        counters["evaluated_simplices"] += 1;
        try {
            if (s1->is_copos(matrix_a, eps, temp_1_n, counters))
                counters["final_simplices"] += 1;
        } catch(BundfussException &e) {
            return std::make_tuple(false, e.what(), counters);
        }

        // Add child 1
        if (lt(s1->get_copos(), eps))
            working_list.push_back(s1);
        else
            counters["Q>-eps"] += 1;

        // Evaluate child 2
        counters["evaluated_simplices"] += 1;
        try {
            if (s2->is_copos(matrix_a, eps, temp_1_n, counters))
                counters["final_simplices"] += 1;
        } catch(BundfussException &e) {
            return std::make_tuple(false, e.what(), counters);
        }

        // Add child 2
        if (lt(s2->get_copos(), eps))
            working_list.push_back(s2);
        else
            counters["Q>-eps"] += 1;
    }

    // Free temp memory
    delete[] temp_1_n;

    return std::make_tuple(true, "A es copositiva (End)", counters);
}

/**
 * Bundfuss bisection with monotonicity
 * 
 * @param matrix_a: Square matrix A to check copositivity
 * @param n: Matrix A size
 * @param eps: Accuracy
 * @return True if A is copositive False if not
 */
std::tuple< bool, std::string, std::map<std::string, int> > bundfuss_bisection_mono(double* matrix_a, int n, double eps) {
    // Temporal arrays
    double* temp_1_n = new double[n];

    // Counters
    std::map<std::string, int> counters; // Counter -> Count
    counters["generated_simplices"] = 0;
    counters["evaluated_simplices"] = 0;
    counters["divided_simplices"] = 0;
    counters["final_simplices"] = 0;
    //counters["max_tree_level"] = 0;
    counters["monotonicity"] = 0;
    counters["Q>-eps"] = 0;
    counters["xAx"] = 0;

    // Initial simplex
    SimplexBundfuss* s0 = new SimplexBundfuss(n);
    counters["generated_simplices"] += 1;
    //counters["max_tree_level"] = 0;

    // Evaluate initial simplex
    counters["evaluated_simplices"] += 1;
    try {
        if (s0->is_copos(matrix_a, eps, temp_1_n, counters))
            counters["final_simplices"] += 1;
    } catch(BundfussException &e) {
        return std::make_tuple(false, e.what(), counters);
    }

    std::vector<SimplexBundfuss*> working_list;

    // Check initial simplex
    if (lt(s0->get_copos(), eps))
        working_list.push_back(s0);
    else
        counters["Q>-eps"] += 1;

    SimplexBundfuss* working_simplex;
    while (!working_list.empty()) {
        // Get simplex
        working_simplex = working_list.back(); // Depth search
        working_list.pop_back();

        // Divide mono
        auto [new_z, new_m] = working_simplex->monotonicity(matrix_a, temp_1_n);
        if (new_m < working_simplex->get_m()) {
            counters["monotonicity"] += 1;

            SimplexBundfuss* reduced_simplex = working_simplex->get_hat(new_z, new_m);

            // Evaluate reduced
            counters["evaluated_simplices"] += 1;
            try {
                if (reduced_simplex->is_copos(matrix_a, eps, temp_1_n, counters))
                    counters["final_simplices"] += 1;
            } catch(BundfussException &e) {
                return std::make_tuple(false, e.what(), counters);
            }

            // Divided reduced
            if (lt(reduced_simplex->get_copos(), eps)) {
                auto [s1, s2] = reduced_simplex->divide();
                counters["divided_simplices"] += 1;
                counters["generated_simplices"] += 2;
                //counters["max_tree_level"] = s1->get_level();

                delete reduced_simplex;
                delete working_simplex;

                // Evaluate child 1
                counters["evaluated_simplices"] += 1;
                try {
                    if (s1->is_copos(matrix_a, eps, temp_1_n, counters))
                        counters["final_simplices"] += 1;
                } catch(BundfussException &e) {
                    return std::make_tuple(false, e.what(), counters);
                }

                // Add child 1
                if (lt(s1->get_copos(), eps))
                    working_list.push_back(s1);
                else
                    counters["Q>-eps"] += 1;

                // Evaluate child 2
                counters["evaluated_simplices"] += 1;
                try {
                    if (s2->is_copos(matrix_a, eps, temp_1_n, counters))
                        counters["final_simplices"] += 1;
                } catch(BundfussException &e) {
                    return std::make_tuple(false, e.what(), counters);
                }

                // Add child 2
                if (lt(s2->get_copos(), eps))
                    working_list.push_back(s2);
                else
                    counters["Q>-eps"] += 1;
            } else {
                counters["Q>-eps"] += 1;
                
                delete reduced_simplex;
                delete working_simplex;
            }
        } else {
            auto [s1, s2] = working_simplex->divide();
            counters["divided_simplices"] += 1;
            counters["generated_simplices"] += 2;
            //counters["max_tree_level"] = s1->get_level();

            delete working_simplex;

            // Evaluate child 1
            counters["evaluated_simplices"] += 1;
            try {
                if (s1->is_copos(matrix_a, eps, temp_1_n, counters))
                    counters["final_simplices"] += 1;
            } catch(BundfussException &e) {
                return std::make_tuple(false, e.what(), counters);
            }

            // Add child 1
            if (lt(s1->get_copos(), eps))
                working_list.push_back(s1);
            else
                counters["Q>-eps"] += 1;

            // Evaluate child 2
            counters["evaluated_simplices"] += 1;
            try {
                if (s2->is_copos(matrix_a, eps, temp_1_n, counters))
                    counters["final_simplices"] += 1;
            } catch(BundfussException &e) {
                return std::make_tuple(false, e.what(), counters);
            }

            // Add child 2
            if (lt(s2->get_copos(), eps))
                working_list.push_back(s2);
            else
                counters["Q>-eps"] += 1;
        }
    }

    // Free temp memory
    delete[] temp_1_n;

    return std::make_tuple(true, "A es copositiva (End)", counters);
}

#endif /* BUNDFUSS_H_ */