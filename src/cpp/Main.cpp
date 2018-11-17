/**
    Main
    Main.cpp
    Purpose: Main for Copositivity detection (C++17 version)

    @author J.M.G. Salmerón
    @version 1.1 11/10/2018
        
    Created on: 18/09/2018
*/

#include <fstream> // Input
#include "clara.hpp" // Argument
#include "Facets.hpp"
#include "Bundfuss.hpp"

// Load matrix clique
double* load_matrix_clique(std::string matrix_a_path, int t, int n, auto parser) {
    double* matrix_a_original = new double[n*n];

    std::ifstream file(matrix_a_path);
    if(!file)
    {
        std::cerr << "Error in command line: " << "Error opening matrix file" << std::endl;
        std::cerr << std::endl;
        parser.writeToStream(std::cout);
        exit(1);
    }

    int i = 0, j = 0;
    std::string elem;
    while(file >> elem) {
        matrix_a_original[i*n+j] = std::stod(elem);
        j++;
        if (j == n) {
            j = 0;
            i++;
            if (i == n)
                i = 0;
        }
    }

    /// E = np.ones((n,n))
    double* matrix_e = new double[n*n];
    for (int i = 0; i < n*n; i++)
        matrix_e[i] = 1;
    /// A1 = (t - 1) * E
    double* matrix_a1 = product_scalar_matrix(t-1, matrix_e, n, n);
    /// A2 = t * A
    double* matrix_a2 = product_scalar_matrix(t, matrix_a_original, n, n);
    /// A = A1 - A2
    double* matrix_a = matrix_minus_matrix(matrix_a1, matrix_a2, n, n);
    /// Free temporal matrices
    delete[] matrix_a_original;
    delete[] matrix_e;
    delete[] matrix_a1;
    delete[] matrix_a2;

    return matrix_a;
}
///////////////////////////////////////////////////////////

// Main
int main(int argc, char* argv[]) {
    chronotime::start_time();

    // Error file
    std::ofstream err("error.txt");
    std::cerr.rdbuf(err.rdbuf());
    std::cerr << "Error log" << std::endl;

    // Arguments
    int n = 0;
    int t = 0;
    std::string division_method = ".";
    std::string matrix_a_path = ".";
    
    bool showhelp = false;

    /// Arg parser
    auto parser =
        clara::Opt( [&n] (int lambda_n) {
            if (lambda_n <= 0) {
                return clara::ParserResult::runtimeError("Dimension must be greater or equal than 1.");
            } else {
                n = lambda_n;
                return clara::ParserResult::ok(clara::ParseResultType::Matched);
            }
        }, "dimension")
        ["-n"]["--dimension"]
        ("Problem dimension (integer) [Required]") |
        clara::Opt(matrix_a_path, "matrix_path")["-M"]["--matrix"]("Matrix file (string) [n!=3,5]") |
        clara::Opt(t, "dimacs_t")["-t"]("Dimacs t value (integer) [n!=3,5]") |
        clara::Opt(division_method, "division_method")["-D"]["--division"]("Division method (string) ['facet', 'bundfuss', 'zbund'] [Required]") |
        clara::Help(showhelp);

    try {
        auto result = parser.parse(clara::Args(argc, argv));
        if (!result) {
            std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
            std::cerr << std::endl;
            parser.writeToStream(std::cout);
            exit(1);
        } else {
            std::cout <<
            "D: " << division_method << std::endl <<
            "n: " << n << std::endl <<
            "M: " << matrix_a_path << std::endl <<
            "t: " << t << std::endl << std::endl;
        }
    } catch (std::exception const & e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        exit(1);
    }

    // Load matrix
    double* matrix_a;
    switch (n) {
        case 3:
            matrix_a = new double[n*n];
            matrix_a[0*n+0] =  2.0; matrix_a[0*n+1] = -1.0; matrix_a[0*n+2] = -1.0;
            matrix_a[1*n+0] = -1.0; matrix_a[1*n+1] =  2.0; matrix_a[1*n+2] = -1.0;
            matrix_a[2*n+0] = -1.0; matrix_a[2*n+1] = -1.0; matrix_a[2*n+2] =  2.0;
            break;
        case 5:
            matrix_a = new double[n*n];
            matrix_a[0*n+0] =  1.0; matrix_a[0*n+1] = -1.0; matrix_a[0*n+2] =  1.0; matrix_a[0*n+3] =  1.0; matrix_a[0*n+4] = -1.0;
            matrix_a[1*n+0] = -1.0; matrix_a[1*n+1] =  1.0; matrix_a[1*n+2] = -1.0; matrix_a[1*n+3] =  1.0; matrix_a[1*n+4] =  1.0;
            matrix_a[2*n+0] =  1.0; matrix_a[2*n+1] = -1.0; matrix_a[2*n+2] =  1.0; matrix_a[2*n+3] = -1.0; matrix_a[2*n+4] =  1.0;
            matrix_a[3*n+0] =  1.0; matrix_a[3*n+1] =  1.0; matrix_a[3*n+2] = -1.0; matrix_a[3*n+3] =  1.0; matrix_a[3*n+4] = -1.0;
            matrix_a[4*n+0] = -1.0; matrix_a[4*n+1] =  1.0; matrix_a[4*n+2] =  1.0; matrix_a[4*n+3] = -1.0; matrix_a[4*n+4] =  1.0;
            break;
        default:
            matrix_a = load_matrix_clique(matrix_a_path, t, n, parser);
    }

    // Let's go
    if (division_method == "facet") {
        auto [ copos, result, counters, chol_evals ] = facetas_copos_sorted_bylevel_live(matrix_a, n);
        
        std::cout   << std::endl << /*termcolor::blink <<*/ termcolor::bold
                    << (copos ? termcolor::green : termcolor::red) << result
                    << termcolor::reset << std::endl;
        
        for ( auto [key, value] : counters)
            std::cout << "\t" << key << ": " << value << std::endl;
        
        std::cout << "\tchol: ";
        for ( auto [key, value] : chol_evals) // key is a must here, sorry for the warning
            std::cout << value << " ";
        std::cout << std::endl;

    } else if (division_method == "bundfuss") {
        auto [ copos, result, counters ] = bundfuss_bisection(matrix_a, n, -accuracy);

        std::cout   << std::endl << /*termcolor::blink <<*/ termcolor::bold
                    << (copos ? termcolor::green : termcolor::red) << result
                    << termcolor::reset << std::endl;
        
        for ( auto [key, value] : counters)
            std::cout << "\t" << key << ": " << value << std::endl;

    } else if (division_method == "zbund") {
        auto [ copos, result, counters ] = bundfuss_bisection_mono(matrix_a, n, -accuracy);

        std::cout   << std::endl << /*termcolor::blink <<*/ termcolor::bold
                    << (copos ? termcolor::green : termcolor::red) << result
                    << termcolor::reset << std::endl;
        
        for ( auto [key, value] : counters)
            std::cout << "\t" << key << ": " << value << std::endl;

    } else {
        std::cout << termcolor::on_red << termcolor::white << "Unknown division method" << termcolor::reset << std::endl;
        parser.writeToStream(std::cout);
    }

    std::cout   << std::endl << "Total_time: "
                << termcolor::bold << chronotime::format_time(chronotime::final_time()) << " seconds" << termcolor::reset
                << std::endl;

    return 0;
}
