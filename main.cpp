#include <iostream>
#include <string>
#include <cstdlib>
#include <gmpxx.h>
#include <list>
#include <chrono>
#include <stdexcept>

#include "trial_division.h"
#include "helper.h"
#include "elliptic_curve.h"


/**
 * @brief Main function for the factorization program.
 *
 * This program factors a number based on user input and command-line parameters. The number
 * can represent a Fermat, Cunningham, RSA challenge number, or a custom (test) number. The
 * factorization can be customized using several optional parameters, including bounds for
 * elliptic curve factorization and toggling trial division.
 *
 * ## Command-line Parameters
 *
 * The program accepts the following optional parameters:
 *
 * 1. `--num_mode` or `-m` (string, optional):
 *    - Specifies the type of number to be factored.
 *    - Supported values:
 *      - `Fermat`: Factors a Fermat number, F_k.
 *      - `Cunningham`: Factors a Cunningham number, C_k.
 *      - `Test`: Factors a custom test number, B_k.
 *      - `RSA`: Factors an RSA challenge number, R_k.
 *    - If not provided, the program defaults to "DirectNumber" mode, that is,
 *      the user will be prompted to directly input the number to be factored.
 *
 * 2. `--stage1_bound` or `-b` (positive integer, optional):
 *    - Specifies the stage-1 bound (B), that is, the bound for prime bases, for elliptic curve factorization.
 *    - Defaults to `DEFAULT_B`.
 *
 * 3. `--stage2_bound` or `-c` (positive integer, optional):
 *    - Specifies the stage-2 bound (C), that is, the bound for the smallest prime factor of N, for elliptic curve factorization.
 *    - Defaults to `DEFAULT_C`.
 *
 * 4. `--num_curves` or `-n` (positive integer, optional):
 *    - Specifies the number of elliptic curves used for factorization.
 *    - Defaults to `DEFAULT_NUM_CURVES`.
 *
 * 5. `--no_trial_division` or `-nt` (boolean flag, optional):
 *    - Disables trial division as a preliminary factorization step.
 *    - Defaults to false (trial division is enabled by default).
 *
 * ## Input Handling
 *
 * - If `--num_mode` or `-m` is provided, the program will ask the user for a specific k-value
 *   corresponding to the selected mode.
 * - If no mode is specified, the program will prompt the user to enter a natural number to be
 *   factored directly.
 * - The input k-value or number must be a positive integer.
 *
 * ## Execution Workflow
 *
 * 1. Parse and validate the command-line arguments.
 * 2. Map the `--num_mode` argument to a `NumberMode` enumeration.
 * 3. Prompt the user to input the number or k-value to factor.
 * 4. Generate the number to be factored based on the input and mode.
 * 5. Perform factorization:
 *    - If trial division is enabled, attempt to factorize the number using trial division up
 *      to a predefined bound (`TRIAL_DIV_BOUND`).
 *    - Use elliptic curve factorization for the remaining composite number, based on the
 *      stage-1 and stage-2 bounds and the specified number of curves.
 * 6. Print the elapsed time for the factorization process as well as the prime factors.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return int Program exit status.
 */
int main(int argc, char *argv[]) {
    // Argument defaults
    std::string mode_param;
    mpz_class B = DEFAULT_B;
    mpz_class C = DEFAULT_C;
    int num_curves = DEFAULT_NUM_CURVES;
    bool no_trial_division = DEFAULT_NO_TRIAL_DIVISION;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--num_mode" || arg == "-m") {
            if (i + 1 < argc) mode_param = argv[++i];  // assign mode_param and increment i
            else {
                std::cerr << "Error: --num_mode/-m requires a value." << std::endl;
                return EXIT_FAILURE;
            }
        } else if (arg == "--stage1_bound" || arg == "-b") {
            if (i + 1 < argc && is_positive_integer(argv[i + 1])) B = std::stoi(argv[++i]);
            else {
                std::cerr << "Error: --stage1_bound/-b requires a positive integer value." << std::endl;
                return EXIT_FAILURE;
            }
        } else if (arg == "--stage2_bound" || arg == "-c") {
            if (i + 1 < argc && is_positive_integer(argv[i + 1])) C = std::stoi(argv[++i]);
            else {
                std::cerr << "Error: --stage2_bound/-c requires a positive integer value." << std::endl;
                return EXIT_FAILURE;
            }
        } else if (arg == "--num_curves" || arg == "-n") {
            if (i + 1 < argc && is_positive_integer(argv[i + 1])) num_curves = std::stoi(argv[++i]);
            else {
                std::cerr << "Error: --num_curves/-n requires a positive integer value." << std::endl;
                return EXIT_FAILURE;
            }
        } else if (arg == "--no_trial_division" || arg == "-nt") {
            no_trial_division = true; // Boolean flag, no value expected
        } else {
            std::cerr << "Error: Unknown parameter " << arg << std::endl;
            return EXIT_FAILURE;
        }
    }

    // Map mode_param to NumberMode
    NumberMode mode = DirectNumber;  // default (in case no `--num_mode` provided)
    std::string question = "natural number k";
    if (!mode_param.empty()) {
        if (mode_param == "Fermat") {
            mode = FermatNumber;
            question = "Fermat-Number F_k";
        } else if (mode_param == "Cunningham") {
            mode = CunninghamNumber;
            question = "Cunningham-Number C_k";
        } else if (mode_param == "Test") {
            mode = TestNumber;
            question = "Test number B_k";
        } else if (mode_param == "RSA") {
            mode = RSANumber;
            question = "RSA-Number R_k";
        } else {
            std::cerr << "Error: Invalid --num_mode value: " << mode_param << std::endl;
            return EXIT_FAILURE;
        }
    }

    // Get user input
    std::cout << "Which " << question << " should be factored?" << std::endl;
    std::cout  << "Input k:  " ;
    std::string number_as_string;
    std::cin >> number_as_string;
    if (!is_positive_integer(number_as_string)) {
        std::cerr << "Error: Input must be a positive integer." << std::endl;
        return EXIT_FAILURE;
    }

    // Generate the number to be factored based on user input
    mpz_class N;
    try {
        N = generate_number_from_user_input(number_as_string, mode);
    } catch (const UnknownNumberModeException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    if (N <= 1) {
        std::cerr << "Error: Please provide an integer > 1." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Valid input received. The program will attempt to factor: " << N.get_str() << std::endl;
    std::list<Factor> factors;

    // Start timing ---------
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();

    // Trial division
    if (not no_trial_division) {  // trial division enabled
        mpz_class trial_division_bound(TRIAL_DIV_BOUND);
        factors = trial_division_bounded(N,trial_division_bound);
    }

    if (N > 1) {
        // Generate starting values for B (stage-1 bound) and C (stage-2 bound) for elliptic curve factorization
        mpz_class sqrt_N;
        mpz_sqrt(sqrt_N.get_mpz_t(), N.get_mpz_t());
        if (C == DEFAULT_C || C > (sqrt_N + 1)) {  // None or unreasonable value for C provided by user
            C = get_smallest_prime_bound(N, mode==RSANumber);
        }
        if (B == DEFAULT_B) {
            B = calculate_base_prime_bound_from_smallest_prime_bound(C);
        }

        //ToDo: Take out
        //////////////////// TEST ///////////////
        std::cout << "N: " << N << std::endl;
        std::cout << "C: " << C << std::endl;
        std::cout << "B: " << B << std::endl;
        std::cout << "num_curves: " << num_curves << std::endl;
        std::cout << "no_trial_division: " << no_trial_division << std::endl;
        //////////////////// TEST ///////////////

        factorize_with_elliptic_curves(N, B, C, num_curves, factors);
    }

    // Stop timing ----------
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    // Calculate elapsed time and print it
    std::chrono::duration<double, std::milli> float_ms = end - start;
    std::cout << "Factorization in " << float_ms.count() << " milliseconds completed" << std::endl;

    std::cout << "Found Factorization:\n";
    std::list<Factor>::iterator it;
    for (it = factors.begin(); it != factors.end(); it++)
        it->printpp();

    return EXIT_SUCCESS; // Program completed successfully
}

// Getestete Funktionen:
// - get_smallest_base_biggest_exponent_for_perfect_power
// - calculate_prime_bound_from_smallest_prime_bound


// ToDo: brauche noch Funktion, die mit gegebenen N einen sinnvollen Startwert für B und C
//  evtl. Idee: lasse mir C von "außen" mitgeben und berechne von dort B (mit Formel (14) - siehe S. 566)
// ToDo: dann alles testen was elliptic Curve class angeht!
// ToDo: next steps:
//  - Recherchieren, überlegen, wie ich C (evtl. aus N) berechne...
//  - Logik Schritt für Schritt durchgehen und prüfen
//  - Verbesserungsvorschläge und prüfen lassen
//  - Testen aller einzelnen funktionen
//  - Testen des Gesamtkonstruktes
//  - Faktorisieren der Testzahlen
// Use elliptic curves factorization for remainder
// ToDo: Funktionen schreiben: 1) get_C_from_N (das C als Funktion von N berechnet - siehe Maple)
// ToDo: dann diesen Wert nutzen, um damit B zu berechnen (diese Funktion gibt es schon).
// ToDo: diese Funktionen (falls B und C nicht mitgebenen wurden; d.h. 0 sind) hier anwenden und weitergeben...
//  Funktionen auch innerhalb des algorithmus für elliptische Kurven nutzen!
//  siehe auch was Prof. Kionke geschrieben hat (in Email)



// Info: wie maple B und C berechnet: https://www.maplesoft.com/applications/Preview.aspx?id=3528
//  - if nargs = 2 then B := boundB; else B := 2*length(n); fi;
//  - if nargs = 3 then C := boundC; else C := ilog[2](n); fi;


// Potential Improvements:
// - Use faster way to identify all primes less or equal to bound B (e.g., sieve of Atkins), das muss aber abgeglichen
//   werden, da wenn ich sieb verwende, dann berechne ich wirklich alle, ohne dass ich u.U. wirklich alle benötige
// - Efficiency can be improved at several points (but often at the cost of readability). U.a.:
//      - Berechnen der Primzahlen <= B
//      - Handling der perfect powers
// - Könnte multiprocessing nutzen, um parallel mehrere durchläufe für Lenstra zu machen
// - Siehe Cohen S489: Montgomerys Methode (Parallel Inverse Modulo N)
// Muss wenn ich durchlaufe mal Memory checken (ob der voll läuft)

