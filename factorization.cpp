#include <iostream>
#include <string>
#include <cstdlib>
#include <gmpxx.h>
#include <list>
#include <chrono>
#include <iterator>
#include <stdexcept>
#include <sstream>

#include "testzahlen.h"


/**
 * @class SingularEllipticCurveException
 * @brief Custom exception thrown when trying to initialize a singular elliptic curve.
 */
class SingularEllipticCurveException : public std::invalid_argument {
public:
    explicit SingularEllipticCurveException(const mpz_class& a, const mpz_class& b, const mpz_class& n)
        : std::invalid_argument(create_message(a, b, n)) {}

private:
    // Helper function to format the exception message
    static std::string create_message(const mpz_class& a, const mpz_class& b, const mpz_class& n) {
        std::ostringstream oss;
        oss << "Invalid inputs a: " << a << ", b: " << b << ", n: " << n
            << "; discriminant (4 * a^3 + 27 * b^2) mod n is not invertible in Z/nZ.";
        return oss.str();
    }
};


class EllipticCurve {
private:
    mpz_class a, b, n;

public:
    EllipticCurve(const mpz_class& a, const mpz_class& b, const mpz_class& n) : a(a), b(b), n(n) {
        // Calculate discriminant with modular reduction at each step (to  minimizes growth of intermediate results)
        mpz_class discriminant = (((4 * a) % n) * (a % n) * (a % n)) % n;  // (4 * a^3) mod n
        const mpz_class tpm_term = (((27 * b) % n) * (b % n)) % n;        // (27 * b^2) mod n
        discriminant = (discriminant + tpm_term) % n;           // (4 * a^3 + 27 * b^2) mod n

        // Check if discriminant is invertible in Z/nZ
        mpz_class gcd_result;
        mpz_gcd(gcd_result.get_mpz_t(), discriminant.get_mpz_t(), n.get_mpz_t());

        if (gcd_result != 1) {
            throw SingularEllipticCurveException(a, b, n);
        }
    }

    const mpz_class& get_a() const { return a; }
    const mpz_class& get_b() const { return b; }
    const mpz_class& get_modulus() const { return n; }

    /**
     * @brief Overloaded equality operator to compare two elliptic curves.
     * @param other The elliptic curve to compare with.
     * @return True if both curves have the same parameters (a, b, and n), otherwise false.
     */
    bool operator==(const EllipticCurve& other) const {
        return (a == other.a && b == other.b && n == other.n);
    }

    /**
     * @brief Overloaded inequality operator to allow convenient comparison.
     * @param other The elliptic curve to compare with.
     * @return True if curves are not the same.
     */
    bool operator!=(const EllipticCurve& other) const {
        return !(*this == other);
    }
};

// auf der Elliptic curve noch eine Methode implementieren, die prüft, ob ein Punkt auf der Kurve liegt
// dann drei Methoden implementieren (jeweils mit Exception, falls invertierung über Z/nZ nicht durchgeführt werden kann):
//    - add points
//    - point doubling
//    - scalar multiplication
// dann alles testen (auch den == Operator)


/**
 * @class Factor
 * @brief Represents a (prime) factor with an associated exponent.
 *
 * Stores information about a number, including whether it is prime,
 * its exponent (used for representing powers such as 2^6), and
 * the factor itself (preferably a prime number).
 */
class Factor {
public:
    unsigned int exponent;  ///< Exponent to represent powers in factorization (e.g., 2^6)
    mpz_class factor; ///< Factor, ideally a prime number
    int is_prime; ///< Indicates if 'factor' is prime (2 = definitely prime, 1 = probably prime, 0 = not prime)

    /**
     * @brief Default constructor initializes the factor with exponent 0, factor 1, and non-prime status.
     */
    Factor()
    {   exponent = 0;
        is_prime = 0;
        factor = mpz_class("1");
    }

    /**
     * @brief Prints the factor in a specific format based on its properties.
     *
     * Outputs the factor in a different format based on whether it is prime
     * or non-prime. Probable primes are suffixed with '_?'.
     */
    void printpp()
    {
        char open = '(';
        char close = ')';

        if(is_prime == 0) // Non-prime factor
        {
            open = '[';
            close = ']';
        }


        std::cout << open << factor;

        if(exponent > 1)
            std::cout << "^" << exponent;
        if(is_prime == 1) // Probable prime
            std::cout << "_?";

        std::cout << close << " ";
    }
};

/**
 * @brief Divides out the maximum power of a (prime) number P from T.
 *
 * Takes a number T and a (prime) number P, returning the maximum exponent `e`
 * such that P^e divides T. After calling this function, T is updated to T / P^e.
 *
 * @param t The number to be divided (will be modified).
 * @param P The prime number divisor.
 * @return unsigned int The maximum exponent `e` such that P^e divides T.
 */
unsigned int divide_out_maximal_power(mpz_class& t, mpz_class P)
{
    unsigned int exponent = 0;
    while(t % P == 0)
    {
        t = t / P;
        exponent++;
    }
    return exponent;
}


/**
 * @brief Factorizes N (and modifies it) using trial division up to a bound B.
 *
 * Factorizes the number N by trial division up to a specified bound B
 * and checks if the remainder is a probable prime. Returns a list of found
 * prime factors. (For large values of N, selecting a large B can significantly
 * increase runtime.)
 *
 * @param N The number to factorize (will be modified).
 * @param B The upper bound for trial division.
 * @return std::list<Factor> A list of found factors with their exponents.
 */
std::list<Factor> trial_division_bounded(mpz_class& N, mpz_class B) {
    mpz_class C = sqrt(N)+1;
    if (B > C) B = C;

    std::list<Factor> factors;
    Factor ppt;

    // Check if N is prime (27 tests: 1 Baillie-PSW and 2 Rabin-Miller tests)
    int a = mpz_probab_prime_p(N.get_mpz_t(), 27);
    if( a > 0) // N is probably prime
    {
        ppt.exponent = 1;
        ppt.factor = N;
        ppt.is_prime = a;
        factors.push_back(ppt);
        N = mpz_class("1");
        return factors;
    }

    // Factorize N by trial division
    // First try 2 and 3, then numbers congruent to 1 or 5 modulo 6
    mpz_class P("2");
    ppt.exponent = divide_out_maximal_power(N,P);
    if(ppt.exponent > 0)
    {
        ppt.factor = P;
        ppt.is_prime = 2;
        factors.push_back(ppt);

    }
    P++; // P = 3
    ppt.exponent = divide_out_maximal_power(N,P);
    if(ppt.exponent > 0)
    {
        ppt.factor = P;
        ppt.is_prime = 2;
        factors.push_back(ppt);

    }

    P+=3; // P = 6
    while(N > 1 and P <= B )
    {
        ppt.exponent = divide_out_maximal_power(N,P-1);
        if(ppt.exponent > 0)
        {
            ppt.factor = P-1;
            ppt.is_prime = 2;
            factors.push_back(ppt);

        }
        ppt.exponent = divide_out_maximal_power(N,P+1);
        if(ppt.exponent > 0)
        {
            ppt.factor = P+1;
            ppt.is_prime = 2;
            factors.push_back(ppt);
        }
        P += 6;
    }

    if( N > 1)
    {
        a = mpz_probab_prime_p(N.get_mpz_t(), 27);
        if( a == 0) // Remaining number is probably prime
        {
            ppt.exponent = 1;
            ppt.factor = N;
            ppt.is_prime = a;
            factors.push_back(ppt);
            N = mpz_class("1");
            return factors;
        }
    }

    return factors;
}


/**
 * @class UnknownNumberModeException
 * @brief Custom exception thrown when an invalid number mode is provided.
 */
class UnknownNumberModeException : public std::invalid_argument {
public:
    explicit UnknownNumberModeException(const std::string& message)
        : std::invalid_argument(message) {}
};


/**
 * @brief Helper function to check if a string represents an integer.
 *
 * @param s Input string.
 * @return true if the string is an integer; otherwise, false.
 */
bool is_integer(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
};


/**
 * @enum NumberMode
 * @brief Modes for generating specific numbers to be factorized.
 */
enum NumberMode : int {
    DirectNumber = 0,
    FermatNumber = 1,
    CunninghamNumber = 2,
    TestNumber = 3,
    RSANumber = 4,
    NotImplemented = 5
};


/**
 * @brief Generates numbers for factorization based on user input and mode.
 *
 * Converts a user-provided input (string and mode) into a number to be factorized.
 * Function implementations are in `testzahlen.h` and `testzahlen.cpp`.
 *
 * @param number_as_string The number provided by the user.
 * @param number_mode Mode to interpret the input number.
 * @return mpz_class The generated number.
 */
mpz_class get_number_from_user_input(const std::string& number_as_string, const NumberMode number_mode) {
    switch (number_mode) {
        case DirectNumber: return mpz_class(number_as_string);
        case FermatNumber: return Fermat(std::stoi(number_as_string));
        case CunninghamNumber: return Cunningham(std::stoi(number_as_string));
        case TestNumber: return TestzahlB(std::stoi(number_as_string));
        case RSANumber: return RSAZahl(std::stoi(number_as_string));
        default: throw UnknownNumberModeException("Unknown mode " + std::to_string(number_mode) + " provided.");
    }
}


/**
 * @brief Main function for the factorization program.
 *
 * Accepts one of the following parameters to specify the number mode:
 *   - `-F` for Fermat numbers
 *   - `-C` for Cunningham numbers
 *   - `-B` for (custom) Test numbers
 *   - `-R` for RSA challenge numbers
 *
 * If no parameter is provided, prompts the user to input a natural number.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return int Program exit status.
 */
int main(int argc, char *argv[]) {
    NumberMode mode = DirectNumber;
    std::string question ("natural number k");

    if(argc == 2)  // one parameter provided
    {
        std::string param (argv[1]);
        if(param == "-F")
        {
            mode = FermatNumber;
            question = "Fermat-Number F_k";
        }
        else if(param == "-C")
        {
            mode = CunninghamNumber;
            question = "Cunningham-Number C_k";
        }
        else if(param == "-B")
        {
            mode = TestNumber;
            question = "Test number B_k";
        }
        else if(param == "-R")
        {
            mode = RSANumber;
            question = "RSA-Number R_k";
        }
        else {  // Invalid parameter provided
            std::cerr << "Error: Unknown Parameter " << param << std::endl;
            return EXIT_FAILURE; // Return non-zero status to indicate an error
        }
    }

    if (argc > 2)  // Too many arguments
    {
        std::cerr << "Error: Too many parameters provided." << std::endl;
        return EXIT_FAILURE;
    }

    // Get user input
    std::cout << "Which " << question << " should be factored?" << std::endl;
    std::cout  << "Input k:  " ;
    std::string number_as_string;
    std::cin >> number_as_string;
    if (!is_integer(number_as_string)) {
        std::cerr << "Error: Input must be an integer." << std::endl;
        return EXIT_FAILURE;
    }

    // Generate the number to be factored based on user input
    mpz_class N;
    try {
        N = get_number_from_user_input(number_as_string, mode);
    } catch (const UnknownNumberModeException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Valid input received. The program will attempt to factor: " << N << std::endl;


    // Factorize up to the limit (trial division up to 10^8, adjustable if necessary)
    mpz_class B("100000000");  // via Variable aus Header file mitgeben
    std::list<Factor> factors = trial_division_bounded(N,B);

    std::cout << "Found the following prime factors:" << std::endl;
    for (std::list<Factor>::iterator it = factors.begin(); it != factors.end(); it++)
        it->printpp();

    std::cout << "Remainder of N: " << N << std::endl;
    // Now use elliptic curve algorithm for remainder
    std::string a_as_string;
    std::string b_as_string;
    std::string n_as_string;

    std::cout << "a: " << std::endl;
    std::cin >> a_as_string;
    mpz_class a = mpz_class(a_as_string);

    std::cout << "b: "<< std::endl;
    std::cin >> b_as_string;
    mpz_class b = mpz_class(b_as_string);

    std::cout << "n: " << std::endl;
    std::cin >> n_as_string;
    mpz_class n = mpz_class(n_as_string);

    std::cout << "Try to instantiate elliptic curve" << std::endl;
    EllipticCurve x = EllipticCurve(a,b,n);
    std::cout << "Done!" << std::endl;

    return EXIT_SUCCESS; // Program completed successfully
}