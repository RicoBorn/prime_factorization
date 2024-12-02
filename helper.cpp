#include <iostream>
#include <string>
#include <list>
#include <stdexcept>
#include <random>
#include <cmath>
#include <gmpxx.h>

#include "helper.h"
#include "test_numbers.h"


/**
 * @brief Helper function to check if a string represents a positive integer (including 0).
 *
 * @param s Input string.
 * @return true if the string is a positive integer (including 0); otherwise, false.
 */
bool is_positive_integer(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}


/**
 * @brief Computes the representation of d as an element in the ring Z/nZ.
 *
 * @param d The integer to be reduced modulo n.
 * @param n The modulus of the ring Z/nZ. Must be greater than 0.
 * @return The canonical representative of d in Z/nZ, satisfying 0 <= result < n.
 * @throws std::invalid_argument if n < 1.
 */
mpz_class mod_ring(const mpz_class& d, const mpz_class& n) {
    if (n < 1) {
        throw std::invalid_argument("Modulus n must be greater than 0.");
    }
    mpz_class result = d % n;
    if (result < 0) {
        result += n;
    }
    return result;
}


/**
 * @brief Calculates the number of decimal digits in a given integer.
 *
 * This function determines the total number of decimal (base-10) digits in the
 * absolute value of the input integer `N`. Special case: if `N` is 0, the
 * function correctly returns 1, as 0 has one decimal digit.
 *
 * @param N An integer of type `mpz_class` whose decimal digit count is to be determined.
 * @return size_t The number of decimal digits in the absolute value of `N`.
 */
size_t get_number_of_decimal_digits(const mpz_class& N) {
    if (N == 0) {
        return 1; // Special case: 0 has 1 decimal digit
    }
    return mpz_sizeinbase(N.get_mpz_t(), 10);
}


/**
 * @brief Divides out the maximum power of a number P from T. The latter will be updated by the function.
 *
 * Takes a number T and a number P, returning the maximum exponent `e`
 * such that P^e divides T. After calling this function, T is updated to T / P^e.
 *
 * @param T The number to be divided (will be modified).
 * @param P The prime number divisor.
 * @return unsigned int The maximum exponent `e` such that P^e divides T.
 */
unsigned int divide_out_maximal_power(mpz_class& T, const mpz_class& P)
{
    unsigned int exponent = 0;
    while(T % P == 0)
    {
        T = T / P;
        exponent++;
    }
    return exponent;
}


/**
 * @brief Default constructor initializes the factor with exponent 0, factor 1, and non-prime status.
 */
Factor::Factor()
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
void Factor::printpp()
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


/**
 * @brief Decomposes a number into a base and exponent if it is a perfect power.
 *
 * Determines the smallest base `k` and largest exponent `m` such that `N = k^m`.
 * That is, for all other bases (a > 1) and exponents (b > 1), such that `N = a^b`, it holds `m >= b`.
 *
 * @param N The number to decompose.
 * @return std::pair<mpz_class, unsigned int> A pair containing the base `k` and exponent `m`.
 * @throw std::invalid_argument if `N` is not a perfect power.
 */
std::pair<mpz_class, unsigned int> get_smallest_base_biggest_exponent_for_perfect_power(const mpz_class& N) {
    // Ensure N is greater than 1 (no meaningful decomposition for N <= 1)
    if (N <= 1) {
        throw std::invalid_argument("N must be greater than 1.");
    }

    // If N = k^m (for natural numbers k,m >1), than m <= ceil(log2(N))
    mpz_class k;
    const unsigned int upper_bound_exponent = ceil(log2(N.get_d()));
    for (unsigned int i = upper_bound_exponent; i >= 2; --i) {  // test potential exponents, starting with upper bound
        // mpz_root returns non-zero if the computation was exact, i.e., if N is k to the i-th power.
        if (mpz_root(k.get_mpz_t(), N.get_mpz_t(), i) > 0) {
            // Verify the result by computing k^i and comparing with N
            mpz_class check;
            mpz_pow_ui(check.get_mpz_t(), k.get_mpz_t(), i);
            if (check == N) {
                return {k, i}; // Return the pair (k, m)
            }
        }
    }
    throw std::invalid_argument("N: " + N.get_str() + " is not a perfect power.");
}


UnknownNumberModeException::UnknownNumberModeException(const std::string& message)
        : std::invalid_argument(message) {}


/**
 * @brief Generates numbers for factorization based on user input and mode.
 *
 * Converts a user-provided input (string and mode) into a number to be factored.
 * Function implementations are in `test_numbers.h` and `test_numbers.cpp`.
 *
 * @param number_as_string The number provided by the user.
 * @param number_mode Mode to interpret the input number.
 * @return mpz_class The generated number.
 */
mpz_class generate_number_from_user_input(const std::string& number_as_string, const NumberMode number_mode) {
    switch (number_mode) {
        case DirectNumber: return mpz_class(number_as_string);
        case FermatNumber: return Fermat(std::stoi(number_as_string));
        case CunninghamNumber: return Cunningham(std::stoi(number_as_string));
        case TestNumber: return TestNumberB(std::stoi(number_as_string));
        case RSANumber: return TestRSANumber(std::stoi(number_as_string));
        default: throw UnknownNumberModeException("Unknown mode " + std::to_string(number_mode) + " provided.");
    }
}
