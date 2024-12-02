#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <stdexcept>
#include <gmpxx.h>


/**
 * @brief Helper function to check if a string represents a positive integer (including 0).
 */
bool is_positive_integer(const std::string& s);


/**
 * @brief Computes the representation of d as an element in the ring Z/nZ.
 */
mpz_class mod_ring(const mpz_class& d, const mpz_class& n);


/**
 * @brief Calculates the number of decimal digits in a given integer.
 */
size_t get_number_of_decimal_digits(const mpz_class& N);


/**
 * @brief Divides out the maximum power of a number P from T. The latter will be updated by the function.
 */
unsigned int divide_out_maximal_power(mpz_class& T, const mpz_class& P);


/**
 * @class Factor
 * @brief Represents a (prime) factor with an associated exponent.
 */
class Factor {
public:
    unsigned int exponent;  ///< Exponent to represent powers in factorization (e.g., 2^6)
    mpz_class factor; ///< Factor, ideally a prime number
    int is_prime; ///< Indicates if 'factor' is prime (2 = definitely prime, 1 = probably prime, 0 = not prime)

    /**
     * @brief Default constructor initializes the factor with exponent 0, factor 1, and non-prime status.
     */
    Factor();

    /**
     * @brief Prints the factor in a specific format based on its properties.
     */
    void printpp();
};


/**
 * @brief Decomposes a number into a base and exponent if it is a perfect power.
 */
std::pair<mpz_class, unsigned int> get_smallest_base_biggest_exponent_for_perfect_power(const mpz_class& N);


/**
 * @class UnknownNumberModeException
 * @brief Custom exception thrown when an invalid number mode is provided.
 */
class UnknownNumberModeException : public std::invalid_argument {
public:
    explicit UnknownNumberModeException(const std::string& message);
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
 */
mpz_class generate_number_from_user_input(const std::string& number_as_string, const NumberMode number_mode);

#endif
