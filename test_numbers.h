#ifndef TEST_NUMBERS_H
#define TEST_NUMBERS_H
#include <gmpxx.h>

/**
 * @brief Helper function for fast exponentiation.
 */
unsigned long int aux_fast_exp(int n, int p);


/**
 * @brief Computes the k-th Fermat number.
 */
mpz_class Fermat(unsigned int k);


/**
 * @brief Generates a Cunningham number with k digits.
 */
mpz_class Cunningham(int k);


/**
 * @brief Generates predefined test numbers B1, ..., B10.
 */
mpz_class TestNumberB(int k);


/**
 * @brief Retrieves the RSA number R_k.
 */
mpz_class TestRSANumber(int k);
#endif
