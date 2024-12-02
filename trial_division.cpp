#include <iostream>
#include <gmpxx.h>
#include <list>

#include "helper.h"
#include "trial_division.h"


const std::string TRIAL_DIV_BOUND = "100000000";
const bool DEFAULT_NO_TRIAL_DIVISION = false;

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
    const mpz_class C = sqrt(N) + 1;
    if (B > C) B = C;

    std::list<Factor> factors;
    Factor ppt;

    // Check if N is prime (27 tests)
    int a = mpz_probab_prime_p(N.get_mpz_t(), 27);
    if ( a > 0) // N is probably prime
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
    if (B < P) {return factors;}
    ppt.exponent = divide_out_maximal_power(N,P);
    if (ppt.exponent > 0)
    {
        ppt.factor = P;
        ppt.is_prime = 2;
        factors.push_back(ppt);

    }
    P++; // P = 3
    if (B < P) {return factors;}
    ppt.exponent = divide_out_maximal_power(N,P);
    if (ppt.exponent > 0)
    {
        ppt.factor = P;
        ppt.is_prime = 2;
        factors.push_back(ppt);

    }

    P+=3; // P = 6
    while (N > 1 and P <= B )
    {
        ppt.exponent = divide_out_maximal_power(N,P-1);
        if (ppt.exponent > 0)
        {
            ppt.factor = P-1;
            ppt.is_prime = 2;
            factors.push_back(ppt);

        }
        ppt.exponent = divide_out_maximal_power(N,P+1);
        if (ppt.exponent > 0)
        {
            ppt.factor = P+1;
            ppt.is_prime = 2;
            factors.push_back(ppt);
        }
        P += 6;
    }

    if (N > 1)
    {
        a = mpz_probab_prime_p(N.get_mpz_t(), 27);
        if (a > 0) // Remaining number is probably prime
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
