#ifndef TRIAL_DIVISION_H
#define TRIAL_DIVISION_H

#include <gmpxx.h>
#include <list>
#include "helper.h"

extern const std::string TRIAL_DIV_BOUND;
extern const bool DEFAULT_NO_TRIAL_DIVISION;

/**
 * @brief Factorizes N (and modifies it) using trial division up to a bound B.
 */
std::list<Factor> trial_division_bounded(mpz_class& N, mpz_class B);

#endif //TRIAL_DIVISION_H
