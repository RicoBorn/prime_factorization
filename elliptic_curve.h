#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include <iostream>
#include <string>
#include <gmpxx.h>
#include <list>
#include <stdexcept>

#include "helper.h"


extern const int DEFAULT_NUM_CURVES;  // Number of elliptic curves to try before increasing bounds
extern const int DEFAULT_B;  // stage-1 bound: if set to 0, it will be calculated at run time (as a function of C)
extern const int DEFAULT_C;  // stage-2 bound: if set to 0, it will be calculated at run time (as a function of N)
extern const bool DEFAULT_IN_PARALLEL;  // Whether to run lenstra algorithm in multiple threads
extern const int DEFAULT_NUM_THREADS;  // Number of threads to use when running lenstra algorithm in parallel


/**
 * @class SingularEllipticCurveException
 * @brief Custom exception thrown when trying to initialize a singular elliptic curve.
 */
class SingularEllipticCurveException : public std::invalid_argument {
public:
    mpz_class a;     ///< a \in Z/nZ (coefficient in the short Weierstrass-Equation).
    mpz_class b;     ///< b \in Z/nZ (coefficient in the short Weierstrass-Equation).
    mpz_class n;     ///< Modulus n.
    mpz_class gcd;   ///< GCD of ((4 * a^3 + 27 * b^2) mod n) and n. To be a singular curve, this must not be 1

    SingularEllipticCurveException(const mpz_class& a, const mpz_class& b, const mpz_class& n, const mpz_class& gcd);
    ~SingularEllipticCurveException() _NOEXCEPT override = default;

private:
    static std::string create_message(const mpz_class& a, const mpz_class& b, const mpz_class& n, const mpz_class& gcd);
};


/**
 * @class NonInvertibleElementException
 * @brief Exception thrown when an element is not invertible modulo n in Z/nZ.
 */
class NonInvertibleElementException : public std::invalid_argument {
public:
    mpz_class n;        ///< Modulus n.
    mpz_class element;  ///< Element for which inverse was attempted.
    mpz_class gcd;      ///< GCD of the element and n.

    NonInvertibleElementException(const mpz_class& n, const mpz_class& element, const mpz_class& gcd);
    ~NonInvertibleElementException() _NOEXCEPT override = default;

private:
    static std::string create_message(const mpz_class& n, const mpz_class& element, const mpz_class& gcd);
};


/**
 * @class PointNotOnCurveException
 * @brief Custom exception thrown when a point is not on the elliptic curve.
 */
class PointNotOnCurveException : public std::invalid_argument {
public:
    explicit PointNotOnCurveException(const std::string& message);
};


/**
 * @class DiscriminantMultipleOfNException
 * @brief Custom exception thrown when discriminant of elliptic curve is multiple of N.
 */
class DiscriminantMultipleOfNException : public std::invalid_argument {
public:
    explicit DiscriminantMultipleOfNException(const std::string& message);
};


/**
 * @class UnsuccessfulLenstraAlgorithmException
 * @brief Custom exception thrown when reaching end of lenstra algorithm, that is, it terminated unsuccessfully.
 */
class UnsuccessfulLenstraAlgorithmException : public std::runtime_error {
public:
    explicit UnsuccessfulLenstraAlgorithmException(const std::string& message);
    explicit UnsuccessfulLenstraAlgorithmException(const char * message);
};


/**
* @brief Computes the stage-1 bound (i.e., prime base bound) `B` from the stage-2 bound (i.e., the smallest prime bound `C`).
 */
mpz_class calculate_base_prime_bound_from_smallest_prime_bound(const mpz_class& C);


/**
* @brief Calculates stage-2 bound (bound for smallest prime) from N.
 */
mpz_class get_smallest_prime_bound(mpz_class& N, const bool is_rsa_number);


/**
 * @brief Computes the exponent `e` such that e = floor(log_p(C + 2*sqrt(C) + 1)).
 */
mpz_class get_exp_for_prime(const mpz_class& C, const mpz_class& p);


/**
 * @brief Validates the input N for the Lenstra elliptic curve factorization algorithm.
 */
void check_input_for_lenstra(const mpz_class& N);


/**
 * @class ECPoint
 * @brief Represents a point on an elliptic curve over ring Z/nZ.
 *
 * This class models a point on an elliptic curve, characterized by coordinates (x, y)
 * It supports two types of points:
 * - Regular points with specific (x, y) coordinates.
 * - The "point at infinity," represented as a special case often used in elliptic curve operations.
 */
class ECPoint {
public:
    mpz_class x, y;
    bool is_infinity;

    ECPoint(const mpz_class& x, const mpz_class& y);
    ECPoint();
    [[nodiscard]] std::string to_string() const;
};


/**
 * @class EllipticCurve
 * @brief Represents an elliptic curve over the ring Z/nZ.
 *
 * This class models an elliptic curve defined over the finite ring Z/nZ, specified by the
 * Weierstrass equation: y^2 = x^3 + ax + b (mod n). The curve is initialized with parameters
 * (a, b, n), and checks are performed to ensure that it is non-singular.
 *
 * Non-singularity is enforced by requiring that the discriminant, (4 * a^3 + 27 * b^2) mod n,
 * is invertible in Z/nZ. If this condition fails, an exception is thrown.
 */
class EllipticCurve {
private:
    mpz_class a, b, n;

public:
    EllipticCurve(const mpz_class& a, const mpz_class& b, const mpz_class& n);

    [[nodiscard]] const mpz_class& get_a() const;
    [[nodiscard]] const mpz_class& get_b() const;
    [[nodiscard]] const mpz_class& get_modulus() const;
    [[nodiscard]] std::string to_string() const;
    [[nodiscard]] bool point_is_on_curve(const ECPoint& point) const;
    [[nodiscard]] bool points_are_equal(const ECPoint& P, const ECPoint& Q) const;
    [[nodiscard]] ECPoint double_point(const ECPoint& P) const;
    [[nodiscard]] ECPoint add_unequal_points(const ECPoint& P, const ECPoint& Q) const;
    [[nodiscard]] ECPoint add_points(const ECPoint& P, const ECPoint& Q) const;
    [[nodiscard]] ECPoint scalar_multiplication(const mpz_class& k, const ECPoint& P) const;
};


/**
 * @brief Runs the Lenstra elliptic curve factorization algorithm in a separate thread.
 */
void run_lenstra_in_thread(const mpz_class& N, const mpz_class& B, const mpz_class& C,  bool& found_divisor, std::mutex& task_mutex, mpz_class& divisor);


/**
 * @brief Executes Lenstra's elliptic curve factorization algorithm multiple times in parallel.
 */
mpz_class run_lenstra_algorithm_multiple_times_in_parallel(const mpz_class& N, const mpz_class& B, const mpz_class& C, const int& m_times, const int& n_threads);


/**
 * @brief Executes Lenstra's elliptic curve factorization algorithm on the given number.
 */
mpz_class run_lenstra_algorithm(const mpz_class& N, const mpz_class& B, const mpz_class& C);


/**
 * @brief Executes the Lenstra elliptic curve factorization algorithm multiple times.
 */
mpz_class run_lenstra_algorithm_multiple_times(const mpz_class& N, const mpz_class& B, const mpz_class& C, const int& m_times);


/**
 * @brief Handles the factorization of a perfect power by decomposing it into its base and exponent.
 */
void handle_perfect_power(mpz_class& N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors, const bool& in_parallel, const int& n_threads);


/**
 * @brief Factorizes a number using elliptic curve factorization.
 */
void factorize_with_elliptic_curves(mpz_class& N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors, const bool& in_parallel, const int& n_threads);


#endif //ELLIPTIC_CURVE_H
