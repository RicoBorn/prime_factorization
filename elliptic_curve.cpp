#include <iostream>
#include <string>
#include <gmpxx.h>
#include <list>
#include <stdexcept>
#include <random>
#include <cmath>
#include <thread>
#include <vector>
#include <mutex>

#include "elliptic_curve.h"
#include "helper.h"


const int DEFAULT_NUM_CURVES = 80;  // Number of elliptic curves to try before increasing bounds
const int DEFAULT_B = 0;  // stage-1 bound: if set to 0, it will be calculated at run time (as a function of C)
const int DEFAULT_C = 0;  // stage-2 bound: if set to 0, it will be calculated at run time (as a function of N)
const bool DEFAULT_IN_PARALLEL = false;  // Whether to run lenstra algorithm in multiple threads
const int DEFAULT_NUM_THREADS = 10;  // Number of threads to use when running lenstra algorithm in parallel


SingularEllipticCurveException::SingularEllipticCurveException(const mpz_class& a, const mpz_class& b, const mpz_class& n, const mpz_class& gcd)
        : std::invalid_argument(create_message(a, b, n, gcd)), a(a), b(b), n(n), gcd(gcd) {
}

std::string SingularEllipticCurveException::create_message(const mpz_class& a, const mpz_class& b, const mpz_class& n, const mpz_class& gcd) {
        std::string message(
            "Singular curve with inputs a: " + a.get_str()
            + ", b: " + b.get_str()
            + ", n: " + n.get_str()
            + "; (4 * a^3 + 27 * b^2) mod n is not invertible in Z/nZ, with GCD: " + gcd.get_str());
        return message;
}


NonInvertibleElementException::NonInvertibleElementException(const mpz_class& n, const mpz_class& element, const mpz_class& gcd)
        : std::invalid_argument(create_message(n, element, gcd)), n(n), element(element), gcd(gcd) {
}

std::string NonInvertibleElementException::create_message(const mpz_class& n, const mpz_class& element, const mpz_class& gcd) {
        std::string message(
            "Element " + element.get_str()
            + " is not invertible modulo " + n.get_str()
            + ". GCD with modulus: " + gcd.get_str());
        return message;
}


PointNotOnCurveException::PointNotOnCurveException(const std::string& message)
        : std::invalid_argument(message) {
}


DiscriminantMultipleOfNException::DiscriminantMultipleOfNException(const std::string& message)
        : std::invalid_argument(message) {
}


UnsuccessfulLenstraAlgorithmException::UnsuccessfulLenstraAlgorithmException(const std::string& message)
        : std::runtime_error(message) {
}

UnsuccessfulLenstraAlgorithmException::UnsuccessfulLenstraAlgorithmException(const char * message)
        : std::runtime_error(message) {
}


/**
 * @brief Computes the stage-1 bound (i.e., prime base bound) `B` from the stage-2 bound (i.e., the smallest prime bound `C`).
 *
 * Uses the formula `B = exp(sqrt(log(C) * log(log(C)) / 2))` to compute an approximate prime bound.
 *
 * @param C The smallest prime bound parameter.
 * @return mpz_class The calculated prime bound as an integer.
 */
mpz_class calculate_base_prime_bound_from_smallest_prime_bound(const mpz_class& C) {
    if (C <= 3) {  // for C < e, ln(ln(C)) not defined
        return C;  // simply return C for those edge cases
    }

    // Convert C to double: there is no log and exp function in GMP library (we might lose precision at this point)
    const double C_double = C.get_d();
    // Calculate B as double
    const double B_double = exp(sqrt((log(C_double)*log(log(C_double)))/2));

    mpz_class B(B_double + 1);  // (B_double + 1) because constructor of mpz_class is flooring

    return B;
}


/**
 * @brief Calculates stage-2 bound (bound for smallest prime) from N.
 *
 * This function calculates the stage-2 bound as follows:
 *   - If: N is known to be RSA-Number there is not much value of setting the bound for smallest prime to much less than sqrt(N)
 *     we calculate bound = floor(sqrt(N))/100 (as a heuristic)
 *   - If: N is a small number (i.e., up to 20 digits), we simply calculate bound = floor(sqrt(N))
 *   - Else: we follow: https://www.maplesoft.com/applications/Preview.aspx?id=3528 and calculate
 *     bound = max(floor(log2(N)), known_viable_bound), with known_viable_bound=100000000000 (is based on empirical values)
 *
 * @param N The (composite) number.
 * @param is_rsa_number Boolean flag whether number is known to be an RSA-Number.
 * @return mpz_class The stage-2 bound.
 */
mpz_class get_smallest_prime_bound(mpz_class& N, const bool is_rsa_number) {
    mpz_class sqrt_N;

    if (is_rsa_number) {
        mpz_sqrt(sqrt_N.get_mpz_t(), N.get_mpz_t());
        return sqrt_N / 100;
    }

    if (get_number_of_decimal_digits(N) <= 20) {
        mpz_sqrt(sqrt_N.get_mpz_t(), N.get_mpz_t());
        return sqrt_N;

    }

    const mpz_class viable_bound("100000000000");  // empirical value; known to result in reasonable run times
    const size_t msb = mpz_sizeinbase(N.get_mpz_t(), 2);
    mpz_class bound = msb - 1;
    if (viable_bound >= bound) {
        return viable_bound;
    }

    return bound;
}


/**
 * @brief Computes the exponent `e` such that `e = floor(log_p(C + 2*sqrt(C) + 1))`.
 *
 * This function is part of the Lenstra algorithm. The term `C + 2*sqrt(C) + 1` stems from the Hasse Bound (1933).
 * Since the GMP library lacks native logarithmic and exponential functions, precision may be lost when converting to and from double.
 *
 * @param C The input integer for which the computation is performed.
 * @param p The prime number to compute the exponent for.
 * @return mpz_class The exponent `e` as an arbitrary-precision integer.
 */
mpz_class get_exp_for_prime(const mpz_class& C, const mpz_class& p) {
    // Convert C to double: there is no log and exp function in GMP library (we might lose precision at this point)
    const double C_double = C.get_d();
    // step 1: calculate C + 2*sqrt(C) + 1
    const double term = C_double + 2 * sqrt(C_double) + 1;
    // calculate log with basis p (using logarithmic laws)
    const double e_double = log(term) / log(p.get_d());
    mpz_class e(e_double);

    return e;
}


/**
 * @brief Validates the input `N` for the Lenstra elliptic curve factorization algorithm.
 *
 * Ensures that `N` is not divisible by 2 and 3, and that `N` is not a perfect power.
 *
 * @param N The number to be validated.
 * @throw std::invalid_argument if the number is divisible by 2, 3, or is a perfect power.
 */
void check_input_for_lenstra(const mpz_class& N) {
    if (N % 2 == 0) {
        throw std::invalid_argument("N: " + N.get_str() + "is divisible by 2.");
    }
    if (N % 3 == 0) {
        throw std::invalid_argument("N: " + N.get_str() + "is divisible by 3.");
    }
    if (mpz_perfect_power_p(N.get_mpz_t()) != 0) {
        throw std::invalid_argument("N: " + N.get_str() + "is a perfect power.");
    }

}


/**
 * @class ECPoint
 * @brief Represents a point on an elliptic curve over ring Z/nZ.
 *
 * This class models a point on an elliptic curve, characterized by coordinates (x, y)
 * It supports two types of points:
 * - Regular points with specific (x, y) coordinates.
 * - The "point at infinity," represented as a special case often used in elliptic curve operations.
 */
ECPoint::ECPoint(const mpz_class& x, const mpz_class& y) : x(x), y(y), is_infinity(false) {}
ECPoint::ECPoint() : is_infinity(true) {}
std::string ECPoint::to_string() const {
        if (is_infinity) return "At Infinity";
        return "(" + x.get_str() + ", " + y.get_str() + ")";
}


/**
 * @brief Constructs a new elliptic curve with specified coefficients and modulus.
 * @param a The coefficient 'a' in the elliptic curve equation (i.e., the short Weierstrass-Equation).
 * @param b The coefficient 'b' in the elliptic curve equation.
 * @param n The modulus defining the ring Z/nZ over which the curve is defined.
 * @throws SingularEllipticCurveException if the curve is singular (discriminant is not invertible).
 *
 * This constructor initializes an elliptic curve using the given parameters, and verifies
 * that the discriminant (4 * a^3 + 27 * b^2) mod n is invertible in Z/nZ.
 */
EllipticCurve::EllipticCurve(const mpz_class& a, const mpz_class& b, const mpz_class& n): a(a), b(b), n(n) {
    // Calculate discriminant with modular reduction for intermediate results (to prevent growth of values)
    mpz_class discriminant = mod_ring(((4 * a * a * a) % n) + ((27 * b * b) % n), n);  // (4 * a^3 + 27 * b^2) mod n

    // Check if discriminant is invertible in Z/nZ
    mpz_class gcd_result;
    mpz_gcd(gcd_result.get_mpz_t(), discriminant.get_mpz_t(), n.get_mpz_t());

    if (gcd_result != 1) {
        throw SingularEllipticCurveException(a, b, n, gcd_result);
    }
}

const mpz_class& EllipticCurve::get_a() const { return a; }
const mpz_class& EllipticCurve::get_b() const { return b; }
const mpz_class& EllipticCurve::get_modulus() const { return n; }

std::string EllipticCurve::to_string() const {
    return "E(" + a.get_str() + "," + b.get_str() + ",Z/" + n.get_str() + "Z)";
}

/**
 * @brief Checks if a given point is on the elliptic curve.
 * @param point The ECPoint to be checked.
 * @return True if the point lies on the elliptic curve, false otherwise.
 *
 * This function verifies that the given point satisfies the curve's equation:
 * y^2 ≡ x^3 + ax + b (mod n). For the point at infinity, it returns true.
 */
bool EllipticCurve::point_is_on_curve(const ECPoint& point) const {
    if (point.is_infinity) {
        // The point at infinity is always considered on the curve.
        return true;
    }

    // Calculate (y^2 - (x^3 + ax + b)) mod n
    mpz_class equation = (point.y * point.y) % n;  // y^3 mod n
    equation -= (point.x * point.x * point.x) % n; // - x^3 mod n
    equation = mod_ring(equation - ((a * point.x) % n) - b, n);   // (y^2 - (x^3 + ax + b)) mod n

    // Check if the equation is zero modulo n
    return equation == 0;
}

/**
 * @brief Checks if two points are equal modulo n.
 * @param P The first ECPoint to compare.
 * @param Q The second ECPoint to compare.
 * @return True if the points are equal modulo n, false otherwise.
 *
 * This method considers two points equal if:
 * - Both are the "point at infinity."
 * - Both have the same x and y coordinates modulo n.
 */
bool EllipticCurve::points_are_equal(const ECPoint& P, const ECPoint& Q) const {
    // Check if both are the point at infinity
    if (P.is_infinity && Q.is_infinity) {
        return true;
    }

    // If one is at infinity and the other is not, they are not equal
    if (P.is_infinity || Q.is_infinity) {
        return false;
    }

    // Compare x and y coordinates modulo n
    return (mod_ring(P.x, n) == mod_ring(Q.x,n)) && (mod_ring(P.y, n) == mod_ring(Q.y, n));
}

/**
 * @brief Performs the point doubling operation on the elliptic curve.
 * @param P The point to double.
 * @return The resulting point after doubling.
 * @throws PointNotOnCurveException if P is not on the curve.
 * @throws NonInvertibleElementException if the slope denominator is not invertible.
 */
ECPoint EllipticCurve::double_point(const ECPoint& P) const {
    if (!point_is_on_curve(P)) {
        throw PointNotOnCurveException("Point " + P.to_string() + " is not on curve " + (this)->to_string());
    }

    if (P.is_infinity || P.y % n == 0) {
        return ECPoint(); // Result is the point at infinity.
    }

    // Calculate the slope (lambda)
    const mpz_class numerator = mod_ring(3 * P.x * P.x + a, n);
    mpz_class denominator = mod_ring(2 * P.y, n);  // will be inverted in next step

    if (mpz_invert(denominator.get_mpz_t(), denominator.get_mpz_t(), n.get_mpz_t()) == 0) {
        mpz_class gcd;
        mpz_gcd(gcd.get_mpz_t(), denominator.get_mpz_t(), n.get_mpz_t());
        throw NonInvertibleElementException(n, denominator, gcd);
    }

    const mpz_class lambda = mod_ring(numerator * denominator, n);

    // Calculate resulting coordinates
    const mpz_class x3 = mod_ring(lambda * lambda - 2 * P.x, n);
    const mpz_class y3 = mod_ring(lambda * (P.x - x3) - P.y, n);

    return ECPoint(x3, y3);
}

/**
 * @brief Adds two unequal points on the elliptic curve.
 * @param P The first point.
 * @param Q The second point.
 * @return The resulting point after addition.
 * @throws PointNotOnCurveException if either P or Q is not on the curve.
 * @throws std::invalid_argument if P and Q are equal.
 * @throws NonInvertibleElementException if the slope denominator is not invertible.
 */
ECPoint EllipticCurve::add_unequal_points(const ECPoint& P, const ECPoint& Q) const {
    if (!point_is_on_curve(P)) {
        throw PointNotOnCurveException("Point " + P.to_string() + " is not on curve " + (this)->to_string());
    }
    if (!point_is_on_curve(Q)) {
        throw PointNotOnCurveException("Point " + Q.to_string() + " is not on curve " + (this)->to_string());
    }
    if (points_are_equal(P, Q)) {
        throw std::invalid_argument("Points " + P.to_string() + " and " + Q.to_string() + " are equal. Use point doubling instead.");
    }

    if (P.is_infinity) return Q;
    if (Q.is_infinity) return P;
    // Check if Q = -P
    const ECPoint inverse_P(P.x, -1*P.y);
    if (points_are_equal(inverse_P, Q)) return ECPoint(); // Result is the point at infinity.

    // Calculate the slope (lambda)
    const mpz_class numerator = mod_ring(Q.y - P.y, n);
    mpz_class denominator = mod_ring(Q.x - P.x, n);

    if (mpz_invert(denominator.get_mpz_t(), denominator.get_mpz_t(), n.get_mpz_t()) == 0) {
        mpz_class gcd;
        mpz_gcd(gcd.get_mpz_t(), denominator.get_mpz_t(), n.get_mpz_t());
        throw NonInvertibleElementException(n, denominator, gcd);
    }

    const mpz_class lambda = mod_ring(numerator * denominator, n);

    // Calculate resulting coordinates
    const mpz_class x3 = mod_ring(lambda * lambda - P.x - Q.x, n);
    const mpz_class y3 = mod_ring(lambda * (P.x - x3) - P.y, n);

    return ECPoint(x3, y3);
}

/**
 * @brief Adds two points on the elliptic curve after verifying they lie on the curve.
 * @param P First point on the elliptic curve.
 * @param Q Second point on the elliptic curve.
 * @return The sum of P and Q as a new ECPoint.
 * @throws PointNotOnCurveException if either P or Q is not on the curve.
 * @throws NonInvertibleElementException if the denominator in slope calculation is not invertible.
 */
ECPoint EllipticCurve::add_points(const ECPoint& P, const ECPoint& Q) const {
    if (points_are_equal(P, Q)) {
        return double_point(P);
    }

    return add_unequal_points(P, Q);
}

/**
 * @brief Performs scalar multiplication of a point on the elliptic curve, using the double-and-add algorithm.
 * @param k The scalar by which to multiply the point.
 * @param P The point to be multiplied.
 * @return The resulting point after scalar multiplication.
 * @throws PointNotOnCurveException if the point P is not on the curve.
 */
ECPoint EllipticCurve::scalar_multiplication(const mpz_class& k, const ECPoint& P) const {
    if (!point_is_on_curve(P)) {
        throw PointNotOnCurveException("Point " + P.to_string() + " is not on curve " + (this)->to_string());
    }

    ECPoint result; // Initialize result as the point at infinity (neutral element for addition)
    ECPoint current = P;
    mpz_class scalar = k;

    while (scalar > 0) {
        // Check the least significant bit
        if (scalar % 2 == 1) {
            if (result.is_infinity) {  // here, no computation needed
                result = current;
            } else {
                result = add_points(result, current);
            }
        }

        // Double the point for the next bit
        current = double_point(current);

        // Right shift the scalar (equivalent to dividing by 2 and flooring the result to next integer)
        scalar /= 2;
    }

    return result;
}


/**
 * @brief Executes Lenstra's elliptic curve factorization algorithm on the given number.
 *
 * This function tries to factorize `N` using an elliptic curve over the field `Z/NZ`.
 * The algorithm iterates through primes up to `B` and applies scalar multiplication
 * on the elliptic curve.
 *
 * @param N The number to be factorized.
 * @param B The bound for the primes used in scalar multiplication.
 * @param C The bound for the smallest prime.
 * @return mpz_class A non-trivial divisor of `N` if found.
 * @throw UnsuccessfulLenstraAlgorithmException if no factor is found after all primes.
 * @throw DiscriminantMultipleOfNException if the curve discriminant is a multiple of `N`.
 * @throw NonInvertibleElementException if a non-invertible element over `Z/NZ` is encountered. This is the `success` case.
 */
mpz_class run_lenstra_algorithm(const mpz_class& N, const mpz_class& B, const mpz_class& C)
{
    check_input_for_lenstra(N);  // check if N satisfies prerequisites for the algorithm

    // Random values for a,u,v
    mpz_class a, u, v;
    gmp_randstate_t state;
    gmp_randinit_default(state); // Initialize the random state
    gmp_randseed_ui(state, std::random_device{}()); // Seed with a random value
    mpz_urandomm(a.get_mpz_t(), state, N.get_mpz_t());
    mpz_urandomm(u.get_mpz_t(), state, N.get_mpz_t());
    mpz_urandomm(v.get_mpz_t(), state, N.get_mpz_t());

    // Calculate b in Z/NZ, such that P(u,v) is on E(a,b,Z/NZ)
    const mpz_class b = mod_ring(((v * v) % N) - ((u * u * u) % N) - ((a * u) % N), N); // b = (v^2 - u^3 - au) mod N

    try {
        const EllipticCurve E(a,b,N);  // Try to initialize curve, may result in SingularEllipticCurveException
        const ECPoint P(u,v);  // Elliptic point initialization
        ECPoint Q(P.x, P.y);  // Q = P

        // Iterate through primes up to B
        mpz_class prime("2"); // Start with the first prime
        while (prime <= B) {
            mpz_class e = get_exp_for_prime(C,prime);
            // For-loop from 0 up to e-1
            for (mpz_class j = 0; j < e; ++j) {
                Q = E.scalar_multiplication(prime,Q);
            }
            // Get the next prime
            mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
        }
    } catch (const SingularEllipticCurveException& exc) {
        if (exc.gcd > 1 && exc.gcd < N) {  // we found non-trivial divisor of N
            return exc.gcd;
        }
        // Now, we know exc.gcd == N
        throw DiscriminantMultipleOfNException("Discriminant is multiple of N: " + N.get_str());
    } catch (const NonInvertibleElementException& exc) {  // we found an element not invertible over Z/nZ
        if (exc.gcd > 1 && exc.gcd < N) {  // we found non-trivial divisor of N
            return exc.gcd;
        }
        throw UnsuccessfulLenstraAlgorithmException("Unlucky event GCD(element, N) = N ");
    }
    // No factor found, algorithm unsuccessful
    throw UnsuccessfulLenstraAlgorithmException("Reached end of lenstra algorithm.");
}


/**
 * @brief Runs the Lenstra elliptic curve factorization algorithm in a separate thread.
 *
 * This function is designed to execute the `run_lenstra_algorithm` in a thread. It attempts
 * to find a non-trivial divisor of the input number `N` using elliptic curve factorization.
 * The function acquires a lock to ensure that only one thread can update the shared variables
 * (`found_divisor` and `divisor`) at a time, preventing race conditions in a multithreaded
 * environment. If a divisor is found, the thread updates the shared `divisor` and sets the
 * `found_divisor` flag to `true` to signal that a factor has been discovered.
 *
 * The function catches any exceptions thrown by the `run_lenstra_algorithm` function and ignores
 * them if no divisor is found, allowing other threads to continue processing.
 *
 * @param N The number to be factorized. This is passed to the Lenstra algorithm.
 * @param B The bound for the primes used in scalar multiplication for elliptic curve factorization.
 * @param C The bound for the smallest prime, used in the calculation of the number of scalar
 *          multiplications per prime.
 * @param found_divisor A Boolean variable that signals if a divisor has been found.
 * @param task_mutex A mutex used for synchronizing access to shared variables.
 * @param divisor A `mpz_class` where the divisor will be stored if found.
 */
void run_lenstra_in_thread(const mpz_class& N, const mpz_class& B, const mpz_class& C,  bool& found_divisor, std::mutex& task_mutex, mpz_class& divisor) {
    try {
        const mpz_class tmp(run_lenstra_algorithm(N, B, C));
        std::lock_guard<std::mutex> lock(task_mutex);  // Lock for other threads
        if (!found_divisor) {  // Check if a divisor has already been found
            divisor = tmp;
            found_divisor = true;
        }
    } catch (const DiscriminantMultipleOfNException&) {
    } catch (const UnsuccessfulLenstraAlgorithmException&) {
    }
}


/**
 * @brief Executes Lenstra's elliptic curve factorization algorithm multiple times in parallel.
 *
 * This function attempts to factorize the number `N` using the Lenstra elliptic curve
 * factorization algorithm, running it up to `m_times` times in parallel across multiple threads.
 * The algorithm is executed concurrently with a maximum of `n_threads` threads running at a time.
 * If a non-trivial divisor of `N` is found by any thread, the result is returned immediately.
 * If all attempts fail, an exception is thrown after all attempts have been made.
 *
 * @param N The number to be factorized. This is passed to the Lenstra algorithm.
 * @param B The bound for the primes used in scalar multiplication for elliptic curve factorization.
 * @param C The bound for the smallest prime, used in the calculation of the number of scalar
 *          multiplications per prime.
 * @param m_times The number of times the Lenstra algorithm should be executed in parallel.
 * @param n_threads The maximum number of threads to run concurrently.
 * @return mpz_class A non-trivial divisor of `N` if found by any of the parallel attempts.
 * @throw UnsuccessfulLenstraAlgorithmException if no divisor is found after `m_times` parallel attempts.
 */
mpz_class run_lenstra_algorithm_multiple_times_in_parallel(const mpz_class& N, const mpz_class& B, const mpz_class& C, const int& m_times, const int& n_threads) {
    // Synchronization variables
    bool found_divisor = false;  // Indicates whether a divisor has been found by one of the threads
    mpz_class divisor("1");
    std::mutex task_mutex;

    int attempts = 0;

    // Repeated execution of threads in blocks
    while (attempts < m_times && !found_divisor) {
        std::vector<std::thread> threads;

        // Start up to n_threads or the remaining tasks
        for (int j = 0; j < n_threads && attempts < m_times; ++j) {
            threads.emplace_back(run_lenstra_in_thread, N, B, C, std::ref(found_divisor), std::ref(task_mutex), std::ref(divisor));
            ++attempts;  // Increment the total number of attempts
        }

        // Wait for all threads to complete
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }
    }

    // Check if valid divisor was found
    if (found_divisor && divisor > 1) {
        return divisor;
    }

    throw UnsuccessfulLenstraAlgorithmException("All " + std::to_string(m_times) + " attempts of Lenstra algorithm failed.");
}


/**
 * @brief Executes the Lenstra elliptic curve factorization algorithm multiple times.
 *
 * Repeats the Lenstra algorithm up to `m_times` in an attempt to find a factor.
 * If all attempts fail, an exception is thrown.
 *
 * @param N The number to be factorized.
 * @param B The bound for primes used in the algorithm.
 * @param C The bound for the smallest prime.
 * @param m_times Number of attempts to perform.
 * @return mpz_class A non-trivial divisor of `N` if found.
 * @throw UnsuccessfulLenstraAlgorithmException if all attempts fail.
 */
mpz_class run_lenstra_algorithm_multiple_times(const mpz_class& N, const mpz_class& B, const mpz_class& C, const int& m_times)
{
    for (int i = 0; i < m_times; ++i) {
        try {
            return run_lenstra_algorithm(N, B, C);
        } catch (const DiscriminantMultipleOfNException&) {
        } catch (const UnsuccessfulLenstraAlgorithmException&) {
        }
    }

    // If all attempts failed, throw exception
    throw UnsuccessfulLenstraAlgorithmException("All " + std::to_string(m_times) + " attempts of Lenstra algorithm failed.");
}


/**
 * @brief Handles the factorization of a perfect power by decomposing it into its base and exponent.
 *
 * This function performs recursive factorization on the base and updates the list of factors.
 *
 * @param N The perfect power to factorize.
 * @param B The prime bound for elliptic curve factorization.
 * @param C Parameter controlling the randomness of curves.
 * @param m_curves Number of elliptic curves to use for factorization.
 * @param factors The list of factors to be updated with results.
 * @param in_parallel Flag whether to run the lenstra algorithm in parallel (multiple threads)
 * @param n_threads The number of threads to use (if the lenstra algorithm run in parallel)
 */
void handle_perfect_power(mpz_class& N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors, const bool& in_parallel, const int& n_threads) {
    // Decompose N into a^b (minimal a, maximal b)
    auto [k, m] = get_smallest_base_biggest_exponent_for_perfect_power(N);

    // Check if k is prime
    const int prime_status = mpz_probab_prime_p(k.get_mpz_t(), 27); // GMP primality test
    if (prime_status > 0) { // k is prime or probably prime; add it to factors (with exponent m) and return
        Factor factor;
        factor.exponent = m;
        factor.factor = k;
        factor.is_prime = prime_status;
        factors.push_back(factor);
        N = mpz_class("1");
        return;
    }

    // k is composite, factorize further; Info: k cannot be a perfect power itself (otherwise m would've not been maximal)
    std::list<Factor> tmp_factors;
    factorize_with_elliptic_curves(k, B, C, m_curves, tmp_factors, in_parallel, n_threads);
    // Now that k is completely factorized, update prime factors (by multiplying all exponents by m)
    for (auto& factor : tmp_factors) {
        factor.exponent *= m;  // Multiply each exponent by m
        factors.push_back(factor);  // Add the updated factor to the main list
    }

}


/**
 * @brief Factorizes a number using elliptic curve factorization.
 *
 * Combines perfect power decomposition, and Lenstra's elliptic curve algorithm
 * to fully factorize `N` into its prime factors.
 *
 * @param N The number to factorize (modified during computation).
 * @param B The prime bound for elliptic curve factorization.
 * @param C Parameter controlling the randomness of curves.
 * @param m_curves Number of elliptic curves to use for factorization.
 * @param factors The list of factors to be updated with results.
 * @param in_parallel Flag whether to run the lenstra algorithm in parallel (multiple threads)
 * @param n_threads The number of threads to use (if the lenstra algorithm run in parallel)
 */
void factorize_with_elliptic_curves(mpz_class& N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors, const bool& in_parallel, const int& n_threads) {
    // First, ensure N is not divisible by 2 or 3 (divide out if necessary)
    Factor factor;
    mpz_class P("2");
    factor.exponent = divide_out_maximal_power(N,P);
    if(factor.exponent > 0)
    {
        factor.factor = P;
        factor.is_prime = 2;
        factors.push_back(factor);

    }
    P++; // P = 3
    factor.exponent = divide_out_maximal_power(N,P);
    if(factor.exponent > 0)
    {
        factor.factor = P;
        factor.is_prime = 2;
        factors.push_back(factor);

    }

    // Base case: stop if N is 1
    if (N == 1)
        return;

    // Check if N is prime and stop
    int prime_status = mpz_probab_prime_p(N.get_mpz_t(), 27); // GMP primality test
    if (prime_status > 0) { // N is prime or probably prime
        factor.exponent = 1;
        factor.factor = N;
        factor.is_prime = prime_status;
        factors.push_back(factor);
        N = mpz_class("1");
        return;
    }

    // Check if N is a perfect power, and handle this case
    if (mpz_perfect_power_p(N.get_mpz_t()) > 0) {
        handle_perfect_power(N, B, C, m_curves, factors, in_parallel, n_threads);
        return; // Stop further factorization since perfect power handling will decompose N
    }

    // Apply Lenstra's algorithm to find a non-trivial divisor
    mpz_class C_tmp = C;
    mpz_class B_tmp = B;
    mpz_class divisor("1");
    while (divisor <= 1) {
        try {
            if (in_parallel) {
                // Run in parallel using multiple threads
                divisor = run_lenstra_algorithm_multiple_times_in_parallel(N, B_tmp, C_tmp, m_curves, n_threads);
            } else {
                // Run sequentially
                divisor = run_lenstra_algorithm_multiple_times(N, B_tmp, C_tmp, m_curves);
            }
        } catch (const UnsuccessfulLenstraAlgorithmException&) {  // Terminated unsuccessfully
            C_tmp = C_tmp * 2; // Increase (i.e, double) C
            B_tmp = calculate_base_prime_bound_from_smallest_prime_bound(C_tmp);  // Update B accordingly
        }
    }

    // Check if the divisor is prime
    prime_status = mpz_probab_prime_p(divisor.get_mpz_t(), 27);
    if (prime_status > 0) { // Divisor is prime or probably prime
        const unsigned int exponent = divide_out_maximal_power(N, divisor);
        // Add the prime factor
        factor.factor = divisor;
        factor.exponent = exponent;
        factor.is_prime = prime_status;
        factors.push_back(factor);
    } else {
        // Divisor is not prime, further factorize it (recursion)
        mpz_class divisor_copy = divisor; // Create a copy to avoid modification (need original divisor below)
        factorize_with_elliptic_curves(divisor_copy, B_tmp, C_tmp, m_curves, factors, in_parallel, n_threads);
        N /= divisor;  // Update N by dividing it by the divisor
    }

    // Factorize the remaining part of N
    factorize_with_elliptic_curves(N, B_tmp, C_tmp, m_curves, factors, in_parallel, n_threads);
}
