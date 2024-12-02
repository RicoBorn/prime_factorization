#include <iostream>
#include <string>
#include <cstdlib>
#include <gmpxx.h>
#include <list>
#include <chrono>
#include <iterator>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <random>
#include <cmath>
#include <utility>

#include "testzahlen.h"

/* compilation: g++ -std=c++11 factorization.cpp testzahlen.cpp -lgmp -lgmpxx */


// Define constants
const std::string TRIAL_DIV_BOUND = "100000000";
const int DEFAULT_NUM_CURVES = 50;  // Number of elliptic curves to try before increasing bounds
const int DEFAULT_B = 0;  // stage-1 bound: if set to 0, it will be calculated at run time (as a function of C)
const int DEFAULT_C = 0;  // stage-2 bound: if set to 0, it will be calculated at run time (as a function of N)
const bool DEFAULT_NO_TRIAL_DIVISION = false;


/**
 * @class SingularEllipticCurveException
 * @brief Custom exception thrown when trying to initialize a singular elliptic curve.
 */
class SingularEllipticCurveException : public std::invalid_argument {
public:
    mpz_class a;     ///< a \in Z/nZ.
    mpz_class b;     ///< b \in Z/nZ.
    mpz_class n;     ///< Modulus n.
    mpz_class gcd;   ///< GCD of ((4 * a^3 + 27 * b^2) mod n) and n. To be a singular curve, this must not be 1

    SingularEllipticCurveException(const mpz_class& a, const mpz_class& b, const mpz_class& n, const mpz_class& gcd)
        : std::invalid_argument(create_message(a, b, n, gcd)), a(a), b(b), n(n), gcd(gcd) {}

    // Explicitly mark the destructor as noexcept to match the base class
    ~SingularEllipticCurveException() _NOEXCEPT override = default;

private:
    // Helper function to format the exception message
    static std::string create_message(const mpz_class& a, const mpz_class& b, const mpz_class& n, const mpz_class& gcd) {
        std::string message(
            "Singular curve with inputs a: " + a.get_str()
            + ", b: " + b.get_str()
            + ", n: " + n.get_str()
            + "; (4 * a^3 + 27 * b^2) mod n is not invertible in Z/nZ, with GCD: " + gcd.get_str());
        return message;
    }
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

    NonInvertibleElementException(const mpz_class& n, const mpz_class& element, const mpz_class& gcd)
        : std::invalid_argument(create_message(n, element, gcd)), n(n), element(element), gcd(gcd) {}

    // Explicitly mark the destructor as noexcept to match the base class
    ~NonInvertibleElementException() _NOEXCEPT override = default;

private:
    // Helper function to format the exception message
    static std::string create_message(const mpz_class& n, const mpz_class& element, const mpz_class& gcd) {
        std::string message(
            "Element " + element.get_str()
            + " is not invertible modulo " + n.get_str()
            + ". GCD with modulus: " + gcd.get_str());
        return message;
    }
};



/**
 * @class PointNotOnCurveException
 * @brief Custom exception thrown when a point is not on the elliptic curve.
 */
class PointNotOnCurveException : public std::invalid_argument {
public:
    explicit PointNotOnCurveException(const std::string& message)
        : std::invalid_argument(message) {}
};

/**
 * @class DiscriminantMultipleOfNException
 * @brief Custom exception thrown when discriminant of elliptic curve is multiple of N.
 */
class DiscriminantMultipleOfNException : public std::invalid_argument {
public:
    explicit DiscriminantMultipleOfNException(const std::string& message)
        : std::invalid_argument(message) {}
};

/**
 * @class UnsuccessfulLenstraAlgorithmException
 * @brief Custom exception thrown when reaching end of lenstra algorithm, that is, it terminated unsuccessfully.
 */
class UnsuccessfulLenstraAlgorithmException : public std::runtime_error {
public:
    explicit UnsuccessfulLenstraAlgorithmException(const std::string& message)
        : std::runtime_error(message) {}

    explicit UnsuccessfulLenstraAlgorithmException(const char * message)
        : std::runtime_error(message) {}
};


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
 * @brief Calculates stage-2 bound (bound for smallest prime) from N.
 *
 * This function calculates the stage-2 bound as follows:
 *   - If: N is known to be RSA-Number there is not value of setting the bound for smallest prime to much less than sqrt(N)
 *     we calculate bound = floor(sqrt(N))/1000
 *   - If: N is a small number (i.e., up to 20 (decimal) digits), we simply calculate bound = floor(sqrt(N))
 *   - Else: it follows: https://www.maplesoft.com/applications/Preview.aspx?id=3528 and calculates bound = floor(log2(N))
 *
 * @param N The (composite) number.
 * @param is_rsa_number Boolean flag whether number is known to be an RSA-Number.
 * @return mpz_class The stage-2 bound.
 */
mpz_class get_smallest_prime_bound(mpz_class& N, const bool is_rsa_number) {
    //ToDo: mit meinem Programm "rumspielen" ich muss sehen, wie lange es dauert, bei verschiedenen C Werten und dass hier mit aufnehmen
    //ToDo: evtl. noch "known lower bound" mit einbinden... (kennen ja durch trial division schon,dass kein kleiner Teiler mehr von N existiert...)

    mpz_class sqrt_N;

    if (is_rsa_number) {
        mpz_sqrt(sqrt_N.get_mpz_t(), N.get_mpz_t());
        return sqrt_N / 100;
    }

    if (get_number_of_decimal_digits(N) < 20) {
        mpz_sqrt(sqrt_N.get_mpz_t(), N.get_mpz_t());
        return sqrt_N;

    }
    const size_t msb = mpz_sizeinbase(N.get_mpz_t(), 2);
    mpz_class bound = msb - 1;

    return bound;
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
class ECPoint {
public:
    mpz_class x, y;
    bool is_infinity;

public:
    ECPoint(const mpz_class& x, const mpz_class& y) : x(x), y(y), is_infinity(false) {}
    ECPoint() : is_infinity(true) {}

    std::string to_string() const {
        if (is_infinity) return "At Infinity";
        return "(" + x.get_str() + ", " + y.get_str() + ")";
    }
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
    EllipticCurve(const mpz_class& a, const mpz_class& b, const mpz_class& n): a(a), b(b), n(n) {
        // Calculate discriminant with modular reduction for intermediate results (to prevent growth of values)
        mpz_class discriminant = mod_ring(((4 * a * a * a) % n) + ((27 * b * b) % n), n);  // (4 * a^3 + 27 * b^2) mod n

        // Check if discriminant is invertible in Z/nZ
        mpz_class gcd_result;
        mpz_gcd(gcd_result.get_mpz_t(), discriminant.get_mpz_t(), n.get_mpz_t());

        if (gcd_result != 1) {
            throw SingularEllipticCurveException(a, b, n, gcd_result);
        }
    }

    const mpz_class& get_a() const { return a; }
    const mpz_class& get_b() const { return b; }
    const mpz_class& get_modulus() const { return n; }

    // to_string method
    std::string to_string() const {
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
    bool point_is_on_curve(const ECPoint& point) const {
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
    bool points_are_equal(const ECPoint& P, const ECPoint& Q) const {
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
    ECPoint double_point(const ECPoint& P) const {
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
    ECPoint add_unequal_points(const ECPoint& P, const ECPoint& Q) const {
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
    ECPoint add_points(const ECPoint& P, const ECPoint& Q) const {
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
    ECPoint scalar_multiplication(const mpz_class& k, const ECPoint& P) const {
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

};


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
                std::cout << "Iteration j = " << j << std::endl;  //ToDo: get out
                Q = E.scalar_multiplication(prime,Q);
            }
            // Dummy statement (replace this with actual logic in the loop)
            std::cout << "Processing prime: " << prime.get_str() << std::endl;  //ToDo: get out
            // Get the next prime
            mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
        }
    } catch (const SingularEllipticCurveException& exc) {
        if (exc.gcd > 1 && exc.gcd < N) {  // we found non-trivial divisor of N; i.e., gdc(((4 * a^3 + 27 * b^2) mod n), n)
            return exc.gcd;
        }
        // Now, we know exc.gcd == N
        throw DiscriminantMultipleOfNException("Discriminant is multiple of N: " + N.get_str());
    } catch (const NonInvertibleElementException& exc) {  // we found an element not invertible over Z/nZ
        return exc.gcd;
    }
    // No factor found, algorithm unsuccessful
    throw UnsuccessfulLenstraAlgorithmException("Reached end of lenstra algorithm.");
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
 * @brief Helper function to check if a string represents a positive integer (including 0).
 *
 * @param s Input string.
 * @return true if the string is a positive integer (including 0); otherwise, false.
 */
bool is_positive_integer(const std::string& s) {
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
 * Converts a user-provided input (string and mode) into a number to be factored.
 * Function implementations are in `testzahlen.h` and `testzahlen.cpp`.
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
        case TestNumber: return TestzahlB(std::stoi(number_as_string));
        case RSANumber: return RSAZahl(std::stoi(number_as_string));
        default: throw UnknownNumberModeException("Unknown mode " + std::to_string(number_mode) + " provided.");
    }
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

// Function prototypes
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
 */
void handle_perfect_power(mpz_class N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors);


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
 */
void factorize_with_elliptic_curves(mpz_class N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors);


void handle_perfect_power(mpz_class N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors) {
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
    factorize_with_elliptic_curves(k, B, C, m_curves, tmp_factors);
    // Now that k is completely factorized, update prime factors (by multiplying all exponents by m)
    for (auto& factor : tmp_factors) {
        factor.exponent *= m;  // Multiply each exponent by m
        factors.push_back(factor);  // Add the updated factor to the main list
    }

};


void factorize_with_elliptic_curves(mpz_class N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors) {
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
        handle_perfect_power(N, B, C, m_curves, factors);
        return; // Stop further factorization since perfect power handling will decompose N
    }

    // Apply Lenstra's algorithm to find a non-trivial divisor
    mpz_class C_tmp = C;
    mpz_class B_tmp = B;
    mpz_class divisor("1");
    while (divisor <= 1) {
        try {
            divisor = run_lenstra_algorithm_multiple_times(N, B_tmp, C_tmp, m_curves);
        } catch (const UnsuccessfulLenstraAlgorithmException&) {  // terminated unsuccessfully
        }
        C_tmp = C_tmp * 2; // Increase (i.e, double) C
        std::cout << "UnsuccessfulLenstraAlgorithmException. Increased C to C_tmp: " << C_tmp << std::endl;  // ToDo: raus
        B_tmp = calculate_base_prime_bound_from_smallest_prime_bound(C_tmp);  // Update B accordingly
        std::cout << "UnsuccessfulLenstraAlgorithmException. Increased B to B_tmp: " << B_tmp << std::endl;  // ToDo: raus

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
        factorize_with_elliptic_curves(divisor_copy, B_tmp, C_tmp, m_curves, factors);  //ToDo: teste mal hier ob mit B und C (statt B_tmp, C_tmp)
        N /= divisor;  // Update N by dividing it by the divisor
    }

    // Factorize the remaining part of N
    factorize_with_elliptic_curves(N, B_tmp, C_tmp, m_curves, factors);
}


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
