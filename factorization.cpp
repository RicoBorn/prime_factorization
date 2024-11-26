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
#include "factorization.h"

/* compilation: g++ -std=c++11 factorization.cpp testzahlen.cpp -lgmp -lgmpxx */


// Define constants
const std::string TRIAL_DIV_BOUND = "1000000";
const int NUM_ELLIPTIC_CURVES = 50;


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
     * @param a The coefficient 'a' in the elliptic curve equation.
     * @param b The coefficient 'b' in the elliptic curve equation.
     * @param n The modulus defining the ring Z/nZ over which the curve is defined.
     * @throws SingularEllipticCurveException if the curve is singular (discriminant is not invertible).
     *
     * This constructor initializes an elliptic curve using the given parameters, and verifies
     * that the discriminant (4 * a^3 + 27 * b^2) mod n is invertible in Z/nZ. The discriminant
     * is computed with modular reduction at each step to avoid overflow in intermediate values.
     */
    EllipticCurve(const mpz_class& a, const mpz_class& b, const mpz_class& n) : a(a), b(b), n(n) {
        // Calculate discriminant with modular reduction for intermediate results (to prevent growth of values)
        mpz_class discriminant = (((4 * a * a * a) % n) + ((27 * b * b) % n)) % n;  // (4 * a^3 + 27 * b^2) mod n

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

        // Calculate y^2 - (x^3 + ax + b) mod n
        mpz_class equation = ((point.y % n) * (point.y % n)) % n;
        equation -= ((((point.x % n) * (point.x % n)) % n) * (point.x % n)) % n; // x^3 mod n
        equation = (equation - (((a % n) * (point.x % n)) % n - (b % n))) % n;   // y^2 - (x^3 + ax + b) mod n

        // Check if the equation is zero modulo n
        return equation % n == 0;
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
        return (P.x % n == Q.x % n) && (P.y % n == Q.y % n);
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
            throw PointNotOnCurveException("Point P (" + P.x.get_str() + ", " + P.y.get_str() + ") is not on the elliptic curve.");
        }

        if (P.is_infinity || P.y == 0) {
            return ECPoint(); // Result is the point at infinity.
        }

        // Calculate the slope (lambda)
        mpz_class numerator = (3 * P.x * P.x + a) % n;
        mpz_class denominator = (2 * P.y) % n;

        mpz_class gcd;
        if (mpz_invert(denominator.get_mpz_t(), denominator.get_mpz_t(), n.get_mpz_t()) == 0) {
            mpz_gcd(gcd.get_mpz_t(), denominator.get_mpz_t(), n.get_mpz_t());
            throw NonInvertibleElementException(n, 2 * P.y, gcd);
        }

        mpz_class lambda = (numerator * denominator) % n;

        // Calculate resulting coordinates
        mpz_class x3 = (lambda * lambda - 2 * P.x) % n;
        if (x3 < 0) x3 += n;

        mpz_class y3 = (lambda * (P.x - x3) - P.y) % n;
        if (y3 < 0) y3 += n;

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
        // ToDo: checke, ob inverse dann point at infinity ausgeben
        if (!point_is_on_curve(P)) {
            throw PointNotOnCurveException("Point P (" + P.x.get_str() + ", " + P.y.get_str() + ") is not on the elliptic curve.");
        }
        if (!point_is_on_curve(Q)) {
            throw PointNotOnCurveException("Point Q (" + Q.x.get_str() + ", " + Q.y.get_str() + ") is not on the elliptic curve.");
        }
        if (points_are_equal(P, Q)) {
            throw std::invalid_argument("Points P and Q are equal. Use point doubling instead.");
        }

        if (P.is_infinity) return Q;
        if (Q.is_infinity) return P;

        // Calculate the slope (lambda)
        mpz_class numerator = (Q.y - P.y) % n;
        mpz_class denominator = (Q.x - P.x) % n;

        mpz_class gcd;
        if (mpz_invert(denominator.get_mpz_t(), denominator.get_mpz_t(), n.get_mpz_t()) == 0) {
            mpz_gcd(gcd.get_mpz_t(), denominator.get_mpz_t(), n.get_mpz_t());
            throw NonInvertibleElementException(n, Q.x - P.x, gcd);
        }

        mpz_class lambda = (numerator * denominator) % n;

        // Calculate resulting coordinates
        mpz_class x3 = (lambda * lambda - P.x - Q.x) % n;
        if (x3 < 0) x3 += n;

        mpz_class y3 = (lambda * (P.x - x3) - P.y) % n;
        if (y3 < 0) y3 += n;

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
            throw PointNotOnCurveException("Point P (" + P.x.get_str() + ", " + P.y.get_str() + ") is not on the elliptic curve.");
        }

        // Initialize result as the point at infinity (neutral element for addition)
        ECPoint result; // Defaults to the point at infinity
        ECPoint current = P;
        mpz_class scalar = k;

        while (scalar > 0) {
            // Check the least significant bit
            if (scalar % 2 == 1) {
                if (result.is_infinity) {
                    result = current; // no computation needed
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
 * @brief Computes the exponent `e` such that `e = log_p(C + 2*sqrt(C) + 1)`.
 *
 * This function is part of the Lenstra algorithm. The term `C + 2*sqrt(C) + 1` arises from the Hasse Bound (1933).
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
void check_input(const mpz_class& N) {
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
    check_input(N);  // check if N satisfies prerequisites

    // Random values for a,u,v
    mpz_class a, u, v;
    gmp_randstate_t state;
    gmp_randinit_default(state); // Initialize the random state
    gmp_randseed_ui(state, std::random_device{}()); // Seed with a random value

    mpz_urandomm(a.get_mpz_t(),state,N.get_mpz_t());
    mpz_urandomm(u.get_mpz_t(),state,N.get_mpz_t());
    mpz_urandomm(v.get_mpz_t(),state,N.get_mpz_t());

    // Calculate b, such that P(u,v) is on E(a,b,Z/NZ)
    mpz_class const b = (((v * v) % N) - ((u * u * u) % N) - ((a * u) % N) % N); // b = (v^2 - u^3 - au) mod N

    try {
        // Try to initialize elliptic curve
        EllipticCurve E(a,b,N);
        // Elliptic point initialization
        const ECPoint P(u,v);
        ECPoint Q(P.x, P.y);  // Q = P

        // Iterate through primes up to B
        mpz_class prime("2"); // Start with the first prime
        while (prime <= B) {
            mpz_class e = get_exp_for_prime(C,prime);
            // For-loop from 0 up to e-1
            for (mpz_class j = 0; j < e; ++j) {
                std::cout << "Iteration j = " << j << std::endl;  //ToDo: get out
                const ECPoint tmp = E.scalar_multiplication(prime,Q);
                Q = tmp;
            }
            // Dummy statement (replace this with actual logic in the loop)
            std::cout << "Processing prime: " << prime.get_str() << std::endl;  //ToDo: get out
            // Get the next prime
            mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
        }
        // No factor found, algorithm unsuccessful
        throw UnsuccessfulLenstraAlgorithmException("Reached end of lenstra algorithm.");
    } catch (const SingularEllipticCurveException& exc) {
        if (exc.gcd > 1 && exc.gcd < N) {  // we found non-trivial divisor of N; i.e., gdc(((4 * a^3 + 27 * b^2) mod n), n)
            return exc.gcd;
        }
        // Now, we know exc.gcd == N
        throw DiscriminantMultipleOfNException("Discriminant is multiple of N: " + N.get_str());
    } catch (const NonInvertibleElementException& exc) {  // we found an element not invertible over Z/nZ
        return exc.gcd;
    }
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
 * @brief Computes the prime bound `B` based on the smallest prime bound parameter `C`.
 *
 * Uses the formula `B = exp(sqrt(log(C) * log(log(C)) / 2))` to compute an approximate prime bound.
 *
 * @param C The smallest prime bound parameter.
 * @return mpz_class The calculated prime bound as an integer.
 * @throw std::invalid_argument if `C` is less than or equal to 1.
 */
mpz_class calculate_prime_bound_from_smallest_prime_bound(const mpz_class& C) {
    // Ensure C > 1 because ln(ln(C)) is undefined for C <= 1
    if (C <= 1) {
        throw std::invalid_argument("C must be greater than 1.");
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
    if (B < P) {return factors;}
    ppt.exponent = divide_out_maximal_power(N,P);
    if(ppt.exponent > 0)
    {
        ppt.factor = P;
        ppt.is_prime = 2;
        factors.push_back(ppt);

    }
    P++; // P = 3
    if (B < P) {return factors;}
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
        if( a > 0) // Remaining number is probably prime
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
mpz_class generate_get_number_from_user_input(const std::string& number_as_string, const NumberMode number_mode) {
    switch (number_mode) {
        case DirectNumber: return mpz_class(number_as_string);
        case FermatNumber: return Fermat(std::stoi(number_as_string));
        case CunninghamNumber: return Cunningham(std::stoi(number_as_string));
        case TestNumber: return TestzahlB(std::stoi(number_as_string));
        case RSANumber: return RSAZahl(std::stoi(number_as_string));
        default: throw UnknownNumberModeException("Unknown mode " + std::to_string(number_mode) + " provided.");
    }
}


// return k,m, such that N=k^m (and m is maximal, meaning for all a,b, such that N = a^b, it holds m>=b)
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

    // If N = k^m (for natural numbers k,m >1), than m <= floor(log2(N))
    mpz_class k;
    unsigned int upper_bound_exponent = floor(log2(N.get_d()));
    for (unsigned int i = upper_bound_exponent; i >= 2; --i) {
        // test potential exponents, starting with upper bound
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
    int prime_status = mpz_probab_prime_p(k.get_mpz_t(), 27); // GMP primality test
    if (prime_status > 0) { // k is prime or probably prime; add it to factors (with exponent m) and return
        Factor factor;
        factor.exponent = m;
        factor.factor = k;
        factor.is_prime = prime_status;
        factors.push_back(factor);
        N = mpz_class("1");
        return;
    }

    // k is composite and must be further factorized; Info: k cannot be a perfect power itself (otherwise m would've not be maximal)
    std::list<Factor> tmp_factors;
    factorize_with_elliptic_curves(N, B, C, m_curves, tmp_factors);
    // Now that k is completely factorized, update factors (by multiplying all exponents by m)
    for (auto& factor : tmp_factors) {
        factor.exponent *= m; // Multiply each exponent by m
        factors.push_back(factor);    // Add the updated factor to the main list
    }

};


void factorize_with_elliptic_curves(mpz_class N, const mpz_class& B, const mpz_class& C, const int& m_curves, std::list<Factor>& factors) {

    // First, make sure N is not divisible by 2 or 3 (divide out if necessary)
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

    // Check if N is a perfect power, and handle this case
    if (mpz_perfect_power_p(N.get_mpz_t()) > 0) {
        handle_perfect_power(N, B, C, m_curves, factors);
        return; // Stop further factorization since perfect power handling will decompose N
    }

    // Check if N is prime
    int prime_status = mpz_probab_prime_p(N.get_mpz_t(), 27); // GMP primality test
    if (prime_status > 0) { // N is prime or probably prime
        Factor factor;
        factor.exponent = 1;
        factor.factor = N;
        factor.is_prime = prime_status;
        factors.push_back(factor);
        N = mpz_class("1");
        return;
    }

    // Apply Lenstra's algorithm to find a non-trivial divisor
    mpz_class C_tmp = C;
    mpz_class B_tmp = calculate_prime_bound_from_smallest_prime_bound(C_tmp);
    mpz_class divisor("1");
    while (divisor <= 1) {
        try {
            divisor = run_lenstra_algorithm_multiple_times(N, B_tmp, C_tmp, m_curves);
        } catch (const UnsuccessfulLenstraAlgorithmException&) {
        }
        C_tmp = C_tmp * 2; // Double C
        B_tmp = calculate_prime_bound_from_smallest_prime_bound(C_tmp);

    }

    // Check if the divisor is prime
    prime_status = mpz_probab_prime_p(divisor.get_mpz_t(), 27);
    if (prime_status > 0) { // Divisor is prime or probably prime
        // Divide out the maximal power of the prime divisor
        unsigned int exponent = divide_out_maximal_power(N, divisor);

        // Add the prime factor
        Factor factor;
        factor.factor = divisor;
        factor.exponent = exponent;
        factor.is_prime = prime_status; // Definitely prime or probable prime
        factors.push_back(factor);
    } else {
        // Divisor is not prime, further factorize it
        mpz_class divisor_copy = divisor; // Create a copy to avoid modification
        factorize_with_elliptic_curves(divisor_copy, B, C, m_curves, factors);

        // Update N by dividing it by the divisor
        N /= divisor;
    }

    // Factorize the remaining part of N
    factorize_with_elliptic_curves(N, B, C, m_curves, factors);
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
    std::string question("natural number k");

    if(argc == 2)  // one parameter provided
    {
        std::string param(argv[1]);
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
        N = generate_get_number_from_user_input(number_as_string, mode);
    } catch (const UnknownNumberModeException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Valid input received. The program will attempt to factor: " << N << std::endl;

    // Start timing ---------
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();

    // Factorize up to a limit (trial division up to 10^8, adjustable if necessary)
    mpz_class trial_division_bound(TRIAL_DIV_BOUND);  // ToDo: via Variable aus Header file mitgeben
    std::list<Factor> factors = trial_division_bounded(N,trial_division_bound);
    // Use elliptic curves factorization for remainder
    factorize_with_elliptic_curves(N, B, C, NUM_ELLIPTIC_CURVES, factors);

    // Stop timing ----------
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    // Calculate elapsed time and print it
    std::chrono::duration<double, std::milli> float_ms = end - start;
    std::cout << "Factorization in " << float_ms.count() << " milliseconds completed" << std::endl;
    
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


// Potential Improvements:
// - Use faster way to identify all primes less or equal to bound B (e.g., sieve of Atkins), das muss aber abgeglichen
//   werden, da wenn ich sieb verwende, dann berechne ich wirklich alle, ohne dass ich u.U. wirklich alle benötige
// - Efficiency can be improved at several points (but often at the cost of readability). U.a.:
//      - Berechnen der Primzahlen <= B
//      - Handling der perfect powers
// - Könnte multiprocessing nutzen, um parallel mehrere durchläufe für Lenstra zu machen
// Muss wenn ich durchlaufe mal Memory checken (ob der voll läuft)
