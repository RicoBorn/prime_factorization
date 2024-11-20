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

#include "testzahlen.h"

/* compilation: g++ -std=c++11 factorization.cpp testzahlen.cpp -lgmp -lgmpxx */

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
 * @class DiscriminantMultipleOfN
 * @brief Custom exception thrown when discriminant of elliptic curve is multiple of N.
 */
class DiscriminantMultipleOfN : public std::invalid_argument {
public:
    explicit DiscriminantMultipleOfN(const std::string& message)
        : std::invalid_argument(message) {}
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
        // Calculate discriminant with modular reduction at each step (to  minimizes growth of intermediate results)
        mpz_class discriminant = (((4 * a) % n) * (a % n) * (a % n)) % n;  // (4 * a^3) mod n
        const mpz_class tpm_term = (((27 * b) % n) * (b % n)) % n;        // (27 * b^2) mod n
        discriminant = (discriminant + tpm_term) % n;           // (4 * a^3 + 27 * b^2) mod n

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

// ToDo: brauche noch Funktion, die mit gegebenen N einen sinnvollen Startwert für B und C
mpz_class lenstra_algorithm(const mpz_class& N, const mpz_class& B, const mpz_class& C)
{
    if (N % 2 == 0) {
        throw std::invalid_argument("N: " + N.get_str() + "is divisible by 2.");
    }
    if (N % 3 == 0) {
        throw std::invalid_argument("N: " + N.get_str() + "is divisible by 3.");
    }
    if (mpz_perfect_power_p(N.get_mpz_t()) != 0) {
        throw std::invalid_argument("N: " + N.get_str() + "is a perfect power.");
    }

    mpz_class a, u, v, g;

    // Create, initialize, and seed a new random number generator.
    gmp_randstate_t state;
    gmp_randinit_default(state); // Initialize the random state
    gmp_randseed_ui(state, std::random_device{}()); // Seed with a random value

    // Random values for a,u,v
    mpz_urandomm(a.get_mpz_t(),state,N.get_mpz_t());
    mpz_urandomm(u.get_mpz_t(),state,N.get_mpz_t());
    mpz_urandomm(v.get_mpz_t(),state,N.get_mpz_t());

    mpz_class const b = (((v * v) % N) - ((u * u * u) % N) - ((a * u) % N) % N); // b = (v^2 - u^3 - au) mod N
    mpz_class discriminant = (((4 * a * a * a) % N) - ((27 * b * b) % N)) % N;
    mpz_gcd(g.get_mpz_t(), discriminant.get_mpz_t(), N.get_mpz_t());
    if (1 < g < N) {
        return g;
    }

    if (g == N) {
        throw DiscriminantMultipleOfN("Discriminant: " + discriminant.get_str() + "is multiple of N: " + N.get_str());
    }

    // Elliptic curve point initialization
    ECPoint P(u,v);
    ECPoint Q(P.x, P.y);
    mpz_class t("1");

    // Iterate through primes up to B
    mpz_class prime("2"); // Start with the first prime
    while (prime <= B) {
        // evtl. mpf_class nutzen, für log Operationen?

        // Dummy statement (replace this with actual logic in the loop)
        std::cout << "Processing prime: " << prime.get_str() << std::endl;

        // Get the next prime
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    }


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


/**
 * @brief Generates all prime numbers <= N using GMP's mpz_nextprime function.
 * @param N The upper bound (inclusive) for the range of primes, of type mpz_class.
 * @return A vector of mpz_class containing all primes <= N.
 */
std::vector<mpz_class> generate_primes_nextprime(const mpz_class& N) {
    std::vector<mpz_class> primes;

    if (N < 2) {
        return primes; // No primes less than 2
    }

    mpz_class current = 2; // Start from the first prime

    while (current <= N) {
        primes.push_back(current);
        mpz_nextprime(current.get_mpz_t(), current.get_mpz_t()); // Get the next prime
    }

    return primes;
}

/**
 * @brief Generates all prime numbers <= N using the Sieve of Eratosthenes.
 * @param N The upper bound (inclusive) for the range of primes, of type mpz_class.
 * @return A vector of mpz_class containing all primes <= N.
 */
std::vector<mpz_class> generate_primes_sieve(const mpz_class& N) {
    if (N < 2) {
        return {}; // No primes less than 2
    }

    long upper_bound = N.get_si(); // Convert N to long
    std::vector<bool> is_prime(upper_bound + 1, true);
    is_prime[0] = is_prime[1] = false; // 0 and 1 are not primes

    // Sieve of Eratosthenes
    for (long i = 2; i * i <= upper_bound; ++i) {
        if (is_prime[i]) {
            for (long j = i * i; j <= upper_bound; j += i) {
                is_prime[j] = false;
            }
        }
    }

    // Collect primes
    std::vector<mpz_class> primes;
    for (long i = 2; i <= upper_bound; ++i) {
        if (is_prime[i]) {
            primes.push_back(mpz_class(i));
        }
    }

    return primes;
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
    /// Start Testing ////
    ///
    /// compare the two prime functions: geben für ein N die gleiche Anzahl an primes aus?
    /// ACHTUNG: eigentlich brauche ich für den lenstra algorithmus gar nicht alle primes direkt berechenen,
    /// sondern nach und nach ok (kann ja sein, dass bei einem direkt abschmiert..), würde für mpz_nextprime sprechen..
    std::cout  << "Input N:  " ;
    std::string number_as_string2;
    std::cin >> number_as_string2;
    mpz_class N2 = mpz_class(number_as_string2);

    lenstra_algorithm(N2, N2, N2);





    std::vector<mpz_class> z,y;

    // Zeitmessung starten---------
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    z = generate_primes_sieve(N2);
    //Zeitmessung stoppen----------
    std::chrono::time_point<std::chrono::high_resolution_clock>
    end = std::chrono::high_resolution_clock::now();
    //Vergangene Zeitspanne bestimmen und ausgeben
    std::chrono::duration<double, std::milli> float_ms = end - start;
    std::cout << "Berechnung in " << float_ms.count() << " Millisekunden abgeschlossen" << std::endl;
    std::cout << "Number of primes " << z.size() << std::endl;

    // Zeitmessung starten---------
    std::chrono::time_point<std::chrono::high_resolution_clock> start2 = std::chrono::high_resolution_clock::now();
    y = generate_primes_nextprime(N2);
    //Zeitmessung stoppen----------
    std::chrono::time_point<std::chrono::high_resolution_clock>
    end2 = std::chrono::high_resolution_clock::now();
    //Vergangene Zeitspanne bestimmen und ausgeben
    std::chrono::duration<double, std::milli> float_ms2 = end2 - start2;
    std::cout << "Berechnung in " << float_ms2.count() << " Millisekunden abgeschlossen" << std::endl;
    std::cout << "Number of primes " << y.size() << std::endl;

    /////// End Testing ///
    ///
    ///
    ///
    ///
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


    // Factorize up to the limit (trial division up to 10^8, adjustable if necessary)
    mpz_class B("100000000");  // via Variable aus Header file mitgeben
    std::list<Factor> factors = trial_division_bounded(N,B);

    std::cout << "Found the following prime factors:" << std::endl;
    for (std::list<Factor>::iterator it = factors.begin(); it != factors.end(); it++)
        it->printpp();

    std::cout << "Remainder of N: " << N << std::endl;
    // Now use elliptic curve algorithm for remainder

    // TEEEEEEEEEEEST
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

    mpz_class denominator("2");
    mpz_class num("5");
    int output;
    output = mpz_invert(denominator.get_mpz_t(), denominator.get_mpz_t(), num.get_mpz_t());
    std::cout << "denominator: " << denominator << std::endl;
    std::cout << "output: " << output << std::endl;

    return EXIT_SUCCESS; // Program completed successfully
}

// ToDo: dann alles testen was elliptic Curve class angeht!
// Potential Improvements:
// - Use faster way to identify all primes less or equal to bound B (e.g., sieve of Atkins), das muss aber abgeglichen
//   werden, da wenn ich sieb verwende, dann berechne ich wirklich alle, ohne dass ich u.U. wirklich alle benötige
// - Efficiency can be improved at several points (but often at the cost of readability)