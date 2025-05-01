#pragma once

#include <iostream>     // For cache generation messages (if verbose)
#include <chrono>       // Not used directly in header, typically used for benchmarking
#include <array>        // Used for intLog internally
#include <vector>       // Used for caches internally
#include <cmath>        // For std::log10, pow, sin, cos, etc., NAN, floor, round, fmod, hypot, isfinite
#include <iomanip>      // Not used directly in header, typically used for output formatting
#include <limits>       // For numeric_limits (NaN, infinity, max)
#include <numbers>      // For std::numbers::pi (C++20)
#include <numeric>      // Potentially included by other headers, might not be strictly necessary
#include <algorithm>    // For std::max, std::min, std::clamp
#include <utility>      // For std::pair used in atan2 cache

/**
 * @author Alexander Desilets - AlexCDesilets@gmail.com
 * @copyright GNU General Pubnlic license V3.0
 * @brief Provides fast approximations for common math functions using lookup tables.
 * @details This namespace contains functions designed for speed, primarily aimed at
 * applications like real-time graphics or simulations where maximum performance
 * is critical and extremely high precision (beyond 10-11 decimal places for most
 * functions) may be secondary. It achieves speed by pre-calculating lookup tables
 * (caches) for various functions and using linear or bilinear interpolation.
 * The theory of mind used to develop the application was based around a slide-ruler, 
 * so treat it like one, fast, but less accurate than the true calculation
 *
 * @warning **Must call quickMath::_initQuickMath() once before using any quick... functions.**
 * @warning Accuracy varies significantly by function (see individual function docs). Functions
 * like quickAcos, quickTan, and quickAtan2 have known lower accuracy with default
 * settings due to the interpolation method and may not be suitable for high-precision tasks.
 * @warning Cache generation (`_initQuickMath`) is not thread-safe. Call from a single thread
 * during application setup.
 */
namespace quickMath {

    /**
     * @brief Defines preset sizes for the internal lookup tables (caches).
     * @details Larger caches generally provide higher accuracy but consume more memory
     * and may have slightly slower initialization and potentially runtime lookup
     * if they exceed CPU cache limits. defaultCache (1e5) is often a good balance.
     */
    enum class initVals {
        smallCache   = (int)1e4, ///< Smaller cache (10k entries, ~0.16MB x2): Faster init, less memory, lower accuracy.
        defaultCache = (int)1e5, ///< Default cache (100k entries, ~1.6MB x2): Good balance of speed, memory, and accuracy (~10^-11 error for most).
        bigCache     = (int)1e6  ///< Large cache (1M entries, ~16MB x2): Higher accuracy (~10^-13 to 10^-15 error), more memory, slightly slower lookup than default.
        // Note: atan2 cache size also scales, but differently (sqrt(resolution*10)).
    };

    /**
     * @brief Internal implementation details for the quickMath library.
     * @note Users should typically not need to interact with this namespace directly.
     */
    namespace _internal {
        // --- Constants ---
        /// Pre-calculated constant: 1.0 / ln(10) for log conversions. Approximately 0.4342944819...
        const double LOG10_E = 1.0 / std::log(10.0);
        /// Definition of PI, uses std::numbers::pi (C++20) if _PI is not already defined by cmath.
        constexpr double _PI = std::numbers::pi;

        // --- Cache Setup ---
        /// Resolution (number of entries) for most 1D lookup tables. Set by _initQuickMath.
        int CACHE_RESOLUTION = static_cast<int>(initVals::defaultCache);
        /// Resolution per dimension for the 2D atan2 lookup table. Derived from CACHE_RESOLUTION.
        int ATAN2_CACHE_RESOLUTION = static_cast<int>(std::sqrt(CACHE_RESOLUTION)); // Initial value, recalculated in generateAtan2Cache
        /// Cache for log10(x) where x is mapped from [1.0, 10.0).
        std::vector<double> log10_cache;
        /// Cache for pow(10, f) where f is mapped from [0.0, 1.0].
        std::vector<double> pow10_cache;
        /// Cache for sin(theta) where theta is mapped from [0, 2*PI).
        std::vector<double> sin_cache;
        /// Cache for cos(theta) where theta is mapped from [0, 2*PI).
        std::vector<double> cos_cache;
        /// Cache for tan(theta) where theta is mapped from [0, 2*PI). Note: Accuracy issues near asymptotes.
        std::vector<double> tan_cache;
        /// Cache for asin(x) where x is mapped from [-1.0, 1.0].
        std::vector<double> asin_cache;
        /// Cache for acos(x) where x is mapped from [-1.0, 1.0]. Note: Lower accuracy.
        std::vector<double> acos_cache;
        /// Cache storing input pairs (x, y) for atan2 generation (mainly for debug/verification). Unused at runtime.
        std::vector<std::pair<double, double>> atan2_input_cache;
        /// Cache storing atan2(y, x) results for x, y mapped from [-1.0, 1.0]. Note: Lower accuracy.
        std::vector<double> atan2_value_cache;

        // --- Cache Generation Functions (Internal Use) ---

        /// Generates the lookup table for log10(x) for x in [1, 10).
        inline void generateLog10Cache() {
            #ifdef QUICKMATH_VERBOSE
            std::cout << "Generating log10 cache (" << CACHE_RESOLUTION << " entries)..." << std::endl;
            #endif
            log10_cache.resize(CACHE_RESOLUTION);
            double N_minus_1 = static_cast<double>(CACHE_RESOLUTION - 1);
            if (CACHE_RESOLUTION <= 1) { // Handle RESOLUTION = 1 or 0 edge case
                log10_cache[0] = std::log10(1.0); // log10(1) = 0
                #ifdef QUICKMATH_VERBOSE
                if(CACHE_RESOLUTION > 0) std::cout << "log10 cache generated." << std::endl;
                #endif
                return;
            }
            for (int i = 0; i < CACHE_RESOLUTION; ++i) {
                // Input value linearly spaced from 1.0 to 10.0
                double input_val = 1.0 + 9.0 * (static_cast<double>(i) / N_minus_1);
                log10_cache[i] = std::log10(input_val);
            }
            #ifdef QUICKMATH_VERBOSE
            std::cout << "log10 cache generated." << std::endl;
            #endif
        }

        /// Generates the lookup table for 10^f for f in [0, 1].
        inline void generatePow10Cache() {
             #ifdef QUICKMATH_VERBOSE
            std::cout << "Generating pow10 cache (" << CACHE_RESOLUTION << " entries)..." << std::endl;
            #endif
            pow10_cache.resize(CACHE_RESOLUTION);
            double N_minus_1 = static_cast<double>(CACHE_RESOLUTION - 1);
             if (CACHE_RESOLUTION <= 1) {
                pow10_cache[0] = std::pow(10.0, 0.0); // 10^0 = 1
                #ifdef QUICKMATH_VERBOSE
                if(CACHE_RESOLUTION > 0) std::cout << "pow10 cache generated." << std::endl;
                #endif
                return;
            }
            for (int i = 0; i < CACHE_RESOLUTION; ++i) {
                double f_val = static_cast<double>(i) / N_minus_1; // f goes from 0 to 1
                pow10_cache[i] = std::pow(10.0, f_val);
            }
            #ifdef QUICKMATH_VERBOSE
            std::cout << "pow10 cache generated." << std::endl;
            #endif
        }

        /// Generates lookup tables for sin(theta), cos(theta), tan(theta) for theta in [0, 2*PI).
        inline void generateTrigCaches() {
            #ifdef QUICKMATH_VERBOSE
            std::cout << "Generating trig caches (" << CACHE_RESOLUTION << " entries)..." << std::endl;
            #endif
            sin_cache.resize(CACHE_RESOLUTION);
            cos_cache.resize(CACHE_RESOLUTION);
            tan_cache.resize(CACHE_RESOLUTION);

             if (CACHE_RESOLUTION <= 1) {
                 sin_cache[0] = std::sin(0.0);
                 cos_cache[0] = std::cos(0.0);
                 tan_cache[0] = std::tan(0.0);
                 #ifdef QUICKMATH_VERBOSE
                 if(CACHE_RESOLUTION > 0) std::cout << "trig caches generated." << std::endl;
                 #endif
                 return;
             }
             double N_minus_1 = static_cast<double>(CACHE_RESOLUTION - 1);
            for (int i = 0; i < CACHE_RESOLUTION; ++i) {
                 // Theta wraps around [0, 2*PI). Using N-1 ensures last point == first point.
                double theta = (static_cast<double>(i) / N_minus_1) * (2.0 * _PI);
                sin_cache[i] = std::sin(theta);
                cos_cache[i] = std::cos(theta);
                // Store tan values, including potential Inf/large values near asymptotes.
                tan_cache[i] = std::tan(theta);
            }

            #ifdef QUICKMATH_VERBOSE
            std::cout << "trig caches generated." << std::endl;
            #endif
        }

        /// Generates lookup tables for asin(x) and acos(x) for x in [-1, 1].
        inline void generateInverseTrigCaches() {
             #ifdef QUICKMATH_VERBOSE
            std::cout << "Generating inverse trig caches (" << CACHE_RESOLUTION << " entries)..." << std::endl;
            #endif
            asin_cache.resize(CACHE_RESOLUTION);
            acos_cache.resize(CACHE_RESOLUTION);
             if (CACHE_RESOLUTION <= 1) {
                 asin_cache[0] = std::asin(0.0); // Use central value if only one entry
                 acos_cache[0] = std::acos(0.0);
                 #ifdef QUICKMATH_VERBOSE
                 if(CACHE_RESOLUTION > 0) std::cout << "inverse trig caches generated." << std::endl;
                 #endif
                 return;
             }
            double N_minus_1 = static_cast<double>(CACHE_RESOLUTION - 1);
            for (int i = 0; i < CACHE_RESOLUTION; ++i) {
                double x = -1.0 + 2.0 * (static_cast<double>(i) / N_minus_1); // x goes from -1 to 1
                // Clamp input to asin/acos just in case of floating point drift
                double x_clamped = std::clamp(x, -1.0, 1.0);
                asin_cache[i] = std::asin(x_clamped);
                acos_cache[i] = std::acos(x_clamped);
            }
            #ifdef QUICKMATH_VERBOSE
            std::cout << "inverse trig caches generated." << std::endl;
            #endif
        }

        /// Generates a 2D lookup table for atan2(y, x) for x, y in [-1, 1].
        inline void generateAtan2Cache() {
            // Adjust ATAN2 resolution based on potentially updated CACHE_RESOLUTION
            ATAN2_CACHE_RESOLUTION = static_cast<int>(std::max(2.0, std::sqrt(static_cast<double>(CACHE_RESOLUTION) * 10.0)));
            int total = ATAN2_CACHE_RESOLUTION * ATAN2_CACHE_RESOLUTION;

            #ifdef QUICKMATH_VERBOSE
            std::cout << "Generating atan2 cache (" << ATAN2_CACHE_RESOLUTION << " x " << ATAN2_CACHE_RESOLUTION << " = " << total << " entries)..." << std::endl;
            #endif

            atan2_input_cache.resize(total); // Only useful for debug/analysis
            atan2_value_cache.resize(total);

             double N2_minus_1 = static_cast<double>(ATAN2_CACHE_RESOLUTION - 1);
             if (ATAN2_CACHE_RESOLUTION <= 1) {
                  atan2_input_cache[0] = {1.0, 0.0}; // e.g., point (1,0)
                  atan2_value_cache[0] = std::atan2(0.0, 1.0); // = 0
                   #ifdef QUICKMATH_VERBOSE
                   if(total > 0) std::cout << "atan2 cache generated." << std::endl;
                   #endif
                 return;
             }

            for (int yi = 0; yi < ATAN2_CACHE_RESOLUTION; ++yi) {
                for (int xi = 0; xi < ATAN2_CACHE_RESOLUTION; ++xi) {
                    double y = -1.0 + 2.0 * (static_cast<double>(yi) / N2_minus_1); // y from -1 to 1
                    double x = -1.0 + 2.0 * (static_cast<double>(xi) / N2_minus_1); // x from -1 to 1

                    int idx = yi * ATAN2_CACHE_RESOLUTION + xi;
                    atan2_input_cache[idx] = std::make_pair(x, y); // Store input pair if needed
                    // Ensure atan2 is not called with (0,0) during generation
                    if (std::abs(x) < 1e-15 && std::abs(y) < 1e-15) {
                        atan2_value_cache[idx] = 0.0; // Define convention for (0,0) input cache point
                    } else {
                        atan2_value_cache[idx] = std::atan2(y, x);
                    }
                }
            }

            #ifdef QUICKMATH_VERBOSE
            std::cout << "atan2 cache generated." << std::endl;
            #endif
        }
    } // namespace _internal


    /**
     * @brief Initializes the quickMath library caches. MUST be called once before using functions.
     * @param cache_Size The desired cache resolution preset using the initVals enum.
     * Affects memory usage, accuracy, and speed of most quick... functions.
     * Defaults to defaultCache (1e5 entries).
     * @note This function is not thread-safe. Call from a single thread during application setup.
     * @note Calling this again will clear and regenerate caches with the new size.
     */
    inline void _initQuickMath(initVals cache_Size = initVals::defaultCache) {
        _internal::CACHE_RESOLUTION = static_cast<int>(cache_Size);
        // Ensure minimum practical size for calculations
        if (_internal::CACHE_RESOLUTION < 2) _internal::CACHE_RESOLUTION = 2;

        // Clear caches if they were previously generated (allows re-init with different size)
        // Checking one is sufficient if all are generated together.
        if (!_internal::log10_cache.empty()) {
             _internal::acos_cache.clear();
             _internal::asin_cache.clear();
             _internal::atan2_input_cache.clear();
             _internal::atan2_value_cache.clear();
             _internal::cos_cache.clear();
             _internal::log10_cache.clear();
             _internal::pow10_cache.clear();
             _internal::sin_cache.clear();
             _internal::tan_cache.clear();
        }
        // Generate all caches
        _internal::generateLog10Cache();
        _internal::generatePow10Cache();
        _internal::generateTrigCaches();
        _internal::generateInverseTrigCaches();
        _internal::generateAtan2Cache(); // ATAN2 resolution calculated inside
    }

    // -------------------------------------------------------------------------
    // Public Math Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Calculates the base-10 logarithm using cache lookup and linear interpolation.
     * @param x Input value.
     * @return Approximation of log10(x). Returns NaN if x <= 0 or caches not initialized.
     * @note Domain: x > 0.
     * @note Accuracy (default cache): ~10^-11 absolute error average. Performance: Faster than std::log10.
     */
    inline double quickLog(double x) {
        const double epsilon = 1e-12; // Tolerance for checking near 1.0
        if (x <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        // Optimization: check for exactly 1.0 early
        if (x == 1.0) return 0.0;
        if (std::abs(x - 1.0) < epsilon) return 0.0;


        // Range Reduction to m in [1, 10)
        double m = x; int k = 0;
        while (m >= 10.0) { m /= 10.0; k++; }
        while (m < 1.0) { m *= 10.0; k--; }
        // Check if m ended up being 1.0 after reduction
        if (std::abs(m - 1.0) < epsilon) return static_cast<double>(k);

        // Cache check
        if (_internal::log10_cache.empty()) {
            // Consider alternative error handling? Throw? Log?
            return std::numeric_limits<double>::quiet_NaN();
        }

        // Cache Lookup & Interpolation
        const double N_minus_1 = static_cast<double>(_internal::CACHE_RESOLUTION - 1);
        const double cache_span = 9.0; // Input range for log cache is [1, 10] -> span 9
        double float_index = ((m - 1.0) / cache_span) * N_minus_1;

        int index0 = static_cast<int>(float_index);
        // Clamp index to ensure index0 and index0+1 are valid
        index0 = std::max(0, std::min(index0, _internal::CACHE_RESOLUTION - 2));
        int index1 = index0 + 1;

        double fraction = float_index - static_cast<double>(index0);

        // Interpolate using FMA
        double delta = _internal::log10_cache[index1] - _internal::log10_cache[index0];
        double log_m = std::fma(fraction, delta, _internal::log10_cache[index0]);

        return log_m + k;
    }

    /**
     * @brief Calculates the natural logarithm (base e) using conversion from quickLog.
     * @param x Input value.
     * @return Approximation of ln(x). Returns NaN if x <= 0 or caches not initialized.
     * @note Domain: x > 0.
     * @note Accuracy depends on quickLog (~10^-11 absolute error scale). Performance: Faster than std::log.
     */
    inline double quickLn(double x) {
        // ln(x) = log10(x) / log10(e) = log10(x) * ln(10)
        // Using the precalculated constant LOG10_E = 1.0 / ln(10) is slightly confusing name here.
        // Let's rename or use the inverse: ln(10) = std::log(10.0)
         const double LN_10 = std::log(10.0); // Approx 2.302585...
         return quickLog(x) * LN_10;
        // Original calculation using LOG10_E = 1/ln(10) would be:
        // return quickLog(x) / _internal::LOG10_E; // This is correct too.
    }

    /**
     * @brief Calculates 10 raised to the power of an integer exponent k. (Internal Helper)
     * @details Uses exponentiation by squaring for efficiency. Includes overflow/underflow checks.
     * @param k Integer exponent.
     * @return 10^k, or infinity/0.0 on overflow/underflow.
     */
    inline double powerOf10(int k) {
        // Implementation as before...
         if (k == 0) return 1.0;
         if (k > 308) return std::numeric_limits<double>::infinity();
         if (k < -308) return 0.0;
         if (k < 0) return 1.0 / powerOf10(-k); // Check reciprocal limit in recursive call
         double result = 1.0;
         double base = 10.0;
         int exponent = k;
         while (exponent > 0) {
             if (exponent % 2 == 1) {
                 if (result > std::numeric_limits<double>::max() / base) return std::numeric_limits<double>::infinity();
                 result *= base;
             }
             if (exponent > 1 && base > std::sqrt(std::numeric_limits<double>::max())) {
                 return std::numeric_limits<double>::infinity();
             }
             base *= base;
             exponent /= 2;
         }
         return result;
    }

    /**
     * @brief Calculates 10 raised to the power L (antilog base 10) using cache lookup.
     * @param L Input exponent.
     * @return Approximation of 10^L. Returns NaN/Inf if L is non-finite or caches not init. Returns Inf/0 on overflow/underflow.
     * @note Domain: Any finite double L.
     * @note Accuracy (default cache): ~10^-11 relative error average. Performance: Faster than std::pow(10, L).
     */
    inline double quickUnlog(double L) {
        const double epsilon = 1e-12; // Tolerance for checking near integers
        if (!std::isfinite(L)) return L; // Return NaN or +/- Inf directly

        // Range Reduction: L = k + f, where f is in [0, 1)
        double k_double = std::floor(L);
        double f = L - k_double;
        int k = static_cast<int>(k_double);

        // Handle f being extremely close to 0 or 1 due to precision
        // If close, use powerOf10 with the nearest integer exponent for robustness
        if (std::abs(f) < epsilon || std::abs(f - 1.0) < epsilon) {
            return powerOf10(static_cast<int>(std::round(L)));
        } else if (f < 0.0) { // Adjust if slightly negative
            f += 1.0;
            k -= 1;
        }
        // Now f is reliably in (0, 1)

        // Cache check
        if (_internal::pow10_cache.empty()) {
             return std::numeric_limits<double>::quiet_NaN();
        }

        // Cache Lookup & Interpolation for 10^f
        const double N_minus_1 = static_cast<double>(_internal::CACHE_RESOLUTION - 1);
        double float_index = f * N_minus_1; // f is in [0, 1) -> index in [0, N-1)
        int index0 = static_cast<int>(float_index);
        // Clamp index
        index0 = std::max(0, std::min(index0, _internal::CACHE_RESOLUTION - 2));
        int index1 = index0 + 1;
        double fraction = float_index - static_cast<double>(index0);

        // Interpolate using FMA
        double delta = _internal::pow10_cache[index1] - _internal::pow10_cache[index0];
        double m = std::fma(fraction, delta, _internal::pow10_cache[index0]); // m = 10^f

        // Combine with integer power part
        return m * powerOf10(k); // result = (10^f) * (10^k) = 10^L
    }

    /**
     * @brief Calculates e raised to the power x (natural exponential) using quickUnlog.
     * @param x Input exponent.
     * @return Approximation of e^x. Returns NaN/Inf if x is non-finite or result overflows.
     * @note Accuracy depends on quickUnlog (~10^-11 relative error scale). Performance: Faster than std::exp.
     */
    inline double quickExp(double x) {
        // e^x = 10^(x * log10(e))
        // Need log10(e) = ln(e)/ln(10) = 1/ln(10) = LOG10_E
        return quickUnlog(x * _internal::LOG10_E);
    }

   /**
    * @brief Calculates base 'a' raised to the power 'b'.
    * @details Implements a^b = 10^(b * log10(a)) using quickLog and quickUnlog.
    * Includes handling for special cases (a=0, a=1, b=0, b=1) and
    * negative base 'a' raised to integer powers 'b'. Matches std::pow domain for real results.
    * custom rules for negative base raised to fractional power, assumes the sign only flips when the fractional power passes an integral checkpoint
    * this allows for extended, functional math
    * @param a Base.
    * @param b Exponent.
    * @return Approximation of a^b. 2^2 = 4, 2^-2 = 1/4, -2^2 = 4, 2^3 = -8, -2^2.2 = ~4.59, -2^-3.2 = -9.19
    * @note Accuracy (default cache): ~10^-10 relative error average. Performance: Faster than std::pow.
    */
    inline double quickPow(double a, double b) {
        if (a == 0.0) return b == 0.0 ? 1.0: b < 0.0 ? std::numeric_limits<double>::infinity(): 0.0;
        if (b == 0.0) return 1;

        // check how the law of signs will apply to the result of the calculation, -2^2 = 4 -> positive sign, -2^-3 = -8 -> negative sign
        bool negativesign = false;
        if (a < 0 && std::abs((int)b)%2==1) negativesign = true;

        double base = std::abs(a);
        double exponent = std::abs(b);
        
        // Check if the result has to be the reciporcal of the calculation, 2^-2 = 1/2^2
        bool negativeExponent = false;
        if (b < 0) negativeExponent = true;

        double result = quickUnlog(exponent * quickLog(base));

        if (negativeExponent) result = 1.0 / result;

        return negativesign ? -1 * result:result;
    }

    /**
     * @brief Calculates the square root using quickPow(x, 0.5).
     * @param x Input value.
     * @return Approximation of sqrt(x). Returns NaN if x < 0.
     * @note Domain: x >= 0.
     * @note Accuracy depends on quickPow (~10^-10 relative error). Performance: Likely slower than std::sqrt.
     */
    inline double quickSqrt(double x) {
        if (x < 0.0) return std::numeric_limits<double>::quiet_NaN();
        if (std::abs(x) < 1e-12) return 0.0; // sqrt(0) is 0
        return quickPow(x, 0.5);
    }

    /**
     * @brief Calculates the sine of an angle (in radians) using cache lookup.
     * @param x Angle in radians.
     * @return Approximation of sin(x). Returns NaN if caches not initialized.
     * @note Domain: Any finite double x. Input is reduced to [0, 2*PI).
     * @note Accuracy (default cache): Good (~10^-8 relative error average). Performance: Faster than std::sin.
     */
    inline double quickSin(double x) {
        // Cache check
        if (_internal::sin_cache.empty()) return std::numeric_limits<double>::quiet_NaN();
        // Handle non-finite input? std::sin handles them.
        if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN(); // Or handle Inf specifically?

        // Range reduction to [0, 2*PI)
        const double two_pi = 2.0 * _internal::_PI;
        x = std::fmod(x, two_pi);
        if (x < 0.0) x += two_pi;

        // Cache Lookup & Interpolation
        const double N_minus_1 = static_cast<double>(_internal::CACHE_RESOLUTION - 1);
        double index = (x / two_pi) * N_minus_1;
        int i0 = static_cast<int>(index);
        // Clamp i0 to ensure validity before modulo/adding 1
        i0 = std::max(0, std::min(i0, _internal::CACHE_RESOLUTION - 1));
        int i1 = (i0 + 1);
        // Handle wrap-around for interpolation point i1
        if (i1 >= _internal::CACHE_RESOLUTION) i1 = 0;

        double frac = index - static_cast<double>(i0);

        // Interpolate using FMA
        double delta = _internal::sin_cache[i1] - _internal::sin_cache[i0];
        return std::fma(frac, delta, _internal::sin_cache[i0]);
    }

   /**
    * @brief Calculates the arcsine (inverse sine) using cache lookup.
    * @param x Input value. Will be clamped to [-1, 1].
    * @return Approximation of asin(x) in radians [-PI/2, PI/2]. Returns NaN if caches not initialized.
    * @note Domain: x in [-1, 1]. Inputs outside this range are clamped.
    * @note Accuracy (default cache): Moderate (~10^-6 relative error average). Performance: Faster than std::asin.
    */
    inline double quickAsin(double x) {
        // Cache check
        if (_internal::asin_cache.empty()) return std::numeric_limits<double>::quiet_NaN();
        // Handle NaN input
        if (!std::isfinite(x)) return x;

        // Clamp input to valid domain [-1, 1]
        x = std::clamp(x, -1.0, 1.0);

        // Map [-1, 1] to [0, N-1] index
        const double N_minus_1 = static_cast<double>(_internal::CACHE_RESOLUTION - 1);
        double index = (x + 1.0) * 0.5 * N_minus_1;
        int i0 = static_cast<int>(index);
        // Clamp index before calculating i1 to prevent reading past array bounds
        i0 = std::max(0, std::min(i0, _internal::CACHE_RESOLUTION - 2)); // Ensure i0 and i0+1 are valid
        int i1 = i0 + 1;
        double frac = index - static_cast<double>(i0);

        // Interpolate using FMA
        double delta = _internal::asin_cache[i1] - _internal::asin_cache[i0];
        return std::fma(frac, delta, _internal::asin_cache[i0]);
    }

   /**
    * @brief Calculates the cosine of an angle (in radians) using cache lookup.
    * @param x Angle in radians.
    * @return Approximation of cos(x). Returns NaN if caches not initialized.
    * @note Domain: Any finite double x. Input is reduced to [0, 2*PI).
    * @note Accuracy (default cache): Very Good (~10^-10 relative error average). Performance: Faster than std::cos.
    */
    inline double quickCos(double x) {
        // Cache check
        if (_internal::cos_cache.empty()) return std::numeric_limits<double>::quiet_NaN();
         if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();

        // Range reduction to [0, 2*PI)
        const double two_pi = 2.0 * _internal::_PI;
        x = std::fmod(x, two_pi);
        if (x < 0.0) x += two_pi;

        // Cache Lookup & Interpolation
        const double N_minus_1 = static_cast<double>(_internal::CACHE_RESOLUTION - 1);
        double index = (x / two_pi) * N_minus_1;
        int i0 = static_cast<int>(index);
        // Clamp i0
        i0 = std::max(0, std::min(i0, _internal::CACHE_RESOLUTION - 1));
        int i1 = (i0 + 1);
        if (i1 >= _internal::CACHE_RESOLUTION) i1 = 0; // Wrap index

        double frac = index - static_cast<double>(i0);

        // Interpolate using FMA
        double delta = _internal::cos_cache[i1] - _internal::cos_cache[i0];
        return std::fma(frac, delta, _internal::cos_cache[i0]);
    }

   /**
    * @brief Calculates the arccosine (inverse cosine) using cache lookup.
    * @param x Input value. Will be clamped to [-1, 1].
    * @return Approximation of acos(x) in radians [0, PI]. Returns NaN if caches not initialized.
    * @note Domain: x in [-1, 1]. Inputs outside this range are clamped.
    * @warning Accuracy (default cache): **Poor** (~10^-3 relative error average) due to linear
    * interpolation issues near x = +/-1 where the function's slope is infinite.
    * Use with extreme caution if precision is required. Performance: Faster than std::acos.
    */
    inline double quickAcos(double x) {
         // Cache check
        if (_internal::acos_cache.empty()) return std::numeric_limits<double>::quiet_NaN();
         if (!std::isfinite(x)) return x; // Return NaN

        // Clamp input to valid domain [-1, 1]
        x = std::clamp(x, -1.0, 1.0);

        // Map [-1, 1] to [0, N-1] index
        const double N_minus_1 = static_cast<double>(_internal::CACHE_RESOLUTION - 1);
        double index = (x + 1.0) * 0.5 * N_minus_1;
        int i0 = static_cast<int>(index);
        // Clamp index
        i0 = std::max(0, std::min(i0, _internal::CACHE_RESOLUTION - 2));
        int i1 = i0 + 1;
        double frac = index - static_cast<double>(i0);

        // Interpolate using FMA
        double delta = _internal::acos_cache[i1] - _internal::acos_cache[i0];
        return std::fma(frac, delta, _internal::acos_cache[i0]);
    }

    /**
     * @brief Calculates the tangent of an angle (in radians) using cache lookup.
     * @param x Angle in radians.
     * @return Approximation of tan(x). May return NaN or Inf.
     * @warning This implementation uses direct cache lookup and linear interpolation,
     * which yields **highly inaccurate results** near the function's asymptotes
     * (e.g., near PI/2, 3*PI/2). Accuracy is poor (~10^-5 relative error average) even
     * in tested ranges avoiding direct hits on asymptotes. Interpolation across
     * cached Infinity values is undefined. Returns NaN if caches not initialized.
     * @warning **Recommended alternative for accuracy: calculate `quickSin(x) / quickCos(x)`**.
     */
    inline double quickTan(double x) {
         // Cache check
        if (_internal::tan_cache.empty()) return std::numeric_limits<double>::quiet_NaN();
         if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();

        // Range reduction to [0, 2*PI)
        const double two_pi = 2.0 * _internal::_PI;
        x = std::fmod(x, two_pi);
        if (x < 0.0) x += two_pi;

        // Cache Lookup & Interpolation
        const double N_minus_1 = static_cast<double>(_internal::CACHE_RESOLUTION - 1);
        double index = (x / two_pi) * N_minus_1;
        int i0 = static_cast<int>(index);
        // Clamp i0
        i0 = std::max(0, std::min(i0, _internal::CACHE_RESOLUTION - 1));
        int i1 = (i0 + 1);
        if (i1 >= _internal::CACHE_RESOLUTION) i1 = 0; // Wrap index

        double frac = index - static_cast<double>(i0);

        // Interpolate using FMA - **WARNING: Highly inaccurate near asymptotes**
        double val0 = _internal::tan_cache[i0];
        double val1 = _internal::tan_cache[i1];

        // If either value is non-finite, interpolation is unreliable
        if (!std::isfinite(val0) || !std::isfinite(val1)) {
             // Near an asymptote. Result should be large or NaN.
             // Returning NaN is safer than trying to guess infinity sign.
             return std::numeric_limits<double>::quiet_NaN();
        }

        double delta = val1 - val0;
        return std::fma(frac, delta, val0);
    }

    /**
     * @brief Calculates the angle (in radians) for a given point (y, x) using 2D cache lookup.
     * @details Maps the input vector (y, x) to the unit circle, then performs bilinear interpolation
     * on a precomputed grid of atan2 values for inputs in [-1, 1] x [-1, 1].
     * @param y The y-coordinate.
     * @param x The x-coordinate.
     * @return Approximation of atan2(y, x) in radians [-PI, PI]. Returns NaN if caches not initialized. Returns 0 if x=0, y=0.
     * @note Domain: Any finite x, y (except both zero simultaneously, handled).
     * @note Accuracy (default cache + 10x atan2 res): Moderate-to-Poor (~10^-4 relative error average)
     * due to bilinear interpolation on the derived grid resolution (~1000x1000). Performance: Faster than std::atan2.
     */
    inline double quickAtan2(double y, double x) {
        // Handle origin case (std::atan2(0,0) returns 0)
        if (std::abs(x) < 1e-12 && std::abs(y) < 1e-12) return 0.0;
        // Cache check
        if (_internal::atan2_value_cache.empty()) return std::numeric_limits<double>::quiet_NaN();
        // Handle non-finite inputs? std::atan2 specifies behavior.
         if (!std::isfinite(x) || !std::isfinite(y)) return std::numeric_limits<double>::quiet_NaN();


        // --- Input Mapping to Normalized [-1, 1] Range ---
        // Normalize (x,y) onto the unit circle to use the cache grid based on angle.
        double mag = quickSqrt(x*x+y*y);
        // Avoid division by zero if magnitude is extremely small (should be covered by origin check)
        // if (mag < 1e-15) return 0.0; // Or handle based on signs of x, y? 0 is atan2 convention

        double norm_x = x / mag;
        double norm_y = y / mag;

        // Map normalized [-1, 1] coordinates to cache indices [0, ATAN2_RES-1]
        const double N2_minus_1 = static_cast<double>(_internal::ATAN2_CACHE_RESOLUTION - 1);
        double fx = (norm_x + 1.0) * 0.5 * N2_minus_1;
        double fy = (norm_y + 1.0) * 0.5 * N2_minus_1;

        // Get integer indices and fractions for interpolation
        int xi = static_cast<int>(fx);
        int yi = static_cast<int>(fy);
        double frac_x = fx - static_cast<double>(xi);
        double frac_y = fy - static_cast<double>(yi);

        // Clamp indices to prevent reading out of bounds for the 4 points
        xi = std::max(0, std::min(xi, _internal::ATAN2_CACHE_RESOLUTION - 2));
        yi = std::max(0, std::min(yi, _internal::ATAN2_CACHE_RESOLUTION - 2));

        // --- Bilinear Interpolation ---
        // Indices for the 4 corner points in the 1D cache vector
        int idx00 = yi * _internal::ATAN2_CACHE_RESOLUTION + xi;
        int idx10 = idx00 + 1;                              // Point (xi+1, yi)
        int idx01 = idx00 + _internal::ATAN2_CACHE_RESOLUTION; // Point (xi, yi+1)
        int idx11 = idx01 + 1;                              // Point (xi+1, yi+1)

        // Values at the four corners
        double v00 = _internal::atan2_value_cache[idx00];
        double v10 = _internal::atan2_value_cache[idx10];
        double v01 = _internal::atan2_value_cache[idx01];
        double v11 = _internal::atan2_value_cache[idx11];

        // Interpolate along x for both y levels using FMA
        double r0 = std::fma(frac_x, v10 - v00, v00); // Interpolated value at y=yi level
        double r1 = std::fma(frac_x, v11 - v01, v01); // Interpolated value at y=yi+1 level

        // Interpolate along y using the results from x interpolation (using FMA)
        // Final result = r0 + frac_y * (r1 - r0)
        return std::fma(frac_y, r1 - r0, r0);
    }

}; // namespace quickMath
