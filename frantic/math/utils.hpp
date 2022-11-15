// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cassert>
#include <stdexcept>
#include <vector>

#ifdef _WIN32
#include <intrin.h>
#endif

#include <boost/cstdint.hpp>

#include <boost/lexical_cast.hpp>

#include <boost/math/special_functions/fpclassify.hpp>

// File: utils.hpp
// Contains a set of useful mathematical values and tools.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>

namespace frantic {
namespace math {

// const: pi_value
// The value of pi to 15 decimal places
const double pi_value = 3.141592653589793;

// clamp a value between an upper and lower bound.
template <typename T>
inline T clamp( T val, T min_val, T max_val ) {
    if( val < min_val )
        return min_val;
    else if( val > max_val )
        return max_val;
    return val;
    // return (std::min)(max_val, (std::max)(val, min_val));
}

// Function: radians_to_degrees
// Converts the given angle in radians to an angle in degrees.
//
// Parameter:
// radians - the angle in radians to convert to degrees.
//
// Returns:
// A float representing the angle in degrees.
const double degrees_per_radian = 180.0 / pi_value;
inline float radians_to_degrees( float radians ) { return radians * static_cast<float>( degrees_per_radian ); }
inline double radians_to_degrees( double radians ) { return radians * degrees_per_radian; }

// Function: degrees_to_radians
// Converts the given angle in degrees to an angle in radians.
//
// Parameter:
// degrees - the angle in degrees to convert to radians.
//
// Returns:
// A float representing the angle in radians.
const double radians_per_degree = pi_value / 180.0;
inline float degrees_to_radians( float degrees ) { return degrees * static_cast<float>( radians_per_degree ); }
inline double degrees_to_radians( double degrees ) { return degrees * radians_per_degree; }

template <typename T>
inline T get_absolute( T value ) {
    if( value < T( 0 ) )
        return -value;
    else
        return value;
}

// Function: lerp
// Finds the point that is linearly interpolated between
// points 'a' and 'b' at time 't'.
//
// Parameters:
// a - the starting point
// b - the ending point
// t - the time to do the interpolation
//
// Returns:
// A value which will be linear interpolation between
// points 'a' and 'b' at time 't'.
template <typename T>
inline T lerp( T a, T b, float t ) {
    return ( 1.0f - t ) * a + t * b;
}

// Mix is the same as above with alpha clamped between zero and one.
template <typename T>
inline T mix( T x, T y, double alpha ) {
    alpha = clamp<double>( alpha, 0.0, 1.0 );
    return T( ( x * alpha ) + ( y * ( 1.0 - alpha ) ) );
}

/**
 * A relative error version of a "less than" comparison.
 * \param lhs The left hand side of the comparison
 * \param rhs The right hand side of the comparison
 * \param tolerance The pseudo-relative error tolerance to use in the comparison.
 * \return True if 'lhs' is less that 'rhs' in a relative sense.
 */
template <typename T>
inline bool relative_less( T lhs, T rhs, T tolerance ) {
    return rhs - lhs > tolerance * ( static_cast<T>( 1 ) + std::abs( rhs ) ); // NOTE: (1 + |rhs|) is a fancy way to get
                                                                              // relative error at large scale and
                                                                              // absolute error at small scale)
}

template <typename T>
inline bool relative_less_equal( T lhs, T rhs, T tolerance ) {
    return lhs - rhs <= tolerance * ( static_cast<T>( 1 ) + std::abs( rhs ) ); // NOTE: (1 + |rhs|) is a fancy way to
                                                                               // get relative error at large scale and
                                                                               // absolute error at small scale)
}

template <typename T>
inline bool relative_greater( T lhs, T rhs, T tolerance ) {
    return lhs - rhs > tolerance * ( static_cast<T>( 1 ) + std::abs( rhs ) ); // NOTE: (1 + |rhs|) is a fancy way to get
                                                                              // relative error at large scale and
                                                                              // absolute error at small scale)
}

template <typename T>
inline bool relative_greater_equal( T lhs, T rhs, T tolerance ) {
    return rhs - lhs <= tolerance * ( static_cast<T>( 1 ) + std::abs( rhs ) ); // NOTE: (1 + |rhs|) is a fancy way to
                                                                               // get relative error at large scale and
                                                                               // absolute error at small scale)
}

template <typename T>
inline bool relative_equal_to( T lhs, T rhs, T tolerance ) {
    return std::abs( lhs - rhs ) <=
           tolerance * ( static_cast<T>( 1 ) + std::abs( rhs ) ); // NOTE: (1 + |rhs|) is a fancy
                                                                  // way to get relative error at
                                                                  // large scale and absolute
                                                                  // error at small scale)
}

template <typename T>
inline bool relative_not_equal_to( T lhs, T rhs, T tolerance ) {
    return std::abs( lhs - rhs ) >
           tolerance * ( static_cast<T>( 1 ) + std::abs( rhs ) ); // NOTE: (1 + |rhs|) is a fancy
                                                                  // way to get relative error at
                                                                  // large scale and absolute
                                                                  // error at small scale)
}

// Function: round
// Rounds the given float to the nearest whole number.
//
// Parameters:
// value - the float to round.
//
// Returns:
// The specified float rounded to the nearest whole number (float).
inline float round( float value ) { return floor( value + 0.5f ); }

// Function: round
// Rounds the given double to the nearest whole number.
//
// Parameters:
// value - the double to round.
//
// Returns:
// The specified double rounded to the nearest whole number (double).
inline double round( double value ) { return floor( value + 0.5 ); }

// converts a value from one linear range to another.
// it is the programmer's responsibility to ensure that
template <typename T, typename U>
inline U linearConvertRange( T value, T fromLB, T fromUB, U toLB, U toUB ) {
    assert( fromUB != fromLB );

    T offset = value - fromLB;
    T fromRange = fromUB - fromLB;

    T toRange = (T)( toUB - toLB );

    return toLB + (U)( ( offset * toRange ) / fromRange );
}

template <typename T>
inline T smoothstep( T val, T lb, T ub ) {
    // assert(lb != ub);
    if( val >= ub )
        return static_cast<T>( 1 );
    else if( val <= lb )
        return static_cast<T>( 0 );
    else {
        T x = ( val - lb ) / ( ub - lb );
        return x * x * ( static_cast<T>( 3 ) - ( static_cast<T>( 2 ) * x ) );
    }
}

template <class T>
inline T square( const T& x ) {
    return x * x;
}

//
// Function: cotangent
//
template <typename T>
T cot( const T& value ) {
    return static_cast<T>( std::tan( pi_value / 2 - value ) );
}

// Calculates (a*b)%c in a manner resistant to overflow. C must be < 2^63 since 2*(a%c)
// must be < 2^64 (The largest unsigned 64bit integer.)
//
// Since (a*b)%c = ((a%c)*(b%c))%c, and (b%c) = (b63*2^63 + b62*2^62 + ... + b0)%c where bN is the Nth binary digit of
// b.
// Using Horner's algorithm:
//   b[0] = "63rd bit of b"
//   b[i+1] = b[i]*2 + "(63-i)'th bit of b"
// Then b[63] == b.  Multiply a*b, you get
//   ab[0] = a * "63rd bit of b"
//   ab[i+1] = ab[i]*2 + a * "(63-i)'th bit of b"
// Then ab[63] == a*b.  Now, we want to compute this modulo c, so
//   ab_c[0] = (a % c) * "63rd bit of b"
//   ab_c[i+1] = (ab_c[i]*2 + a * "(63-i)'th bit of b") % c
// Then ab_c[63] == (a*b) % c.  More conservatively, the iterative step could be
//   ab_c[i+1] = ((ab_c[i]*2) % c + (a % c) * "(63-i)'th bit of b") % c
inline boost::int64_t safe_mul_mod( boost::uint64_t a, boost::uint64_t b, boost::uint64_t c ) {
    if( c >= ( static_cast<boost::uint64_t>( 1 ) << 63 ) )
        throw std::runtime_error( "safe_mul_mod(): The divisor " + boost::lexical_cast<std::string>( c ) +
                                  " is too large" );

    boost::uint64_t result = 0;
    boost::uint64_t aModC = a % c;
    boost::uint64_t mask = ( static_cast<boost::uint64_t>( 1 ) << 63 );
    // Find the first non-zero bit
    while( !( b & mask ) && mask != 0 )
        mask >>= 1;
    // Use Horner's rule/algorithm/scheme/method to evaluate the expression
    for( ; mask != 0; mask >>= 1 ) {
        result = ( result << 1 ) % c;
        if( b & mask )
            result = ( result + aModC ) % c;
    }

    return result;
}

inline void byte_swap_4bytes( char* data ) {
    std::swap( data[0], data[3] );
    std::swap( data[1], data[2] );
}

inline void byte_swap_2bytes( char* data ) { std::swap( data[0], data[1] ); }

inline void byte_swap( float& data ) { byte_swap_4bytes( (char*)&data ); }

inline void byte_swap( boost::uint32_t& data ) { byte_swap_4bytes( (char*)&data ); }

inline void byte_swap( boost::uint16_t& data ) { byte_swap_2bytes( (char*)&data ); }

inline float mm_to_inch( float v ) { return v * 0.0393700787f; }
inline float inch_to_mm( float v ) { return v * 25.4f; }

enum LengthUnit { MILLIMETERS, METERS, INCHES, FEET };

inline float convert_to_mm( float dist, const LengthUnit& unitType ) {
    switch( unitType ) {
    case FEET:
        return frantic::math::inch_to_mm( dist * 12.0f );
    case INCHES:
        return frantic::math::inch_to_mm( dist );
    case METERS:
        return dist * 1000.0f;
    case MILLIMETERS:
        return dist;
    default:
        throw std::runtime_error( "convert_to_mm failed: Attempted to convert unknown LengthUnit" );
    }
}

inline bool is_finite( float v ) { return ( boost::math::isfinite )( v ); }

#ifdef _MSC_VER
inline bool is_finite( double v ) { return _finite( v ) != 0; }
inline bool is_infinite( double v ) {
    double inf = std::numeric_limits<double>::infinity();
    return v == inf || v == -inf;
}
inline bool is_nan( double v ) { return _isnan( v ) != 0; }
#else
inline bool is_finite( double v ) { return !std::isinf( v ) && !std::isnan( v ); }
inline bool is_infinite( double v ) { return std::isinf( v ); }
inline bool is_nan( double v ) { return std::isnan( v ); }
#endif

template <class T>
void write_scilab_vector( std::ostream& out, const std::vector<T>& vector, const std::string& varName,
                          bool columnVector ) {
    out << varName << " = [";
    for( unsigned i = 0; i < vector.size(); ++i ) {
        out << vector[i];
        if( columnVector )
            out << ";";
        else
            out << ",";
    }
    out << "];\n" << std::endl;
}

template <typename T, int N>
bool solve_by_elimination( T mat[N * N], T result[N] ) {
    T* rows[N]; // Used for row swapping

    for( int i = 0; i < N; i++ )
        rows[i] = &mat[i * N];

    for( int i = 0; i < N; i++ ) {
        // Find the largest possible pivot point
        T maxValue = fabs( rows[i][i] );
        int maxIndex = i;

        for( int j = i + 1; j < N; j++ ) {
            T val = fabs( rows[j][i] );
            if( val > maxValue ) {
                maxValue = val;
                maxIndex = j;
            }
        }

        // Set the largest to the pivot point, or return false if the matrix is degenerate
        if( maxValue == 0 )
            return false;

        if( maxIndex != i ) {
            T* m = rows[i];
            T r = result[i];

            rows[i] = rows[maxIndex];
            result[i] = result[maxIndex];

            rows[maxIndex] = m;
            result[maxIndex] = r;
        }

        result[i] /= rows[i][i];
        for( int k = N - 1; k >= i; k-- )
            rows[i][k] /= rows[i][i];

        // Inspect each pivot
        for( int j = 0; j < N; j++ ) {
            if( j != i ) {
                T factor = rows[j][i];
                result[j] -= result[i] * factor;
                for( int k = i; k < N; k++ )
                    rows[j][k] -= rows[i][k] * factor;
            }
        }
    }

    return true;
}

inline unsigned int round_up_to_power_of_two( unsigned int x ) {
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x++;

    return x;
}

/**
 * Return the log2 of an integer.  This is equivalent to the position of the
 * most significant bit set in the integer.
 *
 * @param n the number to process.  If this is zero, then an exception is
 *          thrown.
 * @return the log2 of n.  If n is zero, then an exception is thrown.
 */
inline boost::uint32_t log2_uint32( boost::uint32_t n ) {
    if( n == 0 ) {
        throw std::runtime_error( "log2_uint32 Error: value is zero" );
    }

#ifdef _WIN32
    unsigned long index;
    _BitScanReverse( &index, n );
    return index;
#else
    const int leadingZeros = __builtin_clz( n );
    return 31 - leadingZeros;
#endif
}

inline bool is_power_of_two( boost::uint32_t n ) {
    if( n > 0 ) {
        return ( n & ( n - 1 ) ) == 0;
    }
    return false;
}
} // namespace math
} // namespace frantic
