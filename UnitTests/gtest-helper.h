// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include "gtest/gtest.h"

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/named_channel_data.hpp>
#include <frantic/channels/property_map.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/quat4f.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

using frantic::channels::channel_general_accessor;
using frantic::channels::channel_map;
using frantic::channels::data_type_t;
using frantic::channels::property_map;
using frantic::graphics::boundbox3t;
using frantic::graphics::quat4t;
using frantic::graphics::transform4t;
using frantic::graphics::vector3t;
using frantic::particles::particle_array;
using frantic::particles::particle_file_metadata;
using frantic::particles::particle_file_stream_factory_object;
using frantic::particles::streams::particle_istream;

#define EXPECT_VECTOR3F_EQ( expected, actual )                                                                         \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpVector3tEQ<float>, expected, actual );

#define EXPECT_VECTOR3FD_EQ( expected, actual )                                                                        \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpVector3tEQ<double>, expected, actual );

#define EXPECT_VECTOR3F_NEAR( expected, actual, abs_error )                                                            \
    EXPECT_PRED_FORMAT3( ::testing::internal::vector3tNear<float>, expected, actual, abs_error )

#define EXPECT_VECTOR3FD_NEAR( expected, actual, abs_error )                                                           \
    EXPECT_PRED_FORMAT3( ::testing::internal::vector3tNear<double>, expected, actual, abs_error )

#define EXPECT_BOUNDBOX3F_EQ( expected, actual )                                                                       \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpBoundbox3tEQ<float>, expected, actual );

#define EXPECT_BOUNDBOX3FD_EQ( expected, actual )                                                                      \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpBoundbox3tEQ<double>, expected, actual );

#define EXPECT_QUAT4F_EQ( expected, actual )                                                                           \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpQuat4tEQ<float>, expected, actual );

#define EXPECT_QUAT4FD_EQ( expected, actual )                                                                          \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpQuat4tEQ<double>, expected, actual );

#define EXPECT_QUAT4F_NEAR( expected, actual, abs_error )                                                              \
    EXPECT_PRED_FORMAT3( ::testing::internal::cmpQuat4tNear<float>, expected, actual, abs_error );

#define EXPECT_QUAT4FD_NEAR( expected, actual, abs_error )                                                             \
    EXPECT_PRED_FORMAT3( ::testing::internal::cmpQuat4tNear<double>, expected, actual, abs_error );

#define EXPECT_SAME_CHANNELS( expected, actual )                                                                       \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpChannelMap, expected, actual );

#define EXPECT_TRANSFORM4F_EQ( expected, actual )                                                                      \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpTransform4tEQ<float>, expected, actual )

#define EXPECT_TRANSFORM4FD_EQ( expected, actual )                                                                     \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpTransform4tEQ<double>, expected, actual )

#define EXPECT_TRANSFORM4F_NEAR( expected, actual, abs_error )                                                         \
    EXPECT_PRED_FORMAT3( ::testing::internal::cmpTransform4tNear<float>, expected, actual, abs_error )

#define EXPECT_TRANSFORM4FD_NEAR( expected, actual, abs_error )                                                        \
    EXPECT_PRED_FORMAT3( ::testing::internal::cmpTransform4tNear<double>, expected, actual, abs_error )

#define EXPECT_PARTICLE_ARRAY_EQ( expected, actual )                                                                   \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpParticleArray, expected, actual );

#define EXPECT_PROPERTY_MAP_EQ( expected, actual )                                                                     \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpPropertyMap, expected, actual );

#define EXPECT_PARTICLE_FILE_EQ( expected, actual )                                                                    \
    EXPECT_PRED_FORMAT2( ::testing::internal::cmpParticleFile, expected, actual );

#define EXPECT_PARTICLE_FILE_EQ_POS_HINT( expected, actual, position_type_hint )                                       \
    EXPECT_PRED_FORMAT3( ::testing::internal::cmpParticleFile, expected, actual, position_type_hint );

namespace testing {
namespace internal {
template <typename RawType>
bool almostEquals( vector3t<RawType> expected, vector3t<RawType> actual ) {
    FloatingPoint<RawType> lhsX( expected.x ), lhsY( expected.y ), lhsZ( expected.z );
    FloatingPoint<RawType> rhsX( actual.x ), rhsY( actual.y ), rhsZ( actual.z );
    if( lhsX.AlmostEquals( rhsX ) && lhsY.AlmostEquals( rhsY ) && lhsZ.AlmostEquals( rhsZ ) ) {
        return true;
    }
    return false;
}

template <typename RawType>
bool almostEquals( quat4t<RawType> expected, quat4t<RawType> actual ) {
    FloatingPoint<RawType> lhsX( expected.x ), lhsY( expected.y ), lhsZ( expected.z ), lhsW( expected.w );
    FloatingPoint<RawType> rhsX( actual.x ), rhsY( actual.y ), rhsZ( actual.z ), rhsW( actual.w );
    if( lhsX.AlmostEquals( rhsX ) && lhsY.AlmostEquals( rhsY ) && lhsZ.AlmostEquals( rhsZ ) &&
        lhsW.AlmostEquals( rhsW ) ) {
        return true;
    }
    return false;
}

template <typename RawType>
bool almostEquals( transform4t<RawType> expected, transform4t<RawType> actual ) {
    for( int i = 0; i < 4 * 4; ++i ) {
        FloatingPoint<RawType> lhs( expected[i] ), rhs( actual[i] );
        if( !lhs.AlmostEquals( rhs ) )
            return false;
    }
    return true;
}

template <typename RawType>
AssertionResult cmpVector3tEQ( const char* expected_expression, const char* actual_expression,
                               vector3t<RawType> expected, vector3t<RawType> actual ) {
    const vector3t<RawType> lhs( expected ), rhs( actual );
    if( almostEquals( lhs, rhs ) ) {
        return AssertionSuccess();
    }
    ::std::stringstream expected_ss;
    expected_ss << std::setprecision( std::numeric_limits<RawType>::digits10 + 2 ) << expected;

    ::std::stringstream actual_ss;
    actual_ss << std::setprecision( std::numeric_limits<RawType>::digits10 + 2 ) << actual;
    return EqFailure( expected_expression, actual_expression, StringStreamToString( &expected_ss ),
                      StringStreamToString( &actual_ss ), false );
}

template <typename RawType>
AssertionResult vector3tNear( const char* expr1, const char* expr2, const char* abs_error_expr,
                              const vector3t<RawType>& val1, const vector3t<RawType>& val2, RawType abs_error ) {
    const RawType diff = vector3t<RawType>::distance( val1, val2 );
    if( diff <= abs_error )
        return AssertionSuccess();

    // TODO(wan): do not print the value of an expression if it's
    // already a literal.
    return AssertionFailure() << "The distance between " << expr1 << " and " << expr2 << " is " << diff
                              << ", which exceeds " << abs_error_expr << ", where\n"
                              << expr1 << " evaluates to " << val1 << ",\n"
                              << expr2 << " evaluates to " << val2 << ", and\n"
                              << abs_error_expr << " evaluates to " << abs_error << ".";
}

template <typename RawType>
AssertionResult cmpBoundbox3tEQ( const char* expected_expression, const char* actual_expression,
                                 boundbox3t<RawType> expected, boundbox3t<RawType> actual ) {
    const vector3t<RawType> lhsMin( expected.minimum() ), lhsMax( expected.maximum() );
    const vector3t<RawType> rhsMin( actual.minimum() ), rhsMax( actual.maximum() );
    if( almostEquals( lhsMin, rhsMin ) && almostEquals( lhsMax, rhsMax ) ) {
        return AssertionSuccess();
    }
    ::std::stringstream expected_ss;
    expected_ss << std::setprecision( std::numeric_limits<RawType>::digits10 + 2 ) << expected;

    ::std::stringstream actual_ss;
    actual_ss << std::setprecision( std::numeric_limits<RawType>::digits10 + 2 ) << actual;
    return EqFailure( expected_expression, actual_expression, StringStreamToString( &expected_ss ),
                      StringStreamToString( &actual_ss ), false );
}

template <typename RawType>
AssertionResult cmpQuat4tEQ( const char* expected_expression, const char* actual_expression, quat4t<RawType> expected,
                             quat4t<RawType> actual ) {
    const quat4t<RawType> lhs( expected ), rhs( actual );
    if( almostEquals( lhs, rhs ) ) {
        return AssertionSuccess();
    }
    ::std::stringstream expected_ss;
    expected_ss << std::setprecision( std::numeric_limits<RawType>::digits10 + 2 ) << expected;

    ::std::stringstream actual_ss;
    actual_ss << std::setprecision( std::numeric_limits<RawType>::digits10 + 2 ) << actual;
    return EqFailure( expected_expression, actual_expression, StringStreamToString( &expected_ss ),
                      StringStreamToString( &actual_ss ), false );
}

template <typename RawType>
AssertionResult cmpQuat4tNear( const char* expr1, const char* expr2, const char* abs_error_expr,
                               const quat4t<RawType>& val1, const quat4t<RawType>& val2, RawType abs_error ) {
    const RawType diff = ( val1 - val2 ).magnitude();
    if( diff <= abs_error )
        return AssertionSuccess();

    // TODO(wan): do not print the value of an expression if it's
    // already a literal.
    return AssertionFailure() << "The distance between " << expr1 << " and " << expr2 << " is " << diff
                              << ", which exceeds " << abs_error_expr << ", where\n"
                              << expr1 << " evaluates to " << val1 << ",\n"
                              << expr2 << " evaluates to " << val2 << ", and\n"
                              << abs_error_expr << " evaluates to " << abs_error << ".";
}

template <typename RawType>
AssertionResult cmpTransform4tEQ( const char* expected_expression, const char* actual_expression,
                                  transform4t<RawType> expected, transform4t<RawType> actual ) {
    const transform4t<RawType> lhs( expected ), rhs( actual );
    if( almostEquals( lhs, rhs ) ) {
        return AssertionSuccess();
    }
    ::std::stringstream expected_ss;
    expected_ss << std::setprecision( std::numeric_limits<RawType>::digits10 + 2 ) << expected;

    ::std::stringstream actual_ss;
    actual_ss << std::setprecision( std::numeric_limits<RawType>::digits10 + 2 ) << actual;
    return EqFailure( expected_expression, actual_expression, StringStreamToString( &expected_ss ),
                      StringStreamToString( &actual_ss ), false );
}

template <typename RawType>
AssertionResult cmpTransform4tNear( const char* expected_expression, const char* actual_expression,
                                    const char* abs_error_expr, transform4t<RawType> expected,
                                    transform4t<RawType> actual, RawType abs_error ) {
    RawType diff = 0;
    for( int i = 0; i < 16; ++i ) {
        diff += frantic::math::square( expected[i] - actual[i] );
    }
    if( diff <= abs_error * abs_error )
        return AssertionSuccess();

    // TODO: do not print the value of an expression if it's
    // already a literal.
    return AssertionFailure() << "The distance between\n"
                              << expected << "\nand\n"
                              << actual << " is " << std::sqrt( diff ) << ", which exceeds " << abs_error_expr
                              << ", where\n"
                              << expected_expression << " evaluates to " << expected << ",\n"
                              << actual_expression << " evaluates to " << actual << ", and\n"
                              << abs_error_expr << " evaluates to " << abs_error << ".";
}

// Test if the two channel maps have channels with the same names, types, and arities.
AssertionResult cmpChannelMap( const char* expected_cm_expr, const char* actual_cm_expr, const channel_map& expected_cm,
                               const channel_map& actual_cm );

// Test if the expected and actual channel data pointers have the same data.
AssertionResult cmpChannelData( const char* expected_expr, const char* actual_expr, const char* channel_name_expr,
                                const char* data_type_expr, const char* arity_expr, const char* expected,
                                const char* actual, const frantic::tstring& channel_name, data_type_t data_type,
                                size_t arity );

// Test if the particle arrays have the same particle count and particle data.
AssertionResult cmpParticleArray( const char* expected_expr, const char* actual_expr, const particle_array& expected,
                                  const particle_array& actual );

// Test if the istreams have equivalent channel map, particle count, and particle data.
AssertionResult cmpParticleData( const char* expected_expr, const char* actual_expr,
                                 boost::shared_ptr<particle_istream> expected,
                                 boost::shared_ptr<particle_istream> actual );

// Test if the property maps have the same property fields and data values.
AssertionResult cmpPropertyMap( const char* expected_expr, const char* actual_expr,
                                const frantic::channels::property_map& expected,
                                const frantic::channels::property_map& actual );

// Test if the metadata have the same general metadata and channel metadata.
AssertionResult cmpMetadata( const char* expected_expr, const char* actual_expr, const particle_file_metadata& expected,
                             const particle_file_metadata& actual );

// Test if the files produce istreams with the same particle data and metadata. Optionally specify position_type_hint.
AssertionResult cmpParticleFile( const char* expected_expr, const char* actual_expr,
                                 const char* position_type_hint_expr, const frantic::tstring& expected_file,
                                 const frantic::tstring& actual_file, data_type_t position_type_hint );

AssertionResult cmpParticleFile( const char* expected_expr, const char* actual_expr,
                                 const frantic::tstring& expected_file, const frantic::tstring& actual_file );
} // namespace internal
} // namespace testing
