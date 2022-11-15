// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/ray3f.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace geometry {

namespace mixed_kdtree_detail {
typedef boost::uint32_t index_t;

enum kdtreeEvent_type { END_ = 0, PARALLEL, BEGIN };

enum kdtreeSide_type {
    LEFT = 0,
    RIGHT,
};

const int COST_TRAVERSAL = 15;
const int COST_INTERSECT = 20;

// the tree now stores heterogeneous primitives
// their intersect cost may vary widely
// so should this hold the primitive's intersect cost as well ?
// then each primitive will need to provide its intersection
// cost
struct kdtreeEvent {
    kdtreeEvent_type type;
    float location;
    index_t index;
    int axis;

    kdtreeEvent() {}
    kdtreeEvent( index_t i, kdtreeEvent_type t, int a, float l )
        : type( t )
        , location( l )
        , index( i )
        , axis( a ) {}
    bool operator<( const kdtreeEvent& rhs ) const {
        return location < rhs.location ||
               ( location == rhs.location && ( axis < rhs.axis || ( axis == rhs.axis && type < rhs.type ) ) );
    }
};

inline float SAH_cost( float voxelSA, float leftSA, float rightSA, std::size_t nLeft, std::size_t nRight ) {
    const float lambda = ( ( nLeft == 0 || nRight == 0 ) ? 0.8f : 1.f );
    return lambda * ( COST_TRAVERSAL + COST_INTERSECT * ( leftSA / voxelSA * nLeft + rightSA / voxelSA * nRight ) );
}

inline void extract_events( std::vector<kdtreeEvent>& e, index_t primitiveIndex,
                            const frantic::graphics::boundbox3f& box ) {
    for( int axis = 0; axis < 3; ++axis ) {
        if( box.size( axis ) > 0 ) {
            e.push_back( kdtreeEvent( primitiveIndex, BEGIN, axis, box.minimum()[axis] ) );
            e.push_back( kdtreeEvent( primitiveIndex, END_, axis, box.maximum()[axis] ) );
        } else {
            e.push_back( kdtreeEvent( primitiveIndex, PARALLEL, axis, box.minimum()[axis] ) );
        }
    }
}
} // namespace mixed_kdtree_detail

namespace detail {

class common_named_channel_setters { //: public output_channel_map_listener {
    frantic::channels::channel_cvt_accessor<float> m_tAccessorFloat;
    frantic::channels::channel_cvt_accessor<double> m_tAccessorDouble;
    frantic::channels::channel_cvt_accessor<float> m_distanceAccessorFloat;
    frantic::channels::channel_cvt_accessor<double> m_distanceAccessorDouble;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_positionAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_normalAccessor;
    frantic::channels::channel_cvt_accessor<double> m_timeAccessor;

    void reset( void );

  public:
    common_named_channel_setters();
    common_named_channel_setters( frantic::channels::channel_map& cm );

    void set_channel_map( const frantic::channels::channel_map& cm );
    void set_ray( char* /*data*/, const frantic::graphics::ray3f& /*ray*/ ) const;
    void set_t( char* data, float t ) const;
    void set_t( char* data, double t ) const;
    void set_distance( char* data, float distance ) const;
    void set_distance( char* data, double distance ) const;
    void set_position( char* data, const frantic::graphics::vector3f& position ) const;
    void set_normal( char* data, const frantic::graphics::vector3f& normal ) const;
    void set_time( char* data, const double time ) const;
};

class convert_and_set_channel {
    std::size_t m_destOffset;
    std::size_t m_arity;
    frantic::channels::channel_type_convertor_function_t m_convertAndSet;

  public:
    convert_and_set_channel();

    convert_and_set_channel( frantic::channels::channel& destCh, std::size_t sourceDataArity,
                             frantic::channels::data_type_t sourceDataType,
                             const frantic::tstring& channelNameForErrorMessage );
    void reset( frantic::channels::channel& destCh, std::size_t sourceDataArity,
                frantic::channels::data_type_t sourceDataType, const frantic::tstring& channelNameForErrorMessage );
    void set( char* dest, const char* value ) const;
    std::size_t arity( void ) const;
};
} // namespace detail
} // namespace geometry
} // namespace frantic
