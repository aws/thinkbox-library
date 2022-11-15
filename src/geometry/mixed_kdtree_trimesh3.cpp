// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/foreach.hpp>

#include <frantic/geometry/mixed_kdtree.hpp>
#include <frantic/geometry/mixed_kdtree_trimesh3.hpp>

#include <frantic/geometry/triangle_utils.hpp>
#include <frantic/math/polynomial_roots.hpp>

using namespace frantic;
using namespace frantic::geometry;
using namespace frantic::channels;

namespace frantic {
namespace geometry {

namespace detail {

//
// trimesh3_vertex_velocity_accessor
//

trimesh3_vertex_velocity_accessor::trimesh3_vertex_velocity_accessor( const frantic::geometry::trimesh3& mesh )
    : m_hasVelocity( false )
    , m_hasConvert( false ) {
    if( mesh.has_vertex_channel( _T("Velocity") ) ) {
        m_ca = mesh.get_vertex_channel_general_accessor( _T("Velocity") );
        if( m_ca.arity() != 3 ) {
            throw std::runtime_error(
                "trimesh3_vertex_velocity_accessor Error: Velocity must have arity 3, but instead it is " +
                boost::lexical_cast<std::string>( m_ca.arity() ) );
        }
        if( m_ca.has_custom_faces() ) {
            throw std::runtime_error( "trimesh3_vertex_velocity_accessor Error: The velocity channel has custom faces, "
                                      "which are not supported." );
        }
        if( m_ca.data_type() != frantic::channels::data_type_float32 ) {
            m_convert = get_channel_type_convertor_function( m_ca.data_type(), frantic::channels::data_type_float32,
                                                             _T("Velocity") );
            m_hasConvert = true;
        }
        m_hasVelocity = true;
    }
}

frantic::graphics::vector3f trimesh3_vertex_velocity_accessor::get( const std::size_t vertexNumber ) const {
    frantic::graphics::vector3f velocity( 0 );
    if( m_hasVelocity ) {
        if( m_hasConvert ) {
            m_convert( reinterpret_cast<char*>( &velocity ), m_ca.data( vertexNumber ), 3 );
        } else {
            velocity = *reinterpret_cast<const frantic::graphics::vector3f*>( m_ca.data( vertexNumber ) );
        }
    }
    return velocity;
}

//
// set_channel_from_trimesh3_vertex_channel
//

set_channel_from_trimesh3_vertex_channel::set_channel_from_trimesh3_vertex_channel(
    channels::channel& outChannel, const frantic::geometry::trimesh3& mesh, const frantic::tstring& channelName )
    : m_inputChannelAccessor( mesh.get_vertex_channel_general_accessor( channelName ) )
    , m_outputOffset( outChannel.offset() ) {
    // m_setter.reset( outChannel, m_inputChannelAccessor.arity(), m_inputChannelAccessor.data_type(), channelName );
    if( outChannel.arity() != m_inputChannelAccessor.arity() ) {
        throw std::runtime_error(
            std::string( "set_channel_from_trimesh3_vertex_channel Error: cannot convert channel " +
                         frantic::strings::to_string( channelName ) +
                         " due to differing arity (mesh input channel arity: " +
                         boost::lexical_cast<std::string>( m_inputChannelAccessor.arity() ) +
                         ", output arity: " + boost::lexical_cast<std::string>( outChannel.arity() ) + ")" )
                .c_str() );
    }
    m_setFromWeightedSum = channel_weighted_sum_combine_and_convert_function( m_inputChannelAccessor.data_type(),
                                                                              outChannel.data_type(), channelName );
}

bool set_channel_from_trimesh3_vertex_channel::is_valid( channels::channel& outChannel, const trimesh3& mesh,
                                                         const frantic::tstring& channelName ) {
    try {
        set_channel_from_trimesh3_vertex_channel test( outChannel, mesh, channelName );
    } catch( const std::exception& /*e*/ ) {
        // std::cout << e.what() << "\n";
        return false;
    }
    return true;
}

void set_channel_from_trimesh3_vertex_channel::set( char* outData, const std::size_t faceNumber,
                                                    const vector3f& barycentricCoord ) const {
    // from get_barycentric in
    // frantic::geometry::trimesh3_vertex_channel_general_accessor_base
    const vector3& face( m_inputChannelAccessor.face( faceNumber ) );
    const char* dataArray[3];
    dataArray[0] = m_inputChannelAccessor.data( face.x );
    dataArray[1] = m_inputChannelAccessor.data( face.y );
    dataArray[2] = m_inputChannelAccessor.data( face.z );
    if( m_setFromWeightedSum ) {
        m_setFromWeightedSum( &barycentricCoord[0], dataArray, 3, m_inputChannelAccessor.arity(),
                              outData + m_outputOffset );
    } else {
        FF_LOG( warning ) << "set_channel_from_trimesh3_vertex_channel::set Warning: called with NULL convertor\n";
    }
}

//
// set_channel_from_trimesh3_face_channel
//

set_channel_from_trimesh3_face_channel::set_channel_from_trimesh3_face_channel( channels::channel& outChannel,
                                                                                const trimesh3& mesh,
                                                                                const frantic::tstring& channelName )
    : m_channelAccessor( mesh.get_face_channel_general_accessor( channelName ) ) {
    m_setter.reset( outChannel, m_channelAccessor.arity(), m_channelAccessor.data_type(), channelName );
}

bool set_channel_from_trimesh3_face_channel::is_valid( channels::channel& outChannel, const trimesh3& mesh,
                                                       const frantic::tstring& channelName ) {
    try {
        set_channel_from_trimesh3_face_channel test( outChannel, mesh, channelName );
    } catch( const std::exception& /*e*/ ) {
        // std::cout << e.what() << "\n";
        return false;
    }
    return true;
}

void set_channel_from_trimesh3_face_channel::set( char* outData, const std::size_t faceNumber ) const {
    m_setter.set( outData, m_channelAccessor.data( faceNumber ) );
}

//
// set_channel_map_data_from_trimesh3
//

void set_channel_map_data_from_trimesh3::reset( void ) {
    m_setFromVertexChannel.clear();
    m_setFromFaceChannel.clear();

    m_setBarycentricCoord.reset();
    m_setGeometricNormal.reset();
    m_setFaceIndex.reset();
    m_setFaceNumber.reset();
}

set_channel_map_data_from_trimesh3::set_channel_map_data_from_trimesh3(
    boost::shared_ptr<frantic::geometry::trimesh3> mesh )
    : m_mesh( mesh ) {}

void set_channel_map_data_from_trimesh3::set_channel_map( const frantic::channels::channel_map& channelMap ) {
    reset();

    // look for matching channels in the output channel map
    for( std::size_t i = 0; i < channelMap.channel_count(); ++i ) {
        frantic::channels::channel ch = channelMap[i];

        if( m_mesh->has_vertex_channel( ch.name() ) ) {
            if( set_channel_from_trimesh3_vertex_channel::is_valid( ch, *m_mesh, ch.name() ) ) {
                m_setFromVertexChannel.push_back( set_channel_from_trimesh3_vertex_channel( ch, *m_mesh, ch.name() ) );
                continue;
            }
        }

        if( m_mesh->has_face_channel( ch.name() ) ) {
            if( set_channel_from_trimesh3_face_channel::is_valid( ch, *m_mesh, ch.name() ) ) {
                m_setFromFaceChannel.push_back( set_channel_from_trimesh3_face_channel( ch, *m_mesh, ch.name() ) );
                continue;
            }
        }

        if( ch.name() == _T("BarycentricCoord") ) {
            m_setBarycentricCoord = channelMap.get_cvt_accessor<frantic::graphics::vector3f>( ch.name() );
        } else if( ch.name() == _T("FaceNumber") ) {
            m_setFaceNumber = channelMap.get_cvt_accessor<boost::int32_t>( ch.name() );
        } else if( ch.name() == _T("FaceIndex") ) {
            m_setFaceIndex = channelMap.get_cvt_accessor<boost::int32_t>( ch.name() );
        } else if( ch.name() == _T("GeometricNormal") ) {
            m_setGeometricNormal = channelMap.get_cvt_accessor<frantic::graphics::vector3f>( ch.name() );
        }
    }
}

void set_channel_map_data_from_trimesh3::set_channel_map_data(
    char* outData, const boost::int32_t faceNumber, const frantic::graphics::vector3f& barycentricCoord,
    const frantic::graphics::vector3f& geometricNormal ) const {
    BOOST_FOREACH( const set_channel_from_trimesh3_vertex_channel& setter, m_setFromVertexChannel ) {
        setter.set( outData, faceNumber, barycentricCoord );
    }

    BOOST_FOREACH( const set_channel_from_trimesh3_face_channel& setter, m_setFromFaceChannel ) {
        setter.set( outData, faceNumber );
    }

    if( m_setFaceIndex.is_valid() )
        m_setFaceIndex.set( outData, faceNumber );
    if( m_setFaceNumber.is_valid() )
        m_setFaceNumber.set( outData, faceNumber );
    if( m_setBarycentricCoord.is_valid() )
        m_setBarycentricCoord.set( outData, barycentricCoord );
    if( m_setGeometricNormal.is_valid() )
        m_setGeometricNormal.set( outData, geometricNormal );
}

} // namespace detail

//
// mixed_kdtree_trimesh3
//
/*
void mixed_kdtree_trimesh3::initialize_from_mesh( boost::shared_ptr<frantic::geometry::trimesh3> mesh ) {
    m_primitives.reserve( mesh->face_count() );

    for( std::size_t faceNumber = 0; faceNumber < mesh->face_count(); ++faceNumber ) {
    boost::shared_ptr<mixed_kdtree_trimesh3_primitive> primitive( new mixed_kdtree_trimesh3_primitive( * this,
faceNumber ) );
        //m_primitives.push_back( mixed_kdtree_trimesh3_primitive( * this, faceNumber ) );
        m_primitives.push_back( primitive );
    }
}
*/
/*
mixed_kdtree_trimesh3::mixed_kdtree_trimesh3( boost::shared_ptr<frantic::geometry::trimesh3> mesh, bool initFromMesh )
  :	m_mesh( mesh ),
    m_setChannelMapDataFromMesh( mesh )
{
  if( initFromMesh ) {
    initialize_from_mesh( mesh );
  }
}
*/

mixed_kdtree_trimesh3::mixed_kdtree_trimesh3( boost::shared_ptr<frantic::geometry::trimesh3> mesh )
    : m_mesh( mesh )
    , m_setChannelMapDataFromMesh( mesh ) {
    if( !mesh ) {
        throw std::runtime_error(
            "mixed_kdtree_trimesh3::mixed_kdtree_trimesh3 Error: attempting to contruct from NULL mesh" );
    }
    // initialize_from_mesh( mesh );
}

mixed_kdtree_trimesh3::~mixed_kdtree_trimesh3( void ) {
    // std::cout << " destroying helper\n";
}
/*
boost::shared_ptr<mixed_kdtree_trimesh3> mixed_kdtree_trimesh3::from_face_subset(
boost::shared_ptr<frantic::geometry::trimesh3> mesh, const std::vector<std::size_t> & faceNumbers ) {
  boost::shared_ptr<mixed_kdtree_trimesh3> kdtreeMesh( new mixed_kdtree_trimesh3( mesh ) );

  kdtreeMesh->m_primitives.reserve( faceNumbers.size() );

  BOOST_FOREACH( const std::size_t faceNumber, faceNumbers ) {
    if( faceNumber > mesh->face_count() ) {
      throw std::runtime_error( "mixed_kdtree_trimesh3::from_face_subset Error: face number is out of range of faces in
the mesh" );
    }
    boost::shared_ptr<mixed_kdtree_trimesh3_primitive> primitive( new mixed_kdtree_trimesh3_primitive( * kdtreeMesh,
faceNumber ) ); kdtreeMesh->m_primitives.push_back( primitive );
  }

  return kdtreeMesh;
}
*/
std::size_t mixed_kdtree_trimesh3::get_primitive_count( void ) { return m_mesh->face_count(); }

void mixed_kdtree_trimesh3::get_primitives( mixed_kdtree_primitive_recorder* recorder ) {
    if( !recorder ) {
        throw std::runtime_error( "mixed_kdtree_trimesh3::get_primitives Error: primitive recorder is NULL" );
    }

    if( m_mesh->face_count() > std::numeric_limits<boost::uint32_t>::max() ) {
        throw std::runtime_error(
            "mixed_kdtree_trimesh3::get_primitives Internal Error: Internal Error: the number of faces (" +
            boost::lexical_cast<std::string>( m_mesh->face_count() ) + ") exceeds the maximum face index (" +
            boost::lexical_cast<std::string>( std::numeric_limits<boost::uint32_t>::max() ) + ")." );
    }

    for( std::size_t faceNumber = 0; faceNumber < m_mesh->face_count(); ++faceNumber ) {
        boost::shared_ptr<mixed_kdtree_trimesh3_primitive> primitive(
            new mixed_kdtree_trimesh3_primitive( *this, static_cast<boost::uint32_t>( faceNumber ) ) );
        recorder->insert( primitive );
    }
    /*
      BOOST_FOREACH( boost::shared_ptr<mixed_kdtree_trimesh3_primitive> & primitive, m_primitives ) {
            recorder->insert( primitive );
        }
      */
    recorder->insert_output_channel_map_listener( this );
}

trimesh3& mixed_kdtree_trimesh3::get_mesh_ref( void ) { return *m_mesh; }

void mixed_kdtree_trimesh3::to_channel_map_data( char* outData, const mixed_kdtree_trimesh3_primitive* primitive,
                                                 const mixed_kdtree_point_data* pointData ) const {
    const char* primitiveData( pointData->get_primitive_data() );
    const boost::int32_t faceNumber( primitive->get_face_number() );
    const vector3f barycentricCoord( primitive->get_barycentric_coord( primitiveData ) );
    const vector3f geometricNormal( primitive->get_geometric_normal( primitiveData ) );

    m_setChannelMapDataFromMesh.set_channel_map_data( outData, faceNumber, barycentricCoord, geometricNormal );
}

void mixed_kdtree_trimesh3::set_output_channel_map( const channels::channel_map& channelMap ) {
    m_setChannelMapDataFromMesh.set_channel_map( channelMap );
}

//
// trimesh3_constant_velocity_kdtree_primitive_helper
//

mixed_kdtree_trimesh3_constant_velocity::mixed_kdtree_trimesh3_constant_velocity( boost::shared_ptr<trimesh3> mesh,
                                                                                  const double t0, const double t1 )
    : m_mesh( mesh )
    , m_t0( t0 )
    , m_t1( t1 )
    , m_setChannelMapDataFromMesh( mesh ) {
    if( !mesh ) {
        throw std::runtime_error(
            "mixed_kdtree_trimesh3_constant_velocity::mixed_kdtree_trimesh3_constant_velocity Error: "
            "attempting to contruct from NULL mesh" );
    }

    init_vertex_position_and_displacement( *mesh, t0, t1 );

    // m_primitives.reserve( mesh->face_count() );
    m_staticGeometry = boost::shared_ptr<mixed_kdtree_trimesh3>( new mixed_kdtree_trimesh3( mesh ) );
}

const trimesh3& mixed_kdtree_trimesh3_constant_velocity::get_mesh_ref( void ) { return *m_mesh; }

std::size_t mixed_kdtree_trimesh3_constant_velocity::get_primitive_count( void ) { return m_mesh->face_count(); }

void mixed_kdtree_trimesh3_constant_velocity::get_primitives( mixed_kdtree_primitive_recorder* recorder ) {
    if( !recorder ) {
        throw std::runtime_error(
            "mixed_kdtree_trimesh3_constant_velocity::get_primitives Error: primitive recorder is NULL" );
    }

    detail::trimesh3_vertex_velocity_accessor velocityAccessor( *m_mesh );

    for( std::size_t faceNumber = 0; faceNumber < m_mesh->face_count(); ++faceNumber ) {
        // only use an animated primitive if the face is actually moving
        const vector3& face( m_mesh->get_face( faceNumber ) );
        if( velocityAccessor.has_velocity() && ( velocityAccessor.get( face.x ).get_magnitude_squared() > 0 ||
                                                 velocityAccessor.get( face.y ).get_magnitude_squared() > 0 ||
                                                 velocityAccessor.get( face.z ).get_magnitude_squared() > 0 ) ) {
            boost::shared_ptr<mixed_kdtree_trimesh3_constant_velocity_primitive> p(
                new mixed_kdtree_trimesh3_constant_velocity_primitive( *this, static_cast<boost::int32_t>( faceNumber ),
                                                                       m_t0, m_t1 ) );
            recorder->insert( p );
        } else {
            boost::shared_ptr<mixed_kdtree_trimesh3_primitive> p(
                new mixed_kdtree_trimesh3_primitive( *m_staticGeometry, static_cast<boost::int32_t>( faceNumber ) ) );
            recorder->insert( p );
        }
    }

    recorder->insert_output_channel_map_listener( this );
    recorder->insert_output_channel_map_listener( m_staticGeometry.get() );
}

vector3f
mixed_kdtree_trimesh3_constant_velocity::get_vertex_initial_position( const boost::int32_t vertexNumber ) const {
    return m_vertexInitialPosition[vertexNumber];
}

vector3f mixed_kdtree_trimesh3_constant_velocity::get_vertex_displacement( const boost::int32_t vertexNumber ) const {
    return m_vertexDisplacement[vertexNumber];
}

void mixed_kdtree_trimesh3_constant_velocity::init_vertex_position_and_displacement( const trimesh3& mesh,
                                                                                     const double t0,
                                                                                     const double t1 ) {
    const std::size_t vertexCount( mesh.vertex_count() );
    bool success = false;

    m_vertexInitialPosition.clear();
    m_vertexInitialPosition.reserve( vertexCount );

    m_vertexDisplacement.clear();
    m_vertexDisplacement.reserve( vertexCount );

    if( mesh.has_vertex_channel( _T("Velocity") ) ) {
        // how should we handle arity() != 3 and non-floating point types ?
        // right now this will throw a runtime error
        const_trimesh3_vertex_channel_cvt_accessor<vector3f> getVelocity(
            mesh.get_vertex_channel_cvt_accessor<vector3f>( _T("Velocity") ) );
        const float timeStep = float( t1 - t0 );
        vector3f velocity;

        for( std::size_t vertexNumber = 0; vertexNumber < vertexCount; ++vertexNumber ) {
            velocity = getVelocity[vertexNumber];

            m_vertexInitialPosition.push_back( mesh.get_vertex( vertexNumber ) + velocity * float( t0 ) );
            m_vertexDisplacement.push_back( velocity * timeStep );
        }

        success = true;
    }

    if( !success ) {
        m_vertexInitialPosition.clear();
        std::copy( mesh.vertices_ref().begin(), mesh.vertices_ref().end(), back_inserter( m_vertexInitialPosition ) );
        m_vertexDisplacement.clear();
        m_vertexDisplacement.resize( vertexCount, vector3f( 0 ) );
    }
}

void mixed_kdtree_trimesh3_constant_velocity::reset_channel_accessors( void ) {
    m_setDisplacementDuringTimeStep.reset();
}

// called by its primitives
// this behaviour is nearly shared with the trimesh3 manager
void mixed_kdtree_trimesh3_constant_velocity::to_channel_map_data(
    char* data, const mixed_kdtree_trimesh3_constant_velocity_primitive* primitive,
    const mixed_kdtree_point_data* pointData ) const {
    const char* primitiveData( pointData->get_primitive_data() );
    const boost::int32_t faceNumber( primitive->get_face_number() );
    const vector3f barycentricCoord( primitive->get_barycentric_coord( primitiveData ) );
    const vector3f geometricNormal( primitive->get_geometric_normal( pointData ) );

    m_setChannelMapDataFromMesh.set_channel_map_data( data, faceNumber, barycentricCoord, geometricNormal );

    if( m_setDisplacementDuringTimeStep.is_valid() )
        m_setDisplacementDuringTimeStep.set( data, primitive->get_motion_during_time_step( primitiveData ) );
}

void mixed_kdtree_trimesh3_constant_velocity::set_output_channel_map( const channels::channel_map& channelMap ) {
    reset_channel_accessors();

    m_setChannelMapDataFromMesh.set_channel_map( channelMap );

    // look for matching channels in the output channel map
    for( std::size_t i = 0; i < channelMap.channel_count(); ++i ) {
        frantic::channels::channel ch = channelMap[i];

        if( ch.name() == _T("DisplacementDuringTimeStep") ) {
            // m_setDisplacementDuringTimeStep = detail::set_channel_from_mem_function1<vector3f,
            // mixed_kdtree_trimesh3_constant_velocity_primitive, const char *>( channelMap, ch.name(),
            // &mixed_kdtree_trimesh3_constant_velocity_primitive::get_motion_during_time_step );
            m_setDisplacementDuringTimeStep = channelMap.get_cvt_accessor<vector3f>( ch.name() );
        }
    }
}

//
// mixed_kdtree_trimesh3_constant_velocity_primitive
//

bool mixed_kdtree_trimesh3_constant_velocity_primitive::is_coplanar_point_in_triangle(
    const frantic::graphics::vector3f& pt, const frantic::graphics::vector3f& va, const frantic::graphics::vector3f& vb,
    const frantic::graphics::vector3f& vc, frantic::graphics::vector3f& outBarycentricCoord ) {
    outBarycentricCoord = graphics::compute_barycentric_coordinates( pt, va, vb, vc );

    if( outBarycentricCoord.x >= 0 && outBarycentricCoord.y >= 0 && outBarycentricCoord.z >= 0 ) {
        return true;
    }
    return false;
}

frantic::graphics::vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_x1( void ) const {
    return m_helper.get_vertex_initial_position( m_face.x );
};
frantic::graphics::vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_x2( void ) const {
    return m_helper.get_vertex_initial_position( m_face.y );
};
frantic::graphics::vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_x3( void ) const {
    return m_helper.get_vertex_initial_position( m_face.z );
};
void mixed_kdtree_trimesh3_constant_velocity_primitive::get_x( frantic::graphics::vector3f& x1,
                                                               frantic::graphics::vector3f& x2,
                                                               frantic::graphics::vector3f& x3 ) {
    x1 = get_x1();
    x2 = get_x2();
    x3 = get_x3();
}

frantic::graphics::vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_v1( void ) const {
    return m_helper.get_vertex_displacement( m_face.x );
};
frantic::graphics::vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_v2( void ) const {
    return m_helper.get_vertex_displacement( m_face.y );
};
frantic::graphics::vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_v3( void ) const {
    return m_helper.get_vertex_displacement( m_face.z );
};
void mixed_kdtree_trimesh3_constant_velocity_primitive::get_v( frantic::graphics::vector3f& v1,
                                                               frantic::graphics::vector3f& v2,
                                                               frantic::graphics::vector3f& v3 ) {
    v1 = get_v1();
    v2 = get_v2();
    v3 = get_v3();
}

// todo: store t0 and t1, and have the object disappear outside of these times?
// or clamp the motion to within these times ?
mixed_kdtree_trimesh3_constant_velocity_primitive::mixed_kdtree_trimesh3_constant_velocity_primitive(
    mixed_kdtree_trimesh3_constant_velocity& helper, const boost::int32_t faceNumber, const double /*t0*/,
    const double /*t1*/ )
    : m_helper( helper )
    , m_face( helper.get_mesh_ref().get_face( faceNumber ) )
    , m_faceNumber( faceNumber ) {}

/**
 * assume the ray travels from origin at time 0
 * to origin + direction at time 1
 */
bool mixed_kdtree_trimesh3_constant_velocity_primitive::intersect_particle( const frantic::graphics::ray3f& ray,
                                                                            double tMin, double tMax,
                                                                            mixed_kdtree_ray_observer* observer ) {
    const vector3f px = ray.origin();
    const vector3f pv = ray.direction();

    double roots[3];
    bool rootFound = false;
    int rootCount = 0;

    const vector3f x1( get_x1() );
    const vector3f x2( get_x2() );
    const vector3f x3( get_x3() );

    const vector3f v1( get_v1() );
    const vector3f v2( get_v2() );
    const vector3f v3( get_v3() );

    // FF_LOG( debug ) << "*\tat face " << m_faceNumber << "\n";
    // FF_LOG( debug ) << "\tx" << x1 << x2 << x3 << "\n";
    // FF_LOG( debug ) << "\tv" << v1 << v2 << v3 << "\n";

    // TODO: check if the ray's travel is parallel to the trangle

    // Case where all triangle vertices have the same velocity.
    // This is handled fine by the general procedure, however this is much
    // simpler (1st order linear vs cubic) and I assume it happens often
    // enough to warrant a special case for it.
    if( v1 == v2 && v2 == v3 ) {
        // all of this is just copied from the corresponding parts of
        // the general case below.
        const vector3f x41 = px - x1;
        const vector3f x31 = x3 - x1;
        const vector3f x21 = x2 - x1;

        const vector3f v41 = pv - v1;

        // TODO: will the face normal be cached anywhere ?
        // if the timestep's start is always 0, or we store the
        // start time, then this can just be the original face's
        // normal.
        const vector3f x21_x_x31( vector3f::cross( x21, x31 ) );

        const double a1 = vector3f::dot_double( v41, x21_x_x31 );
        const double a0 = vector3f::dot_double( x41, x21_x_x31 );

        if( fabs( a1 ) > 0 ) {
            roots[0] = -a0 / a1;
            rootCount = 1;
        }
    } else {

        const vector3f x41 = px - x1;
        const vector3f x31 = x3 - x1;
        const vector3f x21 = x2 - x1;

        const vector3f v41 = pv - v1;
        const vector3f v31 = v3 - v1;
        const vector3f v21 = v2 - v1;

        // Two of these cross products are constant.
        // It may be worthwhile to pre-calculate them.
        const vector3f v31_x_v41( vector3f::cross( v31, v41 ) );
        const vector3f x31_x_v21( vector3f::cross( x31, v21 ) ); // constant for all rays
        const vector3f x41_x_v31( vector3f::cross( x41, v31 ) );
        const vector3f x21_x_x31( vector3f::cross( x21, x31 ) ); // constant for all rays

        const double a3 = vector3f::dot_double( v21, v31_x_v41 );
        const double a2 = vector3f::dot_double( x21, v31_x_v41 ) - vector3f::dot_double( v41, x31_x_v21 ) -
                          vector3f::dot_double( v21, x41_x_v31 );
        const double a1 = vector3f::dot_double( v41, x21_x_x31 ) - vector3f::dot_double( x21, x41_x_v31 ) -
                          vector3f::dot_double( x41, x31_x_v21 );
        const double a0 = vector3f::dot_double( x41, x21_x_x31 );

        // FF_LOG( debug ) << "*\tat face " << m_faceNumber << "\n";
        // FF_LOG( debug ) << "\txdiff " << x41 << x31 << x21 << "\n";
        // FF_LOG( debug ) << "\tvdiff " << v41 << v31 << v21 << "\n";
        // FF_LOG( debug ) << "\tcrosses " << v31_x_v41 << x31_x_v21 << x41_x_v31 << x21_x_x31 << "\n";
        // FF_LOG( debug ) << "\ta " << a3 <<" " <<  a2 << " " << a1 << " " << a0 << "\n";
        //  TODO: this uses an epsilon because the get_cubic_roots doesn't
        //  work well for large coefficients.
        //  It may be worthwhile to re-write get_cubic_roots and use a solver
        // based on numerical methods, but the epsilon here covers the
        // problematic case I encountered.
        // TODO: the coefficients input to get_polynomial_roots should at
        // least be conditioned

        // TODO: should also handle non-trivial zeros
        if( fabs( a3 ) > 1.e-5 ) {
            // x^3 + ax^2 + bx + c = 0
            // std::cout << "x^3\n";
            rootCount = frantic::math::polynomial_roots::get_cubic_roots( a2 / a3, a1 / a3, a0 / a3, roots[0], roots[1],
                                                                          roots[2] );
        } else if( fabs( a2 ) > 0 ) {
            // ax^2 + bx + c = 0
            // std::cout << "x^2\n";
            rootCount = frantic::math::polynomial_roots::get_quadratic_roots<double>( a2, a1, a0, roots[0], roots[1] );
        } else if( fabs( a1 ) > 0 ) {
            // std::cout << "x^1\n";
            roots[0] = -a0 / a1;
            rootCount = 1;
        } else {
            // std::cout << "x^0\n";
        }
    }

    vector3f normal;
    vector3f barycentricCoords;
    boundbox3f intersectionBBox;
    double root = std::numeric_limits<double>::max();
    for( int i = 0; i < rootCount; ++i ) {
        // FF_LOG( debug ) << "*\tat face " << m_faceNumber << "\n";
        // FF_LOG( debug ) << "\troot: " << roots[i] << " (" << tMin << "," << tMax << ")\n";

        // std::cout << "root: " << roots[i] << " (" << tMin << "," << tMax << ")\n";// << get_bounds() << "\n";
        if( roots[i] >= tMin && roots[i] <= tMax && roots[i] < root ) {
            // need to make sure the collision point is within the triangle
            const vector3f pt = px + pv * float( roots[i] );
            const vector3f va = x1 + v1 * float( roots[i] );
            const vector3f vb = x2 + v2 * float( roots[i] );
            const vector3f vc = x3 + v3 * float( roots[i] );

            // nb: the first call sets barycentricCoords even if no
            // intersection is found.
            // If it is swapped for a function that does not, then
            // you must calculate the barycentric coords for the
            // intersection elsewhere.
            bool inside = false;
            const float eta = -0.001f;
            if( is_coplanar_point_in_triangle( pt, va, vb, vc, barycentricCoords ) ) {
                inside = true;
            } else if( barycentricCoords.x >= eta && barycentricCoords.y >= eta && barycentricCoords.z >= eta ) {
                if( intersect_triangle( ray, tMin, tMax, va, vb, vc, intersectionBBox ) ) {
                    inside = true;
                }
            }
            if( inside ) {
                rootFound = true;
                root = roots[i];
                normal = vector3f::cross( vb - va, vc - va ).to_normalized();
            }
        }
    }

    // if( rootFound && observer->wants_hit_details( root ) ) {
    if( rootFound ) {
        mixed_kdtree_trimesh3_point_data intersection( &m_helper.get_mesh_ref(), m_faceNumber, barycentricCoords );
        observer->insert( ray, root, root, normal, this, reinterpret_cast<char*>( &intersection ),
                          sizeof( intersection ) );
    }

    return rootFound;
}
// test for intersection with an instantaneous ray
// the geometry is offset to the specified time
bool mixed_kdtree_trimesh3_constant_velocity_primitive::intersect_ray( const frantic::graphics::ray3f& ray, double tMin,
                                                                       double tMax, mixed_kdtree_ray_observer* observer,
                                                                       double time ) {
    bool found = false;

    const vector3f x1( get_x1() );
    const vector3f x2( get_x2() );
    const vector3f x3( get_x3() );

    const vector3f v1( get_v1() );
    const vector3f v2( get_v2() );
    const vector3f v3( get_v3() );

    const vector3f va = x1 + v1 * static_cast<float>( time );
    const vector3f vb = x2 + v2 * static_cast<float>( time );
    const vector3f vc = x3 + v3 * static_cast<float>( time );

    const plane3f facePlane( plane3f::from_triangle( va, vb, vc ) );

    const double distance = facePlane.get_distance_to_intersection( ray );

    if( distance >= tMin && distance <= tMax ) {
        const vector3f pt( ray.at( distance ) );
        vector3f barycentricCoord;
        boundbox3f intersectionBBox;

        // NB: first function call must always set the barycentric coords,
        // even if no intersection is found.
        // If you use a function that does not do so, then you must determine
        // them elsewhere in case intersect_triangle finds an intersection
        // which is_coplanar_point_in_triangle does not.
        bool inside = false;
        const float eta = -0.001f;
        if( is_coplanar_point_in_triangle( pt, va, vb, vc, barycentricCoord ) ) {
            inside = true;
        } else if( barycentricCoord.x >= eta && barycentricCoord.y >= eta && barycentricCoord.z >= eta ) {
            if( intersect_triangle( ray, tMin, tMax, va, vb, vc, intersectionBBox ) ) {
                inside = true;
            }
        }
        if( inside ) {
            // if( observer->wants_hit_details( distance ) ) {
            mixed_kdtree_trimesh3_point_data intersection( &m_helper.get_mesh_ref(), m_faceNumber, barycentricCoord );

            observer->insert( ray, distance, time, facePlane.normal(), this, reinterpret_cast<char*>( &intersection ),
                              sizeof( mixed_kdtree_trimesh3_point_data ) );

            found = true;
            //}
        }
    }

    return found;
}

void mixed_kdtree_trimesh3_constant_velocity_primitive::get_any_point( vector3f& outPoint, vector3f& outNormal,
                                                                       double time ) {
    const vector3f x1( get_x1() );
    const vector3f x2( get_x2() );
    const vector3f x3( get_x3() );

    const vector3f v1( get_v1() );
    const vector3f v2( get_v2() );
    const vector3f v3( get_v3() );

    const vector3f va = x1 + v1 * static_cast<float>( time );
    const vector3f vb = x2 + v2 * static_cast<float>( time );
    const vector3f vc = x3 + v3 * static_cast<float>( time );

    outPoint = 1.0f / 3 * ( va + vb + vc );
    outNormal = vector3f::cross( vb - va, vc - va ).to_normalized();
}

// todo: should the output be restricted to nodebounds ?
bool mixed_kdtree_trimesh3_constant_velocity_primitive::find_nearest_point( const vector3f& point,
                                                                            boundbox3f& /*nodeBounds*/,
                                                                            mixed_kdtree_distance_observer* observer,
                                                                            double time ) {
    bool found = false;

    const vector3f x1( get_x1() );
    const vector3f x2( get_x2() );
    const vector3f x3( get_x3() );

    const vector3f v1( get_v1() );
    const vector3f v2( get_v2() );
    const vector3f v3( get_v3() );

    const vector3f va = x1 + v1 * static_cast<float>( time );
    const vector3f vb = x2 + v2 * static_cast<float>( time );
    const vector3f vc = x3 + v3 * static_cast<float>( time );

    const plane3f facePlane( plane3f::from_triangle( va, vb, vc ) );

    // Project the point onto the triangle's plane
    double distance = facePlane.get_signed_distance_to_plane_double( point );

    if( observer->is_in_range( fabs( distance ) ) ) {

        // Check that the intersection is inside the bounding box (we're intersecting the ray with the intersection of
        // the triangle and the bounding box)
        vector3f projected = point - static_cast<float>( distance ) * facePlane.normal();

        // Now get the barycentric coordinates of this projection
        // vector3f barycentricCoord = m_mesh.compute_barycentric_coordinates( m_faceNumber, projected );
        vector3f barycentricCoord = graphics::compute_barycentric_coordinates( projected, va, vb, vc );

        if( barycentricCoord.is_inf() )
            return false;

        if( barycentricCoord.x >= 0 && barycentricCoord.y >= 0 && barycentricCoord.z >= 0 ) {
            distance = fabs( distance ); // Projected and barycentricCoord are both accurate
        } else {
            detail::nearest_point_on_triangle( projected, barycentricCoord, va, vb, vc );
            distance = vector3f::distance_double( point, projected );
        }

        // if( observer->wants_hit_details( distance ) ) {
        mixed_kdtree_trimesh3_point_data pointData( &m_helper.get_mesh_ref(), m_faceNumber, barycentricCoord );

        observer->insert( projected, distance, time, facePlane.normal(), this, reinterpret_cast<char*>( &pointData ),
                          sizeof( mixed_kdtree_trimesh3_point_data ) );

        found = true;
        //}
    }

    return found;
}

boundbox3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_bounds() const {
    const vector3f x1( get_x1() );
    const vector3f x2( get_x2() );
    const vector3f x3( get_x3() );

    const vector3f v1( get_v1() );
    const vector3f v2( get_v2() );
    const vector3f v3( get_v3() );

    boundbox3f result( x1 );
    result += x2;
    result += x3;

    result += x1 + v1;
    result += x2 + v2;
    result += x3 + v3;

    return result;
}

boundbox3f mixed_kdtree_trimesh3_constant_velocity_primitive::intersect_with( const boundbox3f& box ) const {
    // TODO: make them less bad
    // note that the shape formed by the sweep over time can be self-
    // intersecting.
    boundbox3f result( get_bounds() );
    result.intersect_with( box );
    return result;
}

std::size_t mixed_kdtree_trimesh3_constant_velocity_primitive::get_data_size( void ) const {
    return sizeof( mixed_kdtree_trimesh3_point_data );
}

boost::int32_t mixed_kdtree_trimesh3_constant_velocity_primitive::get_face_number( void ) const { return m_faceNumber; }

vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_geometric_normal(
    const mixed_kdtree_point_data* pointData ) const {
    // need time which data doesn't currently have
    const double time = static_cast<float>( pointData->get_time() );

    const vector3f x1( get_x1() );
    const vector3f x2( get_x2() );
    const vector3f x3( get_x3() );

    const vector3f v1( get_v1() );
    const vector3f v2( get_v2() );
    const vector3f v3( get_v3() );

    const vector3f va = x1 + v1 * static_cast<float>( time );
    const vector3f vb = x2 + v2 * static_cast<float>( time );
    const vector3f vc = x3 + v3 * static_cast<float>( time );

    return vector3f::cross( vb - va, vc - va ).to_normalized();
}

vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_barycentric_coord( const char* data ) const {
    const mixed_kdtree_trimesh3_point_data* pointData(
        reinterpret_cast<const mixed_kdtree_trimesh3_point_data*>( data ) );
    return pointData->barycentricCoords;
}

vector3f mixed_kdtree_trimesh3_constant_velocity_primitive::get_motion_during_time_step( const char* data ) const {
    const mixed_kdtree_trimesh3_point_data* pointData(
        reinterpret_cast<const mixed_kdtree_trimesh3_point_data*>( data ) );
    return ( pointData->barycentricCoords.x * get_v1() ) + ( pointData->barycentricCoords.y * get_v2() ) +
           ( pointData->barycentricCoords.z * get_v3() );
}

void mixed_kdtree_trimesh3_constant_velocity_primitive::to_raytrace_intersection(
    raytrace_intersection& raytraceIntersection, const char* data ) const {
    const mixed_kdtree_trimesh3_point_data* pointData(
        reinterpret_cast<const mixed_kdtree_trimesh3_point_data*>( data ) );
    raytraceIntersection.faceIndex = static_cast<int>( pointData->faceNumber );
    raytraceIntersection.barycentricCoords = pointData->barycentricCoords;
    raytraceIntersection.motionDuringTimeStep = get_motion_during_time_step( data );
}

void mixed_kdtree_trimesh3_constant_velocity_primitive::to_nearest_point_search_result(
    nearest_point_search_result& searchResult, const char* data ) const {
    const mixed_kdtree_trimesh3_point_data* pointData(
        reinterpret_cast<const mixed_kdtree_trimesh3_point_data*>( data ) );
    searchResult.faceIndex = static_cast<int>( pointData->faceNumber ); // should be m_faceNumber too
    searchResult.barycentricCoords = pointData->barycentricCoords;
}

void mixed_kdtree_trimesh3_constant_velocity_primitive::to_channel_map_data(
    char* outData, const mixed_kdtree_point_data* pointData ) const {
    m_helper.to_channel_map_data( outData, this, pointData );
}

//
// mixed_kdtree_trimesh3_primitive
//

vector3f mixed_kdtree_trimesh3_primitive::compute_barycentric_coordinates( const vector3f& pt ) {
    return compute_barycentric_coordinates_with_helpers(
        pt, m_helper.get_mesh_ref().get_vertex( m_face.x ), m_helper.get_mesh_ref().get_vertex( m_face.y ),
        m_helper.get_mesh_ref().get_vertex( m_face.z ), m_barycentric0Axis, m_barycentric1Axis,
        m_barycentricInverseDeterminant );
}

mixed_kdtree_trimesh3_primitive::mixed_kdtree_trimesh3_primitive( mixed_kdtree_trimesh3& helper,
                                                                  const boost::uint32_t faceNumber )
    : m_helper( helper )
    , m_faceNumber( faceNumber )
    , m_facePlane( helper.get_mesh_ref().calculate_face_plane( faceNumber ) ) {
    trimesh3& mesh = helper.get_mesh_ref();

    const vector3& face( mesh.get_face( faceNumber ) );
    m_face = face;
    // m_vertices[0] = mesh.get_vertex(face.x);
    // m_vertices[1] = mesh.get_vertex(face.y);
    // m_vertices[2] = mesh.get_vertex(face.z);
    compute_barycentric_helpers( mesh.get_vertex( face.x ), mesh.get_vertex( face.y ), mesh.get_vertex( face.z ),
                                 m_barycentric0Axis, m_barycentric1Axis, m_barycentricInverseDeterminant );
    // get_plucker( m_vertices[0], m_vertices[1], m_pluckerA );
    // get_plucker( m_vertices[1], m_vertices[2], m_pluckerB );
    // get_plucker( m_vertices[2], m_vertices[0], m_pluckerC );
}

bool is_intersecting( const boundbox3f& box, const frantic::graphics::ray3f& ray, const double start,
                      const double end ) {
    const vector3f startPoint( ray.at( start ) );
    const vector3f endPoint( ray.at( end ) );
    if( box.is_empty() ) {
        return false;
    } else if( box.contains( startPoint ) || box.contains( endPoint ) ) {
        return true;
    } else if( ( startPoint.x < box.minimum().x && endPoint.x < box.minimum().x ) ||
               ( startPoint.y < box.minimum().y && endPoint.y < box.minimum().y ) ||
               ( startPoint.z < box.minimum().z && endPoint.z < box.minimum().z ) ||
               ( startPoint.x > box.maximum().x && endPoint.x > box.maximum().x ) ||
               ( startPoint.y > box.maximum().y && endPoint.y > box.maximum().y ) ||
               ( startPoint.z > box.maximum().z && endPoint.z > box.maximum().z ) ) {
        // A number of special cases of guaranteed non-intersection, which were very numerically
        // unstable in the test below.
        return false;
    } else {
        double minimumOfFarDistances = ( std::numeric_limits<double>::max )();
        double maximumOfNearDistances = ( std::numeric_limits<double>::min )();

        if( ray.direction().x != 0 ) {
            double invDirX = 1.0 / ray.direction().x;
            if( invDirX >= 0 ) {
                maximumOfNearDistances = ( box.minimum().x - ray.origin().x ) * invDirX;
                minimumOfFarDistances = ( box.maximum().x - ray.origin().x ) * invDirX;
            } else {
                maximumOfNearDistances = ( box.maximum().x - ray.origin().x ) * invDirX;
                minimumOfFarDistances = ( box.minimum().x - ray.origin().x ) * invDirX;
            }
        }

        if( ray.direction().y != 0 ) {
            double invDirY = 1.0 / ray.direction().y;
            if( invDirY >= 0 ) {
                maximumOfNearDistances =
                    ( std::max )( maximumOfNearDistances, ( box.minimum().y - ray.origin().y ) * invDirY );
                minimumOfFarDistances =
                    ( std::min )( minimumOfFarDistances, ( box.maximum().y - ray.origin().y ) * invDirY );
            } else {
                maximumOfNearDistances =
                    ( std::max )( maximumOfNearDistances, ( box.maximum().y - ray.origin().y ) * invDirY );
                minimumOfFarDistances =
                    ( std::min )( minimumOfFarDistances, ( box.minimum().y - ray.origin().y ) * invDirY );
            }
        }

        if( ray.direction().z != 0 ) {
            double invDirZ = 1.0 / ray.direction().z;
            if( invDirZ >= 0 ) {
                maximumOfNearDistances =
                    ( std::max )( maximumOfNearDistances, ( box.minimum().z - ray.origin().z ) * invDirZ );
                minimumOfFarDistances =
                    ( std::min )( minimumOfFarDistances, ( box.maximum().z - ray.origin().z ) * invDirZ );
            } else {
                maximumOfNearDistances =
                    ( std::max )( maximumOfNearDistances, ( box.maximum().z - ray.origin().z ) * invDirZ );
                minimumOfFarDistances =
                    ( std::min )( minimumOfFarDistances, ( box.minimum().z - ray.origin().z ) * invDirZ );
            }
        }

        if( maximumOfNearDistances >= end || minimumOfFarDistances <= start ) {
            return false;
        }

        return ( maximumOfNearDistances <= minimumOfFarDistances );
    }
}

/*
bool is_intersecting( const frantic::graphics::boundbox3f & bbox, const frantic::graphics::ray3f & ray, const double t0,
const double t1 ) { float tmin, tmax, tymin, tymax, tzmin, tzmax;

    const float divx = 1 / ray.direction().x;
    if (divx >= 0) {
        tmin = divx * (bbox.minimum().x - ray.origin().x);
        tmax = divx * (bbox.maximum().x - ray.origin().x);
    }
    else {
        tmin = divx * (bbox.maximum().x - ray.origin().x);
        tmax = divx * (bbox.minimum().x - ray.origin().x);
    }

    const float divy = 1 / ray.direction().y;
    if (divy >= 0) {
        tymin = divy * (bbox.minimum().y - ray.origin().y);
        tymax = divy * (bbox.maximum().y - ray.origin().y);
    }
    else {
        tymin = divy * (bbox.maximum().y - ray.origin().y);
        tymax = divy * (bbox.minimum().y - ray.origin().y);
    }
  if ( (tmin > tymax) || (tymin > tmax) ) {
        return false;
  }
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    const float divz = 1 / ray.direction().z;
    if (divz >= 0) {
        tzmin = divz * (bbox.minimum().z - ray.origin().z);
        tzmax = divz * (bbox.maximum().z - ray.origin().z);
    }
    else {
        tzmin = divz * (bbox.maximum().z - ray.origin().z);
        tzmax = divz * (bbox.minimum().z - ray.origin().z);
    }
  if ( (tmin > tzmax) || (tzmin > tmax) ) {
        return false;
  }
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return ( (tmin < static_cast<float>( t1 )) && (tmax > static_cast<float>( t0 )) );
}
*/

/**
 *  Determine whether the ray intersects the triangle within the available
 * precision.  If an intersection is found, then this function produces
 * a bounding box which contains the intersection.
 *
 * @param ray the ray to check for intersections.
 * @param t0 the first value of the ray parameter.
 * @param t1 the last value of the ray parameter.
 * @param va the first vertex of the triangle to find intersections with.
 * @param vb the second vertex.
 * @param vc the third vertex.
 * @param outBBox if the function returns true, then this boundbox contains
 *		the intersection between the ray and the triangle.
 * @param lastDiagonalLength this is set during recursion; you should simply
 *		use the default value.
 * @return true if the ray intersects the triangle within the available
 *		precision, and false otherwise.
 *
 * Recursive algorithm for detecting ray-volume intersections as described
 * in:
 *  Dammertz and Keller.  "Improving Ray Tracing Precision by Object Space
 * Intersection Computation."
 */
bool intersect_triangle( const frantic::graphics::ray3f& ray, const double t0, const double t1,
                         const frantic::graphics::vector3f& va, const frantic::graphics::vector3f& vb,
                         const frantic::graphics::vector3f& vc, frantic::graphics::boundbox3f& outBBox,
                         const double lastDiagonalLength ) {
    frantic::graphics::boundbox3f bbox( va );
    bbox += vb;
    bbox += vc;

    const double diagonalLength = bbox.xsize() + bbox.ysize() + bbox.zsize();

    // TODO: the ray-box volume intersection catches some cases that the
    // local implementation was missing.  This may be for cases in which
    // the ray is on the box ?
    // This gets run a lot, so it's probably worth looking into for
    // performance reasons.  We probably want to change the ray parameter
    // to use doubles, too.
    // if( ! ray.is_intersecting_box_volume( bbox, float( t0 ), float( t1 ) ) ) {
    if( !is_intersecting( bbox, ray, t0, t1 ) ) {
        return false;
    } else if( diagonalLength >= lastDiagonalLength ) {
        outBBox = bbox;
        return true;
    }

    /*

    Triangle is recursively subdivided into four new triangles at every step:
         vc
         /\
      vac /__\ vbc
       /\  /\
     va	/__\/__\ vb
        vab

    */

    const vector3f vab = 0.5f * ( va + vb );
    const vector3f vac = 0.5f * ( va + vc );
    const vector3f vbc = 0.5f * ( vb + vc );

    bool hit = intersect_triangle( ray, t0, t1, va, vab, vac, outBBox, diagonalLength );

    if( !hit ) {
        hit = intersect_triangle( ray, t0, t1, vac, vab, vbc, outBBox, diagonalLength );
    }
    if( !hit ) {
        hit = intersect_triangle( ray, t0, t1, vab, vb, vbc, outBBox, diagonalLength );
    }
    if( !hit ) {
        hit = intersect_triangle( ray, t0, t1, vac, vbc, vc, outBBox, diagonalLength );
    }

    return hit;
}

bool mixed_kdtree_trimesh3_primitive::intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                                     mixed_kdtree_ray_observer* observer, bool fixedTime,
                                                     double time ) {
    bool found = false;

    // frantic::logging::set_logging_level( frantic::logging::level::debug );
    // FF_LOG( debug ) << "checking for intersection with face " << m_faceNumber << std::endl;
    // FF_LOG( debug ) << "- current ray: " << ray.str() << " (" << tMin << "," << tMax << ")" << std::endl;

    // Get the intersection with the plane
    // const double distance = m_facePlane.get_distance_to_intersection( ray );
    bool isCoplanar;
    double distance = m_facePlane.get_distance_to_intersection( ray, isCoplanar );
    // FF_LOG( debug ) << "- got distance " << distance << (isCoplanar ? ", is coplanar" : "" ) << std::endl;
    if( isCoplanar ) {
        distance = tMax;
    }

    if( distance >= tMin && distance <= tMax ) {
        // Check that the intersection is inside the bounding box (we're intersecting the ray with the intersection of
        // the triangle and the bounding box)
        const vector3f pt = ray.at( distance );
        // FF_LOG( debug ) << "- point: " << pt << std::endl;

        vector3f barycentricCoord = compute_barycentric_coordinates( pt );
        // FF_LOG( debug ) << "- barycentric coords: " << barycentricCoord << std::endl;
        //  If all the coordinates are positive, it's inside the triangle
        //  TODO: find a better way to do this
        //  For now I've just loosened the barycentric coord tolerances,
        //  to help avoid rays going through the cracks between
        //  adjacent faces.

        if( barycentricCoord.x >= 0 && barycentricCoord.y >= 0 && barycentricCoord.z >= 0 ) {
            // FF_LOG( debug ) << "- got hit" << std::endl;
            /*
            bool inside = true;
            float rayPlucker[6];
            get_plucker( ray.at( tMin ), ray.at( tMax ), rayPlucker );
            bool side0 = plucker_side( m_pluckerA, rayPlucker ) >= 0;
            bool side1 = plucker_side( m_pluckerB, rayPlucker ) >= 0;
            if( side0 != side1 ) {
                inside = false;
            }
            bool side2 = plucker_side( m_pluckerC, rayPlucker ) >= 0;
            if( side2 != side0 ) {
                inside = false;
            }
            */

            // if( observer->wants_hit_details( distance ) ) {
            mixed_kdtree_trimesh3_point_data intersection( &m_helper.get_mesh_ref(), m_faceNumber, barycentricCoord );

            if( fixedTime ) {
                observer->insert( ray, distance, time, m_facePlane.normal(), this,
                                  reinterpret_cast<char*>( &intersection ),
                                  sizeof( mixed_kdtree_trimesh3_point_data ) );
            } else {
                observer->insert( ray, distance, distance, m_facePlane.normal(), this,
                                  reinterpret_cast<char*>( &intersection ),
                                  sizeof( mixed_kdtree_trimesh3_point_data ) );
            }

            // outIntersection.reset( ray, distance, time, m_facePlane.normal(), this, reinterpret_cast<char *>(&
            // intersection), sizeof( mixed_kdtree_trimesh3_point_data ) );

            found = true;
        } else if( barycentricCoord.x >= -0.001f && barycentricCoord.y >= -0.001f && barycentricCoord.z >= -0.001f ) {
            // FF_LOG( debug ) << "- robust intersection check" << std::endl;
            //  If it's outside the triangle according to the barycentric
            //  test, but it's within an error margin, then try a more
            //  robust check.
            //
            //  Before, I was using a moderate epsilon ( 1.e-4 ) without
            //  a robust check.  This got rid of the false negatives I
            //  was seeing, but added false positives.  I'd imagine such small
            //  false positives aren't much of a problem for rendering, but
            //  they correspond to overhanging faces which had some weird
            //  effects during particle collisions.
            boundbox3f outBBox;
            if( intersect_triangle( ray, tMin, tMax, m_helper.get_mesh_ref().get_vertex( m_face.x ),
                                    m_helper.get_mesh_ref().get_vertex( m_face.y ),
                                    m_helper.get_mesh_ref().get_vertex( m_face.z ), outBBox ) ) {
                // FF_LOG( debug ) << "- got hit" << std::endl;
                mixed_kdtree_trimesh3_point_data intersection( &m_helper.get_mesh_ref(), m_faceNumber,
                                                               barycentricCoord );

                if( fixedTime ) {
                    observer->insert( ray, distance, time, m_facePlane.normal(), this,
                                      reinterpret_cast<char*>( &intersection ),
                                      sizeof( mixed_kdtree_trimesh3_point_data ) );
                } else {
                    observer->insert( ray, distance, distance, m_facePlane.normal(), this,
                                      reinterpret_cast<char*>( &intersection ),
                                      sizeof( mixed_kdtree_trimesh3_point_data ) );
                }

                found = true;
            }
        }
    }
    return found;
}

bool mixed_kdtree_trimesh3_primitive::intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                                     mixed_kdtree_ray_observer* observer, double time ) {
    return intersect_ray( ray, tMin, tMax, observer, true, time );
}

bool mixed_kdtree_trimesh3_primitive::intersect_particle( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                                          mixed_kdtree_ray_observer* observer ) {
    return intersect_ray( ray, tMin, tMax, observer, false, 0 );
}

void mixed_kdtree_trimesh3_primitive::get_any_point( vector3f& outPoint, vector3f& outNormal, double /*time*/ ) {
    outPoint = ( 1.0f / 3 *
                 ( m_helper.get_mesh_ref().get_vertex( m_face.x ) + m_helper.get_mesh_ref().get_vertex( m_face.y ) +
                   m_helper.get_mesh_ref().get_vertex( m_face.z ) ) );
    outNormal = m_facePlane.normal();
}

// todo: should the output be restricted to points that are
// within nodebounds ?
bool mixed_kdtree_trimesh3_primitive::find_nearest_point( const frantic::graphics::vector3f& point,
                                                          frantic::graphics::boundbox3f& /*nodeBounds*/,
                                                          mixed_kdtree_distance_observer* distanceObserver,
                                                          double time ) {
    bool found = false;

    const plane3f plane = m_facePlane; // m_mesh.get_face_plane( m_faceNumber );

    // Project the point onto the triangle's plane
    double distance = plane.get_signed_distance_to_plane_double( point );

    if( distanceObserver->is_in_range( fabs( distance ) ) ) {

        // Check that the intersection is inside the bounding box (we're intersecting the ray with the intersection of
        // the triangle and the bounding box)
        vector3f projected = point - static_cast<float>( distance ) * plane.normal();

        // Now get the barycentric coordinates of this projection
        vector3f barycentricCoord = compute_barycentric_coordinates( projected );

        if( barycentricCoord.is_inf() )
            return false;

        if( barycentricCoord.x >= 0 && barycentricCoord.y >= 0 && barycentricCoord.z >= 0 ) {
            distance = fabs( distance ); // Projected and barycentricCoord are both accurate
        } else {
            const vector3f& A = m_helper.get_mesh_ref().get_vertex( m_face.x );
            const vector3f& B = m_helper.get_mesh_ref().get_vertex( m_face.y );
            const vector3f& C = m_helper.get_mesh_ref().get_vertex( m_face.z );

            detail::nearest_point_on_triangle( projected, barycentricCoord, A, B, C );
            distance = vector3f::distance_double( point, projected );
        }

        // if( distanceObserver->wants_hit_details( distance ) ) {
        mixed_kdtree_trimesh3_point_data pointData( &m_helper.get_mesh_ref(), m_faceNumber, barycentricCoord );

        distanceObserver->insert( projected, distance, time, plane.normal(), this,
                                  reinterpret_cast<char*>( &pointData ), sizeof( mixed_kdtree_trimesh3_point_data ) );

        found = true;
        //}
    }
    return found;
}

frantic::graphics::boundbox3f mixed_kdtree_trimesh3_primitive::get_bounds() const {
    const vector3& face( m_face /*m_helper.get_mesh_ref().get_face(m_faceNumber)*/ );

    frantic::graphics::boundbox3f result;
    result.set_to_point( m_helper.get_mesh_ref().get_vertex( face.x ) );
    result += m_helper.get_mesh_ref().get_vertex( face.y );
    result += m_helper.get_mesh_ref().get_vertex( face.z );

    return result;
}

boundbox3f mixed_kdtree_trimesh3_primitive::intersect_with( const boundbox3f& box ) const {
    const frantic::graphics::vector3& face( m_face /*m_helper.get_mesh_ref().get_face(m_faceNumber)*/ );

    const vector3f& A( m_helper.get_mesh_ref().get_vertex( face.x ) );
    const vector3f& B( m_helper.get_mesh_ref().get_vertex( face.y ) );
    const vector3f& C( m_helper.get_mesh_ref().get_vertex( face.z ) );

    boundbox3f result;
    box.intersect_with_triangle( A, B, C, result );
    return result;
}

// I'm not sure if this makes sense, but I am implementing it for testing purposes
boost::int32_t mixed_kdtree_trimesh3_primitive::get_face_number( void ) const {
    return static_cast<boost::int32_t>( m_faceNumber );
}

vector3f mixed_kdtree_trimesh3_primitive::get_barycentric_coord( const char* data ) const {
    const mixed_kdtree_trimesh3_point_data* pointData(
        reinterpret_cast<const mixed_kdtree_trimesh3_point_data*>( data ) );
    return pointData->barycentricCoords;
}

vector3f mixed_kdtree_trimesh3_primitive::get_geometric_normal( const char* ) const { return m_facePlane.normal(); }

std::size_t mixed_kdtree_trimesh3_primitive::get_data_size( void ) const {
    return sizeof( mixed_kdtree_trimesh3_point_data );
}

void mixed_kdtree_trimesh3_primitive::to_raytrace_intersection( raytrace_intersection& raytraceIntersection,
                                                                const char* data ) const {
    const mixed_kdtree_trimesh3_point_data* pointData(
        reinterpret_cast<const mixed_kdtree_trimesh3_point_data*>( data ) );
    raytraceIntersection.faceIndex = static_cast<int>( pointData->faceNumber );
    raytraceIntersection.barycentricCoords = pointData->barycentricCoords;
    raytraceIntersection.motionDuringTimeStep = vector3f( 0 );
}

void mixed_kdtree_trimesh3_primitive::to_nearest_point_search_result( nearest_point_search_result& searchResult,
                                                                      const char* data ) const {
    const mixed_kdtree_trimesh3_point_data* pointData(
        reinterpret_cast<const mixed_kdtree_trimesh3_point_data*>( data ) );
    searchResult.faceIndex = static_cast<int>( pointData->faceNumber ); // should be m_faceNumber too
    searchResult.barycentricCoords = pointData->barycentricCoords;
}

void mixed_kdtree_trimesh3_primitive::to_channel_map_data( char* outData,
                                                           const mixed_kdtree_point_data* pointData ) const {
    m_helper.to_channel_map_data( outData, this, pointData );
}

} // namespace geometry
} // namespace frantic
