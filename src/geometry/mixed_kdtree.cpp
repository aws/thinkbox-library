// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/mixed_kdtree.hpp>
#include <frantic/geometry/mixed_kdtree_trimesh3.hpp>

#include <frantic/geometry/raytracing.hpp>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include <algorithm>

#include <math.h>

using namespace frantic;
using namespace frantic::graphics;
using namespace frantic::channels;
using namespace frantic::geometry;

namespace frantic {
namespace geometry {

//
// mixed_kdtree
//

void mixed_kdtree::insert_default_channels( frantic::channels::channel_map& cm ) {
    cm.define_channel<float>( _T("Distance") );
    cm.define_channel<vector3f>( _T("Position") );
    cm.define_channel<vector3f>( _T("Normal") );
}

void mixed_kdtree::build_default_channel_map( frantic::channels::channel_map& cm ) {
    cm.reset();
    insert_default_channels( cm );
    cm.end_channel_definition();
}

/*
std::size_t mixed_kdtree::get_max_internal_result_size( void ) const {
    if( ! m_final ) {
        throw std::runtime_error( "mixed_kdtree Error: cannot get the result size until the tree is finalized." );
    }
    return m_maxDataSize + std::max( sizeof( mixed_kdtree_nearest_point ), sizeof( mixed_kdtree_ray_intersection ) );
}
*/
void mixed_kdtree::initialize( void ) {
    m_final = false;
    // m_rootNode = new mixed_kdtree_node<mixed_kdtree_primitive>();
    m_rootNode = new mixed_kdtree_node();
    // m_setChannelMapListeners.push_back( & m_commonNamedChannelSetters );
    /*
      m_logStream.open( "c:/temp/kdtree_log.txt", std::ios_base::out + std::ios_base::app );
    m_outFile.open( "c:/temp/kdtree_log_file.txt", std::ios_base::out + std::ios_base::app );
      frantic::logging::redirect_all_streams( & m_logStream );
      frantic::logging::set_logging_level( frantic::logging::level::debug );
    m_outFile << "Starting log...\n";
    */
}

mixed_kdtree::mixed_kdtree() {
    initialize();

    build_default_channel_map( m_channelMap );
    set_channel_map( m_channelMap );
}

mixed_kdtree::mixed_kdtree( const frantic::channels::channel_map& cm ) {
    initialize();

    set_channel_map( cm );
}

mixed_kdtree::mixed_kdtree( boost::shared_ptr<frantic::geometry::trimesh3> mesh ) {
    // m_outFile << "starting mixed_kdtree from trimesh3...\n";
    initialize();

    build_default_channel_map( m_channelMap );
    set_channel_map( m_channelMap );

    insert( mesh );

    finalize();
    // m_outFile << "done mixed_kdtree from trimesh3.\n";
}

mixed_kdtree::~mixed_kdtree() {
    // m_logStream.close();
    if( m_rootNode != 0 ) {
        delete m_rootNode;
        m_rootNode = 0;
    }
}

void mixed_kdtree::set_channel_map( const frantic::channels::channel_map& channelMap ) {
    // should this be set at construction time instead so we know it before
    // inserting primitives ?
    // if( m_final ) {
    // throw std::runtime_error( "mixed_kdtree Error: cannot set the channel map of a finalized tree." );
    //}
    m_channelMap = channelMap;
    m_commonNamedChannelSetters.set_channel_map( channelMap );
    BOOST_FOREACH( output_channel_map_listener* listener, m_outputChannelMapListeners ) {
        listener->set_output_channel_map( channelMap );
    }
    // need to reset the channel mappings used by primitives
    // so we need to hold an object which controls the conversion
    // from a primitive's plain data to the channel data

    // for now this is being done by each geometry object

    // reset the default output channel data to all zeros
    m_defaultChannelMapData.clear();
    m_defaultChannelMapData.resize( m_channelMap.structure_size(), 0 );
}

const frantic::channels::channel_map& mixed_kdtree::get_channel_map( void ) const { return m_channelMap; }

void mixed_kdtree::set_default_channel_map_data( const std::vector<char>& defaultData ) {
    if( defaultData.size() != m_channelMap.structure_size() ) {
        throw std::runtime_error(
            "mixed_kdtree::set_default_channel_map_data Error: the specified default data does not "
            "match the size of the channel map." );
    }
    m_defaultChannelMapData = defaultData;
}

void mixed_kdtree::insert( boost::shared_ptr<mixed_kdtree_primitive_creator> object ) {
    if( m_final ) {
        throw std::runtime_error( "mixed_kdtree Error: cannot insert into finalized tree" );
    }

    m_objects.push_back( object );
}

void mixed_kdtree::insert( boost::shared_ptr<frantic::geometry::trimesh3> mesh ) {
    boost::shared_ptr<mixed_kdtree_trimesh3> meshObject( new mixed_kdtree_trimesh3( mesh ) );
    insert( meshObject );
}

void mixed_kdtree::insert_with_constant_velocity( boost::shared_ptr<frantic::geometry::trimesh3> mesh, const double t0,
                                                  const double t1 ) {
    boost::shared_ptr<mixed_kdtree_trimesh3_constant_velocity> meshObject(
        new mixed_kdtree_trimesh3_constant_velocity( mesh, t0, t1 ) );
    insert( meshObject );
}

void mixed_kdtree::finalize( void ) {
    // m_outFile << "starting finalize...\n";
    // m_outFile.flush();
    std::size_t primitiveCount = 0;

    BOOST_FOREACH( boost::shared_ptr<mixed_kdtree_primitive_creator>& object, m_objects ) {
        primitiveCount += object->get_primitive_count();
    }

    // Get the primitives and output_channel_map_listeners from
    // all of the primitive_creators.
    mixed_kdtree_primitive_recorder recorder( m_primitives, m_outputChannelMapListeners );
    recorder.reserve( primitiveCount );

    BOOST_FOREACH( boost::shared_ptr<mixed_kdtree_primitive_creator>& object, m_objects ) {
        object->get_primitives( &recorder );
    }

    // Now m_primitives and m_outputChannelMapListeners are
    // filled.

    m_rootBounds = recorder.get_bounds();
    m_rootBounds.expand_fractional( 0.00001f );

    m_maxDataSize = recorder.get_data_size();

    build_kdtree_greedy_SAH_nlogn( *m_rootNode, m_primitives, m_rootBounds );

    m_final = true;

    BOOST_FOREACH( output_channel_map_listener* listener, m_outputChannelMapListeners ) {
        listener->set_output_channel_map( m_channelMap );
    }

    // possible alternative to vectors for primitives' search result data ?
    // m_internalResultPool = boost::shared_ptr<boost::pool<> >( new boost::pool<>( get_max_internal_result_size() ) );
    // m_outFile << "done finalize.\n";
    // m_outFile.flush();
}

boundbox3f mixed_kdtree::get_bounds() const {
    if( !m_final ) {
        throw std::runtime_error( "mixed_kdtree::get_bounds Error: cannot get bounds before finalizing tree" );
    }
    return m_rootBounds;
}

std::size_t mixed_kdtree::get_node_count() const { return m_rootNode->get_node_count(); }

std::size_t mixed_kdtree::get_maximum_depth() const { return m_rootNode->get_maximum_depth(); }

std::size_t mixed_kdtree::get_largest_leaf_size() const { return m_rootNode->get_largest_leaf_size(); }

std::size_t mixed_kdtree::get_leaf_count() const { return m_rootNode->get_leaf_count(); }

std::size_t mixed_kdtree::get_populated_leaf_count() const { return m_rootNode->get_populated_leaf_count(); }

void mixed_kdtree::dump_tree( std::ostream& out ) const {
    out << "----- Dumping mixed KD-Tree\n";
    out << " boundbox: " << m_rootBounds << "\n";
    m_rootNode->dump_tree( out, 0 );
    out << "----- Finished dumping KD-Tree" << std::endl;
}

bool mixed_kdtree::intersect_particle( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                       raytrace_intersection& outIntersection ) const {
    assert( m_rootNode );

    // FF_LOG( debug ) << "intersect_particle " << ray << " tMin " << tMin << " tMax " << tMax << "\n";
    // if( ! m_rootBounds.is_volume_intersecting_line_segment( ray.at( tMin ), ray.at( tMax ) ) )
    // return false;

    mixed_kdtree_ray_intersection intersection;
    find_nearest_intersection_mixed_kdtree_ray_observer observer( intersection );
    // m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_particle<mixed_kdtree_primitive>() );
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_particle() );
    bool found = observer.got_hit();
    if( found ) {
        intersection.to_raytrace_intersection( outIntersection );
        // FF_LOG( debug ) << "got hit\n";
    }
    return found;
}

bool mixed_kdtree::intersect_particle( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                       char* outData ) const {
    assert( m_rootNode );

    // if( ! m_rootBounds.is_volume_intersecting_line_segment( ray.at( tMin ), ray.at( tMax ) ) )
    // return false;

    mixed_kdtree_ray_intersection intersection;
    find_nearest_intersection_mixed_kdtree_ray_observer observer( intersection );
    // m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_particle<mixed_kdtree_primitive>() );
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_particle() );
    bool found = observer.got_hit();
    if( found ) {
        intersection.to_channel_map_data( outData, m_commonNamedChannelSetters, m_defaultChannelMapData );
    }
    return found;
}

// This returns the first intersection between the ray and one of the primitives, or false if none is found.
bool mixed_kdtree::intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                  raytrace_intersection& outIntersection, double time ) const {
    assert( m_rootNode );

    // if( ! m_rootBounds.is_volume_intersecting_line_segment( ray.at( tMin ), ray.at( tMax ) ) )
    // return false;

    mixed_kdtree_ray_intersection intersection;
    find_nearest_intersection_mixed_kdtree_ray_observer observer( intersection );
    // m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray<mixed_kdtree_primitive>( time ) );
    // FF_LOG( debug ) << "intersect_ray called with ray " << ray.str() << " (" << tMin << "," << tMax << ")" <<
    // std::endl;
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray( time ) );
    bool found = observer.got_hit();
    if( found ) {
        intersection.to_raytrace_intersection( outIntersection );
        // FF_LOG( debug ) << "found intersection at " << outIntersection.position.str() << std::endl;
    }
    return found;
}

bool mixed_kdtree::intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax, char* outData,
                                  double time ) const {
    assert( m_rootNode );

    // if( ! m_rootBounds.is_volume_intersecting_line_segment( ray.at( tMin ), ray.at( tMax ) ) )
    // return false;

    mixed_kdtree_ray_intersection intersection;
    find_nearest_intersection_mixed_kdtree_ray_observer observer( intersection );
    // m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray<mixed_kdtree_primitive>( time ) );
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray( time ) );
    bool found = observer.got_hit();
    if( found ) {
        intersection.to_channel_map_data( outData, m_commonNamedChannelSetters, m_defaultChannelMapData );
    }
    return found;
}

bool mixed_kdtree::intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                  std::vector<char>& outData, double time ) const {
    assert( m_rootNode );

    const std::size_t dataSize = get_channel_map().structure_size();
    if( outData.size() != dataSize ) {
        outData.resize( get_channel_map().structure_size() );
    }

    mixed_kdtree_ray_intersection intersection;
    find_nearest_intersection_mixed_kdtree_ray_observer observer( intersection );
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray( time ) );
    bool found = observer.got_hit();
    if( found && dataSize != 0 ) {
        intersection.to_channel_map_data( &outData[0], m_commonNamedChannelSetters, m_defaultChannelMapData );
    }
    return found;
}

// This returns true if there is an intersection between the ray and one of the primitives, or false if none is found.
bool mixed_kdtree::intersects_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax, double time ) const {
    assert( m_rootNode );

    // if( ! m_rootBounds.is_volume_intersecting_line_segment( ray.at( tMin ), ray.at( tMax ) ) )
    // return false;

    // std::cout << "called intersects_ray\n";
    find_any_intersection_mixed_kdtree_ray_observer observer;
    // m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray<mixed_kdtree_primitive>( time ) );
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray( time ) );
    return observer.got_hit();
}

bool mixed_kdtree::intersects_ray_segment( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                           double time ) const {
    return intersects_ray( ray, tMin, tMax, time );
}

// This returns true if there is an intersection between the ray and one of the primitives, or false if none is found.
bool mixed_kdtree::intersects_ray_segment( const frantic::graphics::vector3f& start,
                                           const frantic::graphics::vector3f& end, double time ) const {
    assert( m_rootNode );

    // if( ! m_rootBounds.is_volume_intersecting_line_segment( start, end ) )
    // return false;

    // std::cout << "called intersects_ray\n";
    ray3f ray( start, end - start );
    find_any_intersection_mixed_kdtree_ray_observer observer;
    // m_rootNode->traverse_ray( ray, 0, 1.f, observer, do_intersect_ray<mixed_kdtree_primitive>( time ) );
    // std::cout << "traversing " << ray << ", from " << start << " to " << end << "\n";
    m_rootNode->traverse_ray( ray, 0, 1.f, observer, do_intersect_ray( time ) );
    return observer.got_hit();
}

// only tested indirectly
void mixed_kdtree::intersect_ray_all( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                      std::vector<mixed_kdtree_ray_intersection>& outIntersections,
                                      double time ) const {
    assert( m_rootNode );

    outIntersections.clear();

    // if( ! m_rootBounds.is_volume_intersecting_line_segment( ray.at( tMin ), ray.at( tMax ) ) )
    // return;

    find_all_intersections_mixed_kdtree_ray_observer observer( outIntersections );
    // m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray<mixed_kdtree_primitive>( time ) );
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray( time ) );
}

void mixed_kdtree::intersect_ray_all( const frantic::graphics::ray3f& ray, double tMax,
                                      std::vector<raytrace_intersection>& outIntersections, double time ) const {
    std::vector<mixed_kdtree_ray_intersection> intersections;

    intersect_ray_all( ray, 0, static_cast<double>( tMax ), intersections, time );
    // outIntersections is NOT cleared in the trimesh3_kdtree method
    // so don't clear it here either
    outIntersections.reserve( outIntersections.size() + intersections.size() );
    BOOST_FOREACH( const mixed_kdtree_ray_intersection& intersection, intersections ) {
        raytrace_intersection outIntersection;
        intersection.to_raytrace_intersection( outIntersection );
        outIntersections.push_back( outIntersection );
    }
    std::sort( outIntersections.begin(), outIntersections.end() );
}

void mixed_kdtree::intersect_ray_all( const frantic::graphics::ray3f& ray, double tMax,
                                      particles::particle_array& outIntersections, double time ) const {
    std::vector<mixed_kdtree_ray_intersection> intersections;

    intersect_ray_all( ray, 0, static_cast<double>( tMax ), intersections, time );

    // std::cout << "got " << intersections.size() << " intersections\n";
    // outIntersections is NOT cleared in the trimesh3_kdtree method
    // so don't clear it here either
    std::size_t outIndex = outIntersections.size();
    outIntersections.resize( outIntersections.size() + intersections.size() );
    BOOST_FOREACH( mixed_kdtree_ray_intersection& intersection, intersections ) {
        intersection.to_channel_map_data( outIntersections.at( outIndex ), m_commonNamedChannelSetters,
                                          m_defaultChannelMapData );
        ++outIndex;
    }
    // should this be sorted..?
    // we don't necessarily have distance information for the former intersections
    // so we can't always guarantee that it's sorted
    // we could force the existence of the distance field in the output cchannel
    // map or something like that if sorting is necessary
    // std::sort( outIntersections.begin(), outIntersections.end() );
}

bool mixed_kdtree::find_nearest_point( const vector3f& point, float maxDistance,
                                       nearest_point_search_result& outNearestPoint, double time ) const {
    assert( m_rootNode );

    // temporary nearest point
    mixed_kdtree_nearest_point nearestPoint;
    nearest_mixed_kdtree_distance_observer observer( nearestPoint, maxDistance );
    boundbox3f bounds( m_rootBounds );
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );
    if( observer.got_hit() ) {
        nearestPoint.to_nearest_point_search_result( outNearestPoint );
        return true;
    }
    return false;
}

bool mixed_kdtree::find_nearest_point( const vector3f& point, float maxDistance, vector3f* pOutNearestPoint,
                                       float* pOutNearestDistance, vector3f* pOutNormal, double time ) const {
    assert( m_rootNode );

    mixed_kdtree_nearest_point nearestPoint;
    boundbox3f bounds( m_rootBounds );
    nearest_mixed_kdtree_distance_observer observer( nearestPoint, maxDistance );
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );
    if( observer.got_hit() ) {
        if( pOutNearestDistance )
            *pOutNearestDistance = static_cast<float>( nearestPoint.get_distance() );
        if( pOutNearestPoint )
            *pOutNearestPoint = nearestPoint.get_position();
        if( pOutNormal )
            *pOutNormal = nearestPoint.get_normal();
        return true;
    }
    if( pOutNearestDistance )
        *pOutNearestDistance = std::numeric_limits<float>::max();
    return false;
}

bool mixed_kdtree::find_nearest_point( const vector3f& point, float maxDistance, char* outData, double time ) const {
    assert( m_rootNode );

    // temporary nearest point
    mixed_kdtree_nearest_point nearestPoint;
    nearest_mixed_kdtree_distance_observer observer( nearestPoint, maxDistance );
    boundbox3f bounds( m_rootBounds );
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );
    if( observer.got_hit() ) {
        nearestPoint.to_channel_map_data( outData, m_commonNamedChannelSetters, m_defaultChannelMapData );
        return true;
    }
    return false;
}

bool mixed_kdtree::find_nearest_point( const vector3f& point, float maxDistance, std::vector<char>& outData,
                                       double time ) const {
    assert( m_rootNode );

    const std::size_t dataSize = get_channel_map().structure_size();
    if( outData.size() != dataSize ) {
        outData.resize( get_channel_map().structure_size() );
    }

    // temporary nearest point
    mixed_kdtree_nearest_point nearestPoint;
    nearest_mixed_kdtree_distance_observer observer( nearestPoint, maxDistance );
    boundbox3f bounds( m_rootBounds );
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );
    if( observer.got_hit() && dataSize > 0 ) {
        nearestPoint.to_channel_map_data( &outData[0], m_commonNamedChannelSetters, m_defaultChannelMapData );
        return true;
    }
    return false;
}

void mixed_kdtree::collect_nearest_primitives( const vector3f& point, std::size_t nFaces,
                                               std::vector<nearest_point_search_result>& outResults,
                                               double time ) const {
    assert( m_rootNode );

    outResults.clear();
    std::vector<mixed_kdtree_nearest_point> nearestPoints;
    within_range_mixed_kdtree_distance_observer observer( nearestPoints, std::numeric_limits<float>::max(), nFaces );

    boundbox3f bounds = m_rootBounds;
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );

    BOOST_FOREACH( mixed_kdtree_nearest_point& nearestPoint, nearestPoints ) {
        nearest_point_search_result result;
        nearestPoint.to_nearest_point_search_result( result );
        outResults.push_back( result );
    }
}

void mixed_kdtree::collect_primitives_within_range( const vector3f& point, float maxDistance,
                                                    std::vector<nearest_point_search_result>& outResults,
                                                    double time ) const {
    assert( m_rootNode );

    outResults.clear();
    std::vector<mixed_kdtree_nearest_point> nearestPoints;
    within_range_mixed_kdtree_distance_observer observer( nearestPoints, maxDistance );

    boundbox3f bounds = m_rootBounds;
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );

    BOOST_FOREACH( mixed_kdtree_nearest_point& nearestPoint, nearestPoints ) {
        nearest_point_search_result result;
        nearestPoint.to_nearest_point_search_result( result );
        outResults.push_back( result );
    }
}

void mixed_kdtree::collect_nearest_primitives( const vector3f& point, std::size_t nFaces,
                                               particles::particle_array& outResults, double time ) const {
    assert( m_rootNode );

    outResults.clear();
    std::vector<mixed_kdtree_nearest_point> nearestPoints;
    within_range_mixed_kdtree_distance_observer observer( nearestPoints, std::numeric_limits<float>::max(), nFaces );

    boundbox3f bounds = m_rootBounds;
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );

    outResults.resize( nearestPoints.size() );

    for( std::size_t i = 0; i < nearestPoints.size(); ++i ) {
        nearestPoints[i].to_channel_map_data( outResults.at( i ), m_commonNamedChannelSetters,
                                              m_defaultChannelMapData );
    }
}

void mixed_kdtree::collect_primitives_within_range( const vector3f& point, float maxDistance,
                                                    particles::particle_array& outResults, double time ) const {
    assert( m_rootNode );

    outResults.clear();
    std::vector<mixed_kdtree_nearest_point> nearestPoints;
    within_range_mixed_kdtree_distance_observer observer( nearestPoints, maxDistance );

    boundbox3f bounds = m_rootBounds;
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );

    outResults.resize( nearestPoints.size() );

    for( std::size_t i = 0; i < nearestPoints.size(); ++i ) {
        nearestPoints[i].to_channel_map_data( outResults.at( i ), m_commonNamedChannelSetters,
                                              m_defaultChannelMapData );
    }
}

// This collects all triangles indices within the sphere with center "point" and radius "maxDistance"
void mixed_kdtree::collect_faces_within_range( const vector3f& point, float maxDistance,
                                               std::vector<int>& outFaceNumbers, double time ) const {
    assert( m_rootNode );

    outFaceNumbers.clear();
    std::vector<mixed_kdtree_nearest_point> nearestPoints;
    within_range_mixed_kdtree_distance_observer observer( nearestPoints, maxDistance );

    boundbox3f bounds = m_rootBounds;
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );

    // the primitives should already be unique and sorted due to the observer
    // If not, we should sort them and remove duplicates here, as in
    // the trimesh3_kdtree procedure.
    BOOST_FOREACH( mixed_kdtree_nearest_point& point, nearestPoints ) {
        nearest_point_search_result result;
        point.to_nearest_point_search_result( result );
        if( result.faceIndex >= 0 ) {
            outFaceNumbers.push_back( result.faceIndex );
        }
    }
}

// This collects all triangles indices within the sphere with center "point" and radius "maxDistance"
void mixed_kdtree::collect_primitives_within_range( const vector3f& point, float maxDistance,
                                                    std::vector<const mixed_kdtree_primitive*>& outPrimitives,
                                                    double time ) const {
    assert( m_rootNode );

    outPrimitives.clear();
    std::vector<mixed_kdtree_nearest_point> nearestPoints;
    within_range_mixed_kdtree_distance_observer observer( nearestPoints, maxDistance );

    boundbox3f bounds = m_rootBounds;
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );

    // the primitives should already be unique and sorted due to the observer
    // If not, we should sort them and remove duplicates here, as in
    // the trimesh3_kdtree procedure.
    BOOST_FOREACH( mixed_kdtree_nearest_point& point, nearestPoints ) {
        outPrimitives.push_back( point.get_primitive() );
    }
}

/*
// I think this is covered fine by find_nearest_point ?
float mixed_kdtree::get_distance_to_surface(const vector3f & point, vector3f * pOutNearestPoint, vector3f * pOutNormal,
double time ) const {
  assert( m_rootNode );

    mixed_kdtree_nearest_point nearestPoint;
    boundbox3f bounds( m_rootBounds );
    nearest_mixed_kdtree_distance_observer observer( nearestPoint, std::numeric_limits<float>::max() );
    m_rootNode->traverse_nearest_points( point, bounds, observer, time );
    if( observer.got_hit() ) {
        if(pOutNearestPoint)
            *pOutNearestPoint = nearestPoint.get_position();
        if(pOutNormal)
            *pOutNormal = nearestPoint.get_normal();
        return static_cast<float>( nearestPoint.get_distance() );
    }
    return std::numeric_limits<float>::max();
}
*/

// from volume_collection
//
// I am using the same voting scheme due to intersections with
// edges / vertices.
//
// Unfortunately even if we get rid of cracks between faces
// we cannot simply use the normal of nearby faces, in case the
// nearest point is an edge or a vertex.
// We could get rid of the voting and use a nearby-primitive search
// instead if we add weighed normal calculations for edges and
// vertices.  We can precalculate such normals easily, but I don't know
// if the overhead would be worth it.
bool mixed_kdtree::is_point_in_volume( const vector3f& p, const double time ) const {
    if( !m_rootBounds.contains( p ) )
        return false;

    int result = 0;

    for( std::size_t i = 0; i < m_primitives.size(); ++i ) {
        vector3f hitPoint, hitNormal;

        // find a point in the geometry we can shoot toward
        m_primitives[i]->get_any_point( hitPoint, hitNormal, time );

        // Verify this point doesn't lie in the triangle's plane.
        // TODO: This should be a relative error
        // vector3f norm = vector3f::normalize( vector3f::cross( v1 - v0, v2 - v0 ) );
        if( std::abs( vector3f::dot( hitNormal, ( p - hitPoint ).to_normalized() ) ) < 0.001f )
            continue;

        raytrace_intersection ri;

        if( intersect_ray( ray3f( p, hitPoint - p ), 0.001, 1.001, ri, time ) ) {
            if( vector3f::dot( ri.geometricNormal, ri.ray.direction() ) > 0 ) {
                if( ++result > 1 )
                    return true;
            } else {
                if( --result < -1 )
                    return false;
            }
        }
    }

    // If all the faces were skipped, then the point is co-planar with the
    // entire mesh (or the mesh is empty).
    return false;
}

// from volume_collection
// I made some changes to try to handle edges better
float mixed_kdtree::get_segment_distance_in_volume( const vector3f& p1, const vector3f& p2, const double time ) const {
    const float eps = 1.0e-5f;

    vector3f rayDir( p2 - p1 );
    vector3f normalizedRayDir( rayDir.to_normalized() );
    float distance = rayDir.get_magnitude();

    // accumulator for distance inside volume
    // this is strictly increasing
    // this should end at <= 1.0, barring precision problems
    double accInside = 0;

    // the intersection.distance at the start of the currentRegion
    double currentRegionStart = 0;

    enum region_type { UNKNOWN, INSIDE, OUTSIDE };

    // the region that the ray is inside now
    region_type currentRegion = UNKNOWN;

    // the region that the ray is entering at the intersection
    region_type newRegion = UNKNOWN;

    std::vector<raytrace_intersection> intersections;
    intersect_ray_all( ray3f( p1, rayDir ), 1, intersections, time );

    // I assume here that the intersections are in order of increasing
    // distance, which they should be.
    BOOST_FOREACH( raytrace_intersection& intersection, intersections ) {
        const float projection( vector3f::dot( intersection.geometricNormal, normalizedRayDir ) );
        if( projection > eps ) {
            // going outside
            newRegion = OUTSIDE;
        } else if( projection < -eps ) {
            // going inside
            newRegion = INSIDE;
        }

        // if the region has changed, add the old region to
        // accInside
        if( newRegion != currentRegion ) {
            switch( currentRegion ) {
            case UNKNOWN:
                switch( newRegion ) {
                // assume the currentRegion was the opposite of the new region
                // alternatively we could verify this with is_point_in_volume()

                // currentRegion == UNKNOWN, newRegion == INSIDE
                case INSIDE:
                    // assume currentRegion == OUTSIDE
                    // accInside -= intersection.distance;
                    break;

                // currentRegion == UNKNOWN, newRegion == OUTSIDE
                case OUTSIDE:
                    // assume currentRegion == INSIDE
                    accInside += intersection.distance;
                    break;
                case UNKNOWN:
                    break;
                }
                break;

            case INSIDE:
                accInside += ( intersection.distance - currentRegionStart );
                break;

            case OUTSIDE:
                // accInside -= intersection.distance;
                break;
            }

            currentRegion = newRegion;
            currentRegionStart = intersection.distance;
        }
    }

    switch( currentRegion ) {
    case UNKNOWN:
        return is_point_in_volume( p1, time ) ? distance : 0;
    case INSIDE:
        accInside += ( 1.0 - currentRegionStart );
        break;
    case OUTSIDE:
        // accInside -= 1.0;
        break;
    }

    // this should be unnecessary..
    accInside = math::clamp( accInside, 0.0, 1.0 );

    return static_cast<float>( accInside * distance );
}

//
// distance traversals
//

// determine the single closest distance
least_distance_mixed_kdtree_distance_observer::least_distance_mixed_kdtree_distance_observer()
    : m_distance( std::numeric_limits<double>::max() ) {}
bool least_distance_mixed_kdtree_distance_observer::is_done( void ) { return false; }
double least_distance_mixed_kdtree_distance_observer::get_max_range( void ) { return m_distance; }
bool least_distance_mixed_kdtree_distance_observer::is_in_range( double d ) { return d < get_max_range(); }
// bool least_distance_mixed_kdtree_distance_observer::wants_hit_details( double d ) {
//    if( d < m_distance ) {
//        m_distance = d;
//    }
//    return false;
//}
void least_distance_mixed_kdtree_distance_observer::insert( const vector3f& /*pt*/, const double d,
                                                            const double /*time*/, const vector3f& /*normal*/,
                                                            mixed_kdtree_primitive* /*mixed*/, const char* /*data = 0*/,
                                                            const std::size_t /*dataSize = 0*/ ) {
    if( d < m_distance ) {
        m_distance = d;
    }
}

double least_distance_mixed_kdtree_distance_observer::get_distance( void ) { return m_distance; }

///

// determine the single closest point and primitive
nearest_mixed_kdtree_distance_observer::nearest_mixed_kdtree_distance_observer( mixed_kdtree_nearest_point& outResult,
                                                                                const double maxDistance )
    : m_maxDistance( maxDistance )
    , m_gotHit( false )
    , m_nearestPoint( outResult ) {}
bool nearest_mixed_kdtree_distance_observer::is_done( void ) { return false; }
bool nearest_mixed_kdtree_distance_observer::got_hit( void ) { return m_gotHit; }
double nearest_mixed_kdtree_distance_observer::get_max_range( void ) { return m_maxDistance; }
bool nearest_mixed_kdtree_distance_observer::is_in_range( double d ) { return d < get_max_range(); }
// bool nearest_mixed_kdtree_distance_observer::wants_hit_details( double d ) { return d < get_max_range(); }
void nearest_mixed_kdtree_distance_observer::insert( const vector3f& pt, const double d, const double time,
                                                     const vector3f& normal, mixed_kdtree_primitive* mixed,
                                                     const char* data, const std::size_t dataSize ) {
    if( d < m_maxDistance ) {
        m_nearestPoint.reset( pt, d, time, normal, mixed, data, dataSize );
        m_maxDistance = d;
        m_gotHit = true;
    }
}

///

// find up to a specified number of primitives within a specified search distance

bool within_range_mixed_kdtree_distance_observer::is_collection_full( void ) {
    return m_nearestPoints.size() >= m_maxCount;
}
double within_range_mixed_kdtree_distance_observer::get_greatest_distance( void ) {
    if( m_nearestPoints.size() ) {
        return m_nearestPoints.back().get_distance();
    }
    return 0;
}
void within_range_mixed_kdtree_distance_observer::bubble_up( std::size_t i ) {
    while( i != 0 && m_nearestPoints[i] < m_nearestPoints[i - 1] ) {
        std::swap( m_nearestPoints[i - 1], m_nearestPoints[i] );
        --i;
    }
}
within_range_mixed_kdtree_distance_observer::within_range_mixed_kdtree_distance_observer(
    std::vector<mixed_kdtree_nearest_point>& outResults, const float maxDistance, const std::size_t maxCount )
    : m_maxDistance( maxDistance )
    , m_maxCount( maxCount )
    , m_nearestPoints( outResults ) {
    m_nearestPoints.clear();
}
bool within_range_mixed_kdtree_distance_observer::is_done( void ) { return false; }
double within_range_mixed_kdtree_distance_observer::get_max_range( void ) {
    return is_collection_full() ? get_greatest_distance() : m_maxDistance;
}
bool within_range_mixed_kdtree_distance_observer::is_in_range( double d ) { return d < get_max_range(); }
// bool within_range_mixed_kdtree_distance_observer::wants_hit_details( double d ) { return d < get_max_range(); }
void within_range_mixed_kdtree_distance_observer::insert( const vector3f& pt, const double d, const double time,
                                                          const vector3f& normal, mixed_kdtree_primitive* primitive,
                                                          const char* data, const std::size_t dataSize ) {
    // in m_nearestPoints:
    // - primitve * are unique
    // - sorted by ascending distance
    // original implementation is in collect_nearest_faces() of trimesh3_kdtree_node.hpp
    // implementation here does not assume that the same primitive is at the same distance.
    // one motivation for this is I am not sure how the boundbox should affect the
    // search.  should the found points be constrained to the boundbox?  or are we
    // interested in points on any object, given that the object intersects with the
    // boundbox ?

    // does the distance belong in  the array ?
    if( d < get_max_range() ) {
        // This point will be inserted.
        // Is the primitive already in the list ?
        std::vector<mixed_kdtree_nearest_point>::iterator oldPos;
        oldPos =
            std::find_if( m_nearestPoints.begin(), m_nearestPoints.end(),
                          boost::bind( std::equal_to<const mixed_kdtree_primitive*>(),
                                       boost::bind( &mixed_kdtree_nearest_point::get_primitive, _1 ), primitive ) );
        if( oldPos == m_nearestPoints.end() ) {
            // The primitive is not in the list, but it should be.
            if( is_collection_full() ) {
                // Replace the last (greatest distance) primitive, and maybe
                // move it earlier in the list.
                m_nearestPoints.back().reset( pt, d, time, normal, primitive, data, dataSize );
                bubble_up( m_nearestPoints.size() - 1 );
            } else {
                // Append the primitive, and maybe move it earlier in the list.
                m_nearestPoints.push_back(
                    mixed_kdtree_nearest_point( pt, d, time, normal, primitive, data, dataSize ) );
                bubble_up( m_nearestPoints.size() - 1 );
            }
        } else {
            // The primitive is already in the list.
            // If its distance is reduced, then update its distance, and
            // maybe move it earlier in the list.
            if( ( *oldPos ).get_distance() > d ) {
                ( *oldPos ).reset( pt, d, time, normal, primitive, data, dataSize );
                bubble_up( std::distance( m_nearestPoints.begin(), oldPos ) );
            }
        }
    }
}

//
// ray traversals
//

find_any_intersection_mixed_kdtree_ray_observer::find_any_intersection_mixed_kdtree_ray_observer()
    : m_gotHit( false ) {}

bool find_any_intersection_mixed_kdtree_ray_observer::is_done( void ) { return m_gotHit; }

double find_any_intersection_mixed_kdtree_ray_observer::get_max_range( void ) {
    return std::numeric_limits<double>::max();
}

bool find_any_intersection_mixed_kdtree_ray_observer::is_in_range( double /*t*/ ) { return true; }

// bool find_any_intersection_mixed_kdtree_ray_observer::wants_hit_details( double t ) {
//    if( t < std::numeric_limits<double>::max() ) {
//        m_gotHit = true;
//    }
//    return false;
//}
bool find_any_intersection_mixed_kdtree_ray_observer::got_hit( void ) const { return m_gotHit; }

// void find_any_intersection_mixed_kdtree_ray_observer::insert( const mixed_kdtree_ray_intersection & intersection ) {
//    if( intersection.get_t() < std::numeric_limits<double>::max() ) {
//        m_gotHit = true;
//    }
//}

void find_any_intersection_mixed_kdtree_ray_observer::insert( const ray3f&, const double t, const double,
                                                              const vector3f&, mixed_kdtree_primitive*, const char*,
                                                              const std::size_t ) {
    if( t < std::numeric_limits<double>::max() ) {
        m_gotHit = true;
    }
}

///

find_nearest_intersection_mixed_kdtree_ray_observer::find_nearest_intersection_mixed_kdtree_ray_observer(
    mixed_kdtree_ray_intersection& intersection )
    : m_intersection( intersection )
    , m_gotHit( false ) {}
bool find_nearest_intersection_mixed_kdtree_ray_observer::is_done( void ) { return false; }
double find_nearest_intersection_mixed_kdtree_ray_observer::get_max_range( void ) { return m_intersection.get_t(); }
bool find_nearest_intersection_mixed_kdtree_ray_observer::is_in_range( double t ) {
    if( t < m_intersection.get_t() )
        return true;
    return false;
}
// bool find_nearest_intersection_mixed_kdtree_ray_observer::wants_hit_details( double t ) {
//    return is_in_range( t );
//}
// void find_nearest_intersection_mixed_kdtree_ray_observer::insert( const mixed_kdtree_ray_intersection & intersection
// ) {
//    if( is_in_range( intersection.get_t() ) ) {
//        m_intersection.reset( intersection );
//        m_gotHit = true;
//    }
//}
void find_nearest_intersection_mixed_kdtree_ray_observer::insert( const ray3f& ray, const double t, const double time,
                                                                  const vector3f& normal, mixed_kdtree_primitive* mixed,
                                                                  const char* data, const std::size_t dataSize ) {
    if( is_in_range( t ) ) {
        m_intersection.reset( ray, t, time, normal, mixed, data, dataSize );
        m_gotHit = true;
    }
}
bool find_nearest_intersection_mixed_kdtree_ray_observer::got_hit( void ) const { return m_gotHit; }

///

find_all_intersections_mixed_kdtree_ray_observer::find_all_intersections_mixed_kdtree_ray_observer(
    std::vector<mixed_kdtree_ray_intersection>& intersections )
    : m_intersections( intersections ) {}
bool find_all_intersections_mixed_kdtree_ray_observer::is_done( void ) { return false; }
bool find_all_intersections_mixed_kdtree_ray_observer::is_in_range( double /*t*/ ) { return true; }
double find_all_intersections_mixed_kdtree_ray_observer::get_max_range( void ) {
    return std::numeric_limits<double>::max();
}
// bool find_all_intersections_mixed_kdtree_ray_observer::wants_hit_details( double /*t*/ ) {
//    return true;
//}
// void find_all_intersections_mixed_kdtree_ray_observer::insert( const mixed_kdtree_ray_intersection & intersection ) {
//    m_intersections.push_back( intersection );
//}
void find_all_intersections_mixed_kdtree_ray_observer::insert( const ray3f& ray, const double t, const double time,
                                                               const vector3f& normal, mixed_kdtree_primitive* mixed,
                                                               const char* data, const std::size_t dataSize ) {
    m_intersections.push_back( mixed_kdtree_ray_intersection( ray, t, time, normal, mixed, data, dataSize ) );
}

//
// mixed_kdtree_primitive_recorder
//

mixed_kdtree_primitive_recorder::mixed_kdtree_primitive_recorder(
    std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& primitives,
    std::vector<output_channel_map_listener*>& outputChannelMapListeners )
    : m_primitives( primitives )
    , m_outputChannelMapListeners( outputChannelMapListeners )
    , m_maxDataSize( 0 ) {}

void mixed_kdtree_primitive_recorder::reserve( const std::size_t primitiveCount ) {
    m_primitives.reserve( primitiveCount );
}

void mixed_kdtree_primitive_recorder::insert( boost::shared_ptr<mixed_kdtree_primitive> primitive ) {
    m_primitives.push_back( primitive );
    m_bounds += primitive->get_bounds();
    m_maxDataSize = std::max( m_maxDataSize, primitive->get_data_size() );
}

void mixed_kdtree_primitive_recorder::insert_output_channel_map_listener(
    output_channel_map_listener* outputChannelMapListener ) {
    m_outputChannelMapListeners.push_back( outputChannelMapListener );
}

/*
std::vector<mixed_kdtree_primitive *> & mixed_kdtree_primitive_recorder::get_primitives_ref( void ) {
    return m_primitives;
}
*/

frantic::graphics::boundbox3f mixed_kdtree_primitive_recorder::get_bounds( void ) const { return m_bounds; }

std::size_t mixed_kdtree_primitive_recorder::get_data_size( void ) const { return m_maxDataSize; }

//
// mixed_kdtree_point_data
//

void mixed_kdtree_point_data::set_my_data_from( const char* data, const std::size_t dataSize ) {
    m_primitiveDataSize = dataSize;
    if( dataSize <= m_staticData.size() ) {
        m_useStaticData = true;
    } else {
        m_useStaticData = false;
        m_dynamicData.resize( dataSize );
    }
    if( dataSize ) {
        if( data == 0 ) {
            throw std::runtime_error( "mixed_kdtree_point_data Error: dataSize > 0, but data is NULL" );
        }
        if( m_useStaticData ) {
            std::copy( data, data + dataSize, m_staticData.begin() );
        } else {
            m_dynamicData.insert( m_dynamicData.begin(), data, data + dataSize );
        }
    }
}

mixed_kdtree_point_data::mixed_kdtree_point_data()
    : m_position( 0 )
    , m_time( 0 )
    , m_normal( vector3f::from_zaxis() )
    , m_primitive( 0 )
    , m_useStaticData( true )
    , m_primitiveDataSize( 0 ) {}

mixed_kdtree_point_data::mixed_kdtree_point_data( const vector3f& position, const double time, const vector3f& normal,
                                                  const mixed_kdtree_primitive* primitive )
    : m_position( position )
    , m_time( time )
    , m_normal( normal )
    , m_primitive( primitive )
    , m_useStaticData( true )
    , m_primitiveDataSize( 0 ) {}

mixed_kdtree_point_data::mixed_kdtree_point_data( const vector3f& position, const double time, const vector3f& normal,
                                                  const mixed_kdtree_primitive* primitive, const char* data,
                                                  const std::size_t dataSize )
    : m_position( position )
    , m_time( time )
    , m_normal( normal )
    , m_primitive( primitive )
    , m_useStaticData( true )
    , m_primitiveDataSize( 0 ) {
    set_my_data_from( data, dataSize );
}

void mixed_kdtree_point_data::reset() {
    m_position = vector3f( 0 );
    m_time = 0;
    m_normal = vector3f::from_zaxis();
    m_primitive = 0;

    m_primitiveDataSize = 0;
    m_useStaticData = true;
}

void mixed_kdtree_point_data::reset( const vector3f& position, const double time, const vector3f& normal,
                                     const mixed_kdtree_primitive* primitive, const char* data,
                                     const std::size_t dataSize ) {
    m_position = position;
    m_time = time;
    m_normal = normal;
    m_primitive = primitive;

    m_primitiveDataSize = dataSize;

    set_my_data_from( data, dataSize );
}

void mixed_kdtree_point_data::reset( const vector3f& position, const double time, const vector3f& normal,
                                     const mixed_kdtree_primitive* primitive ) {
    reset( position, time, normal, primitive, 0, 0 );
}

const mixed_kdtree_primitive* mixed_kdtree_point_data::get_primitive( void ) const { return m_primitive; }

vector3f mixed_kdtree_point_data::get_position( void ) const { return m_position; }

double mixed_kdtree_point_data::get_time( void ) const { return m_time; }

void mixed_kdtree_point_data::set_time( double time ) { m_time = time; }

vector3f mixed_kdtree_point_data::get_normal( void ) const { return m_normal; }

const char* mixed_kdtree_point_data::get_primitive_data( void ) const {
    return m_useStaticData ? ( &m_staticData[0] ) : ( &m_dynamicData[0] );
}

std::size_t mixed_kdtree_point_data::get_primitive_data_size( void ) const { return m_primitiveDataSize; }

void mixed_kdtree_point_data::to_raytrace_intersection( raytrace_intersection& raytraceIntersection ) const {
    raytraceIntersection.position = m_position;
    raytraceIntersection.geometricNormal = m_normal;
    // this probably isn't necessary
    // it would let us get things like the face number and
    // barycentric coordinates
    if( m_primitive ) {
        m_primitive->to_raytrace_intersection( raytraceIntersection, get_primitive_data() );
    }
}

void mixed_kdtree_point_data::to_nearest_point_search_result( nearest_point_search_result& searchResult ) const {
    searchResult.position = m_position;
    searchResult.geometricNormal = m_normal;
    // this probably isn't necessary
    if( m_primitive ) {
        m_primitive->to_nearest_point_search_result( searchResult, get_primitive_data() );
    }
}

//
// mixed_kdtree_nearest_point_data
//

mixed_kdtree_nearest_point_data::mixed_kdtree_nearest_point_data()
    : m_distance( std::numeric_limits<double>::max() ) {}

mixed_kdtree_nearest_point_data::mixed_kdtree_nearest_point_data( const double distance )
    : m_distance( distance ) {}

void mixed_kdtree_nearest_point_data::reset( const double distance ) { m_distance = distance; }

double mixed_kdtree_nearest_point_data::get_distance( void ) const { return m_distance; }

bool mixed_kdtree_nearest_point_data::operator<( const mixed_kdtree_nearest_point_data& rhs ) const {
    return m_distance < rhs.m_distance;
}

void mixed_kdtree_nearest_point_data::to_nearest_point_search_result(
    nearest_point_search_result& searchResult ) const {
    searchResult.distance = static_cast<float>( m_distance );
}

//
// mixed_kdtree_nearest_point
//

mixed_kdtree_nearest_point::mixed_kdtree_nearest_point() {}

mixed_kdtree_nearest_point::mixed_kdtree_nearest_point( const vector3f& pt, const double distance, const double time,
                                                        const vector3f& normal, mixed_kdtree_primitive* mixed )
    : mixed_kdtree_nearest_point_data( distance )
    , mixed_kdtree_point_data( pt, time, normal, mixed ) {}

mixed_kdtree_nearest_point::mixed_kdtree_nearest_point( const vector3f& pt, const double distance, const double time,
                                                        const vector3f& normal, mixed_kdtree_primitive* mixed,
                                                        const char* data, const std::size_t dataSize )
    : mixed_kdtree_nearest_point_data( distance )
    , mixed_kdtree_point_data( pt, time, normal, mixed, data, dataSize ) {}

void mixed_kdtree_nearest_point::reset( const vector3f& pt, const double distance, const double time,
                                        const vector3f& normal, mixed_kdtree_primitive* mixed ) {
    mixed_kdtree_nearest_point_data::reset( distance );
    mixed_kdtree_point_data::reset( pt, time, normal, mixed );
}

void mixed_kdtree_nearest_point::reset( const vector3f& pt, const double distance, const double time,
                                        const vector3f& normal, mixed_kdtree_primitive* mixed, const char* data,
                                        const std::size_t dataSize ) {
    mixed_kdtree_nearest_point_data::reset( distance );
    mixed_kdtree_point_data::reset( pt, time, normal, mixed, data, dataSize );
}

void mixed_kdtree_nearest_point::to_nearest_point_search_result( nearest_point_search_result& searchResult ) const {
    mixed_kdtree_nearest_point_data::to_nearest_point_search_result( searchResult );
    mixed_kdtree_point_data::to_nearest_point_search_result( searchResult );
}

//
// mixed_kdtree_ray_intersection_data
//

mixed_kdtree_ray_intersection_data::mixed_kdtree_ray_intersection_data()
    : m_t( std::numeric_limits<double>::max() )
    , m_distance( std::numeric_limits<float>::max() ) //,
{}

mixed_kdtree_ray_intersection_data::mixed_kdtree_ray_intersection_data( const ray3f& ray, const double t )
    : m_t( t )
    , m_ray( ray )
    , m_distance( ray.direction().get_magnitude() * t ) {}

void mixed_kdtree_ray_intersection_data::reset( const ray3f& ray, const double t ) {
    m_t = t;
    m_ray = ray;
    m_distance = m_ray.direction().get_magnitude() * m_t;
}

ray3f mixed_kdtree_ray_intersection_data::get_ray( void ) const { return m_ray; }

double mixed_kdtree_ray_intersection_data::get_t( void ) const { return m_t; }

double mixed_kdtree_ray_intersection_data::get_distance( void ) const {
    return m_distance;
    // return m_ray.direction().get_magnitude() * m_t;
}

bool mixed_kdtree_ray_intersection_data::operator<( const mixed_kdtree_ray_intersection_data& rhs ) const {
    return get_distance() < rhs.get_distance();
}

void mixed_kdtree_ray_intersection_data::to_raytrace_intersection( raytrace_intersection& raytraceIntersection ) const {
    raytraceIntersection.ray = m_ray;
    raytraceIntersection.distance = m_t;
}

//
// mixed_kdtree_ray_intersection
//

mixed_kdtree_ray_intersection::mixed_kdtree_ray_intersection() {}

mixed_kdtree_ray_intersection::mixed_kdtree_ray_intersection( const ray3f& ray, const double t, const double time,
                                                              const vector3f& normal, mixed_kdtree_primitive* mixed )
    : mixed_kdtree_ray_intersection_data( ray, t )
    , mixed_kdtree_point_data( ray.at( t ), time, normal, mixed ) {}

mixed_kdtree_ray_intersection::mixed_kdtree_ray_intersection( const ray3f& ray, const double t, const double time,
                                                              const vector3f& normal, mixed_kdtree_primitive* mixed,
                                                              const char* data, const std::size_t dataSize )
    : mixed_kdtree_ray_intersection_data( ray, t )
    , mixed_kdtree_point_data( ray.at( t ), time, normal, mixed, data, dataSize ) {}

void mixed_kdtree_ray_intersection::reset( const mixed_kdtree_ray_intersection& intersection ) {
    mixed_kdtree_ray_intersection_data::reset( intersection.get_ray(), intersection.get_t() );
    mixed_kdtree_point_data::reset( intersection.get_position(), intersection.get_time(), intersection.get_normal(),
                                    intersection.get_primitive(), intersection.get_primitive_data(),
                                    intersection.get_primitive_data_size() );
}

void mixed_kdtree_ray_intersection::reset( const ray3f& ray, const double t, const double time, const vector3f& normal,
                                           mixed_kdtree_primitive* mixed ) {
    mixed_kdtree_ray_intersection_data::reset( ray, t );
    mixed_kdtree_point_data::reset( ray.at( t ), time, normal, mixed );
}

void mixed_kdtree_ray_intersection::reset( const ray3f& ray, const double t, const double time, const vector3f& normal,
                                           mixed_kdtree_primitive* mixed, const char* data,
                                           const std::size_t dataSize ) {
    mixed_kdtree_ray_intersection_data::reset( ray, t );
    mixed_kdtree_point_data::reset( ray.at( t ), time, normal, mixed, data, dataSize );
}

void mixed_kdtree_ray_intersection::to_raytrace_intersection( raytrace_intersection& raytraceIntersection ) const {
    mixed_kdtree_ray_intersection_data::to_raytrace_intersection( raytraceIntersection );
    mixed_kdtree_point_data::to_raytrace_intersection( raytraceIntersection );
}

namespace detail {

//
// common_named_channel_setters
//

void common_named_channel_setters::reset( void ) {
    m_tAccessorFloat.reset();
    m_tAccessorDouble.reset();
    m_distanceAccessorFloat.reset();
    m_distanceAccessorDouble.reset();
    m_positionAccessor.reset();
    m_normalAccessor.reset();
    m_timeAccessor.reset();
}

common_named_channel_setters::common_named_channel_setters() { reset(); }

common_named_channel_setters::common_named_channel_setters( channels::channel_map& cm ) { set_channel_map( cm ); }

void common_named_channel_setters::set_channel_map( const channels::channel_map& cm ) {
    reset();
    if( cm.has_channel( _T("RayT") ) ) {
        m_tAccessorFloat = cm.get_cvt_accessor<float>( _T("RayT") );
        m_tAccessorDouble = cm.get_cvt_accessor<double>( _T("RayT") );
    }
    if( cm.has_channel( _T("Distance") ) ) {
        m_distanceAccessorFloat = cm.get_cvt_accessor<float>( _T("Distance") );
        m_distanceAccessorDouble = cm.get_cvt_accessor<double>( _T("Distance") );
    }
    if( cm.has_channel( _T("Position") ) ) {
        m_positionAccessor = cm.get_cvt_accessor<vector3f>( _T("Position") );
    }
    if( cm.has_channel( _T("Normal") ) ) {
        m_normalAccessor = cm.get_cvt_accessor<vector3f>( _T("Normal") );
    }
    if( cm.has_channel( _T("Time") ) ) {
        m_timeAccessor = cm.get_cvt_accessor<double>( _T("Time") );
    }
}

void common_named_channel_setters::set_ray( char* /*data*/, const ray3f& /*ray*/ ) const {
    // should this be split into separate Origin and Direction channels or something similar ?
}

void common_named_channel_setters::set_t( char* data, float t ) const { m_tAccessorFloat.set( data, t ); }

void common_named_channel_setters::set_t( char* data, double t ) const { m_tAccessorDouble.set( data, t ); }

void common_named_channel_setters::set_distance( char* data, float distance ) const {
    m_distanceAccessorFloat.set( data, distance );
}

void common_named_channel_setters::set_distance( char* data, double distance ) const {
    m_distanceAccessorDouble.set( data, distance );
}

void common_named_channel_setters::set_position( char* data, const vector3f& position ) const {
    m_positionAccessor.set( data, position );
}

void common_named_channel_setters::set_normal( char* data, const vector3f& normal ) const {
    m_normalAccessor.set( data, normal );
}

void common_named_channel_setters::set_time( char* data, const double time ) const { m_timeAccessor.set( data, time ); }

//
// convert_and_set_channel
//

convert_and_set_channel::convert_and_set_channel()
    : m_destOffset( 0 )
    , m_arity( 0 )
    , m_convertAndSet( 0 ) {}

convert_and_set_channel::convert_and_set_channel( channels::channel& destCh, std::size_t sourceDataArity,
                                                  channels::data_type_t sourceDataType,
                                                  const frantic::tstring& channelNameForErrorMessage )
    : m_destOffset( destCh.offset() )
    , m_arity( destCh.arity() ) {
    const channels::data_type_t destDataType = destCh.data_type();

    if( sourceDataArity != m_arity )
        throw std::runtime_error(
            std::string() + "convert_and_set_channel() - Cannot convert channel: \"" +
            frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
            frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
            frantic::strings::to_string( channel_data_type_str( m_arity, destDataType ) ) + " due to differing arity" );

    m_convertAndSet = get_channel_type_convertor_function( sourceDataType, destDataType, channelNameForErrorMessage );
}

void convert_and_set_channel::reset( channels::channel& destCh, std::size_t sourceDataArity,
                                     channels::data_type_t sourceDataType,
                                     const frantic::tstring& channelNameForErrorMessage ) {
    m_destOffset = destCh.offset();
    m_arity = destCh.arity();

    const channels::data_type_t destDataType = destCh.data_type();

    if( sourceDataArity != m_arity )
        throw std::runtime_error(
            std::string() + "convert_and_set_channel() - Cannot convert channel: \"" +
            frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
            frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
            frantic::strings::to_string( channel_data_type_str( m_arity, destDataType ) ) + " due to differing arity" );

    m_convertAndSet = get_channel_type_convertor_function( sourceDataType, destDataType, channelNameForErrorMessage );
}

void convert_and_set_channel::set( char* dest, const char* value ) const {
    if( m_convertAndSet ) {
        m_convertAndSet( dest + m_destOffset, value, m_arity );
    }
}

std::size_t convert_and_set_channel::arity( void ) const { return m_arity; }

} // namespace detail

//
// mixed_kdtree_node
//

mixed_kdtree_node::mixed_kdtree_node() {
    // Initialize it indicating it's not a leaf node.
    m_axisAndLeafFlag = 0;
    m_children = 0;
}

mixed_kdtree_node::~mixed_kdtree_node() {
    if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
        if( m_primitivePointers != 0 ) {
            delete m_primitivePointers;
            m_primitivePointers = 0;
        }
    } else {
        if( m_children != 0 ) {
            delete[] m_children;
            m_children = 0;
        }
    }
}

// Constructor for making an interior node
// NOTE: Only one initialize function should be called once on a node.  No checks are done to make sure things are
// valid.
// void initialize( mixed_kdtree_node<mixed_kdtree_primitive>* children, int axis, float split )
void mixed_kdtree_node::initialize( mixed_kdtree_node* children, int axis, float split ) {
    m_split = split;
    m_children = children;
    // Set the axis and indicate it's not a leaf
    m_axisAndLeafFlag = axis | 0x00000000;
}

void mixed_kdtree_node::initialize( std::vector<mixed_kdtree_detail::index_t>& primitiveIndices,
                                    const std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& primitives ) {
    if( primitiveIndices.size() > 0 ) {
        m_primitivePointers = new std::vector<mixed_kdtree_primitive*>();
        m_primitivePointers->reserve( primitiveIndices.size() );
        BOOST_FOREACH( const mixed_kdtree_detail::index_t primitiveIndex, primitiveIndices ) {
            m_primitivePointers->push_back( primitives[primitiveIndex].get() );
        }
        std::vector<mixed_kdtree_detail::index_t> empty;
        primitiveIndices.swap( empty );
    } else {
        m_primitivePointers = 0;
    }

    // Indicate it's a leaf
    m_axisAndLeafFlag = 0x80000000;
}

mixed_kdtree_node* mixed_kdtree_node::left_child() { return m_children; }
mixed_kdtree_node* mixed_kdtree_node::right_child() { return m_children + 1; }

const mixed_kdtree_node* mixed_kdtree_node::left_child() const { return m_children; }
const mixed_kdtree_node* mixed_kdtree_node::right_child() const { return m_children + 1; }

int mixed_kdtree_node::get_node_count() const {
    if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
        return 1;
    } else {
        return 1 + left_child()->get_node_count() + right_child()->get_node_count();
    }
}

int mixed_kdtree_node::get_maximum_depth() const {
    if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
        return 0;
    } else {
        return 1 + ( std::max )( left_child()->get_maximum_depth(), right_child()->get_maximum_depth() );
    }
}

int mixed_kdtree_node::get_largest_leaf_size() const {
    if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
        if( m_primitivePointers != 0 )
            return (int)m_primitivePointers->size();
        else
            return 0;
    } else {
        return ( std::max )( left_child()->get_largest_leaf_size(), right_child()->get_largest_leaf_size() );
    }
}

std::size_t mixed_kdtree_node::get_populated_leaf_count() const {
    if( m_axisAndLeafFlag < 0 ) {
        if( m_primitivePointers != 0 ) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return left_child()->get_populated_leaf_count() + right_child()->get_populated_leaf_count();
    }
}

std::size_t mixed_kdtree_node::get_leaf_count() const {
    if( m_axisAndLeafFlag < 0 ) {
        return 1;
    } else {
        return left_child()->get_leaf_count() + right_child()->get_leaf_count();
    }
}

void mixed_kdtree_node::dump_tree( std::ostream& out, int depth ) const {
    if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
        if( m_primitivePointers != 0 ) {
            out << depth << ": child node with " << m_primitivePointers->size() << " children\n";
            // for( unsigned i = 0; i < m_primitiveIndices->size(); ++i ) {
            // vector3 face = mesh.get_face((*m_primitiveIndices)[i]);
            // out << depth << ": face " << (*m_primitiveIndices)[i] << " / " << face << " / " <<
            // mesh.get_vertex(face.x) << " " << mesh.get_vertex(face.y) << " " << mesh.get_vertex(face.z) << "\n";
            //}
        } else {
            out << depth << ": empty leaf node\n";
        }
    } else {
        char dimensions[] = { 'X', 'Y', 'Z' };
        out << depth << ": interior node " << dimensions[m_axisAndLeafFlag & 0x00000003] << " " << m_split << "\n";
        out << depth << ": left child\n";
        left_child()->dump_tree( out, depth + 1 );
        out << depth << ": right child\n";
        right_child()->dump_tree( out, depth + 1 );
    }
}

//
// build_kdtree*
//
#define USE_BUILD_KDTREE_CPP
#ifdef USE_BUILD_KDTREE_CPP
// template<class mixed_kdtree_primitive>
// void build_kdtree_greedy_SAH_nlogn(mixed_kdtree_node<mixed_kdtree_primitive>& node, const
// std::vector<boost::shared_ptr<mixed_kdtree_primitive> > & primitives, const boundbox3f& bounds,
// std::vector<detail::kdtreeEvent>& events, std::vector<int>& indices, std::vector<int>& objFlags, int depth ) {
void build_kdtree_greedy_SAH_nlogn( mixed_kdtree_node& node,
                                    const std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& primitives,
                                    const boundbox3f& bounds, std::vector<mixed_kdtree_detail::kdtreeEvent>& events,
                                    std::vector<mixed_kdtree_detail::index_t>& indices,
                                    std::vector<boost::int8_t>& objFlags, const int maximumDepth, const int depth ) {
    // Could use Wald's stopping condition, instead of stopping at 8
    // triangles.
    //
    // Wald's stopping condition is part of the initial bestCost below:
    // if we don't find anything better, then the recursion stops.
    // from Wald:
    // Terminate( triangles T, voxel V) =
    // {
    //   true if argmin_p C_V(p) > K_I|T|
    //   false otherwise
    // }
    // here we have T : indices( need anything else? ), V: bounds
    // that is: true if (bestCost > COST_INTERSECT * indices.size())
    // or set initial bestCost = COST_INTERSECT * indices.size()

    // if(indices.size() <= 8){
    // node.initialize(indices);
    // return;
    //}
    bool doneSplitting = false;

    const float voxelSA = bounds.get_surface_area();

    // TODO: I think this used to be done indirectly by the
    // event detection code, which counted an AABB as parallel
    // once it got smaller than a threshold.
    //
    // Since this is gone now, I'm instead checking whether
    // the clipped bounds are greater than some fraction
    // of their original extents.  Maybe this should be
    // depend on the average size of features in the scene,
    // or a function of the size of features in the current
    // box ?
    // if( bounds.get_max_dimension() < minimumVoxelLength ) {

    // I've switched to limiting the depth for now

    const std::size_t nMin = 1;

    if( depth > maximumDepth || indices.size() <= nMin ) {
        doneSplitting = true;
    }

    int bestAxis = -1;
    // float bestCost = 0.9f * static_cast<float>( mixed_kdtree_detail::COST_INTERSECT ) * indices.size();
    const float intersectionCost = static_cast<float>( mixed_kdtree_detail::COST_INTERSECT ) * indices.size();
    float bestCost = intersectionCost;
    float bestSplit = std::numeric_limits<float>::max();

    // The side that triangles on the splitting plane should go to.
    // This is \hat{p}_{side} in Wald's paper.
    mixed_kdtree_detail::kdtreeSide_type bestPlanarPartitionSide = mixed_kdtree_detail::LEFT;

    std::size_t nLeft[3], nRight[3];
    nLeft[0] = nLeft[1] = nLeft[2] = 0;
    nRight[0] = nRight[1] = nRight[2] = indices.size();

    if( !doneSplitting ) {
        for( std::size_t i = 0; i < events.size(); /*do nothing*/ ) {
            int axis = events[i].axis;
            float split = events[i].location;

            std::size_t counters[] = { 0, 0, 0 };
            for( ; i < events.size() && events[i].location == split && events[i].axis == axis; ++i )
                ++counters[events[i].type];

            nRight[axis] -= ( counters[mixed_kdtree_detail::END_] + counters[mixed_kdtree_detail::PARALLEL] );

            boundbox3f left( bounds ), right( bounds );
            left.maximum()[axis] = right.minimum()[axis] = split;
            const float leftSA = left.get_surface_area();
            const float rightSA = right.get_surface_area();

            // todo: I'm not sure if this is a reasonable criteria for
            // accepting splits:  each half must have either a primitive
            // in it, or it must have some volume.
            if( ( ( nLeft[axis] + counters[mixed_kdtree_detail::PARALLEL] ) > 0 || left.size( axis ) > 0 ) &&
                ( nRight[axis] > 0 || right.size( axis ) > 0 ) ) {
                float cost = mixed_kdtree_detail::SAH_cost(
                    voxelSA, leftSA, rightSA, nLeft[axis] + counters[mixed_kdtree_detail::PARALLEL], nRight[axis] );
                if( cost < bestCost ) {
                    bestCost = cost;
                    bestSplit = split;
                    bestAxis = axis;
                    bestPlanarPartitionSide = mixed_kdtree_detail::LEFT;
                }
            }
            // If there are no parallel triangles, then there is no
            // need to check for their cost on the right-hand side.
            if( counters[mixed_kdtree_detail::PARALLEL] > 0 ) {
                if( ( nLeft[axis] > 0 || left.size( axis ) > 0 ) &&
                    ( ( nRight[axis] + counters[mixed_kdtree_detail::PARALLEL] ) > 0 || right.size( axis ) > 0 ) ) {
                    float cost = mixed_kdtree_detail::SAH_cost(
                        voxelSA, leftSA, rightSA, nLeft[axis], nRight[axis] + counters[mixed_kdtree_detail::PARALLEL] );
                    if( cost < bestCost ) {
                        bestCost = cost;
                        bestSplit = split;
                        bestAxis = axis;
                        bestPlanarPartitionSide = mixed_kdtree_detail::RIGHT;
                    }
                }
            }

            nLeft[axis] += ( counters[mixed_kdtree_detail::BEGIN] + counters[mixed_kdtree_detail::PARALLEL] );
        }
    }

    // this count may be unnecessary
    // I added it to avoid reallocations.  This seemed beneficial
    // in test cases but it may not be in general.
    std::size_t nBestLeft = 0;
    std::size_t nBestRight = 0;

    // if( bestAxis != -1 && indices.size() > 2 && ! doneSplitting ){
    if( bestAxis != -1 && !doneSplitting ) {
        // std::cout << "splitting axis " << bestAxis << " at " << bestSplit << "\n";
        // We need to classify triangles as being on either side of the split, or crossing
        for( std::vector<mixed_kdtree_detail::index_t>::const_iterator index = indices.begin(), end = indices.end();
             index != end; ++index )
            objFlags[*index] = 0;

        for( std::vector<mixed_kdtree_detail::kdtreeEvent>::const_iterator it = events.begin(), end = events.end();
             it != end; ++it ) {
            if( it->axis != bestAxis )
                continue;

            if( it->location < bestSplit ) {
                if( it->type != mixed_kdtree_detail::BEGIN ) {
                    ++nBestLeft;
                    objFlags[it->index] = -1;
                }
            } else if( it->location > bestSplit ) {
                if( it->type != mixed_kdtree_detail::END_ ) {
                    ++nBestRight;
                    objFlags[it->index] = 1;
                }
            } else {
                if( it->type == mixed_kdtree_detail::END_ ) {
                    ++nBestLeft;
                    objFlags[it->index] = -1;
                } else if( it->type == mixed_kdtree_detail::BEGIN ) {
                    ++nBestRight;
                    objFlags[it->index] = 1;
                } else if( it->type == mixed_kdtree_detail::PARALLEL && it->location == bestSplit ) {
                    if( bestPlanarPartitionSide == mixed_kdtree_detail::LEFT ) {
                        ++nBestLeft;
                        objFlags[it->index] = -1;
                    } else {
                        ++nBestRight;
                        objFlags[it->index] = +1;
                    }
                }
            }
        }

        const std::size_t nBestIntersecting = indices.size() - nBestLeft - nBestRight;

        const std::size_t reserveLeftIndices = nBestIntersecting + nBestLeft;
        const std::size_t reserveRightIndices = nBestIntersecting + nBestRight;
        const std::size_t reserveLeftEvents0 = 6 * nBestLeft;
        const std::size_t reserveRightEvents0 = 6 * nBestRight;
        const std::size_t reserveLeftEvents1 = 6 * nBestIntersecting;
        const std::size_t reserveRightEvents1 = 6 * nBestIntersecting;

        std::vector<mixed_kdtree_detail::kdtreeEvent> leftEvents[3], rightEvents[3];

        // Most events should be split into one side xor the other
        // leftEvents[0].reserve(events.size()/2);
        // rightEvents[0].reserve(events.size()/2);
        leftEvents[0].reserve( reserveLeftEvents0 );
        rightEvents[0].reserve( reserveRightEvents0 );

        for( std::vector<mixed_kdtree_detail::kdtreeEvent>::const_iterator it = events.begin(), end = events.end();
             it != end; ++it ) {
            if( objFlags[it->index] < 0 )
                leftEvents[0].push_back( *it );
            else if( objFlags[it->index] > 0 )
                rightEvents[0].push_back( *it );
        }

        frantic::clear_with_swap( events );

        std::vector<mixed_kdtree_detail::index_t> leftIndices, rightIndices;
        // reserve more than this?
        // leftIndices.reserve(indices.size()/2);
        // rightIndices.reserve(indices.size()/2);
        leftIndices.reserve( reserveLeftIndices );
        rightIndices.reserve( reserveRightIndices );

        boundbox3f leftBounds( bounds ), rightBounds( bounds );
        leftBounds.maximum()[bestAxis] = rightBounds.minimum()[bestAxis] = bestSplit;

        leftEvents[1].reserve( reserveLeftEvents1 );
        rightEvents[1].reserve( reserveRightEvents1 );

        BOOST_FOREACH( const mixed_kdtree_detail::index_t i, indices ) {
            if( objFlags[i] < 0 )
                leftIndices.push_back( i );
            else if( objFlags[i] > 0 )
                rightIndices.push_back( i );
            else {
                // const vector3& f = mesh.get_face(i);
                // const boundbox3f unclippedBounds( primitives[i]->get_bounds() );
                // const float unclippedVolume = unclippedBounds.get_volume();

                boundbox3f clippedBounds;
                clippedBounds = primitives[i]->intersect_with( leftBounds );
                // if( leftBounds.intersect_with_triangle(mesh.get_vertex(f.x), mesh.get_vertex(f.y),
                // mesh.get_vertex(f.z), clippedBounds) ){
                if( !clippedBounds.is_empty() ) {
                    // const int maxAxis = clippedBounds.get_max_dimension_axis();
                    // if ( unclippedBounds.size( maxAxis ) < 16.f * clippedBounds.size( maxAxis ) ) {
                    leftIndices.push_back( i );
                    extract_events( leftEvents[1], i, clippedBounds );
                }

                clippedBounds = primitives[i]->intersect_with( rightBounds );
                // if( rightBounds.intersect_with_triangle(mesh.get_vertex(f.x), mesh.get_vertex(f.y),
                // mesh.get_vertex(f.z), clippedBounds) ){
                if( !clippedBounds.is_empty() ) {
                    // const int maxAxis = clippedBounds.get_max_dimension_axis();
                    // if ( unclippedBounds.size( maxAxis ) < 16.f * clippedBounds.size( maxAxis ) ) {
                    rightIndices.push_back( i );
                    extract_events( rightEvents[1], i, clippedBounds );
                }
            }
        }

        // std::cout << "0: " << leftEvents[0].size() << "/" << leftEvents[0].capacity() << " " << rightEvents[0].size()
        // <<
        // "/" << rightEvents[0].capacity() << std::endl;
        // std::cout << "1: " << leftEvents[1].size() << "/" << leftEvents[1].capacity() << " " << rightEvents[1].size()
        // <<
        // "/" << rightEvents[1].capacity() << std::endl;
        // std::cout << "i: " << leftIndices.size() << "/" << leftIndices.capacity() << " " << rightIndices.size() <<
        // "/" << rightIndices.capacity() << std::endl;

        frantic::clear_with_swap( indices );

        node.initialize( new mixed_kdtree_node[2], bestAxis, bestSplit );

        leftEvents[2].resize( leftEvents[0].size() + leftEvents[1].size() );
        std::sort( leftEvents[1].begin(), leftEvents[1].end() );
        std::merge( leftEvents[0].begin(), leftEvents[0].end(), leftEvents[1].begin(), leftEvents[1].end(),
                    leftEvents[2].begin() );

        // leftEvents[2] and leftIndices
        // will be deallocated inside the recursive function calls.
        // Clear the others now.
        for( int i = 0; i < 2; ++i ) {
            frantic::clear_with_swap( leftEvents[i] );
        }

        rightEvents[2].resize( rightEvents[0].size() + rightEvents[1].size() );
        std::sort( rightEvents[1].begin(), rightEvents[1].end() );
        std::merge( rightEvents[0].begin(), rightEvents[0].end(), rightEvents[1].begin(), rightEvents[1].end(),
                    rightEvents[2].begin() );

        // rightEvents[2] and rightIndices
        // will be deallocated inside the recursive function calls.
        // Clear the others now.
        for( int i = 0; i < 2; ++i ) {
            frantic::clear_with_swap( rightEvents[i] );
        }

        build_kdtree_greedy_SAH_nlogn( *node.left_child(), primitives, leftBounds, leftEvents[2], leftIndices, objFlags,
                                       maximumDepth, depth + 1 );
        build_kdtree_greedy_SAH_nlogn( *node.right_child(), primitives, rightBounds, rightEvents[2], rightIndices,
                                       objFlags, maximumDepth, depth + 1 );
    } else {
        node.initialize( indices, primitives );
    }
}

// template<class mixed_kdtree_primitive>
// void build_kdtree_greedy_SAH_nlogn(mixed_kdtree_node<mixed_kdtree_primitive>& node,
// std::vector<boost::shared_ptr<mixed_kdtree_primitive> >& primitives, const boundbox3f& bounds) {
void build_kdtree_greedy_SAH_nlogn( mixed_kdtree_node& node,
                                    std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& primitives,
                                    const boundbox3f& bounds ) {
    const std::size_t primitiveCount = primitives.size();

    if( primitiveCount > std::numeric_limits<mixed_kdtree_detail::index_t>::max() ) {
        throw std::runtime_error(
            "build_kdtree_greedy_SAH_nlogn Internal Error: the number of primitives (" +
            boost::lexical_cast<std::string>( primitiveCount ) + ") exceeds the maximum primitive index (" +
            boost::lexical_cast<std::string>( std::numeric_limits<mixed_kdtree_detail::index_t>::max() ) + ")." );
    }

    std::vector<mixed_kdtree_detail::index_t> indices;
    std::vector<boost::int8_t> objFlags( primitiveCount, 0 );
    std::vector<mixed_kdtree_detail::kdtreeEvent> events;

    indices.reserve( primitiveCount );
    events.reserve( 6 * primitiveCount );

    for( std::size_t i = 0; i < primitiveCount; ++i ) {
        boundbox3f clippedBounds;
        clippedBounds = primitives[i]->intersect_with( bounds );

        if( !clippedBounds.is_empty() ) {
            indices.push_back( static_cast<mixed_kdtree_detail::index_t>( i ) );
            extract_events( events, static_cast<mixed_kdtree_detail::index_t>( i ), clippedBounds );
        }
    }

    // not sure if this is a good idea
    // This is an attempt to avoid degenerately small voxels.
    // Maybe I should instead adapt this to what's in the current voxel,
    // or set relative to the size of features in the scene ?
    // that I need to do this may also be a symptom of a problem with
    // the splitting planes that I accept.
    // const float minimumVoxelLength = std::min<float>( 0.00001f, 10.f * std::numeric_limits<float>::epsilon() *
    // bounds.get_max_dimension() );

    // heuristic suggested in Havran's PhD thesis, "Heuristic Ray Shooting Algorithms", 2001
    const int maximumDepth =
        static_cast<int>( 1.2 * log( static_cast<double>( primitives.size() ) ) / log( 2.0 ) + 2.0 );

    std::sort( events.begin(), events.end() );

    build_kdtree_greedy_SAH_nlogn( node, primitives, bounds, events, indices, objFlags, maximumDepth, 0 );
}
//*/
#endif
} // namespace geometry
} // namespace frantic
