// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <vector>

#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/particles/particle_grid_tree.hpp>

namespace frantic {
namespace graphics {

/**
 * This function will use dart throwing techniques to pack sample points into a poisson
 * sphere corner region with minimum spacing 'packRadius'.
 *
 * @tparam RandomGenerator A class supporting the boost::random::variate_generator<> concept.
 * @param outPoints A vector which will be filled with the generated points.
 * @param packRadius The minimum spacing between sample points. Currently it is invalid to specify
 *                    a radius > 0.1. Please see the wiki documentation (link above) for justification.
 * @param relativePacking Sets the target percentage of the optimal packing to try and obtain.
 * @param maxAttempts The maximum number of trys to associate with a region generation before failing.
 */
template <class RandomGenerator>
void build_poisson_sphere_corner_region( RandomGenerator& rng, std::vector<frantic::graphics::vector3f>& outPoints,
                                         float packRadius, float relativePacking, std::size_t maxAttempts ) {
    using frantic::graphics::plane3f;
    using frantic::graphics::vector3f;

    float r = packRadius;
    float rt2 = std::sqrt( 2.0f );
    float rt3 = std::sqrt( 3.0f );

    float v0 = r * ( 1 + 2.f * rt2 );
    float v1 = r * ( 3 + rt2 ) / rt2;
    float v2 = r * ( 1 + rt2 );

    float d1 = 1.f / rt2;
    float d2 = 1.f / rt3;

    // The edge region boundaries are covered by the seeding boundbox
    boundbox3f b( vector3f( -v0, -v0, -v0 ), vector3f( v0, v0, v0 ) );

    plane3f cornerPlanes[] = { // The face region boundaries
                               plane3f::from_normal_and_point( vector3f( d1, d1, 0 ), vector3f( v1, v1, 0 ) ),
                               plane3f::from_normal_and_point( vector3f( d1, 0, d1 ), vector3f( v1, 0, v1 ) ),
                               plane3f::from_normal_and_point( vector3f( 0, d1, d1 ), vector3f( 0, v1, v1 ) ),

                               // The center region boundary
                               plane3f::from_normal_and_point( vector3f( d2, d2, d2 ), vector3f( v2, v2, v2 ) ) };

    frantic::channels::channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.end_channel_definition();

    frantic::particles::particle_grid_tree pgt( pcm, 4 * r );

    // The volume of the region being packed by spheres is the seeding bounds + the radius on each
    // side.
    float packingVolume = 2 * ( v0 + r );
    packingVolume = packingVolume * packingVolume * packingVolume;

    std::size_t estParticles =
        static_cast<std::size_t>( ( relativePacking * packingVolume ) / ( 4.f * rt2 * ( r * r * r ) ) );

    outPoints.reserve( estParticles );
    for( std::size_t i = 0, count = 0; i < maxAttempts; ++i ) {
        // TODO: Use a more sophisicated dart-throwing concept
        vector3f p = b.random_vector( rng );
        if( !pgt.has_particle_near( p, 2 * r ) ) {
            pgt.insert( (char*)&p.x );

            // The corner region is symmetric about all axes, so only check the all positive
            // octant.
            vector3f pClip( std::abs( p.x ), std::abs( p.y ), std::abs( p.z ) );
            bool isValid = true;
            for( std::size_t i = 0; i < 4 && isValid; ++i ) {
                if( cornerPlanes[i].get_signed_distance_to_plane( pClip ) >= 0 )
                    isValid = false;
            }

            if( isValid )
                outPoints.push_back( p );

            if( ++count == estParticles )
                break;
        }
    }
}

namespace detail {
/**
 * This function will add a series of poinrs to a particle grid tree, offset with a given
 * vector. This is used to add corner, edge, face, etc. constraints when seeding new
 * samples.
 *
 * @param pgt The particle grid tree to add the points to.
 * @param pts The points to add to the grid tree.
 * @param offset The offset vector to add to each point before insertion into the tree.
 */
inline void add_constraint_points( frantic::particles::particle_grid_tree& pgt,
                                   const std::vector<frantic::graphics::vector3f>& pts, const vector3f& offset ) {
    using frantic::graphics::vector3f;

    std::vector<vector3f>::const_iterator it = pts.begin();
    std::vector<vector3f>::const_iterator itEnd = pts.end();
    for( ; it != itEnd; ++it ) {
        vector3f p = ( *it ) + offset;
        pgt.insert( (char*)&p.x );
    }
}

/**
 * This function will decode an index used in build_poisson_sphere_tiles() for encoding the colors
 * used in edges, faces, and cubes. The returned value is an orientation which is meaningless
 * for cubes, but describes the axis of edge direction or face normal.
 *
 * @param index The index to decode.
 * @param numColors The number of colors that this encoded index can represent.
 * @param numEncoded The number of elements that were encoded in this index. This value should
 *                    be 2 for edges, 4 for faces, and 8 for cubes.
 * @param outCorners An array with 'numEncoded' elements. The decoded colors will be written here.
 * @return The encoded orientation of the edge or face. This is meaningless for cubes.
 */
inline std::size_t decode_index( std::size_t index, std::size_t numColors, std::size_t numEncoded,
                                 std::size_t outCorners[] ) {
    for( std::size_t i = 0; i < numEncoded; ++i, index /= numColors )
        outCorners[i] = ( index % numColors );
    return index; // The orientation is left over.
}

/**
 * This function will encode an edge with two colored endpoints and an orientation into a
 * single index.
 *
 * @param numColors The number of distinct possible colors for each endpoint.
 * @param orientation The orientation of the edge.
 * @param color1 The color on the first end of the edge.
 * @param color2 The color on the second end of the edge.
 * @return The encoded edge index.
 */
inline std::size_t encode_edge_index( std::size_t numColors, std::size_t orientation, std::size_t color1,
                                      std::size_t color2 ) {
    return color1 + ( color2 + orientation * numColors ) * numColors;
}

/**
 * This function will encode a face with four colored corners and an orientation into a
 * single index.
 *
 * @param numColors The number of distinct possible colors for each corner.
 * @param orientation The orientation of the face.
 * @param color1 The color on the first corner of the face.
 * @param color2 The color on the second corner of the face.
 * @param color3 The color on the third corner of the face.
 * @param color4 The color on the fourth corner of the face.
 * @return The encoded face index.
 */
inline std::size_t encode_face_index( std::size_t numColors, std::size_t orientation, std::size_t color1,
                                      std::size_t color2, std::size_t color3, std::size_t color4 ) {
    return color1 + ( color2 + ( color3 + ( color4 + orientation * numColors ) * numColors ) * numColors ) * numColors;
}

/**
 * This function will encode a cube with eight colored corners into a
 * single index.
 *
 * @param numColors The number of distinct possible colors for each corner.
 * @param colors An eight element array of colors, corresponding to the eight cube corners.
 * @return The encoded face index.
 */
inline std::size_t encode_cube_index( std::size_t numColors, std::size_t colors[8] ) {
    // Horner's scheme for computing polynomials.
    std::size_t result = colors[7];
    for( std::size_t i = 6; i >= 0; --i )
        result = result * numColors + colors[i];
    return result;
}
} // namespace detail

/**
 * This function will use dart throwing techniques to pack sample points into a poisson
 * sphere edge region with minimum spacing 'packRadius'.
 *
 * @tparam RandomGenerator A class supporting the boost::random::variate_generator<> concept.
 * @param outPoints A vector which will be filled with the generated points.
 * @param corners A two-element array of vectors which contain the sample points of the corner
 *                 regions adjacent on the ends of the edge region we are seeding.
 * @param orientation The axis of orientation for this edge. 0 is the x-axis, 1 the y-axis, 2 the z-axis.
 * @param packRadius The minimum spacing between sample points. Currently it is invalid to specify
 *                    a radius > 0.1. Please see the wiki documentation (link above) for justification.
 * @param relativePacking Sets the target percentage of the optimal packing to try and obtain.
 * @param maxAttempts The maximum number of trys to associate with a region generation before failing.
 */
template <class RandomGenerator>
void build_poisson_sphere_edge_region( RandomGenerator& rng, std::vector<frantic::graphics::vector3f>& outPoints,
                                       const std::vector<frantic::graphics::vector3f>* corners[2],
                                       std::size_t orientation, float packRadius, float relativePacking,
                                       std::size_t maxAttempts ) {
    using frantic::graphics::plane3f;
    using frantic::graphics::vector3f;

    float r = packRadius;
    float rt2 = std::sqrt( 2.0f );
    float d = 1.f / rt2;

    float v0 = r * ( 1 + 2 * rt2 );
    float v1 = r * ( 1 + rt2 );
    float v2 = r * ( 1 + rt2 ) / rt2;

    transform4f tm, tmInv;
    if( orientation == 1 ) {
        //[0,1,0] [-1,0,0] [0,0,1]
        tm = transform4f( vector3f( 0, 1, 0 ), vector3f( -1, 0, 0 ), vector3f( 0, 0, 1 ) );
        tmInv = transform4f( vector3f( 0, -1, 0 ), vector3f( 1, 0, 0 ), vector3f( 0, 0, 1 ) );
    } else if( orientation == 2 ) {
        //[0,0,1] [0,1,0] [-1,0,0]
        tm = transform4f( vector3f( 0, 0, 1 ), vector3f( 0, 1, 0 ), vector3f( -1, 0, 0 ) );
        tmInv = transform4f( vector3f( 0, 0, -1 ), vector3f( 0, 1, 0 ), vector3f( 1, 0, 0 ) );
    }

    // The corner region boundary and the face region boundaries are are included in
    // the seeding boundbox
    boundbox3f b;
    b += tm.transform_no_translation( vector3f( -0.5f + v0, -v1, -v1 ) );
    b += tm.transform_no_translation( vector3f( 0.5f - v0, v1, v1 ) );

    // The center region boundary
    vector3f p = tm.transform_no_translation( vector3f( 0, v2, v2 ) );
    vector3f n = tmInv.transpose_transform_no_translation( vector3f( 0, d, d ) );
    plane3f edgePlane = plane3f::from_normal_and_point( n, p );

    frantic::channels::channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.end_channel_definition();

    frantic::particles::particle_grid_tree pgt( pcm, 4 * r );

    // Pre-seed the particle grid tree, such that this edge has the
    // proper corner constraints.
    detail::add_constraint_points( pgt, *corners[0], -0.5f * tm.get_column( 0 ) );
    detail::add_constraint_points( pgt, *corners[1], 0.5f * tm.get_column( 0 ) );

    // The volume of the region being packed by spheres is the seeding bounds + the radius on each
    // side.
    float packLength = 2.f * ( 0.5f - v0 + r );
    float packWidth = 2.f * ( v1 + r );
    float packingVolume = packLength * packWidth * packWidth;

    std::size_t estParticles =
        static_cast<std::size_t>( ( relativePacking * packingVolume ) / ( 4.f * rt2 * ( r * r * r ) ) );

    outPoints.reserve( estParticles );
    for( std::size_t i = 0, count = 0; i < maxAttempts; ++i ) {
        // TODO: Use a more sophisicated dart-throwing concept
        vector3f p = b.random_vector( rng );
        if( !pgt.has_particle_near( p, 2 * r ) ) {
            pgt.insert( (char*)&p.x );

            // The region is symmetric about all axes, so only check the all positive
            // octant.
            vector3f pClip( std::abs( p.x ), std::abs( p.y ), std::abs( p.z ) );
            if( edgePlane.get_signed_distance_to_plane( pClip ) < 0 )
                outPoints.push_back( p );

            if( ++count == estParticles )
                break;
        }
    }
}

/**
 * This function will use dart throwing techniques to pack sample points into a poisson
 * sphere face region with minimum spacing 'packRadius'.
 *
 * @tparam RandomGenerator A class supporting the boost::random::variate_generator<> concept.
 * @param outPoints A vector which will be filled with the generated points.
 * @param corners A four-element array of vectors which contain the sample points of the corner
 *                 regions adjacent on the corners of the face region we are seeding.They should
 *                 be sorted such that lower indexed corners come first. Please see above for corner
 *                 ordering.
 * @param edges A four-element array of vectors which contain the sample points of the edge
 *                 regions adjacent on the sides of the face region we are seeding. They should
 *                 be sorted such that lower indexed edges come first. Please see above for edge
 *                 ordering.
 * @param orientation The axis of orientation for this edge. 0 is the x-axis, 1 the y-axis, 2 the z-axis.
 * @param packRadius The minimum spacing between sample points. Currently it is invalid to specify
 *                    a radius > 0.1. Please see the wiki documentation (link above) for justification.
 * @param relativePacking Sets the target percentage of the optimal packing to try and obtain.
 * @param maxAttempts The maximum number of trys to associate with a region generation before failing.
 */
template <class RandomGenerator>
void build_poisson_sphere_face_region( RandomGenerator& rng, std::vector<frantic::graphics::vector3f>& outPoints,
                                       const std::vector<frantic::graphics::vector3f>* corners[4],
                                       const std::vector<frantic::graphics::vector3f>* edges[4],
                                       std::size_t orientation, float packRadius, float relativePacking,
                                       std::size_t maxAttempts ) {
    using frantic::graphics::plane3f;
    using frantic::graphics::vector3f;

    float r = packRadius;
    float rt2 = std::sqrt( 2.0f );
    float d = 1.f / rt2;

    float v0 = r * ( 1 + rt2 );
    float v1 = r * ( 3 + rt2 ) / rt2;

    // Default is the y-oriented faces
    transform4f tm, tmInv;
    if( orientation == 0 ) {
        //[0,1,0] [-1,0,0] [0,0,1]
        tm = transform4f( vector3f( 0, 1, 0 ), vector3f( -1, 0, 0 ), vector3f( 0, 0, 1 ) );
        tmInv = transform4f( vector3f( 0, -1, 0 ), vector3f( 1, 0, 0 ), vector3f( 0, 0, 1 ) );
    } else if( orientation == 2 ) {
        //[1,0,0] [0,0,-1] [0,1,0]
        tm = transform4f( vector3f( 1, 0, 0 ), vector3f( 0, 0, -1 ), vector3f( 0, 1, 0 ) );
        tmInv = transform4f( vector3f( 1, 0, 0 ), vector3f( 0, 0, 1 ), vector3f( 0, -1, 0 ) );
    }

    // The center region boundary and the edge region boundaries are handled by
    // the seeding bounds
    boundbox3f b;
    b += tm.transform_no_translation( vector3f( -0.5f + v0, -r, -0.5f + v0 ) );
    b += tm.transform_no_translation( vector3f( 0.5f - v0, r, 0.5f - v0 ) );

    // The corner region boundary
    vector3f p = tm.transform_no_translation( vector3f( 0.5f - v1, 0, 0.5f - v1 ) );
    vector3f n = tmInv.transpose_transform_no_translation( vector3f( d, 0, d ) );
    plane3f facePlane = plane3f::from_normal_and_point( n, p );

    frantic::channels::channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.end_channel_definition();

    frantic::particles::particle_grid_tree pgt( pcm, 4 * r );

    // Pre-seed the particle grid tree, such that this face has the
    // proper corner and edge constraints.
    detail::add_constraint_points( pgt, *corners[0], -0.5f * ( tm.get_column( 0 ) + tm.get_column( 3 ) ) );
    detail::add_constraint_points( pgt, *corners[1], 0.5f * tm.get_column( 0 ) );
    detail::add_constraint_points( pgt, *corners[2], 0.5f * tm.get_column( 2 ) );
    detail::add_constraint_points( pgt, *corners[3], -0.5f * ( tm.get_column( 0 ) + tm.get_column( 3 ) ) );

    detail::add_constraint_points( pgt, *edges[0], -0.5f * tm.get_column( 3 ) );
    detail::add_constraint_points( pgt, *edges[1], 0.5f * tm.get_column( 3 ) );
    detail::add_constraint_points( pgt, *edges[2], -0.5f * tm.get_column( 0 ) );
    detail::add_constraint_points( pgt, *edges[3], 0.5f * tm.get_column( 0 ) );

    // The volume of the region being packed by spheres is the seeding bounds + the radius on each
    // side.
    float packLength = 2.f * ( 0.5f - v0 + r );
    float packDepth = 4 * r;
    float packingVolume = packLength * packLength * packDepth;

    std::size_t estParticles =
        static_cast<std::size_t>( ( relativePacking * packingVolume ) / ( 4.f * rt2 * ( r * r * r ) ) );

    outPoints.reserve( estParticles );
    for( std::size_t i = 0, count = 0; i < maxAttempts; ++i ) {
        // TODO: Use a more sophisicated dart-throwing concept
        vector3f p = b.random_vector( rng );
        if( !pgt.has_particle_near( p, 2 * r ) ) {
            pgt.insert( (char*)&p.x );

            // The region is symmetric about all axes, so only check the all positive
            // octant.
            vector3f pClip( std::abs( p.x ), std::abs( p.y ), std::abs( p.z ) );
            if( facePlane.get_signed_distance_to_plane( pClip ) < 0 )
                outPoints.push_back( p );

            if( ++count == estParticles )
                break;
        }
    }
}

/**
 * This function will use dart throwing techniques to pack sample points with a poisson
 * sphere distribution into tileable, unit cube with minimum spacing 'packRadius'.
 *
 * @tparam RandomGenerator A class supporting the boost::random::variate_generator<> concept.
 * @param outPoints A vector which will be filled with the generated points.
 * @param corners An eight-element array of vectors which contain the sample points of the corner
 *                 regions adjacent on the corners of the cube we are seeding.They should
 *                 be sorted such that lower indexed corners come first. Please see above for corner
 *                 ordering.
 * @param edges A twelve-element array of vectors which contain the sample points of the edge
 *                 regions adjacent on the edges of the cube we are seeding. They should
 *                 be sorted such that lower indexed edges come first. Please see above for edge
 *                 ordering.
 * @param faces A siz element array of vectors which contain the sample points of the face
 *                 regions adjacent on the sides of the cube we are seeding. They should
 *                 be sorted such that lower indexed faces come first. Please see above for face
 *                 ordering.
 * @param orientation The axis of orientation for this edge. 0 is the x-axis, 1 the y-axis, 2 the z-axis.
 * @param packRadius The minimum spacing between sample points. Currently it is invalid to specify
 *                    a radius > 0.1. Please see the wiki documentation (link above) for justification.
 * @param relativePacking Sets the target percentage of the optimal packing to try and obtain.
 * @param maxAttempts The maximum number of trys to associate with a region generation before failing.
 */
template <class RandomGenerator>
void build_poisson_sphere_cube( RandomGenerator& rng, std::vector<frantic::graphics::vector3f>& outPoints,
                                const std::vector<frantic::graphics::vector3f>* corners[8],
                                const std::vector<frantic::graphics::vector3f>* edges[12],
                                const std::vector<frantic::graphics::vector3f>* faces[6], float packRadius,
                                float relativePacking, std::size_t maxAttempts ) {
    using frantic::graphics::vector3f;

    float r = packRadius;
    float rt2 = std::sqrt( 2.f );

    frantic::channels::channel_map pcm;
    pcm.define_channel<vector3f>( _T("Position") );
    pcm.end_channel_definition();

    frantic::particles::particle_grid_tree pgt( pcm, 4 * r );

    // Pre-seed the particle grid tree, such that this cube has the
    // proper corner, edge, and face constraints.
    detail::add_constraint_points( pgt, *corners[0], vector3f( 0 ) );
    detail::add_constraint_points( pgt, *corners[1], vector3f( 1, 0, 0 ) );
    detail::add_constraint_points( pgt, *corners[2], vector3f( 0, 1, 0 ) );
    detail::add_constraint_points( pgt, *corners[3], vector3f( 1, 1, 0 ) );
    detail::add_constraint_points( pgt, *corners[4], vector3f( 0, 0, 1 ) );
    detail::add_constraint_points( pgt, *corners[5], vector3f( 1, 0, 1 ) );
    detail::add_constraint_points( pgt, *corners[6], vector3f( 0, 1, 1 ) );
    detail::add_constraint_points( pgt, *corners[7], vector3f( 1, 1, 1 ) );

    detail::add_constraint_points( pgt, *edges[0], vector3f( 0.5f, 0, 0 ) );
    detail::add_constraint_points( pgt, *edges[1], vector3f( 0.5f, 1, 0 ) );
    detail::add_constraint_points( pgt, *edges[2], vector3f( 0.5f, 0, 1 ) );
    detail::add_constraint_points( pgt, *edges[3], vector3f( 0.5f, 1, 1 ) );
    detail::add_constraint_points( pgt, *edges[4], vector3f( 0, 0.5f, 0 ) );
    detail::add_constraint_points( pgt, *edges[5], vector3f( 1, 0.5f, 0 ) );
    detail::add_constraint_points( pgt, *edges[6], vector3f( 0, 0.5f, 1 ) );
    detail::add_constraint_points( pgt, *edges[7], vector3f( 1, 0.5f, 1 ) );
    detail::add_constraint_points( pgt, *edges[8], vector3f( 0, 0, 0.5f ) );
    detail::add_constraint_points( pgt, *edges[9], vector3f( 1, 0, 0.5f ) );
    detail::add_constraint_points( pgt, *edges[10], vector3f( 0, 1, 0.5f ) );
    detail::add_constraint_points( pgt, *edges[11], vector3f( 1, 1, 0.5f ) );

    detail::add_constraint_points( pgt, *faces[0], vector3f( 0, 0.5f, 0.5f ) );
    detail::add_constraint_points( pgt, *faces[1], vector3f( 1, 0.5f, 0.5f ) );
    detail::add_constraint_points( pgt, *faces[2], vector3f( 0.5f, 0, 0.5f ) );
    detail::add_constraint_points( pgt, *faces[3], vector3f( 0.5f, 1, 0.5f ) );
    detail::add_constraint_points( pgt, *faces[4], vector3f( 0.5f, 0.5f, 0 ) );
    detail::add_constraint_points( pgt, *faces[5], vector3f( 0.5f, 0.5f, 1 ) );

    boundbox3f b( vector3f( 0, 0, 0 ), vector3f( 1, 1, 1 ) );

    // The volume of the region being packed by spheres is the seeding bounds + the radius on each
    // side.
    float packLength = 1.f + 2.f * r;
    float packingVolume = packLength * packLength * packLength;

    std::size_t estParticles =
        static_cast<std::size_t>( ( relativePacking * packingVolume ) / ( 4.f * rt2 * ( r * r * r ) ) );

    std::size_t count = 0;
    for( std::size_t i = 0; i < maxAttempts; ++i ) {
        vector3f p = b.random_vector( rng );
        if( !pgt.has_particle_near( p, 2 * r ) ) {
            pgt.insert( (char*)&p.x );
            if( ++count == estParticles )
                break;
        }
    }

    // We can make a conservative reservation size, because we know that at least 'count'
    // particles were placed in the center region of the cube.
    outPoints.reserve( count );

    frantic::particles::particle_grid_tree::const_iterator it = pgt.begin();
    frantic::particles::particle_grid_tree::const_iterator itEnd = pgt.end();
    for( ; it != itEnd; ++it ) {
        vector3f p = *reinterpret_cast<const vector3f*>( *it );
        if( p.x >= 0 && p.x < 1.f && p.y >= 0 && p.y < 1.f && p.z >= 0 && p.z < 1.f ) {
            outPoints.push_back( p );
        }
    }
}

/**
 * This function will generate 'numColors'^8 tiles of sample points which are tileable, according
 * to the colors of the 8 corners. Each tile has points generated with a poisson sphere distribution,
 * and a minimum spacing of 'packRadius'.
 *
 * @tparam RandomGenerator A class supporting the boost::random::variate_generator<> concept.
 * @tparam ForwardIterator A class supporting the STL forward iterator concept
 * @param rng A random number generator to use for seeding sample positions.
 * @param itOutTiles A forward iterator over a collection of std::vector<vector3f> objects.
 * @param numColors The number of unique colors to support.
 * @param packRadius The minimum spacing between sample points. Currently it is invalid to specify
 *                    a radius > 0.1. Please see the wiki documentation (link above) for justification.
 * @param relativePacking Sets the target percentage of the optimal packing to try and obtain.
 * @param maxAttempts The maximum number of trys to associate with a region generation before failing.
 */
template <class RandomGenerator, class ForwardIterator>
void build_poisson_sphere_tiles( RandomGenerator& rng, ForwardIterator itOutTiles, std::size_t numColors,
                                 float packRadius, float relativePacking = 0.5f, std::size_t maxAttempts = 100000 ) {
    using frantic::graphics::vector3f;

    // There are two corner colors per edge, and 3 orientations.
    // There are four corner colors per face, and 3 orientations.
    std::size_t numCorners = numColors;
    std::size_t numEdges = 3 * numColors * numColors;
    std::size_t numFaces = 3 * numColors * numColors * numColors * numColors;

    boost::scoped_array<std::vector<vector3f>> pCorners( new std::vector<vector3f>[numCorners] );
    boost::scoped_array<std::vector<vector3f>> pEdges( new std::vector<vector3f>[numEdges] );
    boost::scoped_array<std::vector<vector3f>> pFaces( new std::vector<vector3f>[numFaces] );

    frantic::diagnostics::profiling_section psCorners( _T("Corner Regions") ), psEdges( _T("Edge Regions") ),
        psFaces( _T("Face Regions") ), psCubes( _T("Entire Cube") ), psTotal( _T("Total") );

    psTotal.enter();

    // Build the corner regions
    psCorners.enter();
    for( std::size_t i = 0; i < numCorners; ++i )
        build_poisson_sphere_corner_region( rng, pCorners[i], packRadius, relativePacking, maxAttempts );
    psCorners.exit();

    // Build the edge regions
    psEdges.enter();
    for( std::size_t i = 0; i < numEdges; ++i ) {
        std::size_t colors[2];
        std::size_t orientation = detail::decode_index( i, numColors, 2, colors );

        const std::vector<vector3f>* curCorners[] = { pCorners.get() + colors[0], pCorners.get() + colors[1] };

        build_poisson_sphere_edge_region( rng, pEdges[i], curCorners, orientation, packRadius, relativePacking,
                                          maxAttempts );
    }
    psEdges.exit();

    // Build the face regions
    psFaces.enter();
    for( std::size_t i = 0; i < numFaces; ++i ) {
        std::size_t colors[4];
        std::size_t orientation = detail::decode_index( i, numColors, 4, colors );

        const std::vector<vector3f>* curCorners[] = { pCorners.get() + colors[0], pCorners.get() + colors[1],
                                                      pCorners.get() + colors[2], pCorners.get() + colors[3] };

        std::size_t o1 = ( orientation == 0 ) ? 1 : 0;
        std::size_t o2 = ( orientation == 2 ) ? 1 : 2;
        const std::vector<vector3f>* curEdges[] = {
            pEdges.get() + detail::encode_edge_index( numColors, o1, colors[0], colors[1] ),
            pEdges.get() + detail::encode_edge_index( numColors, o1, colors[2], colors[3] ),
            pEdges.get() + detail::encode_edge_index( numColors, o2, colors[0], colors[2] ),
            pEdges.get() + detail::encode_edge_index( numColors, o2, colors[1], colors[3] ) };

        build_poisson_sphere_face_region( rng, pFaces[i], curCorners, curEdges, orientation, packRadius,
                                          relativePacking, maxAttempts );
    }
    psFaces.exit();

    // There are 'numColors'^8 different cubes based on every possible combination of
    // colors assigned to the eight cube corners.
    std::size_t numCubes;
    numCubes = numColors * numColors;
    numCubes = numCubes * numCubes;
    numCubes = numCubes * numCubes;

    // Build each cube combination
    psCubes.enter();
    for( std::size_t i = 0; i < numCubes; ++i, ++itOutTiles ) {
        std::size_t colors[8];
        detail::decode_index( i, numColors, 8, colors );

        const std::vector<vector3f>* curCorners[] = { pCorners.get() + colors[0], pCorners.get() + colors[1],
                                                      pCorners.get() + colors[2], pCorners.get() + colors[3],
                                                      pCorners.get() + colors[4], pCorners.get() + colors[5],
                                                      pCorners.get() + colors[6], pCorners.get() + colors[7] };

        const std::vector<vector3f>* curEdges[] = {
            pEdges.get() + detail::encode_edge_index( numColors, 0, colors[0], colors[1] ),
            pEdges.get() + detail::encode_edge_index( numColors, 0, colors[2], colors[3] ),
            pEdges.get() + detail::encode_edge_index( numColors, 0, colors[4], colors[5] ),
            pEdges.get() + detail::encode_edge_index( numColors, 0, colors[6], colors[7] ),

            pEdges.get() + detail::encode_edge_index( numColors, 1, colors[0], colors[2] ),
            pEdges.get() + detail::encode_edge_index( numColors, 1, colors[1], colors[3] ),
            pEdges.get() + detail::encode_edge_index( numColors, 1, colors[4], colors[6] ),
            pEdges.get() + detail::encode_edge_index( numColors, 1, colors[5], colors[7] ),

            pEdges.get() + detail::encode_edge_index( numColors, 2, colors[0], colors[4] ),
            pEdges.get() + detail::encode_edge_index( numColors, 2, colors[1], colors[5] ),
            pEdges.get() + detail::encode_edge_index( numColors, 2, colors[2], colors[6] ),
            pEdges.get() + detail::encode_edge_index( numColors, 2, colors[3], colors[7] ) };

        const std::vector<vector3f>* curFaces[] = {
            pFaces.get() + detail::encode_face_index( numColors, 0, colors[0], colors[2], colors[4], colors[6] ),
            pFaces.get() + detail::encode_face_index( numColors, 0, colors[1], colors[3], colors[5], colors[7] ),

            pFaces.get() + detail::encode_face_index( numColors, 1, colors[0], colors[1], colors[4], colors[5] ),
            pFaces.get() + detail::encode_face_index( numColors, 1, colors[2], colors[3], colors[6], colors[7] ),

            pFaces.get() + detail::encode_face_index( numColors, 2, colors[0], colors[1], colors[2], colors[3] ),
            pFaces.get() + detail::encode_face_index( numColors, 2, colors[4], colors[5], colors[6], colors[7] ) };

        build_poisson_sphere_cube( rng, *itOutTiles, curCorners, curEdges, curFaces, packRadius, relativePacking,
                                   maxAttempts );
    }
    psCubes.exit();

    psTotal.exit();

    FF_LOG( stats ) << _T("Poisson Sphere sample generation timing:\n") << '\t' << psCorners << '\n'
                    << '\t' << psEdges << '\n'
                    << '\t' << psFaces << '\n'
                    << '\t' << psCubes << '\n'
                    << '\t' << psTotal << std::endl;
}

} // namespace graphics
} // namespace frantic
