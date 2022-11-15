// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/*
Equations for smoothing taken from

[1] A. Belyaev, I. Bogaevski, and Y. Ohtake. Polyhedral Surface Smoothing
with Simultaneous Mesh Regularization. Geometric Modeling and Processing,
pages 229-237, 2000.

[2] M. Desbrun, M. Meyer, P. Schr�oder, and A. H. Barr. Implicit Fairing of
Irregular Meshes using Diffusion and Curvature Flow.
Computer Graphics (SIGGRAPH 99 Proceedings), pages 317�324, 1999.

[3] A. Belyaev, Y. Ohtake, and A. Pasko. Dynamic Mesh Optimization for Polygonized
Implicit Surfaces with Sharp Features. The Visual Computer, 2003, Volume 19,
Numbers 2-3, Pages 115-126.
*/

// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/geometry/relaxation.hpp>

#include <frantic/geometry/dcel.hpp>
#include <frantic/geometry/dcel_construction.hpp>
#include <frantic/geometry/dcel_iterators.hpp>
#include <frantic/geometry/mesh_interface_utils.hpp>

#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/vector2.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#include <frantic/math/utils.hpp>

#include <boost/foreach.hpp>

namespace frantic {
namespace geometry {
namespace relaxation {

namespace detail {

/**
 * Calculate the mean curvature of the surface at point vertices[pos]. Equation 15 and 16 from [2]
 */
frantic::graphics::vector3f get_mean_curvature( std::vector<frantic::graphics::vector3>& adjFaces,
                                                std::vector<frantic::graphics::vector3f>& vertices, int pos ) {
    using frantic::graphics::vector3f;

    // formula:
    // - 1/2A sum( 1/4Ai ( (PQ dot QR)PR + ( PQ dot RQ )PQ ) )
    // Equation 15 and 16 from M. Desbrun, M. Meyer, P. Schr�oder, and A. H. Barr. Implicit fairing of irregular meshes
    // using diffusion and curvature flow.

    vector3f sum, *Q, *R;
    const vector3f& P = vertices[pos];

    float area = 0;

    for( std::size_t i = 0; i < adjFaces.size(); ++i ) {
        int q = 0, r = 0;
        if( adjFaces[i].x == pos ) {
            q = adjFaces[i].y;
            r = adjFaces[i].z;
        } else if( adjFaces[i].y == pos ) {
            q = adjFaces[i].x;
            r = adjFaces[i].z;
        } else if( adjFaces[i].z == pos ) {
            q = adjFaces[i].x;
            r = adjFaces[i].y;
        }

        Q = &vertices[q];
        R = &vertices[r];

        vector3f PQ = *Q - P, PR = *R - P, QR = *R - *Q, RQ = *Q - *R;

        float faceArea = vector3f::cross( PQ, PR ).get_magnitude() / 2;

        float PQ_dot_QR = vector3f::dot( PQ, QR );
        float PR_dot_RQ = vector3f::dot( PR, RQ );

        vector3f curvature = ( PQ_dot_QR * PR + PR_dot_RQ * PQ ) / ( 4 * faceArea );

        if( PQ_dot_QR != 0 && PR_dot_RQ != 0 )
            sum += curvature;
        area += faceArea;
    }
    return -sum / ( 2 * area );
}

/**
 * Find the average point from an original point and all the neighbouring vertices.
 */
frantic::graphics::vector3f get_umbrella_operator( const std::vector<frantic::graphics::vector3f>& vertices,
                                                   const std::vector<int>& adjacentVertices, int pos ) {
    frantic::graphics::vector3f averagePoint;

    // Fix generation of NaN vertices
    if( adjacentVertices.size() <= 0 )
        return averagePoint;

    for( std::vector<int>::const_iterator it = adjacentVertices.begin(); it != adjacentVertices.end(); ++it )
        averagePoint += vertices[*it];

    averagePoint /= (float)adjacentVertices.size();
    return averagePoint - vertices[pos];
}

/**
 * Uses the discrete mean curvature flow for smoothing and the Laplacian flow for improving mesh sampling rate.
 * Equation 11 from [1]. Gives undesired results for mesh with 'bubbles' within the surface.
 */
frantic::graphics::vector3f laplacian_mean_curvature_smooth( std::vector<frantic::graphics::vector3f>& vertices,
                                                             std::vector<int>& adjVertices,
                                                             std::vector<frantic::graphics::vector3>& adjFaces,
                                                             frantic::graphics::vector3f& normal, int pos, float C ) {
    frantic::graphics::vector3f umbrella = get_umbrella_operator( vertices, adjVertices, pos ); // laplacian
    frantic::graphics::vector3f Hn = get_mean_curvature( adjFaces, vertices, pos );             // mean curvature

    return Hn + C * ( umbrella - frantic::graphics::vector3f::dot( umbrella, normal ) * umbrella );
}

/**
 * Smoothing that examines the angle between the mean curvature vector and normal speed. Smoothing as suggested in
 * equation 13 from [1]. Epsilon should be small, .1 is good.
 * Gives undesired results for mesh with 'bubbles' within the surface.
 */
frantic::graphics::vector3f adaptive_regularizing_mean_curvature_smooth(
    std::vector<frantic::graphics::vector3f>& vertices, std::vector<int>& adjVertices,
    std::vector<frantic::graphics::vector3>& adjFaces, int pos, float elipson ) {
    frantic::graphics::vector3f F;

    // m is U0 / ||U0||
    // Hn is mean curvature

    frantic::graphics::vector3f umbrella = get_umbrella_operator( vertices, adjVertices, pos );
    frantic::graphics::vector3f m = umbrella / umbrella.get_magnitude();

    frantic::graphics::vector3f Hn = get_mean_curvature( adjFaces, vertices, pos );
    float cosTheta = frantic::graphics::vector3f::dot( m, Hn ) / Hn.get_magnitude();

    /*  (13)
      F = |H|m/cosTheta		if		cosTheta > e
      F = 2Hn - |H|m/cosTheta	if		cosTheta < -e
      F = 0					if		|cosTheta| <= e
     */

    if( cosTheta > elipson ) {
        F = Hn.get_magnitude() * m / cosTheta;
    } else if( cosTheta < -elipson ) {
        F = 2 * Hn - ( Hn.get_magnitude() * m ) / cosTheta;
    } else {
        F = vector3f( 0 );
    }

    return F;
}

/**
 * For every vertex in the mesh, find all its neighbouring faces.
 */
void find_adjacent_faces( std::vector<std::vector<frantic::graphics::vector3>>& adjFaces,
                          const frantic::geometry::trimesh3& mesh ) {
    size_t numVertices = mesh.vertex_count();

    adjFaces.clear();
    adjFaces.resize( numVertices );

    for( size_t i = 0; i < numVertices; ++i )
        adjFaces[i].reserve( 6 );

    for( std::vector<frantic::graphics::vector3>::const_iterator it = mesh.faces_ref().begin();
         it != mesh.faces_ref().end(); ++it ) {
        adjFaces[it->x].push_back( *it );
        adjFaces[it->y].push_back( *it );
        adjFaces[it->z].push_back( *it );
    }
}

/**
 * For every vertex in the mesh, find all its neighbouring faces.
 */
void find_adjacent_faces( std::vector<std::vector<int>>& adjFaces, const frantic::geometry::trimesh3& mesh ) {
    size_t numVertices = mesh.vertex_count();
    const std::vector<frantic::graphics::vector3>& faces = mesh.faces_ref();

    adjFaces.clear();
    adjFaces.resize( numVertices );

    for( size_t i = 0; i < numVertices; ++i )
        adjFaces[i].reserve( 6 );

    for( size_t i = 0; i < faces.size(); ++i ) {
        adjFaces[faces[i].x].push_back( (int)i );
        adjFaces[faces[i].y].push_back( (int)i );
        adjFaces[faces[i].z].push_back( (int)i );
    }
}

/**
 * For every vertex in the mesh, find all its neighbouring vertices.
 */
void find_adjacent_vertices( std::vector<std::vector<int>>& connections, frantic::geometry::trimesh3& mesh ) {
    size_t numVertices = mesh.vertex_count();
    connections.resize( numVertices );

    for( size_t i = 0; i < numVertices; i++ )
        connections[i].reserve( 6 );

    // for each face, add the permutations of the vertices into the vector
    std::vector<frantic::graphics::vector3>::iterator it;
    for( it = mesh.faces_ref().begin(); it != mesh.faces_ref().end(); it++ ) {
        int a = it->x;
        int b = it->y;
        int c = it->z;

        bool insertA = true, insertB = true, insertC = true;

        std::vector<int>::iterator end = connections[a].end();
        for( std::vector<int>::iterator num = connections[a].begin(); num != end; num++ ) {
            if( *num == b )
                insertB = false;
            else if( *num == c )
                insertC = false;
        }

        if( insertB )
            connections[a].push_back( b );
        if( insertC )
            connections[a].push_back( c );

        insertA = insertB = insertC = true;

        // same thing for connections[b]
        end = connections[b].end();
        for( std::vector<int>::iterator num = connections[b].begin(); num != end; num++ ) {
            if( *num == a )
                insertA = false;
            else if( *num == c )
                insertC = false;
        }

        if( insertA )
            connections[b].push_back( a );
        if( insertC )
            connections[b].push_back( c );

        insertA = insertB = insertC = true;

        // same thing for connections[c]
        end = connections[c].end();
        for( std::vector<int>::iterator num = connections[c].begin(); num != end; num++ ) {
            if( *num == a )
                insertA = false;
            else if( *num == b )
                insertB = false;
        }

        if( insertA )
            connections[c].push_back( a );
        if( insertB )
            connections[c].push_back( b );
    }
}

/**
 * The N operator fits the normals of the mesh to the normals of the implicit surface. Equation 5.5 (not labeled, in
 * between 5 and 6) in [3].
 */
frantic::graphics::vector3f n_force( std::vector<frantic::graphics::vector3f>& vertices,
                                     std::vector<frantic::graphics::vector3>& faces, std::vector<int>& adjFaces,
                                     int pos, std::vector<frantic::graphics::vector3f>& gradients ) {
    // N(P, v(P), f) = [1 / sum(A(T))] * sum( A(T)v(T) )
    // v(T) = [ PC dot m(T) ]m(T) ---> the projection of vector PC on the m(T) direction, A(T) is area of T
    // m(T) = grad f( C ) / || gradf(C) ||
    // C is centroid of T --> C = (P + P1 + P2 ) / 3
    float sumAreas = 0;
    frantic::graphics::vector3f centroid, m, PC, v, total;

    // need adjacent faces
    for( size_t f = 0; f < adjFaces.size(); ++f ) {

        frantic::graphics::vector3& face = faces[adjFaces[f]];
        frantic::graphics::vector3f& gradient = gradients[adjFaces[f]];

        float area = triangle_area( vertices[face.x], vertices[face.y], vertices[face.z] );
        sumAreas += area;

        centroid = ( vertices[face.x] + vertices[face.y] + vertices[face.z] ) / 3;

        if( gradient.get_magnitude() == 0 )
            m = gradient;
        else
            m = gradient / gradient.get_magnitude();

        PC = centroid - vertices[pos];
        v = frantic::graphics::vector3f::dot( PC, m ) * m;

        total += ( area * v );
    }

    if( sumAreas == 0 )
        return vector3f( 0 );

    return total / sumAreas;
}

/**
 * The R operator equalizes the mesh sampling rate with equation (6) from [3]. Used for the mesh optimization to an
 * implicit surface.
 */
frantic::graphics::vector3f r_force( std::vector<frantic::graphics::vector3f>& vertices, std::vector<int>& adjVertices,
                                     frantic::graphics::vector3f& normal, int pos, float C ) {
    // R(P, v(P)) = C[ U - ( U dot n ) n ]
    // n is the mesh normal at vertex P, C is a positive constant (0.1 used in paper)
    // float C = 0.1f;
    frantic::graphics::vector3f U = get_umbrella_operator( vertices, adjVertices, pos );
    return C * ( U - frantic::graphics::vector3f::dot( U, normal ) * normal );
}

/**********************************************************************************/
/*********************** Classes for the tbb parallel loops ***********************/
/**********************************************************************************/

class apply_lap_smooth {
  private:
    std::vector<frantic::graphics::vector3f>& m_currVertices;
    std::vector<frantic::graphics::vector3f>& m_prevVertices;
    std::vector<std::vector<int>>& m_adjVertices;
    const std::vector<size_t>& m_subset;
    float m_scale;

    apply_lap_smooth& operator=( const apply_lap_smooth& ); // not implemented

  public:
    apply_lap_smooth( std::vector<frantic::graphics::vector3f>& currVertices,
                      std::vector<frantic::graphics::vector3f>& prevVertices,
                      std::vector<std::vector<int>>& adjVertices, const float scale, const std::vector<size_t>& subset )
        : m_currVertices( currVertices )
        , m_prevVertices( prevVertices )
        , m_adjVertices( adjVertices )
        , m_subset( subset )
        , m_scale( scale ) {}

    size_t get_vertex( size_t opIndex ) const {
        if( m_subset.empty() ) {
            return opIndex;
        } else {
            return m_subset[opIndex];
        }
    }

    // equation 2 from [1]
    void operator()( const tbb::blocked_range<size_t>& r ) const {
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            size_t vertexId = get_vertex( i );
            m_currVertices[vertexId] =
                m_prevVertices[vertexId] +
                m_scale * get_umbrella_operator( m_prevVertices, m_adjVertices[vertexId], (int)vertexId );
        }
    }
};

class apply_lap_mean_curve_smooth {
  private:
    std::vector<frantic::graphics::vector3f>& m_currVertices;
    std::vector<frantic::graphics::vector3f>& m_prevVertices;
    std::vector<std::vector<int>>& m_adjVertices;
    std::vector<std::vector<frantic::graphics::vector3>>& m_adjFaces;
    trimesh3_vertex_channel_accessor<vector3f>& m_normals;
    float m_scale;
    float m_constant;

    apply_lap_mean_curve_smooth& operator=( const apply_lap_mean_curve_smooth& ); // not implemented

  public:
    apply_lap_mean_curve_smooth( std::vector<frantic::graphics::vector3f>& currVertices,
                                 std::vector<frantic::graphics::vector3f>& prevVertices,
                                 std::vector<std::vector<int>>& adjVertices,
                                 std::vector<std::vector<frantic::graphics::vector3>>& adjFaces,
                                 trimesh3_vertex_channel_accessor<vector3f>& normals, float scale, float constant )
        : m_currVertices( currVertices )
        , m_prevVertices( prevVertices )
        , m_adjVertices( adjVertices )
        , m_adjFaces( adjFaces )
        , m_normals( normals )
        , m_scale( scale )
        , m_constant( constant ) {}

    // equation 11 from [1]
    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            m_currVertices[i] = m_prevVertices[i] + m_scale * laplacian_mean_curvature_smooth(
                                                                  m_prevVertices, m_adjVertices[i], m_adjFaces[i],
                                                                  m_normals[i], (int)i, m_constant );
        }
    }
};

class apply_adaptive_reg_smooth {
  private:
    std::vector<frantic::graphics::vector3f>& m_currVertices;
    std::vector<frantic::graphics::vector3f>& m_prevVertices;
    std::vector<std::vector<int>>& m_adjVertices;
    std::vector<std::vector<frantic::graphics::vector3>>& m_adjFaces;
    float m_scale;
    float m_epsilon;

    apply_adaptive_reg_smooth& operator=( const apply_adaptive_reg_smooth& ); // not implemented

  public:
    apply_adaptive_reg_smooth( std::vector<frantic::graphics::vector3f>& currVertices,
                               std::vector<frantic::graphics::vector3f>& prevVertices,
                               std::vector<std::vector<int>>& adjVertices,
                               std::vector<std::vector<frantic::graphics::vector3>>& adjFaces, float scale,
                               float epsilon )
        : m_currVertices( currVertices )
        , m_prevVertices( prevVertices )
        , m_adjVertices( adjVertices )
        , m_adjFaces( adjFaces )
        , m_scale( scale )
        , m_epsilon( epsilon ) {}

    // equation 12 from [1]
    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            m_currVertices[i] =
                m_prevVertices[i] + m_scale * adaptive_regularizing_mean_curvature_smooth(
                                                  m_prevVertices, m_adjVertices[i], m_adjFaces[i], (int)i, m_epsilon );
        }
    }
};

// equation between 5 and 6 in [3]
void apply_n_force::operator()( const tbb::blocked_range<size_t>& r ) const {
    for( size_t i = r.begin(); i != r.end(); ++i ) {
        m_currVertices[i] = m_prevVertices[i] + n_force( m_prevVertices, m_faces, m_adjFaces[i], (int)i, m_gradients );
    }
}

template <class VectorType>
class apply_laplacian_smooth_mesh_interface {
    typedef VectorType vector_type;
    typedef typename VectorType::float_type float_type;

    const std::vector<vector_type>& m_inputPoints;
    mesh_channel_cvt<vector_type>& m_meshChannel;
    mesh_interface* m_mesh;
    const dcel& m_edgeStructure;

    const std::vector<size_t>& m_vertexSubset;
    float_type m_scale;

    apply_laplacian_smooth_mesh_interface& operator=( const apply_laplacian_smooth_mesh_interface& ); // not implemented
  public:
    apply_laplacian_smooth_mesh_interface( const std::vector<vector_type>& inputPoints,
                                           mesh_channel_cvt<vector_type>& channel, mesh_interface* mesh,
                                           const dcel& edgeStructure, const std::vector<size_t>& vertexSubset,
                                           float_type scale )
        : m_inputPoints( inputPoints )
        , m_meshChannel( channel )
        , m_mesh( mesh )
        , m_edgeStructure( edgeStructure )
        , m_vertexSubset( vertexSubset )
        , m_scale( scale ) {}

    size_t get_vertex( size_t index ) const {
        if( m_vertexSubset.empty() ) {
            return index;
        } else {
            return m_vertexSubset[index];
        }
    }

    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        // TODO: fix the mesh interface to allow adjacency support on arbitrary channels
        if( m_edgeStructure.vertex_count() > 0 ) {
            run_dcel( r );
        } else {
            run_mesh( r );
        }
    }

    void run_mesh( const tbb::blocked_range<std::size_t>& r ) const {
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            vertex_iterator vIt;
            const size_t vertexId = get_vertex( i );

            vector_type averagePoint;
            size_t numAdjacentVerts = 0;

            m_mesh->init_vertex_iterator( vIt, vertexId );
            do {
                // this is correct, the input points are only of the subset range
                averagePoint += m_inputPoints[m_mesh->get_edge_endpoint( vIt )];
                ++numAdjacentVerts;
            } while( m_mesh->advance_vertex_iterator( vIt ) );

            if( numAdjacentVerts > 0 ) {
                averagePoint /= (float_type)numAdjacentVerts;
                averagePoint -= m_inputPoints[i];
                m_meshChannel.set_value( vertexId, m_inputPoints[i] + m_scale * averagePoint );
            }
        }
    }

    void run_dcel( const tbb::blocked_range<std::size_t>& r ) const {
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            const size_t vertexId = get_vertex( i );

            vector_type averagePoint;
            size_t numAdjacentVerts = 0;

            BOOST_FOREACH( dcel::index_t adjVertex,
                           dcel_vertex_adjacency_cycle_range( m_edgeStructure.get_vertex_halfedge( vertexId ) ) ) {
                // this is correct, the input points are only of the subset range
                averagePoint += m_inputPoints[adjVertex];
                ++numAdjacentVerts;
            }

            if( numAdjacentVerts > 0 ) {
                averagePoint /= (float_type)numAdjacentVerts;
                averagePoint -= m_inputPoints[i];
                m_meshChannel.set_value( vertexId, m_inputPoints[i] + m_scale * averagePoint );
            }
        }
    }
};

template <class VectorType>
void run_apply_laplacian_smooth( const std::vector<VectorType>& inputPoints, mesh_channel_cvt<VectorType>& channel,
                                 mesh_interface* mesh, const dcel& edgeStructure,
                                 const std::vector<size_t>& vertexSubset, typename VectorType::float_type scale ) {
    apply_laplacian_smooth_mesh_interface<VectorType> op( inputPoints, channel, mesh, edgeStructure, vertexSubset,
                                                          scale );
    tbb::parallel_for( tbb::blocked_range<size_t>( 0, channel.get_num_elements() ), op );
}

template <class VectorType>
class copy_mesh_channel_vertices {
    typedef VectorType vector_type;
    std::vector<vector_type>& m_inputPoints;
    const mesh_channel_cvt<vector_type>& m_meshChannel;

    copy_mesh_channel_vertices& operator=( const copy_mesh_channel_vertices& );

  public:
    copy_mesh_channel_vertices( std::vector<vector_type>& inputPoints, const mesh_channel_cvt<VectorType>& meshChannel )
        : m_inputPoints( inputPoints )
        , m_meshChannel( meshChannel ) {}

    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        for( size_t i = r.begin(); i != r.end(); ++i ) {
            m_inputPoints[i] = m_meshChannel.get_value( i );
        }
    }
};

template <class VectorType>
void run_copy_mesh_channel_vertices( std::vector<VectorType>& inputPoints,
                                     const mesh_channel_cvt<VectorType>& meshChannel ) {
    copy_mesh_channel_vertices<VectorType> op( inputPoints, meshChannel );
    tbb::parallel_for( tbb::blocked_range<size_t>( 0, meshChannel.get_num_elements() ), op );
}

template <class VectorType>
void run_laplacian_smooth( frantic::geometry::mesh_interface* mesh, frantic::geometry::mesh_channel* channel,
                           std::size_t iterationCount, typename VectorType::float_type scale ) {
    std::vector<VectorType> intermediatePoints( channel->get_num_elements() );
    mesh_channel_cvt<VectorType> channelCvt( channel );

    dcel edgeStructure;

    if( channelCvt.get_name() != _T( "verts" ) && channelCvt.get_channel_type() == mesh_channel::face_vertex ) {
        mesh_channel_to_dcel( channel, edgeStructure );
    } else {
        mesh->init_adjacency();
    }

    std::vector<size_t> vertexSubset;

    for( size_t i = 0; i < iterationCount; ++i ) {
        run_copy_mesh_channel_vertices<VectorType>( intermediatePoints, channelCvt );
        run_apply_laplacian_smooth( intermediatePoints, channelCvt, mesh, edgeStructure, vertexSubset, scale );
    }
}

} // namespace detail

/**
 * Laplacian smoothing. Equation 2 from [1]. Points move towards the median point of all neighbours multiplied by a
 * weighted amount.
 * Reduces high frequency surface information and tends to flatten the surface. Works on the mesh with holes.
 */
void laplacian_smooth( frantic::geometry::trimesh3& mesh, std::size_t iterationCount, float scale ) {
    laplacian_smooth( mesh, iterationCount, scale, std::vector<size_t>() );
}

void laplacian_smooth( frantic::geometry::trimesh3& mesh, std::size_t iterationCount, float scale,
                       const std::vector<size_t>& subset ) {
    tbb::task_scheduler_init taskScheduleInit;

    // clamp weight to [0..1]
    if( scale > 1 )
        scale = 1.f;
    else if( scale < 0 )
        scale = 0.f;

    std::vector<std::vector<int>> adjVertices; // list of vertices and all their adjacent vertices
    detail::find_adjacent_vertices( adjVertices, mesh );

    std::vector<frantic::graphics::vector3f>& currVertices = mesh.vertices_ref(); // actual vector from mesh object
    std::vector<frantic::graphics::vector3f> prevVertices = currVertices;         // copy of currVertices

    size_t numVertices = subset.size() > 0 ? subset.size() : currVertices.size();

    for( size_t iteration = 0; iteration < iterationCount; ++iteration ) {
        tbb::parallel_for( tbb::blocked_range<size_t>( 0, numVertices ),
                           detail::apply_lap_smooth( currVertices, prevVertices, adjVertices, scale, subset ) );
        prevVertices = currVertices; // copy the new vertices into the prev vector
    }
}

void laplacian_smooth( frantic::geometry::mesh_interface* mesh, std::size_t iterationCount, double scale ) {
    mesh_interface_geom_channel geomChannel( mesh );
    laplacian_smooth( mesh, &geomChannel, iterationCount, scale );
}

void laplacian_smooth( frantic::geometry::mesh_interface* mesh, const frantic::geometry::mesh_channel* channel,
                       std::size_t iterationCount, double scale ) {
    tbb::task_scheduler_init taskScheduleInit;

    frantic::geometry::mesh_channel* unconst = const_cast<frantic::geometry::mesh_channel*>( channel );

    scale = frantic::math::clamp( scale, 0.0, 1.0 );

    bool accepted = false;

    switch( channel->get_data_type() ) {
    case frantic::channels::data_type_float32:
        switch( channel->get_data_arity() ) {
        case 2:
            detail::run_laplacian_smooth<frantic::graphics2d::vector2f>( mesh, unconst, iterationCount,
                                                                         float( scale ) );
            accepted = true;
            break;
        case 3:
            detail::run_laplacian_smooth<frantic::graphics::vector3f>( mesh, unconst, iterationCount, float( scale ) );
            accepted = true;
            break;
        default:
            accepted = false;
        }
        break;
    case frantic::channels::data_type_float64:
        switch( channel->get_data_arity() ) {
        case 3:
            detail::run_laplacian_smooth<frantic::graphics::vector3fd>( mesh, unconst, iterationCount, scale );
            accepted = true;
            break;
        default:
            accepted = false;
        }
        break;
    default:
        accepted = false;
    }

    if( !accepted ) {
        throw std::runtime_error( "laplacian_smooth: Error, not implemented for type: " +
                                  frantic::strings::to_string(
                                      channel_data_type_str( channel->get_data_arity(), channel->get_data_type() ) ) );
    }
}

/**
 * Uses the discrete mean curvature flow for smoothing and the Laplacian flow for improving mesh sampling rate. Equation
 * 11 from [1].
 * Scale / step should also be a small positive number.
 * Gives undesired results for mesh with 'bubbles' within the surface.
 */
void laplacian_mean_curvature_smooth( frantic::geometry::trimesh3& mesh, std::size_t iterationCount, float scale,
                                      float constant ) {
    tbb::task_scheduler_init taskScheduleInit;

    std::vector<std::vector<frantic::graphics::vector3>> adjFaces;
    std::vector<std::vector<int>> adjVertices;

    detail::find_adjacent_faces( adjFaces, mesh );
    detail::find_adjacent_vertices( adjVertices, mesh );

    mesh.build_vertex_normals();
    trimesh3_vertex_channel_accessor<vector3f> normals = mesh.get_vertex_channel_accessor<vector3f>( _T("Normal") );

    // clamp scale to [0..1]
    if( scale > 1 )
        scale = 1.f;
    else if( scale < 0 )
        scale = 0.f;

    std::vector<frantic::graphics::vector3f>& currVertices = mesh.vertices_ref();
    std::vector<frantic::graphics::vector3f> prevVertices = currVertices;

    size_t numVertices = currVertices.size();

    for( size_t iteration = 0; iteration < iterationCount; ++iteration ) {
        tbb::parallel_for( tbb::blocked_range<size_t>( 0, numVertices ),
                           detail::apply_lap_mean_curve_smooth( currVertices, prevVertices, adjVertices, adjFaces,
                                                                normals, scale, constant ) );
        prevVertices = currVertices;
    }
}

/**
 * Smoothing that examines the angle between the mean curvature vector and normal speed. Smoothing as suggested in
 * equation 13 from [1]. Epsilon should be small, .1 is good.
 * Scale / step should also be a small positive number.
 * Gives undesired results for mesh with 'bubbles' within the surface.
 */
void adaptive_regularizing_mean_curvature_smooth( frantic::geometry::trimesh3& mesh, std::size_t iterationCount,
                                                  float scale, float epsilon ) {
    tbb::task_scheduler_init taskScheduleInit;

    std::vector<std::vector<frantic::graphics::vector3>> adjFaces;
    std::vector<std::vector<int>> adjVertices;

    detail::find_adjacent_faces( adjFaces, mesh );
    detail::find_adjacent_vertices( adjVertices, mesh );

    // clamp scale to [0..1]
    if( scale > 1 )
        scale = 1.f;
    else if( scale < 0 )
        scale = 0.f;

    std::vector<frantic::graphics::vector3f>& currVertices = mesh.vertices_ref();
    std::vector<frantic::graphics::vector3f> prevVertices = currVertices;

    size_t numVertices = currVertices.size();

    for( size_t iteration = 0; iteration < iterationCount; ++iteration ) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>( 0, numVertices ),
            detail::apply_adaptive_reg_smooth( currVertices, prevVertices, adjVertices, adjFaces, scale, epsilon ) );
        prevVertices = currVertices;
    }
}
} // namespace relaxation
} // namespace geometry
} // namespace frantic
