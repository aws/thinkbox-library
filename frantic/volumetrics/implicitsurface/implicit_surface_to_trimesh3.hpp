// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
//#define FRANTIC_DISABLE_THREADS

#pragma once

#include <vector>

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics2d/boundrect2.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_grid_tree.hpp>
#include <frantic/volumetrics/implicitsurface/marching_cubes_table.hpp>
#include <frantic/volumetrics/levelset/rle_defined_and_adj_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/rle_plane.hpp>
//#include <frantic/diagnostics/profiling_manager.hpp>
#include <frantic/volumetrics/implicitsurface/level_set_implicit_surface_policies.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>

#pragma warning( push )
#pragma warning( disable : 4512 4100 4244 4245 )
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/queuing_mutex.h>
#include <tbb/spin_mutex.h>
#include <tbb/spin_rw_mutex.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_thread.h>
#pragma warning( pop )

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

// extern frantic::diagnostics::profiling_manager sparsePM;

namespace detail {

template <class T>
void cumulative_sum( std::vector<T>& v ) {
    for( std::size_t i = 1; i < v.size(); ++i ) {
        v[i] += v[i - 1];
    }
}

frantic::graphics::vector3 get_direction_neighbor_offset( int boundaryCase );

struct exposed_voxel_vertices {
    frantic::graphics::vector3 voxelCoord;
    frantic::graphics::vector3 vertices;
    exposed_voxel_vertices( const frantic::graphics::vector3& voxelCoord, const frantic::graphics::vector3& vertices )
        : voxelCoord( voxelCoord )
        , vertices( vertices ) {}
};

struct exposed_block_vertices {
    frantic::graphics::vector3
        blockCoord; // I think this should be implicit given that you are writing; TODO : remove this ?
    frantic::graphics::vector3 voxelCoord;
    frantic::graphics::vector3 vertices;
    exposed_block_vertices( const frantic::graphics::vector3& blockCoord, const frantic::graphics::vector3& voxelCoord,
                            const frantic::graphics::vector3& vertices )
        : blockCoord( blockCoord )
        , voxelCoord( voxelCoord )
        , vertices( vertices ) {}
};

inline int get_block_color( const frantic::graphics::vector3& v ) {
    return ( ( v.z % 2 ) == 0 ? 0 : 4 ) | ( ( v.y % 2 ) == 0 ? 0 : 2 ) | ( ( v.x % 2 ) == 0 ? 0 : 1 );
}

struct block_color_comparison {
    // true : a comes before b
    bool operator()( const frantic::graphics::vector3& a, const frantic::graphics::vector3& b ) {
        return get_block_color( a ) < get_block_color( b );
        const int az = ( a.z % 2 ) != 0;
        const int bz = ( b.z % 2 ) != 0;
        if( az != bz ) {
            return az < bz;
        }

        const int ay = ( a.y % 2 ) != 0;
        const int by = ( b.y % 2 ) != 0;
        if( ay != by ) {
            return ay < by;
        }

        const int ax = ( a.x % 2 ) != 0;
        const int bx = ( b.x % 2 ) != 0;
        return ax < bx;
    }
};

inline boost::int32_t get_invalid_vertex_number() { return -1; }

inline bool is_internal_vertex_number( boost::int32_t i ) { return i >= 0; }

inline bool is_external_vertex_number( boost::int32_t i ) { return i < -1; }

inline bool is_vertex_number( boost::int32_t i ) { return i != get_invalid_vertex_number(); }

inline boost::int32_t encode_external_vertex_number( boost::int32_t vertexNumber ) { return ( -vertexNumber ) - 2; }

inline boost::int32_t decode_external_vertex_number( boost::int32_t encodedVertexNumber ) {
    return ( -encodedVertexNumber ) - 2;
}

struct deferred_center_vertex {
    boost::int32_t vertexNumber;
    boost::int32_t edgeVertexNumbers[12];
    deferred_center_vertex( const boost::int32_t vertexNumber, const boost::int32_t edgeVertexNumbers[12] )
        : vertexNumber( vertexNumber ) {
        memcpy( this->edgeVertexNumbers, edgeVertexNumbers, 12 * sizeof( boost::int32_t ) );
    }
    void to_global_vertices_and_holes( boost::int32_t firstVertexNumber, boost::int32_t firstHoleNumber ) {
        vertexNumber += firstVertexNumber;
        for( int edge = 0; edge < 12; ++edge ) {
            if( is_vertex_number( edgeVertexNumbers[edge] ) ) {
                if( is_external_vertex_number( edgeVertexNumbers[edge] ) ) {
                    edgeVertexNumbers[edge] = encode_external_vertex_number(
                        decode_external_vertex_number( edgeVertexNumbers[edge] ) + firstHoleNumber );
                } else {
                    edgeVertexNumbers[edge] += firstVertexNumber;
                }
            }
        }
    }
};

typedef boost::unordered_map<frantic::graphics::vector3, frantic::graphics::vector3, voxel_coord_hasher> voxel_vertex_t;

void fill_initial_cube_case_values( frantic::graphics2d::size2 size, const std::vector<float>& voxelCornerDensities,
                                    std::vector<unsigned char>& outCubeCases );

void fill_initial_cube_case_values( frantic::graphics2d::size2 size, const std::vector<float>& voxelCornerDensities,
                                    std::vector<unsigned char>& outCubeCases,
                                    std::vector<boost::int32_t>& outNewVertexCount );

void fill_sparse_initial_vertex_indices( const boost::shared_array<float>& voxelCornerDensities,
                                         const rle_plane& currentRLP,
                                         boost::shared_array<frantic::graphics::vector3>& outVertexIndices );

void fill_cube_case_values( frantic::graphics2d::size2 size, const std::vector<unsigned char>& previousCubeCases,
                            const std::vector<float>& voxelCornerDensities, std::vector<unsigned char>& outCubeCases );

void fill_cube_case_values( const frantic::graphics::size3& size, const std::vector<float>& voxelCornerDensities,
                            std::vector<boost::uint8_t>& outCubeCases );

void fill_cube_case_values_mt( tbb::affinity_partitioner& partitioner, frantic::graphics2d::size2 size,
                               const std::vector<unsigned char>& previousCubeCases,
                               const std::vector<float>& voxelCornerDensities, std::vector<unsigned char>& outCubeCases,
                               std::vector<boost::int32_t>& outNewVertexCount );

void fill_cube_case_values_mt( tbb::affinity_partitioner& partitioner, frantic::graphics2d::size2 size,
                               const std::vector<unsigned char>& previousCubeCases,
                               const std::vector<float>& voxelCornerDensities,
                               std::vector<unsigned char>& outCubeCases );

void fill_sparse_cube_case_values( float defaultOutsideDensity,
                                   const boost::shared_array<float>& previousVoxelCornerDensities,
                                   const rle_plane& prevRLP, const boost::shared_array<float>& voxelCornerDensities,
                                   const rle_plane& currentRLP, boost::shared_array<unsigned char>& outDefinedCubeCases,
                                   rle_plane& outDefinedCubeCasesRLP,
                                   boost::shared_array<frantic::graphics::vector3>& outVertexIndices );

void generate_faces_for_plane( const frantic::volumetrics::marching_cubes_table& mct, frantic::graphics2d::size2 size,
                               const std::vector<unsigned char>& previousCubeCases,
                               const std::vector<unsigned char>& currentCubeCases,
                               const std::vector<int>& previousVertexIndices,
                               const std::vector<int>& currentVertexIndices, frantic::geometry::trimesh3& outMesh );

void generate_faces_for_sparse_plane( const marching_cubes_table& mct,
                                      const boost::shared_array<unsigned char> currentCubeCases,
                                      const frantic::volumetrics::rle_plane& currentCubeCasesRLP,
                                      const boost::shared_array<frantic::graphics::vector3>& previousVertexIndices,
                                      const boost::shared_array<frantic::graphics::vector3>& currentVertexIndices,
                                      frantic::geometry::trimesh3& outMesh );

inline boost::int32_t
flush_mesh_buffer( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& inputChannels,
                   frantic::geometry::trimesh3& buffer, tbb::spin_mutex& meshMutex, boost::int32_t& nextVertexNumber,
                   boost::int32_t& nextFaceNumber,
                   std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                   frantic::geometry::trimesh3& outMesh );

inline void flush_exposed_vertex_buffer( std::vector<exposed_voxel_vertices>& exposedVoxelVerticesBuffer,
                                         boost::int32_t firstVertexNumber, tbb::spin_mutex& mutex,
                                         std::vector<exposed_voxel_vertices>& outExposedVoxelVertices );

template <class ImplicitSurfacePolicy>
void generate_vertices_for_box(
    ImplicitSurfacePolicy& isp, const frantic::volumetrics::marching_cubes_table& mct,
    const frantic::graphics::vector3& blockCoord, const frantic::graphics::boundbox3& xyzExtents, const int blockColor,
    const std::vector<float>& voxelCornerDensities, const std::vector<boost::uint8_t>& cubeCases_,
    // const std::vector<exposed_voxel_vertices> & exposedVoxelVertices,
    // const stdext::hash_map<frantic::graphics::vector3,std::pair<size_t,size_t>,voxel_coord_hasher> &
    // blockCoordToExposedVoxelVertices,
    const voxel_vertex_t& exposedVoxelVertices, std::vector<frantic::graphics::vector3>& outVertexIndices,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
    frantic::geometry::trimesh3& outMesh, std::vector<exposed_voxel_vertices>& outExposedBoxVertices,
    typename ImplicitSurfacePolicy::vertex_workspace_t& vertexWorkspace );

template< class ImplicitSurfacePolicy >
void generate_disambiguated_faces_for_box(
	ImplicitSurfacePolicy& isp,
	const frantic::volumetrics::marching_cubes_table& mct,
	const frantic::graphics::boundbox3& xyzExtents,
	const std::vector<float>& voxelCornerDensities_,
	const std::vector<boost::uint8_t>& cubeCases_,
	const std::vector<frantic::graphics::vector3>& vertexIndices_,
	const std::vector<frantic::graphics::vector3f>& externalVertices,
	std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
	frantic::geometry::trimesh3& outMesh,
	typename ImplicitSurfacePolicy::vertex_workspace_t & vertexWorkspace/*,
	frantic::diagnostics::profiling_manager & prof*/ );

template <class ImplicitSurfacePolicy>
class MarchingCubesVertexGenerationBody {

    ImplicitSurfacePolicy& m_isp;
    const frantic::volumetrics::marching_cubes_table& m_mct;
    const frantic::graphics2d::boundrect2& m_xyExtents;
    int m_z;
    const std::vector<float>& m_previousVoxelCornerDensities;
    const std::vector<float>& m_currentVoxelCornerDensities;
    const std::vector<boost::uint8_t>& m_cubeCases;
    const std::size_t m_firstVertexNumber;
    const std::vector<boost::int32_t>& m_newVertexCount;
    std::vector<boost::int32_t>& m_outVertexIndices;
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& m_outNamedVertexChannels;
    std::vector<boost::int32_t>& m_outNewFaceCount;
    frantic::geometry::trimesh3& m_outMesh;
    // Set to write when you are reallocating memory.
    // Set to read when you are reading or writing so you want to prevent reallocations.
    // A r/w lock may no longer be useful here..
    tbb::spin_rw_mutex& m_outMeshMutex;

    MarchingCubesVertexGenerationBody& operator=( const MarchingCubesVertexGenerationBody& );

  public:
    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        ImplicitSurfacePolicy& isp = m_isp;
        const frantic::graphics2d::boundrect2& xyExtents = m_xyExtents;
        int z = m_z;
        const std::vector<float>& previousVoxelCornerDensities = m_previousVoxelCornerDensities;
        const std::vector<float>& currentVoxelCornerDensities = m_currentVoxelCornerDensities;
        // avoid secure scl checks
        // TODO: get rid of this when you're no longer using msvc8/9
        const boost::uint8_t* cubeCases = m_cubeCases.size() > 0 ? &m_cubeCases[0] : 0;
        // const std::vector<unsigned char>& cubeCases = m_cubeCases;
        std::vector<boost::int32_t>& outVertexIndices = m_outVertexIndices;
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels =
            m_outNamedVertexChannels;
        std::vector<boost::int32_t>& outNewFaceCount = m_outNewFaceCount;
        frantic::geometry::trimesh3& outMesh = m_outMesh;
        typename ImplicitSurfacePolicy::vertex_workspace_t vertexWorkspace;

        float cubeCorners[8];
        std::vector<frantic::graphics::vector3> faces( 12 );

        int yStart = static_cast<boost::int32_t>( r.begin() ) + xyExtents.minimum().y,
            yEnd = static_cast<boost::int32_t>( r.end() ) + xyExtents.minimum().y;
        int xStart = xyExtents.minimum().x, xEnd = xyExtents.maximum().x + 1;

        const int xOffset = xStart - xyExtents.minimum().x;
        const int xsize = xyExtents.xsize();

        // Loop over all the values except the initial ones.  This is because the values at the beginning of an axis
        // don't generate any vertices by construction.  Note that the places where out-of-bounds index access could
        // occur won't, because the duplication of density values for the initial row, column, and layer prevent it.
        for( int y = yStart; y < yEnd; ++y ) {
            const int rowIndex = y - xyExtents.minimum().y;
            int i = rowIndex * xyExtents.xsize() + xOffset;
            std::size_t vertexNumber = m_firstVertexNumber + ( rowIndex > 0 ? m_newVertexCount[rowIndex - 1] : 0 );
            boost::int32_t newFaceCount = 0;

            for( int x = xStart; x < xEnd; ++x ) {

                const unsigned char cubeCase = cubeCases[i];

                if( cubeCase != 0x00 && cubeCase != 0xff ) {
                    // I need to lock down the mesh as I increase the size to accomodate new verts.
                    // Once the space is allocated, no locking is needed, all writing should not
                    // conflict.
                    // tbb::spin_rw_mutex::scoped_lock meshLock( m_outMeshMutex );
                    // int vertNumber = (int)outMesh.vertex_count();

                    // Get the new vertices created by this cube case
                    frantic::graphics::vector3 newVerts;
                    marching_cubes_table::get_new_verts( cubeCase, 0, newVerts );

                    if( i >= xsize + 1 ) {
                        cubeCorners[0] = previousVoxelCornerDensities[i - xsize - 1];
                        cubeCorners[1] = previousVoxelCornerDensities[i - xsize];
                        cubeCorners[2] = previousVoxelCornerDensities[i - 1];
                        cubeCorners[3] = previousVoxelCornerDensities[i];
                        cubeCorners[4] = currentVoxelCornerDensities[i - xsize - 1];
                        cubeCorners[5] = currentVoxelCornerDensities[i - xsize];
                        cubeCorners[6] = currentVoxelCornerDensities[i - 1];
                        cubeCorners[7] = currentVoxelCornerDensities[i];
                        m_mct.get_cubecase_faces( cubeCase, cubeCorners, faces );
                        newFaceCount += static_cast<boost::int32_t>( faces.size() );
                    }

                    // meshLock.release();

                    outVertexIndices[i] = static_cast<int>( vertexNumber );

                    // Actually add the vertices to the mesh
                    if( newVerts.x >= 0 ) { // Edge 67, along the X axis
                        float density0 = currentVoxelCornerDensities[i - 1], density1 = currentVoxelCornerDensities[i];
                        frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x - 1, y, z ),
                                                   voxelCorner1 = frantic::graphics::vector3( x, y, z );
                        isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                                vertexNumber++, outNamedVertexChannels, outMesh,
                                                                vertexWorkspace );
                    }
                    if( newVerts.y >= 0 ) { // Edge 57, along the Y axis
                        float density0 = currentVoxelCornerDensities[i - xyExtents.xsize()],
                              density1 = currentVoxelCornerDensities[i];
                        frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x, y - 1, z ),
                                                   voxelCorner1 = frantic::graphics::vector3( x, y, z );
                        isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                                vertexNumber++, outNamedVertexChannels, outMesh,
                                                                vertexWorkspace );
                    }
                    if( newVerts.z >= 0 ) { // Edge 37, along the Z axis
                        float density0 = previousVoxelCornerDensities[i], density1 = currentVoxelCornerDensities[i];
                        frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x, y, z - 1 ),
                                                   voxelCorner1 = frantic::graphics::vector3( x, y, z );
                        isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                                vertexNumber++, outNamedVertexChannels, outMesh,
                                                                vertexWorkspace );
                    }
                }
                ++i;
            }
            outNewFaceCount[rowIndex] = newFaceCount;
        }
    }

    MarchingCubesVertexGenerationBody(
        ImplicitSurfacePolicy& isp, const frantic::volumetrics::marching_cubes_table& mct,
        const frantic::graphics2d::boundrect2& xyExtents, int z, const std::vector<float>& previousVoxelCornerDensities,
        const std::vector<float>& currentVoxelCornerDensities, const std::vector<boost::uint8_t>& cubeCases,
        const std::size_t firstVertexNumber, const std::vector<boost::int32_t>& newVertexCount,
        std::vector<boost::int32_t>& outVertexIndices,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        std::vector<boost::int32_t>& outNewFaceCount, frantic::geometry::trimesh3& outMesh,
        tbb::spin_rw_mutex& outMeshMutex )
        : m_isp( isp )
        , m_mct( mct )
        , m_xyExtents( xyExtents )
        , m_z( z )
        , m_previousVoxelCornerDensities( previousVoxelCornerDensities )
        , m_currentVoxelCornerDensities( currentVoxelCornerDensities )
        , m_cubeCases( cubeCases )
        , m_firstVertexNumber( firstVertexNumber )
        , m_newVertexCount( newVertexCount )
        , m_outVertexIndices( outVertexIndices )
        , m_outNamedVertexChannels( outNamedVertexChannels )
        , m_outNewFaceCount( outNewFaceCount )
        , m_outMesh( outMesh )
        , m_outMeshMutex( outMeshMutex ) {}
};

// This function adds the new vertices generated by the addition of another layer of cubes
template <class ImplicitSurfacePolicy>
void generate_vertices_for_plane_mt(
    tbb::affinity_partitioner& partitioner, ImplicitSurfacePolicy& isp,
    const frantic::volumetrics::marching_cubes_table& mct, const frantic::graphics2d::boundrect2& xyExtents, int z,
    const std::vector<float>& previousVoxelCornerDensities, const std::vector<float>& currentVoxelCornerDensities,
    const std::vector<boost::uint8_t>& cubeCases, const std::vector<boost::int32_t>& newVertexCount,
    std::vector<boost::int32_t>& outVertexIndices,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
    std::vector<boost::int32_t>& outNewFaceCount, frantic::geometry::trimesh3& outMesh ) {
    outVertexIndices.resize( xyExtents.get_area() );
    tbb::spin_rw_mutex outMeshMutex;

    const std::size_t firstVertexNumber = outMesh.vertex_count();
    if( newVertexCount.size() ) {
        const std::size_t addedVertexCount = newVertexCount.back();
        outMesh.add_vertices( addedVertexCount );
        for( std::size_t c = 0; c < outNamedVertexChannels.size(); ++c ) {
            outNamedVertexChannels[c].add_vertices( addedVertexCount );
        }
    }

    if( outNewFaceCount.size() != (size_t)xyExtents.ysize() ) {
        throw std::runtime_error( "generate_vertices_for_plane_mt Error: face count array size does not match ysize (" +
                                  boost::lexical_cast<std::string>( outNewFaceCount.size() ) +
                                  " != " + boost::lexical_cast<std::string>( xyExtents.ysize() ) + ")" );
    }

    MarchingCubesVertexGenerationBody<ImplicitSurfacePolicy> body(
        isp, mct, xyExtents, z, previousVoxelCornerDensities, currentVoxelCornerDensities, cubeCases, firstVertexNumber,
        newVertexCount, outVertexIndices, outNamedVertexChannels, outNewFaceCount, outMesh, outMeshMutex );

    tbb::blocked_range<std::size_t> range( 0, xyExtents.ysize() );
#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body, partitioner );
#else
#pragma message( "Threads are disabled" )
    body( range );
#endif
    cumulative_sum( outNewFaceCount );
}

template <class ImplicitSurfacePolicy>
struct generate_disambiguated_faces_for_plane_body {
    ImplicitSurfacePolicy& m_isp;
    const frantic::volumetrics::marching_cubes_table& m_mct;
    const frantic::graphics2d::size2& m_size;
    const std::vector<boost::uint8_t>& m_previousCubeCases;
    const std::vector<boost::uint8_t>& m_currentCubeCases;
    const std::vector<boost::int32_t>& m_previousVertexIndices;
    const std::vector<boost::int32_t>& m_currentVertexIndices;
    const std::vector<float>& m_previousVoxelCornerDensities;
    const std::vector<float>& m_voxelCornerDensities;
    const std::size_t m_firstFaceNumber;
    const std::vector<boost::int32_t> m_newFaceCount;
    tbb::spin_mutex& m_meshMutex;
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& m_outNamedVertexChannels;
    frantic::geometry::trimesh3& m_outMesh;

    generate_disambiguated_faces_for_plane_body<ImplicitSurfacePolicy>&
    operator=( const generate_disambiguated_faces_for_plane_body<ImplicitSurfacePolicy>& ); // not implemented

    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        const std::size_t yStart = r.begin();
        const std::size_t yEnd = r.end();

        const std::size_t xStart = 1;
        const std::size_t xEnd = m_size.xsize;

        std::vector<boost::int32_t> cubeVerts( 12 );
        std::vector<frantic::graphics::vector3> cubeFaces;
        float cubeCorners[8];

        // avoid secure scl checks
        // TODO: get rid of this when you're no longer using msvc8/9
        const boost::uint8_t* previousCubeCases = m_previousCubeCases.size() > 0 ? &m_previousCubeCases[0] : 0;
        const boost::uint8_t* currentCubeCases = m_currentCubeCases.size() > 0 ? &m_currentCubeCases[0] : 0;
        const boost::int32_t* previousVertexIndices =
            m_previousVertexIndices.size() > 0 ? &m_previousVertexIndices[0] : 0;
        const boost::int32_t* currentVertexIndices = m_currentVertexIndices.size() > 0 ? &m_currentVertexIndices[0] : 0;
        const float* previousVoxelCornerDensities =
            m_previousVoxelCornerDensities.size() > 0 ? &m_previousVoxelCornerDensities[0] : 0;
        const float* voxelCornerDensities = m_voxelCornerDensities.size() > 0 ? &m_voxelCornerDensities[0] : 0;

        typename ImplicitSurfacePolicy::vertex_workspace_t vertexWorkspace;

        for( std::size_t y = yStart; y != yEnd; ++y ) {
            std::size_t index = xStart + y * m_size.xsize;
            std::size_t faceNumber = ( y == 0 ? m_firstFaceNumber : ( m_firstFaceNumber + m_newFaceCount[y - 1] ) );
            const std::size_t yNewFaceCount =
                ( y == 0 ? m_newFaceCount[0] : m_newFaceCount[y] - m_newFaceCount[y - 1] );
            if( yNewFaceCount == 0 ) {
                continue;
            }
            for( std::size_t x = xStart; x != xEnd; ++x ) {
                const unsigned char cubeCase = currentCubeCases[index];
                if( cubeCase != 0x00 && cubeCase != 0xff ) {
                    // Retrieve all the vertices from the edges of this cube
                    marching_cubes_table::fill_verts(
                        cubeCase, m_currentVertexIndices[index], currentCubeCases[index - 1],
                        currentCubeCases[index - m_size.xsize], previousCubeCases[index],
                        currentCubeCases[index - m_size.xsize - 1], previousCubeCases[index - 1],
                        previousCubeCases[index - m_size.xsize], currentVertexIndices[index - 1],
                        currentVertexIndices[index - m_size.xsize], previousVertexIndices[index],
                        currentVertexIndices[index - m_size.xsize - 1], previousVertexIndices[index - 1],
                        previousVertexIndices[index - m_size.xsize], cubeVerts );

                    // retrieve the corner values for disambiguation
                    cubeCorners[0] = previousVoxelCornerDensities[index - m_size.xsize - 1];
                    cubeCorners[1] = previousVoxelCornerDensities[index - m_size.xsize];
                    cubeCorners[2] = previousVoxelCornerDensities[index - 1];
                    cubeCorners[3] = previousVoxelCornerDensities[index];
                    cubeCorners[4] = voxelCornerDensities[index - m_size.xsize - 1];
                    cubeCorners[5] = voxelCornerDensities[index - m_size.xsize];
                    cubeCorners[6] = voxelCornerDensities[index - 1];
                    cubeCorners[7] = voxelCornerDensities[index];
                    const bool newVert = m_mct.get_cubecase_faces( cubeCase, cubeCorners, cubeFaces );
                    if( newVert ) {
                        tbb::spin_mutex::scoped_lock meshLock( m_meshMutex );

                        // average the existing verts for the new vert location
                        // TODO:  a policy based vertex addition would probably be more appropriate with a
                        // solve to find a vert location
                        frantic::graphics::vector3f avg;
                        int count = 0;
                        for( size_t i = 0; i < cubeVerts.size(); ++i ) {
                            if( cubeVerts[i] >= 0 ) {
                                avg += m_outMesh.get_vertex( cubeVerts[i] );
                                ++count;
                            }
                        }
                        avg = avg * ( 1.f / count );

                        const std::size_t vertIndex = m_outMesh.vertex_count();
                        cubeVerts.push_back( static_cast<boost::int32_t>( vertIndex ) );
                        m_outMesh.add_vertex_and_resize_channels( avg );
                        if( !m_outNamedVertexChannels.empty() ) {
                            m_isp.populate_vertex_channels( m_outNamedVertexChannels, vertIndex, avg, vertexWorkspace );
                            // m_isp.set_isosurface_location_world( avg );
                            // m_isp.add_vertex_to_channels( m_outNamedVertexChannels );
                        }
                    }

                    {
                        // tbb::spin_mutex::scoped_lock meshLock( m_meshMutex );

                        for( unsigned i = 0; i < cubeFaces.size(); ++i ) {
                            const frantic::graphics::vector3 face = cubeFaces[i];
                            m_outMesh.faces_ref()[faceNumber].set( cubeVerts[face.x], cubeVerts[face.y],
                                                                   cubeVerts[face.z] );
                            ++faceNumber;
                        }
                    }
                }
                ++index;
            }
        }
    }

    generate_disambiguated_faces_for_plane_body(
        ImplicitSurfacePolicy& isp, const frantic::volumetrics::marching_cubes_table& mct,
        const frantic::graphics2d::size2& size, const std::vector<boost::uint8_t>& previousCubeCases,
        const std::vector<boost::uint8_t>& currentCubeCases, const std::vector<boost::int32_t>& previousVertexIndices,
        const std::vector<boost::int32_t>& currentVertexIndices, const std::vector<float>& previousVoxelCornerDensities,
        const std::vector<float>& voxelCornerDensities, const std::size_t firstFaceNumber,
        const std::vector<boost::int32_t>& newFaceCount, tbb::spin_mutex& meshMutex,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        frantic::geometry::trimesh3& outMesh )
        : m_isp( isp )
        , m_mct( mct )
        , m_size( size )
        , m_previousCubeCases( previousCubeCases )
        , m_currentCubeCases( currentCubeCases )
        , m_previousVertexIndices( previousVertexIndices )
        , m_currentVertexIndices( currentVertexIndices )
        , m_previousVoxelCornerDensities( previousVoxelCornerDensities )
        , m_voxelCornerDensities( voxelCornerDensities )
        , m_firstFaceNumber( firstFaceNumber )
        , m_newFaceCount( newFaceCount )
        , m_meshMutex( meshMutex )
        , m_outNamedVertexChannels( outNamedVertexChannels )
        , m_outMesh( outMesh ) {}
};

template <class ImplicitSurfacePolicy>
void generate_disambiguated_faces_for_plane_mt(
    tbb::affinity_partitioner& partitioner, ImplicitSurfacePolicy& isp,
    const frantic::volumetrics::marching_cubes_table& mct, const frantic::graphics2d::size2& size,
    const std::vector<boost::uint8_t>& previousCubeCases, const std::vector<boost::uint8_t>& currentCubeCases,
    const std::vector<boost::int32_t>& previousVertexIndices, const std::vector<boost::int32_t>& currentVertexIndices,
    const std::vector<float>& previousVoxelCornerDensities, const std::vector<float>& voxelCornerDensities,
    const std::vector<boost::int32_t>& newFaceCount,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
    frantic::geometry::trimesh3& outMesh ) {
    tbb::spin_mutex meshMutex;

    const std::size_t firstFaceNumber = outMesh.face_count();

    if( newFaceCount.size() ) {
        outMesh.add_faces( newFaceCount.back() );
    }

    generate_disambiguated_faces_for_plane_body<ImplicitSurfacePolicy> body(
        isp, mct, size, previousCubeCases, currentCubeCases, previousVertexIndices, currentVertexIndices,
        previousVoxelCornerDensities, voxelCornerDensities, firstFaceNumber, newFaceCount, meshMutex,
        outNamedVertexChannels, outMesh );

    tbb::blocked_range<std::size_t> range( 1, size.ysize );
#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body, partitioner );
#else
#pragma message( "Threads are disabled" )
    body( range );
#endif
}

template <class ImplicitSurfacePolicy>
void generate_disambiguated_faces_for_plane(
    ImplicitSurfacePolicy& isp, const frantic::volumetrics::marching_cubes_table& mct, frantic::graphics2d::size2 size,
    const std::vector<boost::uint8_t>& previousCubeCases, const std::vector<boost::uint8_t>& currentCubeCases,
    const std::vector<boost::int32_t>& previousVertexIndices, const std::vector<boost::int32_t>& currentVertexIndices,
    const std::vector<float>& previousVoxelCornerDensities, const std::vector<float>& voxelCornerDensities,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
    frantic::geometry::trimesh3& outMesh ) {
    // Allocate a vector for the 12 cube vertices, and the faces
    std::vector<boost::int32_t> cubeVerts( 12 );
    std::vector<frantic::graphics::vector3> cubeFaces;

    // DEBUG
    // a face channel accessor for indicating which faces were produced by which cube cases
    // trimesh3_face_channel_accessor<frantic::graphics::vector3> configAcc =
    // outMesh.get_face_channel_accessor<frantic::graphics::vector3>("CubeConfig"); trimesh3_face_channel_accessor<int>
    // caseAcc = outMesh.get_face_channel_accessor<int>("CubeCase");

    // avoid secure scl checks
    // TODO: get rid of this when you're no longer using msvc8/9
    const boost::uint8_t* pPreviousCubeCases = previousCubeCases.size() > 0 ? &previousCubeCases[0] : 0;
    const boost::uint8_t* pCurrentCubeCases = currentCubeCases.size() > 0 ? &currentCubeCases[0] : 0;
    const boost::int32_t* pPreviousVertexIndices = previousVertexIndices.size() > 0 ? &previousVertexIndices[0] : 0;
    const boost::int32_t* pCurrentVertexIndices = currentVertexIndices.size() > 0 ? &currentVertexIndices[0] : 0;
    const float* pPreviousVoxelCornerDensities =
        previousVoxelCornerDensities.size() > 0 ? &previousVoxelCornerDensities[0] : 0;
    const float* pVoxelCornerDensities = voxelCornerDensities.size() > 0 ? &voxelCornerDensities[0] : 0;

    // Skip the first row
    int index = size.xsize;
    float cubeCorners[8];
    for( int y = 1; y < size.ysize; ++y ) {
        // Skip the first column
        ++index;
        for( int x = 1; x < size.xsize; ++x ) {
            unsigned char cubeCase = pCurrentCubeCases[index];
            if( cubeCase != 0x00 && cubeCase != 0xff ) {

                // Retrieve all the vertices from the edges of this cube
                mct.fill_verts( cubeCase, pCurrentVertexIndices[index], pCurrentCubeCases[index - 1],
                                pCurrentCubeCases[index - size.xsize], pPreviousCubeCases[index],
                                pCurrentCubeCases[index - size.xsize - 1], pPreviousCubeCases[index - 1],
                                pPreviousCubeCases[index - size.xsize], pCurrentVertexIndices[index - 1],
                                pCurrentVertexIndices[index - size.xsize], pPreviousVertexIndices[index],
                                pCurrentVertexIndices[index - size.xsize - 1], pPreviousVertexIndices[index - 1],
                                pPreviousVertexIndices[index - size.xsize], cubeVerts );

                // retrieve the corner values for disambiguation
                cubeCorners[0] = pPreviousVoxelCornerDensities[index - size.xsize - 1];
                cubeCorners[1] = pPreviousVoxelCornerDensities[index - size.xsize];
                cubeCorners[2] = pPreviousVoxelCornerDensities[index - 1];
                cubeCorners[3] = pPreviousVoxelCornerDensities[index];
                cubeCorners[4] = pVoxelCornerDensities[index - size.xsize - 1];
                cubeCorners[5] = pVoxelCornerDensities[index - size.xsize];
                cubeCorners[6] = pVoxelCornerDensities[index - 1];
                cubeCorners[7] = pVoxelCornerDensities[index];

                // get the faces for this cube case.
                // it may be necessary to add another vert to the mesh.
                //				frantic::graphics::vector3 cubeCaseDebug;
                bool newVert = mct.get_cubecase_faces( cubeCase, cubeCorners, cubeFaces );
                //				bool newVert = mct.get_cubecase_faces(cubeCase, cubeCorners, cubeFaces,
                // cubeCaseDebug);
                if( newVert ) {

                    // average the existing verts for the new vert location
                    // TODO:  a policy based vertex addition would probably be more appropriate with a
                    // solve to find a vert location
                    frantic::graphics::vector3f avg;
                    int count = 0;
                    for( size_t i = 0; i < cubeVerts.size(); ++i ) {
                        if( cubeVerts[i] >= 0 ) {
                            avg += outMesh.get_vertex( cubeVerts[i] );
                            ++count;
                        }
                    }
                    avg = avg * ( 1.f / count );
                    cubeVerts.push_back( (int)outMesh.vertex_count() );
                    outMesh.add_vertex( avg );
                    if( !outNamedVertexChannels.empty() )
                        isp.add_vertex_to_channels( outNamedVertexChannels );
                }

                frantic::graphics::vector3 face;
                for( unsigned i = 0; i < cubeFaces.size(); ++i ) {
                    face = cubeFaces[i];
                    // std::cout << "face: " << face << std::endl;
                    outMesh.add_face(
                        frantic::graphics::vector3( cubeVerts[face.x], cubeVerts[face.y], cubeVerts[face.z] ) );
                    //					caseAcc.add_face(); caseAcc[caseAcc.size()-1] = cubeCase;
                    //					configAcc.add_face(); configAcc[configAcc.size()-1] =
                    // cubeCaseDebug;
                }
            }
            ++index;
        }
    }
}

// This function adds the new vertices generated by the addition of another layer of cubes
template <class ImplicitSurfacePolicy>
void generate_vertices_for_plane(
    ImplicitSurfacePolicy& isp, const frantic::graphics2d::boundrect2& xyExtents, int z,
    const std::vector<float>& previousVoxelCornerDensities, const std::vector<float>& currentVoxelCornerDensities,
    const std::vector<unsigned char>& cubeCases, std::vector<boost::int32_t>& outVertexIndices,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
    frantic::geometry::trimesh3& outMesh ) {
    outVertexIndices.resize( xyExtents.get_area() );
    int i = 0;

    // Loop over all the values except the initial ones.  This is because the values at the beginning of an axis don't
    // generate any vertices by construction.  Note that the places where out-of-bounds index access could occur won't,
    // because the duplication of density values for the initial row, column, and layer prevent it.
    for( int y = xyExtents.minimum().y; y <= xyExtents.maximum().y; ++y ) {
        for( int x = xyExtents.minimum().x; x <= xyExtents.maximum().x; ++x ) {
            unsigned char cubeCase = cubeCases[i];
            if( cubeCase != 0x00 && cubeCase != 0xff ) {
                outVertexIndices[i] = (boost::int32_t)outMesh.vertex_count();
                // Get the new vertices created by this cube case
                frantic::graphics::vector3 newVerts =
                    marching_cubes_table::get_new_verts( cubeCase, (int)outMesh.vertex_count() );
                // Actually add the vertices to the mesh
                if( newVerts.x >= 0 ) { // Edge 67, along the X axis
                    float density0 = currentVoxelCornerDensities[i - 1], density1 = currentVoxelCornerDensities[i];
                    frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x - 1, y, z ),
                                               voxelCorner1 = frantic::graphics::vector3( x, y, z );
                    isp.find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1 );
                    outMesh.add_vertex( isp.get_isosurface_location_world() );
                    if( !outNamedVertexChannels.empty() )
                        isp.add_vertex_to_channels( outNamedVertexChannels );
                }
                if( newVerts.y >= 0 ) { // Edge 57, along the Y axis
                    float density0 = currentVoxelCornerDensities[i - xyExtents.xsize()],
                          density1 = currentVoxelCornerDensities[i];
                    frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x, y - 1, z ),
                                               voxelCorner1 = frantic::graphics::vector3( x, y, z );
                    isp.find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1 );
                    outMesh.add_vertex( isp.get_isosurface_location_world() );
                    if( !outNamedVertexChannels.empty() )
                        isp.add_vertex_to_channels( outNamedVertexChannels );
                }
                if( newVerts.z >= 0 ) { // Edge 37, along the Z axis
                    float density0 = previousVoxelCornerDensities[i], density1 = currentVoxelCornerDensities[i];
                    frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x, y, z - 1 ),
                                               voxelCorner1 = frantic::graphics::vector3( x, y, z );
                    isp.find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1 );
                    outMesh.add_vertex( isp.get_isosurface_location_world() );
                    if( !outNamedVertexChannels.empty() )
                        isp.add_vertex_to_channels( outNamedVertexChannels );
                }
            }
            ++i;
        }
    }
}

/**
 *	This function adds the new vertices generated by the addition of another layer of cubes.
 *	It does so using sparse plane information.
 *
 *	This is a little more complicated than the dense case, because I dont know if required verts have
 *	been generated yet.  Previous passes may have generated verts that I can reuse, the only way
 *	to know is to check the previous plane vertex indices.  These should already be initialized by previous
 *	function calls and may or may not already have a valid vertex index for a given location.
 *
 *	The verts in the mesh will be generated in a different order from the dense version, but the index information
 *	is kept track of in the outVertexIndices frantic::graphics::vector3s.  Any time a "previous" edge is encountered
 *that requires a vert,	if a valid index into the mesh verts is not already present, a vert is generated, added to the
 *mesh and the appropriate vertex index is stored.
 *
 *	@param	isp	The implicit surface policy used to solve for vert locations when they are required.
 *	@param	z	The current plane.
 *	@param	previousVoxelCornerDensities	corner samples from the previous plane
 *	@param	currentVoxelCornerDensities		corner samples from the current plane
 *	@param	definedCubeCases				defined cube cases for this plane
 *	@param	definedCubeCasesRLP
 *	@param	prevVertexIndices				vert index flags/values from the previous plane
 *	@param	currentVertexIndices			vert index flags/values from the current plane
 *	@param	outNamedVertexChannels			vert channels for the output mesh
 *	@param	outMesh							the output mesh to add the verts to
 */
template <class ImplicitSurfacePolicy>
void generate_vertices_for_sparse_plane(
    ImplicitSurfacePolicy& isp, int z, const boost::shared_array<float> previousVoxelCornerDensities,
    const boost::shared_array<float> currentVoxelCornerDensities,
    const boost::shared_array<unsigned char> definedCubeCases,
    const frantic::volumetrics::rle_plane& definedCubeCasesRLP,
    boost::shared_array<frantic::graphics::vector3>& prevVertexIndices,
    boost::shared_array<frantic::graphics::vector3>& outVertexIndices,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
    frantic::geometry::trimesh3& outMesh ) {
    // std::cout << std::endl << "generate_vertices_for_sparse_plane" << std::endl;
    // definedCubeCasesRLP.dump(std::cout);

    if( definedCubeCasesRLP.is_empty() )
        return;

    // Loop over all the defined cube cases using the rle plane structure.
    int xsize = definedCubeCasesRLP.get_xy_extents().xsize();

    char vertFlags[12];

    for( int yExtent = 1; yExtent < definedCubeCasesRLP.get_xy_extents().ysize(); ++yExtent ) {
        int run = definedCubeCasesRLP.get_y_extent_first_run( yExtent );
        if( run < 0 ) // skip empty extents
            continue;
        int y = yExtent + definedCubeCasesRLP.get_xy_extents().minimum().y;
        int xOffset = yExtent * xsize - definedCubeCasesRLP.get_xy_extents().minimum().x;
        while( run <= definedCubeCasesRLP.get_y_extent_last_run( yExtent ) ) {

            for( int i = definedCubeCasesRLP[run].first; i <= definedCubeCasesRLP[run].second; ++i ) {
                unsigned char cubeCase = definedCubeCases[i];
                if( cubeCase != 0x00 && cubeCase != 0xff ) {
                    // std::cout << "x: " << i << "  cubecase: " << std::hex << (int)cubeCase << std::dec << std::endl;
                    int x = i - xOffset;

                    const int edgeForAxis[3] = { 11, 10, 7 };
                    for( int axis = 0; axis < 3; ++axis ) {
                        // Determine which verts are needed for this case
                        marching_cubes_table::get_required_vert_flags( cubeCase, vertFlags );

                        const int edge = edgeForAxis[axis];

                        if( vertFlags[edge] ) {
                            if( outVertexIndices[i][axis] < 0 ) {
                                int vertCount = (int)outMesh.vertex_count();
                                outVertexIndices[i][axis] = vertCount;

                                // Actually add the vertex to the mesh
                                float density0, density1;
                                frantic::graphics::vector3 voxelCorner0, voxelCorner1;
                                switch( axis ) {
                                case 0:
                                    density0 = currentVoxelCornerDensities[i - 1];
                                    density1 = currentVoxelCornerDensities[i];
                                    voxelCorner0.set( x - 1, y, z );
                                    voxelCorner1.set( x, y, z );
                                    break;
                                case 1:
                                    density0 = currentVoxelCornerDensities[i - xsize];
                                    density1 = currentVoxelCornerDensities[i];
                                    voxelCorner0.set( x, y - 1, z );
                                    voxelCorner1.set( x, y, z );
                                    break;
                                case 2:
                                    density0 = previousVoxelCornerDensities[i];
                                    density1 = currentVoxelCornerDensities[i];
                                    voxelCorner0.set( x, y, z - 1 );
                                    voxelCorner1.set( x, y, z );
                                    break;
                                }

                                isp.find_edge_isosurface_location( density0, density1, voxelCorner0, voxelCorner1 );
                                outMesh.add_vertex( isp.get_isosurface_location_world() );
                                if( !outNamedVertexChannels.empty() )
                                    isp.add_vertex_to_channels( outNamedVertexChannels );
                            }
                        }
                    }
                }
            }
            run++;
        }
    }
}

void copy_color_vertices( const std::vector<exposed_voxel_vertices>& colorExposedVertices,
                          voxel_vertex_t& exposedVoxelVertices );

template <class ImplicitSurfacePolicy>
struct generate_block_geometry_impl {
    ImplicitSurfacePolicy& m_isp;
    const frantic::volumetrics::marching_cubes_table& m_mct;

    const std::vector<frantic::graphics::vector3>& m_meshingBlocks;
    const frantic::graphics::size3 m_meshingBlockSize;
    const std::size_t m_bufferVertsSoftLimit;
    const std::size_t m_bufferFacesSoftLimit;

    const std::vector<frantic::tstring>& m_outVertexChannelNames;

    std::vector<frantic::graphics::vector3f>& m_externalVertices;

    const voxel_vertex_t& m_exposedVoxelVertices;

    tbb::spin_mutex& m_outMeshMutex;
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& m_outVertexChannels;
    frantic::geometry::trimesh3& m_outMesh;
    boost::int32_t& m_nextVertexNumber;
    boost::int32_t& m_nextFaceNumber;
    std::vector<exposed_voxel_vertices>& m_outExposedVertices;

    shared_progress_logger_proxy& m_progressLogger;

    int m_blockColor;

    generate_block_geometry_impl<ImplicitSurfacePolicy>&
    operator=( const generate_block_geometry_impl<ImplicitSurfacePolicy>& ); // not implemented

    void operator()( const tbb::blocked_range<std::size_t>& r ) const {
        frantic::geometry::trimesh3 bufferMesh;
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor> bufferChannels;

        // frantic::diagnostics::profiling_manager prof;
        // prof.new_profiling_section( "total" );
        // prof.new_profiling_section( "density" );
        // prof.new_profiling_section( "cube cases" );
        // prof.new_profiling_section( "vertices" );
        // prof.new_profiling_section( "faces" );
        // prof.new_profiling_section( "init" );
        // prof.new_profiling_section( "face vert lookup" );
        // prof.new_profiling_section( "face vert" );
        // prof.new_profiling_section( "face output" );
        // prof.new_profiling_section( "face lookup" );

        // prof.enter_section( 0 );

        // TODO: change nextVertexNumber / nextFaceNumber to atomics ?
        for( std::size_t i = 0; i < m_outVertexChannels.size(); ++i ) {
            bufferMesh.add_vertex_channel_raw( m_outVertexChannelNames[i], m_outVertexChannels[i].arity(),
                                               m_outVertexChannels[i].data_type() );
            bufferChannels.push_back( bufferMesh.get_vertex_channel_general_accessor( m_outVertexChannelNames[i] ) );
        }

        const frantic::graphics::size3 densityBlockSize( m_meshingBlockSize.xsize() + 1, m_meshingBlockSize.ysize() + 1,
                                                         m_meshingBlockSize.zsize() + 1 );
        const int densityBlockVolume = densityBlockSize.volume();

        std::vector<float> voxelCornerDensities( densityBlockVolume );
        std::vector<frantic::graphics::vector3> vertexIndices( densityBlockVolume );
        std::vector<boost::uint8_t> cubeCases( densityBlockVolume );

        typename ImplicitSurfacePolicy::sparse_voxel_corner_density_workspace_t sparseVoxelCornerDensityWorkspace;
        typename ImplicitSurfacePolicy::vertex_workspace_t vertexWorkspace;

        std::vector<exposed_voxel_vertices> exposedVertexBuffer;

        for( std::size_t i = r.begin(); i != r.end(); ++i ) {
            if( m_progressLogger.is_cancelled() ) {
                break;
            }
            if( bufferMesh.vertex_count() > m_bufferVertsSoftLimit ||
                bufferMesh.face_count() > m_bufferFacesSoftLimit ) {
                const boost::int32_t firstVertexNumber =
                    detail::flush_mesh_buffer( bufferChannels, bufferMesh, m_outMeshMutex, m_nextVertexNumber,
                                               m_nextFaceNumber, m_outVertexChannels, m_outMesh );
                detail::flush_exposed_vertex_buffer( exposedVertexBuffer, firstVertexNumber, m_outMeshMutex,
                                                     m_outExposedVertices );
            }
            const frantic::graphics::vector3& meshingBlock = m_meshingBlocks[i];
            frantic::graphics::boundbox3 blockExtents(
                frantic::graphics::vector3( meshingBlock.x * m_meshingBlockSize.xsize(),
                                            meshingBlock.y * m_meshingBlockSize.ysize(),
                                            meshingBlock.z * m_meshingBlockSize.zsize() ),
                m_meshingBlockSize );
            frantic::graphics::boundbox3 densityExtents( blockExtents.minimum() - frantic::graphics::vector3( 1 ),
                                                         densityBlockSize );

            // prof.enter_section( 1 );
            m_isp.fill_voxel_corner_densities( densityExtents, &voxelCornerDensities[0],
                                               sparseVoxelCornerDensityWorkspace );
            // prof.exit_section( 1 );
            // prof.enter_section( 2 );
            detail::fill_cube_case_values( densityBlockSize, voxelCornerDensities, cubeCases );
            // prof.exit_section( 2 );
            // prof.enter_section( 3 );
            detail::generate_vertices_for_box( m_isp, m_mct, meshingBlock, densityExtents, m_blockColor,
                                               voxelCornerDensities, cubeCases, m_exposedVoxelVertices,
                                               /*m_blockToExposedVoxelVertices,*/ vertexIndices, bufferChannels,
                                               bufferMesh, exposedVertexBuffer, vertexWorkspace );
            // prof.exit_section( 3 );
            // prof.enter_section( 4 );
            detail::generate_disambiguated_faces_for_box( m_isp, m_mct, densityExtents, voxelCornerDensities, cubeCases,
                                                          vertexIndices, m_externalVertices, bufferChannels, bufferMesh,
                                                          vertexWorkspace /*, prof*/ );
            // prof.exit_section( 4 );

            m_progressLogger.inc_progress();
        }

        const boost::int32_t firstVertexNumber =
            detail::flush_mesh_buffer( bufferChannels, bufferMesh, m_outMeshMutex, m_nextVertexNumber, m_nextFaceNumber,
                                       m_outVertexChannels, m_outMesh );
        detail::flush_exposed_vertex_buffer( exposedVertexBuffer, firstVertexNumber, m_outMeshMutex,
                                             m_outExposedVertices );

        // prof.exit_section( 0 );
    }

    generate_block_geometry_impl(
        ImplicitSurfacePolicy& isp, frantic::volumetrics::marching_cubes_table& mct, const int blockColor,
        const std::vector<frantic::graphics::vector3>& meshingBlocks, const frantic::graphics::size3 meshingBlockSize,
        const std::size_t bufferVertsSoftLimit, const std::size_t bufferFacesSoftLimit,
        // const std::vector<exposed_voxel_vertices> & exposedVoxelVertices,
        // const stdext::hash_map<vector3,std::pair<size_t,size_t>,voxel_coord_hasher> & blockToExposedVoxelVertices,
        const voxel_vertex_t& exposedVoxelVertices, std::vector<frantic::graphics::vector3f>& externalVertices,
        std::vector<frantic::tstring>& outVertexChannelNames,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outVertexChannels,
        frantic::geometry::trimesh3& outMesh, boost::int32_t& nextVertexNumber, boost::int32_t& nextFaceNumber,
        std::vector<exposed_voxel_vertices>& outExposedVertices, tbb::spin_mutex& outMeshMutex,
        // tbb::atomic<std::size_t> & progress,
        // frantic::logging::progress_logger & progressLogger,
        shared_progress_logger_proxy& progressLogger
        // tbb::queuing_mutex & progressMutex,
        // tbb::tbb_thread::id guiThreadID/*,
        /*tbb::task_group_context & taskGroupContext*/ )
        : m_isp( isp )
        , m_mct( mct )
        , m_blockColor( blockColor )
        , m_meshingBlocks( meshingBlocks )
        , m_meshingBlockSize( meshingBlockSize )
        , m_bufferVertsSoftLimit( bufferVertsSoftLimit )
        , m_bufferFacesSoftLimit( bufferFacesSoftLimit )
        , m_exposedVoxelVertices( exposedVoxelVertices )
        ,
        // m_blockToExposedVoxelVertices( blockToExposedVoxelVertices ),
        m_externalVertices( externalVertices )
        , m_outVertexChannelNames( outVertexChannelNames )
        , m_outVertexChannels( outVertexChannels )
        , m_outMesh( outMesh )
        , m_nextVertexNumber( nextVertexNumber )
        , m_nextFaceNumber( nextFaceNumber )
        , m_outExposedVertices( outExposedVertices )
        , m_outMeshMutex( outMeshMutex )
        ,
        // m_progress( progress ),
        m_progressLogger( progressLogger ) //,
    // m_progressMutex( progressMutex ),
    // m_guiThreadID( guiThreadID )/*,
    // m_taskGroupContext( taskGroupContext )*/
    {}
};

template <class ImplicitSurfacePolicy>
void generate_block_geometry_mt(
    ImplicitSurfacePolicy& isp, frantic::volumetrics::marching_cubes_table& mct,
    std::vector<frantic::graphics::vector3>& meshingBlocks, const frantic::graphics::size3 meshingBlockSize,
    const std::size_t bufferVertsSoftLimit, const std::size_t bufferFacesSoftLimit,
    std::vector<frantic::tstring>& outVertexChannelNames,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outVertexChannels,
    frantic::geometry::trimesh3& outMesh, shared_progress_logger_proxy& progressLogger ) {
    // tbb::task_group_context taskGroupContext( tbb::task_group_context::isolated );
    tbb::spin_mutex outMeshMutex;
    tbb::auto_partitioner partitioner;

    // frantic::diagnostics::profiling_manager prof;
    // prof.new_profiling_section( "geometry - init" );
    // prof.new_profiling_section( "geometry - mt" );
    // prof.new_profiling_section( "geometry - st" );

    // tbb::atomic<boost::int32_t> nextVertexNumber;
    // tbb::atomic<boost::int32_t> nextFaceNumber;

    // prof.enter_section( 0 );

    std::vector<frantic::graphics::vector3f> externalVertices;

    // std::vector<exposed_block_vertices> colorExposedVertices;
    std::vector<exposed_voxel_vertices> colorExposedVertices;

    // sort the blocks by their corner position in non-overlapping 2x2x2 blocks
    // all blocks in the same corner can be processed in parallel
    // cornerStartIndex holds start and end index for each corner position

    const std::size_t colorCount = 8;
    boost::array<std::size_t, 8> blockColorCount;
    blockColorCount.assign( 0 );

    BOOST_FOREACH( const frantic::graphics::vector3& coord, meshingBlocks ) {
        ++blockColorCount[get_block_color( coord )];
    }

    boost::array<std::size_t, 9> colorStartIndex;
    colorStartIndex[0] = 0;
    for( int corner = 0; corner < 8; ++corner ) {
        colorStartIndex[corner + 1] = colorStartIndex[corner] + blockColorCount[corner];
    }

    {
        std::vector<frantic::graphics::vector3> sortedMeshingBlocks( meshingBlocks.size() );
        boost::array<std::size_t, 8> currentCornerIndex;
        std::copy( colorStartIndex.begin(), colorStartIndex.begin() + 8, currentCornerIndex.begin() );
        BOOST_FOREACH( const frantic::graphics::vector3& coord, meshingBlocks ) {
            const int corner = get_block_color( coord );
            sortedMeshingBlocks[currentCornerIndex[corner]] = coord;
            ++currentCornerIndex[corner];
        }
        std::swap( meshingBlocks, sortedMeshingBlocks );
    }

    boost::int32_t nextVertexNumber = 0;
    boost::int32_t nextFaceNumber = 0;

    // std::vector<exposed_voxel_vertices> exposedVoxelVertices;
    //  todo : separate block colors ?
    // stdext::hash_map<vector3,std::pair<size_t,size_t>,voxel_coord_hasher> blockToExposedVoxelVertices;

    voxel_vertex_t exposedVoxelVertices;
    // exposedVoxelVertices.set_empty_key( frantic::graphics::vector3( std::numeric_limits<boost::int32_t>::max() ) );

    // prof.exit_section( 0 );

    progressLogger.set_progress_end( meshingBlocks.size() );

    for( int color = 0; color < colorCount; ++color ) {
        if( progressLogger.is_cancelled() ) {
            break;
        }

        // FF_LOG( stats ) << "color: " << boost::lexical_cast<std::string>( color ) << " blocks: " <<
        // boost::lexical_cast<std::string>( blockColorCount[color] ) << std::endl; prof.enter_section( 1 );
        tbb::blocked_range<std::size_t> range( colorStartIndex[color], colorStartIndex[color + 1] );

        colorExposedVertices.clear();

        generate_block_geometry_impl<ImplicitSurfacePolicy> body(
            isp, mct, color, meshingBlocks, meshingBlockSize, bufferVertsSoftLimit, bufferFacesSoftLimit,
            exposedVoxelVertices, /*blockToExposedVoxelVertices,*/ externalVertices, outVertexChannelNames,
            outVertexChannels, outMesh, nextVertexNumber, nextFaceNumber, colorExposedVertices, outMeshMutex,
            progressLogger );

#ifndef FRANTIC_DISABLE_THREADS
        tbb::parallel_for( range, body, partitioner /*, taskGroupContext*/ );
#else
#pragma message( "Threads are disabled" )
        body( range );
#endif
        // prof.exit_section( 1 );
        // prof.enter_section( 2 );
        outMesh.set_vertex_count( nextVertexNumber );
        BOOST_FOREACH( frantic::geometry::trimesh3_vertex_channel_general_accessor& outputChannel, outVertexChannels ) {
            outputChannel.set_vertex_count( nextVertexNumber );
        }
        outMesh.set_face_count( nextFaceNumber );

        copy_color_vertices( colorExposedVertices, exposedVoxelVertices /*, blockToExposedVoxelVertices*/ );
        std::copy( outMesh.vertices_ref().begin() + externalVertices.size(), outMesh.vertices_ref().end(),
                   std::back_insert_iterator<std::vector<vector3f>>( externalVertices ) );
        // prof.exit_section( 2 );

        // progressLogger.update_progress( progress, meshingBlocks.size() );
    }

    // FF_LOG( stats ) << prof << std::endl;

    // if( taskGroupContext.is_group_execution_cancelled() ) {
    // throw frantic::logging::progress_cancel_exception( "Cancelled Operation: Frost" );
    //}
};

inline boost::int32_t round_up_to_multiple_of_float_v_size( boost::int32_t i ) {
    const std::size_t width = frantic::simd::float_v::static_size;
    return ( i + width - 1 ) / width * width;
}

} // namespace detail

/**
 * Converts an implicit surface into a triangle mesh.  It handles the proper copying of named channels,
 * whose input is specified by the implicit surface policy.
 *
 * @param  isp              The implicit surface policy evaluats the implicit surface function for the algorithm.
 * @param  outMesh          The output parameter where the resulting mesh will go.
 */
template <class ImplicitSurfacePolicy>
void convert_dense_implicit_surface_to_trimesh3( frantic::channels::channel_propagation_policy& cpp,
                                                 ImplicitSurfacePolicy& mcp, frantic::geometry::trimesh3& outMesh ) {
    static marching_cubes_table mct;

    tbb::task_scheduler_init taskSchedulerInit;
    tbb::affinity_partitioner partitioner;

    // Clear the mesh to start
    outMesh.clear();

    // Initialize the named channels we'll be creating
    std::vector<frantic::tstring> channelNames;
    std::vector<std::pair<size_t, frantic::channels::data_type_t>> channelDataTypes;
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor> outputChannels;
    mcp.get_channel_names( cpp, channelNames );
    // This prepares the output channel general accessors contained within the mc policy, and
    // returns the primitive data types of all the channels.
    mcp.prepare_channel_accessors( channelNames, channelDataTypes );
    // Create all the channels in the output mesh, as well as the vector of writable channel accessors used
    // to fill in the data.
    for( unsigned i = 0; i < channelNames.size(); ++i ) {
        outMesh.add_vertex_channel_raw( channelNames[i], channelDataTypes[i].first, channelDataTypes[i].second );
        outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( channelNames[i] ) );
    }

    // DEBUG
    // a face channel to indicate which cubecase/configuration produced this face.
    // for tracking down bad cubecases that produce bad meshes
    // outMesh.add_face_channel<frantic::graphics::vector3>("CubeConfig");
    // outMesh.add_face_channel<int>("CubeCase");

    // These are the bound we use to scan.
    frantic::graphics::boundbox3 voxelBounds = mcp.get_voxel_bounds();
    frantic::graphics2d::boundrect2 xyExtents( voxelBounds.minimum().x, voxelBounds.maximum().x,
                                               voxelBounds.minimum().y, voxelBounds.maximum().y );
    int zMin = voxelBounds.minimum().z, zMax = voxelBounds.maximum().z;

    frantic::graphics2d::size2 densitySampleDimensions = xyExtents.size();

    // If this is empty, then let's skip all the computation below.
    if( densitySampleDimensions.get_area() <= 0 )
        return;

    std::vector<float> previousVoxelCornerDensities( densitySampleDimensions.get_area() ),
        currentVoxelCornerDensities( densitySampleDimensions.get_area() );

    // Though normally the cube case dimensions should be one smaller, we're adding 1 along the x and y dimensions to
    // hold the cube cases for the duplicated samples at the beginning.  The first row and column are used to create
    // vertices, but are not used to generate faces.
    std::vector<boost::uint8_t> previousCubeCases( densitySampleDimensions.get_area() ),
        currentCubeCases( densitySampleDimensions.get_area() );
    std::vector<boost::int32_t> previousVertexIndices( densitySampleDimensions.get_area() ),
        currentVertexIndices( densitySampleDimensions.get_area() );

    // Get the voxel corner densities at the smallest Z
    mcp.fill_voxel_corner_densities( xyExtents, zMin, previousVoxelCornerDensities );

    // Compute the cube case values based on the voxel corner densities.  Subsequent calculations of the cube cases
    // use previously calculated cube cases, to optimize the computation
    detail::fill_initial_cube_case_values( densitySampleDimensions, previousVoxelCornerDensities, previousCubeCases );

    // Generate the new vertices that would be generate for this set of cubes.  This primes the algorithm for the next
    // step, and we don't generate faces for these cubes.  Rather, the vertices generated are used by the faces in the
    // next plane. For this case, the same plane of voxel corner densities is provided twice because of the boundary
    // conditions we use.
    // detail::generate_vertices_for_plane_mt( partitioner, mcp, xyExtents, zMin, previousVoxelCornerDensities,
    // previousVoxelCornerDensities, previousCubeCases, previousVertexIndices, outputChannels, outMesh );
    detail::generate_vertices_for_plane( mcp, xyExtents, zMin, previousVoxelCornerDensities,
                                         previousVoxelCornerDensities, previousCubeCases, previousVertexIndices,
                                         outputChannels, outMesh );

    // Now go through all the slices of cubes
    for( int z = zMin + 1; z <= zMax; ++z ) {
        // Get the new slice of voxel density samples, and compute the cube cases that arise
        mcp.fill_voxel_corner_densities( xyExtents, z, currentVoxelCornerDensities );
        detail::fill_cube_case_values_mt( partitioner, densitySampleDimensions, previousCubeCases,
                                          currentVoxelCornerDensities, currentCubeCases );

        // Generate the vertices and then the faces for this slice plane
        // detail::generate_vertices_for_plane_mt( partitioner, mcp, xyExtents, z, previousVoxelCornerDensities,
        // currentVoxelCornerDensities, currentCubeCases, currentVertexIndices, outputChannels, outMesh );
        detail::generate_vertices_for_plane( mcp, xyExtents, z, previousVoxelCornerDensities,
                                             currentVoxelCornerDensities, currentCubeCases, currentVertexIndices,
                                             outputChannels, outMesh );
        // detail::generate_faces_for_plane( mct, densitySampleDimensions, previousCubeCases, currentCubeCases,
        // previousVertexIndices, currentVertexIndices, outMesh );
        detail::generate_disambiguated_faces_for_plane(
            mcp, mct, densitySampleDimensions, previousCubeCases, currentCubeCases, previousVertexIndices,
            currentVertexIndices, previousVoxelCornerDensities, currentVoxelCornerDensities, outputChannels, outMesh );

        // Swap the previous and current values for the next iteration of the loop
        previousVertexIndices.swap( currentVertexIndices );
        previousCubeCases.swap( currentCubeCases );
        previousVoxelCornerDensities.swap( currentVoxelCornerDensities );
    }
}
//*/

/**
 * Converts an implicit surface into a triangle mesh.  It handles the proper copying of named channels,
 * whose input is specified by the implicit surface policy, and limited by the channel propagation policy.
 * This function uses the "sparse" implicit surface policy functions whose computational complexity is data
 * driven, rather than volume driven.  The "sparse" functions are also multithreaded, and this function
 * also makes use of the multithreaded generate_vertices_for_plane_mt function, which is a major bottle
 * neck in mesh generation. A null progress logger is used.
 *
 * @param  mcp              The implicit surface policy evaluats the implicit surface function for the algorithm.
 * @param  outMesh          The output parameter where the resulting mesh will go.
 */
template <class ImplicitSurfacePolicy>
void convert_implicit_surface_to_trimesh3( const frantic::channels::channel_propagation_policy& cpp,
                                           ImplicitSurfacePolicy& mcp, frantic::geometry::trimesh3& outMesh ) {
    frantic::logging::null_progress_logger nullLogger;
    convert_implicit_surface_to_trimesh3( cpp, mcp, outMesh, nullLogger );
}

/**
 * Converts an implicit surface into a triangle mesh.  It handles the proper copying of named channels,
 * whose input is specified by the implicit surface policy, and limited by the channel propagation policy.
 * This function uses the "sparse" implicit surface policy functions whose computational complexity is data
 * driven, rather than volume driven.  The "sparse" functions are also multithreaded, and this function
 * also makes use of the multithreaded generate_vertices_for_plane_mt function, which is a major bottle
 * neck in mesh generation.
 *
 * @param  mcp              The implicit surface policy evaluats the implicit surface function for the algorithm.
 * @param  outMesh          The output parameter where the resulting mesh will go.
 * @param  progressLogger	The object to accept progress on this funciton
 */
template <class ImplicitSurfacePolicy>
void convert_implicit_surface_to_trimesh3( const frantic::channels::channel_propagation_policy& cpp,
                                           ImplicitSurfacePolicy& mcp, frantic::geometry::trimesh3& outMesh,
                                           frantic::logging::progress_logger& progressLogger ) {
    static marching_cubes_table mct;

    tbb::task_scheduler_init taskSchedulerInit;
    tbb::affinity_partitioner partitioner;

    // frantic::diagnostics::profiling_manager prof;
    // prof.new_profiling_section( "total" );
    // prof.new_profiling_section( "density" );
    // prof.new_profiling_section( "cube cases" );
    // prof.new_profiling_section( "vertices" );
    // prof.new_profiling_section( "faces" );

    // prof.enter_section( 0 );

    // Clear the mesh to start
    outMesh.clear();

    // Initialize the named channels we'll be creating
    std::vector<frantic::tstring> channelNames;

    std::vector<std::pair<size_t, frantic::channels::data_type_t>> channelDataTypes;

    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor> outputChannels;

    mcp.get_channel_names( cpp, channelNames );
    // This prepares the output channel general accessors contained within the mc policy, and
    // returns the primitive data types of all the channels.
    mcp.prepare_channel_accessors( channelNames, channelDataTypes );

    // Create all the channels in the output mesh, as well as the vector of writable channel accessors used
    // to fill in the data.

    for( unsigned i = 0; i < channelNames.size(); ++i ) {
        outMesh.add_vertex_channel_raw( channelNames[i], channelDataTypes[i].first, channelDataTypes[i].second );
        outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( channelNames[i] ) );
    }

    // These are the bounds we use to scan.
    frantic::graphics::boundbox3 voxelBounds = mcp.get_voxel_bounds();
    // Round the x size up to allow for SIMD processing
    voxelBounds.maximum().x =
        voxelBounds.minimum().x + detail::round_up_to_multiple_of_float_v_size( voxelBounds.xsize() ) - 1;
    frantic::graphics2d::boundrect2 xyExtents( voxelBounds.minimum().x, voxelBounds.maximum().x,
                                               voxelBounds.minimum().y, voxelBounds.maximum().y );
    int zMin = voxelBounds.minimum().z, zMax = voxelBounds.maximum().z;
    const boost::int64_t zSliceCount = zMax - zMin + 1;

    frantic::graphics2d::size2 densitySampleDimensions = xyExtents.size();

    // If this is empty, then let's skip all the computation below.
    if( densitySampleDimensions.get_area() <= 0 || xyExtents.is_empty() )
        return;

    std::vector<float> previousVoxelCornerDensities( densitySampleDimensions.get_area() ),
        currentVoxelCornerDensities( densitySampleDimensions.get_area() );
    progressLogger.update_progress(
        100.f * 0.33f /
        ( zSliceCount + 1 ) ); // give a small amount of progress in case the memory allocations are taking a long time

    // Though normally the cube case dimensions should be one smaller, we're adding 1 along the x and y dimensions to
    // hold the cube cases for the duplicated samples at the beginning.  The first row and column are used to create
    // vertices, but are not used to generate faces.
    std::vector<boost::uint8_t> previousCubeCases( densitySampleDimensions.get_area() ),
        currentCubeCases( densitySampleDimensions.get_area() );
    progressLogger.update_progress( 100.f * 0.67f / ( zSliceCount + 1 ) );
    std::vector<boost::int32_t> previousVertexIndices( densitySampleDimensions.get_area() ),
        currentVertexIndices( densitySampleDimensions.get_area() );
    progressLogger.update_progress( 100.f / ( zSliceCount + 1 ) );

    // cumulative sum of the new vertices and faces in the current plane
    std::vector<boost::int32_t> newVertexCount( densitySampleDimensions.ysize );
    std::vector<boost::int32_t> newFaceCount( densitySampleDimensions.ysize );

    typename ImplicitSurfacePolicy::sparse_voxel_corner_density_workspace_t sparseVoxelCornerDensityWorkspace;

    // Get the voxel corner densities at the smallest Z
    // mcp.fill_voxel_corner_densities( xyExtents, zMin, previousVoxelCornerDensities );
    mcp.fill_sparse_voxel_corner_densities( xyExtents, zMin, &previousVoxelCornerDensities[0],
                                            sparseVoxelCornerDensityWorkspace );
    progressLogger.update_progress( 100.f * 1.33f / ( zSliceCount + 1 ) );

    // Compute the cube case values based on the voxel corner densities.  Subsequent calculations of the cube cases
    // use previously calculated cube cases, to optimize the computation
    detail::fill_initial_cube_case_values( densitySampleDimensions, previousVoxelCornerDensities, previousCubeCases,
                                           newVertexCount );
    progressLogger.update_progress( 100.f * 1.67f / ( zSliceCount + 1 ) );

    // Generate the new vertices that would be generate for this set of cubes.  This primes the algorithm for the next
    // step, and we don't generate faces for these cubes.  Rather, the vertices generated are used by the faces in the
    // next plane. For this case, the same plane of voxel corner densities is provided twice because of the boundary
    // conditions we use.
    detail::generate_vertices_for_plane_mt( partitioner, mcp, mct, xyExtents, zMin, previousVoxelCornerDensities,
                                            previousVoxelCornerDensities, previousCubeCases, newVertexCount,
                                            previousVertexIndices, outputChannels, newFaceCount, outMesh );
    progressLogger.update_progress( 2 * 100.f / ( zSliceCount + 1 ) );

    // Now go through all the slices of cubes

    // Use the progress logger to send update
    for( int z = zMin + 1; z <= zMax; ++z ) {
        // Get the new slice of voxel density samples, and compute the cube cases that arise
        // prof.enter_section( 1 );
        mcp.fill_sparse_voxel_corner_densities( xyExtents, z, &currentVoxelCornerDensities[0],
                                                sparseVoxelCornerDensityWorkspace );
        // prof.exit_section( 1 );
        // prof.enter_section( 2 );
        detail::fill_cube_case_values_mt( partitioner, densitySampleDimensions, previousCubeCases,
                                          currentVoxelCornerDensities, currentCubeCases, newVertexCount );
        // prof.exit_section( 2 );

        // Generate the vertices and then the faces for this slice plane
        // prof.enter_section( 3 );
        detail::generate_vertices_for_plane_mt( partitioner, mcp, mct, xyExtents, z, previousVoxelCornerDensities,
                                                currentVoxelCornerDensities, currentCubeCases, newVertexCount,
                                                currentVertexIndices, outputChannels, newFaceCount, outMesh );
        // prof.exit_section( 3 );
        // detail::generate_faces_for_plane( mct, densitySampleDimensions, previousCubeCases, currentCubeCases,
        // previousVertexIndices, currentVertexIndices, outMesh ); prof.enter_section( 4 );
        detail::generate_disambiguated_faces_for_plane_mt(
            partitioner, mcp, mct, densitySampleDimensions, previousCubeCases, currentCubeCases, previousVertexIndices,
            currentVertexIndices, previousVoxelCornerDensities, currentVoxelCornerDensities, newFaceCount,
            outputChannels, outMesh );
        // prof.exit_section( 4 );

        // Swap the previous and current values for the next iteration of the loop
        previousVertexIndices.swap( currentVertexIndices );
        previousCubeCases.swap( currentCubeCases );
        previousVoxelCornerDensities.swap( currentVoxelCornerDensities );
        progressLogger.update_progress( 100.f * ( 2 + z - zMin ) / ( zSliceCount + 1 ) );
    }
    // prof.exit_section( 0 );
    // FF_LOG( stats ) << prof << std::endl;
}

namespace detail {

inline unsigned char get_cube_case( float a, float b, float c, float d, float e, float f, float g, float h ) {
    return ( a < 0 ? 0x01 : 0x00 ) | ( b < 0 ? 0x02 : 0x00 ) | ( c < 0 ? 0x04 : 0x00 ) | ( d < 0 ? 0x08 : 0x00 ) |
           ( e < 0 ? 0x10 : 0x00 ) | ( f < 0 ? 0x20 : 0x00 ) | ( g < 0 ? 0x40 : 0x00 ) | ( h < 0 ? 0x80 : 0x00 );
}

struct voxel_vertex_cache {
    struct voxel_vertex_cache_entry {
        bool m_set;
        std::pair<std::size_t, std::size_t> range;
        std::size_t hint;
        voxel_vertex_cache_entry()
            : m_set( false ) {}
        bool is_set() { return m_set; }
        void reset( const std::pair<std::size_t, std::size_t>& newRange ) {
            m_set = true;
            hint = newRange.first;
            range = newRange;
        }
        void reset() {
            m_set = true;
            hint = 0;
            range = std::make_pair<std::size_t, std::size_t>( 0, 0 );
        }
    };

    boost::array<voxel_vertex_cache_entry, 26> m_cache;

    typedef const boost::unordered_map<frantic::graphics::vector3, std::pair<size_t, size_t>, voxel_coord_hasher>
        block_lut_t;

    const block_lut_t& m_blockToBoundaryVertices;
    const std::vector<exposed_voxel_vertices>& m_boundaryVertices;

    frantic::graphics::vector3 m_blockCoord;

    void fill_cache_entry( int direction ) {
        if( !m_cache[direction].is_set() ) {
            block_lut_t::const_iterator i =
                m_blockToBoundaryVertices.find( m_blockCoord + get_direction_neighbor_offset( direction ) );
            if( i == m_blockToBoundaryVertices.end() ) {
                m_cache[direction].reset();
            } else {
                m_cache[direction].reset( i->second );
            }
        }
    }
    voxel_vertex_cache& operator=( const voxel_vertex_cache& ); // not implemented
  public:
    voxel_vertex_cache( const frantic::graphics::vector3& blockCoord,
                        const std::vector<exposed_voxel_vertices>& boundaryVertices,
                        block_lut_t& blockToBoundaryVertices )
        : m_blockCoord( blockCoord )
        , m_blockToBoundaryVertices( blockToBoundaryVertices )
        , m_boundaryVertices( boundaryVertices ) {}
    bool try_get( int direction, const frantic::graphics::vector3& voxelCoord,
                  frantic::graphics::vector3& outVertices ) {
        fill_cache_entry( direction );
        for( std::size_t i = m_cache[direction].hint; i != m_cache[direction].range.second; ++i ) {
            if( m_boundaryVertices[i].voxelCoord == voxelCoord ) {
                m_cache[direction].hint = i; // todo : inc
                outVertices = m_boundaryVertices[i].vertices;
                return true;
            }
        }
        for( std::size_t i = m_cache[direction].range.first; i != m_cache[direction].hint; ++i ) {
            if( m_boundaryVertices[i].voxelCoord == voxelCoord ) {
                m_cache[direction].hint = i;
                outVertices = m_boundaryVertices[i].vertices;
                return true;
            }
        }
        return false;
    }
};

int get_boundary_owner_direction_code( int blockColor, int boundaryCase );

template <class ImplicitSurfacePolicy>
class block_vertex_manager {
    const ImplicitSurfacePolicy& m_isp;
    const frantic::volumetrics::marching_cubes_table& m_mct;
    void get_corner_vert_flags( boost::uint8_t cubeCase, boost::uint8_t flags[3] ) {
        flags[0] = ( 0x40 & ( ( cubeCase >> 1 ) ^ cubeCase ) ) != 0; // edge 67
        flags[1] = ( 0x20 & ( ( cubeCase >> 2 ) ^ cubeCase ) ) != 0; // edge 57
        flags[2] = ( 0x08 & ( ( cubeCase >> 4 ) ^ cubeCase ) ) != 0; // edge 37
    }
    const std::vector<float>& m_voxelCornerDensities;
    tbb::spin_rw_mutex m_outMeshMutex;
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& m_outNamedVertexChannels;
    frantic::geometry::trimesh3& m_outMesh;

    // const std::vector<exposed_voxel_vertices> & m_inVoxelVertices;
    const voxel_vertex_t& m_exposedVoxelVertices;

    std::vector<exposed_voxel_vertices>& m_outVoxelVertices;

    std::vector<frantic::graphics::vector3>& m_outVertexIndices;

    // voxel_vertex_cache m_voxelVertexCache;

    typename ImplicitSurfacePolicy::vertex_workspace_t& m_vertexWorkspace;

    int m_xsize;
    int m_xyArea;

    int m_blockColor;

    frantic::graphics::vector3 m_blockCoord;

    block_vertex_manager& operator=( const block_vertex_manager& ); // not implemented

    // init
    // "nonempty" function requires ( cubeCase != 0x00 && cubeCase != 0xff )
    // I moved this out of add_interior_voxel_vertices() for the sake of inlining
    // on an old version of GCC.
    void add_interior_voxel_vertices_nonempty(
        std::size_t i, boost::uint8_t cubeCase,
        const frantic::graphics::vector3&
            voxelCoord /*, boost::int32_t nextVertexNumber, frantic::graphics::vector3 & outVertexNumbers*/ ) {
        const boost::int32_t x = voxelCoord.x;
        const boost::int32_t y = voxelCoord.y;
        const boost::int32_t z = voxelCoord.z;
        boost::uint8_t vertFlag[3];
        // get required verts
        get_corner_vert_flags( cubeCase, vertFlag );
        std::size_t newVertCount = 0;
        for( int axis = 0; axis < 3; ++axis ) {
            if( vertFlag[axis] ) {
                ++newVertCount;
            }
        }
        boost::int32_t nextVertexNumber = static_cast<boost::int32_t>( m_outMesh.vertex_count() );
        m_outMesh.add_vertices( newVertCount );
        BOOST_FOREACH( frantic::geometry::trimesh3_vertex_channel_general_accessor& vertexChannel,
                       m_outNamedVertexChannels ) {
            vertexChannel.add_vertices( newVertCount );
        }
        frantic::graphics::vector3 voxelVertexIndices( get_invalid_vertex_number() );
        if( vertFlag[0] ) {
            float density0 = m_voxelCornerDensities[i - 1], density1 = m_voxelCornerDensities[i];
            if( density0 * density1 > 0 ) {
                throw std::runtime_error( "density mismatch" );
            }
            voxelVertexIndices.x = nextVertexNumber;
            frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x - 1, y, z ),
                                       voxelCorner1 = frantic::graphics::vector3( x, y, z );
            m_isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                      nextVertexNumber++, m_outNamedVertexChannels, m_outMesh,
                                                      m_vertexWorkspace );
        }
        if( vertFlag[1] ) {
            float density0 = m_voxelCornerDensities[i - m_xsize], density1 = m_voxelCornerDensities[i];
            if( density0 * density1 > 0 ) {
                throw std::runtime_error( "density mismatch" );
            }
            voxelVertexIndices.y = nextVertexNumber;
            frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x, y - 1, z ),
                                       voxelCorner1 = frantic::graphics::vector3( x, y, z );
            m_isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                      nextVertexNumber++, m_outNamedVertexChannels, m_outMesh,
                                                      m_vertexWorkspace );
        }
        if( vertFlag[2] ) {
            float density0 = m_voxelCornerDensities[i - m_xyArea], density1 = m_voxelCornerDensities[i];
            if( density0 * density1 > 0 ) {
                throw std::runtime_error( "density mismatch" );
            }
            voxelVertexIndices.z = nextVertexNumber;
            frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x, y, z - 1 ),
                                       voxelCorner1 = frantic::graphics::vector3( x, y, z );
            m_isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                      nextVertexNumber++, m_outNamedVertexChannels, m_outMesh,
                                                      m_vertexWorkspace );
        }
        m_outVertexIndices[i] = voxelVertexIndices;
    }

    // process a voxel
    // - retrieve vertices if external
    //   - cache hash lookup
    // - expose vertices
    //
    // "nonempty" function requires ( cubeCase != 0x00 && cubeCase != 0xff )
    // I moved this out of add_interior_voxel_vertices() for the sake of inlining
    // on an old version of GCC.
    void add_boundary_voxel_vertices_nonempty(
        std::size_t i, boost::uint8_t cubeCase, int /*boundaryCase*/,
        const frantic::graphics::vector3&
            voxelCoord /*, boost::int32_t nextVertexNumber, frantic::graphics::vector3 & outVertexNumbers*/ ) {
        const boost::int32_t x = voxelCoord.x;
        const boost::int32_t y = voxelCoord.y;
        const boost::int32_t z = voxelCoord.z;
        boost::uint8_t vertFlag[3];
        // get required verts
        get_corner_vert_flags( cubeCase, vertFlag );
        boost::uint8_t vertUninitialized[3] = { vertFlag[0], vertFlag[1], vertFlag[2] };
        // try to fill external verts
        frantic::graphics::vector3 voxelVertexIndices( get_invalid_vertex_number() );
        // const int directionCode = get_boundary_owner_direction_code( m_blockColor, boundaryCase );
        // const bool fetchBoundary = true; //directionCode >= 0;
        // if( fetchBoundary ) {
        {
            // vector3 foundVertices;
            voxel_vertex_t::const_iterator iter = m_exposedVoxelVertices.find( voxelCoord );
            const bool hit = iter != m_exposedVoxelVertices.end();
            // const bool hit = m_voxelVertexCache.try_get( directionCode, voxelCoord, foundVertices );
            if( hit ) {
                const frantic::graphics::vector3 foundVertices = iter->second;
                for( int axis = 0; axis < 3; ++axis ) {
                    if( vertUninitialized[axis] && is_vertex_number( foundVertices[axis] ) ) {
                        voxelVertexIndices[axis] = detail::encode_external_vertex_number( foundVertices[axis] );
                        vertUninitialized[axis] = 0;
                    }
                }
                /*
                for( int axis = 0; axis < 3; ++axis ) {
                  if( vertUninitialized[axis] ) {
                    FF_LOG( warning ) << "no connection found" << std::endl;
                    //throw std::runtime_error( "no connection found" );
                  }
                }
                */
            }
        }
        std::size_t newVertCount = 0;
        for( int axis = 0; axis < 3; ++axis ) {
            if( vertUninitialized[axis] ) {
                ++newVertCount;
            }
        }
        boost::int32_t nextVertexNumber = static_cast<boost::int32_t>( m_outMesh.vertex_count() );
        m_outMesh.add_vertices( newVertCount );
        BOOST_FOREACH( frantic::geometry::trimesh3_vertex_channel_general_accessor& vertexChannel,
                       m_outNamedVertexChannels ) {
            vertexChannel.add_vertices( newVertCount );
        }
        // fill internal or remaining verts
        if( vertUninitialized[0] ) {
            float density0 = m_voxelCornerDensities[i - 1], density1 = m_voxelCornerDensities[i];
            if( density0 * density1 > 0 ) {
                throw std::runtime_error( "density mismatch" );
            }
            voxelVertexIndices.x = nextVertexNumber;
            frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x - 1, y, z ),
                                       voxelCorner1 = frantic::graphics::vector3( x, y, z );
            m_isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                      nextVertexNumber++, m_outNamedVertexChannels, m_outMesh,
                                                      m_vertexWorkspace );
        }
        if( vertUninitialized[1] ) {
            float density0 = m_voxelCornerDensities[i - m_xsize], density1 = m_voxelCornerDensities[i];
            if( density0 * density1 > 0 ) {
                throw std::runtime_error( "density mismatch" );
            }
            voxelVertexIndices.y = nextVertexNumber;
            frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x, y - 1, z ),
                                       voxelCorner1 = frantic::graphics::vector3( x, y, z );
            m_isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                      nextVertexNumber++, m_outNamedVertexChannels, m_outMesh,
                                                      m_vertexWorkspace );
        }
        if( vertUninitialized[2] ) {
            float density0 = m_voxelCornerDensities[i - m_xyArea], density1 = m_voxelCornerDensities[i];
            if( density0 * density1 > 0 ) {
                throw std::runtime_error( "density mismatch" );
            }
            voxelVertexIndices.z = nextVertexNumber;
            frantic::graphics::vector3 voxelCorner0 = frantic::graphics::vector3( x, y, z - 1 ),
                                       voxelCorner1 = frantic::graphics::vector3( x, y, z );
            m_isp.add_edge_isosurface_vertex_to_mesh( density0, density1, voxelCorner0, voxelCorner1,
                                                      nextVertexNumber++, m_outNamedVertexChannels, m_outMesh,
                                                      m_vertexWorkspace );
        }
        if( newVertCount > 0 ) {
            m_outVoxelVertices.push_back( exposed_voxel_vertices( voxelCoord, voxelVertexIndices ) );
        }
        /*
        if( ! fetchBoundary ) {
          for( int axis = 0; axis < 3; ++axis ) {
            if( voxelVertexIndices[axis] != detail::get_invalid_vertex_number() && detail::is_external_vertex_number(
        voxelVertexIndices[axis] ) ) { throw std::runtime_error( "publishing external vertex" );
            }
          }
          m_outExposedBlockVertices.push_back( exposed_block_vertices( m_blockCoord, voxelCoord, voxelVertexIndices ) );
        }
        */
        m_outVertexIndices[i] = voxelVertexIndices;
    }

  public:
    block_vertex_manager(
        ImplicitSurfacePolicy& isp, const frantic::volumetrics::marching_cubes_table& mct,
        const frantic::graphics::vector3& blockCoord, const frantic::graphics::boundbox3& xyzExtents,
        const int blockColor, const std::vector<float>& voxelCornerDensities,
        // const std::vector<exposed_voxel_vertices> & inVoxelVertices,
        // const stdext::hash_map<vector3,std::pair<std::size_t,std::size_t>,voxel_coord_hasher> &
        // blockCoordToVoxelVertices,
        const voxel_vertex_t& exposedVoxelVertices, std::vector<frantic::graphics::vector3>& outVertexIndices,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        frantic::geometry::trimesh3& outMesh,
        // std::vector<detail::exposed_block_vertices> & outExposedBlockVertices,
        std::vector<detail::exposed_voxel_vertices>& outExposedBlockVertices,
        typename ImplicitSurfacePolicy::vertex_workspace_t& vertexWorkspace )
        : m_isp( isp )
        , m_mct( mct )
        , m_voxelCornerDensities( voxelCornerDensities )
        ,
        // m_blockCoordToVoxelVertices( blockCoordToVoxelVertices ),
        m_outNamedVertexChannels( outNamedVertexChannels )
        , m_outMesh( outMesh )
        ,
        // m_outExposedBlockVertices( outExposedBlockVertices ),
        m_exposedVoxelVertices( exposedVoxelVertices )
        , m_vertexWorkspace( vertexWorkspace )
        ,
        // m_voxelVertexCache( blockCoord, inVoxelVertices, blockCoordToVoxelVertices ),
        m_xsize( xyzExtents.xsize() )
        , m_xyArea( xyzExtents.xsize() * xyzExtents.ysize() )
        , m_blockColor( blockColor )
        , m_blockCoord( blockCoord )
        , m_outVertexIndices( outVertexIndices )
        , m_outVoxelVertices( outExposedBlockVertices ) {}
    // init
    inline void add_interior_voxel_vertices(
        std::size_t i, boost::uint8_t cubeCase,
        const frantic::graphics::vector3&
            voxelCoord /*, boost::int32_t nextVertexNumber, frantic::graphics::vector3 & outVertexNumbers*/ ) {
        if( cubeCase != 0x00 && cubeCase != 0xff ) {
            add_interior_voxel_vertices_nonempty( i, cubeCase, voxelCoord );
        }
    }
    // process a voxel
    // - retrieve vertices if external
    //   - cache hash lookup
    // - expose vertices
    inline void add_boundary_voxel_vertices(
        std::size_t i, boost::uint8_t cubeCase, int boundaryCase,
        const frantic::graphics::vector3&
            voxelCoord /*, boost::int32_t nextVertexNumber, frantic::graphics::vector3 & outVertexNumbers*/ ) {
        if( cubeCase != 0x00 && cubeCase != 0xff ) {
            add_boundary_voxel_vertices_nonempty( i, cubeCase, boundaryCase, voxelCoord );
        }
    }
};

template <class ImplicitSurfacePolicy>
void generate_vertices_for_box(
    ImplicitSurfacePolicy& isp, const frantic::volumetrics::marching_cubes_table& mct,
    const frantic::graphics::vector3& blockCoord, const frantic::graphics::boundbox3& xyzExtents, const int blockColor,
    const std::vector<float>& voxelCornerDensities, const std::vector<boost::uint8_t>& cubeCases_,
    // const std::vector<exposed_voxel_vertices> & exposedVoxelVertices,
    // const stdext::hash_map<frantic::graphics::vector3,std::pair<size_t,size_t>,voxel_coord_hasher> &
    // blockCoordToExposedVoxelVertices,
    const voxel_vertex_t& exposedVoxelVertices, std::vector<frantic::graphics::vector3>& outVertexIndices,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
    frantic::geometry::trimesh3& outMesh, std::vector<exposed_voxel_vertices>& outExposedBoxVertices,
    typename ImplicitSurfacePolicy::vertex_workspace_t& vertexWorkspace ) {
    using std::size_t;

    const int zMin = xyzExtents.minimum().z;
    const int yMin = xyzExtents.minimum().y;
    const int xMin = xyzExtents.minimum().x;

    const int xMax = xyzExtents.maximum().x;
    const int yMax = xyzExtents.maximum().y;
    const int zMax = xyzExtents.maximum().z;

    const boost::uint8_t* cubeCases = cubeCases_.size() > 0 ? &cubeCases_[0] : 0;

    block_vertex_manager<ImplicitSurfacePolicy> boxVertexManager(
        isp, mct, blockCoord, xyzExtents, blockColor, voxelCornerDensities, exposedVoxelVertices,
        /*blockCoordToExposedVoxelVertices,*/ outVertexIndices, outNamedVertexChannels, outMesh, outExposedBoxVertices,
        vertexWorkspace );

    size_t i = 0;
    // zMin
    //   yMin
    boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 0, frantic::graphics::vector3( xMin, yMin, zMin ) );
    ++i;
    for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 1, frantic::graphics::vector3( x, yMin, zMin ) );
    }
    boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 2, frantic::graphics::vector3( xMax, yMin, zMin ) );
    ++i;
    //   y interior
    for( int32_t y = yMin + 1; y < yMax; ++y ) {
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 3, frantic::graphics::vector3( xMin, y, zMin ) );
        ++i;
        for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
            boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 4,
                                                          frantic::graphics::vector3( x, y, zMin ) );
        }
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 5, frantic::graphics::vector3( xMax, y, zMin ) );
        ++i;
    }
    //   yMax
    boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 6, frantic::graphics::vector3( xMin, yMax, zMin ) );
    ++i;
    for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 7, frantic::graphics::vector3( x, yMax, zMin ) );
    }
    boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 8, frantic::graphics::vector3( xMax, yMax, zMin ) );
    ++i;

    // z interior
    for( int32_t z = zMin + 1; z < zMax; ++z ) {
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 9, frantic::graphics::vector3( xMin, yMin, z ) );
        ++i;
        for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
            boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 10,
                                                          frantic::graphics::vector3( x, yMin, z ) );
        }
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 11,
                                                      frantic::graphics::vector3( xMax, yMin, z ) );
        ++i;
        //   y interior
        for( int32_t y = yMin + 1; y < yMax; ++y ) {
            boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 12,
                                                          frantic::graphics::vector3( xMin, y, z ) );
            ++i;
            for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
                boxVertexManager.add_interior_voxel_vertices( i, cubeCases[i], frantic::graphics::vector3( x, y, z ) );
            }
            boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 13,
                                                          frantic::graphics::vector3( xMax, y, z ) );
            ++i;
        }
        //   yMax
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 14,
                                                      frantic::graphics::vector3( xMin, yMax, z ) );
        ++i;
        for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
            boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 15,
                                                          frantic::graphics::vector3( x, yMax, z ) );
        }
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 16,
                                                      frantic::graphics::vector3( xMax, yMax, z ) );
        ++i;
    }

    // zMax
    boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 17, frantic::graphics::vector3( xMin, yMin, zMax ) );
    ++i;
    for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 18,
                                                      frantic::graphics::vector3( x, yMin, zMax ) );
    }
    boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 19, frantic::graphics::vector3( xMax, yMin, zMax ) );
    ++i;
    //   y interior
    for( int32_t y = yMin + 1; y < yMax; ++y ) {
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 20,
                                                      frantic::graphics::vector3( xMin, y, zMax ) );
        ++i;
        for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
            boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 21,
                                                          frantic::graphics::vector3( x, y, zMax ) );
        }
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 22,
                                                      frantic::graphics::vector3( xMax, y, zMax ) );
        ++i;
    }
    //   yMax
    boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 23, frantic::graphics::vector3( xMin, yMax, zMax ) );
    ++i;
    for( int32_t x = xMin + 1; x < xMax; ++x, ++i ) {
        boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 24,
                                                      frantic::graphics::vector3( x, yMax, zMax ) );
    }
    boxVertexManager.add_boundary_voxel_vertices( i, cubeCases[i], 25, frantic::graphics::vector3( xMax, yMax, zMax ) );
    ++i;
}

template< class ImplicitSurfacePolicy >
void generate_disambiguated_faces_for_box(
	ImplicitSurfacePolicy& isp,
	const frantic::volumetrics::marching_cubes_table& mct,
	const frantic::graphics::boundbox3& xyzExtents,
	const std::vector<float>& voxelCornerDensities_,
	const std::vector<boost::uint8_t>& cubeCases_,
	const std::vector<frantic::graphics::vector3>& vertexIndices_,
	const std::vector<frantic::graphics::vector3f>& externalVertices,
	std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
	frantic::geometry::trimesh3& outMesh,
	typename ImplicitSurfacePolicy::vertex_workspace_t & vertexWorkspace/*,
	frantic::diagnostics::profiling_manager & prof*/ )
{
    const int zEnd = xyzExtents.maximum().z - xyzExtents.minimum().z + 1;
    const int yEnd = xyzExtents.maximum().y - xyzExtents.minimum().y + 1;
    const int xEnd = xyzExtents.maximum().x - xyzExtents.minimum().x + 1;

    const int xyArea = xyzExtents.xsize() * xyzExtents.ysize();
    const int xsize = xyzExtents.xsize();

    const boost::uint8_t* cubeCases = cubeCases_.size() > 0 ? &cubeCases_[0] : 0;
    const float* voxelCornerDensities = voxelCornerDensities_.size() > 0 ? &voxelCornerDensities_[0] : 0;
    const frantic::graphics::vector3* vertexIndices = vertexIndices_.size() > 0 ? &vertexIndices_[0] : 0;

    std::vector<frantic::graphics::vector3> cubeFaces;

    for( int z = 1; z != zEnd; ++z ) {
        for( int y = 1; y != yEnd; ++y ) {
            std::size_t i = z * xyArea + y * xsize + 1;
            for( int x = 1; x != xEnd; ++x ) {
                const unsigned char cubeCase = cubeCases[i];
                if( cubeCase != 0x00 && cubeCase != 0xff ) {
                    // Retrieve all the vertices from the edges of this cube
                    // prof.enter_section( 6 );
                    char vertFlags[12];
                    marching_cubes_table::get_required_vert_flags( cubeCase, vertFlags );

                    boost::int32_t cubeVerts[13];
                    cubeVerts[0] =
                        vertFlags[0] ? vertexIndices[i - xyArea - xsize].x : detail::get_invalid_vertex_number();
                    cubeVerts[1] = vertFlags[1] ? vertexIndices[i - xyArea - 1].y : detail::get_invalid_vertex_number();
                    cubeVerts[2] = vertFlags[2] ? vertexIndices[i - xsize - 1].z : detail::get_invalid_vertex_number();
                    cubeVerts[3] = vertFlags[3] ? vertexIndices[i - xyArea].y : detail::get_invalid_vertex_number();
                    cubeVerts[4] = vertFlags[4] ? vertexIndices[i - xsize].z : detail::get_invalid_vertex_number();
                    cubeVerts[5] = vertFlags[5] ? vertexIndices[i - xyArea].x : detail::get_invalid_vertex_number();
                    cubeVerts[6] = vertFlags[6] ? vertexIndices[i - 1].z : detail::get_invalid_vertex_number();
                    cubeVerts[7] = vertFlags[7] ? vertexIndices[i].z : detail::get_invalid_vertex_number();
                    cubeVerts[8] = vertFlags[8] ? vertexIndices[i - xsize].x : detail::get_invalid_vertex_number();
                    cubeVerts[9] = vertFlags[9] ? vertexIndices[i - 1].y : detail::get_invalid_vertex_number();
                    cubeVerts[10] = vertFlags[10] ? vertexIndices[i].y : detail::get_invalid_vertex_number();
                    cubeVerts[11] = vertFlags[11] ? vertexIndices[i].x : detail::get_invalid_vertex_number();
                    cubeVerts[12] = detail::get_invalid_vertex_number();
                    // prof.exit_section( 6 );

                    // retrieve the corner values for disambiguation
                    // prof.enter_section( 9 );
                    float cubeCorners[8];
                    cubeCorners[0] = voxelCornerDensities[i - xyArea - xsize - 1];
                    cubeCorners[1] = voxelCornerDensities[i - xyArea - xsize];
                    cubeCorners[2] = voxelCornerDensities[i - xyArea - 1];
                    cubeCorners[3] = voxelCornerDensities[i - xyArea];
                    cubeCorners[4] = voxelCornerDensities[i - xsize - 1];
                    cubeCorners[5] = voxelCornerDensities[i - xsize];
                    cubeCorners[6] = voxelCornerDensities[i - 1];
                    cubeCorners[7] = voxelCornerDensities[i];
                    const bool newVert = mct.get_cubecase_faces( cubeCase, cubeCorners, cubeFaces );
                    // prof.exit_section( 9 );
                    // prof.enter_section( 7 );
                    if( newVert ) {
                        // tbb::spin_mutex::scoped_lock meshLock( m_meshMutex );

                        // average the existing verts for the new vert location
                        // TODO:  a policy based vertex addition would probably be more appropriate with a
                        // solve to find a vert location
                        frantic::graphics::vector3f avg;
                        int count = 0;
                        BOOST_FOREACH( boost::int32_t cubeVert, cubeVerts ) {
                            if( is_vertex_number( cubeVert ) ) {
                                if( is_internal_vertex_number( cubeVert ) ) {
                                    avg += outMesh.get_vertex( cubeVert );
                                } else {
                                    avg += externalVertices[decode_external_vertex_number( cubeVert )];
                                }
                                ++count;
                            }
                        }
                        avg *= 1.f / count;

                        const boost::int32_t vertexNumber = static_cast<boost::int32_t>( outMesh.vertex_count() );
                        cubeVerts[12] = vertexNumber;
                        outMesh.add_vertex( avg );
                        if( !outNamedVertexChannels.empty() ) {
                            BOOST_FOREACH( frantic::geometry::trimesh3_vertex_channel_general_accessor& vertexChannel,
                                           outNamedVertexChannels ) {
                                vertexChannel.add_vertex();
                            }
                            isp.populate_vertex_channels( outNamedVertexChannels, vertexNumber, avg, vertexWorkspace );
                        }
                    }
                    // prof.exit_section( 7 );
                    // prof.enter_section( 8 );

                    const frantic::graphics::vector3* const cubeFacesBegin = cubeFaces.size() > 0 ? &cubeFaces[0] : 0;
                    const frantic::graphics::vector3* const cubeFacesEnd = cubeFacesBegin + cubeFaces.size();

                    for( const frantic::graphics::vector3* iter = cubeFacesBegin; iter != cubeFacesEnd; ++iter ) {
                        // same check is run in convert_particle_implicit_surface_to_trimesh3 later
                        /*
                        if( cubeVerts[iter->x] == detail::get_invalid_vertex_number() || cubeVerts[iter->y] ==
                        detail::get_invalid_vertex_number() || cubeVerts[iter->z] == detail::get_invalid_vertex_number()
                        ) { throw std::runtime_error( "invalid vertex number" );
                        }
                        */
                        outMesh.add_face( cubeVerts[iter->x], cubeVerts[iter->y], cubeVerts[iter->z] );
                    }
                }
                ++i;
            }
        }
    }
}

inline void flush_exposed_vertex_buffer( std::vector<exposed_voxel_vertices>& exposedVoxelVerticesBuffer,
                                         boost::int32_t firstVertexNumber, tbb::spin_mutex& mutex,
                                         std::vector<exposed_voxel_vertices>& outExposedVoxelVertices ) {
    if( exposedVoxelVerticesBuffer.size() ) {
        tbb::spin_mutex::scoped_lock lock( mutex );

        BOOST_FOREACH( exposed_voxel_vertices& exposedVoxelVertices, exposedVoxelVerticesBuffer ) {
            for( int corner = 0; corner < 3; ++corner ) {
                if( is_vertex_number( exposedVoxelVertices.vertices[corner] ) ) {
                    if( is_external_vertex_number( exposedVoxelVertices.vertices[corner] ) ) {
                        exposedVoxelVertices.vertices[corner] =
                            decode_external_vertex_number( exposedVoxelVertices.vertices[corner] );
                    } else {
                        exposedVoxelVertices.vertices[corner] += firstVertexNumber;
                    }
                }
            }
        }

        std::copy( exposedVoxelVerticesBuffer.begin(), exposedVoxelVerticesBuffer.end(),
                   std::back_insert_iterator<std::vector<exposed_voxel_vertices>>( outExposedVoxelVertices ) );
        exposedVoxelVerticesBuffer.clear();
    }
}

inline boost::int32_t
flush_mesh_buffer( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& inputChannels,
                   frantic::geometry::trimesh3& buffer, tbb::spin_mutex& meshMutex, boost::int32_t& nextVertexNumber,
                   boost::int32_t& nextFaceNumber,
                   std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                   frantic::geometry::trimesh3& outMesh ) {
    if( buffer.vertex_count() || buffer.face_count() ) {
        tbb::spin_mutex::scoped_lock meshLock( meshMutex );

        const boost::int32_t firstVertexNumber = nextVertexNumber;
        nextVertexNumber += static_cast<boost::int32_t>( buffer.vertex_count() );

        const boost::int32_t firstFaceNumber = nextFaceNumber;
        nextFaceNumber += static_cast<boost::int32_t>( buffer.face_count() );

        if( firstVertexNumber + buffer.vertex_count() > outMesh.vertex_count() ) {
            const boost::int32_t newVertexCount =
                std::max<boost::int32_t>( nextVertexNumber, 2 * static_cast<boost::int32_t>( outMesh.vertex_count() ) );
            outMesh.set_vertex_count( newVertexCount );
            BOOST_FOREACH( frantic::geometry::trimesh3_vertex_channel_general_accessor& outputChannel,
                           outputChannels ) {
                outputChannel.set_vertex_count( newVertexCount );
            }
        }
        if( firstFaceNumber + buffer.face_count() > outMesh.face_count() ) {
            outMesh.set_face_count(
                std::max<boost::int32_t>( nextFaceNumber, 2 * static_cast<boost::int32_t>( outMesh.face_count() ) ) );
        }

        // meshLock.release();

        BOOST_FOREACH( frantic::graphics::vector3& face, buffer.faces_ref() ) {
            for( int corner = 0; corner < 3; ++corner ) {
                const boost::int32_t vertexNumber = face[corner];
                if( is_external_vertex_number( vertexNumber ) ) {
                    face[corner] = decode_external_vertex_number( vertexNumber );
                } else {
                    face[corner] = vertexNumber + firstVertexNumber;
                }
            }
        }

        // meshLock.acquire( meshMutex, false );

        if( buffer.vertex_count() ) {
            memcpy( &outMesh.vertices_ref()[firstVertexNumber], &buffer.vertices_ref()[0],
                    buffer.vertex_count() * sizeof( frantic::graphics::vector3f ) );
            for( std::size_t i = 0; i < outputChannels.size(); ++i ) {
                frantic::geometry::trimesh3_vertex_channel_general_accessor& outputChannel = outputChannels[i];
                frantic::geometry::trimesh3_vertex_channel_general_accessor& inputChannel = inputChannels[i];
                memcpy( outputChannel.data( firstVertexNumber ), inputChannel.data( 0 ),
                        buffer.vertex_count() * outputChannel.primitive_size() );
            }
        }
        if( buffer.face_count() ) {
            memcpy( &outMesh.faces_ref()[firstFaceNumber], &buffer.faces_ref()[0],
                    buffer.face_count() * sizeof( frantic::graphics::vector3 ) );
        }

        // meshLock.release();

        buffer.set_vertex_count( 0 );
        BOOST_FOREACH( frantic::geometry::trimesh3_vertex_channel_general_accessor& bufferChannel, inputChannels ) {
            bufferChannel.set_vertex_count( 0 );
        }
        buffer.set_face_count( 0 );

        return firstVertexNumber;
    }
    return -1;
}

} // namespace detail

template <class ImplicitSurfacePolicy>
void convert_particle_implicit_surface_to_trimesh3( ImplicitSurfacePolicy& mcp,
                                                    const frantic::channels::channel_propagation_policy& cpp,
                                                    frantic::geometry::trimesh3& outMesh,
                                                    shared_progress_logger_proxy& progressLogger ) {
    static marching_cubes_table mct;

    tbb::task_scheduler_init taskSchedulerInit;

    // Clear the mesh to start
    outMesh.clear();

    // frantic::diagnostics::profiling_manager prof;
    // prof.new_profiling_section( "total" );
    // prof.new_profiling_section( "init" );
    // prof.new_profiling_section( "geometry generation" );

    // prof.enter_section( 0 );
    // prof.enter_section( 1 );

    // Initialize the named channels we'll be creating
    std::vector<frantic::tstring> channelNames;

    std::vector<std::pair<size_t, frantic::channels::data_type_t>> channelDataTypes;

    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor> outputChannels;

    mcp.get_channel_names( cpp, channelNames );
    // This prepares the output channel general accessors contained within the mc policy, and
    // returns the primitive data types of all the channels.
    mcp.prepare_channel_accessors( channelNames, channelDataTypes );

    // Create all the channels in the output mesh, as well as the vector of writable channel accessors used
    // to fill in the data.

    for( unsigned i = 0; i < channelNames.size(); ++i ) {
        outMesh.add_vertex_channel_raw( channelNames[i], channelDataTypes[i].first, channelDataTypes[i].second );
        outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( channelNames[i] ) );
    }

    // These are the bounds we use to scan.
    frantic::graphics::boundbox3 voxelBounds = mcp.get_voxel_bounds();
    frantic::graphics::boundbox3 xyzExtents( voxelBounds.minimum().x, voxelBounds.maximum().x, voxelBounds.minimum().y,
                                             voxelBounds.maximum().y, voxelBounds.minimum().z,
                                             voxelBounds.maximum().z );

    if( xyzExtents.is_empty() )
        return;

    const frantic::graphics::size3 meshingBlockSize( 15, 15, 15 ); // meshing block length
    const frantic::graphics::size3 densityBlockSize( meshingBlockSize.xsize() + 1, meshingBlockSize.ysize() + 1,
                                                     meshingBlockSize.zsize() + 1 );
    // const int densityBlockVolume = densityBlockSize.volume();

    std::vector<frantic::graphics::vector3> meshingBlocks;

    const float meshingVoxelLength = mcp.get_voxel_coord_system().voxel_length();
    frantic::volumetrics::voxel_coord_system meshingVCS = mcp.get_voxel_coord_system();

    const frantic::graphics::size3f worldMeshingBlockSize( meshingBlockSize.xsize() * meshingVoxelLength,
                                                           meshingBlockSize.ysize() * meshingVoxelLength,
                                                           meshingBlockSize.zsize() * meshingVoxelLength );
    mcp.get_affected_blocks( worldMeshingBlockSize, meshingBlocks, progressLogger );

    if( progressLogger.is_cancelled() ) {
        outMesh.clear_and_deallocate();
        throw frantic::logging::progress_cancel_exception( "Meshing" );
    }

    // const int maxVertsPerCube = 4;
    // const int maxFacesPerCube = 36;

    // const int maxVertsPerBlock = maxVertsPerCube * densityBlockVolume;
    // const int maxFacesPerBlock = maxFacesPerCube * densityBlockVolume;

    const std::size_t bufferVertsSoftLimit = 4096;
    const std::size_t bufferFacesSoftLimit = 8192;

    // const std::size_t bufferVerts = bufferVertsSoftLimit + maxVertsPerBlock;
    // const std::size_t bufferFaces = bufferFacesSoftLimit + maxFacesPerBlock;

    // prof.exit_section( 1 );

    // prof.enter_section( 2 );

    detail::generate_block_geometry_mt( mcp, mct, meshingBlocks, meshingBlockSize, bufferVertsSoftLimit,
                                        bufferFacesSoftLimit, channelNames, outputChannels, outMesh, progressLogger );

    // prof.exit_section( 2 );

    if( progressLogger.is_cancelled() ) {
        outMesh.clear_and_deallocate();
        throw frantic::logging::progress_cancel_exception( "Meshing" );
    }

    // temporary sanity check
    for( std::size_t i = 0; i < outMesh.face_count(); ++i ) {
        const frantic::graphics::vector3 face = outMesh.get_face( i );
        for( int corner = 0; corner < 3; ++corner ) {
            if( face[corner] < 0 || static_cast<std::size_t>( face[corner] ) >= outMesh.vertex_count() ) {
                throw std::runtime_error( " vertex number out of bounds" );
            }
        }
    }

    // prof.exit_section( 0 );
    // FF_LOG( stats ) << prof << std::endl;
}

/**
 * Converts an implicit surface into a triangle mesh.  It handles the proper copying of named channels,
 * whose input is specified by the implicit surface policy.
 *
 * @param  isp              The implicit surface policy evaluats the implicit surface function for the algorithm.
 * @param  outMesh          The output parameter where the resulting mesh will go.
 *
template< class ImplicitSurfacePolicy >
void convert_sparse_implicit_surface_to_trimesh3( frantic::channels::channel_propagation_policy& cpp,
ImplicitSurfacePolicy& isp, frantic::geometry::trimesh3& outMesh)
{
  static marching_cubes_table mct;

  tbb::task_scheduler_init taskSchedulerInit;

  //sparsePM = frantic::diagnostics::profiling_manager();
  //sparsePM.new_profiling_section( "corner densities" );
  //sparsePM.new_profiling_section( "cube cases" );
  //sparsePM.new_profiling_section( "vertices" );
  //sparsePM.new_profiling_section( "faces" );
  //sparsePM.new_profiling_section( "density init" );
  //sparsePM.new_profiling_section( "density particle search" );
  //sparsePM.new_profiling_section( "density eval" );
  //sparsePM.new_profiling_section( "density end" );

  ImplicitSurfacePolicy::sparse_voxel_corner_density_workspace_t sparseVoxelCornerDensityWorkspace;

  // Clear the mesh to start
  outMesh.clear();

  // Initialize the named channels we'll be creating
  std::vector<std::string> channelNames;
  std::vector< std::pair<size_t,frantic::channels::data_type_t> > channelDataTypes;
  std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor> outputChannels;
  isp.get_channel_names( cpp, channelNames );

  // This prepares the output channel general accessors contained within the mc policy, and
  // returns the primitive data types of all the channels.
  isp.prepare_channel_accessors( channelNames, channelDataTypes );

  // Create all the channels in the output mesh, as well as the vector of writable channel accessors used
  // to fill in the data.
  for( unsigned i = 0; i < channelNames.size(); ++i ) {
    outMesh.add_vertex_channel_raw( channelNames[i], channelDataTypes[i].first, channelDataTypes[i].second );
    outputChannels.push_back( outMesh.get_vertex_channel_general_accessor(channelNames[i]) );
  }

  // These are the bounds we use to scan for corner samples
  frantic::graphics::boundbox3 voxelBounds = isp.get_voxel_bounds();
  frantic::graphics2d::boundrect2 xyExtents( voxelBounds.minimum().x, voxelBounds.maximum().x, voxelBounds.minimum().y,
voxelBounds.maximum().y ); int zMin = voxelBounds.minimum().z, zMax = voxelBounds.maximum().z;

  frantic::graphics2d::size2 densitySampleDimensions = xyExtents.size();

  // If this is empty, then let's skip all the computation below.
  if( densitySampleDimensions.get_area() <= 0 )
    return;

  // arrays for the density samples
  boost::shared_array<float> previousVoxelCornerDensities( new float[densitySampleDimensions.get_area()] ),
                 currentVoxelCornerDensities( new float[densitySampleDimensions.get_area()] );

  // Though normally the cube case dimensions should be one smaller, we're adding 1 along the x and y
  // dimensions to keep the rle plane structures and vertex index arrays consistent.
  boost::shared_array<unsigned char> previousDefinedCubeCases( new unsigned char[densitySampleDimensions.get_area()] ),
                     currentDefinedCubeCases( new unsigned char[densitySampleDimensions.get_area()] );

  // arrays of vertex indices into the trimesh
  boost::shared_array<frantic::graphics::vector3> previousVertexIndices, currentVertexIndices;
  previousVertexIndices.reset( new frantic::graphics::vector3[densitySampleDimensions.get_area()]);
  currentVertexIndices.reset( new frantic::graphics::vector3[densitySampleDimensions.get_area()]);

  // Get the voxel corner densities at the smallest Z
  rle_plane previousVoxelCornersRLP(xyExtents), currentVoxelCornersRLP(xyExtents);
  rle_plane previousDefinedCubeCasesRLP(xyExtents), currentDefinedCubeCasesRLP(xyExtents);
  isp.fill_sparse_voxel_corner_densities( xyExtents, zMin, previousVoxelCornerDensities.get(), previousVoxelCornersRLP,
sparseVoxelCornerDensityWorkspace );

  // initialize the previous vert indices for any adjacent defined corners on the initial plane
  detail::fill_sparse_initial_vertex_indices( previousVoxelCornerDensities,
                        previousVoxelCornersRLP,
                        previousVertexIndices);

  // Go through all the slices of cubes
  for( int z = zMin+1; z <= zMax; ++z ) {

    // Get the new slice of voxel density samples
    //sparsePM.enter_section(0);
    isp.fill_sparse_voxel_corner_densities( xyExtents, z, currentVoxelCornerDensities.get(), currentVoxelCornersRLP,
sparseVoxelCornerDensityWorkspace );
    //sparsePM.exit_section(0);

    // compute the cube cases that arise
    //sparsePM.enter_section(1);
    detail::fill_sparse_cube_case_values( isp.get_default_outside_distance(),
                        previousVoxelCornerDensities,
                        previousVoxelCornersRLP,
                        currentVoxelCornerDensities,
                        currentVoxelCornersRLP,
                        currentDefinedCubeCases,
                        currentDefinedCubeCasesRLP,
                        currentVertexIndices);
    //sparsePM.exit_section(1);

    // Generate the vertices
    //sparsePM.enter_section(2);
    detail::generate_vertices_for_sparse_plane( isp, z, previousVoxelCornerDensities, currentVoxelCornerDensities,
currentDefinedCubeCases, currentDefinedCubeCasesRLP,previousVertexIndices, currentVertexIndices, outputChannels, outMesh
);
    //detail::generate_vertices_for_sparse_plane_mt( isp, z, previousVoxelCornerDensities, currentVoxelCornerDensities,
currentDefinedCubeCases, currentDefinedCubeCasesRLP,previousVertexIndices, currentVertexIndices, outputChannels, outMesh
);
    //sparsePM.exit_section(2);

    // Generate the new faces
    //sparsePM.enter_section(3);
    detail::generate_faces_for_sparse_plane( mct, currentDefinedCubeCases, currentDefinedCubeCasesRLP,
previousVertexIndices, currentVertexIndices, outMesh );
    //detail::generate_disambiguated_faces_for_sparse_plane( isp, mct, previousVoxelCornerDensities,
currentVoxelCornerDensities, currentDefinedCubeCases, currentDefinedCubeCasesRLP, previousVertexIndices,
currentVertexIndices, outputChannels, outMesh );
    //sparsePM.exit_section(3);

    // Swap the previous and current values for the next iteration of the loop
    previousVertexIndices.swap( currentVertexIndices );
    previousDefinedCubeCases.swap( currentDefinedCubeCases );
    previousDefinedCubeCasesRLP.swap( currentDefinedCubeCasesRLP );
    previousVoxelCornerDensities.swap( currentVoxelCornerDensities );
    previousVoxelCornersRLP.swap(currentVoxelCornersRLP);
  }
  //FF_LOG( stats ) << sparsePM << std::endl;
}
*/

/**
 * Converts a set of particles into a triangle mesh using a union of spheres implicit
 * surface.  Ensure that the voxel length of the provided particle_grid_tree
 * is at least equal to the particle effect radius.  The conversion is done sparsely.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as particleEffectRadius.
 * @param  maximumParticleRadius This as big or larger than the radius of the largest particle in the system.
 * @param  particleRadiusToEffectRadiusScale  Multiplied by a particle radius, this gives the metaball field effect
 * radius around a particle.
 * @param  implicitThreshold     The isosurface value which defines the implicit surface.
 * @param  vertexRefinement      The number of iterations used to refine the vertex positions,
 *                               solving for the actual isosurface position.
 * @param  meshingVCS            The voxel coordinate system which specifies the grid used
 *                               for the marching cubes meshing.
 * @param  outMesh               The output parameter where the resulting mesh will go.
 */
void union_of_spheres_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                     const frantic::channels::channel_propagation_policy& cpp,
                                                     float maximumParticleRadius,
                                                     float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                     const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                     frantic::geometry::trimesh3& outMesh );

/**
 * Converts a set of particles into a triangle mesh using a union of spheres implicit
 * surface.  Ensure that the voxel length of the provided particle_grid_tree
 * is at least equal to the particle effect radius.  The conversion is done sparsely.
 * Passes a progress logger to the conversion method.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as particleEffectRadius.
 * @param  maximumParticleRadius This as big or larger than the radius of the largest particle in the system.
 * @param  particleRadiusToEffectRadiusScale  Multiplied by a particle radius, this gives the metaball field effect
 * radius around a particle.
 * @param  implicitThreshold     The isosurface value which defines the implicit surface.
 * @param  vertexRefinement      The number of iterations used to refine the vertex positions,
 *                               solving for the actual isosurface position.
 * @param  meshingVCS            The voxel coordinate system which specifies the grid used
 *                               for the marching cubes meshing.
 * @param  outMesh               The output parameter where the resulting mesh will go.
 * @param  progressLogger		 The object to accept progress on the conversion method.
 */
void union_of_spheres_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                     const frantic::channels::channel_propagation_policy& cpp,
                                                     float maximumParticleRadius,
                                                     float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                     const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                     frantic::geometry::trimesh3& outMesh,
                                                     frantic::logging::progress_logger& progressLogger );

void union_of_spheres_convert_sparse_particles_to_trimesh3(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maximumParticleRadius, float particleRadiusToEffectRadiusScale, float implicitThreshold,
    const voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh,
    frantic::logging::progress_logger& progressLogger );

void union_of_spheres_convert_sparse_particles_to_trimesh3(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maximumParticleRadius, float particleRadiusToEffectRadiusScale, float implicitThreshold,
    const voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh,
    shared_progress_logger_proxy& progressLogger );

struct union_of_spheres_sparse_params {
    frantic::particles::particle_grid_tree& particles;
    const frantic::channels::channel_propagation_policy& cpp;
    float maximumParticleRadius;
    float particleRadiusToEffectRadiusScale;
    float implicitThreshold;
    const voxel_coord_system& meshingVCS;
    int vertexRefinement;
    frantic::geometry::trimesh3& outMesh;
    shared_progress_logger_proxy& progressLogger;

    union_of_spheres_sparse_params& operator=( const union_of_spheres_sparse_params& ); // not implemented

    union_of_spheres_sparse_params( frantic::particles::particle_grid_tree& particles,
                                    const frantic::channels::channel_propagation_policy& cpp,
                                    float maximumParticleRadius, float particleRadiusToEffectRadiusScale,
                                    float implicitThreshold, const voxel_coord_system& meshingVCS, int vertexRefinement,
                                    frantic::geometry::trimesh3& outMesh,
                                    shared_progress_logger_proxy& progressLogger );
};

void union_of_spheres_convert_sparse_particles_to_trimesh3( union_of_spheres_sparse_params& params );

/**
 * Converts a set of particles into a triangle mesh using a metaballs implicit
 * surface.  Ensure that the voxel length of the provided particle_grid_tree
 * is at least equal to the particle effect radius.  The conversion is done sparsely.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as particleEffectRadius.
 * @param  maximumParticleRadius This as big or larger than the radius of the largest particle in the system.
 * @param  particleRadiusToEffectRadiusScale  Multiplied by a particle radius, this gives the metaball field effect
 * radius around a particle.
 * @param  implicitThreshold     The isosurface value which defines the implicit surface.
 * @param  vertexRefinement      The number of iterations used to refine the vertex positions,
 *                               solving for the actual isosurface position.
 * @param  meshingVCS            The voxel coordinate system which specifies the grid used
 *                               for the marching cubes meshing.
 * @param  outMesh               The output parameter where the resulting mesh will go.
 */
void metaball_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                             const frantic::channels::channel_propagation_policy& cpp,
                                             float maximumParticleRadius, float particleRadiusToEffectRadiusScale,
                                             float implicitThreshold, const voxel_coord_system& meshingVCS,
                                             int vertexRefinement, frantic::geometry::trimesh3& outMesh );

/**
 * Converts a set of particles into a triangle mesh using a metaballs implicit
 * surface.  Ensure that the voxel length of the provided particle_grid_tree
 * is at least equal to the particle effect radius.  The conversion is done sparsely.
 * Passes a progress logger to the conversion method.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as particleEffectRadius.
 * @param  maximumParticleRadius This as big or larger than the radius of the largest particle in the system.
 * @param  particleRadiusToEffectRadiusScale  Multiplied by a particle radius, this gives the metaball field effect
 * radius around a particle.
 * @param  implicitThreshold     The isosurface value which defines the implicit surface.
 * @param  vertexRefinement      The number of iterations used to refine the vertex positions,
 *                               solving for the actual isosurface position.
 * @param  meshingVCS            The voxel coordinate system which specifies the grid used
 *                               for the marching cubes meshing.
 * @param  outMesh               The output parameter where the resulting mesh will go.
 * @param  progressLogger		 The object to accept progress on the conversion method.
 */
void metaball_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                             const frantic::channels::channel_propagation_policy& cpp,
                                             float maximumParticleRadius, float particleRadiusToEffectRadiusScale,
                                             float implicitThreshold, const voxel_coord_system& meshingVCS,
                                             int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                             frantic::logging::progress_logger& progressLogger );

void metaball_convert_sparse_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                    const frantic::channels::channel_propagation_policy& cpp,
                                                    float maximumParticleRadius,
                                                    float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                    const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                    frantic::geometry::trimesh3& outMesh,
                                                    frantic::logging::progress_logger& progressLogger );

void metaball_convert_sparse_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                    const frantic::channels::channel_propagation_policy& cpp,
                                                    float maximumParticleRadius,
                                                    float particleRadiusToEffectRadiusScale, float implicitThreshold,
                                                    const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                    frantic::geometry::trimesh3& outMesh,
                                                    shared_progress_logger_proxy& progressLogger );

struct metaball_sparse_params {
    frantic::particles::particle_grid_tree& particles;
    const frantic::channels::channel_propagation_policy& cpp;
    float maximumParticleRadius;
    float particleRadiusToEffectRadiusScale;
    float implicitThreshold;
    const voxel_coord_system& meshingVCS;
    int vertexRefinement;
    frantic::geometry::trimesh3& outMesh;
    shared_progress_logger_proxy& progressLogger;

    metaball_sparse_params& operator=( const metaball_sparse_params& ); // not implemented

    metaball_sparse_params( frantic::particles::particle_grid_tree& particles,
                            const frantic::channels::channel_propagation_policy& cpp, float maximumParticleRadius,
                            float particleRadiusToEffectRadiusScale, float implicitThreshold,
                            const voxel_coord_system& meshingVCS, int vertexRefinement,
                            frantic::geometry::trimesh3& outMesh, shared_progress_logger_proxy& progressLogger );
};

void metaball_convert_sparse_particles_to_trimesh3( metaball_sparse_params& params );

/**
 * Converts a set of particles into a triangle mesh using the implicit
 * surface described in the paper by Zhu and Bridson, "Animating Sand as a Fluid" from SIGGRAPH 2005.
 * Ensure that the voxel length of the provided particle_grid_tree is at least equal to the kernel compact support.
 * The radius of the particle must be stored in a float "Radius" channel. The conversion is done sparsely.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as particleEffectRadius.
 * @param  maxParticleRadius	 The max radius of particles in the system.
 * @param  effectRadius			 The radius of effect used in conjunction with the maxParticleRadius to
 *								 determine the sphere about each particle, which is the
 *compact support of the interpolation kernel.
 * @param  lowDensityTrimmingDensity    An attempt to trim away low density regions which would still be inside the
 *surface due to the formulation.  This parameter defines the density at which trimming begins.
 * @param  lowDensityTrimmingStrength   This is a multiplier which scales the magnitude of the low density trimming.
 * @param  vertexRefinement      The number of iterations used to refine the vertex positions,
 *                               solving for the actual isosurface position.
 * @param  meshingVCS            The voxel coordinate system which specifies the grid used
 *                               for the marching cubes meshing.
 * @param  outMesh               The output parameter where the resulting mesh will go.
 */
void zhu_bridson_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float maxParticleRadius, float effectRadius,
                                                float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
                                                const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                frantic::geometry::trimesh3& outMesh );

/**
 * Converts a set of particles into a triangle mesh using the implicit
 * surface described in the paper by Zhu and Bridson, "Animating Sand as a Fluid" from SIGGRAPH 2005.
 * Ensure that the voxel length of the provided particle_grid_tree is at least equal to the kernel compact support.
 * The radius of the particle must be stored in a float "Radius" channel. The conversion is done sparsely.
 * Passes a progress logger to the conversion method.
 *
 * @param  particles             The input particles, in a particle_grid_tree.  The voxel length of this
 *                               particle_grid_tree should be at least as big as particleEffectRadius.
 * @param  maxParticleRadius	 The max radius of particles in the system.
 * @param  effectRadius			 The radius of effect used in conjunction with the maxParticleRadius to
 *								 determine the sphere about each particle, which is the
 *compact support of the interpolation kernel.
 * @param  lowDensityTrimmingDensity    An attempt to trim away low density regions which would still be inside the
 *surface due to the formulation.  This parameter defines the density at which trimming begins.
 * @param  lowDensityTrimmingStrength   This is a multiplier which scales the magnitude of the low density trimming.
 * @param  vertexRefinement      The number of iterations used to refine the vertex positions,
 *                               solving for the actual isosurface position.
 * @param  meshingVCS            The voxel coordinate system which specifies the grid used
 *                               for the marching cubes meshing.
 * @param  outMesh               The output parameter where the resulting mesh will go.
 * @param  progressLogger		 The object to accept progress on the conversion method.
 */
void zhu_bridson_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float maxParticleRadius, float effectRadius,
                                                float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
                                                const voxel_coord_system& meshingVCS, int vertexRefinement,
                                                frantic::geometry::trimesh3& outMesh,
                                                frantic::logging::progress_logger& progressLogger );

void zhu_bridson_convert_sparse_particles_to_trimesh3(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maxParticleRadius, float effectRadiusScale, float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
    const voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh,
    frantic::logging::progress_logger& progressLogger );

void zhu_bridson_convert_sparse_particles_to_trimesh3(
    frantic::particles::particle_grid_tree& particles, const frantic::channels::channel_propagation_policy& cpp,
    float maxParticleRadius, float effectRadius, float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
    const voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh,
    shared_progress_logger_proxy& progressLogger );

struct zhu_bridson_sparse_params {
    frantic::particles::particle_grid_tree& particles;
    const frantic::channels::channel_propagation_policy& cpp;
    float maximumParticleRadius;
    float particleRadiusToEffectRadiusScale;
    float lowDensityTrimmingDensity;
    float lowDensityTrimmingStrength;
    const voxel_coord_system& meshingVCS;
    int vertexRefinement;
    frantic::geometry::trimesh3& outMesh;
    shared_progress_logger_proxy& progressLogger;

    zhu_bridson_sparse_params& operator=( const zhu_bridson_sparse_params& ); // not implemented

    zhu_bridson_sparse_params( frantic::particles::particle_grid_tree& particles,
                               const frantic::channels::channel_propagation_policy& cpp, float maximumParticleRadius,
                               float particleRadiusToEffectRadiusScale, float lowDensityTrimmingDensity,
                               float lowDensityTrimmingStrength, const voxel_coord_system& meshingVCS,
                               int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                               shared_progress_logger_proxy& progressLogger );
};

void zhu_bridson_convert_sparse_particles_to_trimesh3( zhu_bridson_sparse_params& params );

void anisotropic_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                int vertexRefinement, frantic::geometry::trimesh3& outMesh );

void anisotropic_convert_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                const frantic::channels::channel_propagation_policy& cpp,
                                                float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                                frantic::logging::progress_logger& progressLogger );

void anisotropic_convert_sparse_particles_to_trimesh3( frantic::particles::particle_grid_tree& particles,
                                                       const frantic::channels::channel_propagation_policy& cpp,
                                                       float implicitThreshold, const voxel_coord_system& meshingVCS,
                                                       int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                                       frantic::logging::progress_logger& progressLogger );

struct anisotropic_sparse_params {
    frantic::particles::particle_grid_tree& particles;
    const frantic::channels::channel_propagation_policy& cpp;
    float implicitThreshold;
    const voxel_coord_system& meshingVCS;
    int vertexRefinement;
    frantic::geometry::trimesh3& outMesh;
    shared_progress_logger_proxy& progressLogger;

    anisotropic_sparse_params& operator=( const anisotropic_sparse_params& ); // not implemented

    anisotropic_sparse_params( frantic::particles::particle_grid_tree& particles,
                               const frantic::channels::channel_propagation_policy& cpp, float implicitThreshold,
                               const voxel_coord_system& meshingVCS, int vertexRefinement,
                               frantic::geometry::trimesh3& outMesh, shared_progress_logger_proxy& progressLogger );
};

void anisotropic_convert_sparse_particles_to_trimesh3( anisotropic_sparse_params& params );

/**
 * This function converts the provided level set into a triangle mesh, feeding the actual sample values
 * into marching cubes as the corner density values.  The vertex positions are linearly interpolated between
 * these sample values.  This is likely to remain the fastest level set to mesh conversion, and is suitable for
 * preview meshes.
 *
 * @param  levelSet  The input level set.
 * @param  outMesh   The output triangle mesh.
 * @param  clipBounds  how many voxels to clip the meshing volume bounds by
 */
void convert_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                   frantic::geometry::trimesh3& outMesh, int clipBounds = 1 );

/**
 * This function converts the provided level set into a triangle mesh, feeding the actual sample values
 * into marching cubes as the corner density values.  The vertex positions are linearly interpolated between
 * these sample values.  This is likely to remain the fastest level set to mesh conversion, and is suitable for
 * preview meshes.
 *
 * @param  levelSet  The input level set.
 * @param  cpp   Specifies which rle channels should be added to the mesh
 * @param  outMesh   The output triangle mesh.
 * @param  clipBounds  how many voxels to clip the meshing volume bounds by
 */
void convert_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                   frantic::channels::channel_propagation_policy& cpp,
                                   frantic::geometry::trimesh3& outMesh, int clipBounds = 1 );

/**
 * This function converts the provided level set into a triangle mesh, feeding the actual sample values
 * into marching cubes as the corner density values.  The vertex positions are linearly interpolated between
 * these sample values.  This is likely to remain the fastest level set to mesh conversion, and is suitable for
 * preview meshes.
 *
 * @param  levelSet  The input level set.
 * @param  cpp   Specifies which rle channels should be added to the mesh
 * @param  outMesh   The output triangle mesh.
 * @param  progressLogger   Receives progress updates during the conversion
 * @param  clipBounds  how many voxels to clip the meshing volume bounds by
 */
void convert_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                   frantic::channels::channel_propagation_policy& cpp,
                                   frantic::geometry::trimesh3& outMesh,
                                   frantic::logging::progress_logger& progressLogger, int clipBounds = 1 );

/*
 * This function converts the provided level set into a triangle mesh sparsely, feeding the actual sample values
 * into marching cubes as the corner density values.  The vertex positions are linearly interpolated between
 * these sample values.  This is likely to remain the fastest level set to mesh conversion, and is suitable for
 * preview meshes.
 *
 * @param  levelSet  The input level set.
 * @param  outMesh   The output triangle mesh.
 * @param  clipBounds  how many voxels to clip the meshing volume bounds by
 */
// void convert_levelset_to_trimesh3_sparse( const frantic::volumetrics::levelset::rle_level_set& levelSet,
// frantic::geometry::trimesh3& outMesh, int clipBounds = 1 );

/*
 * This function converts the provided level set into a triangle mesh sparsely, feeding the actual sample values
 * into marching cubes as the corner density values.  The vertex positions are linearly interpolated between
 * these sample values.  This is likely to remain the fastest level set to mesh conversion, and is suitable for
 * preview meshes.
 *
 * @param  levelSet  The input level set.
 * @param  cpp   Specifies which rle channels should be added to the mesh
 * @param  outMesh   The output triangle mesh.
 * @param  clipBounds  how many voxels to clip the meshing volume bounds by
 */
// void convert_levelset_to_trimesh3_sparse( const frantic::volumetrics::levelset::rle_level_set& levelSet,
// frantic::channels::channel_propagation_policy& cpp, frantic::geometry::trimesh3& outMesh, int clipBounds = 1 );

template <class FilterType>
void convert_filtered_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                            frantic::channels::channel_propagation_policy& cpp,
                                            const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                            int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                            frantic::logging::progress_logger& progressLogger, int clipBounds ) {
    // Create the policy, and call the implicit surface conversion function
    reconstruction_filtered_rle_level_set_is_policy<FilterType> mcp( levelSet, meshingVCS, vertexRefinement,
                                                                     FilterType(), clipBounds );
    frantic::volumetrics::implicitsurface::convert_implicit_surface_to_trimesh3( cpp, mcp, outMesh, progressLogger );
}

template <class FilterType>
void convert_filtered_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                            frantic::channels::channel_propagation_policy& cpp,
                                            const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                            int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                            int clipBounds ) {
    frantic::logging::null_progress_logger nullLogger;
    convert_filtered_levelset_to_trimesh3<FilterType>( levelSet, cpp, meshingVCS, vertexRefinement, outMesh, nullLogger,
                                                       clipBounds );
}

template <class FilterType>
void convert_filtered_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                            const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                            int vertexRefinement, frantic::geometry::trimesh3& outMesh,
                                            int clipBounds ) {
    frantic::channels::channel_propagation_policy cpp;
    convert_filtered_levelset_to_trimesh3<FilterType>( levelSet, cpp, meshingVCS, vertexRefinement, outMesh,
                                                       clipBounds );
}

/**
 * This function converts the provided level set into a triangle mesh, using a tricubic filter for reconstructing
 * the function value.  Currently it uses a bspline filter for the reconstruction.
 *
 * @todo  This function should accept a parameter to specify what reconstruction filter to use.
 *
 * @param  levelSet  The input level set.
 * @param  meshingVCS  The voxel coord system for the meshing volume
 * @param  vertexRefinement  The number of vertex refinements to be performed by the solver.
 * @param  outMesh   The output triangle mesh.
 * @param  clipBounds  how many voxels to clip the meshing volume bounds by
 */
void convert_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                   const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement,
                                   frantic::geometry::trimesh3& outMesh, int clipBounds = 1 );

// with channel propagation policy
void convert_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                   frantic::channels::channel_propagation_policy& cpp,
                                   const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement,
                                   frantic::geometry::trimesh3& outMesh, int clipBounds = 1 );

// with channel propagation policy and progress logger
void convert_levelset_to_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                   frantic::channels::channel_propagation_policy& cpp,
                                   const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement,
                                   frantic::geometry::trimesh3& outMesh,
                                   frantic::logging::progress_logger& progressLogger, int clipBounds = 1 );

/**
 * This function converts the provided level set into a triangle mesh sparsely, using a tricubic filter for
reconstructing
 * the function value.  Currently it uses a bspline filter for the reconstruction.
 *
 * @todo  This is a little buggy and needs to be looked at before use.
 * @todo  This function should accept a parameter to specify what reconstruction filter to use.
 *
 * @param  levelSet  The input level set.
 * @param  meshingVCS  The voxel coord system for the meshing volume
 * @param  vertexRefinement  The number of vertex refinements to be performed by the solver.
 * @param  outMesh   The output triangle mesh.
 * @param  clipBounds  Whether or not to clip the bounds of the meshing volume (e.g. to remove the box that would
otherwise get tacked onto the bottom of an open liquid surface) void convert_levelset_to_trimesh3_sparse( const
frantic::volumetrics::levelset::rle_level_set& levelSet, frantic::channels::channel_propagation_policy& cpp, const
frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement, frantic::geometry::trimesh3& outMesh, int
clipBounds = 1);
 */

/**
 *  Fix topology issues that may exist in a mesh output by our marching
 * cubes implementation.
 *
 *  In some cases, our marching cubes implementation outputs a mesh
 * that is not a closed manifold.  This function fixes the mesh topology
 * so that it becomes a closed manifold.
 *
 *  The specific feature that we correct is a double-sided quad connecting
 * to surfaces together -- imaging a flattened drinking straw connecting
 * two inflated balloons.
 *
 * @param mesh the mesh to process.  The procedure assumes that the mesh was
 *    produced by one of the
 *    frantic::volumetrics::implicitsurface::*_to_trimesh3() functions --
 *    it is not suitable for correcting general topology issues.
 * @param progressLogger the progress logger to use while processing the mesh.
 */
void fix_marching_cubes_topology( frantic::geometry::trimesh3& mesh,
                                  frantic::logging::progress_logger& progressLogger );

/**
 * This function creates a channel_map that is optimized for Frost and SIMD calculations.
 * Channels are ordered starting with "Position" and "Radius" and followed in descending order of
 * size. All channels of the same type will be grouped together.
 *
 * @param  channelMap  The list of channels the use to create the optimized channel_map.
 * @param  alignmentPadding  Ensure that the final structure of the channel_map is alignmed
 *    to this many bytes, default value is 4.
 */
frantic::channels::channel_map create_optimized_channel_map( const frantic::channels::channel_map& channelMap,
                                                             size_t alignmentPadding = 4 );

} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * This function converts the provided level set into a debugging mesh, based on the rle_index_spec and the level
 * set signed distance function.   It uses the "Color" channel of the output mesh to mark parts of the mesh.
 *
 * @param  levelSet  The input level set.
 * @param  cubeSize  The size of debugging cubes to make.
 * @param  outMesh   The output triangle mesh.
 */
void convert_levelset_to_ris_debug_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                             float cubeSize, frantic::geometry::trimesh3& outMesh );

void convert_levelset_to_ris_debug_trimesh3( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                             float cubeSize, frantic::geometry::trimesh3& outMesh,
                                             frantic::logging::progress_logger& progressLogger );

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
