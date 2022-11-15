// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/mem_fn.hpp>

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_named_channels.hpp>
#include <frantic/math/reconstruction_filters.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

#pragma warning( push )
#pragma warning( disable : 4512 4100 4245 )
#include <tbb/spin_rw_mutex.h>
#pragma warning( pop )

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

namespace detail {
struct empty_workspace {};
} // namespace detail

/**
 * This is a policy class which provides the implicit surface by directly copying the level set distance function
 * values. Also known as "Preview Trilinear" in the GUI.
 *
 */
class direct_linear_rle_level_set_is_policy {
  public:
    typedef detail::empty_workspace sparse_voxel_corner_density_workspace_t;
    typedef detail::empty_workspace vertex_workspace_t;

  private:
    const frantic::volumetrics::levelset::rle_level_set& m_levelSet;
    int m_boundsClip;

    std::vector<frantic::volumetrics::levelset::const_rle_channel_general_accessor> m_inputChannels;

    // When an isosurface location is found, these parameters specify where it is.  In other implementations of the is
    // policy, the array of filter weights and data indices would be stored to reuse for all the different named
    // channels.
    frantic::graphics::vector3f m_isosurfaceLocationWorld;
    // isosurfaceLocationWorld = (1-alpha) * corner_sample_coord_to_world(voxelCorner0) + alpha *
    // corner_sample_coord_to_world(voxelCorner1);
    frantic::graphics::vector3 m_isosurfaceLocationVoxelCorner0, m_isosurfaceLocationVoxelCorner1;
    float m_isosurfaceLocationAlpha;

    void populate_vertex_channels(
        const std::vector<frantic::volumetrics::levelset::const_rle_channel_general_accessor>& inputChannels,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
        const frantic::graphics::vector3f& vert ) const;

    // Disable copy constructor
    direct_linear_rle_level_set_is_policy& operator=( const direct_linear_rle_level_set_is_policy& ); // not implemented
  public:
    // function docs in cpp file
    direct_linear_rle_level_set_is_policy( const frantic::volumetrics::levelset::rle_level_set& ls,
                                           int boundsClip = 1 );
    float get_default_outside_distance() const { return m_levelSet.get_outside_distance(); };
    const frantic::volumetrics::voxel_coord_system& get_voxel_coord_system() const;
    frantic::graphics::boundbox3 get_voxel_bounds() const;
    frantic::graphics::vector3f
    corner_sample_coord_to_world( const frantic::graphics::vector3& voxelCornerCoord ) const;
    void find_edge_isosurface_location( float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
                                        const frantic::graphics::vector3& voxelCorner1 );
    frantic::graphics::vector3f find_internal_isosurface_location( float density0, float density1,
                                                                   const frantic::graphics::vector3& voxelCorner0,
                                                                   const frantic::graphics::vector3& voxelCorner1 );
    const frantic::graphics::vector3f& get_isosurface_location_world() const;
    void get_channel_names( std::vector<frantic::tstring>& outNames ) const;
    void get_channel_names( const frantic::channels::channel_propagation_policy& cpp,
                            std::vector<frantic::tstring>& outNames ) const;
    void get_channel_info( frantic::channels::channel_propagation_policy cpp, std::vector<frantic::tstring>& outNames,
                           std::vector<frantic::channels::data_type_t>& outChannelType,
                           std::vector<size_t>& outChannelArity, std::vector<size_t>& outChannelPrimitiveSize ) const;
    void prepare_channel_accessors( const std::vector<frantic::tstring>& names,
                                    std::vector<std::pair<std::size_t, frantic::channels::data_type_t>>& outDataTypes );
    void set_isosurface_location_world( const frantic::graphics::vector3f& /*worldLocation*/ ){}; // not yet implemented
    void
    add_vertex_to_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels );
    void get_interface_widths( float& interfaceWidthInside, float& interfaceWidthOutside ) const;
    int get_exterior_region_code() const;

    void fill_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                      std::vector<float>& outVoxelCornerValues ) const;
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, frantic::volumetrics::rle_plane& rlp );
    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, detail::empty_workspace& );
    void fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                   std::vector<frantic::tstring>& channelNames, std::vector<char*>& channelData,
                                   frantic::volumetrics::rle_plane& rlp );

    void add_edge_isosurface_vertex_to_mesh(
        float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
        const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        frantic::geometry::trimesh3& outMesh, detail::empty_workspace& ) const;
    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 const frantic::channels::channel_propagation_policy& cpp ) const;
    void
    populate_vertex_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              vertex_workspace_t& workspace ) const;
};

namespace detail {
struct reconstruction_filtered_rle_level_set_vertex_workspace {
    std::vector<boost::int32_t> inputDataIndexes;
    std::vector<float> collapsedSamples;
    std::vector<const char*> inputDataPointers;
    std::vector<int> namedChannelReconstructionIndex;
    std::vector<float> namedChannelReconstructionCoeff;
};
} // namespace detail

/**
 * This is a policy class which provides the implicit surface using a 4x4x4 reconstruction filter.
 *
 */
template <class ReconstructionFilter>
class reconstruction_filtered_rle_level_set_is_policy {
  public:
    typedef detail::empty_workspace sparse_voxel_corner_density_workspace_t;
    typedef detail::reconstruction_filtered_rle_level_set_vertex_workspace vertex_workspace_t;

  private:
    const frantic::volumetrics::levelset::rle_level_set& m_levelSet;

    // How much to refine the vertex positions
    int m_vertexRefinement;

    // The voxel grid used to sample the implicit surface.
    frantic::volumetrics::voxel_coord_system m_meshingVCS;
    // The bounding box to mesh
    frantic::graphics::boundbox3 m_meshingVCSBounds;

    std::vector<frantic::volumetrics::levelset::const_rle_channel_general_accessor> m_inputChannels;

    // This is an array of indexes and coefficients to be used for reconstructing the values for the
    // named channels
    std::vector<int> m_namedChannelReconstructionIndex;
    std::vector<float> m_namedChannelReconstructionCoeff;
    // This is the world location currently being processed on the isosurface
    frantic::graphics::vector3f m_isosurfaceLocationWorld;

    ReconstructionFilter m_reconFilter;

    // Disable copy constructor
    reconstruction_filtered_rle_level_set_is_policy&
    operator=( const reconstruction_filtered_rle_level_set_is_policy& ) {}

    void find_edge_isosurface_location_impl(
        float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
        const frantic::graphics::vector3& voxelCorner1, frantic::graphics::vector3f& outIsosurfaceLocationWorld,
        std::vector<int>& outNamedChannelReconstructionIndex, std::vector<float>& outNamedChannelReconstructionCoeff,
        detail::reconstruction_filtered_rle_level_set_vertex_workspace& workspace ) const {
        using namespace frantic::graphics;

        const frantic::graphics::vector3f lsVoxelCoord0 =
            m_levelSet.get_voxel_coord_system().get_voxel_coord( corner_sample_coord_to_world( voxelCorner0 ) ) -
            frantic::graphics::vector3f( 0.5f );
        const frantic::graphics::vector3f lsVoxelCoord1 =
            m_levelSet.get_voxel_coord_system().get_voxel_coord( corner_sample_coord_to_world( voxelCorner1 ) ) -
            frantic::graphics::vector3f( 0.5f );

        boundbox3 inputSampleBounds;
        // Compute the input sample bounds, by including the right 4x4x4 grid around both voxelCorner0 and voxelCorner1.
        vector3 minCoord( (int)floorf( lsVoxelCoord0.x - 1.f ), (int)floorf( lsVoxelCoord0.y - 1.f ),
                          (int)floorf( lsVoxelCoord0.z - 1.f ) );
        inputSampleBounds += minCoord;
        inputSampleBounds += minCoord + vector3( 3, 3, 3 );
        minCoord.set( (int)floorf( lsVoxelCoord1.x - 1.f ), (int)floorf( lsVoxelCoord1.y - 1.f ),
                      (int)floorf( lsVoxelCoord1.z - 1.f ) );
        inputSampleBounds += minCoord;
        inputSampleBounds += minCoord + vector3( 3, 3, 3 );

        size3 inputSampleBoundsSize = inputSampleBounds.size();

        std::vector<boost::int32_t>& inputDataIndexes = workspace.inputDataIndexes;
        inputDataIndexes.resize( inputSampleBounds.get_volume() );
        m_levelSet.get_rle_index_spec().fill_data_index_box( inputSampleBounds, &inputDataIndexes[0] );

        /////////////////
        // Collapse the samples into one line of samples along the root-finding axis
        /////////////////

        float perpFilterPlane[16];
        float collapsedFloatCoord0, collapsedFloatCoord1;
        std::vector<float>& collapsedSamples = workspace.collapsedSamples;

        int axis;
        if( voxelCorner0.x != voxelCorner1.x ) {
            axis = 0;

            float yFloat = lsVoxelCoord0.y - inputSampleBounds.yminimum(),
                  zFloat = lsVoxelCoord0.z - inputSampleBounds.zminimum();

            // Create the YZ filter plane which will be used to collapse the samples to a line in the X direction,
            // and to build the 64 reconstruction weights for named channel reconstruction.
            float yFilterLine[4];
            for( int y = 0; y < 4; ++y )
                yFilterLine[y] = (float)m_reconFilter( yFloat - y );
            int index = 0;
            for( int z = 0; z < 4; ++z ) {
                float zFilterCoeff = (float)m_reconFilter( zFloat - z );
                for( int y = 0; y < 4; ++y, ++index ) {
                    perpFilterPlane[index] = yFilterLine[y] * zFilterCoeff;
                }
            }

            // Figure out the coordinate extents for the iterative solve.
            collapsedFloatCoord0 = lsVoxelCoord0.x - inputSampleBounds.xminimum();
            collapsedFloatCoord1 = lsVoxelCoord1.x - inputSampleBounds.xminimum();

            // Collapse the YZ plane
            collapsedSamples.resize( inputSampleBoundsSize.xsize() );
            for( size_t x = 0; x < collapsedSamples.size(); ++x ) {
                float sample = 0;
                int index = 0;
                for( int z = 0; z < 4; ++z ) {
                    for( int y = 0; y < 4; ++y, ++index ) {
                        sample += perpFilterPlane[index] *
                                  m_levelSet.get_using_data_index(
                                      inputDataIndexes[x + inputSampleBoundsSize.xsize() *
                                                               ( y + inputSampleBoundsSize.ysize() * z )] );
                    }
                }
                collapsedSamples[x] = sample;
            }

            // fout << "\n";
            // for( size_t x = 0; x < collapsedSamples.size(); ++x )
            //	fout << collapsedSamples[x] << " ";
            // fout << "\n";
            // fout << "lsVoxelCoord0 == " << lsVoxelCoord0 << ", lsVoxelCoord1 == " << lsVoxelCoord0 << "\n";
            // fout << "density0 == " << density0 << " at coordinate " << collapsedFloatCoord0 << "\n";
            // fout << "density1 == " << density1 << " at coordinate " << collapsedFloatCoord1 << "\n";
            // fout << "\n";

            // float density = 0;
            // int collapsedIndexBegin = max( (int)floorf(collapsedFloatCoord0 - 1.f), 0 );
            // int collapsedIndexEnd = min( collapsedIndexBegin + 4, (int)collapsedSamples.size() );
            // for( int i = collapsedIndexBegin; i < collapsedIndexEnd; ++i ) {
            //	fout << collapsedSamples[i] << " * " << (float)m_reconFilter(collapsedFloatCoord0 - i);
            //	if( i < collapsedIndexEnd-1 )
            //		fout << " + ";
            //	density += collapsedSamples[i] * (float)m_reconFilter(collapsedFloatCoord0 - i);
            // }
            // fout << " = " << density << "\n";

            // density = 0;
            // collapsedIndexBegin = max( (int)floorf(collapsedFloatCoord1 - 1.f), 0 );
            // collapsedIndexEnd = min( collapsedIndexBegin + 4, (int)collapsedSamples.size() );
            // for( int i = collapsedIndexBegin; i < collapsedIndexEnd; ++i ) {
            //	fout << collapsedSamples[i] << " * " << (float)m_reconFilter(collapsedFloatCoord1 - i);
            //	if( i < collapsedIndexEnd-1 )
            //		fout << " + ";
            //	density += collapsedSamples[i] * (float)m_reconFilter(collapsedFloatCoord1 - i);
            // }
            // fout << " = " << density << "\n";

        } else if( voxelCorner0.y != voxelCorner1.y ) {
            axis = 1;

            float xFloat = lsVoxelCoord0.x - inputSampleBounds.xminimum(),
                  zFloat = lsVoxelCoord0.z - inputSampleBounds.zminimum();

            // Create the XZ filter plane which will be used to collapse the samples to a line in the Y direction,
            // and to build the 64 reconstruction weights for named channel reconstruction.
            float xFilterLine[4];
            for( int x = 0; x < 4; ++x )
                xFilterLine[x] = (float)m_reconFilter( xFloat - x );
            int index = 0;
            for( int z = 0; z < 4; ++z ) {
                float zFilterCoeff = (float)m_reconFilter( zFloat - z );
                for( int x = 0; x < 4; ++x, ++index ) {
                    perpFilterPlane[index] = xFilterLine[x] * zFilterCoeff;
                }
            }

            // Figure out the coordinate extents for the iterative solve.
            collapsedFloatCoord0 = lsVoxelCoord0.y - inputSampleBounds.yminimum();
            collapsedFloatCoord1 = lsVoxelCoord1.y - inputSampleBounds.yminimum();

            // Collapse the XZ plane
            collapsedSamples.resize( inputSampleBoundsSize.ysize() );
            for( size_t y = 0; y < collapsedSamples.size(); ++y ) {
                float sample = 0;
                int index = 0;
                for( int z = 0; z < 4; ++z ) {
                    for( int x = 0; x < 4; ++x, ++index ) {
                        sample += perpFilterPlane[index] *
                                  m_levelSet.get_using_data_index(
                                      inputDataIndexes[x + inputSampleBoundsSize.xsize() *
                                                               ( y + inputSampleBoundsSize.ysize() * z )] );
                    }
                }
                collapsedSamples[y] = sample;
            }
        }
        // else if( voxelCorner0.z != voxelCorner1.z )
        else {
            axis = 2;

            float xFloat = lsVoxelCoord0.x - inputSampleBounds.xminimum(),
                  yFloat = lsVoxelCoord0.y - inputSampleBounds.yminimum();

            // Create the XY filter plane which will be used to collapse the samples to a line in the Z direction,
            // and to build the 64 reconstruction weights for named channel reconstruction.
            float xFilterLine[4];
            for( int x = 0; x < 4; ++x )
                xFilterLine[x] = (float)m_reconFilter( xFloat - x );
            int index = 0;
            for( int y = 0; y < 4; ++y ) {
                float yFilterCoeff = (float)m_reconFilter( yFloat - y );
                for( int x = 0; x < 4; ++x, ++index ) {
                    perpFilterPlane[index] = xFilterLine[x] * yFilterCoeff;
                }
            }

            // Figure out the coordinate extents for the iterative solve.
            collapsedFloatCoord0 = lsVoxelCoord0.z - inputSampleBounds.zminimum();
            collapsedFloatCoord1 = lsVoxelCoord1.z - inputSampleBounds.zminimum();

            // Collapse the XY plane
            collapsedSamples.resize( inputSampleBoundsSize.zsize() );
            for( size_t z = 0; z < collapsedSamples.size(); ++z ) {
                float sample = 0;
                int index = 0;
                for( int y = 0; y < 4; ++y ) {
                    for( int x = 0; x < 4; ++x, ++index ) {
                        sample += perpFilterPlane[index] *
                                  m_levelSet.get_using_data_index(
                                      inputDataIndexes[x + inputSampleBoundsSize.xsize() *
                                                               ( y + inputSampleBoundsSize.ysize() * z )] );
                    }
                }
                collapsedSamples[z] = sample;
            }
        }

        /////////////////
        // Find the root in the collapsed samples
        /////////////////

        // Make sure coordinate 0 has the negative value (this is assumed in the iterative solve below)
        if( density0 > density1 ) {
            std::swap( density0, density1 );
            std::swap( collapsedFloatCoord0, collapsedFloatCoord1 );
        }

        float alpha = 0;
        for( int vr = 0; vr < m_vertexRefinement; ++vr ) {

            // Special case when we hit the isosurface exactly
            if( density1 == 0 )
                break;

            // Compute the isosurface intersection, assuming a linear function
            alpha = density0 / ( density0 - density1 );
            float collapsedFloatCoord = ( 1 - alpha ) * collapsedFloatCoord0 + alpha * collapsedFloatCoord1;

            // evaluate the density at this location
            float density = 0;
            int collapsedIndexBegin = std::max( (int)floorf( collapsedFloatCoord - 1.f ), 0 );
            int collapsedIndexEnd = std::min( collapsedIndexBegin + 4, (int)collapsedSamples.size() );
            for( int i = collapsedIndexBegin; i < collapsedIndexEnd; ++i ) {
                density += collapsedSamples[i] * (float)m_reconFilter( collapsedFloatCoord - i );
            }

            // Replace either the negative or the positive coordinate with the new refined density
            if( density < 0 ) {
                density0 = density;
                collapsedFloatCoord0 = collapsedFloatCoord;
            } else {
                density1 = density;
                collapsedFloatCoord1 = collapsedFloatCoord;
            }
        }

        // Compute the final refined position
        alpha = density0 / ( density0 - density1 );
        float collapsedFloatCoord = ( 1 - alpha ) * collapsedFloatCoord0 + alpha * collapsedFloatCoord1;

        // Build the actual isosurface coordinate in world space, and set up the reconstruction indexes
        // and coefficients for the named channels.
        outNamedChannelReconstructionIndex.clear();
        outNamedChannelReconstructionCoeff.clear();
        switch( axis ) {
        case 0: {
            outIsosurfaceLocationWorld = m_levelSet.get_voxel_coord_system().get_world_coord(
                frantic::graphics::vector3f( collapsedFloatCoord + inputSampleBounds.xminimum() + 0.5f,
                                             lsVoxelCoord0.y + 0.5f, lsVoxelCoord0.z + 0.5f ) );

            // Figure out the range to use along the search axis
            int collapsedIndexBegin = std::max( (int)floorf( collapsedFloatCoord - 1.f ), 0 );
            int collapsedIndexEnd = std::min( collapsedIndexBegin + 4, (int)collapsedSamples.size() );
            float collapsedAxisReconFilter[4];
            for( int i = collapsedIndexBegin; i < collapsedIndexEnd; ++i )
                collapsedAxisReconFilter[i - collapsedIndexBegin] = (float)m_reconFilter( collapsedFloatCoord - i );

            // Go through all the example points, and produce the dataIndex/coefficient pairs.  Also sum
            // up all the coefficients so we can normalize them
            float coefficientSum = 0;
            for( int z = 0; z < 4; ++z ) {
                for( int y = 0; y < 4; ++y ) {
                    for( int x = collapsedIndexBegin; x < collapsedIndexEnd; ++x ) {
                        int dataIndex = inputDataIndexes[x + inputSampleBoundsSize.xsize() *
                                                                 ( y + inputSampleBoundsSize.ysize() * z )];
                        if( dataIndex >= 0 ) {
                            float coefficient =
                                collapsedAxisReconFilter[x - collapsedIndexBegin] * perpFilterPlane[y + 4 * z];
                            if( coefficient != 0 ) {
                                outNamedChannelReconstructionIndex.push_back( dataIndex );
                                outNamedChannelReconstructionCoeff.push_back( coefficient );
                                /*
                                if( !_finite(coefficient) ) {
                                  ofstream fout( "c:\\maxtests\\levelsetmesher\\dbg.txt", ios::app|ios::out );
                                  fout << "Bad coefficient " << coefficient << "\n";
                                  fout << "collapsedAxisReconFilter[x-collapsedIndexBegin]: " <<
                                collapsedAxisReconFilter[x-collapsedIndexBegin] << "\n"; fout << "perpFilterPlane[y + 4
                                * z]: " << perpFilterPlane[y + 4 * z] << "\n";
                                }
                                */
                                coefficientSum += coefficient;
                            }
                        }
                    }
                }
            }
            // Normalize the coefficients
            if( coefficientSum != 0 ) {
                for( unsigned i = 0; i < outNamedChannelReconstructionCoeff.size(); ++i )
                    outNamedChannelReconstructionCoeff[i] /= coefficientSum;
            }

            break;
        }
        case 1: {
            outIsosurfaceLocationWorld =
                m_levelSet.get_voxel_coord_system().get_world_coord( frantic::graphics::vector3f(
                    lsVoxelCoord0.x + 0.5f, collapsedFloatCoord + inputSampleBounds.yminimum() + 0.5f,
                    lsVoxelCoord0.z + 0.5f ) );

            // Figure out the range to use along the search axis
            int collapsedIndexBegin = std::max( (int)floorf( collapsedFloatCoord - 1.f ), 0 );
            int collapsedIndexEnd = std::min( collapsedIndexBegin + 4, (int)collapsedSamples.size() );
            float collapsedAxisReconFilter[4];
            for( int i = collapsedIndexBegin; i < collapsedIndexEnd; ++i )
                collapsedAxisReconFilter[i - collapsedIndexBegin] = (float)m_reconFilter( collapsedFloatCoord - i );

            // Go through all the xample points, and produce the dataIndex/coefficient pairs.  Also sum
            // up all the coefficients so we can normalize them
            float coefficientSum = 0;
            for( int z = 0; z < 4; ++z ) {
                for( int y = collapsedIndexBegin; y < collapsedIndexEnd; ++y ) {
                    for( int x = 0; x < 4; ++x ) {
                        int dataIndex = inputDataIndexes[x + inputSampleBoundsSize.xsize() *
                                                                 ( y + inputSampleBoundsSize.ysize() * z )];
                        if( dataIndex >= 0 ) {
                            float coefficient =
                                collapsedAxisReconFilter[y - collapsedIndexBegin] * perpFilterPlane[x + 4 * z];
                            if( coefficient != 0 ) {
                                outNamedChannelReconstructionIndex.push_back( dataIndex );
                                outNamedChannelReconstructionCoeff.push_back( coefficient );
                                coefficientSum += coefficient;
                            }
                        }
                    }
                }
            }
            // Normalize the coefficients
            if( coefficientSum != 0 ) {
                for( unsigned i = 0; i < outNamedChannelReconstructionCoeff.size(); ++i )
                    outNamedChannelReconstructionCoeff[i] /= coefficientSum;
            }

            break;
        }
        case 2: {
            outIsosurfaceLocationWorld = m_levelSet.get_voxel_coord_system().get_world_coord(
                frantic::graphics::vector3f( lsVoxelCoord0.x + 0.5f, lsVoxelCoord0.y + 0.5f,
                                             collapsedFloatCoord + inputSampleBounds.zminimum() + 0.5f ) );

            // Figure out the range to use along the search axis
            int collapsedIndexBegin = std::max( (int)floorf( collapsedFloatCoord - 1.f ), 0 );
            int collapsedIndexEnd = std::min( collapsedIndexBegin + 4, (int)collapsedSamples.size() );
            float collapsedAxisReconFilter[4];
            for( int i = collapsedIndexBegin; i < collapsedIndexEnd; ++i )
                collapsedAxisReconFilter[i - collapsedIndexBegin] = (float)m_reconFilter( collapsedFloatCoord - i );

            // Go through all the example points, and produce the dataIndex/coefficient pairs.  Also sum
            // up all the coefficients so we can normalize them
            float coefficientSum = 0;
            for( int z = collapsedIndexBegin; z < collapsedIndexEnd; ++z ) {
                for( int y = 0; y < 4; ++y ) {
                    for( int x = 0; x < 4; ++x ) {
                        int dataIndex = inputDataIndexes[x + inputSampleBoundsSize.xsize() *
                                                                 ( y + inputSampleBoundsSize.ysize() * z )];
                        if( dataIndex >= 0 ) {
                            float coefficient =
                                collapsedAxisReconFilter[z - collapsedIndexBegin] * perpFilterPlane[x + 4 * y];
                            if( coefficient != 0 ) {
                                outNamedChannelReconstructionIndex.push_back( dataIndex );
                                outNamedChannelReconstructionCoeff.push_back( coefficient );
                                coefficientSum += coefficient;
                            }
                        }
                    }
                }
            }
            // Normalize the coefficients
            if( coefficientSum != 0 ) {
                for( unsigned i = 0; i < outNamedChannelReconstructionCoeff.size(); ++i )
                    outNamedChannelReconstructionCoeff[i] /= coefficientSum;
            }

            break;
        }
        default:
            throw std::runtime_error(
                "reconstruction_filtered_rle_level_set_is_policy.find_edge_isosurface_location_impl: "
                "Internal error, unexpected axis." );
        }
    }

    // copied from add_vertex_to_channels()
    void set_vertex_channels( std::vector<geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              std::size_t vertexNumber, const std::vector<int>& namedChannelReconstructionIndex,
                              const std::vector<float>& namedChannelReconstructionCoeff,
                              detail::reconstruction_filtered_rle_level_set_vertex_workspace& workspace ) const {
        if( namedChannelReconstructionIndex.empty() ) {
            // If we got no coefficients, set all the values to 0
            for( size_t i = 0, ie = outputChannels.size(); i != ie; ++i )
                memset( outputChannels[i].data( vertexNumber ), 0, outputChannels[i].primitive_size() );
        } else {
            std::vector<const char*>& dataPointers = workspace.inputDataPointers;
            dataPointers.resize( namedChannelReconstructionIndex.size() );
            for( size_t i = 0, ie = outputChannels.size(); i != ie; ++i ) {
                const frantic::volumetrics::levelset::const_rle_channel_general_accessor& inputChannel =
                    m_inputChannels[i];
                frantic::geometry::trimesh3_vertex_channel_general_accessor& outputChannel = outputChannels[i];
                // Initialize all the input data pointers
                // for( size_t j = 0, je = dataPointers.size(); j != je; ++j )
                // dataPointers[j] = inputChannel.data(namedChannelReconstructionIndex[j]);
                std::transform( namedChannelReconstructionIndex.begin(), namedChannelReconstructionIndex.end(),
                                dataPointers.begin(),
                                std::bind( &frantic::volumetrics::levelset::const_rle_channel_general_accessor::data,
                                           &inputChannel, std::placeholders::_1 ) );
                // Call the weighted sum function to reconstruct the value
                inputChannel.get_weighted_sum_combine_function()(
                    &namedChannelReconstructionCoeff[0], &dataPointers[0], namedChannelReconstructionCoeff.size(),
                    outputChannel.arity(), outputChannel.data( vertexNumber ) );

                /*
                // DEBUGGING
                if( m_inputChannels[i].data_type() == channels::data_type_float32 &&
                  !_finite(*reinterpret_cast<float*>(outputChannels[i].data(outputChannels[i].size()-1))) )
                {
                  ofstream fout( "c:\\temp\\dbg_cubic.txt", ios::app|ios::out );
                  fout << "Vertex " << outputChannels[i].size()-1 << " is bad\n";
                  fout << "Float weights:\n";
                  for( size_t j = 0, je = dataPointers.size(); j != je; ++j )
                    fout << m_namedChannelReconstructionCoeff[j] << " ";
                  fout << "\n";
                  fout << "Data indexes:\n";
                  for( size_t j = 0, je = dataPointers.size(); j != je; ++j )
                    fout << m_namedChannelReconstructionIndex[j] << " ";
                  fout << "\n";
                }
                //*/
            }
        }
    }

    /**
     *	This function populates all vertex channels using data from the implicit
     *	surface policy for the input vertex.
     *
     *	@param	inputChannels	Accessors to particle channels used for population.
     *	@param	outputChannels	Accessors to each channel that should be populated.
     *	@param	vert			Location of the vertex.
     *	@param	vertIndex		Index of the vertex.
     */
    void populate_vertex_channels(
        const std::vector<frantic::volumetrics::levelset::const_rle_channel_general_accessor>& inputChannels,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, std::size_t vertIndex,
        const frantic::graphics::vector3f& vert ) const {
        if( inputChannels.size() != outputChannels.size() )
            throw std::runtime_error( "reconstruction_filtered_rle_level_set_is_policy.populate_vertex_channels: "
                                      "inputChannels and outputChannels are not the same size" );

        frantic::graphics::vector3f voxelCoordinate( m_levelSet.get_voxel_coord_system().get_voxel_coord( vert ) );
        voxelCoordinate -= frantic::graphics::vector3f( 0.5f, 0.5f, 0.5f );
        frantic::graphics::vector3 minCoordinate( (int)floor( voxelCoordinate.x - 1.f ),
                                                  (int)floor( voxelCoordinate.y - 1.f ),
                                                  (int)floor( voxelCoordinate.z - 1.f ) );
        frantic::graphics::vector3 maxCoordinate( minCoordinate.x + 3, minCoordinate.y + 3, minCoordinate.z + 3 );
        frantic::graphics::boundbox3 sampleBounds( minCoordinate, maxCoordinate );
        std::vector<boost::int32_t> dataIndex( 64 );
        std::size_t dataIndexSize = 0;

        float totalWeights = 0.f;
        std::vector<float> dx( 4 );
        std::vector<float> dy( 4 );
        std::vector<float> dz( 4 );
        std::vector<float> multipliers( 64 );
        std::vector<const char*> dataPointers( 64 );

        dx[0] = (float)m_reconFilter( voxelCoordinate.x - minCoordinate.x );
        dx[1] = (float)m_reconFilter( voxelCoordinate.x - minCoordinate.x - 1 );
        dx[2] = (float)m_reconFilter( maxCoordinate.x - voxelCoordinate.x - 1 );
        dx[3] = (float)m_reconFilter( maxCoordinate.x - voxelCoordinate.x );

        dy[0] = (float)m_reconFilter( voxelCoordinate.y - minCoordinate.y );
        dy[1] = (float)m_reconFilter( voxelCoordinate.y - minCoordinate.y - 1 );
        dy[2] = (float)m_reconFilter( maxCoordinate.y - voxelCoordinate.y - 1 );
        dy[3] = (float)m_reconFilter( maxCoordinate.y - voxelCoordinate.y );

        dz[0] = (float)m_reconFilter( voxelCoordinate.z - minCoordinate.z );
        dz[1] = (float)m_reconFilter( voxelCoordinate.z - minCoordinate.z - 1 );
        dz[2] = (float)m_reconFilter( maxCoordinate.z - voxelCoordinate.z - 1 );
        dz[3] = (float)m_reconFilter( maxCoordinate.z - voxelCoordinate.z );

        m_levelSet.get_rle_index_spec().fill_data_index_box( sampleBounds, &dataIndex[0] );

        for( std::size_t z = 0; z < 4; z++ )
            for( std::size_t y = 0; y < 4; y++ )
                for( std::size_t x = 0; x < 4; x++ )
                    if( dataIndex[16 * z + 4 * y + x] >= 0 ) {
                        multipliers[dataIndexSize] = ( dx[x] * dy[y] * dz[z] );
                        totalWeights += multipliers[dataIndexSize];
                        dataIndex[dataIndexSize++] = dataIndex[16 * z + 4 * y + x];
                    }

        // normalize the weights
        if( totalWeights < 1 && totalWeights != 0 )
            for( std::size_t i = 0; i < dataIndexSize; i++ )
                multipliers[i] /= totalWeights;

        for( std::size_t i = 0; i < inputChannels.size(); ++i ) {
            for( std::size_t j = 0; j < dataIndexSize; j++ ) {
                dataPointers[j] = inputChannels[i].data( dataIndex[j] );
            }
            inputChannels[i].get_weighted_sum_combine_function()( &multipliers[0], &dataPointers[0], dataIndexSize,
                                                                  inputChannels[i].arity(),
                                                                  outputChannels[i].data( vertIndex ) );
        }
    }

  public:
    void add_edge_isosurface_vertex_to_mesh(
        float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
        const frantic::graphics::vector3& voxelCorner1, std::size_t vertexNumber,
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outNamedVertexChannels,
        frantic::geometry::trimesh3& outMesh,
        detail::reconstruction_filtered_rle_level_set_vertex_workspace& workspace ) const {
        frantic::graphics::vector3f isosurfaceLocationWorld;
        find_edge_isosurface_location_impl( density0, density1, voxelCorner0, voxelCorner1, isosurfaceLocationWorld,
                                            workspace.namedChannelReconstructionIndex,
                                            workspace.namedChannelReconstructionCoeff, workspace );

        { // scope for mesh lock
            // add a vert to the mesh at the surface location
            // tbb::spin_rw_mutex::scoped_lock meshLock( meshMutex, false );
            outMesh.vertices_ref()[vertexNumber] = isosurfaceLocationWorld;
            set_vertex_channels( outNamedVertexChannels, vertexNumber, workspace.namedChannelReconstructionIndex,
                                 workspace.namedChannelReconstructionCoeff, workspace );
        }
    }

    reconstruction_filtered_rle_level_set_is_policy( const frantic::volumetrics::levelset::rle_level_set& levelSet,
                                                     const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                     int vertexRefinement, const ReconstructionFilter& reconFilter,
                                                     int boundsClip = 1 )
        : m_levelSet( levelSet )
        , m_meshingVCS( meshingVCS )
        , m_vertexRefinement( vertexRefinement )
        , m_reconFilter( reconFilter ) {
        frantic::graphics::boundbox3 voxelBounds = m_levelSet.get_rle_index_spec().outer_bounds();
        if( !boundsClip )
            voxelBounds.expand( 1 );

        m_meshingVCSBounds =
            m_meshingVCS.get_voxel_bounds( m_levelSet.get_voxel_coord_system().get_world_bounds( voxelBounds ) );
        m_meshingVCSBounds = frantic::graphics::boundbox3(
            m_meshingVCSBounds.minimum(), m_meshingVCSBounds.maximum() - frantic::graphics::vector3( 1, 1, 1 ) );

        if( !boundsClip )
            m_meshingVCSBounds.expand( 1 );
        else
            m_meshingVCSBounds.expand( -boundsClip + 1 );
    }

    const voxel_coord_system& get_voxel_coord_system() const { return m_meshingVCS; }

    // This returns the XY bounds within which to operate
    frantic::graphics::boundbox3 get_voxel_bounds() const { return m_meshingVCSBounds; }

    float get_default_outside_distance() const { return m_levelSet.get_outside_distance(); }

    void fill_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                      std::vector<float>& outVoxelCornerValues ) const {
        if( xyExtents.get_area() != outVoxelCornerValues.size() ) {
            throw std::runtime_error(
                "reconstruction_filtered_rle_level_set_is_policy.fill_voxel_corner_densities Error: extents area (" +
                boost::lexical_cast<std::string>( xyExtents.get_area() ) + ") does not match the array size (" +
                boost::lexical_cast<std::string>( outVoxelCornerValues.size() ) + ")" );
        }
        if( outVoxelCornerValues.size() > 0 ) {
#ifndef FRANTIC_DISABLE_THREADS
            m_levelSet.fill_plane_mt( m_meshingVCS, m_reconFilter, xyExtents, z, &outVoxelCornerValues[0] );
#else
#pragma message( "Threads are disabled" )
            m_levelSet.fill_plane( m_meshingVCS, m_reconFilter, xyExtents, z, outVoxelCornerValues );
#endif
        }
    }

    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, detail::empty_workspace& ) {
        m_levelSet.fill_plane_mt( m_meshingVCS, m_reconFilter, xyExtents, z, outVoxelCornerValues );
    }

    frantic::graphics::vector3f
    corner_sample_coord_to_world( const frantic::graphics::vector3& voxelCornerCoord ) const {
        return m_meshingVCS.get_world_voxel_center( voxelCornerCoord );
    }

    /**
     * This function does a solve to find the actual isosurface location between two voxel centers.  It
     * is assumed that the signs of density0 and density1 are different, and that voxelCorner0
     * has a smaller coordinate along one of X, Y or Z.
     *
     * @param  density0  The evaluated density at voxelCorner0.
     * @param  density1  The evaluated density at voxelCorner1.
     * @param  voxelCorner0  The first voxel coordinate.
     * @param  voxelCorner1  The second voxel coordinate, which is equal to voxelCorner0 except for one coordinate which
     *                       is one greater.
     */
    void find_edge_isosurface_location( float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
                                        const frantic::graphics::vector3& voxelCorner1 ) {
        detail::reconstruction_filtered_rle_level_set_vertex_workspace workspace;
        find_edge_isosurface_location_impl( density0, density1, voxelCorner0, voxelCorner1, m_isosurfaceLocationWorld,
                                            m_namedChannelReconstructionIndex, m_namedChannelReconstructionCoeff,
                                            workspace );
    }

    /**
     * This function computes an additional vertex location internal to a voxel given two
     * opposing voxel corners.  It doesn't modify the internal state of the is policy.
     *
     * @param  density0  The density function at voxelCorner0.
     * @param  density1  The density function at voxelCorner1.
     * @param  voxelCorner0  The first of two corners between which the isosurface is contained.
     * @param  voxelCorner1  The second of two corners between which the isosurface is contained.
     * @return vector3f	 The 3D location of the vertex on the isosurface internal to the voxel.
     *
     * TODO: This is just trilinear at the moment.  Need the filtered version of this solve.
     */
    frantic::graphics::vector3f find_internal_isosurface_location( float density0, float density1,
                                                                   const frantic::graphics::vector3& voxelCorner0,
                                                                   const frantic::graphics::vector3& voxelCorner1 ) {
        using namespace frantic::graphics;

        float densityDelta = density1 - density0;
        float isosurfaceLocationAlpha;

        if( densityDelta == 0 )
            isosurfaceLocationAlpha = 0.5f;
        else
            isosurfaceLocationAlpha = fabsf( density0 / densityDelta );
        vector3 isosurfaceLocationVoxelCorner0 = voxelCorner0;
        vector3 isosurfaceLocationVoxelCorner1 = voxelCorner1;
        return ( 1 - isosurfaceLocationAlpha ) * corner_sample_coord_to_world( voxelCorner0 ) +
               isosurfaceLocationAlpha * corner_sample_coord_to_world( voxelCorner1 );
    }

    const frantic::graphics::vector3f& get_isosurface_location_world() const { return m_isosurfaceLocationWorld; }

    void set_isosurface_location_world( const frantic::graphics::vector3f& /*worldLocation*/ ){}; // not yet implemented

    ////////
    // Functions for dealing with the named channels
    ////////

    /**
     * This retrieves the named channels in the input data structure.
     *
     * @param  outNames  A vector of strings into which the channel names are added.
     */
    void get_channel_names( std::vector<frantic::tstring>& outNames ) const {
        m_levelSet.get_channel_names( outNames );
    }

    /**
     *	This retrieves the named channels in the level set that the policy is operating on.
     *	It adds the names to the existing vector, so if you want just the names this function provides,
     *	make sure to initialize the output vector to empty yourself.  It uses the channel propagation
     *	policy to determine which channels should be returned.
     *
     *	@param	outNames	A vector of channel name strings.
     */
    void get_channel_names( const frantic::channels::channel_propagation_policy& cpp,
                            std::vector<frantic::tstring>& outNames ) const {
        std::vector<frantic::tstring> channelNames;
        m_levelSet.get_channel_names( channelNames );
        size_t channelCount = channelNames.size();

        // if ( cpp.is_channel_included("SignedDistance") ) {
        //	outNames.push_back("SignedDistance");
        // }

        for( size_t i = 0; i < channelCount; ++i ) {
            if( cpp.is_channel_included( channelNames[i] ) ) {
                outNames.push_back( channelNames[i] );
            }
        }
    }

    /**
     * This retrieves all the channel information associated with the named channels
     * in isp.
     *
     * @param	cpp						A channel propagation policy.
     * @param	outName					A vector of strings into which the channel names are added.
     * @param	outChannelType			A vector of data_type_t for each channel
     * @param	outChannelArity			A vector of size_t arities for each channel
     * @param	outChannelPrimitiveSize	A vector of size_t primitive sizes for each channel.
     */
    void get_channel_info( frantic::channels::channel_propagation_policy cpp, std::vector<frantic::tstring>& outName,
                           std::vector<frantic::channels::data_type_t>& outChannelType,
                           std::vector<size_t>& outChannelArity, std::vector<size_t>& outChannelPrimitiveSize ) const {

        std::vector<std::string> channelNames;
        m_levelSet.get_channel_names( channelNames );

        if( cpp.is_channel_included( _T("SignedDistance") ) ) {
            outName.push_back( _T("SignedDistance") );
            outChannelType.push_back( frantic::channels::data_type_float32 );
            outChannelArity.push_back( 1 );
            outChannelPrimitiveSize.push_back( 4 );
        }

        for( size_t i = 0; i < channelNames.size(); ++i ) {

            if( cpp.is_channel_included( channelNames[i] ) ) {
                frantic::volumetrics::levelset::const_rle_channel_general_accessor acc =
                    m_levelSet.get_channel_general_accessor( channelNames[i] );
                outName.push_back( channelNames[i] );
                outChannelType.push_back( acc.data_type() );
                outChannelArity.push_back( acc.arity() );
                outChannelPrimitiveSize.push_back( acc.primitive_size() );
            }
        }
    }

    /**
     * This function is called so that the channel accessors can be prepared,
     * and to report the data types for the algorithm to create corresponding
     * channels in the output data structure.
     *
     * @param  names  The channel names for which to prepare accessors.
     * @param  outDataTypes  The data types that the source data structure has for the corresponding names.
     */
    void
    prepare_channel_accessors( const std::vector<frantic::tstring>& names,
                               std::vector<std::pair<std::size_t, frantic::channels::data_type_t>>& outDataTypes ) {
        outDataTypes.clear();
        outDataTypes.reserve( names.size() );
        m_inputChannels.clear();
        m_inputChannels.reserve( names.size() );
        for( size_t i = 0, ie = names.size(); i != ie; ++i ) {
            m_inputChannels.push_back( m_levelSet.get_channel_general_accessor( names[i] ) );
            outDataTypes.push_back(
                std::make_pair( m_inputChannels.back().arity(), m_inputChannels.back().data_type() ) );
        }
    }

    /**
     * Given a set of writable trimesh3 accessors, this should add one data element to each of the channels,
     * which corresponds the the channel names given in prepare_channel_accessors.  The function uses the internal state
     * of the is policy for the location at which to look up the data
     *
     * @param  outputChannels  An array of trimesh3 vertex channel accessors.  One data element gets added to each.
     */
    void add_vertex_to_channels( std::vector<geometry::trimesh3_vertex_channel_general_accessor>& outputChannels ) {
        if( m_namedChannelReconstructionIndex.empty() ) {
            // If we got no coefficients, set all the values to 0
            for( size_t i = 0, ie = outputChannels.size(); i != ie; ++i )
                memset( outputChannels[i].add_vertex(), 0, outputChannels[i].primitive_size() );
        } else {
            std::vector<const char*> dataPointers( m_namedChannelReconstructionIndex.size() );
            for( size_t i = 0, ie = outputChannels.size(); i != ie; ++i ) {
                // Initialize all the input data pointers
                for( size_t j = 0, je = dataPointers.size(); j != je; ++j )
                    dataPointers[j] = m_inputChannels[i].data( m_namedChannelReconstructionIndex[j] );
                // Call the weighted sum function to reconstruct the value
                m_inputChannels[i].get_weighted_sum_combine_function()(
                    &m_namedChannelReconstructionCoeff[0], &dataPointers[0], m_namedChannelReconstructionCoeff.size(),
                    outputChannels[i].arity(), outputChannels[i].add_vertex() );

                /*
                // DEBUGGING
                if( m_inputChannels[i].data_type() == channels::data_type_float32 &&
                  !_finite(*reinterpret_cast<float*>(outputChannels[i].data(outputChannels[i].size()-1))) )
                {
                  ofstream fout( "c:\\temp\\dbg_cubic.txt", ios::app|ios::out );
                  fout << "Vertex " << outputChannels[i].size()-1 << " is bad\n";
                  fout << "Float weights:\n";
                  for( size_t j = 0, je = dataPointers.size(); j != je; ++j )
                    fout << m_namedChannelReconstructionCoeff[j] << " ";
                  fout << "\n";
                  fout << "Data indexes:\n";
                  for( size_t j = 0, je = dataPointers.size(); j != je; ++j )
                    fout << m_namedChannelReconstructionIndex[j] << " ";
                  fout << "\n";
                }
                //*/
            }
        }
    }

    void get_interface_widths( float& interfaceWidthInside, float& interfaceWidthOutside ) const {
        interfaceWidthInside = m_levelSet.get_interface_voxel_width_inside();
        interfaceWidthOutside = m_levelSet.get_interface_voxel_width_outside();
    }

    int get_exterior_region_code() const { return m_levelSet.get_rle_index_spec().get_exterior_region_code(); }

    void fill_sparse_voxel_corner_densities( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                             float* outVoxelCornerValues, frantic::volumetrics::rle_plane& rlp ) {
        std::vector<std::string> channelNames;
        channelNames.push_back( "SignedDistance" );
        std::vector<char*> channelData;
        channelData.push_back( (char*)outVoxelCornerValues );
        m_levelSet.fill_sparse_plane_channel_data( m_meshingVCS, frantic::math::b_spline_filter(), xyExtents, z,
                                                   channelNames, channelData, rlp );
        // m_levelSet.fill_sparse_plane(m_meshingVCS, frantic::math::b_spline_filter(), xyExtents, z,
        // outVoxelCornerValues, rlp);
    }

    void fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                   std::vector<frantic::tstring>& channelNames, std::vector<char*>& channelData,
                                   frantic::volumetrics::rle_plane& rlp ) {

        m_levelSet.fill_sparse_plane_channel_data( m_meshingVCS, frantic::math::b_spline_filter(), xyExtents, z,
                                                   channelNames, channelData, rlp );
    }

    /**
     *	This function populates all vertex channels specified in the channel propagation
     *	policy using info from the surface policy.
     *	<p>
     *	All other channel info is preserved.
     *
     *	@param	outMesh		Mesh containing vertex and channel info.
     *	@param	cpp			A channel propagation policy.
     */
    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 const frantic::channels::channel_propagation_policy& cpp ) const {
        std::vector<frantic::tstring> previousNames; // all channels names in the mesh
        std::vector<frantic::tstring> names;         // all channels that are going to be populated
        std::vector<std::pair<std::size_t, frantic::channels::data_type_t>>
            dataTypes; // types and arity of each channel to be populated
        std::vector<frantic::volumetrics::levelset::const_rle_channel_general_accessor>
            inputChannels; // channel accessors to channel_map
        std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>
            outputChannels; // channel accessors to vertex_channel

        // If channels already exist throw error.
        outMesh.get_vertex_channel_names( previousNames );
        for( std::size_t i = 0; i < previousNames.size(); ++i )
            if( cpp.is_channel_included( previousNames[i] ) )
                throw std::runtime_error(
                    "reconstruction_filtered_rle_level_set_is_policy.populate_mesh_channels: Attempt to "
                    "populate a channel that already exists." );

        // Find the names of all channels found in both cpp and the particles.
        get_channel_names( cpp, names );

        // Prepare the input channel accessors and get the data types which
        // will be used to create new channels in the mesh.
        dataTypes.reserve( names.size() );
        inputChannels.reserve( names.size() );
        for( std::size_t i = 0; i < names.size(); ++i ) {
            inputChannels.push_back( m_levelSet.get_channel_general_accessor( names[i] ) );
            dataTypes.push_back( std::make_pair( inputChannels.back().arity(), inputChannels.back().data_type() ) );
        }

        // Create the new vertex and get accessors to those vertex channels.
        for( std::size_t i = 0; i < names.size(); ++i ) {
            outMesh.add_vertex_channel_raw( names[i], dataTypes[i].first, dataTypes[i].second );
            outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( names[i] ) );
        }

        // Populate each vertex channel.
        for( std::size_t i = 0; i < outMesh.vertex_count(); ++i ) {
            populate_vertex_channels( inputChannels, outputChannels, i, outMesh.get_vertex( i ) );
        }
    }
    void
    populate_vertex_channels( std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
                              size_t vertIndex, const frantic::graphics::vector3f& vert,
                              vertex_workspace_t& /*workspace*/ ) const {
        populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, vert );
    }
};

} // namespace implicitsurface
} // namespace volumetrics
}; // namespace frantic
