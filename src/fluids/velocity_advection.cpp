// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/fluids/velocity_advection.hpp>
#include <frantic/volumetrics/levelset/rle_defined_and_adj_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_defined_box_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>
#include <frantic/volumetrics/rle_weno_interpolation.hpp>
/*
// for debug
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_file_io.hpp>
#include <frantic/volumetrics/rle_weno_interpolation.hpp>
#include <frantic/volumetrics/voxel_field_visualization.hpp>
*/

#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

namespace frantic {
namespace fluids {

using namespace std;
using namespace frantic::graphics;
using namespace frantic::logging;
using namespace frantic::volumetrics;
using namespace frantic::volumetrics::levelset;

struct profiling_group {
    frantic::diagnostics::profiling_section psAdvectionMain;
    frantic::diagnostics::profiling_section psAdvectionNoIter;
    frantic::diagnostics::profiling_section psTraceRK;
    frantic::diagnostics::profiling_section psLerp;
    frantic::diagnostics::profiling_section psInit;
    frantic::diagnostics::profiling_section psLookups;
    frantic::diagnostics::profiling_section psGetIndices;
};

profiling_group velAdvect;
profiling_group oldVelAdvect;

void ready_profiling() {
    velAdvect.psAdvectionMain.set_name( _T("Main Advection Loop") );
    velAdvect.psAdvectionNoIter.set_name( _T("Main Advection w/o Iterator") );
    velAdvect.psTraceRK.set_name( _T("RK Trace") );
    velAdvect.psLerp.set_name( _T("Trilinear Total") );
    velAdvect.psInit.set_name( _T("Trilinear Init") );
    velAdvect.psLookups.set_name( _T("Trilinear Lookup") );
    velAdvect.psGetIndices.set_name( _T("Trilinear Get Indices") );

    oldVelAdvect.psAdvectionMain.set_name( _T("Old Main Advection Loop") );
    oldVelAdvect.psAdvectionNoIter.set_name( _T("Old Main Advection w/o Iterator") );
    oldVelAdvect.psTraceRK.set_name( _T("Old RK Trace") );
    oldVelAdvect.psLerp.set_name( _T("Old Trilinear Total") );
    oldVelAdvect.psInit.set_name( _T("Old Trilinear Init") );
    oldVelAdvect.psLookups.set_name( _T("Old Trilinear Lookup") );
    oldVelAdvect.psGetIndices.set_name( _T("Old Trilinear Get Indices") );
}

void dump_profiling( std::basic_ostream<frantic::tchar>& out ) {
    out << velAdvect.psGetIndices << std::endl;
    out << velAdvect.psLookups << std::endl;
    out << velAdvect.psInit << std::endl;
    out << velAdvect.psLerp << std::endl;
    out << velAdvect.psTraceRK << std::endl;
    out << velAdvect.psAdvectionNoIter << std::endl;
    out << velAdvect.psAdvectionMain << std::endl;
    out << _T("\n\t****\n\n");
    out << oldVelAdvect.psGetIndices << std::endl;
    out << oldVelAdvect.psLookups << std::endl;
    out << oldVelAdvect.psInit << std::endl;
    out << oldVelAdvect.psLerp << std::endl;
    out << oldVelAdvect.psTraceRK << std::endl;
    out << oldVelAdvect.psAdvectionNoIter << std::endl;
    out << oldVelAdvect.psAdvectionMain << std::endl;
    out << _T("\n\t****\n\n");
}

namespace detail {

float get_maximum_voxel_motion_from_staggered_velocity( const rle_voxel_field& field, const float dt ) {
    const_rle_channel_accessor<vector3f> velAcc = field.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    float maxVelocityMagnitude = 0;
    for( size_t i = 0; i < velAcc.size(); ++i ) {
        for( int axis = 0; axis < 3; ++axis ) {
            maxVelocityMagnitude = std::max<float>( maxVelocityMagnitude, fabsf( velAcc[i][axis] ) );
        }
    }

    return maxVelocityMagnitude * fabsf( dt ) / field.get_voxel_coord_system().voxel_length();
}

} // namespace detail

/**
 *  Perform velocity advection using RLE lookups to sample the velocity field.
 * rle_scanline_window_advect is faster for lesser advection motion.
 */
class rle_scanline_lookup_advect {
    const boundbox3& m_bounds;
    const rle_index_spec &m_velRIS, m_srcRIS;
    const voxel_coord_system &m_srcVcs, m_velVcs;

    const_rle_channel_accessor<vector3f>& m_velocityAcc;
    const_rle_channel_accessor<vector3f>& m_srcVelAcc;
    rle_channel_accessor<vector3f>& m_resultVelAcc;

    float m_dt;

    rle_scanline_lookup_advect& operator=( const rle_scanline_lookup_advect& ); // not implemented

  public:
    rle_scanline_lookup_advect( const boundbox3& bounds, const rle_index_spec& velRIS, const rle_index_spec& srcRIS,
                                const voxel_coord_system& srcVcs, const voxel_coord_system& velVcs,
                                const_rle_channel_accessor<vector3f>& velocityAcc,
                                const_rle_channel_accessor<vector3f>& srcVelAcc,
                                rle_channel_accessor<vector3f>& resultVelAcc, float dt, bool isSameRleIndexSpecs )
        : m_bounds( bounds )
        , m_velRIS( velRIS )
        , m_srcRIS( srcRIS )
        , m_srcVcs( srcVcs )
        , m_velVcs( velVcs )
        , m_velocityAcc( velocityAcc )
        , m_srcVelAcc( srcVelAcc )
        , m_resultVelAcc( resultVelAcc )
        , m_dt( dt ) {
        if( !isSameRleIndexSpecs ) {
            throw std::runtime_error( "rle_scanline_lookup_advect Error: rle_index_specs are not equal!" );
        }
    }

    void operator()( const tbb::blocked_range2d<size_t>& r ) const {

        const vector3& boundsMin = m_bounds.minimum();

        cached_trilerp_vector3f_staggered velLookup( m_velRIS );

        for( size_t c = r.rows().begin(); c != r.rows().end(); ++c ) {
            int z = static_cast<int>( c ) + boundsMin.z;
            for( size_t b = r.cols().begin(); b != r.cols().end(); ++b ) {
                int y = static_cast<int>( b ) + boundsMin.y;

                for( rle_run_iterator i( m_velRIS, static_cast<int>( b ), static_cast<int>( c ) ), ie; i != ie; ++i ) {
                    boost::int32_t dataIndex = i.get_data_index();
                    if( dataIndex < 0 ) {
                        continue;
                    }
                    for( boost::int32_t x = i.get_xmin(); x <= i.get_xmax(); ++x, ++dataIndex ) {
                        const vector3 cellVoxelCoord( x, y, z );
                        vector3f newVelocity;

                        // get the current lookup velocity
                        for( int comp = 0; comp < 3; ++comp ) {
                            // shift to the real world location of the cell sample
                            // we are using only the x,y,z faces so we need to double the index
                            const vector3f faceVoxelCoord = m_velVcs.get_voxel_face_center( cellVoxelCoord, comp * 2 );
                            const vector3f faceWorldCoord = m_velVcs.get_world_face_center( cellVoxelCoord, comp * 2 );

                            vector3f vel;
                            // TODO: don't silently ignore failed lookups
                            velLookup.get( m_velocityAcc, faceVoxelCoord, vel );

                            const vector3f midpoint = faceWorldCoord - 0.5f * vel * m_dt;
                            vector3f vel2;
                            // TODO: don't silently ignore failed lookups
                            velLookup.get( m_velocityAcc, m_velVcs.get_voxel_coord( midpoint ), vel2 );

                            const vector3f quarterpoint = faceWorldCoord - 0.75f * vel2 * m_dt;
                            vector3f vel3;
                            // TODO: don't silently ignore failed lookups
                            velLookup.get( m_velocityAcc, m_velVcs.get_voxel_coord( quarterpoint ), vel3 );

                            const vector3f finalPosition =
                                faceWorldCoord - ( 2 * vel + 3 * vel2 + 4 * vel3 ) * m_dt / 9.f;
                            vector3f finalPosVelocity;
                            // TODO: don't silently ignore failed lookups
                            velLookup.get( m_srcVelAcc, m_srcVcs.get_voxel_coord( finalPosition ), finalPosVelocity );
                            newVelocity[comp] = finalPosVelocity[comp];
                        }
                        m_resultVelAcc[dataIndex] = newVelocity;
                    }
                }
            }
        }
    }
};

/**
 *  Perform velocity advection using a sliding window of data indices to
 * sample the velocity field.
 */
class rle_scanline_window_advect {
    const boundbox3& m_bounds;
    const size3& m_boxSize;
    const rle_index_spec &m_velRIS, m_srcRIS;
    const voxel_coord_system &m_srcVcs, m_velVcs; //, m_resultVcs;
    vector3 m_centerOffset;
    int cellBoxIndex;
    bool m_sameRIS;

    const_rle_channel_accessor<vector3f>& m_velocityAcc;
    const_rle_channel_accessor<vector3f>& m_srcVelAcc;
    rle_channel_accessor<vector3f>& m_resultVelAcc;

    float m_dt;

    rle_scanline_window_advect& operator=( const rle_scanline_window_advect& ); // not implemented

  public:
    rle_scanline_window_advect( const boundbox3& bounds, const size3& boxSize, const rle_index_spec& velRIS,
                                const rle_index_spec& srcRIS, const voxel_coord_system& srcVcs,
                                const voxel_coord_system& velVcs, const_rle_channel_accessor<vector3f>& velocityAcc,
                                const_rle_channel_accessor<vector3f>& srcVelAcc,
                                rle_channel_accessor<vector3f>& resultVelAcc, float dt, bool isSameRleIndexSpecs )
        : m_bounds( bounds )
        , m_boxSize( boxSize )
        , m_velRIS( velRIS )
        , m_srcRIS( srcRIS )
        , m_srcVcs( srcVcs )
        , m_velVcs( velVcs )
        , m_sameRIS( isSameRleIndexSpecs )
        , m_velocityAcc( velocityAcc )
        , m_srcVelAcc( srcVelAcc )
        , m_resultVelAcc( resultVelAcc )
        , m_dt( dt ) {
        m_centerOffset = vector3( m_boxSize.xsize() / 2 );
        cellBoxIndex = m_centerOffset.x + m_centerOffset.y * m_boxSize.xsize() +
                       m_centerOffset.z * m_boxSize.xsize() * m_boxSize.ysize();
    }

    static size3 get_index_box_size( float maxVoxelMotion ) {
        const int sideLength = 3 + 2 * static_cast<int>( max<float>( 0, ceilf( maxVoxelMotion - 0.49f ) ) );
        return size3( sideLength );
    }

    void operator()( const tbb::blocked_range2d<size_t>& r ) const {

        const std::vector<boost::int32_t>& bcToRunIndex = m_velRIS.get_bc_to_run_index_vector();
        const std::vector<run_data>& runIndexData = m_velRIS.get_run_index_data_vector();
        const vector3& boundsMin = m_bounds.minimum();

        vector3f faceVoxelCoord, faceWorldCoord, finalPosition;

        int ysize = m_bounds.ysize(); //, zsize = m_bounds.zsize();

        for( size_t c = r.rows().begin(); c != r.rows().end(); ++c ) {
            int z = static_cast<int>( c ) + boundsMin.z;
            size_t cOffset = c * ysize;
            for( size_t b = r.cols().begin(); b != r.cols().end(); ++b ) {
                int y = static_cast<int>( b ) + boundsMin.y;
                size_t bcIndex = b + cOffset;

                int runRangeStart = bcToRunIndex[bcIndex];
                int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

                // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the
                // runs with the iterator
                if( runRangeStart != runRangeEnd ||
                    runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                    // construct a block iterator for this scanline

                    // for velocity advection we need one voxel in the negative direction, and two in hte positive
                    // to make sure we have all the data indices that we will need to process the lookup. This means we
                    // will end up with more overall indices, but should be faster than carefully choosing and finding
                    // only the ones we need int xmin=boundsMin.x-centerOffset.x, xmax = boundsMin.x+centerOffset.x;

                    // construct our iterator box
                    boundbox3 box(
                        vector3( boundsMin.x - m_centerOffset.x, y - m_centerOffset.y, z - m_centerOffset.z ),
                        m_boxSize );

                    for( levelset::rle_defined_box_iterator boxIter( m_velRIS, box );
                         !boxIter.is_xplane_finished( m_centerOffset.x ); ++boxIter ) {

                        const boost::int32_t* const dataIndices = boxIter.get_indices();

                        const boundbox3& currentBox = boxIter.current_box();
                        const vector3& currentMin = currentBox.minimum();

                        boost::int32_t cellIndex = dataIndices[cellBoxIndex];
                        // logging::error << "\ncell index: " << cellIndex  << "\n" << endl;
                        // logging::error << "Center " << centerOffset.x << " is_finished: " <<
                        // boxIter.is_xplane_finished(centerOffset.x) << endl;
                        // for( int t = 0; t < boxSize.xsize() ; ++t ){
                        //	logging::error << "Offset " << t << " is_finished: " << boxIter.is_xplane_finished(t) <<
                        //endl;
                        //}

                        // if the center cell is not actually defined then we can skip
                        if( cellIndex < 0 ) {
                            // double check against a naive lookup
                            // int realIndex = ris.XYZtoDataIndex( vector3( currentMin + centerOffset ) );
                            // if (realIndex != cellIndex ) {
                            //	logging::error << "Real Index: " << realIndex << endl;
                            //	logging::error << "Iter Index: " << cellIndex << endl;
                            //	throw std::runtime_error( "Real Index does not line up with the Iterator Index!");
                            //}
                            continue;
                        }

                        // logging::error << "currentmin: " << currentMin << endl;
                        vector3 cellVoxelCoord = currentMin + m_centerOffset;
                        // logging::error << "Cell Coord: " << cellVoxelCoord << endl;
                        vector3f newVelocity;

                        // get the current lookup velocity
                        for( int comp = 0; comp < 3; ++comp ) {
                            // logging::error << "Computing Component: " << comp << endl;

                            // shift to the real world location of the cell sample
                            // we are using only the x,y,z faces so we need to double the index
                            faceVoxelCoord = m_velVcs.get_voxel_face_center( cellVoxelCoord, comp * 2 );

                            faceWorldCoord = m_velVcs.get_world_face_center( cellVoxelCoord, comp * 2 );

                            // trace the final position using the velocity field provided
                            vector3f vel;
                            trilerp_vector3f_staggered_from_indices( m_velocityAcc, faceVoxelCoord, dataIndices,
                                                                     currentBox, vel );

                            vector3f midpoint = faceWorldCoord - 0.5f * vel * m_dt;

                            vector3f vel2;
                            trilerp_vector3f_staggered_from_indices(
                                m_velocityAcc, m_velVcs.get_voxel_coord( midpoint ), dataIndices, currentBox, vel2 );

                            vector3f quarterpoint = faceWorldCoord - 0.75f * vel2 * m_dt;

                            vector3f vel3;
                            trilerp_vector3f_staggered_from_indices( m_velocityAcc,
                                                                     m_velVcs.get_voxel_coord( quarterpoint ),
                                                                     dataIndices, currentBox, vel3 );

                            finalPosition = faceWorldCoord - ( 2 * vel + 3 * vel2 + 4 * vel3 ) * m_dt / 9.f;

                            // look up the value to copy from the source field
                            if( m_sameRIS ) {
                                // since the index specs are the same we can use the same indices we have in the block
                                // iterator to perform the final lookup
                                vector3f finalPosVelocity;
                                trilerp_vector3f_staggered_from_indices( m_srcVelAcc,
                                                                         m_srcVcs.get_voxel_coord( finalPosition ),
                                                                         dataIndices, currentBox, finalPosVelocity );
                                newVelocity[comp] = finalPosVelocity[comp];
                            } else {
                                // we have to perform a full lookup in the source velocity field
                                throw std::runtime_error( "rle_index_specs are not equal! wtf?" );
                            }
                        }
                        // logging::error << "At Index: " << cellIndex << " new vel: " << 	newVelocity << endl;
                        m_resultVelAcc[cellIndex] = newVelocity;
                    }
                }
            }
        }
    }
};

void velocity_advect( const rle_voxel_field& source, rle_voxel_field& result, const rle_voxel_field& velocityField,
                      float dt ) {
    const float maxVoxelMotion = detail::get_maximum_voxel_motion_from_staggered_velocity( velocityField, dt );
    velocity_advect( source, result, velocityField, maxVoxelMotion, dt );
}

void velocity_advect( const rle_voxel_field& source, rle_voxel_field& result, const rle_voxel_field& velocityField,
                      float maxVoxelMotion, float dt ) {
    // TODO : The velocity advection routines can use undefined source
    // regions, which will result in erroneous 0 velocity.  We should
    // use the nearest defined voxels or extrapolated values instead.
    //
    // This same change should be made in velocity lookup part of
    // level set advection.
    tbb::task_scheduler_init taskScheduleInit;

    const frantic::volumetrics::voxel_coord_system& srcVcs = source.get_voxel_coord_system();
    const rle_index_spec& srcRleIndex = source.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> srcVelAcc = source.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    rle_voxel_field temp( srcVcs, srcRleIndex );
    temp.add_channel<vector3f>( _T("StaggeredVelocity") );
    // temp.zero_channel( "StaggeredVelocity" );
    rle_channel_accessor<vector3f> velNext = temp.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    // logging::set_logging_level( 3 );

    FF_LOG( error ) << "**** (vel advect )staggered velocity channel size=" << srcVelAcc.size() << " ****" << endl;
    FF_LOG( debug ) << "Starting iteration" << std::endl;

    const frantic::volumetrics::voxel_coord_system& velocityVcs = velocityField.get_voxel_coord_system();
    const rle_index_spec& velocityRleIndex = velocityField.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> velocityAcc =
        velocityField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    bool sameRleIndexSpecs = velocityRleIndex == srcRleIndex;

    const rle_index_spec& ris = temp.get_rle_index_spec();
    const boundbox3 bounds = ris.outer_bounds();

    // run the tbb:: advector object in parallel (in theory)

    tbb::blocked_range2d<size_t> range( 0, bounds.zsize(), 0, bounds.ysize() );

    logging::error << "threading block range: " << bounds.ysize() << ", " << bounds.zsize() << endl;

    // Use rle lookups for large steps, and sliding windows for smaller
    // steps.
    // TODO : we'd probably be better off using a small window (5^3 or 7^3)
    // all the time, and switching to rle lookups if we go out of bounds.
    if( maxVoxelMotion > 3.5f ) {
        typedef rle_scanline_lookup_advect advector_t;

        advector_t parallelAdvector( bounds, velocityRleIndex, srcRleIndex, srcVcs, velocityVcs, velocityAcc, srcVelAcc,
                                     velNext, dt, sameRleIndexSpecs );
#ifndef FRANTIC_DISABLE_THREADS
        tbb::parallel_for( range, parallelAdvector, tbb::auto_partitioner() );
#else
#pragma message( "Threads are disabled" )
        parallelAdvector( range );
#endif
    } else { // maxVoxelMotion <= 3.5f
        typedef rle_scanline_window_advect advector_t;

        if( maxVoxelMotion > 12.f ) {
            throw std::runtime_error( "velocity_advect Error: the maximum velocity in the advected velocity field (" +
                                      boost::lexical_cast<std::string>( maxVoxelMotion ) +
                                      " voxel lengths) exceeds the maximum limit (12 voxel lengths)." );
        }

        const int indexBoxLength = maxVoxelMotion >= 0 ? advector_t::get_index_box_size( maxVoxelMotion ).xsize() : 5;
        const size3 boxSize( indexBoxLength );
        advector_t parallelAdvector( bounds, boxSize, velocityRleIndex, srcRleIndex, srcVcs, velocityVcs, // resultVcs,
                                     velocityAcc, srcVelAcc, velNext, dt, sameRleIndexSpecs );

#ifndef FRANTIC_DISABLE_THREADS
        tbb::parallel_for( range, parallelAdvector, tbb::auto_partitioner() );
#else
#pragma message( "Threads are disabled" )
        parallelAdvector( range );
#endif
    }

    result.swap( temp );

    // logging::stats << "Advected Field Size=" << result.size() << " timing:" <<
    // velAdvect.psAdvectionMain.last_timing_seconds() << std::endl;
}

void serial_velocity_advect( const rle_voxel_field& source, rle_voxel_field& result,
                             const rle_voxel_field& velocityField, float dt ) {

    const frantic::volumetrics::voxel_coord_system& srcVcs = source.get_voxel_coord_system();
    const rle_index_spec& srcRleIndex = source.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> srcVelAcc = source.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    rle_voxel_field temp( srcVcs, srcRleIndex );
    temp.add_channel<vector3f>( _T("StaggeredVelocity") );
    // temp.zero_channel( "StaggeredVelocity" );
    rle_channel_accessor<vector3f> velNext = temp.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    // logging::set_logging_level( 3 );

    FF_LOG( error ) << "**** (vel advect )staggered velocity channel size=" << srcVelAcc.size() << " ****" << endl;
    FF_LOG( debug ) << "Starting iteration" << std::endl;

    const frantic::volumetrics::voxel_coord_system& velocityVcs = velocityField.get_voxel_coord_system();
    const rle_index_spec& velocityRleIndex = velocityField.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> velocityAcc =
        velocityField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    bool sameRleIndexSpecs = velocityRleIndex == srcRleIndex;

    const rle_index_spec& ris = temp.get_rle_index_spec();
    boundbox3 bounds = ris.outer_bounds();

    vector3 boundsMin = bounds.minimum();
    vector3 boundsMax = bounds.maximum();

    const vector<boost::int32_t>& bcToRunIndex = ris.get_bc_to_run_index_vector();

    int ysize = bounds.ysize();

    // std::vector<int> dataIndices(bounds.get_volume()) ;
    // const boost::int32_t const* dataIndices;

    int count = 0;

    const std::vector<run_data>& runIndexData = ris.get_run_index_data_vector();

    // precompute the index to the center cell (2,2,2) in terms of the local 5x5x5 box
    // this is written out for clarity, but could be collapsed to save a few cyclces...
    size3 boxSize( 5, 5, 5 );
    vector3 centerOffset( 2 );

    vector3f finalPosition;

    int cellBoxIndex =
        centerOffset.x + centerOffset.y * boxSize.xsize() + centerOffset.z * boxSize.xsize() * boxSize.ysize();

    velAdvect.psAdvectionMain.enter();
    for( int z = boundsMin.z; z <= boundsMax.z; ++z ) {
        int cIndex = ris.z_to_c( z );
        for( int y = boundsMin.y; y <= boundsMax.y; ++y ) {
            // compute the bcIndex of the scanline
            int bcIndex = ris.y_to_b( y ) + cIndex * ysize;

            // logging::error << "Computing Scanline: [" << y << ", " << z << "]: bcIndex:" << bcIndex << endl;
            // get the run range of the scanline
            int runRangeStart = bcToRunIndex[bcIndex];
            int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

            // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the runs
            // with the iterator
            if( runRangeStart != runRangeEnd || runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                // construct a block iterator for this scanline

                // for velocity advection we need one voxel in the negative direction, and two in hte positive
                // to make sure we have all the data indices that we will need to process the lookup. This means we will
                // end up with more overall indices, but should be faster than carefully choosing and finding only the
                // ones we need int xmin=boundsMin.x-centerOffset.x, xmax = boundsMin.x+centerOffset.x;

                // construct our iterator box
                // boundbox3 box2(xmin,xmax, y-2,y+2,z-2,z+2);

                boundbox3 box( vector3( boundsMin.x - centerOffset.x, y - centerOffset.y, z - centerOffset.z ),
                               boxSize );
                // if( box2 != box )
                //	throw std::runtime_error("OUr Boxes don't match: " + box.str() + " != " + box2.str() );

                for( levelset::rle_defined_box_iterator boxIter( ris, box );
                     !boxIter.is_xplane_finished( centerOffset.x ); ++boxIter ) {

                    velAdvect.psAdvectionNoIter.enter();

                    // logging::error << "Current Box: " << boxIter.current_box() << endl;
                    //++steps;

                    // velAdvect.psGetIndices.enter();
                    // boxIter.get_indices( &dataIndices[0] );
                    const boost::int32_t* const dataIndices = boxIter.get_indices();

                    // for( int i=0;i<velocityCache.size(); ++i){
                    //	velocityCache[i] = velocityAcc[dataIndices[i]];
                    //}
                    // velAdvect.psGetIndices.exit();

                    const boundbox3& currentBox = boxIter.current_box();
                    const vector3& currentMin = currentBox.minimum();

                    boost::int32_t cellIndex = dataIndices[cellBoxIndex];
                    // logging::error << "\ncell index: " << cellIndex  << "\n" << endl;
                    // logging::error << "Center " << centerOffset.x << " is_finished: " <<
                    // boxIter.is_xplane_finished(centerOffset.x) << endl;
                    // for( int t = 0; t < boxSize.xsize() ; ++t ){
                    //	logging::error << "Offset " << t << " is_finished: " << boxIter.is_xplane_finished(t) << endl;
                    //}

                    // if the center cell is not actually defined then we can skip
                    if( cellIndex < 0 ) {
                        // double check against a naive lookup
                        // int realIndex = ris.XYZtoDataIndex( vector3( currentMin + centerOffset ) );
                        // if (realIndex != cellIndex ) {
                        //	logging::error << "Real Index: " << realIndex << endl;
                        //	logging::error << "Iter Index: " << cellIndex << endl;
                        //	throw std::runtime_error( "Real Index does not line up with the Iterator Index!");
                        //}
                        velAdvect.psAdvectionNoIter.exit();
                        continue;
                    }

                    // logging::error << "currentmin: " << currentMin << endl;
                    // debugVisitedIndices.push_back( cellIndex );

                    vector3 cellVoxelCoord = currentMin + centerOffset;
                    // logging::error << "Cell Coord: " << cellVoxelCoord << endl;
                    vector3f newVelocity;

                    // get the current lookup velocity
                    for( int comp = 0; comp < 3; ++comp ) {
                        // logging::error << "Computing Component: " << comp << endl;

                        // shift to the real world location of the cell sample
                        // we are using only the x,y,z faces so we need to double the index
                        vector3f faceVoxelCoord = srcVcs.get_voxel_face_center( cellVoxelCoord, comp * 2 );

                        vector3f faceWorldCoord = srcVcs.get_world_face_center( cellVoxelCoord, comp * 2 );
                        // if( velocityVcs.get_voxel_coord( faceWorldCoord) != faceVoxelCoord )
                        //	throw std::runtime_error("Face Coords do not line up");

                        // trace the final position using the velocity field provided

                        velAdvect.psTraceRK.enter();

                        velAdvect.psLerp.enter();
                        vector3f vel;
                        trilerp_vector3f_staggered_from_indices( velocityAcc, faceVoxelCoord, dataIndices, currentBox,
                                                                 vel );
                        velAdvect.psLerp.exit();

                        // FF_LOG(debug) << "\tvel 1: " << vel << std::endl;
                        // peform the trace

                        vector3f midpoint = faceWorldCoord - 0.5f * vel * dt;
                        // FF_LOG(debug)<< "\tmid: " << midpoint << std::endl;

                        velAdvect.psLerp.enter();
                        trilerp_vector3f_staggered_from_indices( velocityAcc, velocityVcs.get_voxel_coord( midpoint ),
                                                                 dataIndices, currentBox, vel );
                        velAdvect.psLerp.exit();

                        // logging::error << "LookUp 2" << endl;

                        finalPosition = faceWorldCoord - vel * dt;
                        velAdvect.psTraceRK.exit();

                        // look up the value to copy from the source field
                        vector3f velocity;
                        if( sameRleIndexSpecs ) {
                            // since the index specs are the same we can use the same indices we have in the block
                            // iterator to perform the final lookup
                            velAdvect.psLerp.enter();
                            trilerp_vector3f_staggered_from_indices( srcVelAcc, srcVcs.get_voxel_coord( finalPosition ),
                                                                     dataIndices, currentBox, velocity );
                            velAdvect.psLerp.exit();
                        } else {

                            // we have to perform a full lookup in the source velocity field
                            throw std::runtime_error( "rle_index_specs are not equal! wtf?" );
                        }

                        newVelocity[comp] = velocity[comp];
                        if( fabsf( newVelocity[comp] ) > 1e28f ) {
                            logging::error << "Velocity=" << newVelocity[comp] << endl;
                            logging::error << "Voxel Coord=" << srcVcs.get_voxel_coord( finalPosition ) << endl;
                            throw std::runtime_error( "serial_velocity_advect() Error: NaN Velocity found!" );
                        }

                        ++count;
                    }

                    // logging::error << "At Index: " << cellIndex << " new vel: " << 	newVelocity << endl;
                    velNext[cellIndex] = newVelocity;
                    velAdvect.psAdvectionNoIter.exit();
                }

                // logging::error << " There were " << steps << " steps taken in this scaline with " << negIndices << "
                // negative indices founds" << endl;
            }
            // else {
            //	int index = ris.XYZtoDataIndex( vector3(boundsMin.x, y, z)  );
            //	if ( index != ris.get_run_index_data_vector()[runRangeStart].dataIndex ) {
            //		logging::error << "RunRange: ( " << runRangeStart << ", " << runRangeEnd << " )" << endl;
            //		logging::error << "Naive DataIndex: " << index << endl;
            //		logging::error << "Run DataIndex: " << ris.get_run_index_data_vector()[runRangeStart].dataIndex  <<
            //endl; 		throw std::runtime_error( "BC Index lookup has failed the generic XYZToDataIndex lookup" );
            //	}

            //	// logging::error << "Empty Scanline, Coord: [" << y << "," << z << "] BC: " << bcIndex << " Run Range:
            //["
            //  //                << runRangeStart << ", " << runRangeEnd << "]" << endl;
            //}
        }
    }
    velAdvect.psAdvectionMain.exit();

    result.swap( temp );

    // if( !result.check_consistency(cout) )
    //	throw std::runtime_error("velocity_advect() - (post) Result Voxel Field Failed Consistency Check");

    logging::stats << "Advected Field Size=" << result.size()
                   << " timing:" << frantic::strings::to_tstring( velAdvect.psAdvectionMain.last_timing_seconds() )
                   << std::endl;
}

void serial_velocity_weno_advect( rle_voxel_field& source, rle_voxel_field& result, rle_voxel_field& velocityField,
                                  float dt ) {

    create_staggered_smoothness_indicator_x_channel( velocityField, _T("StaggeredVelocity"), _T("WenoIndicatorOne"),
                                                     _T("WenoIndicatorTwo") );
    create_staggered_smoothness_indicator_x_channel( source, _T("StaggeredVelocity"), _T("WenoIndicatorOne"),
                                                     _T("WenoIndicatorTwo") );

    const frantic::volumetrics::voxel_coord_system& srcVcs = source.get_voxel_coord_system();
    const rle_index_spec& srcRleIndex = source.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> srcVelAcc = source.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    rle_voxel_field temp( srcVcs, srcRleIndex );
    temp.add_channel<vector3f>( _T("StaggeredVelocity") );
    // temp.zero_channel( "StaggeredVelocity" );
    rle_channel_accessor<vector3f> velNext = temp.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    // logging::set_logging_level( 3 );

    FF_LOG( error ) << "**** (vel advect )staggered velocity channel size=" << srcVelAcc.size() << " ****" << endl;
    FF_LOG( debug ) << "Starting iteration" << std::endl;

    const frantic::volumetrics::voxel_coord_system& velocityVcs = velocityField.get_voxel_coord_system();
    const rle_index_spec& velocityRleIndex = velocityField.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> velocityAcc =
        velocityField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    bool sameRleIndexSpecs = velocityRleIndex == srcRleIndex;

    const rle_index_spec& ris = temp.get_rle_index_spec();
    boundbox3 bounds = ris.outer_bounds();

    const vector3& boundsMin = bounds.minimum();
    const vector3& boundsMax = bounds.maximum();

    const vector<boost::int32_t>& bcToRunIndex = ris.get_bc_to_run_index_vector();

    int ysize = bounds.ysize();
    int count = 0;

    const std::vector<run_data>& runIndexData = ris.get_run_index_data_vector();

    const_rle_channel_accessor<float> indicator1Acc =
        velocityField.get_channel_accessor<float>( _T("WenoIndicatorOne" ) );
    const_rle_channel_accessor<float> indicator2Acc =
        velocityField.get_channel_accessor<float>( _T("WenoIndicatorTwo") );

    const_rle_channel_accessor<float> srcIndicator1Acc = source.get_channel_accessor<float>( _T("WenoIndicatorOne") );
    const_rle_channel_accessor<float> srcIndicator2Acc = source.get_channel_accessor<float>( _T("WenoIndicatorTwo") );

    // precompute the index to the center cell (2,2,2) in terms of the local 5x5x5 box
    // this is written out for clarity, but could be collapsed to save a few cyclces...
    size3 boxSize( 7, 7, 7 );
    vector3 centerOffset( 3 );

    vector3f finalPosition;

    int cellBoxIndex =
        centerOffset.x + centerOffset.y * boxSize.xsize() + centerOffset.z * boxSize.xsize() * boxSize.ysize();

    for( int z = boundsMin.z; z <= boundsMax.z; ++z ) {
        int cIndex = ris.z_to_c( z );
        for( int y = boundsMin.y; y <= boundsMax.y; ++y ) {
            // compute the bcIndex of the scanline
            int bcIndex = ris.y_to_b( y ) + cIndex * ysize;

            // logging::error << "Computing Scanline: [" << y << ", " << z << "]: bcIndex:" << bcIndex << endl;
            // get the run range of the scanline
            int runRangeStart = bcToRunIndex[bcIndex];
            int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

            // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the runs
            // with the iterator
            if( runRangeStart != runRangeEnd || runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                // construct a block iterator for this scanline

                // for velocity advection we need one voxel in the negative direction, and two in hte positive
                // to make sure we have all the data indices that we will need to process the lookup. This means we will
                // end up with more overall indices, but should be faster than carefully choosing and finding only the
                // ones we need
                boundbox3 box( vector3( boundsMin.x - centerOffset.x, y - centerOffset.y, z - centerOffset.z ),
                               boxSize );

                for( levelset::rle_defined_box_iterator boxIter( ris, box );
                     !boxIter.is_xplane_finished( centerOffset.x ); ++boxIter ) {

                    velAdvect.psAdvectionNoIter.enter();

                    const boost::int32_t* const dataIndices = boxIter.get_indices();

                    const boundbox3& currentBox = boxIter.current_box();
                    const vector3& currentMin = currentBox.minimum();

                    boost::int32_t cellIndex = dataIndices[cellBoxIndex];
                    // if the center cell is not actually defined then we can skip
                    if( cellIndex < 0 ) {
                        // logging::error << "undef: "  << currentMin + centerOffset << std::endl;
                        continue;
                    }

                    vector3 cellVoxelCoord = currentMin + centerOffset;

                    vector3f newVelocity;

                    // get the current lookup velocity
                    for( int comp = 0; comp < 3; ++comp ) {

                        // logging::error << "Computing Component: " << comp << endl;

                        vector3f faceWorldCoord, midpoint, finalPosition, velocity;
                        vector3f vel1, vel2;

                        try {
                            // shift to the real world location of the cell sample
                            // we are using only the x,y,z faces so we need to double the index
                            vector3f faceVoxelCoord = srcVcs.get_voxel_face_center( cellVoxelCoord, comp * 2 );

                            faceWorldCoord = srcVcs.get_world_face_center( cellVoxelCoord, comp * 2 );

                            // trace the final position using the velocity field provided
                            vector3f vel = staggered_weno3_lookup( dataIndices, boxSize, currentMin, indicator1Acc,
                                                                   indicator2Acc, velocityAcc, faceVoxelCoord );
                            midpoint = faceWorldCoord - 0.5f * vel * dt;
                            // FF_LOG(debug)<< "\tmid: " << midpoint << std::endl;

                            vel =
                                staggered_weno3_lookup( dataIndices, boxSize, currentMin, indicator1Acc, indicator2Acc,
                                                        velocityAcc, velocityVcs.get_voxel_coord( midpoint ) );
                            //
                            // logging::error << "LookUp 2" << endl;
                            finalPosition = faceWorldCoord - vel * dt;

                            // look up the value to copy from the source field

                            if( sameRleIndexSpecs ) {
                                /*if( comp==2 && cellVoxelCoord.x == boundsMin.x && cellVoxelCoord.y == (boundsMin.y )
                                && cellVoxelCoord.z == boundsMin.z+6 ) { logging::set_logging_level(5);
                                }*/
                                // since the index specs are the same we can use the same indices we have in the block
                                // iterator to perform the final lookup
                                velocity = staggered_weno3_lookup( dataIndices, boxSize, currentMin, srcIndicator1Acc,
                                                                   srcIndicator2Acc, srcVelAcc,
                                                                   srcVcs.get_voxel_coord( finalPosition ) );
                                /*if( comp ==2 && cellVoxelCoord.x == boundsMin.x&& cellVoxelCoord.y == (boundsMin.y )
                                && cellVoxelCoord.z == boundsMin.z+6 ) { logging::set_logging_level(3);
                                }*/

                            } else {

                                // we have to perform a full lookup in the source velocity field
                                throw std::runtime_error( "rle_index_specs are not equal! wtf?" );
                            }
                        } catch( std::exception& e ) {

                            staggered_weno3_debug_dump( dataIndices, boxSize, currentMin, indicator1Acc, indicator2Acc,
                                                        velocityAcc, velocityVcs.get_voxel_coord( midpoint ) );

                            logging::error << "v1: " << vel1 << std::endl;
                            logging::error << "v2: " << vel2 << std::endl;

                            vector3f lerp1, lerp2;
                            trilerp_vector3f_staggered( velocityRleIndex, velocityAcc,
                                                        velocityVcs.get_voxel_face_center( cellVoxelCoord, comp * 2 ),
                                                        lerp1 );
                            trilerp_vector3f_staggered( velocityRleIndex, velocityAcc,
                                                        velocityVcs.get_voxel_coord( midpoint ), lerp2 );
                            logging::error << "lerp1: " << lerp1 << std::endl;
                            logging::error << "lerp2: " << lerp2 << std::endl;

                            logging::error << "face:  " << velocityVcs.get_voxel_coord( faceWorldCoord ) << std::endl;
                            logging::error << "mid:   " << velocityVcs.get_voxel_coord( midpoint ) << std::endl;
                            logging::error << "final: " << velocityVcs.get_voxel_coord( finalPosition ) << std::endl;
                            throw std::runtime_error( e.what() );
                        }

                        newVelocity[comp] = velocity[comp];
                        if( fabsf( newVelocity[comp] ) > 1e28f ) {
                            logging::error << "Velocity=" << newVelocity[comp] << endl;
                            logging::error << "Voxel Coord=" << srcVcs.get_voxel_coord( finalPosition ) << endl;
                            throw std::runtime_error( "serial_velocity_weno_advect Error: NaN Velocity found!" );
                        }

                        ++count;
                    }

                    // logging::error << "At Index: " << cellIndex << " new vel: " << 	newVelocity << endl;
                    velNext[cellIndex] = newVelocity;
                }
            }
        }
    }

    result.swap( temp );

    // clean up the indicator channels
    velocityField.erase_channel( _T("WenoIndicatorOne") );
    velocityField.erase_channel( _T("WenoIndicatorTwo") );

    source.erase_channel( _T("WenoIndicatorOne") );
    source.erase_channel( _T("WenoIndicatorTwo") );

    logging::stats << "Advected Field Size=" << result.size()
                   << " timing:" << frantic::strings::to_tstring( velAdvect.psAdvectionMain.last_timing_seconds() )
                   << std::endl;
}

class rle_scanline_weno3_advect {
    const boundbox3& m_bounds;
    const size3& m_boxSize;
    const rle_index_spec &m_velRIS, m_srcRIS;
    const voxel_coord_system &m_srcVcs, m_velVcs;
    vector3 m_centerOffset;
    int cellBoxIndex;
    bool m_sameRIS;

    // precompute the index to the center cell (2,2,2) in terms of the local 5x5x5 box
    // this is written out for clarity, but could be collapsed to save a few cyclces...

    const_rle_channel_accessor<vector3f>& m_velocityAcc;
    const_rle_channel_accessor<vector3f>& m_srcVelAcc;
    rle_channel_accessor<vector3f>& m_resultVelAcc;

    const_rle_channel_accessor<float>&m_velIndicatorX1Acc, &m_velIndicatorX2Acc;
    const_rle_channel_accessor<float>&m_srcIndicatorX1Acc, &m_srcIndicatorX2Acc;

    float m_dt;

    rle_scanline_weno3_advect& operator=( const rle_scanline_weno3_advect& ); // not implemented

  public:
    rle_scanline_weno3_advect(
        const boundbox3& bounds, const size3& boxSize, const rle_index_spec& velRIS, const rle_index_spec& srcRIS,
        const voxel_coord_system& srcVcs, const voxel_coord_system& velVcs,
        const_rle_channel_accessor<vector3f>& velocityAcc, const_rle_channel_accessor<vector3f>& srcVelAcc,
        rle_channel_accessor<vector3f>& resultVelAcc, const_rle_channel_accessor<float>& velIndicatorX1Acc,
        const_rle_channel_accessor<float>& velIndicatorX2Acc, const_rle_channel_accessor<float>& srcIndicatorX1Acc,
        const_rle_channel_accessor<float>& srcIndicatorX2Acc, float dt, bool isSameRleIndexSpecs )
        : m_bounds( bounds )
        , m_boxSize( boxSize )
        , m_velRIS( velRIS )
        , m_srcRIS( srcRIS )
        , m_srcVcs( srcVcs )
        , m_velVcs( velVcs )
        , m_sameRIS( isSameRleIndexSpecs )
        , m_velocityAcc( velocityAcc )
        , m_srcVelAcc( srcVelAcc )
        , m_resultVelAcc( resultVelAcc )
        , m_velIndicatorX1Acc( velIndicatorX1Acc )
        , m_velIndicatorX2Acc( velIndicatorX2Acc )
        , m_srcIndicatorX1Acc( srcIndicatorX1Acc )
        , m_srcIndicatorX2Acc( srcIndicatorX2Acc )
        , m_dt( dt ) {
        m_centerOffset = vector3( m_boxSize.xsize() / 2 );
        cellBoxIndex = m_centerOffset.x + m_centerOffset.y * m_boxSize.xsize() +
                       m_centerOffset.z * m_boxSize.xsize() * m_boxSize.ysize();
    }

    static size3 get_index_box_size( float maxVoxelMotion ) {
        const int boxLength = 5 + 2 * static_cast<int>( max<float>( 0, ceilf( maxVoxelMotion - 0.49f ) ) );
        return size3( boxLength );
    }

    void operator()( const tbb::blocked_range2d<size_t>& r ) const {

        const std::vector<boost::int32_t>& bcToRunIndex = m_velRIS.get_bc_to_run_index_vector();
        const std::vector<run_data>& runIndexData = m_velRIS.get_run_index_data_vector();
        const vector3& boundsMin = m_bounds.minimum();

        vector3f faceVoxelCoord, faceWorldCoord, vel, finalPosition;

        int ysize = m_bounds.ysize(); //, zsize = m_bounds.zsize();

        for( size_t c = r.rows().begin(); c != r.rows().end(); ++c ) {
            int z = static_cast<int>( c ) + boundsMin.z;
            size_t cOffset = c * ysize;
            for( size_t b = r.cols().begin(); b != r.cols().end(); ++b ) {
                int y = static_cast<int>( b ) + boundsMin.y;
                size_t bcIndex = b + cOffset;

                // logging::error << "Computing Scanline: [" << y << ", " << z << "]: bcIndex:" << bcIndex << endl;
                // get the run range of the scanline
                int runRangeStart = bcToRunIndex[bcIndex];
                int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

                // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the
                // runs with the iterator
                if( runRangeStart != runRangeEnd ||
                    runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                    // construct a block iterator for this scanline

                    // for velocity advection we need one voxel in the negative direction, and two in hte positive
                    // to make sure we have all the data indices that we will need to process the lookup. This means we
                    // will end up with more overall indices, but should be faster than carefully choosing and finding
                    // only the ones we need
                    boundbox3 box(
                        vector3( boundsMin.x - m_centerOffset.x, y - m_centerOffset.y, z - m_centerOffset.z ),
                        m_boxSize );

                    for( levelset::rle_defined_box_iterator boxIter( m_velRIS, box );
                         !boxIter.is_xplane_finished( m_centerOffset.x ); ++boxIter ) {

                        const boost::int32_t* const dataIndices = boxIter.get_indices();

                        const boundbox3& currentBox = boxIter.current_box();
                        const vector3& currentMin = currentBox.minimum();

                        boost::int32_t cellIndex = dataIndices[cellBoxIndex];
                        // if the center cell is not actually defined then we can skip
                        if( cellIndex < 0 ) {
                            // logging::error << "undef: "  << currentMin + centerOffset << std::endl;
                            continue;
                        }

                        vector3 cellVoxelCoord = currentMin + m_centerOffset;

                        vector3f newVelocity;

                        // get the current lookup velocity
                        for( int comp = 0; comp < 3; ++comp ) {
                            vector3f faceWorldCoord, midpoint, quarterpoint, finalPosition, velocity;
                            try {
                                // shift to the real world location of the cell sample
                                // we are using only the x,y,z faces so we need to double the index
                                vector3f faceVoxelCoord = m_srcVcs.get_voxel_face_center( cellVoxelCoord, comp * 2 );

                                faceWorldCoord = m_srcVcs.get_world_face_center( cellVoxelCoord, comp * 2 );

                                vector3f vel1 =
                                    staggered_weno3_lookup( dataIndices, m_boxSize, currentMin, m_velIndicatorX1Acc,
                                                            m_velIndicatorX2Acc, m_velocityAcc, faceVoxelCoord );
                                midpoint = faceWorldCoord - 0.5f * vel1 * m_dt;

                                vector3f vel2 = staggered_weno3_lookup(
                                    dataIndices, m_boxSize, currentMin, m_velIndicatorX1Acc, m_velIndicatorX2Acc,
                                    m_velocityAcc, m_velVcs.get_voxel_coord( midpoint ) );

                                quarterpoint = faceWorldCoord - 0.75f * vel2 * m_dt;

                                vector3f vel3 = staggered_weno3_lookup(
                                    dataIndices, m_boxSize, currentMin, m_velIndicatorX1Acc, m_velIndicatorX2Acc,
                                    m_velocityAcc, m_velVcs.get_voxel_coord( quarterpoint ) );

                                finalPosition = faceWorldCoord - ( 2 * vel1 + 3 * vel2 + 4 * vel3 ) * m_dt / 9.f;

                                // look up the value to copy from the source field
                                if( m_sameRIS ) {
                                    velocity = staggered_weno3_lookup(
                                        dataIndices, m_boxSize, currentMin, m_srcIndicatorX1Acc, m_srcIndicatorX2Acc,
                                        m_srcVelAcc, m_srcVcs.get_voxel_coord( finalPosition ) );
                                } else {

                                    // we have to perform a full lookup in the source velocity field
                                    throw std::runtime_error( "rle_index_specs are not equal! wtf?" );
                                }
                            } catch( std::exception& e ) {

                                staggered_weno3_debug_dump( dataIndices, m_boxSize, currentMin, m_velIndicatorX2Acc,
                                                            m_velIndicatorX2Acc, m_velocityAcc,
                                                            m_velVcs.get_voxel_coord( midpoint ) );

                                vector3f lerp1, lerp2;
                                trilerp_vector3f_staggered( m_velRIS, m_velocityAcc,
                                                            m_velVcs.get_voxel_face_center( cellVoxelCoord, comp * 2 ),
                                                            lerp1 );
                                trilerp_vector3f_staggered( m_velRIS, m_velocityAcc,
                                                            m_velVcs.get_voxel_coord( midpoint ), lerp2 );
                                logging::error << "lerp1: " << lerp1 << std::endl;
                                logging::error << "lerp2: " << lerp2 << std::endl;

                                logging::error << "face:  " << m_velVcs.get_voxel_coord( faceWorldCoord ) << std::endl;
                                logging::error << "mid:   " << m_velVcs.get_voxel_coord( midpoint ) << std::endl;
                                logging::error << "final: " << m_velVcs.get_voxel_coord( finalPosition ) << std::endl;
                                throw std::runtime_error( e.what() );
                            }

                            newVelocity[comp] = velocity[comp];

                            if( fabsf( newVelocity[comp] ) > 1e28f ) {
                                logging::error << "Velocity=" << newVelocity[comp] << endl;
                                logging::error << "Voxel Coord=" << m_srcVcs.get_voxel_coord( finalPosition ) << endl;
                                throw std::runtime_error( "rle_scanline_weno3_advect() Error: NaN Velocity found!" );
                            }
                        }

                        // logging::error << "At Index: " << cellIndex << " new vel: " << 	newVelocity << endl;
                        m_resultVelAcc[cellIndex] = newVelocity;
                    }
                }
            }
        }
    }
};

void velocity_weno_advect( rle_voxel_field& source, rle_voxel_field& result, rle_voxel_field& velocityField,
                           float maxVoxelMotion, float dt ) {
    tbb::task_scheduler_init taskSchedulerInit;
    typedef rle_scanline_weno3_advect advector_t;

    create_staggered_smoothness_indicator_x_channel( velocityField, _T("StaggeredVelocity"), _T("WenoIndicatorOne"),
                                                     _T("WenoIndicatorTwo") );
    create_staggered_smoothness_indicator_x_channel( source, _T("StaggeredVelocity"), _T("WenoIndicatorOne"),
                                                     _T("WenoIndicatorTwo") );

    const frantic::volumetrics::voxel_coord_system& srcVcs = source.get_voxel_coord_system();
    const rle_index_spec& srcRleIndex = source.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> srcVelAcc = source.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    rle_voxel_field temp( srcVcs, srcRleIndex );
    temp.add_channel<vector3f>( _T("StaggeredVelocity") );
    // temp.zero_channel( "StaggeredVelocity" );
    rle_channel_accessor<vector3f> velNext = temp.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    // logging::set_logging_level( 3 );

    FF_LOG( error ) << "**** (vel advect )staggered velocity channel size=" << srcVelAcc.size() << " ****" << endl;
    FF_LOG( debug ) << "Starting iteration" << std::endl;

    const frantic::volumetrics::voxel_coord_system& velocityVcs = velocityField.get_voxel_coord_system();
    const rle_index_spec& velocityRleIndex = velocityField.get_rle_index_spec();
    const_rle_channel_accessor<vector3f> velocityAcc =
        velocityField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    bool sameRleIndexSpecs = velocityRleIndex == srcRleIndex;

    const rle_index_spec& ris = temp.get_rle_index_spec();
    boundbox3 bounds = ris.outer_bounds();

    // const vector<boost::int32_t>& bcToRunIndex = ris.get_bc_to_run_index_vector();

    // int ysize = bounds.ysize();
    // int count =0;

    // const std::vector<run_data>& runIndexData = ris.get_run_index_data_vector();

    const_rle_channel_accessor<float> indicator1Acc =
        velocityField.get_channel_accessor<float>( _T("WenoIndicatorOne" ) );
    const_rle_channel_accessor<float> indicator2Acc =
        velocityField.get_channel_accessor<float>( _T("WenoIndicatorTwo") );

    const_rle_channel_accessor<float> srcIndicator1Acc = source.get_channel_accessor<float>( _T("WenoIndicatorOne") );
    const_rle_channel_accessor<float> srcIndicator2Acc = source.get_channel_accessor<float>( _T("WenoIndicatorTwo") );

    // precompute the index to the center cell (2,2,2) in terms of the local 5x5x5 box
    // this is written out for clarity, but could be collapsed to save a few cyclces...
    const size3 boxSize = advector_t::get_index_box_size( maxVoxelMotion );

    advector_t parallelAdvector( bounds, boxSize, velocityRleIndex, srcRleIndex, srcVcs, velocityVcs, velocityAcc,
                                 srcVelAcc, velNext, indicator1Acc, indicator2Acc, srcIndicator1Acc, srcIndicator2Acc,
                                 dt, sameRleIndexSpecs );

    // run the tbb:: advector object in parallel (in theory)
    tbb::blocked_range2d<size_t> range( 0, bounds.zsize(), 0, bounds.ysize() );

    logging::error << "threading block range: " << bounds.ysize() << ", " << bounds.zsize() << endl;
#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, parallelAdvector, tbb::auto_partitioner() );
#else
#pragma message( "Threads are disabled" )
    parallelAdvector( range );
#endif

    result.swap( temp );

    // clean up the indicator channels
    velocityField.erase_channel( _T("WenoIndicatorOne") );
    velocityField.erase_channel( _T("WenoIndicatorTwo") );

    source.erase_channel( _T("WenoIndicatorOne") );
    source.erase_channel( _T("WenoIndicatorTwo") );

    logging::stats << "Advected Field Size=" << result.size()
                   << " timing:" << frantic::strings::to_tstring( velAdvect.psAdvectionMain.last_timing_seconds() )
                   << std::endl;
}

void velocity_weno_advect( rle_voxel_field& source, rle_voxel_field& result, rle_voxel_field& velocityField,
                           float dt ) {
    const float maxVoxelMotion = detail::get_maximum_voxel_motion_from_staggered_velocity( velocityField, dt );
    velocity_weno_advect( source, result, velocityField, maxVoxelMotion, dt );
}

// const rle_voxel_field& source,
void trace_runge_kutta_2( const voxel_coord_system& vcs, const rle_index_spec& ris,
                          const_rle_channel_accessor<vector3f>& velAcc, const vector3f& worldLocation, float dt,
                          vector3f& finalPosition ) {
    cached_trilerp_vector3f_staggered velocityTrilerp( ris );

    vector3f vel;
    velocityTrilerp.get( velAcc, vcs.get_voxel_coord( worldLocation ), vel );

    const vector3f midpoint = worldLocation + 0.5f * vel * dt;
    velocityTrilerp.get( velAcc, vcs.get_voxel_coord( midpoint ), vel );

    finalPosition = worldLocation + vel * dt;
}

void trace_runge_kutta_2_unstaggered( const voxel_coord_system& vcs, const rle_index_spec& ris,
                                      const_rle_channel_accessor<vector3f>& velAcc, const vector3f& worldLocation,
                                      float dt, vector3f& finalPosition ) {
    vector3f vel;
    trilerp_vector3f( ris, velAcc, vcs.get_voxel_coord( worldLocation ), vel );

    if( logging::is_logging_debug() ) {
        std::cout << "lookup voxel coord: " << vcs.get_voxel_coord( worldLocation ) << std::endl;
        std::cout << "vel: " << vel << std::endl;
    }

    // peform the trace
    vector3f midpoint = worldLocation + 0.5f * vel * dt;
    trilerp_vector3f( ris, velAcc, vcs.get_voxel_coord( midpoint ), vel );

    if( logging::is_logging_debug() ) {
        std::cout << "lookup mid voxel coord: " << vcs.get_voxel_coord( midpoint ) << std::endl;
        std::cout << "vel: " << vel << std::endl;
    }

    finalPosition = worldLocation + vel * dt;

    if( logging::is_logging_debug() ) {
        std::cout << "final voxel coord: " << vcs.get_voxel_coord( finalPosition ) << std::endl;
    }
}

void trace_runge_kutta_3( const voxel_coord_system& vcs, const rle_index_spec& ris,
                          const_rle_channel_accessor<vector3f>& velAcc, const vector3f& worldLocation, float dt,
                          vector3f& finalPosition ) {
    cached_trilerp_vector3f_staggered velocityTrilerp( ris );
    vector3f vel, vel2, vel3;

    velocityTrilerp.get( velAcc, vcs.get_voxel_coord( worldLocation ), vel );

    const vector3f midpoint = worldLocation + 0.5f * dt * vel;
    velocityTrilerp.get( velAcc, vcs.get_voxel_coord( midpoint ), vel2 );

    const vector3f quarterpoint = worldLocation + 0.75f * dt * vel2;
    velocityTrilerp.get( velAcc, vcs.get_voxel_coord( quarterpoint ), vel3 );

    finalPosition = worldLocation + dt * ( ( 2.f * vel ) + ( 3.f * vel2 ) + ( 4.f * vel3 ) ) / 9.f;
}

void trace_runge_kutta_4( const voxel_coord_system& vcs, const rle_index_spec& ris,
                          const_rle_channel_accessor<vector3f>& velAcc, const vector3f& worldLocation, float dt,
                          vector3f& finalPosition ) {
    // TODO: this should use the "from_indices" version and reuse indices if it can
    vector3f vel1, vel2, vel3, vel4;
    trilerp_vector3f_staggered( ris, velAcc, vcs.get_voxel_coord( worldLocation ), vel1 );
    trilerp_vector3f_staggered( ris, velAcc, vcs.get_voxel_coord( worldLocation + 0.5f * dt * vel1 ), vel2 );
    trilerp_vector3f_staggered( ris, velAcc, vcs.get_voxel_coord( worldLocation + 0.5f * dt * vel2 ), vel3 );
    trilerp_vector3f_staggered( ris, velAcc, vcs.get_voxel_coord( worldLocation + dt * vel3 ), vel4 );

    finalPosition = worldLocation + dt * ( vel1 + 2 * ( vel2 + vel3 ) + vel4 ) / 6.0f;
}

void compute_curl_channel( rle_voxel_field& field, const frantic::tstring staggeredVectorChannel,
                           const frantic::tstring& curlChannel ) {

    field.add_channel<vector3f>( curlChannel );

    const_rle_channel_accessor<vector3f> vector = field.get_channel_accessor<vector3f>( staggeredVectorChannel );
    rle_channel_accessor<vector3f> curl = field.get_channel_accessor<vector3f>( curlChannel );

    rle_defined_and_adj_iterator i( field.get_rle_index_spec() ), ie( field.get_rle_index_spec(), true );

    float dx = field.get_voxel_coord_system().voxel_length();

    for( ; i != ie; ++i ) {

        boost::int32_t xpos = i.get_x_pos_data_index(), ypos = i.get_y_pos_data_index(),
                       zpos = i.get_z_pos_data_index();

        boost::int32_t dataIndex = i.get_center_data_index();

        const vector3f& vectorNeg = vector[dataIndex];
        vector3f& curlValue = curl[dataIndex];

        vector3f vectorXpos = ( xpos >= 0 ) ? vector[xpos] : vector3f();
        vector3f vectorYpos = ( ypos >= 0 ) ? vector[ypos] : vector3f();
        vector3f vectorZpos = ( zpos >= 0 ) ? vector[zpos] : vector3f();

        curlValue.x = ( ( vectorYpos.z - vectorNeg.z ) - ( vectorZpos.y - vectorNeg.y ) ) / dx;
        curlValue.y = ( ( vectorZpos.x - vectorNeg.x ) - ( vectorXpos.z - vectorNeg.z ) ) / dx;
        curlValue.z = ( ( vectorXpos.y - vectorNeg.y ) - ( vectorYpos.x - vectorNeg.x ) ) / dx;
    }
}

void compute_convective_derivative_channel( rle_voxel_field& field, const frantic::tstring staggeredVectorChannel,
                                            const frantic::tstring& derivativeChannel ) {

    const rle_index_spec& ris = field.get_rle_index_spec();

    field.add_channel<vector3f>( derivativeChannel );

    const_rle_channel_accessor<vector3f> vector = field.get_channel_accessor<vector3f>( staggeredVectorChannel );
    rle_channel_accessor<vector3f> derivative = field.get_channel_accessor<vector3f>( derivativeChannel );

    float dx = field.get_voxel_coord_system().voxel_length();

    // rle_defined_and_adj_iterator i( field.get_rle_index_spec() ), ie( field.get_rle_index_spec(),true);

    size3 boxSize( 3 );

    boundbox3 bounds = ris.outer_bounds();
    const int ysize = bounds.ysize();

    const vector3& boundsMin = bounds.minimum();
    const vector3& boundsMax = bounds.maximum();

    const std::vector<boost::int32_t>& bcToRunIndex = ris.get_bc_to_run_index_vector();
    const std::vector<run_data>& runIndexData = ris.get_run_index_data_vector();

    int zIndexOffset = boxSize.xsize() * boxSize.ysize();
    int yIndexOffset = boxSize.xsize();

    boost::int32_t centerBoxIndex = 1 + yIndexOffset + zIndexOffset;

    tbb::blocked_range2d<int> r( boundsMin.z, boundsMax.z + 1, boundsMin.y, boundsMax.y + 1 );

    for( int z = r.rows().begin(); z != r.rows().end(); ++z ) {
        int cIndex = ris.z_to_c( z );
        for( int y = r.cols().begin(); y != r.cols().end(); ++y ) {
            // compute the bcIndex of the scanline
            int bcIndex = ris.y_to_b( y ) + cIndex * ysize;

            // logging::error << "Computing Scanline: [" << y << ", " << z << "]: bcIndex:" << bcIndex << endl;
            // get the run range of the scanline
            int runRangeStart = bcToRunIndex[bcIndex];
            int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

            // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the runs
            // with the iterator
            if( runRangeStart != runRangeEnd || runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {

                boundbox3 box( vector3( boundsMin.x - 1, y - 1, z - 1 ), boxSize );

                // std::cout << "Starting Box: " << box << endl;

                for( rle_defined_box_iterator boxIter( ris, box ); !boxIter.is_xplane_finished( 1 ); ++boxIter ) {
                    const boost::int32_t* dataIndices = boxIter.get_indices();

                    if( dataIndices[centerBoxIndex] < 0 ) {
                        continue;
                    }

                    vector3f zeroVector;

                    const vector3f& centerVoxel = vector[dataIndices[centerBoxIndex]];

                    const vector3f& centerXposVoxel =
                        ( dataIndices[centerBoxIndex + 1] >= 0 ) ? vector[dataIndices[centerBoxIndex + 1]] : zeroVector;
                    const vector3f& centerYposVoxel = ( dataIndices[centerBoxIndex + yIndexOffset] >= 0 )
                                                          ? vector[dataIndices[centerBoxIndex + yIndexOffset]]
                                                          : zeroVector;
                    const vector3f& centerZposVoxel = ( dataIndices[centerBoxIndex + zIndexOffset] >= 0 )
                                                          ? vector[dataIndices[centerBoxIndex + zIndexOffset]]
                                                          : zeroVector;

                    const vector3f& centerXnegVoxel =
                        ( dataIndices[centerBoxIndex - 1] >= 0 ) ? vector[dataIndices[centerBoxIndex - 1]] : zeroVector;
                    const vector3f& centerYnegVoxel = ( dataIndices[centerBoxIndex - yIndexOffset] >= 0 )
                                                          ? vector[dataIndices[centerBoxIndex - yIndexOffset]]
                                                          : zeroVector;
                    const vector3f& centerZnegVoxel = ( dataIndices[centerBoxIndex - zIndexOffset] >= 0 )
                                                          ? vector[dataIndices[centerBoxIndex - zIndexOffset]]
                                                          : zeroVector;

                    // const vector3f& centerXnegYnegVoxel = (dataIndices[centerBoxIndex - yIndexOffset-1] >= 0 ) ?
                    // vector[dataIndices[centerBoxIndex - yIndexOffset-1]] : zeroVector;
                    // const vector3f& centerXnegZnegVoxel = (dataIndices[centerBoxIndex - zIndexOffset-1] >= 0 ) ?
                    // vector[dataIndices[centerBoxIndex - zIndexOffset-1]] : zeroVector;
                    const vector3f& centerXnegYposVoxel = ( dataIndices[centerBoxIndex + yIndexOffset - 1] >= 0 )
                                                              ? vector[dataIndices[centerBoxIndex + yIndexOffset - 1]]
                                                              : zeroVector;
                    const vector3f& centerXnegZposVoxel = ( dataIndices[centerBoxIndex + zIndexOffset - 1] >= 0 )
                                                              ? vector[dataIndices[centerBoxIndex + zIndexOffset - 1]]
                                                              : zeroVector;

                    const vector3f& centerXposYnegVoxel = ( dataIndices[centerBoxIndex - yIndexOffset + 1] >= 0 )
                                                              ? vector[dataIndices[centerBoxIndex - yIndexOffset + 1]]
                                                              : zeroVector;
                    const vector3f& centerXposZnegVoxel = ( dataIndices[centerBoxIndex - zIndexOffset + 1] >= 0 )
                                                              ? vector[dataIndices[centerBoxIndex - zIndexOffset + 1]]
                                                              : zeroVector;
                    const vector3f& centerXposYposVoxel = ( dataIndices[centerBoxIndex + yIndexOffset + 1] >= 0 )
                                                              ? vector[dataIndices[centerBoxIndex + yIndexOffset + 1]]
                                                              : zeroVector;
                    const vector3f& centerXposZposVoxel = ( dataIndices[centerBoxIndex + zIndexOffset + 1] >= 0 )
                                                              ? vector[dataIndices[centerBoxIndex + zIndexOffset + 1]]
                                                              : zeroVector;

                    // const vector3f& centerYnegZnegVoxel = (dataIndices[centerBoxIndex - yIndexOffset - zIndexOffset]
                    // >= 0 ) ? vector[dataIndices[centerBoxIndex - yIndexOffset - zIndexOffset]] : zeroVector;
                    const vector3f& centerYnegZposVoxel =
                        ( dataIndices[centerBoxIndex - yIndexOffset + zIndexOffset] >= 0 )
                            ? vector[dataIndices[centerBoxIndex - yIndexOffset + zIndexOffset]]
                            : zeroVector;
                    const vector3f& centerYposZnegVoxel =
                        ( dataIndices[centerBoxIndex + yIndexOffset - zIndexOffset] >= 0 )
                            ? vector[dataIndices[centerBoxIndex + yIndexOffset - zIndexOffset]]
                            : zeroVector;
                    const vector3f& centerYposZposVoxel =
                        ( dataIndices[centerBoxIndex + yIndexOffset + zIndexOffset] >= 0 )
                            ? vector[dataIndices[centerBoxIndex + yIndexOffset + zIndexOffset]]
                            : zeroVector;

                    // const vector3f& centerXnegYnegZnegVoxel = (dataIndices[centerBoxIndex - yIndexOffset -
                    // zIndexOffset-1] >= 0 ) ? vector[dataIndices[centerBoxIndex - yIndexOffset - zIndexOffset-1]] :
                    // zeroVector; const vector3f& centerXnegYposZnegVoxel = (dataIndices[centerBoxIndex + yIndexOffset
                    // - zIndexOffset-1] >= 0 ) ? vector[dataIndices[centerBoxIndex + yIndexOffset - zIndexOffset-1]] :
                    // zeroVector; const vector3f& centerXnegYnegZposVoxel = (dataIndices[centerBoxIndex - yIndexOffset
                    // + zIndexOffset-1] >= 0 ) ? vector[dataIndices[centerBoxIndex - yIndexOffset + zIndexOffset-1]] :
                    // zeroVector; const vector3f& centerXnegYposZposVoxel = (dataIndices[centerBoxIndex + yIndexOffset
                    // + zIndexOffset-1] >= 0 ) ? vector[dataIndices[centerBoxIndex + yIndexOffset + zIndexOffset-1]] :
                    // zeroVector;

                    // const vector3f& centerXposYnegZnegVoxel = (dataIndices[centerBoxIndex - yIndexOffset -
                    // zIndexOffset+1] >= 0 ) ? vector[dataIndices[centerBoxIndex - yIndexOffset - zIndexOffset+1]] :
                    // zeroVector; const vector3f& centerXposYposZnegVoxel = (dataIndices[centerBoxIndex + yIndexOffset
                    // - zIndexOffset+1] >= 0 ) ? vector[dataIndices[centerBoxIndex + yIndexOffset - zIndexOffset+1]] :
                    // zeroVector; const vector3f& centerXposYnegZposVoxel = (dataIndices[centerBoxIndex - yIndexOffset
                    // + zIndexOffset+1] >= 0 ) ? vector[dataIndices[centerBoxIndex - yIndexOffset + zIndexOffset+1]] :
                    // zeroVector; const vector3f& centerXposYposZposVoxel = (dataIndices[centerBoxIndex + yIndexOffset
                    // + zIndexOffset+1] >= 0 ) ? vector[dataIndices[centerBoxIndex + yIndexOffset + zIndexOffset+1]] :
                    // zeroVector;

                    vector3f& Du = derivative[dataIndices[centerBoxIndex]];

                    float dudy =
                        ( centerYposVoxel.x + centerXposYposVoxel.x - centerYnegVoxel.x - centerXposYnegVoxel.x );
                    float dudz =
                        ( centerZposVoxel.x + centerXposZposVoxel.x - centerZnegVoxel.x - centerXposZnegVoxel.x );

                    Du.x = ( ( centerXposVoxel.x - centerVoxel.x ) + 0.25f * ( dudy + dudz ) ) / dx;

                    float dvdx =
                        ( centerXposVoxel.y + centerXposYposVoxel.y - centerXnegVoxel.y - centerXnegYposVoxel.y );
                    float dvdz =
                        ( centerZposVoxel.y + centerYposZposVoxel.y - centerZnegVoxel.y - centerYposZnegVoxel.y );

                    Du.y = ( ( centerYposVoxel.y - centerVoxel.y ) + 0.25f * ( dvdx + dvdz ) ) / dx;

                    float dwdx =
                        ( centerXposVoxel.z + centerXposZposVoxel.z - centerXnegVoxel.z - centerXnegZposVoxel.z );
                    float dwdy =
                        ( centerYposVoxel.z + centerYposZposVoxel.z - centerYnegVoxel.z - centerYnegZposVoxel.z );

                    Du.z = ( ( centerZposVoxel.z - centerVoxel.z ) + 0.25f * ( dwdx + dwdy ) ) / dx;

                    // if( !Du.is_zero() )
                    //	std::cout << "deriv: " << Du << endl;
                }
            }
        }
    }
}

class rle_scanline_velocity_correct_and_clamp {
    const rle_index_spec& m_ris;
    const const_rle_channel_accessor<vector3f>& m_resultVel;
    const const_rle_channel_accessor<vector3f>& m_origVel;
    const const_rle_channel_accessor<vector3f>& m_errorEstimateVel;
    rle_channel_accessor<vector3f>& m_outVel;
    const vector3f* m_advectedVelChannel;

    rle_scanline_velocity_correct_and_clamp&
    operator=( const rle_scanline_velocity_correct_and_clamp& ); // not implemented

  public:
    rle_scanline_velocity_correct_and_clamp( const rle_index_spec& ris,
                                             const const_rle_channel_accessor<vector3f>& resultVel,
                                             const const_rle_channel_accessor<vector3f>& origVel,
                                             const const_rle_channel_accessor<vector3f>& errorEstimateVel,
                                             rle_channel_accessor<vector3f>& outVel,
                                             const vector3f* advectedVelChannel )
        : m_ris( ris )
        , m_resultVel( resultVel )
        , m_origVel( origVel )
        , m_errorEstimateVel( errorEstimateVel )
        , m_outVel( outVel )
        , m_advectedVelChannel( advectedVelChannel ) {}

    // nb : using int xyz range
    void operator()( const tbb::blocked_range2d<int>& r ) const {
        const boundbox3 bounds = m_ris.outer_bounds();
        const vector3& boundsMin = bounds.minimum();

        const int ysize = bounds.ysize();

        const size3 boxSize( 3 );
        const vector3 centerOffset( 1 );
        const int cellBoxIndex =
            centerOffset.x + centerOffset.y * boxSize.xsize() + centerOffset.z * boxSize.xsize() * boxSize.ysize();

        const std::vector<boost::int32_t>& bcToRunIndex = m_ris.get_bc_to_run_index_vector();
        const std::vector<run_data>& runIndexData = m_ris.get_run_index_data_vector();

        for( int z = r.rows().begin(); z != r.rows().end(); ++z ) {
            int cIndex = m_ris.z_to_c( z );
            for( int y = r.cols().begin(); y != r.cols().end(); ++y ) {
                // compute the bcIndex of the scanline
                int bcIndex = m_ris.y_to_b( y ) + cIndex * ysize;

                // logging::error << "Computing Scanline: [" << y << ", " << z << "]: bcIndex:" << bcIndex << endl;
                // get the run range of the scanline
                int runRangeStart = bcToRunIndex[bcIndex];
                int runRangeEnd = bcToRunIndex[bcIndex + 1] - 2;

                // provided this scanline has at least 2 real runs (and therefore likely some defined data) process the
                // runs with the iterator
                if( runRangeStart != runRangeEnd ||
                    runIndexData[runRangeStart].x != runIndexData[runRangeStart + 1].x ) {
                    // construct a block iterator for this scanline

                    // for velocity advection we need one voxel in the negative direction, and two in hte positive
                    // to make sure we have all the data indices that we will need to process the lookup. This means we
                    // will end up with more overall indices, but should be faster than carefully choosing and finding
                    // only the ones we need
                    boundbox3 box( vector3( boundsMin.x - centerOffset.x, y - centerOffset.y, z - centerOffset.z ),
                                   boxSize );

                    for( levelset::rle_defined_box_iterator boxIter( m_ris, box );
                         !boxIter.is_xplane_finished( centerOffset.x ); ++boxIter ) {

                        const boost::int32_t* const dataIndices = boxIter.get_indices();
                        const boost::int32_t cellIndex = dataIndices[cellBoxIndex];

                        if( cellIndex < 0 ) {
                            continue;
                        }

                        // It is generally better to clamp the velocity to the local neighborhood to avoid overshooting
                        // and creating suprious spikes in the velocity field.
                        const vector3f vel =
                            m_resultVel[cellIndex] - ( m_errorEstimateVel[cellIndex] - m_origVel[cellIndex] ) * 0.5f;

                        boundbox3f clampBounds;
                        vector3f& min = clampBounds.minimum();
                        vector3f& max = clampBounds.maximum();

                        for( int face = 0; face < 3; ++face ) {
                            // clamp the velocity the original neighborhood if we have an overshoot
                            const vector3 abc = centerOffset - vector3( 1 ) + vector3::from_axis( face );

                            for( int k = 0; k < 2; ++k ) {
                                int zIndexOffset = ( abc.z + k ) * boxSize.xsize() * boxSize.ysize();
                                for( int j = 0; j < 2; ++j ) {
                                    int yzIndexOffset = ( abc.y + j ) * boxSize.xsize() + zIndexOffset;
                                    for( int i = 0; i < 2; ++i ) {
                                        int boxIndex = ( abc.x + i ) + yzIndexOffset;
                                        boost::int32_t dataIndex = dataIndices[boxIndex];

                                        if( dataIndex < 0 ) {
                                            continue;
                                        }

                                        if( m_origVel[dataIndex][face] < min[face] )
                                            min[face] = m_origVel[dataIndex][face];
                                        if( m_origVel[dataIndex][face] > max[face] )
                                            max[face] = m_origVel[dataIndex][face];
                                    }
                                }
                            }
                        }

                        if( m_advectedVelChannel ) {
                            vector3f newVelocity;
                            for( int axis = 0; axis < 3; ++axis ) {
                                if( min[axis] <= vel[axis] && vel[axis] <= max[axis] ) {
                                    // use error-corrected velocity
                                    newVelocity[axis] = vel[axis];
                                } else {
                                    // clamp to semi-lagrangian advection
                                    newVelocity[axis] = m_advectedVelChannel[cellIndex][axis];
                                }
                            }
                            m_outVel[cellIndex] = newVelocity;
                        } else {
                            m_outVel[cellIndex] = clampBounds.clamp( vel );
                        }
                    }
                }
            }
        }
    }
};

/**
 *  Perform the velocity correction and clamping used in
 * velocity_maccormack_advect and velocity_bfecc_advect.
 *
 *  Note: the parameter names are based on those used in
 * velocity_maccormack_advect -- they do not match those used
 * in BFECC !
 *
 *  The velocity correction calculates:
 *
 *    outVel = resultVel + 0.5 * ( origVel - errorEstimateVel )
 *
 *  For BFECC, we pass in resultVel=origVel to get:
 *
 *    outVel = origVel + 0.5 * ( origVel - errorEstimateVel )
 *
 *  New extrema in outVel are detected by comparing each staggered sample
 * of outVel with the 8 surrounding staggered samples in origVel.
 * (TODO: should this be the advection's source point instead of the
 * destination point?)
 *
 *  How extrema are corrected depends on whether you provide
 * advectedVelChannel.  If advectedVelChannel is given, then the extrema
 * are replaced with the staggered sample value from that field.
 * Usually advectedVelChannel is the result of semi-lagrangian
 * advection, which is currently named resultVel in the calling functions.
 * If advectedVelAcc is NULL, then new extrema are clamped to lie within
 * the range of values found in the 8 surrounding staggered samples.
 *
 */
void velocity_correct_and_clamp( const rle_index_spec& ris, const const_rle_channel_accessor<vector3f>& resultVel,
                                 const const_rle_channel_accessor<vector3f>& origVel,
                                 const const_rle_channel_accessor<vector3f>& errorEstimateVel,
                                 rle_channel_accessor<vector3f>& outVel, const vector3f* advectedVelChannel = 0 ) {
    tbb::task_scheduler_init taskScheduleInit;

    const boundbox3 bounds = ris.outer_bounds();
    const vector3& boundsMin = bounds.minimum();
    const vector3& boundsMax = bounds.maximum();

    rle_scanline_velocity_correct_and_clamp body( ris, resultVel, origVel, errorEstimateVel, outVel,
                                                  advectedVelChannel );

    tbb::blocked_range2d<int> range( boundsMin.z, boundsMax.z + 1, boundsMin.y, boundsMax.y + 1 );

#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, body, tbb::auto_partitioner() );
#else
#pragma message( "Threads are disabled" )
    body( range );
#endif
}

void velocity_maccormack_advect( const rle_voxel_field& advectedField, rle_voxel_field& result,
                                 const rle_voxel_field& srcField, float dt ) {
    const float maxSrcVoxelMotion = detail::get_maximum_voxel_motion_from_staggered_velocity( srcField, dt );

    rle_voxel_field errorEstimateField( advectedField.get_voxel_coord_system() );

    frantic::fluids::velocity_advect( advectedField, result, srcField, maxSrcVoxelMotion, dt );

    frantic::fluids::velocity_advect( result, errorEstimateField, srcField, maxSrcVoxelMotion, -dt );

    rle_channel_accessor<vector3f> resultVel = result.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );
    const_rle_channel_accessor<vector3f> origVel =
        advectedField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );
    rle_channel_accessor<vector3f> errorEstimateVel =
        errorEstimateField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    const rle_index_spec& ris = result.get_rle_index_spec();

    velocity_correct_and_clamp( ris, resultVel, origVel, errorEstimateVel, resultVel,
                                reinterpret_cast<vector3f*>( resultVel.data( 0 ) ) );
}

/**
 *  Perform velocity advection using RLE lookups to sample the velocity field.
 * rle_scanline_window_advect is faster for lesser advection motion.
 */
class rle_scanline_lookup_advect_and_revert {
    const boundbox3& m_bounds;
    const rle_index_spec& m_ris;
    const voxel_coord_system& m_vcs;

    const_rle_channel_accessor<vector3f>& m_sourceAcc;
    const_rle_channel_accessor<vector3f>& m_revertAcc;
    const_rle_channel_accessor<vector3f>& m_velocityAcc;
    rle_channel_accessor<vector3f>& m_resultAcc;

    float m_dt;

    rle_scanline_lookup_advect_and_revert& operator=( const rle_scanline_lookup_advect_and_revert& ); // not implemented

  public:
    rle_scanline_lookup_advect_and_revert( const boundbox3& bounds, const rle_index_spec& ris,
                                           const voxel_coord_system& vcs,
                                           const_rle_channel_accessor<vector3f>& sourceAcc,
                                           const_rle_channel_accessor<vector3f>& revertAcc,
                                           const_rle_channel_accessor<vector3f>& velocityAcc,
                                           rle_channel_accessor<vector3f>& resultAcc, float dt )
        : m_bounds( bounds )
        , m_ris( ris )
        , m_vcs( vcs )
        , m_sourceAcc( sourceAcc )
        , m_revertAcc( revertAcc )
        , m_velocityAcc( velocityAcc )
        , m_resultAcc( resultAcc )
        , m_dt( dt ) {}

    void operator()( const tbb::blocked_range2d<size_t>& r ) const {

        const vector3& boundsMin = m_bounds.minimum();

        cached_trilerp_vector3f_staggered velLookup( m_ris );

        for( size_t c = r.rows().begin(); c != r.rows().end(); ++c ) {
            int z = static_cast<int>( c ) + boundsMin.z;
            for( size_t b = r.cols().begin(); b != r.cols().end(); ++b ) {
                int y = static_cast<int>( b ) + boundsMin.y;

                for( rle_run_iterator i( m_ris, static_cast<int>( b ), static_cast<int>( c ) ), ie; i != ie; ++i ) {
                    boost::int32_t dataIndex = i.get_data_index();
                    if( dataIndex < 0 ) {
                        continue;
                    }
                    for( boost::int32_t x = i.get_xmin(); x <= i.get_xmax(); ++x, ++dataIndex ) {
                        const vector3 cellVoxelCoord( x, y, z );
                        vector3f newVelocity;

                        // get the current lookup velocity
                        for( int comp = 0; comp < 3; ++comp ) {
                            // shift to the real world location of the cell sample
                            // we are using only the x,y,z faces so we need to double the index
                            const vector3f faceVoxelCoord = m_vcs.get_voxel_face_center( cellVoxelCoord, comp * 2 );
                            const vector3f faceWorldCoord = m_vcs.get_world_face_center( cellVoxelCoord, comp * 2 );

                            vector3f vel;
                            // TODO: don't silently ignore failed lookups
                            velLookup.get( m_velocityAcc, faceVoxelCoord, vel );

                            const vector3f midpoint = faceWorldCoord - 0.5f * vel * m_dt;
                            vector3f vel2;
                            // TODO: don't silently ignore failed lookups
                            velLookup.get( m_velocityAcc, m_vcs.get_voxel_coord( midpoint ), vel2 );

                            const vector3f quarterpoint = faceWorldCoord - 0.75f * vel2 * m_dt;
                            vector3f vel3;
                            // TODO: don't silently ignore failed lookups
                            velLookup.get( m_velocityAcc, m_vcs.get_voxel_coord( quarterpoint ), vel3 );

                            const vector3f finalPosition =
                                faceWorldCoord - ( 2 * vel + 3 * vel2 + 4 * vel3 ) * m_dt / 9.f;
                            vector3f finalPosVelocity;
                            // TODO: don't silently ignore failed lookups
                            velLookup.get( m_sourceAcc, m_vcs.get_voxel_coord( finalPosition ), finalPosVelocity );

                            const boundbox3f clampRange =
                                velLookup.get_sample_range( m_revertAcc, m_vcs.get_voxel_coord( finalPosition ) );
                            if( clampRange.minimum()[comp] <= clampRange.maximum()[comp] ) {
                                if( finalPosVelocity[comp] < clampRange.minimum()[comp] ||
                                    finalPosVelocity[comp] > clampRange.maximum()[comp] ) {
                                    // TODO: don't silently ignore failed lookups
                                    velLookup.get( m_revertAcc, m_vcs.get_voxel_coord( finalPosition ),
                                                   finalPosVelocity );
                                }
                            }
                            newVelocity[comp] = finalPosVelocity[comp];
                        }
                        m_resultAcc[dataIndex] = newVelocity;
                    }
                }
            }
        }
    }
};

void velocity_advect_and_revert_extrema( const rle_voxel_field& source, const rle_voxel_field& revert,
                                         rle_voxel_field& result, const rle_voxel_field& velocityField, float dt ) {
    // TODO : The velocity advection routines can use undefined source
    // regions, which will result in erroneous 0 velocity.  We should
    // use the nearest defined voxels or extrapolated values instead.
    //
    // This same change should be made in velocity lookup part of
    // level set advection.
    tbb::task_scheduler_init taskScheduleInit;

    const frantic::volumetrics::voxel_coord_system& vcs = source.get_voxel_coord_system();
    const rle_index_spec& ris = source.get_rle_index_spec();
    const rle_index_spec& revertRIS = revert.get_rle_index_spec();
    const rle_index_spec& velocityFieldRIS = velocityField.get_rle_index_spec();

    const bool sameRleIndexSpecs = ( ris == revertRIS ) && ( ris == velocityFieldRIS );
    if( !sameRleIndexSpecs ) {
        throw std::runtime_error( "velocity_advect_and_revert_extrema Error: the rle index specs are different." );
    }

    if( vcs != revert.get_voxel_coord_system() ) {
        throw std::runtime_error(
            "velocity_advect_and_revert_extrema Error: the source and revert voxel coordinate systems are different." );
    }
    if( vcs != velocityField.get_voxel_coord_system() ) {
        throw std::runtime_error(
            "velocity_advect_and_revert_extrema Error: the source and velocityField voxel coordinate "
            "systems are different." );
    }

    const_rle_channel_accessor<vector3f> sourceAcc = source.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );
    const_rle_channel_accessor<vector3f> revertAcc = revert.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );
    const_rle_channel_accessor<vector3f> velocityAcc =
        velocityField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    rle_voxel_field temp( vcs, ris );
    temp.add_channel<vector3f>( _T("StaggeredVelocity") );
    // temp.zero_channel( "StaggeredVelocity" );
    rle_channel_accessor<vector3f> resultAcc = temp.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    // logging::set_logging_level( 3 );

    FF_LOG( error ) << "**** (vel advect )staggered velocity channel size=" << sourceAcc.size() << " ****" << endl;
    FF_LOG( debug ) << "Starting iteration" << std::endl;

    const boundbox3 bounds = ris.outer_bounds();

    // run the tbb:: advector object in parallel (in theory)

    tbb::blocked_range2d<size_t> range( 0, bounds.zsize(), 0, bounds.ysize() );

    logging::error << "threading block range: " << bounds.ysize() << ", " << bounds.zsize() << endl;

    typedef rle_scanline_lookup_advect_and_revert advector_t;

    advector_t parallelAdvector( bounds, ris, vcs, sourceAcc, revertAcc, velocityAcc, resultAcc, dt );

#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_for( range, parallelAdvector, tbb::auto_partitioner() );
#else
#pragma message( "Threads are disabled" )
    parallelAdvector( range );
#endif

    result.swap( temp );
}

void velocity_bfecc_advect( const rle_voxel_field& advectedField, rle_voxel_field& result,
                            const rle_voxel_field& srcField, float dt ) {
    // TODO: revert when advecting from occlusions or sim bounds
    const float maxSrcVoxelMotion = detail::get_maximum_voxel_motion_from_staggered_velocity( srcField, dt );

    rle_voxel_field errorEstimateField( advectedField.get_voxel_coord_system() );

    frantic::fluids::velocity_advect( advectedField, result, srcField, maxSrcVoxelMotion, dt );

    frantic::fluids::velocity_advect( result, errorEstimateField, srcField, maxSrcVoxelMotion, -dt );

    const_rle_channel_accessor<vector3f> origVel =
        advectedField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );
    rle_channel_accessor<vector3f> errorEstimateVel =
        errorEstimateField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    const rle_index_spec& ris = result.get_rle_index_spec();

    // \overline{\phi}^n = ( 3 \phi^n - \hat{\phi}^n)/2
    for( std::size_t i = 0; i < ris.data_size(); ++i ) {
        errorEstimateVel[i] = 1.5f * origVel[i] - 0.5f * errorEstimateVel[i];
    }

    frantic::fluids::velocity_advect_and_revert_extrema( errorEstimateField, advectedField, result, srcField, dt );
}

void velocity_modified_bfecc_advect( const rle_voxel_field& advectedField, rle_voxel_field& result,
                                     const rle_voxel_field& srcField, float dt ) {
    const float maxSrcVoxelMotion = detail::get_maximum_voxel_motion_from_staggered_velocity( srcField, dt );

    rle_voxel_field errorEstimateField( advectedField.get_voxel_coord_system() );

    frantic::fluids::velocity_advect( advectedField, result, srcField, maxSrcVoxelMotion, dt );

    frantic::fluids::velocity_advect( result, errorEstimateField, srcField, maxSrcVoxelMotion, -dt );

    rle_channel_accessor<vector3f> resultVel = result.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );
    const_rle_channel_accessor<vector3f> origVel =
        advectedField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );
    rle_channel_accessor<vector3f> errorEstimateVel =
        errorEstimateField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    const rle_index_spec& ris = result.get_rle_index_spec();

    velocity_correct_and_clamp( ris, resultVel, origVel, errorEstimateVel, errorEstimateVel );

    frantic::fluids::velocity_advect( errorEstimateField, result, srcField, maxSrcVoxelMotion, dt );
}
} // namespace fluids
} // namespace frantic
