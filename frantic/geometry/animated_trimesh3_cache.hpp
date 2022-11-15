// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/files/files.hpp>
#include <frantic/files/paths.hpp>
#include <frantic/geometry/trimesh3_file_io.hpp>

#include <frantic/files/file_sequence.hpp>

#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace geometry {

//  This class is designed to handle a number of use cases, and should satisfy the following requirements for now.
//
// API Support:
// - Be able to use and generate 3ds max Mesh objects for fast operations within a max plugin.
// - Be able to use and generate trimesh3 objects for use outside of max.
//
// FILE Format Support:
// - Support sequences of mesh files with consistent vertex counts.
// - Support a single mesh file that remains consistent across time.
// - Support file numbering based on frames or subframe ticks.
//
// What it doesn't do:
// - support animated transformations.  This is only caching the underlying triangle mesh.  The animation of
// transformations
//   should happen outside of this.  For example, in 3ds max, transform animations are done at the inode level rather
//   than the geometry object level.

class animated_trimesh3_cache {
    // This parameter is 1 for normal caches, and will vary based on the framerate for ticks-based caches.
    int m_countsPerFrame;
    // The frames that are cached to a file should be at times
    // m_countIntervalStart, m_countIntervalStart + m_countsPerFrame, m_countIntervalStart + 2*m_countsPerFrame, ...,
    // m_countIntervalEnd
    int m_countIntervalStart, m_countIntervalEnd;
    // This offset is added to file lookups
    int m_countOffset;
    // Adjusts the speed of the mesh animation.  This will cause the system to interpolate frames.
    float m_timeSpeedFactor;
    // This is the padding of the numbered frames
    int m_framePadding;
    // Whether the mesh cache is animated.  This may be determined by whether the filename has numbers just before the
    // extension, for example.
    bool m_animated;
    // The files always look like m_filePrefix + strings::zero_pad(count, m_framePadding) + m_filePostfix
    frantic::tstring m_filePrefix, m_filePostfix;

  public:
    ///////////////
    // Constructors
    ///////////////

    animated_trimesh3_cache() {
        m_countsPerFrame = 1;
        m_countOffset = 0;
        m_framePadding = 4;
        m_animated = false;
        m_countIntervalStart = 0;
        m_countIntervalEnd = 0;
        m_timeSpeedFactor = 1.f;
    }

    // This constructs a mesh cache from an input file.  Whether the sequence is animated is determined by the existence
    // of a sequence number in the file name.
    animated_trimesh3_cache( const frantic::tstring& filename, int countsPerFrame = 1, int countOffset = 0 ) {
        set( filename, countsPerFrame, countOffset );
    }

    ///////////////
    // Setters
    ///////////////

    void set( const frantic::tstring& filename, int countsPerFrame = 1, int countOffset = 0 ) {
        // Set the initial parameters
        m_countsPerFrame = countsPerFrame;
        m_countOffset = countOffset;
        m_timeSpeedFactor = 1.f;

        // Split apart the prototype filename
        frantic::files::split_sequence_path( filename, m_filePrefix, m_framePadding, m_countIntervalStart,
                                             m_filePostfix );
        if( m_framePadding > 0 ) {
            // Adjust the number which was in the file to account for the count offset
            std::pair<int, int> ivalid = frantic::files::get_sequence_range( frantic::strings::to_tstring( filename ) );
            m_countIntervalStart = ivalid.first;
            m_countIntervalEnd = ivalid.second;

            // m_countIntervalStart += countOffset;
            // m_countIntervalEnd = m_countIntervalStart;
            // expand_count_interval();	<--Deprecated

            m_animated = true;
        } else {
            m_animated = false;
            m_countIntervalStart = 0;
            m_countIntervalEnd = 0;
        }
    }

    void set( const frantic::tstring& filePrefix, int framePadding, const frantic::tstring& filePostfix,
              int countsPerFrame, int countOffset, bool animated, float timeSpeedFactor, int countIntervalStart,
              int countIntervalEnd ) {
        m_filePrefix = filePrefix;
        m_framePadding = framePadding;
        m_filePostfix = filePostfix;
        m_countsPerFrame = countsPerFrame;
        m_countOffset = countOffset;
        m_animated = animated;
        m_timeSpeedFactor = timeSpeedFactor;
        m_countIntervalStart = countIntervalStart;
        m_countIntervalEnd = countIntervalEnd;
    }

    void set_filePrefix( const frantic::tstring& filePrefix ) { m_filePrefix = filePrefix; }

    void set_filePostfix( const frantic::tstring& filePostfix ) { m_filePostfix = filePostfix; }

    void set_framePadding( int framePadding ) { m_framePadding = framePadding; }

    void set_countOffset( int countOffset ) { m_countOffset = countOffset; }

    void set_timeSpeedFactor( float timeSpeedFactor ) { m_timeSpeedFactor = timeSpeedFactor; }

    void set_animated( bool animated ) { m_animated = animated; }

    void set_countInterval( int countStart, int countEnd ) {
        m_countIntervalStart = countStart;
        m_countIntervalEnd = countEnd;
    }

    ///////////////
    // Getters
    ///////////////

    const frantic::tstring& get_filePrefix() const { return m_filePrefix; }

    const frantic::tstring& get_filePostfix() const { return m_filePostfix; }

    int get_framePadding() const { return m_framePadding; }

    int get_countsPerFrame() const { return m_countsPerFrame; }

    int get_countOffset() const { return m_countOffset; }

    bool get_animated() const { return m_animated; }

    float get_timeSpeedFactor() const { return m_timeSpeedFactor; }

    std::pair<int, int> get_countInterval() const { return std::make_pair( m_countIntervalStart, m_countIntervalEnd ); }

    frantic::tstring get_nearest_mesh_filename( int count ) const {
        float cacheCountNumber = to_cache_count_number( (float)count );
        int nearestCacheCountNumber = get_nearest_cache_count_number( cacheCountNumber );
        return get_mesh_filename( nearestCacheCountNumber );
    }

    frantic::tstring get_save_mesh_filename( int count ) const {
        if( m_timeSpeedFactor != 1.f )
            throw std::runtime_error(
                "animated_trimesh3_cache.get_save_mesh_filename: Can only save meshes when the time "
                "speed is set to 1.0.  It is currently set to " +
                boost::lexical_cast<std::string>( m_timeSpeedFactor ) );

        return get_mesh_filename( to_cache_count_number( count ) );
    }

    ///////////////
    // Operations
    ///////////////

    // Looks at the files on disk, adjusting the interval to fit the largest sequence of files that match the file
    // pattern
    void adjust_count_interval() {
        if( get_animated() ) {
            std::pair<int, int> ivalid =
                frantic::files::get_sequence_range( frantic::strings::to_tstring( get_mesh_filename( 0 ) ) );
            m_countIntervalStart = ivalid.first;
            m_countIntervalEnd = ivalid.second;
        }
    }

    // This will scan the cache's data directory and return the interval of the lowest frame found to the highest
    void adjust_to_largest_interval() {
        const std::pair<int, int> cInterval = files::get_sequence_range(
            frantic::strings::to_tstring( m_filePrefix ) + strings::zero_pad( 0, m_framePadding ) +
            frantic::strings::to_tstring( m_filePostfix ) );
        m_countIntervalStart = cInterval.first;
        m_countIntervalEnd = cInterval.second;
    }

    void expand_count_interval() {
        while( files::file_exists(
            get_mesh_filename( to_cache_count_number( m_countIntervalStart - m_countsPerFrame ) ) ) )
            m_countIntervalStart -= m_countsPerFrame;

        while(
            files::file_exists( get_mesh_filename( to_cache_count_number( m_countIntervalEnd + m_countsPerFrame ) ) ) )
            m_countIntervalEnd += m_countsPerFrame;
    }

    void load_mesh_nearest( int count, trimesh3& outMesh ) { return load_mesh_nearest( (float)count, outMesh ); }

    void load_mesh_nearest( float count, trimesh3& outMesh ) {
        float cacheCountNumber = to_cache_count_number( count );
        int nearestCacheCountNumber = get_nearest_cache_count_number( cacheCountNumber );
        load_mesh_internal( nearestCacheCountNumber, outMesh );
    }

    void load_mesh_interpolated( int count, trimesh3& outMesh ) {
        return load_mesh_interpolated( (float)count, outMesh );
    }

    void load_mesh_interpolated( float count, trimesh3& outMesh ) {
        float cacheCountNumber = to_cache_count_number( count );
        int belowCacheCountNumber = 0, aboveCacheCountNumber = 0;
        get_nearest_cache_count_interval( cacheCountNumber, belowCacheCountNumber, aboveCacheCountNumber );
        if( belowCacheCountNumber == aboveCacheCountNumber ) {
            load_mesh_internal( belowCacheCountNumber, outMesh );
        } else {
            trimesh3 m1, m2;
            load_mesh_internal( belowCacheCountNumber, m1 );
            load_mesh_internal( aboveCacheCountNumber, m2 );
            // TODO: Write this from_linear_interpolation function
            // float alpha = (float)(cacheCountNumber - belowCacheCountNumber) / (float)(aboveCacheCountNumber -
            // belowCacheCountNumber); outMesh.from_linear_interpolation( m1, m2, alpha );
            throw std::runtime_error( "TODO: trimesh3::load_mesh_interpolated isn't implemented yet." );
        }
    }

    void write_mesh( int count, trimesh3& mesh ) {
        if( m_timeSpeedFactor != 1.f )
            throw std::runtime_error(
                "animated_trimesh3_cache.write_mesh: Can only save meshes when the time speed is set "
                "to 1.0.  It is currently set to " +
                boost::lexical_cast<std::string>( m_timeSpeedFactor ) );

        int cacheCountNumber = to_cache_count_number( count );
        write_mesh_internal( cacheCountNumber, mesh );
    }

  private:
    int to_cache_count_number( int count ) const {
        float cacheCountNumber =
            m_countIntervalStart + ( count - m_countIntervalStart ) * m_timeSpeedFactor - m_countOffset;
        if( cacheCountNumber > m_countIntervalEnd )
            cacheCountNumber = (float)m_countIntervalEnd;

        if( cacheCountNumber < m_countIntervalStart )
            cacheCountNumber = (float)m_countIntervalStart;

        return (int)( frantic::math::round( cacheCountNumber ) );
    }

    // Previous algorithm was BUGGZILLA. It mostly just scaled the whole time segment, not the active cache segment.
    float to_cache_count_number( float count ) const {
        float cacheCountNumber =
            m_countIntervalStart + ( count - m_countIntervalStart ) * m_timeSpeedFactor - m_countOffset;
        if( cacheCountNumber > m_countIntervalEnd ) {
            // TODO: Commented out these lines because it seems stupid to have side effects in this function!
            //       Should figure out why it's here, and fix any consequences removing the side effects causes
            // So this frame isn't in our currently known interval. I'll check if its file exists
            // then assign it to the interval endpoint if it does.
            //	if( files::file_exists( get_mesh_filename( (int)cacheCountNumber ) ) )
            //		m_countIntervalEnd = (int)cacheCountNumber;

            cacheCountNumber = (float)m_countIntervalEnd;
        }

        if( cacheCountNumber < m_countIntervalStart ) {
            // TODO: Commented out these lines because it seems stupid to have side effects in this function!
            //       Should figure out why it's here, and fix any consequences removing the side effects causes
            // So this frame isn't in our currently known interval. I'll check if its file exists
            // then assign it to the interval endpoint if it does.
            // if( files::file_exists( get_mesh_filename( (int)cacheCountNumber ) ) )
            //	m_countIntervalStart = (int)cacheCountNumber;

            cacheCountNumber = (float)m_countIntervalStart;
        }

        return cacheCountNumber;
    }

    void shrink_count_interval() {
        if( m_timeSpeedFactor != 1.f )
            throw std::runtime_error( "animated_trimesh3_cache.expand_count_interval: Can only shrink the source file "
                                      "interval when the time factor is 1.0.  It is currently set to " +
                                      boost::lexical_cast<std::string>( m_timeSpeedFactor ) );

        int savedStart = m_countIntervalStart, savedEnd = m_countIntervalEnd;

        while( m_countIntervalStart <= m_countIntervalEnd &&
               !files::file_exists( get_mesh_filename( to_cache_count_number( m_countIntervalStart ) ) ) )
            m_countIntervalStart += m_countsPerFrame;

        while( m_countIntervalStart <= m_countIntervalEnd &&
               !files::file_exists( get_mesh_filename( to_cache_count_number( m_countIntervalEnd ) ) ) )
            m_countIntervalEnd -= m_countsPerFrame;

        // If we shrank it to empty, restore the interval
        if( m_countIntervalStart > m_countIntervalEnd ) {
            m_countIntervalStart = savedStart;
            m_countIntervalEnd = savedEnd;
        }
    }

    int get_nearest_cache_count_number( float cacheCountNumber ) const {
        int startAbove = (int)ceil( cacheCountNumber );
        int startBelow = (int)floor( cacheCountNumber );
        if( startAbove == startBelow ) { // If the number is an integer, take the potential faster route of just a
                                         // single file existence test
            if( files::file_exists( get_mesh_filename( startBelow ) ) )
                return startBelow;
            startAbove++;
            startBelow--;
        } else if( startAbove - cacheCountNumber <
                   cacheCountNumber - startBelow ) { // If the interval above is smaller, start the search above.
            if( files::file_exists( get_mesh_filename( startAbove ) ) )
                return startAbove;
            startAbove++;
        }
        // search with a radius of m_countsPerFrame for the nearest numbered file
        for( int i = 0; i < m_countsPerFrame; ++i ) {
            // First look below
            if( files::file_exists( get_mesh_filename( startBelow ) ) )
                return startBelow;
            startBelow--;
            // Then look above
            if( files::file_exists( get_mesh_filename( startAbove ) ) )
                return startAbove;
            startAbove++;
        }
        // If we found nothing, just round to the nearest integer
        return (int)frantic::math::round( cacheCountNumber );
    }

    // This finds the nearest files surrounding the requested count number
    void get_nearest_cache_count_interval( float cacheCountNumber, int& outBeforeCountNumber,
                                           int& outAfterCountNumber ) {
        // Only check within stepsize for speed sake
        int endTime = ( std::min )( (int)ceil( cacheCountNumber ) + m_countsPerFrame, m_countIntervalEnd );
        int startTime = ( std::max )( (int)floor( cacheCountNumber ) - m_countsPerFrame, m_countIntervalStart );

        outBeforeCountNumber = (int)floor( cacheCountNumber );
        while( outBeforeCountNumber >= startTime &&
               !frantic::files::file_exists( get_mesh_filename( outBeforeCountNumber ) ) ) {
            outBeforeCountNumber--;
        }

        outAfterCountNumber = (int)ceil( cacheCountNumber );
        while( outAfterCountNumber <= endTime &&
               !frantic::files::file_exists( get_mesh_filename( outAfterCountNumber ) ) ) {
            outAfterCountNumber++;
        }
    }

    // This returns the filename of the specified count.  The file count number is already within the filename numbering
    // system
    frantic::tstring get_mesh_filename( int cacheCountNumber ) const {
        frantic::tstring filename = m_filePrefix;
        if( m_animated )
            filename += strings::zero_pad( cacheCountNumber, m_framePadding );
        filename += m_filePostfix;
        return filename;
    }

    void load_mesh_internal( int cacheCountNumber, trimesh3& mesh ) {
        frantic::tstring filename = get_mesh_filename( cacheCountNumber );

        if( !frantic::files::file_exists( filename ) )
            throw std::runtime_error( "animated_trimesh3_cache.load_mesh: Requested mesh file \"" +
                                      frantic::strings::to_string( filename ) + "\" does not exist" );

        load_mesh_file( filename, mesh );
    }

    void write_mesh_internal( int cacheCountNumber, trimesh3& mesh ) {
        frantic::tstring filename = get_mesh_filename( cacheCountNumber );

        write_mesh_file( filename, mesh );
    }
};

} // namespace geometry
} // namespace frantic
