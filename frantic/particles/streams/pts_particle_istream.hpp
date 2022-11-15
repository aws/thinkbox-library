// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <boost/algorithm/string/join.hpp>

#include <frantic/particles/streams/particle_istream.hpp>

#include <frantic/files/csv_files.hpp>
#include <frantic/files/files.hpp>

#include <frantic/channels/channel_map_adaptor.hpp>

namespace frantic {
namespace particles {
namespace streams {

/**
 * The PTS file format from Leica Cyclone, while simple, is not clearly documented anywhere. Perhaps the closest
 * to reasonable documentation is here (though it claims intensity is required, and explains how to use Excel to
 * add an intensity channel if it is missing from the scanner output):
 * http://helpdesk.microsurvey.com/index.php?/Knowledgebase/Article/View/513/0/cannot-access-data-error-when-opening-a-pts-format-pointcloud
 *
 * The general pattern is as follows:
 *
 * <integer particle count>
 * X Y Z I
 * ... (as many as the count specified
 * <integer particle count>
 * X Y Z I
 * ... (as many as the count specified
 *
 * The individual line formats can be as follows, with distinction between RGB and N via whether the numbers are
 * float or integer. TODO: Our implementation does not support the case with normals presently.
 *
 * X Y Z
 * X Y Z I
 * X Y Z R G B
 * X Y Z I R G B
 * X Y Z I NX NY NZ
 * X Y Z I R G B NZ NY NZ
 */
class pts_particle_istream : public particle_istream {
    frantic::tstring m_name;
    boost::int64_t m_lineNumber;
    boost::int64_t m_currentParticleIndex;
    boost::int64_t m_badLineCount, m_maxBadLineLogMessages;

    // This channel map is deduced from the first few lines of the PTS file
    frantic::channels::channel_map m_nativeChannelMap;
    // This channel map is what was requested, and is used for all processing
    frantic::channels::channel_map m_particleChannelMap;
    std::vector<char> m_defaultParticleBuffer;

    // Accessors and channel info into `m_particleChannelMap` for processing one particle
    const frantic::channels::channel* m_positionChannel;
    const frantic::channels::channel* m_intensityChannel;
    frantic::channels::channel_cvt_accessor<float> m_intensityAcc;
    const frantic::channels::channel* m_colorChannel;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_colorAcc;
    const frantic::channels::channel* m_normalChannel;

    files::file_ptr m_fin;
    bool m_finReopened;
    boost::shared_ptr<frantic::files::line_reader_interface> m_lineReader;
    boost::int64_t m_fileSize;

    // Reuse the same string and vector memory for each particle access
    std::string m_line;
    std::vector<std::string> m_fields;

    static std::string get_file_error_message( const frantic::tstring& filename, const std::string& message ) {
        return "pts_particle_istream: In file '" + frantic::strings::to_string( filename ) + "': " + message;
    }

    std::string get_line_error_message( const frantic::tstring& filename, const std::string& message ) {
        return "pts_particle_istream: In file '" + frantic::strings::to_string( filename ) + "', line " +
               boost::lexical_cast<std::string>( m_lineNumber ) + ": " + message;
    }

    /** Reads the first few lines from the file to figure out the native channel map, then resets the file pointer. This
     * function also validates that it looks sort of like a PTS file */
    void deduce_native_channel_map( const frantic::tstring& file, files::file_ptr& fin,
                                    channels::data_type_t positionTypeHint );

    void initialize_stream( const frantic::tstring& file, channels::data_type_t positionTypeHint );

    /** Prints a warning log message for the currently processed line */
    void log_bad_line( const std::string& msg );
    /** Prints a warning log message about all the warnings, if there were any */
    void log_bad_line_summary();

  public:
    pts_particle_istream( const frantic::tstring& filename, channels::data_type_t positionTypeHint )
        : m_name( filename )
        , m_lineNumber( 0 )
        , m_currentParticleIndex( -1 )
        , m_badLineCount( 0 )
        , m_maxBadLineLogMessages( 20 )
        , m_finReopened( false ) {
        initialize_stream( m_name, positionTypeHint );
        set_channel_map( m_nativeChannelMap );
    }

    pts_particle_istream( const frantic::tstring& filename, const frantic::channels::channel_map& channelMap )
        : m_name( filename )
        , m_lineNumber( 0 )
        , m_currentParticleIndex( -1 )
        , m_badLineCount( 0 )
        , m_maxBadLineLogMessages( 20 )
        , m_finReopened( false ) {
        initialize_stream( m_name, channelMap[_T("Position")].data_type() );
        set_channel_map( channelMap );
    }

    virtual ~pts_particle_istream() { close(); }

    void close() { m_fin.close(); }

    std::size_t particle_size() const { return m_particleChannelMap.structure_size(); }

    frantic::tstring name() const { return m_name; }

    boost::int64_t particle_count() const {
        // .pts files have a count for a block, but there may be multiple blocks in the file, so we can't know the whole
        // count. Due to high occurrence of inconsistent/corrupt .pts files, we're also doing more liberal parsing, thus
        // don't trust the count anyway.
        return -1;
    }

    boost::int64_t particle_index() const { return m_currentParticleIndex; }

    boost::int64_t particle_count_left() const { return -1; }

    boost::int64_t particle_progress_count() const { return m_fileSize; }

    boost::int64_t particle_progress_index() const {
        // Use the # bytes read in the file for progress
        if( m_lineReader ) {
            return m_lineReader->get_file_progress();
        } else {
            return m_finReopened ? m_fileSize : 0;
        }
    }

    const channels::channel_map& get_channel_map() const { return m_particleChannelMap; }

    // Access to the channel_map in the file
    const channels::channel_map& get_native_channel_map() const { return m_nativeChannelMap; }

    void set_default_particle( char* buffer ) {
        m_particleChannelMap.copy_structure( &m_defaultParticleBuffer[0], buffer );
    }

    //////////////////////
    // Manipulation and retrieval functions
    //////////////////////

    // This allows you to change the particle layout that's being loaded on the fly, in case it couldn't
    // be set correctly at creation time.
    void set_channel_map( const channels::channel_map& particleChannelMap );

    bool get_particle( char* rawParticleBuffer );

    bool get_particles( char* particleBuffer, std::size_t& numParticles );
};

} // namespace streams
} // namespace particles
} // namespace frantic
