// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/pts_particle_istream.hpp>

using frantic::graphics::vector3f;
using frantic::graphics::vector3fd;
using std::string;
using std::vector;
using namespace frantic;
using namespace frantic::particles::streams;

namespace {
struct pts_field_locations {
    // Each int is -1 or its starting index within the fields.
    // If position is -1, the line was invalid.
    int position, intensity, color, normal;
    pts_field_locations()
        : position( -1 )
        , intensity( -1 )
        , color( -1 )
        , normal( -1 ) {}
};
} // anonymous namespace

/** Detects where the PTS fields are using some simple heuristics. Sets `out.X` to the found index if X is detected,
 * otherwise leaves it untouched */
static void detect_pts_fields( const std::vector<std::string>& fields, pts_field_locations& out ) {
    switch( fields.size() ) {
    case 3:
        // Just position
        out.position = 0;
        break;
    case 4:
        // Position and intensity
        out.position = 0;
        out.intensity = 3;
        break;
    case 6:
        // Position and (color or normal)
        out.position = 0;
        if( fields[3].find( '.' ) == string::npos ) {
            out.color = 3;
        } else {
            out.normal = 3;
        }
        break;
    case 7:
        // Position, intensity, and (color or normal)
        out.position = 0;
        out.intensity = 3;
        if( fields[4].find( '.' ) == string::npos ) {
            out.color = 4;
        } else {
            out.normal = 4;
        }
        break;
    case 10:
        // Position, intensity, color, and normal
        out.position = 0;
        out.intensity = 3;
        out.color = 4;
        out.normal = 7;
        break;
    }
}

/** Splits a line by spaces. Tries to avoid reallocation in `out`. */
static void split_line( std::string& line, std::vector<std::string>& out ) {
    if( line.empty() ) {
        out.clear();
        return;
    }
    size_t outPos = 0;
    const char *begin = &line[0], *end = begin + line.size();
    // Skip spaces
    while( begin < end && isspace( *begin ) )
        ++begin;
    while( begin < end ) {
        // Find the end of the token
        const char* pos = begin;
        while( pos < end && !isspace( *pos ) )
            ++pos;
        // Copy the string, avoiding reallocation in `out`
        if( outPos < out.size() ) {
            out[outPos].resize( pos - begin );
            memcpy( &out[outPos][0], begin, pos - begin );
        } else {
            out.push_back( string( begin, pos ) );
        }
        ++outPos;
        begin = pos;
        // Skip spaces
        while( begin < end && isspace( *begin ) )
            ++begin;
    }
    if( outPos < out.size() )
        out.resize( outPos );
}

/** Reads a single line from the file, and splits it apart by whitespace. Returns true if a line was read */
static bool read_line( files::line_reader_interface& lr, std::vector<std::string>& out ) {
    string line;
    if( !lr.getline( line ) )
        return false;
    split_line( line, out );
    return true;
}

void pts_particle_istream::log_bad_line( const std::string& msg ) {
    ++m_badLineCount;
    if( m_badLineCount < m_maxBadLineLogMessages ) {
        FF_LOG( warning ) << "Bad PTS line " << m_lineNumber << ": " << msg.c_str() << std::endl;
    }
}

void pts_particle_istream::log_bad_line_summary() {
    if( m_badLineCount > 0 ) {
        FF_LOG( warning ) << "Ignored " << m_badLineCount << " bad lines in PTS file \"" << m_name << "\"" << std::endl;
    }
}

void pts_particle_istream::deduce_native_channel_map( const frantic::tstring& file, files::file_ptr& fin,
                                                      channels::data_type_t positionTypeHint ) {
    boost::shared_ptr<frantic::files::line_reader_interface> lr = frantic::files::create_line_reader( fin.get(), true );
    // First, make sure the file seems valid.
    // We expect to find an integer (and a newline) in the first line
    // allow whitespace following the integer, but limit the number of characters we examine
    if( !lr->getline( m_line ) || m_line.empty() ) {
        throw frantic::invalid_particle_file_exception(
            get_file_error_message( file, "Expected a count in the first line" ) );
    }
    if( isspace( m_line[0] ) ) {
        throw frantic::invalid_particle_file_exception(
            get_file_error_message( file, "Unexpected white space before the count in first line." ) );
    }
    size_t i = 0;
    // Expect a point count
    for( ; i != m_line.size(); ++i ) {
        if( i != 0 && isspace( m_line[i] ) )
            break;
        if( !isdigit( m_line[i] ) ) {
            throw frantic::invalid_particle_file_exception( get_file_error_message(
                file, "Unexpected characters in first line.  Expected only a count in first line." ) );
        }
    }
    // We allow space, but nothing else after the point count
    for( ; i != m_line.size(); ++i ) {
        if( !isspace( m_line[i] ) ) {
            throw frantic::invalid_particle_file_exception(
                get_file_error_message( file, "The first line should only contain a count" ) );
        }
    }

    // With the first line out of the way, and seemingly a valid point count, read a few lines
    // to deduce the native channel map
    pts_field_locations ptsFields;
    vector<string> fields;
    for( i = 0; i < 20; ++i ) {
        if( !read_line( *lr.get(), fields ) )
            break;
        // All the fields should be parsable as double
        try {
            for( size_t j = 0; j != fields.size(); ++j ) {
                boost::lexical_cast<double>( fields[j] );
            }
        } catch( const boost::bad_lexical_cast& ) {
            throw frantic::invalid_particle_file_exception(
                get_file_error_message( file, "Non-numerical data found near beginning of file" ) );
        }
        detect_pts_fields( fields, ptsFields );
    }

    // Define the native channel map based on what we saw
    m_nativeChannelMap.define_channel( _T("Position"), 3, positionTypeHint );
    if( ptsFields.intensity >= 0 )
        m_nativeChannelMap.define_channel<float>( _T("Intensity") );
    if( ptsFields.color >= 0 )
        m_nativeChannelMap.define_channel( _T("Color"), 3, channels::data_type_float16 );
    if( ptsFields.normal >= 0 )
        m_nativeChannelMap.define_channel( _T("Normal"), 3, channels::data_type_float32 );
    m_nativeChannelMap.end_channel_definition();

    // Reset the file
    lr.reset();
    std::clearerr( fin );
    std::fseek( fin, 0, SEEK_SET );
}

void pts_particle_istream::initialize_stream( const frantic::tstring& file, channels::data_type_t positionTypeHint ) {
    if( positionTypeHint == channels::data_type_invalid ) {
        positionTypeHint = channels::data_type_float32;
    }

    m_finReopened = false;

    files::file_ptr fin( frantic::files::tfopen( file.c_str(), _T("rb") ) );
    if( !fin ) {
        throw std::runtime_error( get_file_error_message( file, "Unable to open file for reading." ) );
    }

    int error = fseek64( fin, 0, std::ios::end );
    if( error == 0 ) { // 0 means successful
        m_fileSize = ftell64( fin );
    } else {
        m_fileSize = -1;
    }
    error = fseek64( fin, 0, std::ios::beg );
    if( error ) {
        throw std::runtime_error( "pts_particle_istream: fseek failed for input file \"" +
                                  frantic::strings::to_string( file ) + "\"" );
    }

    deduce_native_channel_map( file, fin, positionTypeHint );
}

void pts_particle_istream::set_channel_map( const channels::channel_map& particleChannelMap ) {
    std::vector<char> newDefaultParticle( particleChannelMap.structure_size() );
    if( m_defaultParticleBuffer.size() > 0 ) {
        frantic::channels::channel_map_adaptor defaultAdaptor( particleChannelMap, m_particleChannelMap );
        defaultAdaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticleBuffer[0] );
    } else {
        memset( &newDefaultParticle[0], 0, particleChannelMap.structure_size() );
    }
    m_defaultParticleBuffer.swap( newDefaultParticle );

    // Set the map and the adaptor
    m_particleChannelMap = particleChannelMap;

    // Reset all the channel info for parsing
    m_positionChannel = &m_particleChannelMap[_T("Position")];
    if( m_particleChannelMap.has_channel( _T("Intensity") ) ) {
        m_intensityChannel = &m_particleChannelMap[_T("Intensity")];
        m_intensityAcc = m_particleChannelMap.get_cvt_accessor<float>( _T("Intensity") );
    } else {
        m_intensityChannel = NULL;
        m_intensityAcc.reset();
    }
    if( m_particleChannelMap.has_channel( _T("Color") ) ) {
        m_colorChannel = &m_particleChannelMap[_T("Color")];
        m_colorAcc = m_particleChannelMap.get_cvt_accessor<vector3f>( _T("Color") );
    } else {
        m_colorChannel = NULL;
        m_colorAcc.reset();
    }
    if( m_particleChannelMap.has_channel( _T("Normal") ) ) {
        m_normalChannel = &m_particleChannelMap[_T("Normal")];
    } else {
        m_normalChannel = NULL;
    }
}

/** Parses a channel from a series of adjacent string fields */
static void parse_channel_value_from_fields( const vector<string>& fields, int fieldIndex,
                                             const frantic::channels::channel& ch, char* outBuffer ) {
    size_t elementSize = channels::sizeof_channel_data_type( ch.data_type() );
    frantic::channels::data_type_t dataType = ch.data_type();
    const string* field = &fields[fieldIndex];
    outBuffer += ch.offset();
    for( size_t i = 0, iend = ch.arity(); i != iend; ++i, outBuffer += elementSize, ++field ) {
        channels::parse_channel_value_from_string( dataType, *field, outBuffer );
    }
}

bool pts_particle_istream::get_particle( char* rawParticleBuffer ) {
    if( !m_finReopened ) {
        m_finReopened = true;
        m_fin.reset( frantic::files::tfopen( m_name.c_str(), _T("rb") ) );
        if( !m_fin ) {
            throw std::runtime_error( get_file_error_message( name(), "Could not reopen the file." ) );
        }

        m_lineReader = frantic::files::create_line_reader( m_fin.get(), true );
    } else if( !m_fin ) {
        return false;
    }

    // Start with the default particle, to avoid complicated logic about which fields are available
    m_particleChannelMap.copy_structure( rawParticleBuffer, m_defaultParticleBuffer );

    while( true ) {
        if( !m_lineReader->getline( m_line ) ) {
            // Reached the end of file, so close things up
            m_lineReader.reset();
            m_fin.close();
            log_bad_line_summary();
            // Free the memory used by the line and fields buffers
            string().swap( m_line );
            vector<string>().swap( m_fields );
            return false;
        }
        ++m_lineNumber;
        // Skip excessively long lines
        if( m_line.size() > 1000 ) {
            log_bad_line( "Excessively long line with " + boost::lexical_cast<string>( m_line.size() ) +
                          " characters" );
            string().swap( m_line );
            continue;
        }

        split_line( m_line, m_fields );

        try {
            // Figure out the field mapping based on the number of fields
            pts_field_locations ptsFields;
            detect_pts_fields( m_fields, ptsFields );
            if( ptsFields.position < 0 ) {
                // Log an error for lines with unexpected field count, except for 1 field which we ignore because it
                // represents a block size (which we ignore due to our relaxed parsing)
                if( m_fields.size() > 1 )
                    log_bad_line( "Line has unexpected number of fields " +
                                  boost::lexical_cast<string>( m_fields.size() ) );
                continue;
            }

            // Set the values in the destination particle
            parse_channel_value_from_fields( m_fields, ptsFields.position, *m_positionChannel, rawParticleBuffer );
            if( ptsFields.intensity >= 0 && m_intensityChannel ) {
                if( m_fields[ptsFields.intensity].find( '.' ) == string::npos ) {
                    // Treat as an integer with -4096 to 4095 range if it has no decimal point
                    int value = 0;
                    channels::parse_channel_value_from_string( channels::data_type_int32, m_fields[ptsFields.intensity],
                                                               reinterpret_cast<char*>( &value ) );
                    m_intensityAcc.set( rawParticleBuffer, ( value + 4096 ) / 8191.f );
                } else {
                    // Treat as a float with 0.0-1.0 range if it has a decimal point
                    parse_channel_value_from_fields( m_fields, ptsFields.intensity, *m_intensityChannel,
                                                     rawParticleBuffer );
                }
            }
            if( ptsFields.color >= 0 && m_colorAcc.is_valid() ) {
                boost::uint8_t color[3] = { 0, 0, 0 };
                channels::parse_channel_value_from_string( channels::data_type_uint8, m_fields[ptsFields.color + 0],
                                                           reinterpret_cast<char*>( &color[0] ) );
                channels::parse_channel_value_from_string( channels::data_type_uint8, m_fields[ptsFields.color + 1],
                                                           reinterpret_cast<char*>( &color[1] ) );
                channels::parse_channel_value_from_string( channels::data_type_uint8, m_fields[ptsFields.color + 2],
                                                           reinterpret_cast<char*>( &color[2] ) );
                m_colorAcc.set( rawParticleBuffer, vector3f( color[0] / 255.f, color[1] / 255.f, color[2] / 255.f ) );
            }
            if( ptsFields.normal >= 0 && m_normalChannel ) {
                parse_channel_value_from_fields( m_fields, ptsFields.normal, *m_normalChannel, rawParticleBuffer );
            }
        } catch( const std::exception& e ) {
            log_bad_line( string() + "Line is invalid: " + e.what() );
            continue;
        }

        return true;
    }
}

bool pts_particle_istream::get_particles( char* particleBuffer, std::size_t& numParticles ) {
    std::size_t particleSize = m_particleChannelMap.structure_size();

    for( std::size_t i = 0; i < numParticles; ++i ) {
        if( !get_particle( particleBuffer + i * particleSize ) ) {
            numParticles = i;
            return false;
        }
    }

    return true;
}