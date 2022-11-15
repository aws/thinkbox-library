// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/assign/list_of.hpp>
#include <functional>
#include <map>

#include <frantic/files/csv_files.hpp>
#include <frantic/misc/string_functions.hpp>

#include <frantic/files/files.hpp>

#include <ctype.h>

using namespace std;

namespace frantic {
namespace files {

/**
 * This function splits a single line of a .csv file into an array of entries.
 * See http://www.rfc-editor.org/rfc/rfc4180.txt for details.
 *
 * You can use feof(in) to check if the end of file was reached by this call.
 *
 * @param  in  The C FILE* pointer.
 * @param  streamName  The stream or file name, for error messages.
 * @param  outColumnData  A vector of strings for the values.
 */
void read_csv_file_line( FILE* in, const frantic::tstring& streamName, std::vector<std::string>& outColumnData ) {
    outColumnData.clear();
    csv_reader csvReader( in, streamName, false );
    csvReader.read_line( outColumnData );
}

//
// csv_reader
//

bool csv_reader::is_delimiter( char c ) { return m_delimiter == c; }

void csv_reader::init( FILE* file, const frantic::tstring& streamName, bool useFileBuffer ) {
    m_streamName = streamName;
    m_delimiter = ',';
    m_useFileBuffer = useFileBuffer;
    m_lineReader = create_line_reader( file, useFileBuffer );
    m_file = file;

    m_fileSize = -1;
    if( file ) {
        boost::int64_t initialPosition = ftell64( file );
        if( initialPosition >= 0 ) {
            int error = fseek64( file, 0, std::ios::end );
            if( error == 0 ) { // 0 means successful
                m_fileSize = ftell64( file );
            }
            error = fseek64( file, initialPosition, std::ios::beg );
            if( error ) {
                throw std::runtime_error( "csv_reader::init fseek failed for input file \"" +
                                          frantic::strings::to_string( streamName ) + "\"" );
            }
        }
    }

    m_emitEmptyColumns = true;
    m_beginRowLineNumber = -1;
    m_lineNumber = -1;
}

bool csv_reader::getline( std::string& outLine ) {
    bool success = false;

    if( m_lineReader ) {
        success = m_lineReader->getline( outLine );
    } else {
        throw std::runtime_error( "csv_reader.getline: File reader is NULL" );
    }

    if( success ) {
        ++m_lineNumber;
    }

    return success;
}

void csv_reader::do_emit( const std::string& s, std::vector<std::string>& outColumnData, std::size_t& outColumn ) {
    if( outColumn >= outColumnData.size() ) {
        outColumnData.resize( outColumn + 1 );
    }
    outColumnData[outColumn++].assign( s );
}

void csv_reader::maybe_emit( const std::string& s, std::vector<std::string>& outColumnData, std::size_t& outColumn ) {
    if( s.size() > 0 || m_emitEmptyColumns ) {
        do_emit( s, outColumnData, outColumn );
    }
}

void csv_reader::maybe_emit( const std::string& s, std::string::size_type startIndex, std::string::size_type endIndex,
                             std::vector<std::string>& outColumnData, std::size_t& outColumn ) {
    if( endIndex > startIndex || m_emitEmptyColumns ) {
        if( outColumn >= outColumnData.size() ) {
            outColumnData.resize( outColumn + 1 );
        }
        outColumnData[outColumn++].assign( s, startIndex, endIndex - startIndex );
    }
}

csv_reader::csv_reader( FILE* file, const frantic::tstring& streamName, bool useBuffer ) {
    init( file, streamName, useBuffer );
}

csv_reader::~csv_reader() {}

void csv_reader::set_delimiter( char c ) { m_delimiter = c; }

void csv_reader::set_emit_empty_columns( bool emitEmptyColumns ) { m_emitEmptyColumns = emitEmptyColumns; }

bool csv_reader::read_line( std::vector<std::string>& outColumnData ) {
    // Read in a line, until the first EOL character
    const bool result = getline( m_line );
    std::size_t outColumn = 0;

    if( result ) {
        m_beginRowLineNumber = m_lineNumber;
    }

    // Parse the line into a sequence of tokens.
    std::string::size_type index = 0, nextDelimIndex, nextQuoteIndex;
    nextDelimIndex = m_line.find( m_delimiter, index );
    nextQuoteIndex = m_line.find( '"', index );
    while( index < m_line.size() ) {
        // If this entry is enclosed in quotes, deal with that specially
        if( nextQuoteIndex != std::string::npos && nextQuoteIndex < nextDelimIndex ) {
            boost::int64_t openQuoteLineNumber = m_lineNumber;

            // Find the quote character, making sure that only whitespace occurs before
            char lineChar = m_line[index];
            while( lineChar != '"' ) {
                if( !isspace( lineChar ) )
                    throw std::runtime_error( "csv_reader.read_line: The input file \"" +
                                              frantic::strings::to_string( m_streamName ) + "\" line " +
                                              boost::lexical_cast<std::string>( openQuoteLineNumber + 1 ) +
                                              " had a quote inside a non-quoted column data entry." );
                lineChar = m_line[++index];
            }

            // Build up the token.  Special cases are "", which gets collapsed to a single quote, and the end of line,
            // which requires us to get another line from the file, resulting in a \n character in the token.
            ++index;
            nextQuoteIndex = m_line.find( '"', index );
            bool done = false;
            m_token.clear();
            while( !done ) {
                if( nextQuoteIndex != std::string::npos ) {
                    // Two double quotes in a row ("") get collapsed into a double quote character
                    if( nextQuoteIndex + 1 < m_line.size() && m_line[nextQuoteIndex + 1] == '"' ) {
                        m_token.append( m_line, index, nextQuoteIndex - index + 1 );
                        index = nextQuoteIndex + 2;
                        nextQuoteIndex = m_line.find( '"', index );
                    }
                    // Otherwise, if there is just one quote, it must be the end of the quoted string
                    else {
                        m_token.append( m_line, index, nextQuoteIndex - index );
                        index = nextQuoteIndex + 1;
                        done = true;
                    }
                } else {
                    m_token.append( m_line, index, m_line.size() - index );
                    m_token += "\n";
                    if( !getline( m_line ) )
                        throw std::runtime_error( "csv_reader.read_line: The input file \"" +
                                                  frantic::strings::to_string( m_streamName ) + "\" line " +
                                                  boost::lexical_cast<std::string>( openQuoteLineNumber + 1 ) +
                                                  " had a non-terminated quoted column data entry." );
                    index = 0;
                    nextQuoteIndex = m_line.find( '"', index );
                }
            }

            // Now ensure that there is just whitespace up to the next ',' character.
            if( index < m_line.size() ) {
                lineChar = m_line[index];
                while( lineChar != m_delimiter && index < m_line.size() ) {
                    if( !isspace( lineChar ) )
                        throw std::runtime_error( "csv_reader.read_line: The input file \"" +
                                                  frantic::strings::to_string( m_streamName ) + "\" line " +
                                                  boost::lexical_cast<std::string>( m_lineNumber + 1 ) +
                                                  " had a quote inside a non-quoted column data entry." );
                    if( ++index < m_line.size() )
                        lineChar = m_line[index];
                }
                // Skip past the ',' character to the start of the next token.
                if( index < m_line.size() )
                    ++index;
            }

            do_emit( m_token, outColumnData, outColumn );
            nextQuoteIndex = m_line.find( '"', index );
        }
        // If there are more delimiters, grab the the token up to the delimiter
        else if( nextDelimIndex != std::string::npos ) {
            maybe_emit( m_line, index, nextDelimIndex, outColumnData, outColumn );
            index = nextDelimIndex + 1;
        }
        // Otherwise, grab this token up till the end of the line, and parsing this line is done
        else {
            maybe_emit( m_line, index, m_line.size(), outColumnData, outColumn );
            index = m_line.size();
        }
        nextDelimIndex = m_line.find( m_delimiter, index );
    }
    // Special case the situation when there's an empty value at the end of the line.
    if( m_line.size() > 0 && m_line[m_line.size() - 1] == m_delimiter ) {
        maybe_emit( "", outColumnData, outColumn );
    }

    outColumnData.resize( outColumn );

    return result;
}

boost::int64_t csv_reader::get_line_number() const { return m_beginRowLineNumber; }

boost::int64_t csv_reader::file_progress_count() const { return m_fileSize; }

boost::int64_t csv_reader::file_progress_index() const { return m_lineReader->get_file_progress(); }

bool guess_csv_delimiter( FILE* in, const frantic::tstring& streamName, char& outDelimiter ) {
    std::string line;
    frantic::strings::getline( in, line );
    bool inQuote = false;

    if( !in ) {
        throw std::runtime_error( "guess_csv_delimiter Error: input file is NULL" );
    }

    // delimiters, starting from top priority
    std::vector<char> delimiters = boost::assign::list_of( ',' )( ';' )( ':' )( '\t' )( ' ' );

    typedef std::map<char, std::size_t> delimiter_map_t;
    delimiter_map_t delimiterMap;
    for( std::vector<char>::iterator i = delimiters.begin(); i != delimiters.end(); ++i ) {
        delimiterMap[*i] = 0;
    }

    boost::int64_t lineNumber = 1;
    boost::int64_t openQuoteLineNumber = 1;

    int index = 0;

    // go through each character in the record
    while( index < static_cast<int>( line.size() ) ) {
        const char c = line[index];
        if( inQuote ) {
            if( c == '\"' ) {
                // inside a quoted block, two quotes in a row are interpreted
                // as a single escaped quote
                if( index + 1 < static_cast<int>( line.size() ) && line[index + 1] == '\"' ) {
                    ++index;
                } else {
                    inQuote = false;
                }
            } else {
                // ignore delimiters inside quotes
            }
            if( index + 1 == static_cast<int>( line.size() ) ) {
                frantic::strings::getline( in, line );
                index = -1;
                ++lineNumber;
            }
        } else {
            if( c == '\"' ) {
                inQuote = true;
                openQuoteLineNumber = lineNumber;
            } else {
                // count the incidence of the delimiters
                delimiter_map_t::iterator i = delimiterMap.find( c );
                if( i != delimiterMap.end() ) {
                    ++( i->second );
                }
            }
        }
        ++index;
    }
    if( inQuote ) {
        throw std::runtime_error( "guess_csv_delimiter: The input file \"" + frantic::strings::to_string( streamName ) +
                                  "\" line " + boost::lexical_cast<std::string>( openQuoteLineNumber ) +
                                  " had a non-terminated quoted column data entry." );
    }
    // pick first delimiter with count > 0
    for( std::vector<char>::iterator i = delimiters.begin(); i != delimiters.end(); ++i ) {
        if( delimiterMap[*i] > 0 ) {
            outDelimiter = *i;
            return true;
        }
    }
    return false;
}

//
// line_reader_interface
//

line_reader_interface::~line_reader_interface() {}

//
// buffered_line_reader
//

/**
 *  Read a FILE into an internal buffer, and return getline and
 * fgetc requests from this buffer.
 *
 *  This is used internally by csv_reader.
 */
class buffered_line_reader : public line_reader_interface {
    std::size_t m_index;
    std::vector<char> m_buffer;
    FILE* m_file;
    std::size_t m_bufferSize;
    char m_prevChar;
    boost::int64_t m_fileProgress;

    void prepare_buffer();

  public:
    buffered_line_reader();
    buffered_line_reader( FILE* file );

    boost::int64_t get_file_progress();
    /**
     *  Get a line of text from the input file.
     *
     * @param line if the function returns true, then this is the next
     *		line of text in the input file.
     * @return true is a line could be read, or false if an error occurred
     *		or the end of file was reached.
     */
    bool getline( std::string& line );
};

void buffered_line_reader::prepare_buffer() {
    if( m_bufferSize <= 0 ) {
        throw std::runtime_error( "bufferSize must be > 0" );
    }
    m_buffer.resize( m_bufferSize );
    size_t result = fread( &m_buffer[0], sizeof( char ), m_bufferSize, m_file );
    m_buffer.resize( result );
    m_index = 0;
    m_fileProgress += result;
}

buffered_line_reader::buffered_line_reader()
    : m_index( std::numeric_limits<std::size_t>::max() )
    , m_file( 0 )
    , m_bufferSize( 0 ) {}

buffered_line_reader::buffered_line_reader( FILE* file )
    : m_index( std::numeric_limits<std::size_t>::max() )
    , m_file( file )
    , m_bufferSize( 32768 )
    , m_prevChar( 0 )
    , m_fileProgress( 0 ) {}

boost::int64_t buffered_line_reader::get_file_progress() { return m_fileProgress; };

bool buffered_line_reader::getline( std::string& line ) {
    line.clear();
    if( m_index >= m_buffer.size() ) {
        prepare_buffer();
    }
    if( m_index < m_buffer.size() ) {
        for( ;; ) {
            { // scope for current m_buffer and its associated variables
                const char* buffer = &m_buffer[0];

                // skip the '\n' in "\r\n"
                if( m_prevChar == '\r' && buffer[m_index] == '\n' ) {
                    m_prevChar = '\n';
                    ++m_index;
                }

                const std::size_t startIndex = m_index;
                const std::size_t bufferSize = m_buffer.size();

                for( ; m_index < bufferSize; ++m_index ) {
                    bool doneLine = false;

                    const char c = buffer[m_index];
                    if( c == '\r' || c == '\n' ) {
                        line.append( buffer + startIndex, m_index - startIndex );
                        ++m_index;
                        doneLine = true;
                    }

                    m_prevChar = c;

                    if( doneLine ) {
                        return true;
                    }
                }
                line.append( buffer + startIndex, bufferSize - startIndex );
            }
            prepare_buffer();
            if( m_index == m_buffer.size() ) {
                if( line.empty() ) {
                    // skip the final '\n' when the file ends with "\r\n"
                    return false;
                } else {
                    return true;
                }
            }
        }
    } else {
        return false;
    }
}

//
// unbuffered_line_reader
//

class unbuffered_line_reader : public line_reader_interface {
    FILE* m_file;
    int m_prevChar;
    boost::int64_t m_fileProgress;

  public:
    unbuffered_line_reader( FILE* file );

    boost::int64_t get_file_progress();
    /**
     *  Get a line of text from the input file.
     *
     * @param line if the function returns true, then this is the next
     *		line of text in the input file.
     * @return true is a line could be read, or false if an error occurred
     *		or the end of file was reached.
     */
    bool getline( std::string& line );
};

unbuffered_line_reader::unbuffered_line_reader( FILE* file )
    : m_file( file )
    , m_prevChar( 0 )
    , m_fileProgress( 0 ) {}

boost::int64_t unbuffered_line_reader::get_file_progress() { return m_fileProgress; };

bool unbuffered_line_reader::getline( std::string& outString ) {
    outString.clear();

    int c = std::fgetc( m_file );
    if( c != EOF ) {
        ++m_fileProgress;
    }

    // skip the '\n' in "\r\n"
    if( m_prevChar == '\r' && c == '\n' ) {
        m_prevChar = c;
        c = std::fgetc( m_file );
        if( c != EOF ) {
            ++m_fileProgress;
        }
    }

    if( c == '\n' || c == '\r' ) {
        m_prevChar = c;
        return true;
    }

    if( c == EOF ) {
        m_prevChar = c;
        return false;
    }

    for( ;; ) {
        outString += static_cast<std::string::value_type>( c );
        c = std::fgetc( m_file );
        m_prevChar = c;
        if( c != EOF ) {
            ++m_fileProgress;
        }
        if( c == '\n' || c == '\r' || c == EOF )
            return true;
    }
}

boost::shared_ptr<frantic::files::line_reader_interface> create_line_reader( FILE* file, bool useFileBuffer ) {
    if( useFileBuffer ) {
        return boost::make_shared<buffered_line_reader>( file );
    } else {
        return boost::make_shared<unbuffered_line_reader>( file );
    }
}

} // namespace files
} // namespace frantic
