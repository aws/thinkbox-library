// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/csv_particle_istream.hpp>

using namespace std;
using namespace frantic;
using frantic::channels::channel_map;
using frantic::channels::data_type_t;

namespace frantic {
namespace particles {
namespace streams {

namespace {

void guess_native_channel_map_from_numeric_headers( const std::vector<std::string>& columnHeaders, bool noError,
                                                    data_type_t positionTypeHint,
                                                    frantic::channels::channel_column_map& outGuess ) {
    frantic::channels::channel_map channelMap;

    if( !noError && columnHeaders.size() < 3 ) {
        throw frantic::invalid_particle_file_exception()
            << "csv_particle_istream : CSV particle file must have at least 3 columns for position information.";
    } else {
        // Always specifies a position channel with xyz components
        channelMap.define_channel( _T("Position"), 3,
                                   ( positionTypeHint == channels::data_type_invalid ) ? channels::data_type_float32
                                                                                       : positionTypeHint );

        if( columnHeaders.size() == 6 ) {
            channelMap.define_channel( _T("Color"), 3, frantic::channels::data_type_float16 );
        } else {
            // Otherwise, remaining columns with Data1, Data2, etc...
            for( size_t i = 3; i < columnHeaders.size(); ++i ) {
                channelMap.define_channel( _T( "Data" ) + boost::lexical_cast<frantic::tstring>( i - 2 ), 1,
                                           frantic::channels::data_type_float32 );
            }
        }
    }

    channelMap.end_channel_definition();

    outGuess.reset( channelMap, columnHeaders.size() );

    for( size_t i = 0; i < 3; ++i ) {
        outGuess.set_column_mapping( i, _T("Position"), i );
    }

    if( columnHeaders.size() == 6 ) {
        for( size_t i = 0; i < 3; ++i ) {
            outGuess.set_column_mapping( 3 + i, _T("Color"), i );
        }
    } else {
        for( size_t i = 3; i < columnHeaders.size(); ++i ) {
            outGuess.set_column_mapping( i, _T( "Data" ) + boost::lexical_cast<frantic::tstring>( i - 2 ), 0 );
        }
    }
}

/**
 * Holds the data type and optionally the index of a numbered column.
 */
struct column_mapping {
    column_mapping() {}
    column_mapping( size_t _column, frantic::channels::data_type_t _type )
        : column( _column )
        , type( _type )
        , hasIndex( false ) {}
    column_mapping( size_t _column, frantic::channels::data_type_t _type, size_t _index )
        : column( _column )
        , type( _type )
        , hasIndex( true )
        , index( _index ) {}
    size_t column;
    frantic::channels::data_type_t type;
    bool hasIndex;
    size_t index;
    bool operator<( const column_mapping& other ) const { return index < other.index; }
};

/**
 * Attempts to construct a channel map and associated column mappings to it. This allows interleaved channels, and
 * 'gaps' in the data
 *
 * The guessing method assumes columns are in the form  "type name[index]"
 * - "type" is optional on all columns, however any columns with the same "name" should match the type of all other
 *   columns with that name, or be blank
 * -- If all columns of the same name are blank, the type will default to float32
 * - "[index]" is optional on all columns
 * -- If indices are omitted for all columns, their order and arity is implicit by the order they appear
 * -- If indices are given for all columns, they must be contiguous, and either 0 or 1-based indices
 * -- It is an error if only some indices are specified
 *
 * These rules are effectively arbitrary, if any changes seem reasonable to you, go for it. Note that these rules are
 * loosened in 'noError' mode to prevent the method from ever failing, at the cost of accepting weird inputs with weird
 * (though still consistent) results.
 *
 * \param columnHeaders  Set of strings at the head of the csv file
 * \param file  Name of the file which the data is loaded from
 * \param noError Loosens the rules to prevent the method ever throwing any errors. Used to present the user with a best
 *                guess and lets them refine it later
 * \param positionTypeHint  The dtype (data_type_float32 or data_type_float64) that the Position channel should have.
 * Use data_type_invalid to let the function choose. \param outGuess The resulting channel column mapping for this csv
 * file
 *
 */
void guess_native_channel_map_from_text_headers( const std::vector<std::string>& columnHeaders,
                                                 const frantic::tstring& file, bool noError,
                                                 frantic::channels::data_type_t positionTypeHint,
                                                 frantic::channels::channel_column_map& outGuess ) {
    typedef std::map<frantic::tstring, std::vector<column_mapping>> column_map;
    column_map columnMappings;
    vector<frantic::tstring> columnNames;

    std::vector<std::string> headerTokens;
    for( unsigned i = 0; i < columnHeaders.size(); ++i ) {
        headerTokens.clear();
        frantic::strings::split( columnHeaders[i], headerTokens );

        // I think its reasonable to say that an empty column can be treated like a gap
        if( headerTokens.empty() ) {
            continue;
        }

        if( !noError && headerTokens.size() > 2 )
            throw frantic::invalid_particle_file_exception()
                << "csv_particle_istream: The input file \"" + frantic::strings::to_string( file ) +
                       "\" has an invalid column header entry, \"" + columnHeaders[i] + "\".";

        frantic::tstring channelName;
        column_mapping mapping;
        mapping.column = i;
        mapping.type = frantic::channels::data_type_invalid;
        mapping.hasIndex = false;

        if( headerTokens.size() == 1 ) {
            channelName = frantic::strings::to_tstring( headerTokens[0] );
            if( mapping.type == frantic::channels::data_type_invalid ) {
                // "If the float32/float16 definition is missing, the type will default to float32."
                mapping.type = frantic::channels::data_type_float32;
            }
        } else {
            mapping.type = frantic::channels::channel_data_type_from_string( headerTokens[0] );

            if( mapping.type == frantic::channels::data_type_invalid ) {
                if( !noError ) {
                    throw frantic::invalid_particle_file_exception(
                        "csv_particle_istream: The input file \"" + frantic::strings::to_string( file ) +
                        "\" has an invalid column header type, \"" + headerTokens[0] +
                        "\".  Some examples of valid type names are \"float32\" and \"int16\"." );
                }
            }

            channelName = frantic::strings::to_tstring( headerTokens[1] );
        }

        size_t startBracketIndex = channelName.find( '[' );
        frantic::tstring indexNumberStr;
        int indexNumber = -1;
        // The indexing bracket can't be at index 0, there must be a channel name in there too
        if( startBracketIndex == 0 ) {
            if( noError ) {
                continue;
            } else {
                throw frantic::invalid_particle_file_exception(
                    "csv_particle_istream: The input file \"" + frantic::strings::to_string( file ) +
                    "\" has an invalid column header name, \"" + columnHeaders[i] +
                    "\".  There must be a channel name, not just an index." );
            }
        }
        // If there are indexing brackets, extract the index number
        else if( startBracketIndex != frantic::tstring::npos ) {
            size_t endBracketIndex = channelName.find( ']', startBracketIndex );
            if( endBracketIndex != channelName.size() - 1 ) {
                if( noError ) {
                    if( endBracketIndex < startBracketIndex || endBracketIndex == frantic::tstring::npos ) {
                        endBracketIndex = channelName.size();
                    }
                } else {
                    throw frantic::invalid_particle_file_exception(
                        "csv_particle_istream: The input file \"" + frantic::strings::to_string( file ) +
                        "\" has an invalid column header name, \"" + columnHeaders[i] +
                        "\".  The index brackets are incorrect, they should be like \"Channel[1]\"." );
                }
            } else {
                // Extract the index
                indexNumberStr = channelName.substr( startBracketIndex + 1, endBracketIndex - startBracketIndex - 1 );
                // Strip off the indexing brackets from the name
                channelName.resize( startBracketIndex );
            }
        }

        // Apply the position data type if provided
        if( positionTypeHint != channels::data_type_invalid && channelName == _T("Position") ) {
            mapping.type = positionTypeHint;
        }

        if( indexNumberStr.empty() ) {
            indexNumber = -1;
        } else {
            try {
                indexNumber = boost::lexical_cast<int>( indexNumberStr );
            } catch( boost::bad_lexical_cast e ) {
                if( !noError ) {
                    throw frantic::invalid_particle_file_exception()
                        << "csv_particle_istream: The input file \"" << frantic::strings::to_string( file )
                        << "\" has an invalid column header name, \"" << columnHeaders[i]
                        << "\".  The number inside the index brackets is incorrect, it must be a valid integer.)";
                }
            }

            if( !noError && indexNumber < 0 ) {
                throw frantic::invalid_particle_file_exception()
                    << "csv_particle_istream: The input file \"" << frantic::strings::to_string( file )
                    << "\" has an invalid column header name, \"" << columnHeaders[i]
                    << "\".  The number inside the index brackets is incorrect, it must be a nonnegative integer.";
            }
        }

        if( indexNumber >= 0 ) {
            mapping.hasIndex = true;
            mapping.index = size_t( indexNumber );
        }

        if( !frantic::channels::is_valid_channel_name( channelName ) ) {
            if( noError ) {
                continue;
            } else {
                throw frantic::invalid_particle_file_exception()
                    << "csv_particle_istream: The input file \"" << frantic::strings::to_string( file )
                    << "\" has an invalid column header name, \"" << columnHeaders[i]
                    << "\".  Names must have alphanumeric characters and \'_\' only.";
            }
        }

        // Keep track of the column names in order so we can preserve the ordering from the CSV (somewhat)
        if( columnMappings.find( channelName ) == columnMappings.end() )
            columnNames.push_back( channelName );
        columnMappings[channelName].push_back( mapping );
    }

    frantic::channels::channel_map channelMap;

    for( vector<frantic::tstring>::const_iterator it = columnNames.begin(); it != columnNames.end(); ++it ) {

        const frantic::tstring& channelName = *it;
        std::vector<column_mapping>& channelColumns = columnMappings[channelName];

        frantic::channels::data_type_t typeGuess = frantic::channels::data_type_invalid;
        std::set<size_t> indexMapping;

        for( size_t j = 0; j < channelColumns.size(); ++j ) {
            if( typeGuess == frantic::channels::data_type_invalid ) {
                typeGuess = channelColumns[j].type;
            } else if( !noError ) {
                if( channelColumns[j].type != frantic::channels::data_type_invalid &&
                    channelColumns[j].type != typeGuess ) {
                    throw frantic::invalid_particle_file_exception()
                        << "csv_particle_istream: Data type for column \"" << columnHeaders[channelColumns[j].column]
                        << "\" does not match previous column for channel \""
                        << frantic::strings::to_string( channelName ) << "\".";
                }
            }

            if( channelColumns[j].hasIndex ) {
                // Duplicate indices on the same channel is an error
                if( !noError && indexMapping.find( channelColumns[j].index ) != indexMapping.end() ) {
                    throw frantic::invalid_particle_file_exception()
                        << "csv_particle_istream: Duplicate index " << channelColumns[j].index << " for channel \""
                        << frantic::strings::to_string( channelName ) << "\" found in column " << j << ".";
                }
                indexMapping.insert( channelColumns[j].index );
            } else if( !noError && !indexMapping.empty() ) {
                throw frantic::invalid_particle_file_exception()
                    << "csv_particle_istream: Mixed index specification found for channel \""
                    << frantic::strings::to_string( channelName )
                    << "\". Indices must be specified either for all columns of a given channel, or none.";
            }
        }

        if( typeGuess == frantic::channels::data_type_invalid ) {
            typeGuess = frantic::channels::data_type_float32;
        }

        // Assign ordering if not given, or attempt to salvage an ordering in noError mode
        if( noError || indexMapping.empty() ) {
            size_t currentIndex = 0;
            for( size_t j = 0; j < channelColumns.size(); ++j ) {
                if( !channelColumns[j].hasIndex ) {
                    while( indexMapping.find( currentIndex ) != indexMapping.end() ) {
                        ++currentIndex;
                    }
                    channelColumns[j].index = currentIndex;
                    indexMapping.insert( currentIndex );
                    ++currentIndex;
                }
            }
        }

        std::stable_sort( channelColumns.begin(), channelColumns.end() );

        // attempt to stably realign duplicates if we don't want to just error out
        if( noError ) {
            for( size_t j = 0; j < channelColumns.size(); ++j ) {
                channelColumns[j].index = j;
            }
        } else {
            bool oneBased = channelColumns[0].index == 1;

            for( size_t j = 0; j < channelColumns.size(); ++j ) {
                if( oneBased ) {
                    --channelColumns[j].index;
                }
                if( channelColumns[j].index != j ) {
                    throw frantic::invalid_particle_file_exception()
                        << "csv_particle_istream: Index specifications for channel \""
                        << frantic::strings::to_string( channelName ) << "\" were not a contiguous 1 or 0-based list.";
                }
            }
        }

        channelMap.define_channel( channelName, channelColumns.size(), typeGuess );
    }

    channelMap.end_channel_definition();
    outGuess.reset( channelMap, columnHeaders.size() );

    for( std::map<frantic::tstring, std::vector<column_mapping>>::iterator it = columnMappings.begin();
         it != columnMappings.end(); ++it ) {

        // Since we already sorted each channel's columns by index, we can assign them in order
        for( size_t j = 0; j < it->second.size(); ++j ) {
            if( it->second[j].type == frantic::channels::data_type_invalid ) {
                if( !noError ) {
                    throw std::logic_error( "csv_particle_istream: an invalid data type was detected in column " +
                                            boost::lexical_cast<std::string>( it->second[j].column ) );
                }
            }
            outGuess.set_column_mapping( it->second[j].column, it->first, it->second[j].index );
        }
    }
}

} // anonymous namespace

csv_particle_istream::csv_format_guess
csv_particle_istream::guess_channel_column_map( const std::vector<std::string>& columnHeaders,
                                                const frantic::tstring& file, bool noError,
                                                data_type_t positionTypeHint, channel_column_map& outMapping ) {
    bool anyNumbers = false, allNumbers = true;
    for( unsigned i = 0; i < columnHeaders.size(); ++i ) {
        std::string entry = frantic::strings::trim( columnHeaders[i] );
        const char* str = entry.c_str();
        char* end = 0;
        frantic::locale::strtod_c( str, &end );
        if( end == str + entry.size() )
            anyNumbers = true;
        else
            allNumbers = false;
    }

    csv_format_guess guessType = CSV_PARTICLE_GUESS_UNKNOWN;

    if( allNumbers ) {
        guess_native_channel_map_from_numeric_headers( columnHeaders, noError, positionTypeHint, outMapping );
        guessType = CSV_PARTICLE_GUESS_NUMBERS;
    } else if( anyNumbers ) {
        if( !noError ) {
            throw frantic::invalid_particle_file_exception()
                << "csv_particle_istream: The input file \"" + frantic::strings::to_string( file ) +
                       "\" is not a valid particle .csv file.  The column header must not contain any numeric values.";
        }
        guessType = CSV_PARTICLE_GUESS_MIXED;
    } else {
        guess_native_channel_map_from_text_headers( columnHeaders, file, noError, positionTypeHint, outMapping );
        guessType = CSV_PARTICLE_GUESS_HEADERS;
    }

    return guessType;
}

void csv_particle_istream::initialize_stream( const frantic::tstring& file, const channel_column_map& expectedFormat,
                                              char delimiter, int headerRowCount, data_type_t positionTypeHint,
                                              bool noError ) {
    m_name = file;
    m_fin.reset( frantic::files::tfopen( file.c_str(), _T("rb") ) );

    if( !m_fin )
        throw std::runtime_error( "csv_particle_istream: Failed to open file \"" + frantic::strings::to_string( file ) +
                                  "\" for reading." );

    if( delimiter == '\0' ) {
        bool foundDelimiter = guess_csv_delimiter( m_fin, file, m_delimiter );
        if( !foundDelimiter ) {
            m_delimiter = ',';
            FF_LOG( debug ) << "csv_particle_istream: In file \"" << name() << "\": Unable to find delimiter.  Using \'"
                            << boost::lexical_cast<frantic::tstring>( m_delimiter ) << "\'." << std::endl;
        } else {
            FF_LOG( debug ) << "csv_particle_istream: In file \"" << name() << "\": Using delimiter \'"
                            << boost::lexical_cast<frantic::tstring>( m_delimiter ) << "\'." << std::endl;
        }
        std::fseek( m_fin, 0, SEEK_SET );
    } else {
        m_delimiter = delimiter;
    }

    // Load in the header
    std::vector<std::string> columnHeaders;
    m_csvReader = get_csv_reader( m_fin, file, false, m_delimiter );
    m_csvReader.read_line( columnHeaders );

    if( columnHeaders.empty() )
        throw frantic::invalid_particle_file_exception()
            << "csv_particle_istream: The input file \"" + frantic::strings::to_string( file ) +
                   "\" is not a valid particle .csv file.  The first line must not be empty.";

    if( expectedFormat.is_empty() ) {
        csv_format_guess guessType =
            guess_channel_column_map( columnHeaders, file, noError, positionTypeHint, m_channelColumnMapping );
        if( guessType == CSV_PARTICLE_GUESS_NUMBERS ) {
            headerRowCount = 0;
        } else {
            headerRowCount = 1;
        }
    } else {
        m_channelColumnMapping = expectedFormat;
        m_channelColumnMapping.remove_trailing_columns();

        if( m_channelColumnMapping.column_count() > columnHeaders.size() ) {
            throw frantic::invalid_particle_file_exception()
                << "csv_particle_istream: Suggested particle format has more columns than in file \"" +
                       frantic::strings::to_string( file ) + "\".  "
                << "Suggested format has "
                << boost::lexical_cast<std::string>( m_channelColumnMapping.column_count() ) + ", "
                << "but file has only " + boost::lexical_cast<std::string>( columnHeaders.size() ) + ".";
        }
    }

    if( headerRowCount <= 0 ) {
        std::clearerr( m_fin );
    }

    m_headerRowCount = headerRowCount;
    m_fin.close();
    m_finReopened = false;
}

void csv_particle_istream::initialize_stream( const frantic::tstring& file, char delimiter,
                                              data_type_t positionTypeHint, bool noError ) {
    initialize_stream( file, channel_column_map(), delimiter, -1, positionTypeHint, noError );
}

} // namespace streams
} // namespace particles
} // namespace frantic
