// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/strings/to_string_classic.hpp>
#include <frantic/strings/tstring.hpp>

#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>

#include <cctype>
#include <cstdio>

#if defined( _MSC_VER )
#pragma warning( push )
// TODO: this is to get rid of all the 'deprecated' warnings about getenv and such
#pragma warning( disable : 4996 )
#endif


namespace frantic {
namespace strings {

/**
 * This function should read a single line from a file, delimited by a newline or EOF. It is intended to match the
 * behaviour of the C++ std::getline template function, only for cstdio.
 *
 * If this function returns false, the caller should check ferror(in) to see whether it was just an EOF or if there was
 * an error.
 *
 * @param  in  The C FILE* pointer.
 * @param  outString  The output string in which to place the line.
 */
inline bool getline( FILE* in, std::string& outString ) {
    int c;

    outString.clear();
    c = std::fgetc( in );
    if( c == '\n' )
        return true;
    if( c == EOF )
        return false;
    do {
        outString += static_cast<std::string::value_type>( c );
        c = std::fgetc( in );
        if( c == '\n' || c == EOF )
            return true;
    } while( 1 );
}

/**
 * Converts an integer that is in byte units to a human-readable string with a unit suffix.
 * @param sizeInBytes An integer to convert.
 */
inline frantic::tstring int_bytes_to_string( size_t sizeInBytes ) {
    std::string results = "";
    std::string suffixes[4];
    suffixes[0] = " bytes";
    suffixes[1] = " KB"; // technically a KiB -Paul.
    suffixes[2] = " MB"; // technically a MiB -Paul.
    suffixes[3] = " GB"; // technically a GiB -Paul.

    double val = (double)sizeInBytes;
    int division = 0;
    bool done = false;

    while( !done ) {
        double divided = val / 1024.0f;
        if( divided < 1 || division == 3 ) {
            std::stringstream stream;
            val = (int)( val * 1000 ) / 1000;
            stream << val << suffixes[division];
            results = stream.str();
            done = true;
        } else {
            val = divided;
            division++;
        }
    }
    return frantic::strings::to_tstring( results );
}

/**
 * This function takes in an int and turns it into a frantic::tstring that is split by commas
 * eg. int_to_comma_seperated_string(1234) = "1,234"
 *
 * @param initialVal the int that will be seperated
 *
 */
inline frantic::tstring int_to_comma_seperated_string( int integerValue ) {
    std::string results = "";

    bool isNeg = false;
    int val = integerValue;
    if( val == 0 ) {
        results = "0";
    } else {
        if( val < 0 ) {
            isNeg = true;
            val = val * -1;
        }
        int precommaVal = val;
        int postcommaVal = 0;
        while( precommaVal > 0 ) {
            std::stringstream stream;
            precommaVal = val / 1000;
            postcommaVal = val - ( 1000 * precommaVal );
            stream << postcommaVal << results;
            results = stream.str();
            if( precommaVal > 0 ) {
                if( postcommaVal < 100 ) {
                    results = "0" + results;
                    if( postcommaVal < 10 ) {
                        results = "0" + results;
                    }
                }
                results = "," + results;
            }
            val = precommaVal;
        }
        if( isNeg ) {
            results = "-" + results;
        }
    }
    return frantic::strings::to_tstring( results );
}

inline bool is_numeric( const std::string& str ) {
    const char* pChar = str.c_str();

    while( *pChar != '\0' ) {
        if( !isdigit( *pChar ) )
            return false;
        pChar++;
    }

    return true;
}

inline std::string to_lower( std::string str ) {
    for( unsigned i = 0; i < str.size(); ++i )
        str[i] = (char)tolower( str[i] );

    return str;
}

// TODO: handle this properly, maybe using something like boost::locale
inline std::wstring to_lower( std::wstring str ) {
    for( unsigned i = 0; i < str.size(); ++i )
        str[i] = (wchar_t)towlower( str[i] );

    return str;
}

inline std::string to_upper( std::string str ) {
    for( unsigned i = 0; i < str.size(); ++i )
        str[i] = (char)toupper( str[i] );

    return str;
}

inline bool starts_with( const std::string& str, const std::string& start ) {
    if( str.size() >= start.size() && str.substr( 0, start.size() ) == start )
        return true;
    return false;
}

inline bool ends_with( const std::string& str, const std::string& end ) {
    if( str.size() >= end.size() && str.substr( str.size() - end.size() ) == end )
        return true;
    return false;
}

inline std::string remove_comment( const std::string& str, const std::string& commentSeparator ) {
    std::string::size_type commentPos = str.find( commentSeparator );

    if( commentPos != std::string::npos )
        return str.substr( 0, commentPos );
    else
        return str;
}

inline char to_hex_char( int i ) {
    if( 0 <= i && i < 10 )
        return char( '0' + i );
    else if( i < 16 )
        return char( 'a' + i - 10 );
    else
        throw std::runtime_error( "Tried to convert out of range value " + boost::lexical_cast<std::string>( i ) +
                                  " to a hexadecimal digit" );
}

inline std::string to_hex_string( char c ) {
    int c1 = ( (unsigned char)c ) / 16, c2 = ( (unsigned char)c ) % 16;
    char result[3];
    result[0] = to_hex_char( c1 );
    result[1] = to_hex_char( c2 );
    result[2] = '\0';
    return result;
}

inline std::string to_hex_string( const char* buffer, int count ) {
    std::string result;
    for( int i = 0; i < count; ++i ) {
        result += to_hex_string( buffer[i] );
        if( ( i + 1 ) % 4 == 0 && i != count - 1 )
            result += " ";
    }
    return result;
}

template <class T>
inline std::string to_hex_string( const T& value ) {
    return to_hex_string( reinterpret_cast<const char*>( &value ), sizeof( T ) );
}

inline bool equals_ignore_case( const std::string& s1, const std::string& s2 ) {
    std::string lower1( to_lower( s1 ) );
    std::string lower2( to_lower( s2 ) );
    return lower1 == lower2;
}

#ifdef _WIN32
inline bool test_environment_variable_name( std::string& name, bool harmonizeCase = false ) {
    // This function tests a variable name to see if it exists in the
    // current environment.  If it does, the method can optionally update
    // the case of the name parameter to match the case of the environment
    // variable

    std::string::size_type equalSign;
    std::string variableName;
    std::string environmentString;
    bool success = false;
    char** search = _environ; // pointer used to traverse environment variable list.

    while( *search && !success ) {
        environmentString = *search;
        equalSign = environmentString.find( '=' );
        if( equalSign == std::string::npos )
            equalSign = environmentString.size();

        variableName = environmentString.substr( 0, equalSign );
        success = equals_ignore_case( variableName, name );
        if( success && harmonizeCase )
            name = variableName;

        ++search;
    }
    return success;
}
#else
inline bool test_environment_variable_name( std::string& name, bool /*harmonizeCase*/ = false ) {
    // In unix, the harmonize case functionality is unneeded because environment variables are
    // case sensitive.  Just return whether the environment variable exists
    return getenv( name.c_str() ) != 0;
}
#endif

/**
 * This function takes in an integer and a desired number of characters to pad to.
 * The function returns a string representing the original integer padded to at least length size
 * with zeros. If size is less than the number of characters needed to represent number as a string
 * the resulting string will be longer than size to account for this.
 * No zeros will be added if the number needs size or smaller characters to be represented.
 *
 * @param number the integer to be represented by the resulting string
 * @param size the minimum number of characters to pad the resulting string to with zeros
 * @return a string representing number padded to at least size with zeros
 *
 */
inline frantic::tstring zero_pad( int number, std::size_t size ) {
    const int absNumber = std::abs( number );
    const frantic::tstring absAsString = frantic::strings::to_string_classic<frantic::tstring>( absNumber );
    const int signLength = number < 0;
    const int numZeros =
        std::max( 0, static_cast<int>( size ) - ( static_cast<int>( absAsString.size() ) + signLength ) );
    const frantic::tstring zeros( numZeros, '0' );
    const frantic::tstring sign = number < 0 ? _T( "-" ) : _T( "" );
    return sign + zeros + absAsString;
}

template <class CharType>
inline std::basic_string<CharType> trim( std::basic_string<CharType> s ) {
    int nonSpaceIndex = 0;

    // left trim
    while( nonSpaceIndex < (int)s.size() && ( s[nonSpaceIndex] == ' ' || s[nonSpaceIndex] == '\t' ) )
        nonSpaceIndex++;
    if( nonSpaceIndex == (int)s.size() )
        return std::basic_string<CharType>();
    if( nonSpaceIndex > 0 )
        s = s.substr( nonSpaceIndex );

    nonSpaceIndex = (int)s.size() - 1;
    while( nonSpaceIndex >= 0 && ( s[nonSpaceIndex] == ' ' || s[nonSpaceIndex] == '\t' ) )
        nonSpaceIndex--;
    if( nonSpaceIndex < (int)s.size() - 1 )
        s = s.substr( 0, nonSpaceIndex + 1 );

    return s;
}

// escapes special characters
template <class CharType>
inline std::basic_string<CharType> get_escaped_string( const std::basic_string<CharType>& s ) {
    std::basic_string<CharType> result;
    int add_from = 0;
    CharType esc_char = 0;
    for( int i = 0; i < (int)s.size(); ++i ) {
        switch( s[i] ) {
        case '\n':
            esc_char = 'n';
            break;
        case '\r':
            esc_char = 'r';
            break;
        case '"':
            esc_char = '"';
            break;
        case '\\':
            esc_char = '\\';
            break;
        case '\t':
            esc_char = 't';
            break;
        default:
            esc_char = '\0';
            break;
        }
        if( esc_char != 0 ) {
            result += s.substr( add_from, i - add_from );
            add_from = i + 1;
            result += '\\';
            result += esc_char;
        }
    }
    result += s.substr( add_from, s.size() - add_from );

    return result;
}

// Puts a string in quotes, and escapes special characters
inline std::string get_quoted_string( const std::string& s ) { return "\"" + get_escaped_string( s ) + "\""; }

// This does the same as what tokenize_string did.  Splits the input string into a series of tokens separated by one or
// more splitChar in between.
template <class CharType>
inline void split( const std::basic_string<CharType>& str, std::vector<std::basic_string<CharType>>& output,
                   CharType splitChar ) {
    typename std::basic_string<CharType>::size_type index = 0, splitIndex;
    while( ( splitIndex = str.find( splitChar, index ) ) != std::basic_string<CharType>::npos ) {
        output.push_back( str.substr( index, splitIndex - index ) );
        index = splitIndex + 1;
    }
    if( index < str.size() )
        output.push_back( str.substr( index ) );
}

namespace detail {
// This does the same as what tokenize_string did.  Splits the input string into a series of tokens separated by one or
// more members
// of splitChars in between.
template <class CharType>
inline void split( const std::basic_string<CharType>& str, std::vector<std::basic_string<CharType>>& output,
                   const std::basic_string<CharType>& splitChars ) {
    typename std::basic_string<CharType>::size_type index = 0, split_index = 0;
    bool inSplitter = true;
    for( index = 0; index < str.size(); ++index ) {
        if( splitChars.find( str[index] ) != std::string::npos ) {
            if( !inSplitter ) {
                output.push_back( str.substr( split_index, index - split_index ) );
                inSplitter = true;
            }
        } else {
            if( inSplitter ) {
                split_index = index;
                inSplitter = false;
            }
        }
    }
    if( !inSplitter ) {
        output.push_back( str.substr( split_index, str.size() - split_index ) );
        inSplitter = true;
    }
}
} // namespace detail

// This does the same as what tokenize_string did.  Splits the input string into a series of tokens separated by one or
// more members
// of splitChars in between.
inline void split( const std::string& str, std::vector<std::string>& output, const std::string& splitChars = "\n\t " ) {
    detail::split<char>( str, output, splitChars );
}

inline void split( const std::wstring& str, std::vector<std::wstring>& output,
                   const std::wstring& splitChars = L"\n\t " ) {
    detail::split<wchar_t>( str, output, splitChars );
}

// Splits apart a search path separated by ; or : characters.  Takes special care to
// navigate around drive letters
inline void split_search_path( const frantic::tstring& str, std::vector<frantic::tstring>& output ) {
    frantic::tstring::size_type index = 0, split_index = 0;

    bool inSplitter = true;
    bool potentialDriveLetter = false;
    for( index = 0; index < str.size(); ++index ) {
        potentialDriveLetter = index - split_index == 1;
        if( ( str[index] == _T( ';' ) || str[index] == _T( ':' ) ) && !potentialDriveLetter ) {
            if( !inSplitter ) {
                output.push_back( str.substr( split_index, index - split_index ) );
                inSplitter = true;
            }
        } else {
            if( inSplitter ) {
                split_index = index;
                inSplitter = false;
            }
        }
    }
    if( !inSplitter ) {
        output.push_back( str.substr( split_index, str.size() - split_index ) );
        inSplitter = true;
    }
}

namespace detail {
template <class CharType>
std::basic_string<CharType> string_replace( std::basic_string<CharType> input,
                                            const std::basic_string<CharType>& target,
                                            const std::basic_string<CharType>& replacement ) {
    typename std::basic_string<CharType>::size_type location = 0;
    while( ( location = input.find( target, location ) ) != std::basic_string<CharType>::npos ) {
        input.replace( location, target.size(), replacement );
        location += replacement.size();
    }
    return input;
}
} // namespace detail

inline std::string string_replace( std::string input, const std::string& target, const std::string& replace ) {
    return detail::string_replace( input, target, replace );
}

inline std::wstring string_replace( std::wstring input, const std::wstring& target, const std::wstring& replace ) {
    return detail::string_replace( input, target, replace );
}

inline std::string char_replace( std::string input, char target, char replacement ) {
    std::string::size_type location = 0;
    while( ( location = input.find( target, location ) ) != std::string::npos ) {
        input[location] = replacement;
        location++;
    }
    return input;
}

// Determines whether or not the character is valid for
// an identifier.
inline bool is_valid_identifier_character( const char& ch ) { return isalnum( ch ) || ch == '_'; }

inline bool is_valid_identifier_character( const wchar_t& ch ) { return iswalnum( ch ) || ch == '_'; }

// makes a valid identifier from an arbitrary string by replacing
// all invalid characters with underscores.
// TODO: handle surrogate characters in wstring
template <class CharType>
inline std::basic_string<CharType> make_valid_identifier( const std::basic_string<CharType>& str ) {
    std::basic_string<CharType> result;
    for( typename std::basic_string<CharType>::const_iterator i( str.begin() ); i != str.end(); ++i ) {
        if( is_valid_identifier_character( *i ) )
            result += *i;
        else
            result += CharType( '_' );
    }
    return result;
}

// inline bool less_ignore_case(const std::string & s1, const std::string & s2)
//{
//	return to_lower(s1) < to_lower(s2);
//}

struct less_ignore_case {
    bool operator()( const std::string& s1, const std::string& s2 ) const { return to_lower( s1 ) < to_lower( s2 ); }
};

// returns the the offset of the last character that is valid for an identifier
// (variable name, environment variable name) starting at the offset id_start
// std::string::npos is returned if the character at start is invalid.
inline std::string::size_type find_end_of_identifier( const std::string& str, std::string::size_type id_start = 0 ) {
    bool valid;
    std::string::size_type end = std::string::npos;

    if( is_valid_identifier_character( str[id_start] ) ) {
        end = id_start;
        valid = true;
        while( valid && ++end < str.size() ) {
            valid = is_valid_identifier_character( str[end] );
        }
        end--;
    }
    return end;
}

// finds the next dos environment variable starting at "start"
// returns the variable name along with the start and end indicies.
// the boolean return value indicates success.
inline bool find_next_dos_environment_variable( std::string input, std::string::size_type& start,
                                                std::string::size_type& end, std::string& variableName ) {
    std::string::size_type v_start = 0, v_end = 0;
    std::string variable;
    bool found = false;
    v_start = input.find( '%', start );
    if( v_start != std::string::npos ) {
        v_end = input.find( '%', v_start + 1 );
        while( !found && v_end != std::string::npos ) {
            // get the variable name from inside the % delimiters
            variable = input.substr( v_start + 1, v_end - v_start - 1 );

            // test the variable for existance in the environment
            // and harmonize the case with the environment.
            if( test_environment_variable_name( variable, true ) ) {
                found = true;
                start = v_start;
                end = v_end;
                variableName = variable;
            } else {
                v_start = v_end;
                v_end = input.find( '%', v_start + 1 );
            }
        }
    }
    return found;
}

inline bool find_next_unix_environment_variable( std::string input, std::string::size_type& start,
                                                 std::string::size_type& end, std::string& variableName ) {
    std::string::size_type ue_start, ue_end;        // unix env variable start and end
    std::string::size_type v_start = 0, v_size = 0; // variable name start and size
    bool found = false;

    // find the first $
    ue_start = input.find( '$', start );
    while( !found && ue_start != std::string::npos && ue_start < ( input.size() - 1 ) ) {
        // determine the start and size of the variable name as well
        // as the start and end of the variable (including special chars)
        if( input[ue_start + 1] == '{' ) {
            if( ( ue_start + 1 ) == ( input.size() - 1 ) ) {
                // if the open brace is at the end
                // of the string, we must fail
                ue_end = std::string::npos;
            } else {
                // handle names in braces
                v_start = ue_start + 2;
                ue_end = input.find( '}', v_start );
                v_size = ue_end - v_start;
            }
        } else {
            // handle names outside of braces.
            v_start = ue_start + 1;
            ue_end = find_end_of_identifier( input, v_start );
            v_size = ue_end - ue_start;
        }

        if( ue_end != std::string::npos && v_size > 0 ) {
            // if we've found a valid variable, replace it
            variableName = input.substr( v_start, v_size );
            start = ue_start;
            end = ue_end;
            found = true;
        } else {
            // no variable found.  Keep searching
            ue_start = input.find( '$', ue_start + 1 );
        }
    }
    return found;
}

inline std::string change_environment_variables_to_dos( std::string input ) {
    std::string::size_type start, end;
    std::string variable;

    start = 0;
    while( find_next_unix_environment_variable( input, start, end, variable ) ) {
        input.replace( start, end - start + 1, "%" + variable + "%" );
    }
    return input;
}

inline std::string change_environment_variables_to_unix( std::string input ) {
    std::string::size_type start, end; // the positions of the start and end '%'
    std::string variable;

    // first test existing unix-style variables to make sure the case is correct
    start = 0;
    while( find_next_unix_environment_variable( input, start, end, variable ) ) {
        if( test_environment_variable_name( variable, true ) ) {
            if( input[start + 1] == '{' )
                input.replace( start + 2, end - start - 2, variable );
            else
                input.replace( start + 1, end - start, variable );
        }
        start = end + 1;
    }

    start = 0;
    while( find_next_dos_environment_variable( input, start, end, variable ) ) {
        if( end != ( input.size() - 1 ) && is_valid_identifier_character( input[end + 1] ) ) {
            // if the next character is valid for an identifier, we must add braces
            input.replace( start, end - start + 1, "${" + variable + "}" );
            start = end + 2;
        } else {
            input.replace( start, end - start + 1, "$" + variable );
            start = end;
        }
    }
    return input;
}

inline std::string substitute_environment_variables( std::string input ) {
    std::string variableName;
    std::string::size_type start, end;
    char* envVar;
    std::string replacement;

    // convert all valid unix environment variables
    start = 0;
    while( find_next_unix_environment_variable( input, start, end, variableName ) ) {
        envVar = getenv( variableName.c_str() );
        if( envVar ) {
            replacement = envVar;
            input.replace( start, end - start + 1, replacement );
            start = start + replacement.size();
        } else {
            start = end + 1;
        }
    }

    // convert all valid dos environment variables
    start = 0;
    while( find_next_dos_environment_variable( input, start, end, variableName ) ) {
        envVar = getenv( variableName.c_str() );
        if( envVar ) {
            replacement = envVar;
            input.replace( start, end - start + 1, replacement );
            start = start + replacement.size();
        } else {
            start = end + 1;
        }
    }

    return input;
}

inline bool contains_environment_variable( const std::string input, std::string::size_type start,
                                           std::string::size_type end ) {
    std::string temp;
    std::string::size_type s2 = start; // the functions below modify start and end
    std::string::size_type e2 = end;   // so we use copies.
    return find_next_dos_environment_variable( input, start, end, temp ) ||
           find_next_unix_environment_variable( input, s2, e2, temp );
}

// Makes sure the command (the first part of the string) is in quotes.
// Example usage:
//  substitute_environment_variables( ensure_quoted_string( "$GELATOHOME/bin/iv.exe myimage.tif" ) )
// Result:
// "\"C:\\Program Files\\NVIDIA Corporation\\Gelato/bin/iv.exe\" myimage.tif"

// will only quote the command if it contains an environment variable
inline std::string ensure_quoted_command( const std::string& input ) {
    using namespace std;
    string result;
    string::size_type i = 0;
    string::size_type start, size;

    // Get to the first non-space
    while( i < input.size() && std::isspace( input[i] ) )
        ++i;
    result += input.substr( 0, i );
    start = i;

    // Add the quote if necessary
    if( i >= input.size() || input[i] == '"' ) {
        return input;
    }

    while( i < input.size() && !std::isspace( input[i] ) )
        ++i;
    size = i - start;

    if( contains_environment_variable( input, start, i ) ) {
        result += "\"";
        result += input.substr( start, size );
        result += "\"";
        result += input.substr( i );
        return result;
    } else
        return input;
}

namespace detail {
inline void substitute_chars_in_command( std::string& input, char searchChar, char replaceChar ) {
    using namespace std;
    // Make sure to use \ characters
    string::size_type i = 0;
    // Get to the first non-space
    while( i < input.size() && std::isspace( input[i] ) )
        ++i;
    if( i < input.size() && input[i] == '"' ) {
        ++i;
        // Process until the closing quote
        while( i < input.size() && input[i] != '"' ) {
            if( input[i] == searchChar )
                input[i] = replaceChar;
            ++i;
        }
    } else {
        // Process until the next space
        while( i < input.size() && !std::isspace( input[i] ) ) {
            if( input[i] == searchChar )
                input[i] = replaceChar;
            ++i;
        }
    }
}
} // namespace detail

inline std::string to_dos_command( const std::string& command ) {
    std::string result = change_environment_variables_to_dos( ensure_quoted_command( command ) );
    // Make it use the '\\' character
    detail::substitute_chars_in_command( result, '/', '\\' );
    return result;
}

inline std::string to_unix_command( const std::string& command ) {
    std::string result = change_environment_variables_to_unix( ensure_quoted_command( command ) );
    // Make it use the '/' character
    detail::substitute_chars_in_command( result, '\\', '/' );
    return result;
}

inline std::string to_command( const std::string& command ) {
    std::string result = substitute_environment_variables( ensure_quoted_command( command ) );
    detail::substitute_chars_in_command( result, '/', '\\' );
    return result;
    //	return detail::add_command_extension( result );
}

// Slightly less accurate than atof, but MUCH FASTER!
inline const char* parse_one_float( const char* str, float& a, bool& succeeded ) {
    succeeded = false;
    const char* s = str;
    bool negative = false;
    float result = 0;
    while( isspace( *s ) )
        ++s;
    // + or -
    if( *s == '+' || *s == '-' )
        negative = ( *s++ == '-' );
    // ####
    if( !isdigit( *s ) ) {
        if( *s != '.' )
            return str;
    } else {
        result += ( ( *s++ ) - '0' );
        while( isdigit( *s ) ) {
            result *= 10;
            result += ( ( *s++ ) - '0' );
        }
    }
    // .#####
    if( *s == '.' ) {
        ++s;
        float fraction = 0.1f;
        while( isdigit( *s ) ) {
            result += fraction * ( ( *s++ ) - '0' );
            fraction *= 0.1f;
        }
        /*
        // The following produced greater differences than the above
        float num = 0, denom = 1;
        while( isdigit(*s) ) {
          num += ((*s++) - '0');
          num *= 10;
          denom *= 10;
        }
        result += num / denom;
        */
    }
    if( negative )
        result = -result;
    // e+-###
    if( *s == 'e' || *s == 'E' ) {
        ++s;
        int exponent = 0;
        bool exponentNegative = false;
        if( *s == '+' || *s == '-' )
            exponentNegative = ( *s++ == '-' );
        if( !isdigit( *s ) )
            return str;
        while( isdigit( *s ) ) {
            exponent *= 10;
            exponent += ( ( *s++ ) - '0' );
        }
        if( exponentNegative )
            exponent = -exponent;

        result *= powf( 10.f, static_cast<float>( exponent ) );
    }

    a = result;
    succeeded = true;
    return s;
}

inline bool parse_n_floats( const std::string& str, float* out, int count ) {
    bool succeeded = true;
    int index = 0;
    const char* strPtr = str.c_str();
    while( succeeded && index < count ) {
        strPtr = parse_one_float( strPtr, out[index], succeeded );
        ++index;
    }
    return succeeded;
}

// converts milliseconds to a pretty, human readable string
template <class CharType>
inline std::basic_string<CharType> ms_to_basic_string( unsigned long ms ) {
    std::basic_stringstream<CharType> str;
    int msr = ms % 1000;
    unsigned long s = ms / 1000;
    unsigned long m = s / 60;
    unsigned long h = m / 60;
    m = m % 60;
    s = s % 60;
    str.width( 2 );
    str.fill( '0' );
    str << h << "h ";
    str.width( 2 );
    str.fill( '0' );
    str << m << "m ";
    str.width( 2 );
    str.fill( '0' );
    str << s << ".";
    str.width( 3 );
    str.fill( '0' );
    str << msr << "s";
    return str.str();
}

// converts milliseconds to a pretty, human readable string
inline frantic::tstring ms_to_string( unsigned long ms ) { return ms_to_basic_string<frantic::tchar>( ms ); }
} // namespace strings
} // namespace frantic

#if defined( _MSC_VER )
#pragma warning( pop )
#endif
