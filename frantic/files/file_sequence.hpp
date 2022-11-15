// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

// TODO: Delete this file, its functionality is taken over by the filename_sequence class.

#include <limits>
#include <string>

#include <frantic/files/paths.hpp>
#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace files {

#ifdef _WIN32
// This function will take a filename pattern of prefix + #### + .extension where #### is any length of digits.
// It will scan the directory designated in the pattern for all files matching the pattern and differing in only
// the digits. It returns an interval consisting of the lowest and highest frame found.
//  For example:
//  If there are files in some directory c:\blah\, like test0004.tif, test0005.tif, test0007.tif
//  get_sequence_range( "c:\\blah\\test0000.tif" ) will return the pair (4,7).
inline std::pair<int, int> get_sequence_range( const frantic::tstring& inFilePattern ) {
    frantic::tstring ext;
    frantic::tstring prefix;
    int numDigits;
    int seqNum;

    int lowestFrame = ( std::numeric_limits<int>::max )(); // Arbitrarily high value
    int highestFrame = -1;

    HANDLE fsHandle;
    WIN32_FIND_DATA fsData;
    frantic::tstring searchPattern;

    split_sequence_path( inFilePattern, prefix, numDigits, seqNum, ext );

    searchPattern = prefix + _T("*") + ext;

    fsHandle = FindFirstFile( searchPattern.c_str(), &fsData );
    if( fsHandle != INVALID_HANDLE_VALUE ) {
        do {
            try {
                seqNum = extract_sequence_number( frantic::tstring( fsData.cFileName ) );
            } catch( const std::exception& ) {
                seqNum = lowestFrame;
            }

            if( seqNum < lowestFrame )
                lowestFrame = seqNum;
            else if( seqNum > highestFrame )
                highestFrame = seqNum;

        } while( FindNextFile( fsHandle, &fsData ) );
    }

    // if(lowestFrame > highestFrame)
    //	lowestFrame = highestFrame;

    return std::pair<int, int>( lowestFrame, highestFrame );
} //~get_sequence_range()

#else
inline std::pair<int, int> get_sequence_range( const frantic::tstring& /*inFilePattern*/ ) {
    throw std::runtime_error( "get_sequence_range has not been implemented for unix yet." );
}
#endif
/*framenumbers::frame_set get_frame_set_from_files(const std::string & inFilePattern){

}*/
} // namespace files
} // namespace frantic
