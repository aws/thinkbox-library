// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <map>
#include <string>

#include <frantic/channels/channel_map.hpp>

namespace frantic {
namespace channels {

namespace detail {
// This represents a block of data that must be transfered using a type converter function.
// For blocks that need to be memcpy'd, the convertor is a memcpy function, and the count
// represents a number of bytes instead of a number of elements.
struct cma_type_conversion_block {
    // The source and destination positions, in bytes.
    unsigned sourcePosition, destPosition;
    // The number of elements to convert.
    unsigned count;
    // The function which knows how to convert such a block.
    channel_type_convertor_function_t convertor;

    cma_type_conversion_block( unsigned sourcePosition_, unsigned destPosition_, unsigned count_,
                               channel_type_convertor_function_t convertor_ )
        : sourcePosition( sourcePosition_ )
        , destPosition( destPosition_ )
        , count( count_ )
        , convertor( convertor_ ) {}
};
} // namespace detail

/**
 * The particle channel map adaptor provides facilities for transfering data from one type of channel_map to another.
 * It analyzes the conversion from the source channel_map to the destination, creating as few block copies and type
 * conversions as possible.
 *
 * Type conversion rules:
 *  * Allowing any conversion between float types
 *  * Allowing any conversion between unsigned int types
 *  * Allowing any conversion between signed int types
 *  * Allowing conversion from an unsigned int to a bigger signed int type
 *  * Disallowing all other conversion types (like signed -> unsigned, int -> float, float->int, etc)
 *  * Strings can only be copied to strings
 *
 */
class channel_map_adaptor {
    std::size_t m_sourceSize, m_destSize;
    // These are all the blocks of data that need to be copied or converted
    std::vector<detail::cma_type_conversion_block> m_typeConversionBlocks;
    std::vector<detail::cma_type_conversion_block> m_defaultCopyBlocks; // first is the offset, second is block length
    bool m_isIdentity;

  public:
    /**
     * Constructs an empty channel_map_adaptor.
     */
    channel_map_adaptor() { clear(); }

    /**
     * The constructor computes the minimal number of memcpy and type conversion blocks required to convert a particle
     * from the source representation to the dest representation.
     *
     * @param  dest             The channel_map specifying the structure of the destination.
     * @param  source           The channel_map specifying the structure of the source.
     * @param  channelRenaming  A map which can be used to copy a channel of one name to a channel of another.
     */
    channel_map_adaptor( const channel_map& dest, const channel_map& source,
                         const std::map<frantic::tstring, frantic::tstring> channelRenaming =
                             ( std::map<frantic::tstring, frantic::tstring>() ) );

    /**
     * This function sets the adaptor, computing the minimal number of memcpy and type conversion blocks required to
     * convert a particle from the source representation to the dest representation.
     *
     * @param  dest             The channel_map specifying the structure of the destination.
     * @param  source           The channel_map specifying the structure of the source.
     * @param  channelRenaming  A map which can be used to copy a channel of one name to a channel of another.
     */
    void set( const channel_map& dest, const channel_map& source,
              const std::map<frantic::tstring, frantic::tstring> channelRenaming =
                  ( std::map<frantic::tstring, frantic::tstring>() ) );

    /**
     * Clears the adaptor to an empty state.
     */
    void clear() {
        m_typeConversionBlocks.clear();
        m_sourceSize = 0;
        m_destSize = 0;
        m_isIdentity = false;
    }

    /**
     * The number of bytes in the source structure.
     */
    std::size_t source_size() const { return m_sourceSize; }

    /**
     * The number of bytes in the destination structure.
     */
    std::size_t dest_size() const { return m_destSize; }

    /**
     * True if this adaptor will simply do a copy of the source structure to the destination.  This can be used
     * to skip the adaptor when an input can be read directly into the destination, instead of reading it into
     * a temporary buffer followed by conversion using the adaptor.
     */
    bool is_identity() const { return m_isIdentity; }

    /**
     * This function copies channel data from a source structure to a destination structure.  For items in the
     * destination structure that didn't have a corresponding item in the source structure,
     * the existing data is untouched.
     *
     * This function will not throw (but could crash if you pass invalid buffers as parameters).
     *
     * @param  destData    The location of the destination structure memory.
     * @param  sourceData  The location of the source structure memory.
     */
    void copy_structure( char* destData, const char* sourceData ) const;

    /**
     * @overload
     */
    void copy_structure( std::vector<char>& destData, const char* sourceData ) const {
        // TODO: Can do checking in a DEBUG mode.
        copy_structure( &destData[0], sourceData );
    }

    /**
     * @overload
     */
    void copy_structure( char* destData, const std::vector<char>& sourceData ) const {
        // TODO: Can do checking in a DEBUG mode.
        copy_structure( destData, &sourceData[0] );
    }

    /**
     * @overload
     */
    void copy_structure( std::vector<char>& destData, const std::vector<char>& sourceData ) const {
        // TODO: Can do checking in a DEBUG mode.
        copy_structure( &destData[0], &sourceData[0] );
    }

    /**
     * Copies the channels in destData not supplied by the src channel map. Takes the value from
     * the supplied default data.
     *
     * @param destData		The location of the destination structure memory.
     * @param defaultData	The location of the default structure memory.
     */
    void copy_structure_gaps_from_default( char* destData, const char* defaultData ) const;

    /**
     * Copies channels from srcData to destData, filling in the "gaps" in destData that was no supplied by
     * srcData with channels from defaultData. This implies defaultData should have the same channel mapping
     * as destData.
     *
     * @param destData		The location of the destination structure memory.
     * @param srcData		The location of the source structure memory.
     * @param defaultData	The location of the default structure memory.
     */
    void copy_structure( char* destData, const char* srcData, const char* defaultData ) const {
        copy_structure( destData, srcData );
        copy_structure_gaps_from_default( destData, defaultData );
    }

    /**
     * @overload
     */
    void copy_structure( std::vector<char>& destData, const std::vector<char>& srcData,
                         const std::vector<char>& defaultData ) const {
        copy_structure( &destData[0], &srcData[0] );
        copy_structure_gaps_from_default( &destData[0], &defaultData[0] );
    }

    /**
     * Dumps the adaptor, for debugging purposes.
     *
     * @param  out  The output stream that the adaptor gets dumped into.
     */
    void dump( std::ostream& out ) const;

    /**
     * This returns the number of block operations which are used to copy the structure.
     * It is used by the unit tests, to ensure that the right number of blocks get created.
     */
    std::size_t block_operation_count() const { return m_typeConversionBlocks.size(); }
};

} // namespace channels
} // namespace frantic
