// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/channels/channel_map_adaptor.hpp>

using namespace std;
using namespace boost;
using namespace frantic::graphics;

namespace {

// This is different from the conversion functions, in that the count parameter represents a
// byte count instead of an element count.
void convert_memcpy( char* out, const char* in, std::size_t byteCount ) { memcpy( out, in, byteCount ); }

} // anonymous namespace

namespace frantic {
namespace channels {

struct individual_element_copy {
    // The source and destination positions, in bytes.
    unsigned sourcePosition, destPosition;
    data_type_t sourceType, destType;
    frantic::tstring name;

    individual_element_copy( unsigned sourcePosition_, unsigned destPosition_, data_type_t sourceType_,
                             data_type_t destType_, const frantic::tstring& name_ )
        : sourcePosition( sourcePosition_ )
        , destPosition( destPosition_ )
        , sourceType( sourceType_ )
        , destType( destType_ )
        , name( name_ ) {}
};

// The constructor computes the minimal number of memcpy and type conversion blocks required to convert a particle from
// the source representation to the dest representation.
channel_map_adaptor::channel_map_adaptor( const channel_map& dest, const channel_map& source,
                                          const map<frantic::tstring, frantic::tstring> channelRenaming ) {
    set( dest, source, channelRenaming );
}

void channel_map_adaptor::set( const channel_map& dest, const channel_map& source,
                               const map<frantic::tstring, frantic::tstring> channelRenaming ) {
    // Reset the existing values just in case
    clear();

    if( !source.channel_definition_complete() )
        throw runtime_error(
            "channel_map_adaptor.set: The source channel_map definition process was not complete.  This "
            "must be finished before an adaptor can be created." );
    if( !dest.channel_definition_complete() )
        throw runtime_error(
            "channel_map_adaptor.set: The destination channel_map definition process was not complete.  "
            "This must be finished before an adaptor can be created." );

    m_sourceSize = source.structure_size();
    m_destSize = dest.structure_size();

    std::set<frantic::tstring> connectedChannels;
    std::vector<individual_element_copy> elementCopies;

    // First convert the mapping into a type-by-type copy representation
    for( unsigned channelNum = 0; channelNum < source.channel_count(); ++channelNum ) {
        const channel& pcSource = source[channelNum];
        frantic::tstring name;
        // Get the name to look for in the destination channel_map, possibly renamed using the provided map<>.
        map<frantic::tstring, frantic::tstring>::const_iterator i = channelRenaming.find( pcSource.name() );
        if( i == channelRenaming.end() )
            name = pcSource.name();
        else
            name = i->second;

        // If this channel has a place to go, then we add the element copies
        if( dest.has_channel( name ) ) {
            connectedChannels.insert( name );

            const channel& pcDest = dest[name];
            if( pcSource.arity() != pcDest.arity() )
                throw runtime_error( "channel_map_adaptor.set: Channel \"" +
                                     frantic::strings::to_string( pcSource.name() ) + "\" has arity " +
                                     boost::lexical_cast<std::string>( pcSource.arity() ) +
                                     " in the source structure and its destination channel, \"" +
                                     frantic::strings::to_string( pcDest.name() ) + "\",  has arity " +
                                     boost::lexical_cast<std::string>( pcDest.arity() ) +
                                     " in the destination structure, but they need to be the same for copying." );
            unsigned dataSizeSource = (unsigned)sizeof_channel_data_type( pcSource.data_type() ),
                     dataSizeDest = (unsigned)sizeof_channel_data_type( pcDest.data_type() );
            for( unsigned element = 0; element < pcSource.arity(); ++element ) {
                elementCopies.push_back(
                    individual_element_copy( unsigned( pcSource.offset() + element * dataSizeSource ),
                                             unsigned( pcDest.offset() + element * dataSizeDest ), pcSource.data_type(),
                                             pcDest.data_type(), pcSource.name() ) );
            }
        }
    }

    // Record the gaps which the source channel map did not fill. The default structure can fill this.
    bool previousDefaultCopyIsMemcpy = false;
    m_defaultCopyBlocks.clear();
    for( std::size_t i = 0; i < dest.channel_count(); ++i ) {
        if( connectedChannels.find( dest[i].name() ) == connectedChannels.end() ) {
            if( is_channel_data_type_pod( dest[i].data_type() ) ) {
                if( previousDefaultCopyIsMemcpy && ( m_defaultCopyBlocks.back().sourcePosition +
                                                     m_defaultCopyBlocks.back().count ) == dest[i].offset() ) {
                    m_defaultCopyBlocks.back().count += (unsigned)dest[i].primitive_size();
                } else {
                    m_defaultCopyBlocks.push_back(
                        detail::cma_type_conversion_block( (unsigned)dest[i].offset(), (unsigned)dest[i].offset(),
                                                           (unsigned)dest[i].primitive_size(), &convert_memcpy ) );
                    previousDefaultCopyIsMemcpy = true;
                }
            } else {
                m_defaultCopyBlocks.push_back( detail::cma_type_conversion_block(
                    (unsigned)dest[i].offset(), (unsigned)dest[i].offset(), (unsigned)dest[i].arity(),
                    get_channel_type_convertor_function( dest[i].data_type(), dest[i].data_type(), dest[i].name() ) ) );
                previousDefaultCopyIsMemcpy = false;
            }
        }
    }

    // Then run-length encode the resulting copies and conversions into a sequence of memcpy and type conversion blocks
    unsigned sourcePosition = 0, destPosition = 0;
    unsigned memcpyBytes = 0, typeConversionElements = 0;
    data_type_t sourceType = data_type_invalid, destType = data_type_invalid;
    std::size_t sourceTypeSize = 0, destTypeSize = 0;
    frantic::tstring typeConversionChannelName;
    for( unsigned i = 0; i < elementCopies.size(); ++i ) {
        // If we can memcpy this one, i.e. the types match and are POD (plain old data)
        if( elementCopies[i].sourceType == elementCopies[i].destType &&
            is_channel_data_type_pod( elementCopies[i].sourceType ) ) {
            // If we can combine this into an existing memcpy block, do it
            if( memcpyBytes > 0 && sourcePosition + memcpyBytes == elementCopies[i].sourcePosition &&
                destPosition + memcpyBytes == elementCopies[i].destPosition ) {
                memcpyBytes += (unsigned)sizeof_channel_data_type( elementCopies[i].sourceType );
            }
            // Otherwise close off any previous blocks, and start a new memcpy block
            else {
                if( memcpyBytes > 0 ) {
                    m_typeConversionBlocks.push_back( detail::cma_type_conversion_block(
                        sourcePosition, destPosition, memcpyBytes, &convert_memcpy ) );
                    memcpyBytes = 0;
                } else if( typeConversionElements > 0 ) {
                    m_typeConversionBlocks.push_back( detail::cma_type_conversion_block(
                        sourcePosition, destPosition, typeConversionElements,
                        get_channel_type_convertor_function( sourceType, destType, typeConversionChannelName ) ) );
                    typeConversionElements = 0;
                }
                // Start the new memcpy block
                sourcePosition = elementCopies[i].sourcePosition;
                destPosition = elementCopies[i].destPosition;
                memcpyBytes = (unsigned)sizeof_channel_data_type( elementCopies[i].sourceType );
            }
        }
        // This one requires a type conversion, or scope-managed copy
        else {
            // If we can combine this into an existing type conversion block, do it
            if( typeConversionElements > 0 && sourceType == elementCopies[i].sourceType &&
                destType == elementCopies[i].destType &&
                sourcePosition + typeConversionElements * sourceTypeSize == elementCopies[i].sourcePosition &&
                destPosition + typeConversionElements * destTypeSize == elementCopies[i].destPosition ) {
                typeConversionElements++;
            }
            // Otherwise close off any previous blocks, and start a new type conversion block
            else {
                if( memcpyBytes > 0 ) {
                    m_typeConversionBlocks.push_back( detail::cma_type_conversion_block(
                        sourcePosition, destPosition, memcpyBytes, &convert_memcpy ) );
                    memcpyBytes = 0;
                } else if( typeConversionElements > 0 ) {
                    m_typeConversionBlocks.push_back( detail::cma_type_conversion_block(
                        sourcePosition, destPosition, typeConversionElements,
                        get_channel_type_convertor_function( sourceType, destType, typeConversionChannelName ) ) );
                    typeConversionElements = 0;
                }
                // Start the new type conversion block
                typeConversionChannelName = elementCopies[i].name;
                sourcePosition = elementCopies[i].sourcePosition;
                destPosition = elementCopies[i].destPosition;
                sourceType = elementCopies[i].sourceType;
                destType = elementCopies[i].destType;
                sourceTypeSize = sizeof_channel_data_type( sourceType );
                destTypeSize = sizeof_channel_data_type( destType );
                typeConversionElements = 1;
            }
        }
    }
    // Close the last block.
    if( memcpyBytes > 0 ) {
        m_typeConversionBlocks.push_back(
            detail::cma_type_conversion_block( sourcePosition, destPosition, memcpyBytes, &convert_memcpy ) );
        memcpyBytes = 0;
    } else if( typeConversionElements > 0 ) {
        //		cout << "Source type is: " << channel_data_type_str(sourceType) << "\n";
        //		cout << "Dest type is: " << channel_data_type_str(destType) << "\n";
        //		cout << "Function is: " << get_type_convertor_string( get_type_convertor_function(sourceType,
        //destType, typeConversionChannelName) ) << "\n";
        m_typeConversionBlocks.push_back( detail::cma_type_conversion_block(
            sourcePosition, destPosition, typeConversionElements,
            get_channel_type_convertor_function( sourceType, destType, typeConversionChannelName ) ) );
        typeConversionElements = 0;
    }

    // Check whether we created an identity map
    if( m_typeConversionBlocks.size() == 1 &&
        m_defaultCopyBlocks.empty() && // Can't be an identity map if we need to insert default values from somewhere.
        m_sourceSize == m_destSize &&
        source.channel_count() >= dest.channel_count() && // Consider it an identity map if we are discarding some
        // m_typeConversionBlocks.front().count <= m_sourceSize && //Made this LEQ since padding can leave some junk
        // bytes that don't matter. BUG BUG BUG! This was making is_identity() true for cases where just the first
        // channel(s) were read.
        m_typeConversionBlocks.front().count == source.structure_size_without_padding() &&
        m_typeConversionBlocks.front().sourcePosition == 0 && m_typeConversionBlocks.front().destPosition == 0 &&
        m_typeConversionBlocks.front().convertor == &convert_memcpy ) {
        m_isIdentity = true;
    } else {
        m_isIdentity = false;
    }
}

void channel_map_adaptor::copy_structure( char* destData, const char* sourceData ) const {
    // Copy all the blocks with type conversion or memcpy
    for( vector<detail::cma_type_conversion_block>::const_iterator i = m_typeConversionBlocks.begin(),
                                                                   ie = m_typeConversionBlocks.end();
         i != ie; ++i )
        i->convertor( destData + i->destPosition, sourceData + i->sourcePosition, i->count );
}

void channel_map_adaptor::copy_structure_gaps_from_default( char* destData, const char* defaultData ) const {
    for( vector<detail::cma_type_conversion_block>::const_iterator it = m_defaultCopyBlocks.begin(),
                                                                   ite = m_defaultCopyBlocks.end();
         it != ite; ++it )
        it->convertor( destData + it->destPosition, defaultData + it->sourcePosition, it->count );
}

void channel_map_adaptor::dump( std::ostream& out ) const {
    out << "Particle Channel Map Adaptor:\n";
    for( vector<detail::cma_type_conversion_block>::const_iterator i = m_typeConversionBlocks.begin();
         i != m_typeConversionBlocks.end(); ++i )
        out << get_channel_type_convertor_debug_string( i->convertor ) << ": sourcePos=" << i->sourcePosition
            << ", destPos=" << i->destPosition << ", count=" << i->count << "\n";
    out << "Default Particle Copy:\n";
    for( vector<detail::cma_type_conversion_block>::const_iterator i = m_defaultCopyBlocks.begin();
         i != m_defaultCopyBlocks.end(); ++i )
        out << get_channel_type_convertor_debug_string( i->convertor ) << ": sourcePos=" << i->sourcePosition
            << ", destPos=" << i->destPosition << ", count=" << i->count << "\n";
}

} // namespace channels
} // namespace frantic
