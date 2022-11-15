// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/channels/channel_operation_compiler.hpp>

namespace frantic {
namespace channels {

channel_operation_compiler::channel_operation_compiler( const frantic::channels::channel_map& pcm,
                                                        const frantic::channels::channel_map& _nativePcm )
    : nativePcm( _nativePcm )
    , internalPcm( pcm )
    , nextOutputOffset( 0 )
    , hasIndexChannel( false ) {}

channel_operation_compiler::~channel_operation_compiler() {
    for( std::size_t i = 0; i < codeSegs.size(); ++i )
        free( codeSegs[i].codeData );
    for( std::size_t i = 0; i < scopedObjs.size(); ++i )
        scopedObjs[i].second( scopedObjs[i].first );
}

void channel_operation_compiler::reset( const frantic::channels::channel_map& pcm,
                                        const frantic::channels::channel_map& _nativePcm ) {
    for( std::size_t i = 0; i < codeSegs.size(); ++i )
        free( codeSegs[i].codeData );
    codeSegs.clear();

    for( std::size_t i = 0; i < scopedObjs.size(); ++i )
        scopedObjs[i].second( scopedObjs[i].first );
    scopedObjs.clear();

    indexAcc.reset();
    tempOffsets.clear();

    hasIndexChannel = false;
    nextOutputOffset = 0;
    internalPcm = pcm;
    nativePcm = _nativePcm;
}

void channel_operation_compiler::set_has_index_channel( bool hasChannel ) {
    hasIndexChannel = hasChannel;
    if( hasChannel )
        indexAcc = internalPcm.get_cvt_accessor<int>( _T("Index") );
}

// void channel_operation_compiler::debug_print( std::ostream& stream ) const {
//	std::map<int,temporary_result>::const_iterator it = tempOffsets.begin(), itEnd = tempOffsets.end();
//	for( ; it != itEnd; ++it )
//		stream << it->first << " : " << frantic::strings::to_string(
// frantic::channels::channel_data_type_str(it->second.arity, it->second.type) ) << " @" << it->second.offset <<
// std::endl; 	stream << "Internal map:\n" << internalPcm << std::endl; 	stream << "Native map:\n" << nativePcm
// << std::endl;
// }

frantic::channels::channel_map& channel_operation_compiler::get_channel_map() { return internalPcm; }

frantic::channels::channel_map& channel_operation_compiler::get_native_channel_map() { return nativePcm; }

const frantic::channels::channel_map& channel_operation_compiler::get_channel_map() const { return internalPcm; }

const frantic::channels::channel_map& channel_operation_compiler::get_native_channel_map() const { return nativePcm; }

temporary_result* channel_operation_compiler::get_node_results( int nodeId ) {
    std::map<int, temporary_result>::iterator it = tempOffsets.find( nodeId );
    if( it != tempOffsets.end() )
        return &it->second;
    return NULL;
}

void channel_operation_compiler::copy_node_results( int destNodeId, int sourceNodeId ) {
    tempOffsets[destNodeId] = tempOffsets[sourceNodeId];
}

temporary_result* channel_operation_compiler::allocate_temporary( int nodeId, data_type_t type, int arity ) {
    int primitiveSize = arity * (int)frantic::channels::sizeof_channel_data_type( type );
    int lastOutputOffset = nextOutputOffset;
    nextOutputOffset += primitiveSize;

    std::map<int, temporary_result>::iterator it = tempOffsets.find( nodeId );
    if( it != tempOffsets.end() )
        throw std::runtime_error( "Recompiled node: " + boost::lexical_cast<std::string>( nodeId ) );

    temporary_result result;
    result.offset = lastOutputOffset;
    result.arity = arity;
    result.type = type;

    it = tempOffsets.insert( std::pair<int, temporary_result>( nodeId, result ) ).first;
    return &it->second;
}

void channel_operation_compiler::append_code_segment( const detail::code_segment& segment ) {
    if( segment.codePtr == NULL )
        throw std::runtime_error( "append_code_segment() - Recieved a NULL fn ptr" );
    codeSegs.push_back( segment );
}

void channel_operation_compiler::eval( char* particle, std::size_t index ) const {
    char* tempStack = (char*)alloca( nextOutputOffset );

    if( hasIndexChannel )
        indexAcc.set( particle, (int)index );

    for( std::size_t i = 0; i < codeSegs.size(); ++i )
        codeSegs[i].codePtr( particle, tempStack, codeSegs[i].codeData );
}

} // namespace channels
} // namespace frantic
