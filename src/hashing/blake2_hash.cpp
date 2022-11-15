// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/hashing/blake2_hash.hpp>

namespace frantic {
namespace hashing {

blake2s_hash::blake2s_hash() {
    int error = blake2s_init( &m_state, BLAKE2S_OUTBYTES );
    if( error ) {
        throw std::runtime_error( "blake2s_hash::blake2s_hash Error initializing hash" );
    }
}

void blake2s_hash::update( const void* in, boost::uint64_t inlen ) {
    int error = blake2s_update( &m_state, reinterpret_cast<const boost::uint8_t*>( in ), inlen );
    if( error ) {
        throw std::runtime_error( "blake2s_hash::update Error updating hash" );
    }
}

void blake2s_hash::final_value( unsigned char* out ) {
    int error = blake2s_final( &m_state, out, BLAKE2S_OUTBYTES );
    if( error ) {
        throw std::runtime_error( "blake2s_hash::final_value Error finalizing hash" );
    }
}

blake2b_hash::blake2b_hash() {
    int error = blake2b_init( &m_state, BLAKE2B_OUTBYTES );
    if( error ) {
        throw std::runtime_error( "blake2b_hash::blake2b_hash Error initializing hash" );
    }
}

void blake2b_hash::update( const void* in, boost::uint64_t inlen ) {
    int error = blake2b_update( &m_state, reinterpret_cast<const boost::uint8_t*>( in ), inlen );
    if( error ) {
        throw std::runtime_error( "blake2b_hash::update Error updating hash" );
    }
}

void blake2b_hash::final_value( unsigned char* out ) {
    int error = blake2b_final( &m_state, out, BLAKE2B_OUTBYTES );
    if( error ) {
        throw std::runtime_error( "blake2b_hash::final_value Error finalizing hash" );
    }
}

blake2sp_hash::blake2sp_hash() {
    int error = blake2sp_init( &m_state, BLAKE2S_OUTBYTES );
    if( error ) {
        throw std::runtime_error( "blake2sp_hash::blake2sp_hash Error initializing hash" );
    }
}

void blake2sp_hash::update( const void* in, boost::uint64_t inlen ) {
    int error = blake2sp_update( &m_state, reinterpret_cast<const boost::uint8_t*>( in ), inlen );
    if( error ) {
        throw std::runtime_error( "blake2sp_hash::update Error updating hash" );
    }
}

void blake2sp_hash::final_value( unsigned char* out ) {
    int error = blake2sp_final( &m_state, out, BLAKE2S_OUTBYTES );
    if( error ) {
        throw std::runtime_error( "blake2sp_hash::final_value Error finalizing hash" );
    }
}

blake2bp_hash::blake2bp_hash() {
    int error = blake2bp_init( &m_state, BLAKE2B_OUTBYTES );
    if( error ) {
        throw std::runtime_error( "blake2bp_hash::blake2bp_hash Error initializing hash" );
    }
}

void blake2bp_hash::update( const void* in, boost::uint64_t inlen ) {
    int error = blake2bp_update( &m_state, reinterpret_cast<const boost::uint8_t*>( in ), inlen );
    if( error ) {
        throw std::runtime_error( "blake2bp_hash::update Error updating hash" );
    }
}

void blake2bp_hash::final_value( unsigned char* out ) {
    int error = blake2bp_final( &m_state, out, BLAKE2B_OUTBYTES );
    if( error ) {
        throw std::runtime_error( "blake2bp_hash::final_value Error finalizing hash" );
    }
}

} // namespace hashing
} // namespace frantic
