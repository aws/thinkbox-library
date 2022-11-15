// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

// These are the functions from Bob Jenkins' public domain hash program
extern "C" {

/*
--------------------------------------------------------------------
 This works on all machines.  To be useful, it requires
 -- that the key be an array of uint32_t's, and
 -- that the length be the number of uint32_t's in the key

 The function hashword() is identical to hashlittle() on little-endian
 machines, and identical to hashbig() on big-endian machines,
 except that the length has to be measured in uint32_ts rather than in
 bytes.  hashlittle() is more complicated than hashword() only because
 hashlittle() has to dance around fitting the key bytes into registers.
--------------------------------------------------------------------
*/
boost::uint32_t hashword( const boost::uint32_t* k,  /* the key, an array of uint32_t values */
                          size_t length,             /* the length of the key, in uint32_ts */
                          boost::uint32_t initval ); /* the previous hash, or an arbitrary value */

void hashword2( const boost::uint32_t* k, /* the key, an array of uint32_t values */
                size_t length,            /* the length of the key, in uint32_ts */
                boost::uint32_t* pc,      /* IN: seed OUT: primary hash value */
                boost::uint32_t* pb );    /* IN: more seed OUT: secondary hash value */

boost::uint32_t hashlittle( const void* key, size_t length, boost::uint32_t initval );

void hashlittle2( const void* key,       /* the key to hash */
                  size_t length,         /* length of the key */
                  boost::uint32_t* pc,   /* IN: primary initval, OUT: primary hash */
                  boost::uint32_t* pb ); /* IN: secondary initval, OUT: secondary hash */

} // extern "C"

// Provide the same hashword function for 64 bit
inline boost::uint64_t hashword64( const boost::uint32_t* key, size_t length, boost::uint64_t initval ) {
    boost::uint32_t pb = ( boost::uint32_t )( initval >> 32 ), pc = (boost::uint32_t)initval;
    hashword2( key, length, &pc, &pb );
    return pc + ( ( (boost::uint64_t)pb ) << 32 );
}

// Provide the same hashlittle function for 64 bit
inline boost::uint64_t hashlittle64( const void* key, size_t length, boost::uint64_t initval ) {
    boost::uint32_t pb = ( boost::uint32_t )( initval >> 32 ), pc = (boost::uint32_t)initval;
    hashlittle2( key, length, &pc, &pb );
    return pc + ( ( (boost::uint64_t)pb ) << 32 );
}
