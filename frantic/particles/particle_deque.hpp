// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/particles/particle_deque_iterator.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

namespace frantic {
namespace particles {

// Returns y = 2^x such that y <= val
inline unsigned floor_power_of_two( unsigned val ) {
    val = val | ( val >> 1 );
    val = val | ( val >> 2 );
    val = val | ( val >> 4 );
    val = val | ( val >> 8 );
    val = val | ( val >> 16 );

    return val - ( val >> 1 );
}

// Returns y = 2^x such that y <= val
/*inline boost::uint32_t floor_power_of_two(boost::uint32_t val)
{
    val = val | (val >> 1);
    val = val | (val >> 2);
    val = val | (val >> 4);
    val = val | (val >> 8);
    val = val | (val >> 16);

    return val - (val >> 1);
}*/

// Returns y = 2^x such that y <= val
inline boost::uint64_t floor_power_of_two( boost::uint64_t val ) {
    val = val | ( val >> 1 );
    val = val | ( val >> 2 );
    val = val | ( val >> 4 );
    val = val | ( val >> 8 );
    val = val | ( val >> 16 );
    val = val | ( val >> 32 );

    return val - ( val >> 1 );
}

class particle_deque {
  private:
    static const int MAX_BLOCK_SIZE = 64 * ( 1 << 10 ); // Maximum block size in bytes (64Kb)
    static const int INDEX_PADDING = 4; // Number of empty blocks to alloc in index array on either end when re-allocing

    typedef char byte;

    std::size_t m_indexBlockCount; // Size of index array
    std::size_t m_indexBlockFront; // Offset into index array for first allocated block
    std::size_t m_indexBlockBack;  // Offset into index array for 1 past the last allocated block
    byte** m_indexBlock;           // Array of ptrs to BLOCK_SIZE particles arrays

    std::size_t m_frontOffset;  // Number of unused particles in front block
    std::size_t m_numParticles; // Number of particles in data structure
    std::size_t m_numPerBlock;  // Number of particles per block

    std::size_t m_particleSize; // Cached locally for good cache consistency
    channels::channel_map m_particleChannelMap;

    inline void grow_front();
    inline void grow_back();

  public:
    typedef particle_deque_iterator iterator;
    typedef particle_deque_const_iterator const_iterator;

    inline particle_deque();
    inline particle_deque( const channels::channel_map& pcm )
        : m_indexBlock( NULL ) {
        reset( pcm );
    }
    inline ~particle_deque();

    inline void reset( const channels::channel_map& pcm );
    inline void clear();

    inline const std::size_t size() const { return m_numParticles; }
    // TODO: Is this correct?
    inline const std::size_t size_in_memory() const { return m_indexBlockCount * m_numPerBlock; }
    inline const bool empty() const { return m_indexBlock == NULL; }

    inline const channels::channel_map& get_channel_map() const { return m_particleChannelMap; }
    inline const std::size_t get_particle_size() const { return m_particleSize; }

    inline void resize( std::size_t size );
    inline void reserve( std::size_t size );

    inline void push_back( const void* buffer );
    inline void push_front( const void* buffer );

    inline void pop_back();
    // TODO: pop_front(), and pop_back() to make a particle queue/stack ... ?

    inline const_iterator begin() const;
    inline const_iterator end() const;
    inline iterator begin();
    inline iterator end();

    inline const char* at( const std::size_t i ) const;
    inline const char* operator[]( const std::size_t i ) const;
    inline const char* back() const;
    inline const char* front() const;

    inline char* at( const std::size_t i );
    inline char* operator[]( const std::size_t i );
    inline char* back();
    inline char* front();

    void insert_particles( boost::shared_ptr<streams::particle_istream> pin );
    void insert_particles( boost::shared_ptr<streams::particle_istream> pin,
                           frantic::logging::progress_logger& progress );

    template <class ForwardIterator>
    void insert_particles( const frantic::channels::channel_map& pcm, ForwardIterator start, ForwardIterator end );
};

particle_deque::particle_deque()
    : m_indexBlock( NULL ) {
    frantic::channels::channel_map pcm;
    pcm.define_channel<frantic::graphics::vector3f>( _T("Position") );
    pcm.end_channel_definition();

    reset( pcm );
}

particle_deque::~particle_deque() { clear(); }

void particle_deque::reset( const channels::channel_map& pcm ) {
    clear();

    m_particleChannelMap = pcm;
    m_particleSize = pcm.structure_size();
    m_numPerBlock = floor_power_of_two( MAX_BLOCK_SIZE / pcm.structure_size() );

    if( m_numPerBlock < 1 )
        throw std::runtime_error( "particle_deque::reset() - the supplied channel map describes a structure too large "
                                  "for the particle_deque" );
}

void particle_deque::clear() {
    if( m_indexBlock ) {
        for( std::size_t i = m_indexBlockFront; i < m_indexBlockBack; ++i )
            free( m_indexBlock[i] );
        free( m_indexBlock );
    }

    m_indexBlock = NULL;
    m_indexBlockCount = m_indexBlockFront = m_indexBlockBack = m_frontOffset = m_numParticles = 0;
}

void particle_deque::grow_back() {
    if( m_indexBlockBack == m_indexBlockCount ) {
        m_indexBlockCount += INDEX_PADDING;
        m_indexBlock = static_cast<byte**>( realloc( m_indexBlock, sizeof( byte* ) * m_indexBlockCount ) );
        if( !m_indexBlock )
            throw std::runtime_error( "particle_deque::grow_back() - failed to reallocate the index block to size: " +
                                      boost::lexical_cast<std::string>( sizeof( byte* ) * m_indexBlockCount ) );
    }

    byte* newBlock = static_cast<byte*>( malloc( m_numPerBlock * m_particleSize ) );
    if( !newBlock )
        throw std::runtime_error( "particle_deque::grow_back() - failed to allocate a new chunk of size: " +
                                  boost::lexical_cast<std::string>( m_numPerBlock * m_particleSize ) );

    m_indexBlock[m_indexBlockBack++] = newBlock;
}

void particle_deque::grow_front() {
    if( m_indexBlockFront == 0 ) {
        byte** newBlock =
            static_cast<byte**>( realloc( m_indexBlock, sizeof( byte* ) * ( m_indexBlockCount + INDEX_PADDING ) ) );
        if( !newBlock )
            throw std::runtime_error(
                "particle_deque::grow_front() - failed to reallocate the index block to size: " +
                boost::lexical_cast<std::string>( sizeof( byte* ) * ( m_indexBlockCount + INDEX_PADDING ) ) );

        memmove( newBlock + INDEX_PADDING, newBlock, sizeof( byte* ) * m_indexBlockCount );

        m_indexBlock = newBlock;
        m_indexBlockFront = INDEX_PADDING; // Shift forward to compensate for new items
        m_indexBlockBack += INDEX_PADDING;
        m_indexBlockCount += INDEX_PADDING;
    }

    byte* newBlock = static_cast<byte*>( malloc( m_numPerBlock * m_particleSize ) );
    if( !newBlock )
        throw std::runtime_error( "particle_deque::grow_front() - failed to allocate a new chunk of size: " +
                                  boost::lexical_cast<std::string>( m_numPerBlock * m_particleSize ) );

    m_indexBlock[--m_indexBlockFront] = newBlock;
    m_frontOffset = m_numPerBlock;
}

void particle_deque::resize( std::size_t size ) {
    if( size == 0 ) {
        m_indexBlockBack = m_indexBlockFront;
        m_frontOffset = 0;
        m_numParticles = 0;
    } else {
        // We need the ability to index 'size' items overall after this method.
        // We have (m_numPerBlock - m_frontOffset) in the front block, so we need
        //'size - (m_numPerBlock - m_frontOffset)' more items. Add m_numPerBlock - 1 in order to round up.
        std::size_t numNeededBlocks = 1 + ( size + m_frontOffset - 1 ) / m_numPerBlock;
        while( numNeededBlocks > ( m_indexBlockBack - m_indexBlockFront ) )
            grow_back();
        m_indexBlockBack = m_indexBlockFront + numNeededBlocks;
        m_numParticles = size;
    } //( size != 0 ){
}

void particle_deque::reserve( std::size_t size ) {
    // There are no constructors in this container, so this is effectively the same result as
    // a normal reserve().
    std::size_t oldSize = m_numParticles;
    if( size > oldSize ) {
        resize( size );
        resize( oldSize );
    }
}

void particle_deque::push_front( const void* buffer ) {
    if( m_frontOffset == 0 )
        grow_front();

    --m_frontOffset;
    ++m_numParticles;
    memcpy( m_indexBlock[m_indexBlockFront] + ( m_frontOffset * m_particleSize ), buffer, m_particleSize );
}

void particle_deque::push_back( const void* buffer ) {
    std::size_t index = ( m_numParticles + m_frontOffset ) / m_numPerBlock + m_indexBlockFront;
    std::size_t offset = ( m_numParticles + m_frontOffset ) & ( m_numPerBlock - 1 );

    assert( index <= m_indexBlockBack );

    if( index == m_indexBlockBack )
        grow_back();

    ++m_numParticles;
    memcpy( m_indexBlock[index] + ( offset * m_particleSize ), buffer, m_particleSize );
}

void particle_deque::pop_back() {
    --m_numParticles;

    // NOTE: This may not leave the deque in the same state as clear() when m_numParticles == 0 after a pop_back().
    //       The reason for that is that after calling push_front(), the first particle may not be at block index 0,
    //       particle index 0.
}

particle_deque::iterator particle_deque::begin() { return iterator( *this, 0 ); }

particle_deque::iterator particle_deque::end() { return iterator( *this, m_numParticles ); }

particle_deque::const_iterator particle_deque::begin() const { return const_iterator( *this, 0 ); }

particle_deque::const_iterator particle_deque::end() const { return const_iterator( *this, m_numParticles ); }

const char* particle_deque::operator[]( const std::size_t i ) const {
    std::size_t index = ( i + m_frontOffset ) / m_numPerBlock + m_indexBlockFront;
    std::size_t offset = ( i + m_frontOffset ) & ( m_numPerBlock - 1 );

    return m_indexBlock[index] + ( offset * m_particleSize );
}

char* particle_deque::operator[]( const std::size_t i ) {
    std::size_t index = ( i + m_frontOffset ) / m_numPerBlock + m_indexBlockFront;
    std::size_t offset = ( i + m_frontOffset ) & ( m_numPerBlock - 1 );

    return m_indexBlock[index] + ( offset * m_particleSize );
}

const char* particle_deque::at( const std::size_t i ) const { return ( *this )[i]; }
char* particle_deque::at( const std::size_t i ) { return ( *this )[i]; }

const char* particle_deque::back() const { return at( m_numParticles - 1 ); }
char* particle_deque::back() { return at( m_numParticles - 1 ); }

const char* particle_deque::front() const { return at( 0 ); }
char* particle_deque::front() { return at( 0 ); }

std::size_t particle_deque_iterator::structure_size() { return m_container->get_channel_map().structure_size(); }
particle_deque_iterator::value_type particle_deque_iterator::operator*() { return m_container->at( m_offset ); }
particle_deque_iterator::value_type particle_deque_iterator::operator[]( std::size_t i ) {
    return m_container->at( m_offset + i );
}

std::size_t particle_deque_const_iterator::structure_size() { return m_container->get_channel_map().structure_size(); }
particle_deque_const_iterator::value_type particle_deque_const_iterator::operator*() {
    return m_container->at( m_offset );
}
particle_deque_const_iterator::value_type particle_deque_const_iterator::operator[]( std::size_t i ) {
    return m_container->at( m_offset + i );
}

inline void particle_deque::insert_particles( boost::shared_ptr<streams::particle_istream> pin ) {
    frantic::logging::null_progress_logger nullProgress;
    insert_particles( pin, nullProgress );
}

inline void particle_deque::insert_particles( boost::shared_ptr<streams::particle_istream> pin,
                                              frantic::logging::progress_logger& progress ) {
#pragma warning( push )
#pragma warning( disable : 4127 )

    // Enforce the stream to have the same map as this container.
    pin->set_channel_map( m_particleChannelMap );

    std::size_t index = ( m_numParticles + m_frontOffset ) / m_numPerBlock + m_indexBlockFront;
    std::size_t offset = ( m_numParticles + m_frontOffset ) & ( m_numPerBlock - 1 );

    boost::int64_t progressIndex;
    boost::int64_t progressCount = pin->particle_progress_count();

    // Loop through, requesting chunks of particles to fill each block of the deque.
    do {
        if( index == m_indexBlockBack )
            grow_back();

        std::size_t expectedCount = m_numPerBlock - offset;
        if( !pin->get_particles( m_indexBlock[index] + ( offset * m_particleSize ), expectedCount ) ) {
            m_numParticles += expectedCount;
            progress.update_progress( m_numParticles, m_numParticles );
            return;
        }

        m_numParticles += expectedCount;
        index = ( m_numParticles + m_frontOffset ) / m_numPerBlock + m_indexBlockFront;
        offset = ( m_numParticles + m_frontOffset ) & ( m_numPerBlock - 1 );

        if( progressCount >= 0 ) {
            progressIndex = pin->particle_progress_index();
            progress.update_progress( progressIndex, progressCount );
        }
    } while( 1 );

#pragma warning( pop )
}

template <class ForwardIterator>
void particle_deque::insert_particles( const frantic::channels::channel_map& pcm, ForwardIterator iter,
                                       ForwardIterator iterEnd ) {
    // make a converter to go to from the incoming channel map to the array's channel map
    frantic::channels::channel_map_adaptor cma( m_particleChannelMap, pcm );
    std::size_t destinationSize = m_particleChannelMap.structure_size();

    if( cma.is_identity() ) {
        for( ; iter != iterEnd; ++iter )
            push_back( *iter );
    } else {
        char* temp = (char*)alloca( m_particleChannelMap.structure_size() );

        m_particleChannelMap.construct_structure( temp );

        for( ; iter != iterEnd; ++iter ) {
            cma.copy_structure( temp, *iter );
            push_back( temp );
        }
    }
}

} // namespace particles
} // namespace frantic
