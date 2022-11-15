// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/particle_array_iterator.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <boost/static_assert.hpp>

#include <new>

#ifdef _MSC_VER
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#define ATOMIC_INCREMENT( val ) InterlockedIncrement( &val )
#define ATOMIC_DECREMENT( val ) InterlockedDecrement( &val )
#else
// This requires GCC 4.1 or so at least. If we need to build on older version, this should be switched to using the
// primitives in <atomic.h> or mutexes
#define ATOMIC_INCREMENT( val ) __sync_add_and_fetch( &val, 1 )
#define ATOMIC_DECREMENT( val ) __sync_add_and_fetch( &val, -1 )
#endif

namespace frantic {
namespace particles {

namespace detail {
class particle_istream_iterator_impl;
};

/**
 * Models an STL input iterator that wraps a particle_istream. The particle_istream is evaluated with repeated calls to
 * particle_istream::get_particles() storing the results in a temporary buffer. Most data is stored in a shared pointer
 * to implementation (which is stored at the beginning of the data buffer) since input iterators can be copy
 * constructed.
 */
class particle_istream_iterator
    : public std::iterator<std::input_iterator_tag, const char*, std::ptrdiff_t, void, void> {
  public:
    /**
     * The default constructor leaves the iterator in a state that will compare equal to an iterator that has advanced
     * past the end of the particle_istream.
     */
    particle_istream_iterator();

    /**
     * Constructs a new iterator that reads particles from the provided stream. The number of buffered particles can be
     * specified, where higher numbers improve parallelism but require more memory. \param particleStream The particle
     * stream to draw particles from while iterating. \param bufferSize The maximum number of particles to buffer at a
     * time.
     */
    particle_istream_iterator( frantic::particles::particle_istream_ptr particleStream,
                               std::size_t bufferSize = 10000 );

    ~particle_istream_iterator();

    /**
     * Copy constructor, required by STL input iterator spec. After copying, both iterators will affect the same stream.
     */
    particle_istream_iterator( const particle_istream_iterator& rhs );

    /**
     * Copy assignment, required by STL input iterator spec. After copying, both iterators will affect the same stream.
     */
    particle_istream_iterator& operator=( const particle_istream_iterator& rhs );

    /**
     * Advances the iterator to the next particle from the stream.
     */
    particle_istream_iterator& operator++();

    /**
     * Advances the iterator to the next particle from the stream. This is the postfix version.
     */
    particle_istream_iterator operator++( int );

    /**
     * Dereferences the iterator by returning a pointer to the current particle.
     * \return A pointer to the beginning of the first particle.
     */
    const char* operator*() const;

    /**
     * \return True if the two iterators are equal. They are equal if they are both past the end or default constructed,
     * or if they are copied from the same original iterator.
     */
    static bool is_equal( const particle_istream_iterator& lhs, const particle_istream_iterator& rhs );

  private:
    frantic::particles::particle_istream_ptr
        m_particleStream; // This could reasonably be moved inside m_pImpl, but I'm not confident about how it would
                          // behave wrt to alignment.

    detail::particle_istream_iterator_impl* m_pImpl;
};

/**
 * This class models a custom iterator scheme that wraps a particle_istream. The particle_istream is evaluated into a
 * buffer in batches at which point you can obtain a random accessor iterator to the beginning and end of the currently
 * buffered batch. Using the advance() member function will discard the current batch and refill the buffer with the
 * next batch of particles.
 */
class particle_istream_block_iterator {
  public:
    typedef frantic::particles::particle_array_iterator buffer_iterator;

  public:
    /**
     * The default constructor leaves the iterator in a state that will compare equal to an iterator that has advanced
     * past the end of the particle_istream.
     */
    particle_istream_block_iterator();

    /**
     * Constructs a new iterator that reads particles from the provided stream. The number of buffered particles can be
     * specified, where higher numbers improve parallelism but require more memory. \param particleStream The particle
     * stream to draw particles from while iterating. \param bufferSize The maximum number of particles to buffer at a
     * time.
     */
    particle_istream_block_iterator( frantic::particles::particle_istream_ptr particleStream,
                                     std::size_t bufferSize = 10000 );

    ~particle_istream_block_iterator();

    /**
     * Copy constructor.After copying, both iterators will affect the same stream.
     */
    particle_istream_block_iterator( const particle_istream_block_iterator& rhs );

    /**
     * Copy assignment. After copying, both iterators will affect the same stream.
     */
    particle_istream_block_iterator& operator=( const particle_istream_block_iterator& rhs );

    /**
     * Discards the current buffered particles and refills the buffer with the next batch of particles. This will
     * invalidate any existing iterators into the buffer. \return False if the stream has been exhausted and this
     * iterator cannot be advanced any further. True otherwise.
     */
    void advance();

    /**
     * Determines if the iterator is still valid, or if it has been exhausted.
     * \return True if the iterator is still valid, and supports calls to begin_buffer() and end_buffer(). False
     * indicates the internal particle_isteam is now empty.
     */
    bool valid() const;

    /**
     * Determines if the current buffer is empty.
     */
    bool empty() const;

    /**
     * Obtains an iterator to the start of the currently buffered particles.
     */
    buffer_iterator begin_buffer();

    /**
     * Obtains an iterator to the end of the currently buffered particles.
     */
    buffer_iterator end_buffer();

    /**
     * \return True if the two iterators are equal. They are equal if they are both past the end or default constructed,
     * or if they are copied from the same original iterator.
     */
    static bool is_equal( const particle_istream_block_iterator& lhs, const particle_istream_block_iterator& rhs );

  private:
    frantic::particles::particle_istream_ptr
        m_particleStream; // This could reasonably be moved inside m_pImpl, but I'm not confident about how it would
                          // behave wrt to alignment.

    detail::particle_istream_iterator_impl* m_pImpl;
};

namespace detail {
/**
 * This class is the internal shared implementation of the read buffer associated with this iterator. Since the iterator
 * is copyable, we need to share the buffer state between all iterator copies. We also need to do some book-keeping on
 * the buffer we are temporarily storing the particle data in, so its all placed inside this impl object.
 *
 * \note The constructors and destructors for this class are private since we can only create them on the heap via
 * impl::allocate() and destroy via. impl::release() when the internal ref count goes to 0.
 */
class particle_istream_iterator_impl {
  public:
    static particle_istream_iterator_impl* allocate( std::size_t stride, std::size_t bufferCount );

    void acquire();

    void release();

    bool refill_buffer( frantic::particles::streams::particle_istream& stream );

    bool advance();

    void* get();

    std::pair<void*, std::size_t> begin_buffer();

    std::pair<void*, std::size_t> end_buffer();

  private:
    particle_istream_iterator_impl( const particle_istream_iterator_impl& rhs );
    particle_istream_iterator_impl& operator=( const particle_istream_iterator_impl& rhs );

    particle_istream_iterator_impl( std::size_t stride, std::size_t bufferCount );

    ~particle_istream_iterator_impl();

  private:
    volatile long m_refCount;  // Number of iterators sharing this data, 4 bytes so the total structure is 32bytes.
    unsigned m_stride;         // Distance in bytes between particles, 4 bytes so the total structure is 32bytes.
    std::size_t m_bufferCount; // Max number of buffered particles, not necessarily the current number.
    void *m_it, *m_itEnd;      // Iterator for current and end of buffered particles

    enum { BUFFER_START_OFFSET = 32 };
};
} // namespace detail

inline particle_istream_iterator::particle_istream_iterator()
    : m_pImpl( NULL ) {}

inline particle_istream_iterator::particle_istream_iterator( frantic::particles::particle_istream_ptr particleStream,
                                                             std::size_t bufferCount )
    : m_particleStream( particleStream )
    , m_pImpl( NULL ) {
    std::size_t stride = particleStream->get_channel_map().structure_size();

    m_pImpl = detail::particle_istream_iterator_impl::allocate( stride, bufferCount );
    m_pImpl->acquire();

    if( !m_pImpl->refill_buffer( *m_particleStream ) ) {
        m_pImpl->release();
        m_pImpl = NULL;
    }
}

inline particle_istream_iterator::~particle_istream_iterator() {
    if( m_pImpl ) {
        m_pImpl->release();
        m_pImpl = NULL;
    }
}

inline particle_istream_iterator::particle_istream_iterator( const particle_istream_iterator& rhs )
    : m_particleStream( rhs.m_particleStream )
    , m_pImpl( rhs.m_pImpl ) {
    if( m_pImpl )
        m_pImpl->acquire();
}

inline particle_istream_iterator& particle_istream_iterator::operator=( const particle_istream_iterator& rhs ) {
    m_particleStream = rhs.m_particleStream;

    if( m_pImpl != rhs.m_pImpl ) { // Prevent self-assignment from causing problems.
        if( m_pImpl )
            m_pImpl->release();

        m_pImpl = rhs.m_pImpl;

        if( m_pImpl )
            m_pImpl->acquire();
    }
}

inline particle_istream_iterator& particle_istream_iterator::operator++() {
    if( !m_pImpl->advance() ) {
        if( !m_pImpl->refill_buffer( *m_particleStream ) ) {
            m_pImpl->release();
            m_pImpl = NULL;
        }
    }

    return *this;
}

inline particle_istream_iterator particle_istream_iterator::operator++( int ) {
    particle_istream_iterator result = *this;

    ++( *this );

    return result;
}

inline const char* particle_istream_iterator::operator*() const { return reinterpret_cast<char*>( m_pImpl->get() ); }

inline bool particle_istream_iterator::is_equal( const particle_istream_iterator& lhs,
                                                 const particle_istream_iterator& rhs ) {
    return ( lhs.m_pImpl == rhs.m_pImpl );
}

inline bool operator==( const particle_istream_iterator& lhs, const particle_istream_iterator& rhs ) {
    return particle_istream_iterator::is_equal( lhs, rhs );
}

inline bool operator!=( const particle_istream_iterator& lhs, const particle_istream_iterator& rhs ) {
    return !particle_istream_iterator::is_equal( lhs, rhs );
}

inline particle_istream_block_iterator::particle_istream_block_iterator()
    : m_pImpl( NULL ) {}

inline particle_istream_block_iterator::particle_istream_block_iterator(
    frantic::particles::particle_istream_ptr particleStream, std::size_t bufferCount )
    : m_particleStream( particleStream )
    , m_pImpl( NULL ) {
    std::size_t stride = particleStream->get_channel_map().structure_size();

    m_pImpl = detail::particle_istream_iterator_impl::allocate( stride, bufferCount );
    m_pImpl->acquire();

    if( !m_pImpl->refill_buffer( *m_particleStream ) ) {
        m_pImpl->release();
        m_pImpl = NULL;
    }
}

inline particle_istream_block_iterator::~particle_istream_block_iterator() {
    if( m_pImpl ) {
        m_pImpl->release();
        m_pImpl = NULL;
    }
}

inline particle_istream_block_iterator::particle_istream_block_iterator( const particle_istream_block_iterator& rhs )
    : m_particleStream( rhs.m_particleStream )
    , m_pImpl( rhs.m_pImpl ) {
    if( m_pImpl )
        m_pImpl->acquire();
}

inline particle_istream_block_iterator&
particle_istream_block_iterator::operator=( const particle_istream_block_iterator& rhs ) {
    m_particleStream = rhs.m_particleStream;

    if( m_pImpl != rhs.m_pImpl ) {
        if( m_pImpl )
            m_pImpl->release();

        m_pImpl = rhs.m_pImpl;

        if( m_pImpl )
            m_pImpl->acquire();
    }
}

inline void particle_istream_block_iterator::advance() {
    if( m_pImpl && !m_pImpl->refill_buffer( *m_particleStream ) ) {
        m_pImpl->release();
        m_pImpl = NULL;
    }
}

inline bool particle_istream_block_iterator::valid() const { return m_pImpl != NULL; }

inline bool particle_istream_block_iterator::empty() const {
    return !m_pImpl || m_pImpl->begin_buffer().first == m_pImpl->end_buffer().first;
}

inline particle_istream_block_iterator::buffer_iterator particle_istream_block_iterator::begin_buffer() {
    std::pair<void*, std::size_t> beginImpl = m_pImpl->begin_buffer();

    return particle_array_iterator( reinterpret_cast<char*>( beginImpl.first ), beginImpl.second );
}

inline particle_istream_block_iterator::buffer_iterator particle_istream_block_iterator::end_buffer() {
    std::pair<void*, std::size_t> endImpl = m_pImpl->end_buffer();

    return particle_array_iterator( reinterpret_cast<char*>( endImpl.first ), endImpl.second );
}

inline bool particle_istream_block_iterator::is_equal( const particle_istream_block_iterator& lhs,
                                                       const particle_istream_block_iterator& rhs ) {
    return lhs.m_pImpl == rhs.m_pImpl;
}

namespace detail {
inline particle_istream_iterator_impl::particle_istream_iterator_impl( std::size_t stride, std::size_t bufferCount )
    : m_refCount( 0 )
    , m_stride( static_cast<unsigned>( stride ) )
    , m_bufferCount( bufferCount ) {
    m_it = m_itEnd = NULL;
}

inline particle_istream_iterator_impl::~particle_istream_iterator_impl() {}

inline particle_istream_iterator_impl* particle_istream_iterator_impl::allocate( std::size_t stride,
                                                                                 std::size_t bufferCount ) {
    // Since sizeof(impl) is different for 32/64bit builds I use 32bytes as the offset and ensure this via a
    // static_assert
    BOOST_STATIC_ASSERT_MSG( sizeof( particle_istream_iterator_impl ) <= BUFFER_START_OFFSET,
                             "sizeof(particle_istream_iterator_impl) has changed and violated an assumption" );

    if( bufferCount < 1u )
        bufferCount = 1u; // Force at least one particle in the buffer.

    // This part is a little funky, but essentially we are allocating raw memory (via operator new, which is kinda like
    // malloc) to store both the impl object AND a data buffer for storing particles. The data buffer begins at pImpl +
    // sizeof(impl). Although alignment is not important for particle data, in this case sizeof(impl) == 32 so the
    // particle buffer begins with appropriate alignment for all numeric types.
    void* pBuffer = operator new( BUFFER_START_OFFSET + stride * bufferCount );

    // Use placement new to construct the impl object at the start of the memory we allocated.
    return new( pBuffer ) particle_istream_iterator_impl( stride, bufferCount );
}

inline void particle_istream_iterator_impl::acquire() { ATOMIC_INCREMENT( m_refCount ); }

inline void particle_istream_iterator_impl::release() {
    if( ATOMIC_DECREMENT( m_refCount ) == 0 ) {
        this->~particle_istream_iterator_impl();

        operator delete( this );
    }
}

inline bool particle_istream_iterator_impl::refill_buffer( frantic::particles::streams::particle_istream& stream ) {
    bool eos;
    std::size_t numParticles;

    char* buffer = reinterpret_cast<char*>( this ) +
                   BUFFER_START_OFFSET; // This is ok since we allocated the memory ourselves via operator new(), so our
                                        // data buffer begins after this + BUFFER_START_OFFSET

    do {
        numParticles = m_bufferCount;
        eos = !stream.get_particles( buffer, numParticles );
    } while( !eos && numParticles == 0 ); // Loop while we didn't get any particles, and haven't hit EOS. Some streams
                                          // could delete all particles that we attempted to read so we end up with 0.

    m_it = buffer;
    m_itEnd = buffer + m_stride * numParticles;

    return !eos || numParticles > 0;
}

inline bool particle_istream_iterator_impl::advance() {
    m_it = reinterpret_cast<char*>( m_it ) + m_stride;
    return ( m_it != m_itEnd );
}

inline void* particle_istream_iterator_impl::get() { return m_it; }

inline std::pair<void*, std::size_t> particle_istream_iterator_impl::begin_buffer() {
    return std::make_pair( m_it, m_stride );
}

inline std::pair<void*, std::size_t> particle_istream_iterator_impl::end_buffer() {
    return std::make_pair( m_itEnd, m_stride );
}
} // namespace detail

} // namespace particles
} // namespace frantic
