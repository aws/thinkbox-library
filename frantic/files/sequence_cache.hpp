// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/caches/lru_cache.hpp>
#include <frantic/files/background_serializer.hpp>
#include <frantic/files/filename_sequence.hpp>
#include <frantic/strings/tstring.hpp>

#include <set>

namespace frantic {
namespace files {

/**
 * \brief Implements a two level cache heirarchy for a data sequence (ex. PRT particle data, XMesh mesh data, etc.)
 * consisting of an least-recently-used (lru) memory cache in conjunction with a file sequence stored on disk (or
 * wherever since the (de)serialization is configured via template).
 *
 * \details Each stored value is associated with a floating point time which identifies it. Values on disk are stored as
 * multiple files in a single directory using a filename_sequence object to create the filenames. A limit can be imposed
 * on the maximum amount of data kept in memory. When the limit is reached, the least-recently-used values will be
 * asynchronously transferred to disk storage to make space.
 *
 * The two main use cases for this class are:
 *  1. A lazy write-through cache for an object which produces a sequence of data. The producer object inserts data into
 * the cache (and possibly reads it back too) and the data is serialized in the background as the cache is filled.
 *  2. A lazy read-back cache for an object which consumes a disk sequence of data. The consumer uses the cache to
 * lazily load data from the disk sequence as it is requested. The least used data is automatically cleared from memory
 * as the cache fills up.
 *
 * \tparam ValueType The type that contains data in the cache and can be serialized.
 * \tparam SerializerType The type that can serialize a 'ValueType' instance. Must implement:
 *                         void serialize( const frantic::tstring& path, const ValueType& val ) const throw(); // NOTE:
 * It must not throw! ValueType deserialize( const frantic::tstring& path ) const; \tparam SizeEstimatorType A functor
 * type that can calculate the size (in generic cache units) of a 'ValueType' instance. Must implement: std::size_t
 * operator()( const ValueType& val ) const;
 *
 * \note This object works with value semantics so if ValueType is an expensive object to copy, consider using a
 * shared_ptr<> instead.
 */
template <class ValueType, class SerializerType, class SizeEstimatorType>
class sequence_cache {
  public:
    typedef typename frantic::caches::lru_cache<double, ValueType>::const_iterator const_memory_iterator;
    typedef typename std::set<double>::const_iterator const_key_iterator;

  public:
    /**
     * Constructor
     * \param cacheLimit The maximum size (in generic cache units) that can be cached in memory before the least
     * recently used item is flushed to disk. \param theSerializer The object responsible for serializing and
     * deserializing a ValueType instance to/from disk. \param theSizeEstimator The object responsible for calculating
     * the size (in generic cache units) of a ValueType instance.
     */
    explicit sequence_cache( std::size_t cacheLimit = 10000000ul,
                             const SerializerType& theSerializer = SerializerType(),
                             const SizeEstimatorType& theSizeEstimator = SizeEstimatorType() );

    enum sync_option { keep_existing, drop_existing, synchronize };

    /**
     * Clears the existing cache and sets the sequenece path pattern that is used when serializing data to disk.
     * \param diskSequencePath The sequence path pattern for storing serialized data. Ex.
     * "C:\Users\Darcy\Test\test_####.prt" If this is empty, the disk cache will be disabled. \param syncOption How the
     * cache contents are handled when assigning the new path. If 'keep_existing' the items held in memory are kept.
     * This option can be weird if some items were serialized to disk since those will disappear. If 'drop_existing' the
     * cache is completely cleared (as if by clear). If 'synchronize' the cache is completely cleared and then
     * repopulated using the new sequence path.
     */
    void set_disk_path( const frantic::tstring& diskSequencePath, sync_option syncOption );

    /**
     * Sets the maximum size (in generic cache units) that can be cached in memory before the least recently used item
     * is flushed to disk. This change may cause data to be flushed from the cache. \param cacheLimit The new limit on
     * the size of the memory cache.
     */
    void set_capacity( std::size_t cacheLimit );

    /**
     * Sets the limit on how many cache units can be waiting for serialization before we wait for the I/O system to
     * catch up. \param pendingQueueLimit The new limit on the size of the items waiting to be serialized.
     */
    void set_pending_capacity( std::size_t pendingQueueLimit );

    /**
     * Sets the number of threads that the contained background_serializer will use to process items.
     * \param numThreads The number of threads to create as background workers.
     */
    void set_num_serializer_threads( std::size_t numThreads );

    /**
     * Adds a new entry to the cache, or replaces an existing entry.
     * \param frame The time (in frames) that theValue is associated with.
     * \param theValue The new value to store in the cache for the given time.
     */
    void insert( double frame, const ValueType& theValue );

    /**
     * Discards the memory and disk caches.
     */
    void clear();

    /**
     * Serializes all items to disk, optionally clearing the memory cache.
     */
    void flush( bool emptyMemoryCache );

    /**
     * Serializes all items to disk asynchronously. The callback function is invoked when serializing is complete.
     */
    template <class Callable>
    void flush_async( const Callable& callback );

    /**
     * Cancels a previously started asynchrous flush. The function returns immediately, but the cancellation may not
     * complete until active serializing items are completed.
     */
    void cancel_flush();

    /**
     * \return A bidirectional iterator to the start of the collection of cache sample times.
     */
    const_key_iterator key_begin() const;

    /**
     * \return A bidirectional iterator past the end of the collection of cache sample times.
     */
    const_key_iterator key_end() const;

    /**
     * Returns an iterator that visits each item currently held in memory. The result is invalidated after calls to
     * insert() or find(). \return An bidirectional iterator to the first memory cached item.
     */
    const_memory_iterator memory_begin() const;

    /**
     * Returns an iterator that visits each item currently held in memory. The result is invalidated after calls to
     * insert() or find(). \return An bidirectional iterator one past the last memory cached item.
     */
    const_memory_iterator memory_end() const;

    /**
     * Returns a copy of the filename sequence pattern used by the disk portion of the cache. Will be an empty string if
     * disk cacheing is disabled.
     */
    frantic::tstring get_disk_path_pattern() const;

    /**
     * \return True if the cache is currently empty
     */
    bool empty() const;

    /**
     * \return The min & max key values in the sequence. Returns NaN,NaN if the sequence is empty.
     */
    std::pair<double, double> get_key_range() const;

    /**
     * \return The current amount of memory in use by the memory cache, in generic cache units.
     */
    std::size_t get_usage() const;

    /**
     * \return The maximum capacity of the memory cache, in generic cache units.
     */
    std::size_t get_capacity() const;

    /**
     * \return The current amount of cache units awaiting serialization.
     */
    std::size_t get_pending_usage() const;

    /**
     * \return The maximum number of cache units that can be awaiting serialization.
     */
    std::size_t get_pending_capacity() const;

    /**
     * If 'frame' is a valid key for an item in the cache (in memory or on disk) then this will return the associated
     * value. \param frame The time (in frames) that we want a cache entry for. \return The value associated with the
     * specified frame, or a default constructed value if there is no associated entry.
     */
    ValueType find( double frame ) const;

    /**
     * Returns the key that is closest in absolute value to the specified time.
     * \param frame The reference time from which we want the closest cache entry key.
     * \return The time of the closest cache entry. Returns NaN if the cache is empty.
     */
    double find_nearest_key( double frame ) const;

    /**
     * Returns the two closest cache keys that are on either side of the specified time. Ex. find_bracketing_keys( 1.3 )
     * return <1, 2> if the cache contained [1, 2, 3, 4, 5]. \param frame The time to find the closest bracketing keys.
     * \return The two closest keys on either side of the specified time. Can return negative or positive infinity if
     * there is not a key on the left or right respectively.
     */
    std::pair<double, double> find_bracketing_keys( double frame ) const;

    /**
     * Blocks the calling thread until all items that are pending serialization have been processed.
     */
    void wait_for_pending() const;

    /**
     * \return Returns the serializer this cache is using to write objects to disk.
     */
    background_serializer<frantic::tstring, ValueType, SerializerType>& get_serializer();

  private:
    class deallocate_callback;

    class clear_callback;

    typedef frantic::caches::lru_cache<double, ValueType> cache_type;

    typedef background_serializer<frantic::tstring, ValueType, SerializerType> serializer_type;

  private:
    /**
     * Called when the memory cache is about to remove an entry. This is where the entry should be serialized to disk if
     * it hasn't been already. \param frame The time associated with 'theVal' \param theVal The data in the cache.
     * \param theReason Indicates why the cache entry is being removed. \return The size of the cache entry being
     * deallocated.
     */
    std::size_t on_cache_entry_flushed( double frame, const ValueType& theVal,
                                        typename cache_type::deallocator_reason theReason ) const;

  private:
    mutable cache_type m_memCache; // Values are held in a limited memory cache here.

    mutable frantic::files::filename_sequence
        m_diskCache; // Values removed from the memory cache are stored on disk governed by this object.

    // Object responsible for asynchronously serializing cache entries to disk. This is a shared_ptr because flush_async
    // might require the serializer to exist beyond the scope of the sequence_cache object (ex. A 3dsMax scene object
    // might use the cache, but reset the scene during a flush).
    boost::shared_ptr<serializer_type> m_serializer;

    std::set<double> m_keyUnion; // Stores the union of the keys for the memory cache and the disk cache so we know what
                                 // keys to look for when querying the caches.

    bool m_diskCacheEnabled;

    SizeEstimatorType m_sizeEstimator;
};

template <class V, class S, class E>
class sequence_cache<V, S, E>::deallocate_callback {
    const sequence_cache* m_pOwner;

  public:
    explicit deallocate_callback( const sequence_cache& owner )
        : m_pOwner( &owner ) {}

    std::size_t operator()( double frame, const V& theVal, typename cache_type::deallocator_reason theReason ) const {
        return m_pOwner->on_cache_entry_flushed( frame, theVal, theReason );
    }
};

template <class ValueType, class SerializerType, class SizeEstimatorType>
inline sequence_cache<ValueType, SerializerType, SizeEstimatorType>::sequence_cache(
    std::size_t cacheLimit, const SerializerType& theSerializer, const SizeEstimatorType& theSizeEstimator )
    : m_memCache( cacheLimit )
    , m_diskCacheEnabled( false )
    , m_sizeEstimator( theSizeEstimator ) {
    m_serializer.reset( new serializer_type( theSerializer, cacheLimit / 2 ) );
}

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::set_disk_path( const frantic::tstring& diskSequencePath, sync_option syncOption ) {
    m_keyUnion.clear();
    m_serializer->clear(); // NOTE: If we add an option to copy existing items to the new path, this cannot be done.

    if( syncOption == keep_existing ) {
        // Repopulate the key union with the items in memory.
        for( typename cache_type::const_iterator it = m_memCache.begin(), itEnd = m_memCache.end(); it != itEnd; ++it )
            m_keyUnion.insert( it.key() );
    } else {
        m_memCache.clear(); // Drop cache data without serializing it.
    }

    if( !diskSequencePath.empty() ) {
        m_diskCacheEnabled = true;
        m_diskCache.get_filename_pattern().set( diskSequencePath );
        m_diskCache.get_frame_set().clear();

        if( syncOption == synchronize ) {
            if( m_diskCache.directory_exists() ) {
                m_diskCache.sync_frame_set();

                m_keyUnion.insert( m_diskCache.get_frame_set().begin(), m_diskCache.get_frame_set().end() );
            }
        } else if( syncOption == keep_existing ) {
            // The pattern might not have a valid directory path if the serializer is going to interpret it.
            if( !m_diskCache.get_filename_pattern().get_directory( false ).empty() )
                m_diskCache.create_directory();

            // Start serializing the items to disk.
            for( typename cache_type::const_iterator it = m_memCache.begin(), itEnd = m_memCache.end(); it != itEnd;
                 ++it ) {
                if( !m_serializer->try_enqueue( m_diskCache[it.key()], it.value(), m_sizeEstimator( it.value() ) ) )
                    break;

                m_diskCache.get_frame_set().add_frame( it.key() );
            }
        }
    } else {
        m_diskCacheEnabled = false;
        m_diskCache.get_filename_pattern().set( _T("") );
        m_diskCache.get_frame_set().clear();
    }
}

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::set_capacity( std::size_t memCacheLimit ) {
    m_memCache.set_capacity( memCacheLimit, deallocate_callback( *this ) );
}

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::set_pending_capacity( std::size_t pendingQueueLimit ) {
    // Set the new capacity of the serialization queue. Only wait for it to fit if we are shrinking the capacity. We
    // might be increasing space but still over-allocated and that would be fine.
    m_serializer->set_capacity( pendingQueueLimit, pendingQueueLimit < m_serializer->get_capacity() );
}

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::set_num_serializer_threads( std::size_t numThreads ) {
    m_serializer->set_num_threads( numThreads );
}

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::insert( double frame, const V& theValue ) {
    std::set<double>::iterator itKeyUnion = m_keyUnion.lower_bound( frame );

    if( itKeyUnion != m_keyUnion.end() && *itKeyUnion == frame ) {
        // We already have an entry for this frame, so just update it.
        m_memCache.insert( std::make_pair( frame, theValue ), m_sizeEstimator( theValue ),
                           deallocate_callback( *this ) );

        // We need to update the disk cache if we are replacing a frame.
        if( m_diskCacheEnabled && m_diskCache.get_frame_set().frame_exists( frame ) ) {
            if( !m_serializer->try_enqueue( m_diskCache[frame], theValue, m_sizeEstimator( theValue ) ) )
                m_diskCache.get_frame_set().remove_frame(
                    frame ); // We didn't queue this sample for writing, so we remove it from the valid set on disk.
        }
    } else {
        // We may need to create the disk cache directory if it doesn't already exist, and isn't empty. We might have a
        // valid pattern without a directory if our serializer is going to modify the path or interpret it as a database
        // key for example.
        if( m_diskCacheEnabled && m_keyUnion.empty() &&
            !m_diskCache.get_filename_pattern().get_directory( false ).empty() )
            m_diskCache.create_directory();

        // We don't have an entry already so add it to the memory cache.
        m_keyUnion.insert( itKeyUnion, frame );
        m_memCache.insert( std::make_pair( frame, theValue ), m_sizeEstimator( theValue ),
                           deallocate_callback( *this ) );

        if( m_diskCacheEnabled ) {
            // This only partially schedules the serialization since it doesn't schedule anything if the queue is full.
            // Use flush() or flush_async() to make sure all items are written to the disk cache.
            // TODO: We might want to reserve some space in the serializer for items that need to be serialized (via
            // caused by find/insert). If we fill the queue
            //       speculatively, we might find that we cannot easily drop an item without waiting on the serializer
            //       queue. We could use a priority to solve this problem too.
            if( m_serializer->try_enqueue( m_diskCache[frame], theValue, m_sizeEstimator( theValue ) ) )
                m_diskCache.get_frame_set().add_frame( frame );
        }
    }
}

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::clear() {
    m_keyUnion.clear();
    m_diskCache.get_frame_set().clear();
    m_memCache.clear();
    m_serializer->clear();
}

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::flush( bool emptyMemoryCache ) {
    for( typename cache_type::const_iterator it = m_memCache.begin(), itEnd = m_memCache.end(); it != itEnd; ++it )
        this->on_cache_entry_flushed( it.key(), it.value(), cache_type::erased );

    if( emptyMemoryCache )
        m_memCache.clear();
}

namespace detail {
struct tie_lifetime {
    boost::shared_ptr<void> m_pManaged;

    boost::function<void( void )> m_impl;

    template <class Callable>
    tie_lifetime( const boost::shared_ptr<void>& pManageable, const Callable& impl )
        : m_pManaged( pManageable )
        , m_impl( impl ) {}

    void operator()() const { m_impl(); }
};
} // namespace detail

template <class V, class S, class E>
template <class Callable>
inline void sequence_cache<V, S, E>::flush_async( const Callable& callback ) {
    if( !m_diskCacheEnabled )
        return;

    std::size_t oldCapacity = m_serializer->get_capacity();

    // Make the serialization queue unbounded temporarily so we can add all items to it.
    m_serializer->set_capacity( ( std::numeric_limits<std::size_t>::max )(), false );

    try {
        for( typename cache_type::const_iterator it = m_memCache.begin(), itEnd = m_memCache.end(); it != itEnd;
             ++it ) {
            if( !m_diskCache.get_frame_set().frame_exists( it.key() ) ) {
                frantic::tstring filePath = m_diskCache[it.key()];

                if( !m_serializer->try_enqueue( filePath, it.value(), m_sizeEstimator( it.value() ) ) )
                    throw std::runtime_error( "Failed to enqueue serialization item" );

                m_diskCache.get_frame_set().add_frame( it.key() );
            }
        }
    } catch( ... ) {
        m_serializer->set_capacity( oldCapacity, false );

        throw;
    }

    // Reset the serializer capacity without waiting for it to actually fit the new capacity. We expect that the queue
    // will be over-filled at this point.
    m_serializer->set_capacity( oldCapacity, false );

    // We need to tie the lifetime of the serializer object to the callback so that the serializer isn't deleted until
    // then.
    m_serializer->invoke_on_idle( detail::tie_lifetime( m_serializer, callback ) );
}

template <class V, class S, class E>
class sequence_cache<V, S, E>::clear_callback {
    sequence_cache<V, S, E>* m_pOwner;

  public:
    clear_callback( sequence_cache<V, S, E>& owner )
        : m_pOwner( &owner ) {}

    void operator()( const frantic::tstring& filePath, const V& ) const {
        bool validFrame = false;
        double frame = 0.f;

        if( m_pOwner->m_diskCache.get_filename_pattern().matches_pattern( filePath, validFrame, frame ) && validFrame )
            m_pOwner->m_diskCache.get_frame_set().remove_frame( frame );
    }
};

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::cancel_flush() {
    m_serializer->clear( clear_callback( *this ) );
}

template <class V, class S, class E>
inline frantic::tstring sequence_cache<V, S, E>::get_disk_path_pattern() const {
    return m_diskCacheEnabled ? m_diskCache.get_filename_pattern().get_pattern() : frantic::tstring();
}

template <class V, class S, class E>
inline typename sequence_cache<V, S, E>::const_key_iterator sequence_cache<V, S, E>::key_begin() const {
    return m_keyUnion.begin();
}

template <class V, class S, class E>
inline typename sequence_cache<V, S, E>::const_key_iterator sequence_cache<V, S, E>::key_end() const {
    return m_keyUnion.end();
}

template <class V, class S, class E>
inline typename sequence_cache<V, S, E>::const_memory_iterator sequence_cache<V, S, E>::memory_begin() const {
    return m_memCache.begin();
}

template <class V, class S, class E>
inline typename sequence_cache<V, S, E>::const_memory_iterator sequence_cache<V, S, E>::memory_end() const {
    return m_memCache.end();
}

template <class V, class S, class E>
inline bool sequence_cache<V, S, E>::empty() const {
    return m_keyUnion.empty();
}

template <class V, class S, class E>
std::pair<double, double> sequence_cache<V, S, E>::get_key_range() const {
    if( m_keyUnion.empty() )
        return std::make_pair( std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() );
    return std::make_pair( *m_keyUnion.begin(), *m_keyUnion.rbegin() );
}

template <class V, class S, class E>
inline std::size_t sequence_cache<V, S, E>::get_usage() const {
    return m_memCache.cache_size();
}

template <class V, class S, class E>
inline std::size_t sequence_cache<V, S, E>::get_capacity() const {
    return m_memCache.cache_capacity();
}

template <class V, class S, class E>
inline std::size_t sequence_cache<V, S, E>::get_pending_usage() const {
    return m_serializer->get_usage();
}

template <class V, class S, class E>
inline std::size_t sequence_cache<V, S, E>::get_pending_capacity() const {
    return m_serializer->get_capacity();
}

template <class V, class S, class E>
inline V sequence_cache<V, S, E>::find( double frame ) const {
    V result;

    std::set<double>::const_iterator itKeyUnion = m_keyUnion.lower_bound( frame );

    if( itKeyUnion != m_keyUnion.end() && *itKeyUnion == frame ) {
        typename cache_type::const_iterator itMem = m_memCache.find( frame );

        if( itMem == m_memCache.end() ) {
            if( m_diskCacheEnabled ) {
                frantic::tstring filePath = m_diskCache[frame];

                // This entry is not in memory, check pending serialized values.
                if( m_serializer->find_pending( filePath, result ) ) {
                    m_memCache.insert( std::make_pair( frame, result ), m_sizeEstimator( result ),
                                       deallocate_callback( *this ) );
                } else {
                    // This entry is not in memory and not pending for serialization, so it must already be on disk.
                    assert( m_diskCache.get_frame_set().frame_exists( frame ) );

                    result = m_serializer->get_serializer().deserialize( filePath );

                    m_memCache.insert( std::make_pair( frame, result ), m_sizeEstimator( result ),
                                       deallocate_callback( *this ) );
                }
            }
        } else {
            result = itMem.value();
        }
    }

    return result;
}

namespace detail {
template <class IteratorType>
inline IteratorType previous( IteratorType it ) {
    return --it;
}
} // namespace detail

template <class V, class S, class E>
inline double sequence_cache<V, S, E>::find_nearest_key( double frame ) const {
    if( m_keyUnion.empty() )
        return std::numeric_limits<double>::infinity();

    std::set<double>::const_iterator itKeyUnion = m_keyUnion.upper_bound( frame );

    if( itKeyUnion == m_keyUnion.end() ) {
        return *detail::previous( itKeyUnion );
    } else if( itKeyUnion == m_keyUnion.begin() ) {
        return *m_keyUnion.begin();
    } else {
        std::set<double>::const_iterator itKeyUnionLeft = detail::previous( itKeyUnion );

        if( ( frame - *itKeyUnionLeft ) < ( *itKeyUnion - frame ) )
            return *itKeyUnionLeft;
        else
            return *itKeyUnion;
    }
}

template <class V, class S, class E>
inline std::pair<double, double> sequence_cache<V, S, E>::find_bracketing_keys( double frame ) const {
    if( m_keyUnion.empty() )
        return std::make_pair( -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() );

    std::set<double>::const_iterator itKeyUnion = m_keyUnion.upper_bound( frame );

    if( itKeyUnion == m_keyUnion.end() ) {
        return std::make_pair( *detail::previous( itKeyUnion ), std::numeric_limits<double>::infinity() );
    } else if( itKeyUnion == m_keyUnion.begin() ) {
        return std::make_pair( -std::numeric_limits<double>::infinity(), *m_keyUnion.begin() );
    } else {
        return std::make_pair( *detail::previous( itKeyUnion ), *itKeyUnion );
    }
}

template <class V, class S, class E>
inline void sequence_cache<V, S, E>::wait_for_pending() const {
    if( m_diskCacheEnabled )
        m_serializer->wait_for_idle();
}

template <class V, class S, class E>
inline background_serializer<frantic::tstring, V, S>& sequence_cache<V, S, E>::get_serializer() {
    return *m_serializer;
}

template <class V, class S, class E>
inline std::size_t
sequence_cache<V, S, E>::on_cache_entry_flushed( double frame, const V& theVal,
                                                 typename cache_type::deallocator_reason theReason ) const {
    assert( m_keyUnion.find( frame ) != m_keyUnion.end() );

    if( theReason == cache_type::erased ) {
        if( m_diskCacheEnabled ) {
            if( !m_diskCache.get_frame_set().frame_exists( frame ) ) {
                m_serializer->enqueue( m_diskCache[frame], theVal, m_sizeEstimator( theVal ) );

                m_diskCache.get_frame_set().add_frame( frame );
            }
        } // else{ m_keyUnion.erase( frame ); } // TODO: Do we want to erase the key in this instance?
    }

    return m_sizeEstimator( theVal );
}

} // namespace files
} // namespace frantic
