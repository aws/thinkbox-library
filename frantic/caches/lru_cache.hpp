// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <boost/function.hpp>
#include <boost/next_prior.hpp>
#include <list>
#include <map>

namespace frantic {
namespace caches {

/**
 * A generic least-recently-used (lru) cache container that is reminiscient of std::map<KeyType, ValueType>.
 * \tparam KeyType The type used as an associative identifier when inserting or finding values in the cache.
 * \tparam ValueType The type stored in the cache. Has value semantics, so consider using a shared_ptr<> for objects
 * with expensive copy constructors or assignment operators. Must support copy construction and assignment.
 *
 * \note This container (like the STL) does not concern itself with thread synchronization issues.
 */
template <class KeyType, class ValueType>
class lru_cache {
    typedef std::list<KeyType> usage_tracker_type;

    // Map the key to the corresponding value and entry in the usage tracker.
    // We keep a reference to the usage tracker entry so that we can update it
    // without performing a linear search.
    typedef std::map<KeyType, std::pair<ValueType, typename usage_tracker_type::iterator>> map_type;

  public:
    typedef KeyType key_type;
    typedef ValueType value_type;

    class const_iterator;

    // NOTE: non-const iterators are not provided because modifying the cache entries might invalidate the cache
    // capacity calculations.

  public:
    /**
     * Indicates the reasons the deallocator is being invoked.
     */
    enum deallocator_reason {
        erased,   // Indicates the cache entry is being deleted from the cache to make space for new items
        replaced, // Indicates the cache entry is being replaced with an updated value.
    };

  public:
    /**
     * Constructor
     * \param maxCapacity The maximum number of generic cache units that may be stored before deallocation will occur.
     */
    explicit lru_cache( std::size_t maxCapacity );

    /**
     * \return True if there are no items currently in the cache.
     */
    bool empty() const;

    /**
     * \return The current number of generic cache units that are occupied.
     * \note It is possible for cache_size() > cache_capacity() if a single item is larger than cache_capactity().
     */
    std::size_t cache_size() const;

    /**
     * \return The maximum number of generic cache units that can be occupied without deallocating entries.
     */
    std::size_t cache_capacity() const;

    /**
     * Removes all items from the cache.
     */
    void clear();

    /**
     * Insert a new entry into the cache, or replace an existing entry. If adding the new entry will cause the
     * cache_size to exceed the cache_capacity, items will be deallocated by invoking the 'deallocator' functor.
     *
     * \tparam DeallocatorType A functor type used when deallocating cache enties. Must implement: 'std::size_t
     * operator()( const KeyType&, const ValueType& theVal, deallocator_reason reason ) const' which returns the number
     * of reclaimed generic cache units due to the deallocation. This must match the value provided when the entry was
     * last inserted or undefined behavior will occur.
     *
     * \param newEntry The new entry to add to the cache.
     * \param newEntrySize The size (in generic cache units) of the entry being added.
     * \param deallocator The functor invoked when deallocating entries.
     */
    template <class DeallocatorType>
    void insert( const std::pair<KeyType, ValueType>& newEntry, std::size_t newEntrySize,
                 const DeallocatorType& deallocator );

    /**
     * Changes the maximum capacity of cache, deallocating items as necessary to fit the new size.
     * \tparam DeallocatorType \see insert() for the type requirements.
     * \param newCapacity The new maximum capacity of the cache.
     * \param deallocator The functor invoked when deallocating entries.
     */
    template <class DeallocatorType>
    void set_capacity( std::size_t newCapacity, const DeallocatorType& deallocator );

    /**
     * \return A bidirectional iterator to the beginning of the cache collection.
     */
    const_iterator begin() const;

    /**
     * \return A bidirectional iterator past the end of the cache collection.
     */
    const_iterator end() const;

    /**
     * Finds the value associated with the key if it is currently in the cache.
     *
     * \note Unlike std::map<>, this function is linear (instead of logarithmic) in the number of cache entries due to
     * tracking the least recently used item.
     *
     * \param key The key to search for.
     * \return An iterator to the found item, or this->end() if the item was not in the cache.
     */
    const_iterator find( const KeyType& key ) const;

  private:
    /**
     * Marks the entry associated with 'key' as recently used.
     * \param key The key for the entry to mark.
     */
    void update_usage( const KeyType& key ) const;

  private:
    struct key_compare {
        bool operator()( const typename map_type::value_type& lhs, const KeyType& rhs ) const;
    };

  private:
    map_type m_cacheData;

    mutable usage_tracker_type m_usageTracker;

    std::size_t m_cacheSize;
    std::size_t m_cacheCapacity;
};

template <class K, class V>
class lru_cache<K, V>::const_iterator {
  public:
    const_iterator( typename map_type::const_iterator iterator )
        : m_iterator( iterator ) {}

    const K& key() const { return m_iterator->first; }

    const V& value() const { return m_iterator->second.first; }

    const_iterator& operator++() {
        ++m_iterator;
        return *this;
    }

    const_iterator operator++( int ) {
        const_iterator result( *this );
        ++m_iterator;
        return result;
    }

    bool operator!=( const const_iterator& other ) const { return m_iterator != other.m_iterator; }

    bool operator==( const const_iterator& other ) const { return m_iterator == other.m_iterator; }

  private:
    typename map_type::const_iterator m_iterator;
};

template <class K, class V>
inline lru_cache<K, V>::lru_cache( std::size_t maxCapacity )
    : m_cacheSize( 0 )
    , m_cacheCapacity( maxCapacity ) {}

template <class K, class V>
inline bool lru_cache<K, V>::empty() const {
    return m_cacheSize == 0;
}

template <class K, class V>
inline std::size_t lru_cache<K, V>::cache_size() const {
    return m_cacheSize;
}

template <class K, class V>
inline std::size_t lru_cache<K, V>::cache_capacity() const {
    return m_cacheCapacity;
}

template <class K, class V>
inline void lru_cache<K, V>::clear() {
    m_cacheSize = 0;
    m_usageTracker.clear();
    m_cacheData.clear();
}

template <class K, class V>
inline bool lru_cache<K, V>::key_compare::operator()( const typename map_type::value_type& lhs, const K& rhs ) const {
    return lhs.first < rhs.first;
}

template <class K, class V>
template <class DeallocatorType>
inline void lru_cache<K, V>::insert( const std::pair<K, V>& newEntry, std::size_t newEntrySize,
                                     const DeallocatorType& cb ) {
    typename map_type::iterator newIt = m_cacheData.lower_bound( newEntry.first );

    if( newIt != m_cacheData.end() && newEntry.first == newIt->first ) {
        m_cacheSize = m_cacheSize + newEntrySize - cb( K(), newIt->second.first, replaced );

        newIt->second.first = newEntry.second;

        typename usage_tracker_type::iterator itUsage = newIt->second.second;

        assert( itUsage != m_usageTracker.end() );

        m_usageTracker.erase( itUsage );
        newIt->second.second = m_usageTracker.end();
    } else {
        m_cacheSize += newEntrySize;

        typedef typename map_type::value_type value_t;
        typedef typename map_type::mapped_type mapped_t;

        newIt =
            m_cacheData.insert( newIt, value_t( newEntry.first, mapped_t( newEntry.second, m_usageTracker.end() ) ) );
    }

    // We might need to deallocate cache entries in order to get back to the desired capacity. We haven't added the new
    // entry to the usage tracker list, so we run no risk of deallocating the object we are trying to add.
    while( m_cacheSize > m_cacheCapacity && !m_usageTracker.empty() ) {
        // We need to deallocate items until our cache fits.
        typename map_type::iterator it = m_cacheData.lower_bound( m_usageTracker.front() );

        assert( it != m_cacheData.end() && it->first == m_usageTracker.front() );

        // Use the callback to mark the deallocation and compute the amount of cache space freed.
        m_cacheSize -= cb( it->first, it->second.first, erased );

        m_cacheData.erase( it );
        m_usageTracker.pop_front();
    }

    // If we removed everything, we can check that the various sizes worked out.
    assert( !m_usageTracker.empty() || m_cacheSize == newEntrySize );

    m_usageTracker.push_back( newEntry.first );
    newIt->second.second = boost::prior( m_usageTracker.end() );
}

template <class K, class V>
template <class DeallocatorType>
inline void lru_cache<K, V>::set_capacity( std::size_t newCapacity, const DeallocatorType& cb ) {
    m_cacheCapacity = newCapacity;

    while( m_cacheSize > m_cacheCapacity && !m_usageTracker.empty() ) {
        // We need to deallocate items until our cache fits.
        typename map_type::iterator it = m_cacheData.lower_bound( m_usageTracker.front() );

        assert( it != m_cacheData.end() && it->first == m_usageTracker.front() );

        // Use the callback to mark the deallocation and compute the amount of cache space freed.
        m_cacheSize -= cb( it->first, it->second.first, erased );

        m_cacheData.erase( it );
        m_usageTracker.pop_front();
    }

    // If we removed everything, we can check that the various sizes worked out.
    assert( !m_usageTracker.empty() || m_cacheSize == 0 );
}

template <class K, class V>
inline typename lru_cache<K, V>::const_iterator lru_cache<K, V>::begin() const {
    return m_cacheData.begin();
}

template <class K, class V>
inline typename lru_cache<K, V>::const_iterator lru_cache<K, V>::end() const {
    return m_cacheData.end();
}

template <class K, class V>
inline typename lru_cache<K, V>::const_iterator lru_cache<K, V>::find( const K& key ) const {
    typename map_type::const_iterator it = m_cacheData.lower_bound( key );

    if( it != m_cacheData.end() && key == it->first ) {
        update_usage( key );

        return it;
    } else {
        return m_cacheData.end();
    }
}

template <class K, class V>
inline void lru_cache<K, V>::update_usage( const K& key ) const {
    typename map_type::const_iterator it = m_cacheData.find( key );

    assert( it != m_cacheData.end() );

    typename usage_tracker_type::iterator itUsage = it->second.second;

    assert( itUsage != m_usageTracker.end() );

    m_usageTracker.splice( m_usageTracker.end(), m_usageTracker, itUsage );
}

} // namespace caches
} // namespace frantic
