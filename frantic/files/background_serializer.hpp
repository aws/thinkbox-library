// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/exception/error_info.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception_ptr.hpp>
#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <boost/scoped_array.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include <tbb/atomic.h> // Replace with Boost.Atomic once we switch to Boost 1.53 (or std::atomic once we're on C++11)
#include <tbb/concurrent_queue.h> // Replace with Boost.Lockfree perhaps?

#include <deque>
#include <limits>

namespace frantic {
namespace files {

// This type is thrown when an exception is handled in the background thread causing it to close.
class background_serializer_critical_error;

/**
 * Maintains a queue of items that is serialized via a background thread. Items are serialized in FIFO manner. The
 * serialization queue may be given a maximum size which will causes attempts to enqueue that exceed the limit to block.
 *
 * \tparam KeyType The information that identifies each serialized item. Typically interpreted as the path to serialize
 * to. \tparam ValueType The data type that is serialized. \tparam SerializeImplType The type of the functor that
 * implements the actual serialization. Must support this interface: void serialize( KeyType, ValueType ) throw();
 *
 * \note SerializeImplType::serialize() is expected to not throw. You should handle serialization exceptions by storing
 * them (using appropriate synchronization) in the instance and extracting them at an appropriate time.
 */
template <class KeyType, class ValueType, class SerializeImplType>
class background_serializer {
  public:
    /**
     * Constructor.
     * \param theSerializeImpl A copy of this is held and used by the background thread to implement the actual
     * serialization process. \param maxPendingSize The number of cache units to allow active in the pending queue
     * before the next attempt to enqueue will block.
     */
    explicit background_serializer( const SerializeImplType& theSerializeImpl,
                                    std::size_t maxPendingSize = ( std::numeric_limits<std::size_t>::max )() );

    /**
     * Destructor. Blocks until the background thread has finished serializing the queued values.
     */
    ~background_serializer();

    /**
     * \return The maximum size of the pending queue buffer, in generic cache units.
     */
    std::size_t get_capacity() const;

    /**
     * \return The current size of the pending queue buffer, in generic cache units. This may exceed get_capacity()
     * since we always accept at least one item.
     */
    std::size_t get_usage() const;

    /**
     * Clears the pending queue without serializing anything.
     */
    void clear();

    /**
     * Clears the pending queue as in 'clear()' but also invokes the provided callable object with the key & value of
     * each item being cleared. \tparam Type supporting 'void operator()( const KeyType&, const ValueType& ) const'
     * \param cb The callback object applied to each object as it is removed from the queue.
     */
    template <class Callable>
    void clear( const Callable& cb );

    /**
     * Sets the number of threads that are created to process pending items.
     * \note Currently, setting the number of threads to a smaller number than already exist will fail.
     * \param numThreads The total number of threads to create.
     */
    void set_num_threads( std::size_t numThreads );

    /**
     * Sets the maximum size of the pending queue buffer, in generic cache units.
     * \param newCapacity The new capacity for the serialization queue.
     * \param waitUntilNotExceeded If true, this function blocks until the current usage does not exceed the specified
     * capacity. Otherwise it returns regardless.
     */
    void set_capacity( std::size_t newCapacity, bool waitUntilNotExceeded );

    /**
     * Adds a value to the queue for serialization on the background thread. This may block if the enqueued data item
     * exceeds the available space of the queue. If there queue is empty and the item still cannot fit, it will be added
     * anyways. \param serializeDest The key value identifying the value to serialize. Typically the path to serialize
     * to. \param newValue The value to serialize. \param dataSize The size in generic cache units of the value to
     * serialize.
     */
    void enqueue( const KeyType& serializeDest, const ValueType& newValue, std::size_t dataSize );

    /**
     * Attempts to enqueue a new value to the pending queue, but does not block if the queue is full.
     * \param serializeDest The key value identifying the value to serialize. Typically the path to serialize to.
     * \param newValue The value to serialize.
     * \param dataSize The size in generic cache units of the value to serialize.
     * \return True if the item was queued, false if the item could not fit into the available queue space.
     */
    bool try_enqueue( const KeyType& serializeDest, const ValueType& newValue, std::size_t dataSize );

    /**
     * Extracts the first pending value that returns true from the predicate. Returns false if no pending values match.
     * \note This will find the most recently queued item if multiple items with the same key were queued.
     * \param serializeDest The key value identifying the value to serialize. Typically the path to serialize to.
     * \param[out] outValue Will be assigned if there is a pending item with the specified key.
     * \return True if an item was found, false if no pending item had the given key.
     */
    bool find_pending( const KeyType& serializeDest, ValueType& outValue ) const;

    /**
     * Blocks the calling thread until all items that are pending serialization have been processed.
     */
    void wait_for_idle() const;

    /**
     * Invokes the callable object when the serialization queue becomes empty.
     * \tparam Callable Must be a compatible type with the constructor for std::function<void(void)>.
     * \param cb The callable object to invoke
     */
    template <class Callable>
    void invoke_on_idle( const Callable& cb );

    /**
     * Assigns a callback that is invoked after each item is serialized.
     * \note The callback is invoked while holding a lock on the background_serializer's mutex so it must not call any
     *       member functions that might need that lock. We will need to switch to a recursive lock to enable this.
     * \tparam Callable A function object that implements: void operator()( const KeyType& ) const;
     * \param callback The functor to invoke after each items is serialized.
     */
    template <class Callable>
    void set_callback( const Callable& callback );

    /**
     * Accesses the serialization implementation.
     */
    SerializeImplType& get_serializer();

    /**
     * \overload
     */
    const SerializeImplType& get_serializer() const;

  private:
    /**
     * Stores information about a single serialization entry.
     */
    struct pending_entry {
        KeyType key;
        ValueType value;
        std::size_t size;

        pending_entry( const KeyType& key, const ValueType& val, std::size_t dataSize );
    };

    /**
     * Executes in the background thread until signalled to complete.
     */
    void thread_impl( boost::optional<pending_entry>& currentItem );

    /**
     * Wait for the serialization queue to empty, then invokes the callback
     */
    void wait_impl( const boost::function<void( void )>& cb );

    /**
     * Performs a check for any exceptions transferred from the background thread. If present the exception will be
     * rethrown. \pre m_mutex is locked by the current thread
     */
    void check_for_exceptions() const;

    /**
     * Fills in a pending_entry object.
     */
    pending_entry make_entry( const KeyType& key, const ValueType& val, std::size_t dataSize );

    /**
     * Functor for comparing the key of a pending_entry.
     */
    struct key_equals {
        bool operator()( const pending_entry& lhs, const KeyType& rhs ) const;
    };

    /**
     * Struct that holds data used by worker thread instances.
     */
    struct thread_data {
        boost::optional<pending_entry> currentItem; // The item the thread is currently processing
    };

    /**
     * \return True if there are multiple worker threads processing items in the queue.
     */
    bool has_multiple_workers() const;

    /**
     * Finds any thread processing an item with the given key.
     * \param serializeKey The key to search for among the items being processed by worker threads.
     * \return A pointer to the thread_data instance for the thread processing an item with the matching key or NULL if
     * none was found.
     */
    const thread_data* find_thread_processing( const KeyType& serializeKey ) const;

  private:
    typedef std::deque<pending_entry> pending_container;

    SerializeImplType m_theImpl; // Object that is responsible for doing the actual serialization.

    std::size_t m_maxCapacity;           // Maximum size of pending value queue, in generic cache units.
    tbb::atomic<std::size_t> m_curUsage; // Current size of '       '     '    , '  '       '     '
    tbb::atomic<bool> m_wantsExit;     // If true, the worker threads should exit as soon as the pending queue empties.
    tbb::atomic<bool> m_hasError;      // If true, one of the worker threads had an unhandled error and closed.
    pending_container m_pendingValues; // Queue of values to serialize

    std::size_t m_numThreads;
    boost::scoped_array<boost::thread> m_threads; // Collection of thread objects.
    std::deque<thread_data> m_threadData;         // The items currently processed by the associated thread.

    mutable boost::mutex m_mutex; // Protects access to the queue. Mutable because even the const member functions must
                                  // lock to inspect the pending queue.
    mutable boost::condition_variable m_next, m_avail; // Conditions for data being available and queue space being
                                                       // available respectively. Mutable due to wait_for_idle().
    mutable boost::exception_ptr
        m_threadException; // The exception associated with the thread. If m_hasError == true then this will be
                           // populated with an exception unhandled in one of the worker threads.

    boost::function<void( const KeyType& )> m_callback; // Callback object invoked after serializing each item.
};

template <class K, class V, class S>
inline background_serializer<K, V, S>::background_serializer( const S& theSerializeImpl, std::size_t maxPendingSize )
    : m_theImpl( theSerializeImpl )
    , m_maxCapacity( maxPendingSize ) {
    m_curUsage = 0;
    m_numThreads = 0;
    m_wantsExit = false;
    m_hasError = false;

    this->set_num_threads( 1u );
}

template <class K, class V, class S>
inline background_serializer<K, V, S>::~background_serializer() {
    m_wantsExit =
        true; // Atomically indicate that the background thread should exit when the pending queue becomes empty.
    m_next.notify_all(); // Wake the background thread if it was waiting for more input.

    for( std::size_t i = 0; i < m_numThreads; ++i )
        m_threads[i].join();
}

template <class K, class V, class S>
inline bool background_serializer<K, V, S>::has_multiple_workers() const {
    return m_numThreads > 1u;
}

template <class K, class V, class S>
inline const typename background_serializer<K, V, S>::thread_data*
background_serializer<K, V, S>::find_thread_processing( const K& serializeKey ) const {
    for( typename std::deque<thread_data>::const_iterator it = m_threadData.begin(), itEnd = m_threadData.end();
         it != itEnd; ++it ) {
        if( it->currentItem && it->currentItem.get().key == serializeKey )
            return &*it;
    }

    return NULL;
}

template <class K, class V, class S>
inline std::size_t background_serializer<K, V, S>::get_capacity() const {
    return m_maxCapacity;
}

template <class K, class V, class S>
inline std::size_t background_serializer<K, V, S>::get_usage() const {
    return m_curUsage;
}

template <class K, class V, class S>
inline void background_serializer<K, V, S>::set_capacity( std::size_t newCapacity, bool waitUntilNotExceeded ) {
    boost::unique_lock<boost::mutex> theLock( m_mutex );

    m_maxCapacity = newCapacity;

    if( waitUntilNotExceeded ) {
        this->check_for_exceptions();

        while( m_curUsage > m_maxCapacity ) {
            m_avail.wait( theLock );

            this->check_for_exceptions();
        }
    }
}

namespace detail {
template <class K, class V>
struct null_op {
    inline void operator()( const K&, const V& ) const {}
};
} // namespace detail

template <class K, class V, class S>
inline void background_serializer<K, V, S>::clear() {
    this->clear( detail::null_op<K, V>() );
}

template <class K, class V, class S>
template <class Callable>
inline void background_serializer<K, V, S>::clear( const Callable& cb ) {
    boost::lock_guard<boost::mutex> theLock( m_mutex );

    while( !m_pendingValues.empty() ) {
        cb( m_pendingValues.front().key, m_pendingValues.front().value );

        m_curUsage -= m_pendingValues.front().size;

        m_pendingValues.pop_front();
    }
}

template <class K, class V, class S>
void background_serializer<K, V, S>::set_num_threads( std::size_t numThreads ) {
    if( numThreads < 1u )
        throw std::invalid_argument( "Invalid argument for set_num_threads()" );

    if( numThreads > m_numThreads ) {
        boost::scoped_array<boost::thread> newThreads( new boost::thread[numThreads] );

        m_threadData.resize( numThreads );

        // Move the old threads
        for( std::size_t i = 0; i < m_numThreads; ++i )
            newThreads[i] = boost::move( m_threads[i] );

        // Initialize the new threads by swapping them into place.
        for( std::size_t i = m_numThreads; i < numThreads; ++i )
            boost::thread( &background_serializer::thread_impl, this, boost::ref( m_threadData[i].currentItem ) )
                .swap( newThreads[i] );

        m_threads.swap( newThreads );
        m_numThreads = numThreads;
    } else if( numThreads < m_numThreads ) {
        // We need to stop some threads.
        for( std::size_t i = numThreads; i < m_numThreads; ++i ) {
            m_threads[i].interrupt();
            m_threads[i].join();
        }

        boost::scoped_array<boost::thread> newThreads( new boost::thread[numThreads] );

        m_threadData.resize( numThreads );

        // Move the old, valid threads
        for( std::size_t i = 0; i < numThreads; ++i )
            newThreads[i] = boost::move( m_threads[i] );

        m_threads.swap( newThreads );
        m_numThreads = numThreads;
    }
}

template <class K, class V, class S>
inline void background_serializer<K, V, S>::enqueue( const K& serializeDest, const V& newValue, std::size_t dataSize ) {
    assert( !m_wantsExit );

    boost::unique_lock<boost::mutex> theLock( m_mutex );

    this->check_for_exceptions();

    // If we have multiple worker threads, we need to make sure that we don't add an item key that already is pending.
    // If we did that would cause the serializer to write to the same target (ie. file) in parallel which would be bad.
    if( this->has_multiple_workers() ) {
        typename pending_container::iterator itPending =
            std::find_if( m_pendingValues.begin(), m_pendingValues.end(),
                          std::bind( key_equals(), std::placeholders::_1, boost::cref( serializeDest ) ) );

        if( itPending != m_pendingValues.end() ) {
            // We need to replace the already pending item, which always succeeds regardless of capacity.
            m_curUsage = m_curUsage - itPending->size + dataSize;

            itPending->value = newValue;
            itPending->size = dataSize;

            return;
        } else {
            // If any thread is currently processing an item with the same key-value, we need to wait until that is not
            // true. We only need to inspect the one thread because the pending queue would never have two items with
            // the same key.
            if( const thread_data* pThreadData = this->find_thread_processing( serializeDest ) ) {
                do {
                    m_avail.wait( theLock );

                    this->check_for_exceptions();
                } while( pThreadData->currentItem && pThreadData->currentItem.get().key == serializeDest );
            }
        }
    }

    // Wait until there is space in the queue.
    while( m_curUsage + dataSize > m_maxCapacity && !m_pendingValues.empty() ) {
        m_avail.wait( theLock );

        this->check_for_exceptions();
    }

    m_curUsage += dataSize;

    // We use push_front so that find_pending will find the most recently pushed value (since it searches front to back)
    m_pendingValues.push_front( make_entry( serializeDest, newValue, dataSize ) );
    m_next.notify_one();
}

template <class K, class V, class S>
inline bool background_serializer<K, V, S>::try_enqueue( const K& serializeDest, const V& newValue,
                                                         std::size_t dataSize ) {
    assert( !m_wantsExit );

    boost::lock_guard<boost::mutex> theLock( m_mutex );

    this->check_for_exceptions();

    // If we have multiple worker threads, we need to make sure that we don't add an item key that already is pending.
    // If we did that would cause the serializer to write to the same target (ie. file) in parallel which would be bad.
    if( this->has_multiple_workers() ) {
        typename pending_container::iterator itPending =
            std::find_if( m_pendingValues.begin(), m_pendingValues.end(),
                          std::bind( key_equals(), std::placeholders::_1, boost::cref( serializeDest ) ) );

        if( itPending != m_pendingValues.end() ) {
            // We need to replace the already pending item, which always succeeds regardless of capacity.
            m_curUsage = m_curUsage - itPending->size + dataSize;

            itPending->value = newValue;
            itPending->size = dataSize;

            return true;
        } else if( this->find_thread_processing( serializeDest ) != NULL ) {
            // If any thread is currently processing an item with the same key-value, we cannot add the new item.
            return false;
        }
    }

    if( m_curUsage + dataSize > m_maxCapacity && !m_pendingValues.empty() )
        return false;

    m_curUsage += dataSize;

    // We use push_front so that find_pending will find the most recently pushed value (since it searches front to back)
    m_pendingValues.push_front( make_entry( serializeDest, newValue, dataSize ) );
    m_next.notify_one();

    return true;
}

template <class K, class V, class S>
inline bool background_serializer<K, V, S>::key_equals::operator()( const pending_entry& lhs, const K& rhs ) const {
    return lhs.key == rhs;
}

template <class K, class V, class S>
inline background_serializer<K, V, S>::pending_entry::pending_entry( const K& _key, const V& _val,
                                                                     std::size_t _dataSize )
    : key( _key )
    , value( _val )
    , size( _dataSize ) {}

template <class K, class V, class S>
inline typename background_serializer<K, V, S>::pending_entry
background_serializer<K, V, S>::make_entry( const K& key, const V& val, std::size_t dataSize ) {
    return pending_entry( key, val, dataSize );
}

template <class K, class V, class S>
inline bool background_serializer<K, V, S>::find_pending( const K& serializeDest, V& outValue ) const {
    boost::lock_guard<boost::mutex> theLock( m_mutex );

    typename pending_container::const_iterator it =
        std::find_if( m_pendingValues.begin(), m_pendingValues.end(),
                      std::bind( key_equals(), std::placeholders::_1, boost::cref( serializeDest ) ) );

    if( it == m_pendingValues.end() ) {
        // Check the worker threads to see if any of the items being processed are the one we are looking for
        if( const thread_data* pThreadData = this->find_thread_processing( serializeDest ) ) {
            outValue = pThreadData->currentItem.get().value;
            return true;
        }

        return false;
    }

    outValue = it->value;

    return true;
}

template <class K, class V, class S>
inline void background_serializer<K, V, S>::wait_for_idle() const {
    boost::unique_lock<boost::mutex> theLock( m_mutex );

    this->check_for_exceptions();

    while( !m_pendingValues.empty() ) {
        m_avail.wait( theLock );

        this->check_for_exceptions();
    }

    // Make sure nothing is being actively processed either.
    for( typename std::deque<thread_data>::const_iterator it = m_threadData.begin(), itEnd = m_threadData.end();
         it != itEnd; ++it ) {
        while( it->currentItem ) {
            m_avail.wait( theLock );

            this->check_for_exceptions();
        }
    }
}

template <class K, class V, class S>
inline void background_serializer<K, V, S>::wait_impl( const boost::function<void( void )>& cb ) {
    this->wait_for_idle();

    cb();
}

template <class K, class V, class S>
template <class Callable>
inline void background_serializer<K, V, S>::invoke_on_idle( const Callable& cb ) {
    boost::thread waitThread( &background_serializer::wait_impl, this, cb );

    waitThread.detach();
}

template <class K, class V, class S>
template <class Callable>
inline void background_serializer<K, V, S>::set_callback( const Callable& callback ) {
    boost::lock_guard<boost::mutex> theLock( m_mutex );

    m_callback = callback;
}

template <class K, class V, class S>
inline S& background_serializer<K, V, S>::get_serializer() {
    return m_theImpl;
}

template <class K, class V, class S>
inline const S& background_serializer<K, V, S>::get_serializer() const {
    return m_theImpl;
}

class background_serializer_critical_error : public virtual boost::exception, public virtual std::exception {
  public:
    typedef boost::error_info<struct tag_nested_exception, boost::exception_ptr> nested_exception_info;

    explicit background_serializer_critical_error( const boost::exception_ptr& pE ) {
        *this << nested_exception_info( pE );
    }

    virtual const char* what() const throw() { return boost::diagnostic_information_what( *this ); }
};

template <class K, class V, class S>
inline void background_serializer<K, V, S>::thread_impl( boost::optional<pending_entry>& currentItem ) {
    boost::unique_lock<boost::mutex> theLock( m_mutex, boost::defer_lock );

    try {
        theLock.lock();

        // We loop until 'm_wantsExit' is set to true AND the queue has been exhausted.
        while( !m_wantsExit ) {
            if( m_pendingValues.empty() )
                m_next.wait( theLock );

            // Process all the items in the queue
            while( !m_pendingValues.empty() ) {
                // Copy the item to serialize into the thread's 'currentItem' location.
                currentItem = m_pendingValues.back();

                // Remove the item from the queue
                m_pendingValues.pop_back();

                // Unlock since we're done with the queue until later and 'serialize' is a potentially long running
                // function.
                theLock.unlock();

                // We specify that serialize DOES NOT THROW. If it does so anyways the thread will close and report the
                // unhandled exception.
                m_theImpl.serialize( currentItem.get().key, currentItem.get().value );

                theLock.lock();

                if( m_callback )
                    m_callback( currentItem.get().key );

                std::size_t theSize = currentItem.get().size; // Make a copy of the size since we want to delete the
                                                              // item before decrementing 'm_curUsage'.

                currentItem = boost::none; // Assigning boost::none will destruct the held object

                m_curUsage -= theSize;

                // We have space to store the next item in the queue. Not ideal to notify a thread that will block on a
                // lock currently held here... We use notify_all() because its possible that multiple threads are
                // waiting to add items to the queue (or for the queue to empty).
                m_avail.notify_all();

                // The thread can be interrupted here to exit early.
                boost::this_thread::interruption_point();
            }
        }
    } catch( const boost::thread_interrupted& ) {
        // We explicitly list thread_interrupted so that the catch( ... ) doesn't process it.
        throw;
    } catch( ... ) {
        // If we didn't have the lock (perhaps due to the exception originating in m_theImpl.serialize().), then we need
        // to hold the lock to prevent a race condition on the exception_ptr.
        if( !theLock.owns_lock() )
            theLock.lock();

        // If this is the first thread to experience a fatal exception, flag it and store the exception_ptr.
        if( m_hasError.compare_and_swap( true, false ) == false )
            // NOTE: boost::current_exception() allows for uniform handling of exceptions. No need to specialize the
            // catch signature unless we actually handle them.
            m_threadException = boost::current_exception();

        // We may need to notify the main thread that we experienced a critical error, otherwise it could block
        // permanently.
        m_avail.notify_all();
    }
}

template <class K, class V, class S>
inline void background_serializer<K, V, S>::check_for_exceptions() const {
    // Precondition: 'm_mutex' is locked by the calling thread.

    // If a background thread has reported an exception, rethrow it here.
    if( m_hasError ) {
        // If we previously reported the exception, report a generic error now.
        if( !m_threadException )
            boost::throw_exception(
                std::runtime_error( "background_serializer is unavailable due to previous errors" ) );

        boost::exception_ptr tempExceptionPtr( m_threadException );
        m_threadException = boost::exception_ptr();

        boost::throw_exception( background_serializer_critical_error( tempExceptionPtr ) );
    }
}

} // namespace files
} // namespace frantic
