// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#ifndef SYNCHRONIZED_QUEUE_HPP
#define SYNCHRONIZED_QUEUE_HPP

#include <boost/thread/mutex.hpp>

// pronounced "deck" !
#include <deque>

namespace frantic {
namespace threads {

template <typename Ttype>
class SynchronizedQueue {
    // TODO: should there be a copy constructor/assignment operator?
    SynchronizedQueue( const SynchronizedQueue& );
    SynchronizedQueue& operator=( const SynchronizedQueue& );

  public:
    typedef Ttype value_type;
    SynchronizedQueue() {}
    ~SynchronizedQueue() {}

    // adds an item to the queue
    void enter( Ttype value ) {
        boost::mutex::scoped_lock lock( m_mutex );
        m_queue.push_back( value );
    }
    // removes an item from the queue.
    // returns false if the queue is empty.
    bool leave( Ttype& result ) {
        boost::mutex::scoped_lock lock( m_mutex );

        if( m_queue.empty() )
            return false;

        result = m_queue.front();
        m_queue.pop_front();

        return true;
    }

    std::size_t size() {
        boost::mutex::scoped_lock lock( m_mutex );
        return m_queue.size();
    }

    bool empty() {
        boost::mutex::scoped_lock lock( m_mutex );
        return m_queue.empty();
    }

    void clear() {
        boost::mutex::scoped_lock lock( m_mutex );
        m_queue.clear();
    }

  protected:
    boost::mutex m_mutex;
    std::deque<Ttype> m_queue;
};

} // namespace threads
} // namespace frantic

#endif
