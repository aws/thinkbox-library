// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#pragma warning( push, 3 )
#pragma warning( disable : 4706 4002 )

#include <boost/exception/all.hpp>
#include <boost/exception_ptr.hpp>
#pragma warning( push, 3 )
#include <boost/thread.hpp>
#pragma warning( pop )

#include <tbb/atomic.h>
#include <tbb/concurrent_queue.h>
#include <tbb/task_scheduler_init.h>

#pragma warning( pop )

namespace frantic {
namespace threads {

#pragma warning( push )
#pragma warning( disable : 4127 )

namespace {

// Copy the exception message into this type because our current compiler does not provide the necessary support to
// capture the original exception across threads. This handler will capture the msg of a std::exception derived
// exception. NOTE: If capturing the exact type of some exceptions is important, make a catch handler above this that
// can handle copying those specific exceptions.
struct thread_exception : virtual boost::exception, virtual std::exception {
    typedef boost::error_info<thread_exception, std::string> info_type;

    virtual const char* what() const throw() {
        if( const std::string* msg = boost::get_error_info<info_type>( *this ) )
            return msg->c_str();
        return "<unknown thread exception>";
    }
};

} // namespace

template <class ProducerConsumerModel>
class buffered_producer_consumer {
    typedef typename ProducerConsumerModel::buffer_type buffer_type;

    /**
     * This class is designed to hold a temporary ptr to an object, deleting it if destroyed while holding an item.
     * I made this instead of using std::unique_ptr because I wanted to be able to set the stored ptr by assigning to
     * a reference (tbb::concurrent_queue::push for example).
     */
    class ptr_guard {
        buffer_type* item;

      public:
        ptr_guard() { item = NULL; }
        ~ptr_guard() {
            if( item )
                delete item;
        }

        inline buffer_type* get() { return item; } // Get a ptr to the currenly owned item
        inline buffer_type*& get_ref() {
            assert( !item );
            return item;
        } // Returns a ref so this pointer can be set. Must not be called with item is non-NULL.
        inline buffer_type* release() {
            buffer_type* result = item;
            item = NULL;
            return result;
        } // Releases ownership of the current item.
    };

    // A counter of the number of active threads. When it becomes zero, all producer threads have finished.
    tbb::atomic<int> numThreadsRemaining;

    // Will be atomically set to true when a thread has registered an exception.
    tbb::atomic<bool> errorOccurred;

    // An exception ptr for transferring the exception (or more likely an approximation) to the main thread.
    boost::exception_ptr error;

    // A pair of queues for passing empty and full buffers between producers and the consumer. The producer threads will
    // be responsible for creating the initial buffers.
    tbb::concurrent_queue<buffer_type*> emptyItems, fullItems;

    template <class T>
    inline static bool concurrent_queue_try_pop( tbb::concurrent_queue<T>& queue, T& outValue ) {
#if TBB_VERSION_MAJOR < 3
        return queue.pop_if_present( outValue );
#else
        return queue.try_pop( outValue );
#endif
    }

    /**
     * This is the code executed by worker threads.
     * @param prod The producer meta-object used to create buffers and Producer::thread_instance objects that do the
     * actual production.
     */
    void producer_fn() {
        try {
            typedef typename ProducerConsumerModel::producer_instance producer_instance;
            producer_instance threadProd( *m_pcModel );

            tbb::task_scheduler_init tsched;

            ptr_guard theItem;

            // Create an extra buffer for later use by this thread.
            theItem.get_ref() = m_pcModel->create_buffer();
            emptyItems.push( theItem.release() );

            // Create the buffer we will be filling.
            theItem.get_ref() = m_pcModel->create_buffer();
            threadProd.init_buffer( theItem.get() );

            while( 1 ) {
                boost::this_thread::interruption_point();

                if( !threadProd.can_produce_more() ) {
                    // This thread should exit since there is no more data to produce. We may need to flush if our
                    // current buffer is non-empty.
                    if( !threadProd.is_buffer_empty( theItem.get() ) ) {
                        threadProd.finish_buffer( theItem.get() );
                        fullItems.push( theItem.release() );
                    } else {
                        emptyItems.push( theItem.release() );
                    }

                    break;
                }

                if( threadProd.is_buffer_full( theItem.get() ) ) {
                    threadProd.finish_buffer( theItem.get() );
                    fullItems.push( theItem.release() );

                    while( !concurrent_queue_try_pop( emptyItems, theItem.get_ref() ) ) {
                        boost::this_thread::interruption_point();
                        boost::this_thread::yield();
                    }

                    threadProd.init_buffer( theItem.get() );
                }

                // We have a non-full buffer so produce data to go into it.
                threadProd.fill_buffer( theItem.get() );
            }
        } catch( const boost::thread_interrupted& ) {
            ; // Do nothing
        } catch( const std::exception& e ) {

            // Atomically determine if another thread has already thrown an exception, and if not store our exception.
            // Otherwise ignore it in favor of the already stored exception.
            if( errorOccurred.compare_and_swap( true, false ) == false )
                error = boost::copy_exception( thread_exception() << thread_exception::info_type( e.what() ) );
        }

        --numThreadsRemaining;
        return;
    }

  private:
    ProducerConsumerModel* m_pcModel;

    boost::thread_group m_threads;

  public:
    buffered_producer_consumer() {
        numThreadsRemaining = 0;
        errorOccurred = false;
    }

    ~buffered_producer_consumer() {
        m_threads.interrupt_all();
        m_threads.join_all();

        buffer_type* item;
        while( concurrent_queue_try_pop( emptyItems, item ) )
            delete item;
        while( concurrent_queue_try_pop( fullItems, item ) )
            delete item;
    }

    void reset( ProducerConsumerModel& pcModel, unsigned int numThreads = 0 ) {
        m_pcModel = &pcModel;

        numThreads = std::max( 1u, numThreads == 0 ? boost::thread::hardware_concurrency() - 1u : numThreads );

        // Set the atomic counter to track the number of outstanding worker threads.
        numThreadsRemaining = numThreads;
        errorOccurred = false;

        for( unsigned int i = 0; i < numThreads; ++i )
            m_threads.create_thread( boost::bind( &buffered_producer_consumer::producer_fn, this ) );
    }

    /**
     * This function will consume data produced by the worker threads until they have all exited.
     * @param cons The consumer implementation object that will process filled buffers.
     */
    void run( bool untilDone = true ) {
        try {
            typedef typename ProducerConsumerModel::consumer_instance consumer_instance;
            consumer_instance threadCons( *m_pcModel );

            ptr_guard theItem;

            do {
                if( !concurrent_queue_try_pop( fullItems, theItem.get_ref() ) ) {
                    threadCons.do_idle_process();

                    while( !concurrent_queue_try_pop( fullItems, theItem.get_ref() ) ) {
                        if( errorOccurred )
                            throw boost::thread_interrupted(); // Can't throw this->error yet, since it may not be
                                                               // finished being created. Throw this instead and sync
                                                               // with all threads.
                        if( numThreadsRemaining == 0 ) {
                            if( concurrent_queue_try_pop(
                                    fullItems,
                                    theItem.get_ref() ) ) // Check again to make sure we didn't get an item before
                                                          // #threads went to 0.
                                break;
                            return;
                        }
                        boost::this_thread::yield();
                    }
                }

                threadCons.consume_buffer( theItem.get() );

                emptyItems.push( theItem.release() );
            } while( untilDone );
        } catch( const boost::thread_interrupted& ) {
            m_threads.interrupt_all();
            m_threads.join_all();

            boost::rethrow_exception( error );
        } catch( ... ) {
            m_threads.interrupt_all();
            m_threads.join_all();

            throw;
        }
    }

    void return_finished_item( std::unique_ptr<buffer_type> item ) { emptyItems.push( item.release() ); }

    std::unique_ptr<buffer_type> steal_finished_item( bool waitForItem = true ) {
        std::unique_ptr<buffer_type> result;

        try {
            buffer_type* theItem = NULL;

            if( waitForItem ) {
                while( !concurrent_queue_try_pop( fullItems, theItem ) ) {
                    if( errorOccurred )
                        throw boost::thread_interrupted(); // Can't throw this->error yet, since it may not be finished
                                                           // being created. Throw this instead and sync with all
                                                           // threads.
                    if( numThreadsRemaining == 0 ) {
                        concurrent_queue_try_pop(
                            fullItems,
                            theItem ); // Check again to make sure we didn't get an item before #threads went to 0.
                        break;
                    }
                    boost::this_thread::yield();
                }
            } else if( !concurrent_queue_try_pop( fullItems, theItem ) && errorOccurred ) {
                throw boost::thread_interrupted(); // Can't throw this->error yet, since it may not be finished being
                                                   // created. Throw this instead and sync with all threads.
            }

            result.reset( theItem );
        } catch( const boost::thread_interrupted& ) {
            m_threads.interrupt_all();
            m_threads.join_all();

            boost::rethrow_exception( error );
        } catch( ... ) {
            m_threads.interrupt_all();
            m_threads.join_all();

            throw;
        }

        return result;
    }

    bool is_done() { return ( numThreadsRemaining == 0 ); }
};

/**
 * This algorithm will operate the Poducer/Consumer pattern for the given template arguments. The construction of the
 * types is defined above.
 *
 * This algorithm will (in parallel) fill buffer objects using ProducerConsumerModel::producer_instance::fill_buffer()
 * and pass them serially to a single Consumer object via. ProducerConsumerModel::consumer_instance::consume_buffer.
 * Once the buffer is consumed it will be made available to another ProducerConsumerModel::producer_instance.
 *
 * @note ProducerConsumerModel::consumer_instance::consume_buffer() will always be called in the context of the thread
 * calling buffered_producer_consumer().
 *
 * @param pcModel An instance of ProducerConsumerModel that supports the interface described above.
 * @param numWorkers The number of worker threads to use, or 0 if the system should decide.
 */
template <class ProducerConsumerModel>
void do_buffered_producer_consumer( ProducerConsumerModel& pcModel, unsigned int numWorkers = 0 ) {
    buffered_producer_consumer<ProducerConsumerModel> theImpl;
    theImpl.reset( pcModel, numWorkers );
    theImpl.run();
}

#pragma warning( pop )

} // namespace threads
} // namespace frantic
