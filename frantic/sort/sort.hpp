// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

//#include <malloc.h>
//#include <process.h>
#include <boost/smart_ptr.hpp>

#pragma warning( push, 3 )
#pragma warning( disable : 4512 4100 )
//#include <tbb/task.h>
#include <tbb/atomic.h>
#include <tbb/parallel_sort.h>
#include <tbb/tbb_thread.h>
#pragma warning( pop )

#define PARALLEL_CUTOFF 10000

// NOTE: You should include windows.h before this file to get threaded_frantic::sort()

namespace frantic {
namespace sort {

namespace detail {
struct DefaultRandIterTraits {
    template <class RandIter>
    static void* get_data_pointer( RandIter p ) {
        return *p;
    }

    template <class RandIter>
    static std::size_t get_data_size( RandIter p ) {
        return p.structure_size();
    }
};
} // namespace detail

/**
 * Swaps the n-byte structure pointed at by iterators p1 and p2. Copies in 32/64 bit chunks if it can.
 */
inline void byte_swap( void* p1, void* p2, std::size_t n ) {
    typedef std::size_t chunk_t; // Use 64-bit on x64, 32-bit on x86
                                 // TODO: Does it make more sense to use 32bit chunks? There is a significant penalty
                                 // for structures n%8 != 0 on 64bit which is fairly common. Most structs are n%4 == 0.

    chunk_t tI, *pA = reinterpret_cast<chunk_t*>( p1 ), *pB = reinterpret_cast<chunk_t*>( p2 );
    for( std::size_t i = 0; i < n / sizeof( chunk_t ); ++i, ++pA, ++pB ) {
        tI = *pA;
        *pA = *pB;
        *pB = tI;
    }

    char tC, *pCA = reinterpret_cast<char*>( pA ), *pCB = reinterpret_cast<char*>( pB );
    for( std::size_t i = 0; i < n % sizeof( chunk_t ); ++i, ++pCA, ++pCB ) {
        tC = *pCA;
        *pCA = *pCB;
        *pCB = tC;
    }
}

/**
 * Swaps the data at the two iterators. get_data_pointer is needed as sometimes this is sorting an array of pointers
 * and we just want the value from the iterator; other times we are sorting an array of objects and we need the
 * addresses of the objects themselves.
 */
template <class RandIter, class RandIterTraits>
inline void iter_swap( RandIter p1, RandIter p2, std::size_t n ) {
    byte_swap( RandIterTraits::get_data_pointer( p1 ), RandIterTraits::get_data_pointer( p2 ), n );
}

/**
 * Function override to get around the inability to pass default template arguments into function templates outside of a
 * class.
 */
template <class RandIter>
inline void iter_swap( RandIter p1, RandIter p2, std::size_t n ) {
    iter_swap<RandIter, detail::DefaultRandIterTraits>( p1, p2, n );
}

/**
 * Returns an iterator that is the median of the three supplied iterators.
 */
template <class RandIter, class Pred>
inline RandIter med3( RandIter a, RandIter b, RandIter c, const Pred& p ) {
    return p( *a, *b ) ? ( p( *b, *c ) ? b : p( *a, *c ) ? c : a ) : ( p( *c, *b ) ? b : p( *c, *a ) ? c : a );
}

/**
 * A O(n^2) sort, of the insertion variety. Use this for very small arrays.
 */
template <class RandIter, class Pred, class RandIterTraits>
inline void insertion_sort( RandIter begin, RandIter end, Pred cmp, std::size_t es ) {
    // TODO: This isn't really an insertion sort. It appears to be a bubble sort.
    for( RandIter pm = begin; pm != end; ++pm )
        for( RandIter pl = pm; pl != begin && cmp( *pl, *( pl - 1 ) ); --pl )
            iter_swap<RandIter, RandIterTraits>( pl, pl - 1, es );
}

template <class RandIter, class Pred>
inline void insertion_sort( RandIter begin, RandIter end, Pred cmp ) {
    insertion_sort<RandIter, Pred, detail::DefaultRandIterTraits>(
        begin, end, cmp, detail::DefaultRandIterTraits::get_data_size( begin ) );
}

/**
 * This function partitions the range bounded by the iterators begin and end, such that
 * it produces the best balanced partition. Returns the number of items less than the pivot,
 * and the number of items greater than the pivot. TODO: Pass in the pivot?
 */
template <class RandIter, class Pred, class ValueType, class RandIterTraits>
inline std::pair<std::size_t, std::size_t> partition( RandIter begin, RandIter end, Pred cmp, std::size_t es ) {
    std::size_t n = ( end - begin );

    if( n <= 1 )
        return std::pair<std::size_t, std::size_t>( 0, 0 );

    { // Choose the pivot and swap it into the leftmost location
        RandIter left, mid, right;
        mid = begin + ( n / 2 );
        if( n > 7 ) {
            left = begin;
            right = end - 1;
            if( n > 40 ) { /* Big arrays, pseudomedian of 9 */
                std::size_t s = ( n / 8 );
                left = med3( left, left + s, left + 2 * s, cmp );
                mid = med3( mid - s, mid, mid + s, cmp );
                right = med3( right - 2 * s, right - s, right, cmp );
            }
            mid = med3( left, mid, right, cmp ); /* Mid-size, med of 3 */
        }

        iter_swap<RandIter, RandIterTraits>( begin, mid, es );
    }

    ValueType pivot = *begin;
    RandIter pa = begin, pb = begin;
    RandIter pc = end - 1, pd = end - 1;

    for( ;; ) {
        while( pb <= pc && !cmp( pivot, *pb ) ) {
            if( !cmp( *pb, pivot ) ) {
                iter_swap<RandIter, RandIterTraits>( pa, pb, es );
                ++pa;
            }
            ++pb;
        }
        while( pc >= pb && !cmp( *pc, pivot ) ) {
            if( !cmp( pivot, *pc ) ) {
                iter_swap<RandIter, RandIterTraits>( pc, pd, es );
                --pd;
            }
            --pc;
        }
        if( pb > pc )
            break;
        iter_swap<RandIter, RandIterTraits>( pb, pc, es );
        ++pb;
        --pc;
    }

    std::size_t s;

    s = std::min<std::size_t>( pa - begin, pb - pa );
    for( RandIter a = begin, b = pb - s; s > 0; --s, ++a, ++b )
        iter_swap<RandIter, RandIterTraits>( a, b, es );

    s = std::min<std::size_t>( pd - pc, end - pd - 1 );
    for( RandIter a = pb, b = end - s; s > 0; --s, ++a, ++b )
        iter_swap<RandIter, RandIterTraits>( a, b, es );

    return std::pair<std::size_t, std::size_t>( pb - pa, pd - pc );
}

template <class RandIter, class Pred>
inline std::pair<std::size_t, std::size_t> partition( RandIter begin, RandIter end, Pred cmp ) {
    return partition<RandIter, Pred, char*, detail::DefaultRandIterTraits>(
        begin, end, cmp, detail::DefaultRandIterTraits::get_data_size( begin ) );
}

/**
 * A single threaded quicksort that degenerates into an insertion sort for small arrays. This sort
 * exists because the std::sort only works on items that can be copied by value. (see std::_Insertion_sort1 for details
 * why)
 */
template <class RandIter, class Pred, class ValueType, class RandIterTraits>
inline void sort( RandIter begin, RandIter end, Pred cmp, std::size_t elementSize ) {
    if( ( end - begin ) < 7 ) { /* Insertion sort on smallest arrays */
        insertion_sort<RandIter, Pred, RandIterTraits>( begin, end, cmp, elementSize );
        return;
    }

    std::pair<std::size_t, std::size_t> partSizes =
        frantic::sort::partition<RandIter, Pred, ValueType, RandIterTraits>( begin, end, cmp, elementSize );

    if( partSizes.first > 1 )
        frantic::sort::sort<RandIter, Pred, ValueType, RandIterTraits>( begin, begin + partSizes.first, cmp,
                                                                        elementSize );
    if( partSizes.second > 1 )
        frantic::sort::sort<RandIter, Pred, ValueType, RandIterTraits>( end - partSizes.second, end, cmp, elementSize );
}

template <class RandIter, class Pred>
inline void sort( RandIter begin, RandIter end, Pred cmp ) {
    sort<RandIter, Pred, char*, detail::DefaultRandIterTraits>( begin, end, cmp,
                                                                detail::DefaultRandIterTraits::get_data_size( begin ) );
}

namespace detail {
struct progress_info {
    tbb::tbb_thread::id m_mainThreadId;
    tbb::atomic<std::size_t> m_currentProgress;
    std::size_t m_totalProgress;

    frantic::logging::progress_logger* m_progress;
};

template <class RandIter, class Pred, class ValueType, class RandIterTraits>
class sort_task : public tbb::task {
    std::pair<RandIter, RandIter> m_range;
    Pred m_pred;

    progress_info& m_progInfo;

    std::size_t m_elementSize;

    void add_progress( std::size_t progress ) {
        std::size_t currentProgress = m_progInfo.m_currentProgress.fetch_and_add( progress ) + progress;
        if( tbb::this_tbb_thread::get_id() == m_progInfo.m_mainThreadId )
            m_progInfo.m_progress->update_progress( static_cast<long>( currentProgress ),
                                                    static_cast<long>( m_progInfo.m_totalProgress ) );
    }

  public:
    sort_task( const std::pair<RandIter, RandIter>& range, const Pred& pred, progress_info& progInfo, std::size_t es )
        : m_range( range )
        , m_pred( pred )
        , m_progInfo( progInfo )
        , m_elementSize( es ) {}

    tbb::task* execute() {
        std::size_t rangeSize = ( m_range.second - m_range.first );
        if( rangeSize < PARALLEL_CUTOFF ) {
            frantic::sort::sort<RandIter, Pred, ValueType, RandIterTraits>( m_range.first, m_range.second, m_pred,
                                                                            m_elementSize );

            add_progress( rangeSize );

            return NULL;
        } else {
            std::pair<std::size_t, std::size_t> sizes =
                frantic::sort::partition<RandIter, Pred, ValueType, RandIterTraits>( m_range.first, m_range.second,
                                                                                     m_pred, m_elementSize );

            std::size_t pivotElementCount = ( m_range.second - sizes.second ) - ( m_range.first + sizes.first );
            add_progress( pivotElementCount );

            tbb::empty_task& c = *new( allocate_continuation() ) tbb::empty_task;
            c.set_ref_count( 2 );

            tbb::task& t2 = *new( c.allocate_child() ) sort_task<RandIter, Pred, ValueType, RandIterTraits>(
                std::make_pair( m_range.second - sizes.second, m_range.second ), m_pred, m_progInfo, m_elementSize );
            spawn( t2 );

            recycle_as_child_of( c );
            m_range.second = m_range.first + sizes.first;

            return this;
        }
    }
};
} // namespace detail

/**
 * This parallel_sort exists because the tbb::parallel_sort only works on items that can be copied by value. (see
 * std::_Insertion_sort1 for details why)
 */
template <class RandIter, class Pred, class ValueType, class RandIterTraits>
inline void parallel_sort( RandIter begin, RandIter end, Pred cmp, frantic::logging::progress_logger& progress,
                           bool forceSingleThreaded = false ) {
    std::size_t elementSize = RandIterTraits::get_data_size( begin );

    if( forceSingleThreaded || ( end - begin ) < PARALLEL_CUTOFF ) {
        frantic::sort::sort<RandIter, Pred, ValueType, RandIterTraits>( begin, end, cmp, elementSize );
        return;
    }

    detail::progress_info theProgressInfo;

    theProgressInfo.m_mainThreadId = tbb::this_tbb_thread::get_id();
    theProgressInfo.m_currentProgress = 0;
    theProgressInfo.m_totalProgress = ( end - begin );
    theProgressInfo.m_progress = &progress;

    tbb::task& t0 = *new( tbb::task::allocate_root() ) detail::sort_task<RandIter, Pred, ValueType, RandIterTraits>(
        std::make_pair( begin, end ), cmp, theProgressInfo, elementSize );
    tbb::task::spawn_root_and_wait( t0 );

    progress.update_progress( theProgressInfo.m_currentProgress, theProgressInfo.m_totalProgress );
}

template <class RandIter, class Pred>
inline void parallel_sort( RandIter begin, RandIter end, Pred cmp, frantic::logging::progress_logger& progress,
                           bool forceSingleThreaded = false ) {
    parallel_sort<RandIter, Pred, char*, detail::DefaultRandIterTraits>( begin, end, cmp, progress,
                                                                         forceSingleThreaded );
}

/**
 * @overload
 */
template <class RandIter, class Pred>
inline void parallel_sort( RandIter begin, RandIter end, Pred cmp ) {
    frantic::logging::null_progress_logger nullLogger;
    parallel_sort( begin, end, cmp, nullLogger );
}

template <class RandIter, class Pred>
inline void threaded_sort( RandIter begin, RandIter end, Pred cmp, std::size_t numThreads ) {
    if( numThreads <= 1 ) {
        frantic::sort::sort( begin, end, cmp );
    } else {
        parallel_sort( begin, end, cmp );
    }
}

} // namespace sort
} // namespace frantic
