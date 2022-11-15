// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file
 * @author Stephen Kiazyk
 *
 * Potentially useful template iterator algorithms that aren't in the c++ standard
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include <tbb/parallel_sort.h>

#include <boost/random/uniform_int.hpp>

namespace frantic {

/**
 * Forms an 'addition pyramid' out of the list, storing the result in-place.
 *  Example: 2 4 7 3 8 3 => 2 6 13 16 24 27
 *
 * @remark Assumes the value type of ForwardIterator is a numeric type that
 *  supports the += operator
 *
 * @param first iterator to the front of the sequence
 * @param last iterator to past the end of the sequence
 */
template <typename ForwardIterator>
void stratify_list( ForwardIterator first, ForwardIterator last ) {
    ForwardIterator next = first;
    ++next;

    while( next != last ) {
        *next += *first;
        ++first;
        ++next;
    }
}

/**
 * Given a set of component labels (just a set of integer labels on each element)
 * group them together into a set of lists, with offset information to denote where
 * each group begins/ends.
 *
 * @param labelsBegin ForwardIterator to the start of the list of component labels
 * @param labelsEnd ForwardIterator to the end of the list of component labels
 * @param outLists OutputIterator to store the sets of components
 * @param outOffsets OutputIterator to store the offset of each component
 *  outComponentLists[outComponentOffsets[c]..outComponentOffsets[c+1]] will contain the faces for component 'c'
 *  (assuming the components are normalized to the range 0..C)
 * @remark the length of outComponentOffsets will be 1 plus the number of components
 */
template <typename ForwardIterator, typename OutputIteratorLists, typename OutputIteratorOffsets>
void group_components( ForwardIterator labelsBegin, ForwardIterator labelsEnd, OutputIteratorLists outLists,
                       OutputIteratorOffsets outOffsets ) {
    typedef typename std::iterator_traits<ForwardIterator>::value_type IntType;

    std::vector<std::pair<IntType, IntType>> listPairs;
    const size_t numElements = std::distance( labelsBegin, labelsEnd );

    if( numElements == 0 ) {
        *outOffsets = 0;
        return;
    }

    listPairs.reserve( numElements );

    IntType currentIndex = 0;
    for( ForwardIterator it = labelsBegin; it != labelsEnd; ++it ) {
        listPairs.push_back( std::make_pair( *it, currentIndex ) );
        ++currentIndex;
    }

#ifndef FRANTIC_DISABLE_THREADS
    tbb::parallel_sort( listPairs.begin(), listPairs.end() );
#else
    std::sort( listPairs.begin(), listPairs.end() );
#endif

    *outOffsets = 0;
    ++outOffsets;

    for( IntType i = 0; i < listPairs.size(); ++i ) {
        *outLists = listPairs[i].second;
        ++outLists;

        if( i > 0 && listPairs[i].first != listPairs[i - 1].first ) {
            *outOffsets = i;
            ++outOffsets;
        }
    }

    *outOffsets = currentIndex;
    ++outOffsets;
}

/**
 * Draw (without replacement) a random sample of size 'sampleSize' from the integer range [0,totalSize-1]
 *
 * @tparam IntType the integer type to draw in the sample
 * @tparam RandomAccessIterator a writable iterator that supports random access
 * @tparam RandomEngine one of the boost random engines
 *
 * @param totalSize the range from which to draw
 * @param outBegin an OutputIterator to the start of the range to write into
 * @param outEnd an OutputIterator to the end of the range to write into
 * @param randomEngine a boost random engine to generate the sample
 * @return a pointer to the end of the written range (will be outEnd unless totalSize is smaller than the output range)
 */
template <typename IntType, typename RandomAccessIterator, typename RandomEngine>
RandomAccessIterator select_index_subset( IntType totalSize, RandomAccessIterator outBegin, RandomAccessIterator outEnd,
                                          RandomEngine randomEngine ) {

    const IntType sampleSize = std::min( IntType( outEnd - outBegin ), totalSize );

    for( IntType i = 0; i < sampleSize; ++i ) {
        *( outBegin + i ) = i;
    }

    boost::uniform_int<size_t> sampleDistro( 0, size_t( sampleSize - 1 ) );

    for( IntType i = sampleSize; i < totalSize; ++i ) {
        boost::uniform_int<size_t> totalDistro( 0, size_t( i ) );

        if( IntType( totalDistro( randomEngine ) ) < sampleSize ) {
            *( outBegin + sampleDistro( randomEngine ) ) = i;
        }
    }

    return outBegin + sampleSize;
}

} // namespace frantic
