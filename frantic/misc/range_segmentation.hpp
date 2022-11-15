// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/misc/algorithm.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <iterator>
#include <vector>

namespace frantic {

/**
 * Provides O(1) access and iteration over a number of subsets which partition a range of integers.
 */
class range_segmentation {
  public:
    class subset_iterator {
      public:
        typedef std::random_access_iterator_tag iterator_category;
        typedef std::size_t value_type;
        typedef ptrdiff_t difference_type;
        typedef size_t* pointer;
        typedef size_t reference;

      private:
        const range_segmentation* m_owner;
        size_t m_subset;
        size_t m_index;

      public:
        // DO NOT USE THE DEFAULT CONSTRUCTOR
        // It exists only to allow std::lower_bound to work on OS X 10.8
        subset_iterator()
            : m_owner( NULL ) {}

        subset_iterator( const range_segmentation* owner, size_t subset, size_t offset = 0 )
            : m_owner( owner )
            , m_subset( subset )
            , m_index( offset ) {}

        subset_iterator( const subset_iterator& other )
            : m_owner( other.m_owner )
            , m_subset( other.m_subset )
            , m_index( other.m_index ) {}

        subset_iterator& operator=( const subset_iterator& other ) {
            m_owner = other.m_owner;
            m_subset = other.m_subset;
            m_index = other.m_index;
            return *this;
        }

        size_t operator*() const { return m_owner->get_subset_face( m_subset, m_index ); }

        subset_iterator& operator++() {
            ++m_index;
            return *this;
        }

        subset_iterator operator++( int ) {
            subset_iterator temp( *this );
            ++( *this );
            return temp;
        }

        subset_iterator& operator--() {
            --m_index;
            return *this;
        }

        subset_iterator operator--( int ) {
            subset_iterator temp( *this );
            --( *this );
            return temp;
        }

        subset_iterator& operator+=( ptrdiff_t offset ) {
            m_index += offset;
            return *this;
        }

        subset_iterator operator+( ptrdiff_t offset ) const {
            subset_iterator temp( *this );
            return temp += offset;
        }

        ptrdiff_t operator-( const subset_iterator& other ) const { return m_index - other.m_index; }

        size_t operator[]( ptrdiff_t offset ) const { return *( ( *this ) + offset ); }

        bool operator==( const subset_iterator& other ) const {
            return m_owner == other.m_owner && m_subset == other.m_subset && m_index == other.m_index;
        }

        bool operator!=( const subset_iterator& other ) const { return !( *this == other ); }

        bool operator<( const subset_iterator& other ) const { return m_index < other.m_index; }
    };

    typedef boost::iterator_range<subset_iterator> subset_range_t;

  private:
    std::vector<size_t> m_subsetElements;
    std::vector<size_t> m_subsetOffsets;
    std::vector<size_t> m_elementAssignments;

  public:
    range_segmentation() {}

    range_segmentation( size_t numFaces ) { reset( numFaces ); }

    range_segmentation( const std::vector<size_t>& elementAssignments, size_t numSubsetsHint = 0 ) {
        assign( elementAssignments, numSubsetsHint );
    }

    template <class InputIterator>
    range_segmentation( InputIterator elementAssignmentsBegin, InputIterator elementAssignmentsEnd,
                        size_t numSubsetsHint = 0 ) {
        assign( elementAssignmentsBegin, elementAssignmentsEnd, numSubsetsHint );
    }

    ~range_segmentation() {}

    void clear() {
        m_elementAssignments.clear();
        m_subsetOffsets.clear();
        m_subsetElements.clear();
    }

    void reset( size_t numFaces ) {
        m_elementAssignments.resize( numFaces );
        m_elementAssignments.assign( numFaces, size_t( -1 ) );
        m_subsetOffsets.resize( 1 );
        m_subsetOffsets[0] = 0;
        m_subsetElements.clear();
    }

    void add_subset( const std::vector<size_t>& list ) { add_subset( list.begin(), list.end() ); }

    template <class InputIterator>
    void add_subset( InputIterator begin, InputIterator end ) {
        const size_t newSubsetId = get_num_subsets();
        size_t subsetSize = 0;

        for( InputIterator it = begin; it != end; ++it ) {
            if( m_elementAssignments.at( *it ) < get_num_subsets() ) {
                throw std::runtime_error( "range_segmentation::add_subset -- face " +
                                          boost::lexical_cast<std::string>( *it ) +
                                          " was already assigned to a subset." );
            }
            m_subsetElements.push_back( *it );
            m_elementAssignments.at( *it ) = newSubsetId;
            ++subsetSize;
        }

        if( subsetSize == 0 ) {
            throw std::runtime_error( "range_segmentation::add_subset -- empty subset" );
        }

        m_subsetOffsets.push_back( m_subsetOffsets.back() + subsetSize );
    }

    template <class InputIterator>
    void merge_from_subrange_segmentation( const range_segmentation& subrangeSeg, InputIterator facesBegin,
                                           InputIterator facesEnd ) {
        const size_t faceListInsertionPoint = m_subsetElements.size();
        const size_t firstNewSubset = get_num_subsets();

        m_subsetElements.resize( faceListInsertionPoint + subrangeSeg.get_assigned_face_count() );
        m_subsetOffsets.resize( m_subsetOffsets.size() + subrangeSeg.get_num_subsets() );

        for( size_t i = 0; i < subrangeSeg.get_num_subsets(); ++i ) {
            m_subsetOffsets[firstNewSubset + 1 + i] =
                faceListInsertionPoint + subrangeSeg.get_subset_offset( i ) + subrangeSeg.get_subset_size( i );
        }

        std::vector<size_t> currentOffsets( subrangeSeg.get_num_subsets(), 0 );

        size_t subrangeFace = 0;
        for( InputIterator it = facesBegin; it != facesEnd; ++it ) {
            const size_t subrangeSubset = subrangeSeg.get_face_subset( subrangeFace );
            ++subrangeFace;

            const size_t newRangeSubset = firstNewSubset + subrangeSubset;
            m_elementAssignments.at( *it ) = newRangeSubset;

            m_subsetElements[get_subset_offset( newRangeSubset ) + currentOffsets[subrangeSubset]] = *it;
            ++currentOffsets[subrangeSubset];
        }

        sort_subsets_range( firstNewSubset, get_num_subsets() );
    }

    void assign( const std::vector<size_t>& assignments, size_t numSubsetsHint = 0 ) {
        assign( assignments.begin(), assignments.end(), numSubsetsHint );
    }

    template <class InputIterator>
    void assign( InputIterator elementAssignmentsBegin, InputIterator elementAssignmentsEnd,
                 size_t numSubsetsHint = 0 ) {
        m_subsetElements.clear();
        m_subsetOffsets.clear();

        if( numSubsetsHint > 0 ) {
            m_subsetOffsets.reserve( numSubsetsHint );
        }

        group_components( elementAssignmentsBegin, elementAssignmentsEnd, std::back_inserter( m_subsetElements ),
                          std::back_inserter( m_subsetOffsets ) );

        m_elementAssignments.resize( m_subsetElements.size() );

        size_t currentSubset = 0;

        for( size_t i = 0; i < m_elementAssignments.size(); ++i ) {
            if( i >= m_subsetOffsets[currentSubset + 1] ) {
                ++currentSubset;
            }

            m_elementAssignments[m_subsetElements[i]] = currentSubset;
        }

        sort_subsets_range( 0, get_num_subsets() );
    }

    size_t get_face_subset( size_t faceId ) const { return m_elementAssignments[faceId]; }

    bool is_face_assigned( size_t faceId ) const { return m_elementAssignments[faceId] < get_num_subsets(); }

    size_t get_num_subsets() const { return m_subsetOffsets.size() - 1; }

    size_t get_subset_size( size_t subsetId ) const {
        return m_subsetOffsets[subsetId + 1] - m_subsetOffsets[subsetId];
    }

    size_t get_subset_offset( size_t subsetId ) const { return m_subsetOffsets[subsetId]; }

    size_t get_assigned_face_count() const { return m_subsetElements.size(); }

    size_t get_assigned_face( size_t i ) const { return m_subsetElements[i]; }

    size_t get_subset_face( size_t subsetId, size_t i ) const {
        const size_t offset = m_subsetOffsets[subsetId];
        return m_subsetElements[offset + i];
    }

    size_t get_face_subset_location( size_t faceId ) const {
        const size_t subsetId = get_face_subset( faceId );
        subset_iterator begin = subset_begin( subsetId );
        subset_iterator end = subset_end( subsetId );
        subset_iterator result = std::lower_bound( begin, end, faceId );
        if( result != end ) {
            return result - begin;
        } else {
            throw std::runtime_error( "range_segmentation::get_face_subset_location: Error, face was not in subset." );
        }
    }

    subset_iterator subset_begin( size_t subsetId ) const { return subset_iterator( this, subsetId ); }

    subset_iterator subset_end( size_t subsetId ) const {
        return subset_iterator( this, subsetId, get_subset_size( subsetId ) );
    }

    subset_range_t subset_range( size_t subsetId ) const {
        return subset_range_t( subset_begin( subsetId ), subset_end( subsetId ) );
    }

  private:
    void sort_subsets_range( size_t startSubset, size_t endSubset ) {
        for( size_t i = startSubset; i < endSubset; ++i ) {
            sort_subset( i );
        }
    }

    void sort_subset( size_t subsetId ) {
        std::sort( m_subsetElements.begin() + m_subsetOffsets[subsetId],
                   m_subsetElements.begin() + m_subsetOffsets[subsetId + 1] );
    }
};

} // namespace frantic
