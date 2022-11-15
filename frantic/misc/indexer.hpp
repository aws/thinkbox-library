// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * Simple class for multi-dimensional indexing
 *
 * The purpose is to provide a simple means to index a linear array as if it was multi-dimensional.
 *
 */
#pragma once

#include <boost/range.hpp>

#include <functional>
#include <numeric>

namespace frantic {

template <class IntType, size_t Dimension>
class indexer {
  public:
    typedef IntType int_type;

  private:
    IntType m_dimensions[Dimension];

  public:
    template <class InputIterator>
    indexer( InputIterator begin, InputIterator end ) {
        initialize_iterator( begin, end );
    }

    template <class InputRange>
    indexer( InputRange range ) {
        initialize_range( range );
    }

    indexer( const indexer& other ) { initialize_iterator( other.m_dimensions, other.m_dimensions + Dimension ); }

    indexer& operator=( const indexer& other ) {
        initialize_iterator( other.m_dimensions, other.m_dimensions + Dimension );
    }

    template <class InputIterator>
    void initialize_iterator( InputIterator begin, InputIterator end ) {
        size_t currentDimension = 0;
        for( InputIterator it = begin; it != end; ++it ) {
            m_dimensions[currentDimension] = *it;
            ++currentDimension;
        }
    }

    template <class InputRange>
    void initialize_range( InputRange range ) {
        initialize_iterator( boost::const_begin( range ), boost::const_end( range ) );
    }

    size_t num_dimensions() const { return Dimension; }

    void set_dimension( IntType dimension, IntType dimensionSize ) { m_dimensions[dimension] = dimensionSize; }

    IntType get_dimension( IntType dimension ) const { return m_dimensions; }

    IntType address_space() const {
        return std::accumulate( m_dimensions, m_dimensions + Dimension, static_cast<IntType>( 1 ),
                                std::multiplies<IntType>() );
    }

    IntType operator()( IntType x0, IntType x1 = 0, IntType x2 = 0, IntType x3 = 0, IntType x4 = 0,
                        IntType x5 = 0 ) const {
        return address( x0, x1, x2, x3, x4, x5 );
    }

    IntType address( IntType x0, IntType x1 = 0, IntType x2 = 0, IntType x3 = 0, IntType x4 = 0,
                     IntType x5 = 0 ) const {
        IntType indices[6] = {
            x0, x1, x2, x3, x4, x5,
        };
        return address( indices );
    }

    IntType address( IntType* indices ) const { return address_iterator( indices, indices + Dimension ); }

    template <class InputIterator>
    IntType address_iterator( InputIterator begin, InputIterator end ) const {
        IntType spaceSize = address_space();

        IntType currentDimension = 0;

        IntType resultingAddress = IntType( 0 );

        for( InputIterator it = begin; it != end; ++it ) {
            spaceSize /= m_dimensions[currentDimension];
            resultingAddress += *it * spaceSize;
            ++currentDimension;
        }

        return resultingAddress;
    }

    template <class InputRange>
    IntType address_range( InputRange range ) const {
        return address_iterator( boost::const_begin( range ), boost::const_end( range ) );
    }
};

template <class IntType>
indexer<IntType, 2> make_indexer( IntType x0, IntType x1 ) {
    IntType indices[2] = {
        x0,
        x1,
    };
    return indexer<IntType, 2>( indices, indices + 2 );
}

template <class IntType>
indexer<IntType, 3> make_indexer( IntType x0, IntType x1, IntType x2 ) {
    IntType indices[3] = {
        x0,
        x1,
        x2,
    };
    return indexer<IntType, 3>( indices, indices + 3 );
}

template <class IntType>
indexer<IntType, 4> make_indexer( IntType x0, IntType x1, IntType x2, IntType x3 ) {
    IntType indices[2] = {
        x0,
        x1,
        x2,
        x3,
    };
    return indexer<IntType, 2>( indices, indices + 4 );
}

template <class IntType>
indexer<IntType, 5> make_indexer( IntType x0, IntType x1, IntType x2, IntType x3, IntType x4 ) {
    IntType indices[5] = {
        x0, x1, x2, x3, x4,
    };
    return indexer<IntType, 3>( indices, indices + 5 );
}

template <class IntType>
indexer<IntType, 6> make_indexer( IntType x0, IntType x1, IntType x2, IntType x3, IntType x4, IntType x5 ) {
    IntType indices[6] = { x0, x1, x2, x3, x4, x5 };
    return indexer<IntType, 3>( indices, indices + 6 );
}

} // namespace frantic
