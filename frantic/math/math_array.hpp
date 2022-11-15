// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "frantic/math/utils.hpp"
#include <iostream>

namespace frantic {
namespace math {

/**
 * A class for performing simple math operations on
 * an array of numbers.  All operations are element
 * by element operations.
 *
 * @author   Brian McKinnon
 * @since    Apr 25, 2007
 */
template <typename DataType, int ArraySize>
class math_array {
  public:
    typedef DataType DType;

    math_array() {
        for( int i = 0; i < ArraySize; i++ )
            this->m_data[i] = 0;
    }
    ~math_array() {}
    int get_array_size() { return ArraySize; }

    DataType get_data( int i ) const { return m_data[i]; }
    void set_data( int i, DataType d ) { m_data[i] = d; }

    const DataType* get_ptr() const { return m_data; }
    DataType* get_ptr() { return m_data; }

    math_array<DataType, ArraySize> operator-() const;
    math_array<DataType, ArraySize> operator+( const math_array<DataType, ArraySize>& ma ) const;
    math_array<DataType, ArraySize> operator-( const math_array<DataType, ArraySize>& ma ) const;
    math_array<DataType, ArraySize> operator*( const math_array<DataType, ArraySize>& ma ) const;
    math_array<DataType, ArraySize> operator/( const math_array<DataType, ArraySize>& ma ) const;

    math_array<DataType, ArraySize>& operator+=( const math_array<DataType, ArraySize>& ma );
    math_array<DataType, ArraySize>& operator-=( const math_array<DataType, ArraySize>& ma );
    math_array<DataType, ArraySize>& operator*=( const math_array<DataType, ArraySize>& ma );
    math_array<DataType, ArraySize>& operator/=( const math_array<DataType, ArraySize>& ma );

    math_array<DataType, ArraySize> operator+( DataType d ) const;
    math_array<DataType, ArraySize> operator-( DataType d ) const;
    math_array<DataType, ArraySize> operator*( DataType d ) const;
    math_array<DataType, ArraySize> operator/( DataType d ) const;

    math_array<DataType, ArraySize>& operator+=( DataType d );
    math_array<DataType, ArraySize>& operator-=( DataType d );
    math_array<DataType, ArraySize>& operator*=( DataType d );
    math_array<DataType, ArraySize>& operator/=( DataType d );

    math_array<DataType, ArraySize> get_abs() const;
    math_array<DataType, ArraySize> get_sqrt() const;

    bool operator==( const math_array<DataType, ArraySize>& ma );
    bool operator!=( const math_array<DataType, ArraySize>& ma );

    DataType get_sum() const;
    DataType get_mean() const;
    DataType get_max() const;
    DataType get_min() const;

  protected:
    /// Stores the data
    DataType m_data[ArraySize];
};

/**
 * A math_array specialization class designed to deal
 * with the case of a zero length array.
 *
 * @author   Brian McKinnon
 * @since    Apr 26, 2007
 */
template <typename DataType>
class math_array<DataType, 0> {
  public:
    typedef DataType DType;

    math_array() {}
    ~math_array() {}
    int get_array_size() { return 0; }

    DataType get_data( int ) const { return 0; }
    void set_data( int, DataType ) {}

    math_array<DataType, 0> get_abs() const { return ( *this ); }
    math_array<DataType, 0> get_sqrt() const { return ( *this ); }

    DataType get_sum() const { return 0; }
    DataType get_mean() const { return 0; }
    DataType get_max() const { return 0; }
    DataType get_min() const { return 0; }

    math_array<DataType, 0> operator-() const { return ( *this ); }
    math_array<DataType, 0> operator+( const math_array<DataType, 0>& ) const { return ( *this ); }
    math_array<DataType, 0> operator-( const math_array<DataType, 0>& ) const { return ( *this ); }
    math_array<DataType, 0> operator*( const math_array<DataType, 0>& ) const { return ( *this ); }
    math_array<DataType, 0> operator/( const math_array<DataType, 0>& ) const { return ( *this ); }

    math_array<DataType, 0>& operator+=( const math_array<DataType, 0>& ) { return ( *this ); }
    math_array<DataType, 0>& operator-=( const math_array<DataType, 0>& ) { return ( *this ); }
    math_array<DataType, 0>& operator*=( const math_array<DataType, 0>& ) { return ( *this ); }
    math_array<DataType, 0>& operator/=( const math_array<DataType, 0>& ) { return ( *this ); }

    math_array<DataType, 0> operator+( DataType ) const { return ( *this ); }
    math_array<DataType, 0> operator-( DataType ) const { return ( *this ); }
    math_array<DataType, 0> operator*( DataType ) const { return ( *this ); }
    math_array<DataType, 0> operator/( DataType ) const { return ( *this ); }

    math_array<DataType, 0>& operator+=( DataType ) { return ( *this ); }
    math_array<DataType, 0>& operator-=( DataType ) { return ( *this ); }
    math_array<DataType, 0>& operator*=( DataType ) { return ( *this ); }
    math_array<DataType, 0>& operator/=( DataType ) { return ( *this ); }

    bool operator==( const math_array<DataType, 0>& ) { return true; }
    bool operator!=( const math_array<DataType, 0>& ) { return true; }
};

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize> math_array<DataType, ArraySize>::operator-() const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = -this->m_data[i];
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>
math_array<DataType, ArraySize>::operator+( const math_array<DataType, ArraySize>& ma ) const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = this->m_data[i] + ma.m_data[i];
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>
math_array<DataType, ArraySize>::operator-( const math_array<DataType, ArraySize>& ma ) const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = this->m_data[i] - ma.m_data[i];
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>
math_array<DataType, ArraySize>::operator*( const math_array<DataType, ArraySize>& ma ) const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = this->m_data[i] * ma.m_data[i];
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>
math_array<DataType, ArraySize>::operator/( const math_array<DataType, ArraySize>& ma ) const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = this->m_data[i] / ma.m_data[i];
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>&
math_array<DataType, ArraySize>::operator+=( const math_array<DataType, ArraySize>& ma ) {
    for( int i = 0; i < ArraySize; i++ ) {
        this->m_data[i] += ma.m_data[i];
    }

    return ( *this );
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>&
math_array<DataType, ArraySize>::operator-=( const math_array<DataType, ArraySize>& ma ) {
    for( int i = 0; i < ArraySize; i++ ) {
        this->m_data[i] -= ma.m_data[i];
    }

    return ( *this );
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>&
math_array<DataType, ArraySize>::operator*=( const math_array<DataType, ArraySize>& ma ) {
    for( int i = 0; i < ArraySize; i++ ) {
        this->m_data[i] *= ma.m_data[i];
    }

    return ( *this );
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>&
math_array<DataType, ArraySize>::operator/=( const math_array<DataType, ArraySize>& ma ) {
    for( int i = 0; i < ArraySize; i++ ) {
        this->m_data[i] /= ma.m_data[i];
    }

    return ( *this );
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize> math_array<DataType, ArraySize>::operator+( DataType d ) const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = this->m_data[i] + d;
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize> math_array<DataType, ArraySize>::operator-( DataType d ) const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = this->m_data[i] - d;
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize> math_array<DataType, ArraySize>::operator*( DataType d ) const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = this->m_data[i] * d;
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize> math_array<DataType, ArraySize>::operator/( DataType d ) const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = this->m_data[i] / d;
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>& math_array<DataType, ArraySize>::operator+=( DataType d ) {
    for( int i = 0; i < ArraySize; i++ ) {
        this->m_data[i] += d;
    }

    return ( *this );
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>& math_array<DataType, ArraySize>::operator-=( DataType d ) {
    for( int i = 0; i < ArraySize; i++ ) {
        this->m_data[i] -= d;
    }

    return ( *this );
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>& math_array<DataType, ArraySize>::operator*=( DataType d ) {
    for( int i = 0; i < ArraySize; i++ ) {
        this->m_data[i] *= d;
    }

    return ( *this );
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize>& math_array<DataType, ArraySize>::operator/=( DataType d ) {
    for( int i = 0; i < ArraySize; i++ ) {
        this->m_data[i] /= d;
    }

    return ( *this );
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize> math_array<DataType, ArraySize>::get_abs() const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = get_absolute( this->m_data[i] );
    }

    return result;
}

template <typename DataType, int ArraySize>
math_array<DataType, ArraySize> math_array<DataType, ArraySize>::get_sqrt() const {
    math_array<DataType, ArraySize> result;

    for( int i = 0; i < ArraySize; i++ ) {
        result.m_data[i] = (DataType)sqrt( (double)this->m_data[i] );
    }

    return result;
}

template <typename DataType, int ArraySize>
bool math_array<DataType, ArraySize>::operator==( const math_array<DataType, ArraySize>& ma ) {
    for( int i = 0; i < ArraySize; i++ ) {
        if( this->m_data[i] != ma.m_data[i] ) {
            return false;
        }
    }
    return true;
}

template <typename DataType, int ArraySize>
bool math_array<DataType, ArraySize>::operator!=( const math_array<DataType, ArraySize>& ma ) {
    for( int i = 0; i < ArraySize; i++ ) {
        if( this->m_data[i] != ma.m_data[i] ) {
            return true;
        }
    }
    return false;
}

template <typename DataType, int ArraySize>
DataType math_array<DataType, ArraySize>::get_sum() const {
    DataType result = 0;

    for( int i = 0; i < ArraySize; i++ ) {
        result += this->m_data[i];
    }

    return result;
}

template <typename DataType, int ArraySize>
DataType math_array<DataType, ArraySize>::get_mean() const {
    return get_sum() / ArraySize;
}

template <typename DataType, int ArraySize>
DataType math_array<DataType, ArraySize>::get_max() const {
    if( ArraySize < 1 ) {
        return 0;
    } else {
        DataType result = this->m_data[0];

        for( int i = 1; i < ArraySize; i++ ) {
            if( this->m_data[i] > result )
                result = this->m_data[i];
        }
    }
}

template <typename DataType, int ArraySize>
DataType math_array<DataType, ArraySize>::get_min() const {
    if( ArraySize < 1 ) {
        return 0;
    } else {
        DataType result = this->m_data[0];

        for( int i = 1; i < ArraySize; i++ ) {
            if( this->m_data[i] < result )
                result = this->m_data[i];
        }
    }
}

} // namespace math
} // namespace frantic
