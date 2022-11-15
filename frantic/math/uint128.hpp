// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

namespace frantic {
namespace math {

class uint128 {
    typedef boost::uint32_t word_t;
    typedef boost::uint64_t dword_t;

    static const int wordCount = 4;

    word_t m_data[4];

    uint128( word_t v3, word_t v2, word_t v1, word_t v0 ) { set( v3, v2, v1, v0 ); }

    word_t& operator[]( std::size_t i ) { return m_data[i]; }

    word_t& at( std::size_t i ) { return m_data[i]; }

    word_t at( std::size_t i ) const { return m_data[i]; }

    static dword_t make_dword( word_t hi, word_t lo ) { return dword_t( hi ) << 32 | dword_t( lo ); }

    static word_t get_hi( dword_t val ) { return static_cast<word_t>( ( val >> 32 ) & 0xffffffff ); }

    static word_t get_lo( dword_t val ) { return static_cast<word_t>( val & 0xffffffff ); }

    word_t top_bit() const { return ( at( wordCount - 1 ) & 0x80000000 ) ? 1 : 0; }

    void set( word_t v3, word_t v2, word_t v1, word_t v0 ) {
        m_data[0] = v0;
        m_data[1] = v1;
        m_data[2] = v2;
        m_data[3] = v3;
    }

    void set( const uint128& other ) {
        for( int i = 0; i < wordCount; ++i ) {
            at( i ) = other.at( i );
        }
    }

    void clear() {
        for( int i = 0; i < wordCount; ++i ) {
            at( i ) = 0;
        }
    }

    void complement() {
        for( int i = 0; i < wordCount; ++i ) {
            at( i ) = ~at( i );
        }
    }

    void inc_nocheck() { add_nocheck( uint128( word_t( 1 ) ) ); }

    void neg() {
        complement();
        inc_nocheck();
    }

    dword_t add_nocheck( const uint128& other ) {
        dword_t sum = 0;
        dword_t carry = 0;
        for( int i = 0; i < wordCount; ++i ) {
            sum = dword_t( at( i ) ) + dword_t( other.at( i ) ) + carry;
            carry = get_hi( sum );
            at( i ) = get_lo( sum );
        }
        return carry;
    }

    void add( const uint128& other ) {
        dword_t carry = add_nocheck( other );
        if( carry ) {
            throw std::overflow_error( "uint128::add overflow" );
        }
    }

    void sub( const uint128& other ) {
        if( is_less_than( other ) ) {
            throw std::overflow_error( "uint128::sub overflow" );
        }

        uint128 temp( other );
        temp.neg();
        add_nocheck( temp );
    }

    void shl_word() {
        if( at( wordCount - 1 ) ) {
            throw std::overflow_error( "uint128::shl_word overflow" );
        }
        for( int i = wordCount - 1; i > 0; --i ) {
            at( i ) = at( i - 1 );
        }
        at( 0 ) = 0;
    }

    bool is_zero() const {
        for( int i = 0; i < wordCount; ++i ) {
            if( at( i ) ) {
                return false;
            }
        }
        return true;
    }

    bool is_eq( const uint128& other ) const {
        for( int i = 0; i < wordCount; ++i ) {
            if( at( i ) != other.at( i ) ) {
                return false;
            }
        }
        return true;
    }

    bool is_less_than( const uint128& other ) const {
        for( int i = wordCount - 1; i >= 0; --i ) {
            if( at( i ) > other.at( i ) ) {
                return false;
            } else if( at( i ) < other.at( i ) ) {
                return true;
            }
        }
        return false;
    }

    void mul( const word_t& rhs ) {
        uint128 acc;

        for( int i = wordCount - 1; i >= 0; --i ) {
            acc.shl_word();
            acc.add( uint128( dword_t( at( i ) ) * rhs ) );
        }

        set( acc );
    }

    void mul( const uint128& other ) {
        uint128 acc;

        for( int i = wordCount - 1; i >= 0; --i ) {
            acc.shl_word();
            uint128 temp( *this );
            temp.mul( other.at( i ) );
            acc += temp;
        }

        set( acc );
    }

    void shl() {
        for( int i = wordCount - 1; i >= 1; --i ) {
            dword_t temp = make_dword( at( i ), at( i - 1 ) ) << 1;
            at( i ) = get_hi( temp );
        }
        at( 0 ) <<= 1;
    }

    void shr() {
        for( int i = 0; i < wordCount - 1; ++i ) {
            dword_t temp = make_dword( at( i + 1 ), at( i ) ) >> 1;
            at( i ) = get_lo( temp );
        }
        at( wordCount - 1 ) >>= 1;
    }

    void div( const uint128& divisor, uint128& remainder ) {
        if( divisor.is_zero() ) {
            throw std::runtime_error( "uint128::div divide by zero" );
        }

        uint128 quotient( *this );
        remainder.clear();

        for( int i = 0; i < 128; ++i ) {
            remainder.shl();
            remainder.at( 0 ) |= quotient.top_bit();
            quotient.shl();

            if( remainder >= divisor ) {
                remainder.sub( divisor );
                quotient.at( 0 ) |= 1;
            }
        }
        set( quotient );
    }

    void div( const uint128& other ) {
        uint128 remainder;
        div( other, remainder );
    }

  public:
    uint128() { clear(); }

    uint128( boost::int32_t val ) {
        if( val < 0 ) {
            throw std::runtime_error( "uint128(int32_t) cannot construct from negative number" );
        }
        const word_t unsignedVal = static_cast<word_t>( val );
        set( 0, 0, 0, unsignedVal );
    }

    uint128( boost::int64_t val ) {
        if( val < 0 ) {
            throw std::runtime_error( "uint128(int64_t) cannot construct from negative number" );
        }
        const dword_t unsignedVal = static_cast<dword_t>( val );
        set( 0, 0, get_hi( unsignedVal ), get_lo( unsignedVal ) );
    }

    uint128( word_t val ) { set( 0, 0, 0, val ); }

    uint128( dword_t val ) { set( 0, 0, get_hi( val ), get_lo( val ) ); }

    template <class IntType>
    inline IntType to_integral() const {
        throw std::logic_error( "uint128_t::to_integral() not defined for this type" );
    }

    // Moved out of class definition due to to_integral<boost::uin64_t> needing to be defined first.
    inline std::string to_string() const;

    bool operator==( const uint128& other ) const { return is_eq( other ); }

    bool operator<( const uint128& other ) const { return is_less_than( other ); }

    bool operator!=( const uint128& other ) const { return !( *this == other ); }

    bool operator>( const uint128& other ) const { return other < *this; }

    bool operator<=( const uint128& other ) const { return !( other < *this ); }

    bool operator>=( const uint128& other ) const { return !( *this < other ); }

    uint128& operator+=( const uint128& other ) {
        add( other );
        return *this;
    }

    const uint128 operator+( const uint128& other ) const {
        uint128 temp( *this );
        temp += other;
        return temp;
    }

    uint128& operator*=( const uint128& other ) {
        mul( other );
        return *this;
    }

    const uint128 operator*( const uint128& other ) const {
        uint128 temp( *this );
        temp *= other;
        return temp;
    }

    uint128 operator/=( const uint128& rhs ) {
        div( rhs );
        return *this;
    }

    const uint128 operator/( const uint128& rhs ) const {
        uint128 temp( *this );
        temp /= rhs;
        return temp;
    }
};

template <>
inline boost::int64_t uint128::to_integral<boost::int64_t>() const {
    if( ( *this ) > uint128( std::numeric_limits<boost::int64_t>::max() ) ) {
        throw std::overflow_error( "uint128::convert_to<int64_t> overflow" );
    }
    return static_cast<boost::int64_t>( make_dword( at( 1 ), at( 0 ) ) );
}

std::string uint128::to_string() const {
    uint128 temp( *this );
    uint128 remainder;
    std::string s;

    while( !temp.is_zero() ) {
        temp.div( boost::uint32_t( 10 ), remainder );
        const std::string charString = boost::lexical_cast<std::string>( remainder.to_integral<boost::int64_t>() );
        s.insert( 0, charString );
    }

    if( s.length() == 0 ) {
        s = "0";
    }

    return s;
}

} // namespace math
} // namespace frantic
