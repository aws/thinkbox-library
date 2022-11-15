// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <frantic/particles/streams/culling_particle_istream.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>

#include <set>

namespace frantic {
namespace particles {
namespace streams {

namespace detail {
/**
 * Traits class for storing an ID list in a std::set<>. This option can store a very sparse set most efficiently, but
 * will use lots of memory for dense sets. Testing is O(log n). Negative IDs are supported correctly. \tparam IDType The
 * type of ID to work with.
 */
template <class IDType>
struct std_set_traits {
    typedef IDType id_type;
    typedef std::set<id_type> set_type;

    inline static bool is_in_set( const set_type& theSet, id_type theVal ) {
        return theSet.find( theVal ) != theSet.end();
    }
};

/**
 * Traits class for storing an ID list in a boost::dynamic_bitset<>. This option uses a bitset to store data, so it
 * optimally stores large, dense sets. It is less efficient for very sparse data stored away from ID 0. Membership
 * testing is O(1). Negative IDs are not supported by this set. \tparam IDType The type of ID to work with.
 */
template <class IDType>
struct bitset_traits {
    typedef IDType id_type;
    typedef boost::dynamic_bitset<> set_type;

    inline static bool is_in_set( const set_type& theSet, id_type theVal ) {
        return theSet.size() > static_cast<std::size_t>( theVal ) && theSet.test( static_cast<std::size_t>( theVal ) );
    }
};
} // namespace detail

/**
 * This culling policy culls particles from a stream if the particle has an id that is either in the provided set, or
 * outside of the provided set.
 * \tparam SetTraits A traits class that determines the type of set to store, the type of ID to work with, and abstracts
 * away the details of determining set membership. Must provide: typedef id_type; // A type for representing particle
 * IDs. Must be supported by frantic::channels::channel_cvt_accessor<>. typedef set_type; // The type of set to use for
 * testing IDs against. We store a boost::shared_ptr<> to an instance of this type. static bool is_in_set( const
 * set_type&, id_type theVal ); // Called to determine if an ID is present in the set.
 */
template <class SetTraits = detail::std_set_traits<boost::int32_t>>
class id_culling_policy : public culling_policy_base<id_culling_policy<SetTraits>> {
    boost::shared_ptr<typename SetTraits::set_type> m_idList;
    frantic::channels::channel_static_cast_const_accessor<typename SetTraits::id_type> m_idAccessor;
    bool m_cullInList;

  public:
    typedef boost::tuples::tuple<boost::shared_ptr<typename SetTraits::set_type>, bool> args_type;

    /**
     * Constructs an id culling policy, whose arguments are a set of ids, and a boolean of whether to cull
     * particles with ids that belong to the set, or do not belong to the set
     *
     * @param args tuple of the policy arguments for the policy
     * @param pcm the channel map to use for the stream
     */
    id_culling_policy( const args_type& args, const frantic::channels::channel_map& pcm )
        : m_idList( args.template get<0>() )
        , m_cullInList( args.template get<1>() ) {
        set_channel_map( pcm );
    }

    // Need a splitting constructor to support TBB
    id_culling_policy( const id_culling_policy& lhs, tbb::split )
        : m_idList( lhs.m_idList )
        , m_cullInList( lhs.m_cullInList )
        , m_idAccessor( lhs.m_idAccessor ) {}

    /**
     * Sets the channel map with which the particles will be read from the stream
     *
     * @param pcm a channel map
     */
    void set_channel_map( const frantic::channels::channel_map& pcm ) {
        m_idAccessor.reset();

        if( pcm.has_channel( _T("ID") ) ) {
            frantic::channels::channel_general_accessor rawAccessor = pcm.get_general_accessor( _T("ID") );
            if( rawAccessor.arity() != 1 || !frantic::channels::is_channel_data_type_int( rawAccessor.data_type() ) ) {
                m_idAccessor.reset();
            } else {
                m_idAccessor.reset( rawAccessor );
            }
        }
    }

    /**
     *
     * Function that tests a particle for culling based on its id and whether it is inside or outside
     * the set of ids provided at construction
     *
     * @param particle - a data pointer to the particle
     * @returns true if the particle should be culled
     */
    bool cull( char* particle ) const {
        return m_idAccessor.is_valid() &&
               m_cullInList == SetTraits::is_in_set( *m_idList, m_idAccessor.get( particle ) );
    }
};

typedef culling_particle_istream<id_culling_policy<>> id_culled_particle_istream;

} // namespace streams
} // namespace particles
} // namespace frantic
