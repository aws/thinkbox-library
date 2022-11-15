// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

// This particle input stream adaptor can be used to duplicate a channel.
// See the constructor's documentation for usage information.

#include <boost/call_traits.hpp>
#include <frantic/channels/channel_map_adaptor.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include <tbb/parallel_for.h>

namespace frantic {
namespace particles {
namespace streams {

class duplicate_channel_particle_istream : public particle_istream {
  private:
    boost::shared_ptr<particle_istream> m_delegate;

    // channel names
    frantic::tstring m_sourceChannelName;
    frantic::tstring m_destChannelName;

    // channel accessors
    frantic::channels::channel_general_accessor m_sourceAcc;
    frantic::channels::channel_general_accessor m_destAcc;
    size_t m_sourceChannelArity;
    bool m_performChannelCopy;

    // our current channel map
    frantic::channels::channel_map m_channelMap;
    frantic::channels::channel_map_adaptor m_cma;
    frantic::channels::channel_type_convertor_function_t m_sourceDestConverter;

    // the native channel map
    frantic::channels::channel_map m_nativeChannelMap;

    // our default particle
    std::vector<char> m_defaultParticle;

    // scratch space for retrieving particles
    std::vector<char> m_sourceParticleBuffer;

    // upgrading arity from arity 1 to arity > 1. this is off by default.
    bool m_allowLowerToHigherArityCopy;
    size_t m_destChannelArity;

  private:
    void init( const frantic::channels::channel_map& inputChannelMap ) {

        // if a default particle was previously set, copy over the old default structure, otherwise, create a new
        // default particle
        if( m_defaultParticle.size() > 0 ) {
            std::vector<char> newDefaultParticle( inputChannelMap.structure_size() );
            frantic::channels::channel_map_adaptor oldToNewChannelMapAdaptor( inputChannelMap, m_channelMap );
            oldToNewChannelMapAdaptor.copy_structure( &newDefaultParticle[0], &m_defaultParticle[0] );
            m_defaultParticle.swap( newDefaultParticle );
        } else {
            m_defaultParticle.resize( inputChannelMap.structure_size() );
            memset( &m_defaultParticle[0], 0, inputChannelMap.structure_size() );
        }

        // this is our new channel map now. the user demands it.
        m_channelMap = inputChannelMap;

        // make sure our delegate stream is giving us all its channels.
        m_delegate->set_channel_map( m_delegate->get_native_channel_map() );
        frantic::channels::channel_map sourceCm = m_delegate->get_channel_map();

        // create native channel map. this will get modified if channels are being copied (below).
        m_nativeChannelMap = m_delegate->get_native_channel_map();

        // decide if channel needs to be copied (source and destination exist).
        m_performChannelCopy = false;
        if( sourceCm.has_channel( m_sourceChannelName ) ) {
            if( m_channelMap.has_channel( m_destChannelName ) ) {
                // the dest channel has been requested. make sure we can actually do the copy, and set up everything for
                // retrieving particles.
                m_performChannelCopy = true;

                m_sourceAcc = sourceCm.get_general_accessor( m_sourceChannelName );
                m_destAcc = m_channelMap.get_general_accessor( m_destChannelName );

                m_sourceChannelArity = m_sourceAcc.arity();
                m_destChannelArity = m_destAcc.arity();

                // can only copy channels of the same arity (or when upgrade is requested)
                if( m_destChannelArity == m_sourceChannelArity ||
                    ( m_allowLowerToHigherArityCopy && m_destChannelArity > m_sourceChannelArity ) ) {

                    // create the converter for data types
                    m_sourceDestConverter = frantic::channels::get_channel_type_convertor_function(
                        m_sourceAcc.data_type(), m_destAcc.data_type(), m_sourceChannelName );

                    // if we've requested the new channel our native channel map, is our delegate's native map, PLUS our
                    // new copied channel.
                    if( m_nativeChannelMap.has_channel( m_destChannelName ) )
                        m_nativeChannelMap.delete_channel(
                            m_destChannelName ); // overwritten channel no longer available to the native map
                    m_nativeChannelMap.append_channel(
                        m_destChannelName, m_destAcc.arity(),
                        m_destAcc.data_type() ); // newly copied channel now available to the native map

                    // as to not mess up the adaptor, get any existing channel with this name out of the way (it's sort
                    // of tricking the adaptor to throw away the channel, but it works great).
                    if( sourceCm.has_channel( m_destChannelName ) )
                        sourceCm.rename_channel( m_destChannelName, m_destChannelName + _T("_OVERWRITTEN") );
                } else {
                    throw std::runtime_error( "Cannot copy channel \"" +
                                              frantic::strings::to_string( m_sourceChannelName ) + "\" to \"" +
                                              frantic::strings::to_string( m_destChannelName ) +
                                              "\" because they have incompatible arities." );
                }
            } else {
                // the dest channel has not actually been requested. still, we need to add it as being available to our
                // native channel map.
                if( !m_nativeChannelMap.has_channel( m_destChannelName ) ) {
                    frantic::channels::channel_general_accessor sourceAcc =
                        sourceCm.get_general_accessor( m_sourceChannelName );
                    m_nativeChannelMap.append_channel( m_destChannelName, sourceAcc.arity(), sourceAcc.data_type() );
                }
            }
        }

        // make adaptor
        // sourceCm should be the stream's native channel map (minus existing destination channel, if it was in the way)
        // m_channelMap is the user's channel map for this stream.
        m_cma.set( m_channelMap, sourceCm );

        // allocate incoming particle buffer
        m_sourceParticleBuffer.resize( sourceCm.structure_size() );
    }

  public:
    /**
     * Creates a stream that duplicates sourceChannelName as destChannelName.
     *
     * @param pin This incoming particle stream
     * @param sourceChannelName The source channel to be duplicated
     * @param destChannelName The channel that the source will be copied into. Note: This stream does not create
     * destChannelName, so if destChannelName does not exist in the source's native channel map, the user will need to
     * call set_channel_map (with a map that has the duplicate channel) on this stream to get that channel when
     * retrieving the particles.
     * @param allowLowerToHigherArityCopy Normally only channels with be duplicated if the source and destination
     * channels have the same arity. If this parameter is set to true, it will ignore arity mismatches when the source
     * channel has lower arity than the destination channel, and just copy the portion needed, and zero the rest.
     */
    duplicate_channel_particle_istream( boost::shared_ptr<particle_istream> pin,
                                        const frantic::tstring& sourceChannelName,
                                        const frantic::tstring& destChannelName,
                                        bool allowLowerToHigherArityCopy = false )
        : m_delegate( pin )
        , m_sourceChannelName( sourceChannelName )
        , m_destChannelName( destChannelName ) {
        m_allowLowerToHigherArityCopy = allowLowerToHigherArityCopy;
        init( m_delegate->get_channel_map() );
    }

    virtual ~duplicate_channel_particle_istream() {}

    // forwards to delegate
    void close() { m_delegate->close(); }
    boost::int64_t particle_count() const { return m_delegate->particle_count(); }
    boost::int64_t particle_index() const { return m_delegate->particle_index(); }
    boost::int64_t particle_count_left() const { return m_delegate->particle_count_left(); }
    boost::int64_t particle_progress_count() const { return m_delegate->particle_progress_count(); }
    boost::int64_t particle_progress_index() const { return m_delegate->particle_progress_index(); }
    boost::int64_t particle_count_guess() const { return m_delegate->particle_count_guess(); }
    frantic::tstring name() const { return m_delegate->name(); }

    std::size_t particle_size() const { return m_channelMap.structure_size(); }

    void set_channel_map( const frantic::channels::channel_map& particleChannelMap ) { init( particleChannelMap ); }

    void set_default_particle( char* rawParticleBuffer ) {
        memcpy( &m_defaultParticle[0], rawParticleBuffer, m_channelMap.structure_size() );
    }

    const frantic::channels::channel_map& get_channel_map() const { return m_channelMap; }

    const frantic::channels::channel_map& get_native_channel_map() const { return m_nativeChannelMap; }

    bool get_particle( char* rawDestParticleBuffer ) {
        char* rawSourceParticleBuffer = &m_sourceParticleBuffer[0];

        // get the delegate's particle (NOTE: it comes in ALWAYS with the delegate's native channel map)
        if( !m_delegate->get_particle( rawSourceParticleBuffer ) )
            return false;

        // convert particle to our current channel map
        m_cma.copy_structure( rawDestParticleBuffer, rawSourceParticleBuffer, &m_defaultParticle[0] );

        // copy the requested channel value
        if( m_performChannelCopy ) {
            if( m_destChannelArity >
                m_sourceChannelArity ) { // not normally used, but if the arity is being upgraded, we
                                         // have to assign it in a loop. kind of a hack for bobo.
                for( std::size_t i = 0; i < m_destChannelArity; ++i ) {
                    // in this case, it only really makes sense if the source arity is one. this function will cause
                    // random memory problems if the arity of the destination does not divide evenly by the arity of the
                    // source.
                    m_sourceDestConverter( m_destAcc.get_channel_data_pointer( rawDestParticleBuffer ) +
                                               i * sizeof_channel_data_type( m_destAcc.data_type() ),
                                           m_sourceAcc.get_channel_data_pointer( rawSourceParticleBuffer ),
                                           m_sourceChannelArity );
                }
            } else {
                m_sourceDestConverter( m_destAcc.get_channel_data_pointer( rawDestParticleBuffer ),
                                       m_sourceAcc.get_channel_data_pointer( rawSourceParticleBuffer ),
                                       m_sourceChannelArity );
            }
        }

        return true;
    }

    bool get_particles( char* buffer, std::size_t& numParticles ) {
        // TODO: threaded?
        size_t offset = 0;
        size_t particleSize = m_channelMap.structure_size();
        while( offset < numParticles * particleSize ) {
            if( !get_particle( buffer + offset ) ) {
                numParticles = offset / particleSize; // am i supposed to modify numParticles?
                return false;
            }
            offset += particleSize;
        }
        return true;
    }
};

} // namespace streams
} // namespace particles
} // namespace frantic
