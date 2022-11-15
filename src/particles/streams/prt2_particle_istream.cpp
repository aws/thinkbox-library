// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/streams/prt2_particle_istream.hpp>

using namespace std;
using namespace frantic;
using namespace frantic::prtfile;
using namespace frantic::particles::streams;

prt2_particle_istream::prt2_particle_istream( const frantic::tstring& file,
                                              frantic::channels::data_type_t positionTypeHint,
                                              const frantic::tstring& chunkName )
    : m_prt2( new prt2_reader() ) {
    m_prt2->open( file );

    set_channel_map_with_position_type( m_prt2->get_channel_map(), positionTypeHint );

    m_chunksList.push_back( chunk_entry::all( chunkName ) );

    init();
}

prt2_particle_istream::prt2_particle_istream( const frantic::tstring& file,
                                              const frantic::channels::channel_map& particleChannelMap,
                                              const frantic::tstring& chunkName )
    : m_prt2( new prt2_reader() ) {
    m_prt2->open( file );

    set_channel_map( particleChannelMap );

    m_chunksList.push_back( chunk_entry::all( chunkName ) );

    init();
}

prt2_particle_istream::prt2_particle_istream( const frantic::tstring& file, const chunk_entry& chunk,
                                              frantic::channels::data_type_t positionTypeHint )
    : m_prt2( new prt2_reader() ) {
    m_prt2->open( file );

    set_channel_map_with_position_type( m_prt2->get_channel_map(), positionTypeHint );

    m_chunksList.push_back( chunk );

    init();
}

prt2_particle_istream::prt2_particle_istream( const frantic::tstring& file,
                                              const frantic::channels::channel_map& particleChannelMap,
                                              const chunk_entry& chunk )
    : m_prt2( new prt2_reader() ) {
    m_prt2->open( file );

    set_channel_map( particleChannelMap );

    m_chunksList.push_back( chunk );

    init();
}

prt2_particle_istream::prt2_particle_istream( const frantic::tstring& file, const chunk_list_t& chunks,
                                              frantic::channels::data_type_t positionTypeHint )
    : m_chunksList( chunks )
    , m_prt2( new prt2_reader() ) {
    m_prt2->open( file );

    set_channel_map_with_position_type( m_prt2->get_channel_map(), positionTypeHint );

    init();
}

prt2_particle_istream::prt2_particle_istream( const frantic::tstring& file,
                                              const frantic::channels::channel_map& particleChannelMap,
                                              const chunk_list_t& chunks )
    : m_chunksList( chunks )
    , m_prt2( new prt2_reader() ) {
    m_prt2->open( file );

    set_channel_map( particleChannelMap );

    init();
}
prt2_particle_istream::prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file,
                                              frantic::channels::data_type_t positionTypeHint,
                                              const frantic::tstring& chunkName )
    : m_prt2( file ) {
    set_channel_map_with_position_type( m_prt2->get_channel_map(), positionTypeHint );

    m_chunksList.push_back( chunk_entry::all( chunkName ) );

    init();
}

prt2_particle_istream::prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file,
                                              const frantic::channels::channel_map& particleChannelMap,
                                              const frantic::tstring& chunkName )
    : m_prt2( file ) {
    set_channel_map( particleChannelMap );

    m_chunksList.push_back( chunk_entry::all( chunkName ) );

    init();
}

prt2_particle_istream::prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file,
                                              const chunk_entry& chunk,
                                              frantic::channels::data_type_t positionTypeHint )
    : m_prt2( file ) {
    set_channel_map_with_position_type( m_prt2->get_channel_map(), positionTypeHint );

    m_chunksList.push_back( chunk );

    init();
}

prt2_particle_istream::prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file,
                                              const frantic::channels::channel_map& particleChannelMap,
                                              const chunk_entry& chunk )
    : m_prt2( file ) {
    set_channel_map( particleChannelMap );

    m_chunksList.push_back( chunk );

    init();
}

prt2_particle_istream::prt2_particle_istream( const boost::shared_ptr<frantic::prtfile::prt2_reader>& file,
                                              const chunk_list_t& chunks,
                                              frantic::channels::data_type_t positionTypeHint )
    : m_chunksList( chunks )
    , m_prt2( file ) {
    set_channel_map_with_position_type( m_prt2->get_channel_map(), positionTypeHint );

    init();
}

void prt2_particle_istream::init() {
    m_particleCount = 0;

    m_nextChunkIndex = 0;
    m_nextParticleIndexInChunk = 0;
    m_chunkParticleCount = 0;

    m_currentChunkIndex = 0;
    m_currentChunkEntryIndex = 0;

    m_currentParticleIndex = -1;

    // Go through the chunks to read, validate it, and calculate the expected particle count
    for( chunk_list_t::iterator iter = m_chunksList.begin(); iter != m_chunksList.end(); ) {
        chunk_entry& entry = *iter;
        for( std::vector<std::pair<boost::int64_t, boost::int64_t>>::iterator entryIter = entry.m_chunkIndices.begin();
             entryIter != entry.m_chunkIndices.end(); ) {
            std::pair<boost::int64_t, boost::int64_t>& entryRange = *entryIter;
            if( entryRange.first < 0 )
                entryRange.first = 0;
            if( entryRange.second < 0 )
                entryRange.second = 0;
            if( entryRange.first > entryRange.second )
                std::swap( entryRange.first, entryRange.second );

            const boost::int64_t totalCount = m_prt2->get_particle_chunk_count( entry.m_chunkName );
            if( entryRange.first >= totalCount ) {
                entryIter = entry.m_chunkIndices.erase( entryIter );
                continue;
            }
            if( entryRange.second >= totalCount ) {
                entryRange.second = totalCount - 1;
            }

            for( boost::int64_t i = entryRange.first; i <= entryRange.second; ++i ) {
                m_particleCount += m_prt2->get_particle_chunk_particle_count( i, entry.m_chunkName );
            }
            ++entryIter;
        }

        if( entry.m_chunkIndices.empty() ) {
            iter = m_chunksList.erase( iter );
        } else {
            ++iter;
        }
    }

    if( m_currentChunkIndex < m_chunksList.size() ) {
        chunk_entry& entry = m_chunksList[m_currentChunkIndex];
        if( m_currentChunkEntryIndex < entry.m_chunkIndices.size() ) {
            m_nextChunkIndex = entry.m_chunkIndices[m_currentChunkEntryIndex].first;
        }
    }

    m_streamName = m_prt2->get_stream_name();
}

prt2_particle_istream::~prt2_particle_istream() { close(); }

void prt2_particle_istream::close() { m_prt2.reset(); }

void prt2_particle_istream::set_channel_map_with_position_type( const frantic::channels::channel_map& pcm,
                                                                frantic::channels::data_type_t positionTypeHint ) {
    const channels::channel& posChan = pcm[_T("Position")];
    if( positionTypeHint == channels::data_type_float64 ) {
        if( m_prt2->is_position_offset() ) {
            // When float64 is requested and the file is Position+Offset, make Position float64 so we can pass along the
            // precision
            channels::channel_map newCm;
            for( size_t i = 0; i != pcm.channel_count(); ++i ) {
                const channels::channel& chan = pcm[i];
                if( chan.name() == _T("Position") ) {
                    newCm.define_channel( chan.name(), chan.arity(), channels::data_type_float64 );
                } else {
                    newCm.define_channel( chan.name(), chan.arity(), chan.data_type() );
                }
            }
            newCm.end_channel_definition();
            set_channel_map( newCm );
        } else {
            // Simply pass through the source precision
            set_channel_map( pcm );
        }
    } else {
        if( posChan.data_type() != channels::data_type_float32 ) {
            // When float32 or invalid is requested and the input is different, downgrade to float32
            channels::channel_map newCm;
            for( size_t i = 0; i != pcm.channel_count(); ++i ) {
                const channels::channel& chan = pcm[i];
                if( chan.name() == _T("Position") ) {
                    newCm.define_channel( chan.name(), chan.arity(), channels::data_type_float32 );
                } else {
                    newCm.define_channel( chan.name(), chan.arity(), chan.data_type() );
                }
            }
            newCm.end_channel_definition();
            set_channel_map( newCm );
        } else {
            // Simply pass through the source precision
            set_channel_map( pcm );
        }
    }
}

void prt2_particle_istream::set_channel_map( const channels::channel_map& particleChannelMap ) {
    boost::scoped_array<char> newDefaultParticle( new char[particleChannelMap.structure_size()] );
    particleChannelMap.construct_structure( newDefaultParticle.get() );

    // If we previously had a default particle set, adapt its channels to the new layout.
    if( m_defaultParticle ) {
        frantic::channels::channel_map_adaptor defaultAdaptor( particleChannelMap, m_particleChannelMap );
        defaultAdaptor.copy_structure( newDefaultParticle.get(), m_defaultParticle.get() );
    }

    m_defaultParticle.swap( newDefaultParticle );

    // Set the map and the adaptor
    m_particleChannelMap = particleChannelMap;
    m_pcmAdaptor.set( m_particleChannelMap, m_prt2->get_channel_map() );

    // In float64 Position+Offset mode, set up for float64 offseting
    m_isFloat64PositionOffset = ( m_particleChannelMap[_T("Position")].data_type() == channels::data_type_float64 &&
                                  m_prt2->is_position_offset() );
    if( m_isFloat64PositionOffset ) {
        m_pos64Accessor = m_particleChannelMap.get_accessor<graphics::vector3fd>( _T("Position") );
        m_pos64Offset.set( 0.0 );
    }
}

void prt2_particle_istream::set_default_particle( char* rawParticleBuffer ) {
    m_particleChannelMap.copy_structure( m_defaultParticle.get(), rawParticleBuffer );
}

bool prt2_particle_istream::chunk_entry_advance() {
    chunk_entry* currentEntry = &m_chunksList[m_currentChunkIndex];
    while( m_nextChunkIndex > currentEntry->m_chunkIndices[m_currentChunkEntryIndex].second ) {
        ++m_currentChunkEntryIndex;
        while( m_currentChunkEntryIndex >= currentEntry->m_chunkIndices.size() ) {
            ++m_currentChunkIndex;
            m_nextChunkIndex = 0;
            m_currentChunkEntryIndex = 0;
            if( m_currentChunkIndex >= m_chunksList.size() )
                return false;
            currentEntry = &m_chunksList[m_currentChunkIndex];
        }
        m_nextChunkIndex = currentEntry->m_chunkIndices[m_currentChunkEntryIndex].first;
    }
    return true;
}

bool prt2_particle_istream::chunk_advance() {
    // A loop to allow for the possibility of zero-sized chunks
    while( m_nextParticleIndexInChunk >= m_chunkParticleCount ) {
        if( m_currentChunkIndex >= m_chunksList.size() ) {
            // Trim the chunk buffer so that its memory usage is minimal.
            // (If this was part of a concatenated_particle_istream, this stream will take up extra memory even if we
            // are done reading from it. The chunk buffer will account for a significant portion of it if we don't do
            // something about it) This is somewhat of a hack since C++ vectors don't directly support trimming.
            std::vector<char>().swap( m_chunkBuffer );
            return false;
        } else {
            while( m_nextChunkIndex >=
                   m_prt2->get_particle_chunk_count( m_chunksList[m_currentChunkIndex].m_chunkName ) ) {
                ++m_nextChunkIndex;
                if( !chunk_entry_advance() ) {
                    std::vector<char>().swap( m_chunkBuffer );
                    return false;
                }
            }
            m_prt2->read_particle_chunk( m_nextChunkIndex, m_chunkBuffer,
                                         m_isFloat64PositionOffset ? &m_pos64Offset : NULL,
                                         m_chunksList[m_currentChunkIndex].m_chunkName );
            ++m_nextChunkIndex;
            chunk_entry_advance();
        }
        m_nextParticleIndexInChunk = 0;
        m_chunkParticleCount = m_chunkBuffer.size() / m_prt2->get_channel_map().structure_size();
    }
    return true;
}

bool prt2_particle_istream::get_particle( char* rawParticleBuffer ) {
    size_t countRead = 1;
    return get_particles( rawParticleBuffer, countRead );
}

bool prt2_particle_istream::get_particles( char* particleBuffer, std::size_t& numParticles ) {
    std::size_t totalRead = 0;
    char* currentParticleBuffer = particleBuffer;
    const size_t chunkParticleSize = m_prt2->get_channel_map().structure_size();
    const size_t outputParticleSize = m_particleChannelMap.structure_size();

    while( totalRead < numParticles && chunk_advance() ) {
        std::size_t leftToRead =
            std::min( numParticles - totalRead, (size_t)m_chunkParticleCount - (size_t)m_nextParticleIndexInChunk );

        const char* chunkParticleBuffer = &m_chunkBuffer[m_nextParticleIndexInChunk * chunkParticleSize];
        if( m_pcmAdaptor.is_identity() && outputParticleSize == chunkParticleSize ) {
            memcpy( currentParticleBuffer, chunkParticleBuffer, leftToRead * outputParticleSize );
        } else if( m_isFloat64PositionOffset ) {
            // Apply the offsets in float64 instead of float32 to keep the precision
            graphics::vector3fd offset( m_pos64Offset.x, m_pos64Offset.y, m_pos64Offset.z );
            for( size_t i = 0; i != leftToRead; ++i ) {
                m_pcmAdaptor.copy_structure( currentParticleBuffer, chunkParticleBuffer, m_defaultParticle.get() );
                m_pos64Accessor( currentParticleBuffer ) += offset;
                chunkParticleBuffer += chunkParticleSize;
                currentParticleBuffer += outputParticleSize;
            }
        } else {
            for( size_t i = 0; i != leftToRead; ++i ) {
                m_pcmAdaptor.copy_structure( currentParticleBuffer, chunkParticleBuffer, m_defaultParticle.get() );
                chunkParticleBuffer += chunkParticleSize;
                currentParticleBuffer += outputParticleSize;
            }
        }

        m_nextParticleIndexInChunk += leftToRead;
        m_currentParticleIndex += leftToRead;
        totalRead += leftToRead;
        currentParticleBuffer = &particleBuffer[outputParticleSize * totalRead];
    }

    numParticles = totalRead;
    return ( totalRead > 0 );
}
