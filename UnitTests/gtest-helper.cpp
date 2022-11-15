// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "UnitTests/gtest-helper.h"

using namespace frantic::channels;
using namespace std;
using frantic::particles::particle_array;
using frantic::particles::particle_file_metadata;
using frantic::particles::particle_file_stream_factory_object;
using frantic::particles::streams::particle_istream;

// Gets name, type, arity, and channel general accessors for each channel in cm
void getAllChannelDefinition( const channel_map& cm, std::vector<frantic::tstring>& outChannelNames,
                              std::vector<data_type_t>& outChannelTypes, std::vector<size_t>& outChannelArities,
                              std::vector<channel_general_accessor>& outAccessors ) {
    data_type_t cm_type;
    size_t cm_arity;
    for( int i = 0; i < cm.channel_count(); i++ ) {
        frantic::tstring channelName;
        cm.get_channel_definition( i, channelName, cm_type, cm_arity );
        outChannelNames.push_back( channelName );
        outChannelTypes.push_back( cm_type );
        outChannelArities.push_back( cm_arity );
        channel_general_accessor accessor = cm.get_general_accessor( channelName );
        outAccessors.push_back( accessor );
    }
}

namespace testing {
namespace internal {
AssertionResult cmpChannelMap( const char* /*expected_cm_expr*/, const char* /*actual_cm_expr*/,
                               const channel_map& expected_cm, const channel_map& actual_cm ) {
    data_type_t expected_type;
    size_t expected_arity;
    data_type_t actual_type;
    size_t actual_arity;
    frantic::tstring channelName;
    basic_stringstream<frantic::tchar> expected_cm_ss;
    expected_cm_ss << "[";
    basic_stringstream<frantic::tchar> actual_cm_ss;
    actual_cm_ss << "[";
    for( int i = 0; i < expected_cm.channel_count(); i++ ) {
        expected_cm.get_channel_definition( i, channelName, expected_type, expected_arity );
        expected_cm_ss << channelName << " ( " << channel_data_type_str( expected_type ) << "[" << expected_arity
                       << "] )";
        if( i < expected_cm.channel_count() - 1 )
            expected_cm_ss << ", ";
    }
    expected_cm_ss << "]";
    for( int i = 0; i < actual_cm.channel_count(); i++ ) {
        actual_cm.get_channel_definition( i, channelName, actual_type, actual_arity );
        actual_cm_ss << channelName << " ( " << channel_data_type_str( actual_type ) << "[" << actual_arity << "] )";
        if( i < actual_cm.channel_count() - 1 )
            actual_cm_ss << ", ";
    }
    actual_cm_ss << "]";

    bool success = true;
    if( expected_cm.channel_count() != actual_cm.channel_count() ) {
        success = false;
    }
    for( int i = 0; i < expected_cm.channel_count(); i++ ) {
        expected_cm.get_channel_definition( i, channelName, expected_type, expected_arity );
        if( !actual_cm.has_channel( channelName ) ) {
            success = false;
        } else {
            actual_cm.get_channel_definition( channelName, actual_type, actual_arity );
            if( expected_type != actual_type ) {
                success = false;
            }
        }
    }
    if( success )
        return AssertionSuccess();
    return AssertionFailure() << "Expected: " << expected_cm_ss.str() << "\nActual: " << actual_cm_ss.str() << "\n";
}

AssertionResult cmpChannelData( const char* /*expected_expr*/, const char* /*actual_expr*/,
                                const char* /*channel_name_expr*/, const char* /*data_type_expr*/,
                                const char* /*arity_expr*/, const char* expected, const char* actual,
                                const frantic::tstring& channel_name, data_type_t data_type, size_t arity ) {
    if( memcmp( expected, actual, arity * sizeof_channel_data_type( data_type ) ) == 0 ) {
        return AssertionSuccess();
    }
    ::std::stringstream expected_ss;
    ::std::stringstream actual_ss;
    expected_ss << "[";
    actual_ss << "[";
    channel_data_type_print( expected_ss, ", ", arity, data_type, expected );
    channel_data_type_print( actual_ss, ", ", arity, data_type, actual );
    expected_ss << "]";
    actual_ss << "]";
    return AssertionFailure() << channel_name << ": " << channel_data_type_str( data_type ) << "[" << arity << "]\n"
                              << "Expected: " << StringStreamToString( &expected_ss ) << "\n"
                              << "Actual: " << StringStreamToString( &actual_ss ) << "\n";
}

AssertionResult const cmpChannelData( const char* expected, const char* actual, const frantic::tstring& channel_name,
                                      data_type_t data_type, size_t arity ) {
    return cmpChannelData( "expected", "actual", "channel name", "data_type_expr", "arity_expr", expected, actual,
                           channel_name, data_type, arity );
}

AssertionResult cmpParticleArray( const char* expected_expr, const char* actual_expr, const particle_array& expected,
                                  const particle_array& actual ) {
    // Put test particles in an array with a matching channel map to ref particles
    particle_array actual_temp = particle_array( expected.get_channel_map() );
    actual_temp.copy_particles_from( actual );
    // Confirm that the arrays have exactly the same channel map, otherwise comparison of data will fail.
    if( expected.get_channel_map() != actual_temp.get_channel_map() ) {
        return AssertionFailure() << "Attempted to set " << actual_expr << " particle array channel map to match "
                                  << expected_expr << " particle array channel map, but conversion failed.\n";
    }

    std::stringstream failMsg;
    // Test particle count
    bool sameParticleCount = ( expected.particle_count() == actual_temp.particle_count() );
    if( !sameParticleCount )
        failMsg << expected_expr << " and " << actual_expr
                << " have different particle counts.\nExpected: " << expected.particle_count()
                << "\nActual: " << actual_temp.particle_count() << "\n";

    // Get channel data
    channel_map expected_cm = expected.get_channel_map();
    std::vector<frantic::tstring> channelNames;
    std::vector<data_type_t> channelTypes;
    std::vector<size_t> channelArities;
    std::vector<channel_general_accessor> accessors;
    getAllChannelDefinition( expected_cm, channelNames, channelTypes, channelArities, accessors );

    // Compare particle arrays
    for( int i = 0; i < min( expected.particle_count(), actual_temp.particle_count() ); i++ ) {
        for( int j = 0; j < expected_cm.channel_count(); j++ ) {
            const char* expected_data = accessors[j].get_channel_data_pointer( expected[i] );
            const char* actual_data = accessors[j].get_channel_data_pointer( actual_temp[i] );
            AssertionResult channelDataCompare =
                cmpChannelData( expected_data, actual_data, channelNames[j], channelTypes[j], channelArities[j] );
            if( !channelDataCompare ) {
                if( !sameParticleCount )
                    failMsg << "\n";
                return AssertionFailure() << StringStreamToString( &failMsg ) << expected_expr << " and " << actual_expr
                                          << " have different particle data at index " << i << ".\n"
                                          << channelDataCompare.message();
            }
        }
    }
    if( !sameParticleCount )
        return AssertionFailure() << StringStreamToString( &failMsg );
    return AssertionSuccess();
}

AssertionResult cmpParticleData( const char* expected_expr, const char* actual_expr,
                                 boost::shared_ptr<particle_istream> expected,
                                 boost::shared_ptr<particle_istream> actual ) {
    const channel_map& expected_cm = expected->get_channel_map();
    const channel_map& actual_cm = actual->get_channel_map();
    std::stringstream failMsg;

    // Test that channel maps are equivalent
    AssertionResult sameChannels = cmpChannelMap( expected_expr, actual_expr, expected_cm, actual_cm );
    if( !sameChannels ) {
        failMsg << expected_expr << " and " << actual_expr << " have different particle channel maps.\n"
                << sameChannels.message() << "\n";

        // Test particle count. If not tested here, it will be by cmpParticleArray.
        bool sameParticleCount = ( expected->particle_count() == actual->particle_count() );
        if( !sameParticleCount )
            failMsg << expected_expr << " and " << actual_expr
                    << " have different particle counts.\nExpected: " << expected->particle_count()
                    << "\nActual: " << actual->particle_count() << "\n";

        return AssertionFailure() << StringStreamToString( &failMsg );
    }

    // Stream particles
    frantic::particles::particle_array expected_particles = particle_array( expected_cm );
    expected_particles.insert_particles( expected );
    particle_array actual_particles = particle_array( actual_cm );
    actual_particles.insert_particles( actual );
    return cmpParticleArray( expected_expr, actual_expr, expected_particles, actual_particles );
}

AssertionResult cmpPropertyMap( const char* expected_expr, const char* actual_expr,
                                const frantic::channels::property_map& expected,
                                const frantic::channels::property_map& actual ) {
    channel_map expected_cm = expected.get_channel_map();
    channel_map actual_cm = actual.get_channel_map();
    AssertionResult cmResult = cmpChannelMap( expected_expr, actual_expr, expected_cm, actual_cm );
    if( !cmResult )
        return AssertionFailure() << expected_expr << " and " << actual_expr
                                  << " have different metadata property fields.\n"
                                  << cmResult.message();

    std::stringstream failMsg;
    std::vector<frantic::tstring> propertyNames;
    std::vector<data_type_t> propertyTypes;
    std::vector<size_t> propertyArities;
    std::vector<channel_general_accessor> accessors;
    getAllChannelDefinition( expected_cm, propertyNames, propertyTypes, propertyArities, accessors );

    bool success = true;
    for( int i = 0; i < expected_cm.channel_count(); i++ ) {
        const char* expected_data = accessors[i].get_channel_data_pointer( expected.get_raw_buffer() );
        std::stringstream expected_ss;
        frantic::channels::channel_data_type_print( expected_ss, ", ", propertyArities[i], propertyTypes[i],
                                                    expected_data );
        const char* actual_data = accessors[i].get_channel_data_pointer( actual.get_raw_buffer() );
        std::stringstream actual_ss;
        frantic::channels::channel_data_type_print( actual_ss, ", ", propertyArities[i], propertyTypes[i],
                                                    actual_data );
        AssertionResult dataResult =
            cmpChannelData( expected_data, actual_data, propertyNames[i], propertyTypes[i], propertyArities[i] );
        if( !dataResult ) {
            success = false;
            failMsg << expected_expr << " and " << actual_expr << " have different metadata property values.\n"
                    << dataResult.message();
        }
    }
    if( success )
        return AssertionSuccess();
    return AssertionFailure() << StringStreamToString( &failMsg );
}

AssertionResult cmpMetadata( const char* expected_expr, const char* actual_expr, const particle_file_metadata& expected,
                             const particle_file_metadata& actual ) {
    const property_map& expected_general_meta = expected.get_general_metadata();
    const property_map& actual_general_meta = actual.get_general_metadata();
    basic_stringstream<frantic::tchar> failMsg;

    // Test general metadata
    AssertionResult generalResult =
        cmpPropertyMap( expected_expr, actual_expr, expected_general_meta, actual_general_meta );
    if( !generalResult )
        failMsg << "\nGeneral Metadata:\n" << generalResult.message();

    // Test channel metadata
    std::vector<frantic::tstring> expected_channels_with_meta;
    std::vector<frantic::tstring> actual_channels_with_meta;
    expected.get_channels_with_metadata( expected_channels_with_meta );
    actual.get_channels_with_metadata( actual_channels_with_meta );
    if( expected_channels_with_meta != actual_channels_with_meta ) {
        basic_stringstream<frantic::tchar> expected_channels_ss;
        basic_stringstream<frantic::tchar> actual_channels_ss;
        for( size_t i = 0; i < expected_channels_with_meta.size(); ++i ) {
            expected_channels_ss << expected_channels_with_meta[i];
        }
        for( size_t i = 0; i < actual_channels_with_meta.size(); ++i ) {
            actual_channels_ss << actual_channels_with_meta[i];
        }
        failMsg << "Streams have different channels with metadata.\nExpected channels: " << expected_channels_ss.str()
                << "\nActual channels: " << actual_channels_ss.str() << std::endl;
        return AssertionFailure() << failMsg.str();
    }
    bool pmSuccess = true;
    basic_stringstream<frantic::tchar> channelFailMsg;
    for( size_t i = 0; i < expected_channels_with_meta.size(); ++i ) {
        frantic::tstring channelName = expected_channels_with_meta[i];
        const property_map* expected_channel_meta = expected.get_channel_metadata( channelName );
        const property_map* actual_channel_meta = actual.get_channel_metadata( channelName );
        AssertionResult pmResult =
            cmpPropertyMap( expected_expr, actual_expr, *expected_channel_meta, *actual_channel_meta );
        if( !pmResult ) {
            pmSuccess = false;
            channelFailMsg << "\n"
                           << channelName << " channel metadata: "
                           << "\n"
                           << pmResult.message();
        }
    }
    if( !pmSuccess )
        failMsg << channelFailMsg.str();
    if( generalResult && pmSuccess )
        return AssertionSuccess();
    return AssertionFailure() << failMsg.str();
}
AssertionResult cmpParticleFile( const char* expected_expr, const char* actual_expr,
                                 const char* /*position_type_hint_expr*/, const frantic::tstring& expected_file,
                                 const frantic::tstring& actual_file, data_type_t position_type_hint ) {
    particle_file_stream_factory_object factory;
    factory.set_position_type_hint( position_type_hint );

    particle_file_metadata expected_meta;
    particle_file_metadata actual_meta;
    // Initialize streams and get their channel_maps
    boost::shared_ptr<particle_istream> expected_stream = factory.create_istream( expected_file, &expected_meta );
    boost::shared_ptr<particle_istream> actual_stream = factory.create_istream( actual_file, &actual_meta );

    AssertionResult particleResult = cmpParticleData( expected_expr, actual_expr, expected_stream, actual_stream );
    AssertionResult metadataResult = cmpMetadata( expected_expr, actual_expr, expected_meta, actual_meta );

    std::stringstream particleFailMsg;
    if( !particleResult )
        particleFailMsg << "\nParticle Data:\n";
    if( particleResult && metadataResult )
        return AssertionSuccess();
    return AssertionFailure() << StringStreamToString( &particleFailMsg ) << particleResult.message()
                              << metadataResult.message();
}

AssertionResult cmpParticleFile( const char* expected_expr, const char* actual_expr,
                                 const frantic::tstring& expected_file, const frantic::tstring& actual_file ) {
    return cmpParticleFile( expected_expr, actual_expr, "no position type hint", expected_file, actual_file,
                            data_type_invalid );
}
} // namespace internal
} // namespace testing
