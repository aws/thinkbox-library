// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"
#include <frantic/files/background_serializer.hpp>
#include <frantic/files/sequence_cache.hpp>
#include <frantic/math/utils.hpp>
#include <fstream>
#include <test_files/files_test_utils.hpp>

using namespace boost::assign;
using namespace frantic::files;
using namespace frantic::files::tests;

/**
 * Testing the behavior of an empty cache here.
 */
TEST( FilesTest, SequenceCacheEmpty ) {

    typedef boost::shared_ptr<serializable> value_type;

    sequence_cache<value_type, test_serializer, size_getter> theCache( 10, test_serializer(), size_getter() );

    ASSERT_TRUE( theCache.empty() );

    ASSERT_EQ( theCache.get_usage(), 0 );

    ASSERT_EQ( theCache.find( 1.0 ), value_type() );

    ASSERT_TRUE( !frantic::math::is_finite( theCache.find_nearest_key( 1.0 ) ) );

    ASSERT_EQ( theCache.find_bracketing_keys( 1.0 ),
               std::make_pair( -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() ) );

    // This should return immediately without doing anything.
    theCache.wait_for_pending();
}

/**
 * Testing the behavior of a single entry in the cache
 */
TEST( FilesTest, SequenceCacheSingle ) {

    typedef boost::shared_ptr<serializable> value_type;

    test_serializer theSerializer;

    sequence_cache<value_type, test_serializer, size_getter> theCache( 10, theSerializer, size_getter() );

    value_type theVal( new serializable );

    theCache.insert( 1.0, theVal );

    ASSERT_TRUE( !theCache.empty() );
    ASSERT_TRUE( theCache.get_usage() > 0 );

    ASSERT_EQ( theCache.find( 1.0 ), theVal );
    ASSERT_EQ( theCache.find( 1.5 ), value_type() );
    ASSERT_EQ( theCache.find( -1000.0 ), value_type() );

    ASSERT_EQ( theCache.find_nearest_key( 1.0 ), 1.0 );
    ASSERT_EQ( theCache.find_nearest_key( 100.0 ), 1.0 );
    ASSERT_EQ( theCache.find_nearest_key( -69.5 ), 1.0 );

    ASSERT_EQ( theCache.find_bracketing_keys( 1.5 ), std::make_pair( 1.0, std::numeric_limits<double>::infinity() ) );
    ASSERT_EQ( theCache.find_bracketing_keys( 0.5 ), std::make_pair( -std::numeric_limits<double>::infinity(), 1.0 ) );
}

/**
 * Testing the behavior of a single entry in the cache being flushed from memory both with and without the disk cache
 * enabled.
 */
TEST( FilesTest, SequenceCacheSingleWithFlush ) {

    typedef boost::shared_ptr<serializable> value_type;

    test_serializer theSerializer;

    sequence_cache<value_type, test_serializer, size_getter> theCache( 10, theSerializer, size_getter() );

    theCache.insert( 1.0, value_type( new serializable ) );
    theCache.flush( false );
    theCache.wait_for_pending();

    std::deque<frantic::tstring> serializedKeys;
    theSerializer.copy_serialized_keys( std::back_inserter( serializedKeys ) );

    // We are asserting that the serializer was never called because we never set the disk path.
    ASSERT_TRUE( serializedKeys.empty() );

    // The item should still be in the cache though, since flushing only ensures the cache entries are available on disk
    // and we didn't specify a path.
    ASSERT_TRUE( !theCache.empty() );
    ASSERT_EQ( theCache.get_usage(), 1 );
    ASSERT_EQ( theCache.find_nearest_key( 1.0 ), 1.0 );
    ASSERT_TRUE( theCache.find( 1.0 ) != value_type() );

    // Clearing the cache will drop everything and leave it empty.
    theCache.clear();

    ASSERT_TRUE( theCache.empty() );
    ASSERT_EQ( theCache.get_usage(), 0 );
    ASSERT_TRUE( !frantic::math::is_finite( theCache.find_nearest_key( 1.0 ) ) );
    ASSERT_TRUE( theCache.find( 1.0 ) == value_type() );

    theCache.set_disk_path( _T("TestPattern_####.ext"),
                            sequence_cache<value_type, test_serializer, size_getter>::keep_existing );

    frantic::files::filename_pattern thePattern( theCache.get_disk_path_pattern() );

    value_type theVal = theSerializer.deserialize( thePattern[1.0] );

    theCache.insert( 1.0, theVal );
    theCache.flush( true ); // We did set a disk path this time, so this will actually serialize. It removes all items
                            // from memory too.
    theCache.wait_for_pending();

    // Get a copy of the keys processed
    theSerializer.copy_serialized_keys( std::back_inserter( serializedKeys ) );

    ASSERT_TRUE( !theCache.empty() );
    ASSERT_EQ( theCache.get_usage(), 0 );
    ASSERT_EQ( theCache.find_nearest_key( 1.0 ), 1.0 );
    ASSERT_NE( theCache.find( 1.0 ), value_type() ); // We should not be getting a default constructed value back.
    ASSERT_NE( theCache.find( 1.0 ),
               theVal ); // We should be getting a different object back (pointing to an equivalent object though).
    ASSERT_EQ( *theCache.find( 1.0 ), *theVal ); // We should be getting an equivalent object back.
    ASSERT_EQ( theCache.find( 1.0 )->m_keyPath, theVal->m_keyPath );
    ASSERT_EQ( serializedKeys.size(), 1 );
    ASSERT_EQ( serializedKeys.front(), thePattern[1.0] );
}

/**
 * Test the cache behavior when inserting more items than there is space for.
 */
TEST( FilesTest, SequenceCacheSerialize ) {
    typedef boost::shared_ptr<serializable> value_type;

    test_serializer theSerializer;

    sequence_cache<value_type, test_serializer, size_getter> theCache( 10, theSerializer, size_getter() );

    theCache.set_disk_path( _T("TestPattern_####.ext"),
                            sequence_cache<value_type, test_serializer, size_getter>::keep_existing );

    frantic::files::filename_pattern thePattern( theCache.get_disk_path_pattern() );

    std::deque<frantic::tstring> expectedPaths;

    for( std::size_t i = 0, iEnd = theCache.get_capacity(); i < iEnd; ++i ) {
        theCache.insert( static_cast<double>( i ), value_type( new serializable ) );

        expectedPaths.push_back( thePattern[static_cast<double>( i )] );
    }

    // There should be nothing serialized at this point, because we merely reached the capacity.
    theCache.wait_for_pending();

    ASSERT_EQ( theCache.get_usage(), theCache.get_capacity() );

    // This next insertion will drop the least recently used item (ie. frame 0)
    theCache.insert( static_cast<double>( theCache.get_capacity() ), value_type( new serializable ) );
    expectedPaths.push_back( thePattern[static_cast<double>( theCache.get_capacity() )] );

    theCache.wait_for_pending(); // We should have serialized 1 item.

    std::deque<frantic::tstring> serializedKeys;

    theSerializer.copy_serialized_keys( std::back_inserter( serializedKeys ) );

    // Assert that we have found the least recently used item in serialized items list.
    // We may have serialized other items though b/c sequence_cache will try to stay busy in background.
    // NOTE: On Linux, gtest's ASSERT_NE tried to print the iterator as if it were a container, so gave a build failure.
    ASSERT_FALSE( std::find( serializedKeys.begin(), serializedKeys.end(), expectedPaths.front() ) ==
                  serializedKeys.end() );

    // Drop everything and check that it was all serialized.
    theCache.flush( true );
    theCache.wait_for_pending();

    serializedKeys.clear();

    theSerializer.copy_serialized_keys( std::back_inserter( serializedKeys ) );

    // Copy the serialization queue and sort it (and the expectedPaths container too) so that we ensure the same files
    // were processed without caring about the order.
    std::stable_sort( serializedKeys.begin(), serializedKeys.end() );
    std::stable_sort( expectedPaths.begin(), expectedPaths.end() );

    ASSERT_EQ( serializedKeys, expectedPaths );
}

// Google fixture for all tests
class FilesTestWithFixture : public ::testing::Test {
  public:
    boost::filesystem::path tempRootFolder, tempFolder;

    virtual void SetUp() {
        tempRootFolder = boost::filesystem::temp_directory_path();
        do {
            tempFolder = boost::filesystem::unique_path( tempRootFolder / "%%%%-%%%%-%%%%-%%%%" );
        } while( !boost::filesystem::create_directory( tempFolder ) ); // Loops while the directory already exists.
    }

    virtual void TearDown() { boost::filesystem::remove_all( tempFolder ); }
};

// Tests the behaviour of sequence_cache::set_disk_path( path, sequence_cache::synchronize )
TEST_F( FilesTestWithFixture, SequenceCacheSynchronize ) {

    std::cout << "testSequenceCacheSynchronize() using folder: " + tempFolder.string() << std::endl;

    frantic::files::filename_sequence theSequence( tempFolder.string<frantic::tstring>(), _T("TestPattern_####.ext") );

    for( std::size_t i = 0; i <= 100ULL; ++i ) {
        std::ofstream theFile;
        theFile.exceptions( std::ios::failbit | std::ios::badbit );

        theFile.open( theSequence[static_cast<double>( i )].c_str(), std::ios::out | std::ios::trunc );
        theFile
            << _T("This is a test file generated by the unit test \"testSequenceCacheSynchronize\". The assoicated ")
               _T("id is: ")
            << i << std::endl;
        theFile.close();
    }

    typedef boost::shared_ptr<serializable> value_type;

    sequence_cache<value_type, test_serializer, size_getter> theCache( 10, test_serializer(), size_getter() );

    theCache.set_disk_path( theSequence.get_filename_pattern().get_pattern(),
                            sequence_cache<value_type, test_serializer, size_getter>::synchronize );

    ASSERT_TRUE( !theCache.empty() );
    ASSERT_EQ( theCache.find_nearest_key( -100.5 ), 0.0 );
    ASSERT_EQ( theCache.find_nearest_key( 231.0 ), 100.0 );

    for( frantic::files::frame_set::const_iterator it = theSequence.get_frame_set().begin(),
                                                   itEnd = theSequence.get_frame_set().end();
         it != itEnd; ++it ) {
        ASSERT_EQ( theCache.find_nearest_key( *it ), *it );
        ASSERT_TRUE( theCache.find( *it ) !=
                     value_type() ); // We should not be getting a default constructed value back.
        ASSERT_TRUE( theCache.find( *it )->m_keyPath ==
                     theSequence[*it] ); // We should be getting an equivalent object back.
    }
}
