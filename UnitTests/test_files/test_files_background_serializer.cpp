// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/files/background_serializer.hpp>
#include <frantic/files/files.hpp>
#include <test_files/files_test_utils.hpp>

using namespace boost::assign;
using namespace frantic::files;
using namespace frantic::files::tests;

// Test that the implementation is actually called. Serialize 1 item and wait.
TEST( FilesTest, BackgroundSerializer1 ) {

    background_serializer<frantic::tstring, boost::shared_ptr<serializable>, test_serializer> theBgSerializer(
        test_serializer(), 10 );

    ASSERT_EQ( theBgSerializer.get_capacity(), 10 );

    std::deque<frantic::tstring> keyList;
    std::deque<boost::shared_ptr<serializable>> valueList;

    keyList += _T("Test1");
    valueList += boost::shared_ptr<serializable>( new serializable( 1 ) );

    theBgSerializer.enqueue( keyList.front(), valueList.front(), valueList.front()->size() );
    theBgSerializer.wait_for_idle();

    ASSERT_EQ( theBgSerializer.get_usage(), 0 );

    std::deque<frantic::tstring> serializedKeys;
    std::deque<boost::shared_ptr<serializable>> serializedValues;

    theBgSerializer.get_serializer().copy_serialized_keys( std::back_inserter( serializedKeys ) );
    theBgSerializer.get_serializer().copy_serialized_values( std::back_inserter( serializedValues ) );

    ASSERT_EQ( keyList, serializedKeys );
    ASSERT_EQ( valueList, serializedValues );
}

// Test that the implementation is actually called and that exceeding the queue capacity is handled. Serialize 20 items
// and wait, then check the results are the same.
TEST( FilesTest, BackgroundSerializer20 ) {
    background_serializer<frantic::tstring, boost::shared_ptr<serializable>, test_serializer> theBgSerializer(
        test_serializer(), 10 );

    std::deque<frantic::tstring> keyList;
    std::deque<boost::shared_ptr<serializable>> valueList;

    for( std::size_t i = 0, iEnd = 20; i < iEnd; ++i ) {
        keyList.push_back( _T("Test") + boost::lexical_cast<frantic::tstring>( i + 1 ) );
        valueList.push_back( boost::shared_ptr<serializable>( new serializable( 1 ) ) );
    }

    std::deque<boost::shared_ptr<serializable>>::const_iterator itValue = valueList.begin();
    for( std::deque<frantic::tstring>::const_iterator it = keyList.begin(), itEnd = keyList.end(); it != itEnd;
         ++it, ++itValue )
        theBgSerializer.enqueue( *it, *itValue, ( *itValue )->size() );

    theBgSerializer.wait_for_idle();

    ASSERT_EQ( theBgSerializer.get_usage(), 0 );

    std::deque<frantic::tstring> serializedKeys;
    std::deque<boost::shared_ptr<serializable>> serializedValues;

    theBgSerializer.get_serializer().copy_serialized_keys( std::back_inserter( serializedKeys ) );
    theBgSerializer.get_serializer().copy_serialized_values( std::back_inserter( serializedValues ) );

    ASSERT_EQ( keyList, serializedKeys );
    ASSERT_EQ( valueList, serializedValues );
}

// Test serializing multiple items with multiple worker threads processing the data.
TEST( FilesTest, BackgroundSerializerMultipleWorkers ) {
    background_serializer<frantic::tstring, boost::shared_ptr<serializable>, test_serializer> theBgSerializer(
        test_serializer(), 10 );

    theBgSerializer.set_num_threads( 4 );

    std::deque<frantic::tstring> keyList;
    std::deque<boost::shared_ptr<serializable>> valueList;

    for( std::size_t i = 0, iEnd = 50; i < iEnd; ++i ) {
        keyList.push_back( _T("Test") + boost::lexical_cast<frantic::tstring>( i + 1 ) );
        valueList.push_back( boost::shared_ptr<serializable>( new serializable( 1 ) ) );
    }

    std::deque<boost::shared_ptr<serializable>>::const_iterator itValue = valueList.begin();
    for( std::deque<frantic::tstring>::const_iterator it = keyList.begin(), itEnd = keyList.end(); it != itEnd;
         ++it, ++itValue )
        theBgSerializer.enqueue( *it, *itValue, ( *itValue )->size() );

    theBgSerializer.wait_for_idle();

    ASSERT_EQ( theBgSerializer.get_usage(), 0 );

    std::deque<frantic::tstring> serializedKeys;
    std::deque<boost::shared_ptr<serializable>> serializedValues;

    theBgSerializer.get_serializer().copy_serialized_keys( std::back_inserter( serializedKeys ) );
    theBgSerializer.get_serializer().copy_serialized_values( std::back_inserter( serializedValues ) );

    // We need to sort the lists before comparing them, because the serialization order is not consistent when using
    // multiple worker threads.
    std::stable_sort( serializedKeys.begin(), serializedKeys.end() );
    std::stable_sort( keyList.begin(), keyList.end() );

    std::stable_sort( serializedValues.begin(), serializedValues.end() );
    std::stable_sort( valueList.begin(), valueList.end() );

    ASSERT_EQ( keyList, serializedKeys );
    ASSERT_EQ( valueList, serializedValues );
}

// Test that inhomogeneous sized items are handled correctly. Serialize 20 items of various sizes and wait, then check
// the results are the same.
TEST( FilesTest, BackgroundSerializerSizes ) {
    background_serializer<frantic::tstring, boost::shared_ptr<serializable>, test_serializer> theBgSerializer(
        test_serializer(), 10 );

    std::deque<frantic::tstring> keyList;
    std::deque<boost::shared_ptr<serializable>> valueList;
    for( std::size_t i = 0, iEnd = 20; i < iEnd; ++i ) {
        keyList.push_back( _T("Test") + boost::lexical_cast<frantic::tstring>( i + 1 ) );
        valueList.push_back( boost::shared_ptr<serializable>( new serializable(
            ( i % 4 ) + 1 ) ) ); // Sequence is 1,2,3,4,1,2,3,4,etc. which will exceed capacity on 5th item
    }

    std::deque<frantic::tstring>::const_iterator it = keyList.begin();
    std::deque<boost::shared_ptr<serializable>>::const_iterator itValue = valueList.begin();

    {
        // Hold the mutex in order to prevent any serialization from ocurring until we specify.
        boost::lock_guard<test_serializer> theLock( theBgSerializer.get_serializer() );

        for( std::deque<frantic::tstring>::const_iterator itEnd = keyList.begin() + 4; it != itEnd; ++it, ++itValue )
            ASSERT_TRUE( theBgSerializer.try_enqueue( *it, *itValue, ( *itValue )->size() ) );

        // The first 4 items sum to usage == 10 (which happens to be the capacity too!)
        ASSERT_EQ( theBgSerializer.get_usage(), 10 );

        // The fifth enqueue will fail because we are at capacity.
        ASSERT_TRUE( !theBgSerializer.try_enqueue( *it, *itValue, ( *itValue )->size() ) );
    }

    // We've released the lock on the serializer, so we can enqueue the rest of the items. Some calls will block until
    // capacity is available.
    for( std::deque<frantic::tstring>::const_iterator itEnd = keyList.end(); it != itEnd; ++it, ++itValue )
        theBgSerializer.enqueue( *it, *itValue, ( *itValue )->size() );

    theBgSerializer.wait_for_idle();

    ASSERT_EQ( theBgSerializer.get_usage(), 0 );

    std::deque<frantic::tstring> serializedKeys;
    std::deque<boost::shared_ptr<serializable>> serializedValues;

    theBgSerializer.get_serializer().copy_serialized_keys( std::back_inserter( serializedKeys ) );
    theBgSerializer.get_serializer().copy_serialized_values( std::back_inserter( serializedValues ) );

    ASSERT_EQ( keyList, serializedKeys );
    ASSERT_EQ( valueList, serializedValues );
}

// Test that the implementation is actually asynchronous. Hold a lock preventing serialization until try_enqueue returns
// false, then wait for completion and check the result.
// This will deadlock (or throw an exception possibly due to ) if the serialization is not called in a separate thread.
TEST( FilesTest, BackgroundSerializerAsync ) {
    background_serializer<frantic::tstring, boost::shared_ptr<serializable>, test_serializer> theBgSerializer(
        test_serializer(), 10 );

    std::deque<frantic::tstring> keyList;
    std::deque<boost::shared_ptr<serializable>> valueList;

    keyList += _T("Test1"), _T("Test2"), _T("Test3"), _T("Test4"), _T("Banana"), _T("Test6"), _T("Test7"),
        _T("Dukes Of Hazard"), _T("Test9"), _T("crev\xe9"), _T("ThisOneShouldBlock");

    for( std::deque<frantic::tstring>::const_iterator it = keyList.begin(), itEnd = keyList.end(); it != itEnd; ++it )
        valueList.push_back( boost::shared_ptr<serializable>( new serializable ) );

    {
        // Hold the mutex in order to prevent any serialization from ocurring until we specify.
        boost::lock_guard<test_serializer> theLock( theBgSerializer.get_serializer() );

        // Insert all but the last item (which should block the thread)
        std::deque<boost::shared_ptr<serializable>>::const_iterator itValue = valueList.begin();
        for( std::deque<frantic::tstring>::const_iterator it = keyList.begin(), itEnd = keyList.end() - 1; it != itEnd;
             ++it, ++itValue )
            ASSERT_TRUE( theBgSerializer.try_enqueue( *it, *itValue, ( *itValue )->size() ) );

        ASSERT_EQ( theBgSerializer.get_usage(), 10 );

        // Assert that adding the last item causes it to block
        ASSERT_TRUE( !theBgSerializer.try_enqueue( keyList.back(), *itValue, ( *itValue )->size() ) );
    }

    // Add the last item now that we've released the lock on the serializer. This should eventually return (ie. not
    // deadlock)
    theBgSerializer.enqueue( keyList.back(), valueList.back(), valueList.back()->size() );

    theBgSerializer.wait_for_idle();

    // Assert that the serializer executed in a different thread than the current one.
    ASSERT_NE( theBgSerializer.get_serializer().get_last_write_id(), boost::this_thread::get_id() );

    std::deque<frantic::tstring> serializedKeys;
    std::deque<boost::shared_ptr<serializable>> serializedValues;

    theBgSerializer.get_serializer().copy_serialized_keys( std::back_inserter( serializedKeys ) );
    theBgSerializer.get_serializer().copy_serialized_values( std::back_inserter( serializedValues ) );

    ASSERT_EQ( keyList, serializedKeys );
    ASSERT_EQ( valueList, serializedValues );
}

// Test that an exception that propogates out of serializer will result in a
// frantic::files::background_serializer_critical_error containing the original exception and
// all further attempts to enqueue to the background_serializer will fail.
TEST( FilesTest, BackgroundSerializerException ) {
    background_serializer<frantic::tstring, boost::shared_ptr<serializable>, test_serializer> theBgSerializer(
        test_serializer(), 10 );

    std::deque<frantic::tstring> keyList;
    std::deque<boost::shared_ptr<serializable>> valueList;

    keyList += _T("Test1"), _T("Test2"), _T("Test3"), _T("Test4"), _T("Banana");

    for( std::deque<frantic::tstring>::const_iterator it = keyList.begin(), itEnd = keyList.end(); it != itEnd; ++it )
        valueList.push_back( boost::shared_ptr<serializable>( new serializable ) );

    std::deque<boost::shared_ptr<serializable>>::const_iterator itValue = valueList.begin();
    for( std::deque<frantic::tstring>::const_iterator it = keyList.begin(), itEnd = keyList.end(); it != itEnd;
         ++it, ++itValue )
        ASSERT_TRUE( theBgSerializer.try_enqueue( *it, *itValue, ( *itValue )->size() ) );

    // Enqueue a special key-value pair that causes an exception when processed.
    theBgSerializer.enqueue( EXCEPTION_REQUEST_KEY, boost::shared_ptr<serializable>( new serializable ), 1 );

    EXPECT_THROW( theBgSerializer.wait_for_idle(), frantic::files::background_serializer_critical_error );
    EXPECT_ANY_THROW(
        theBgSerializer.enqueue( _T("ThisShouldFail"), boost::shared_ptr<serializable>( new serializable ), 1 ) );

    std::deque<frantic::tstring> serializedKeys;
    std::deque<boost::shared_ptr<serializable>> serializedValues;

    theBgSerializer.get_serializer().copy_serialized_keys( std::back_inserter( serializedKeys ) );
    theBgSerializer.get_serializer().copy_serialized_values( std::back_inserter( serializedValues ) );

    ASSERT_EQ( keyList, serializedKeys );
    ASSERT_EQ( valueList, serializedValues );
}

// Test that nothing bad occurs when submitting a really large item that overfills the serializer's capacity. I expect
// it to allow at least one item even if it exceeds the capacity.
TEST( FilesTest, BackgroundSerializerLargeItem ) {
    background_serializer<frantic::tstring, boost::shared_ptr<serializable>, test_serializer> theBgSerializer(
        test_serializer(), 1 );

    std::deque<frantic::tstring> keyList;
    std::deque<boost::shared_ptr<serializable>> valueList;

    keyList.push_back( _T("ThisIsWayBiggerThanTheQueueCapacity") );
    valueList.push_back(
        boost::shared_ptr<serializable>( new serializable( std::numeric_limits<unsigned>::max() - 100u ) ) );

    ASSERT_TRUE( theBgSerializer.try_enqueue( keyList.front(), valueList.front(), valueList.front()->size() ) );

    theBgSerializer.wait_for_idle();

    std::deque<frantic::tstring> serializedKeys;
    std::deque<boost::shared_ptr<serializable>> serializedValues;

    theBgSerializer.get_serializer().copy_serialized_keys( std::back_inserter( serializedKeys ) );
    theBgSerializer.get_serializer().copy_serialized_values( std::back_inserter( serializedValues ) );

    ASSERT_EQ( keyList, serializedKeys );
    ASSERT_EQ( valueList, serializedValues );
}
