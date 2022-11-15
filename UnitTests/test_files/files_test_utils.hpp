// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#ifdef WIN32
#pragma warning( push, 3 )
#endif
#include <boost/assign/std/deque.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#ifdef WIN32
#pragma warning( pop )
#endif

namespace frantic {
namespace files {
namespace tests {
// Represents a serializable item for purposes of testing.
struct serializable {
    frantic::tstring m_keyPath;
    std::size_t m_size;

    serializable()
        : m_size( 1u ) {}

    serializable( std::size_t dataSize )
        : m_size( dataSize ) {}

    std::size_t size() const { return m_size; }
};

inline bool operator==( const serializable& lhs, const serializable& rhs ) {
    return lhs.m_keyPath == rhs.m_keyPath && lhs.m_size == rhs.m_size;
}

template <class T>
struct size_getter_impl {
    static std::size_t apply( const T& val ) { return val.size(); }
};

template <class T>
struct size_getter_impl<boost::shared_ptr<T>> {
    static std::size_t apply( const boost::shared_ptr<T>& pVal ) { return pVal->size(); }
};

// A functor that can get the size of a serializable object (or a shared_ptr to one).
struct size_getter {
    template <class T>
    std::size_t operator()( const T& val ) const {
        // Forward to size_getter_impl<>::apply so we can do template specialization.
        return size_getter_impl<T>::apply( val );
    }
};

// A special serialization path that forces an exception
static const frantic::tchar* EXCEPTION_REQUEST_KEY = _T("Please Throw An Exception Now");

// A serializer object that uses pimpl so that we can treat it with value semantics
class test_serializer {
    struct test_serializer_impl {
        boost::mutex m_mutex;
        boost::thread::id m_lastWriteID;
        std::deque<frantic::tstring> m_serializedKeys;
        std::deque<boost::shared_ptr<serializable>> m_serializedValues;

        test_serializer_impl() {}

        void serialize( const frantic::tstring& path, boost::shared_ptr<serializable> pSerializable ) {
            boost::unique_lock<boost::mutex> theLock( m_mutex );
            boost::this_thread::sleep( boost::posix_time::milliseconds( 100 ) ); // Simulate work.

            if( path == EXCEPTION_REQUEST_KEY )
                boost::throw_exception( std::runtime_error( "Serializing failure" ) );

            m_serializedKeys.push_back( path );
            m_serializedValues.push_back( pSerializable );

            m_lastWriteID = boost::this_thread::get_id();
        }

        boost::shared_ptr<serializable> deserialize( const frantic::tstring& path ) {
            boost::shared_ptr<serializable> pResult( new serializable );

            pResult->m_keyPath = path;

            return pResult;
        }

        boost::thread::id get_last_write_id() const { return m_lastWriteID; }

        template <class OutputIterator>
        void copy_serialized_keys( OutputIterator it ) const {
            std::copy( m_serializedKeys.begin(), m_serializedKeys.end(), it );
        }

        template <class OutputIterator>
        void copy_serialized_values( OutputIterator it ) const {
            std::copy( m_serializedValues.begin(), m_serializedValues.end(), it );
        }

        void lock() { m_mutex.lock(); }

        void unlock() { m_mutex.unlock(); }
    };

    boost::shared_ptr<test_serializer_impl> m_pImpl;

  public:
    test_serializer()
        : m_pImpl( new test_serializer_impl ) {}

    void serialize( const frantic::tstring& path, boost::shared_ptr<serializable> pSerializable ) {
        m_pImpl->serialize( path, pSerializable );
    }

    boost::shared_ptr<serializable> deserialize( const frantic::tstring& path ) { return m_pImpl->deserialize( path ); }

    boost::thread::id get_last_write_id() const { return m_pImpl->get_last_write_id(); }

    template <class OutputIterator>
    void copy_serialized_keys( OutputIterator it ) const {
        m_pImpl->copy_serialized_keys( it );
    }

    template <class OutputIterator>
    void copy_serialized_values( OutputIterator it ) const {
        m_pImpl->copy_serialized_values( it );
    }

    void lock() { m_pImpl->lock(); }

    void unlock() { m_pImpl->unlock(); }
};
} // namespace tests
} // namespace files
} // namespace frantic
