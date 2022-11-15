// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// TODO: Much of this stuff is shared with the rle level set named channels.  We'll probably want to refactor the shared
// code into a common
//       base class once this has settled down a bit.

#pragma once

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/graphics/vector3.hpp>

namespace frantic {
namespace geometry {

using frantic::channels::data_type_t;
using frantic::graphics::vector3;

class trimesh3;
class trimesh3_vertex_channel;
class trimesh3_face_channel;

///////////////////////
// VERTEX CHANNEL GENERAL ACCESSOR BASE
///////////////////////

// This class is used to provide faster access within a named channel.
template <class FaceVector, class DataVector>
class trimesh3_vertex_channel_general_accessor_base {
  protected:
    FaceVector* m_faces;
    DataVector* m_data;

    std::size_t m_arity;         // # of data values for this channel per particle
    data_type_t m_dataType;      // primitive data type
    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    // Function which knows how to do weighted sum combinations of this data type
    frantic::channels::channel_weighted_sum_combine_function_t m_weightedSumFunction;

    bool m_hasCustomFaces;

    trimesh3_vertex_channel_general_accessor_base( FaceVector* faces, DataVector* data, std::size_t arity,
                                                   data_type_t dataType, bool hasCustomFaces )
        : m_faces( faces )
        , m_data( data )
        , m_arity( arity )
        , m_dataType( dataType )
        , m_hasCustomFaces( hasCustomFaces ) {
        if( m_faces == 0 && m_data != 0 )
            throw std::runtime_error(
                "trimesh3_vertex_channel_accessor: Tried to construct an accessor with invalid faces "
                "but valid vertex data.  The faces must always be valid if the vertex data is." );

        m_weightedSumFunction = channel_weighted_sum_combine_function( dataType );
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    friend class trimesh3;

  public:
    trimesh3_vertex_channel_general_accessor_base()
        : m_faces( 0 )
        , m_data( 0 )
        , m_arity( 0 )
        , m_dataType( channels::data_type_int8 )
        , m_primitiveSize( 0 )
        , m_weightedSumFunction( 0 )
        , m_hasCustomFaces( false ) {}

    std::size_t size() const { return m_data ? ( m_data->size() / m_primitiveSize ) : 0; }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    bool valid() const { return m_data != 0; }

    operator bool() const { return m_data != 0; }

    frantic::channels::channel_weighted_sum_combine_function_t get_weighted_sum_combine_function() const {
        return m_weightedSumFunction;
    }

    /**
     * This returns the size of the primitive of this channel.  That's equivalent to arity() * (data type size).
     */
    std::size_t primitive_size() const { return m_primitiveSize; }

    void dump( std::ostream& out ) const {
        if( m_faces ) {
            out << "m_faces size: " << m_faces->size() << "\n";
        } else {
            out << "m_faces is null\n";
        }
        if( m_data ) {
            out << "m_data size: " << m_data->size() << "\n";
        } else {
            out << "m_data is null\n";
        }
        out << "m_arity: " << m_arity << "\n";
        out << "m_dataType: " << m_dataType << "\n";
        out << "m_dataType str: "
            << frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) << "\n";
        out << "m_primitiveSize: " << m_primitiveSize << "\n";
        out << "m_hasCustomFaces: " << m_hasCustomFaces << "\n";
        out << "size(): " << size() << "\n";
    }

    ///////////
    // Access to the faces
    ///////////

    bool has_custom_faces() const { return m_hasCustomFaces; }

    std::size_t face_count() const { return m_faces ? m_faces->size() : 0; }

    const vector3& face( std::size_t f ) const { return ( *m_faces )[f]; }

    ///////////
    // Access to the data
    ///////////

    void print( std::ostream& out, std::size_t i ) const {
        channels::channel_data_type_print( out, m_arity, m_dataType, data( i ) );
    }

    const char* data( std::size_t i ) const { return m_data->ptr_at( i * m_primitiveSize ); }

    void get_barycentric( int faceIndex, const frantic::graphics::vector3f& barycentricCoords, char* outResult ) const {
        vector3 face = ( *m_faces )[faceIndex];
        const char* dataArray[3];
        dataArray[0] = m_data->ptr_at( face.x * m_primitiveSize );
        dataArray[1] = m_data->ptr_at( face.y * m_primitiveSize );
        dataArray[2] = m_data->ptr_at( face.z * m_primitiveSize );

        m_weightedSumFunction( &barycentricCoords[0], dataArray, 3, m_arity, outResult );
    }
};

///////////////////////
// VERTEX CHANNEL GENERAL ACCESSOR
///////////////////////

// Forward declaration of the const version of this class, so it has access to the members
class const_trimesh3_vertex_channel_general_accessor;

// The non-const accessor adds additional methods for writing to the channel
class trimesh3_vertex_channel_general_accessor
    : public trimesh3_vertex_channel_general_accessor_base<std::vector<vector3>, frantic::graphics::raw_byte_buffer> {
  protected:
    trimesh3_vertex_channel_general_accessor( std::vector<vector3>* faces, frantic::graphics::raw_byte_buffer* data,
                                              std::size_t arity, data_type_t dataType, bool hasCustomFaces )
        : trimesh3_vertex_channel_general_accessor_base<std::vector<vector3>, frantic::graphics::raw_byte_buffer>(
              faces, data, arity, dataType, hasCustomFaces ) {}

    friend class const_trimesh3_vertex_channel_general_accessor;
    friend class trimesh3_vertex_channel;

  public:
    trimesh3_vertex_channel_general_accessor() {}

    ///////////
    // Access to the faces
    ///////////

    vector3& face( std::size_t f ) { return ( *m_faces )[f]; }

    // NOTE: When code adds a face here, it is its responsibility to also add the corresponding face to the main mesh
    // faces and
    //       any other named channels with custom faces.
    //       Failure to do so will result in an exception when the accessor for a channel is requested.

    void add_face( const vector3& face ) {
        if( has_custom_faces() )
            m_faces->push_back( face );
        else
            throw std::runtime_error(
                "trimesh3_vertex_channel_accessor.add_face: Tried to add a face to a vertex channel "
                "which doesn't have custom faces." );
    }

    void add_face( int v0, int v1, int v2 ) {
        if( has_custom_faces() )
            m_faces->push_back( vector3( v0, v1, v2 ) );
        else
            throw std::runtime_error(
                "trimesh3_vertex_channel_accessor.add_face: Tried to add a face to a vertex channel "
                "which doesn't have custom faces." );
    }

    // Adds an arbitrary face to the mesh by triangulating it
    void add_face( const std::vector<int>& vertexIndices ) {
        if( has_custom_faces() ) {
            if( vertexIndices.size() > 2 ) {
                int iForward = 0, iBackward = (int)vertexIndices.size() - 1;
                while( iBackward - iForward > 1 ) {
                    // Alternate incrementing iForward and decrementing iBackward to
                    // avoid having all triangles emanate from the same vertex.
                    if( ( iBackward - iForward ) % 2 == 0 ) {
                        m_faces->push_back(
                            vector3( vertexIndices[iForward], vertexIndices[iForward + 1], vertexIndices[iBackward] ) );
                        iForward++;
                    } else {
                        m_faces->push_back( vector3( vertexIndices[iBackward - 1], vertexIndices[iBackward],
                                                     vertexIndices[iForward] ) );
                        iBackward--;
                    }
                }
            }
        } else {
            throw std::runtime_error(
                "trimesh3_vertex_channel_accessor.add_face: Tried to add a face to a vertex channel "
                "which doesn't have custom faces." );
        }
    }

    ///////////
    // Access to the data
    ///////////

    char* data( std::size_t i ) { return m_data->ptr_at( i * m_primitiveSize ); }

    char* add_vertex() { return m_data->add_element( m_primitiveSize ); }

    void add_vertices( size_t numVertsToAdd ) {
        m_data->resize_with_exponential_growth( m_data->size() + numVertsToAdd * m_primitiveSize );
    }

    void set_vertex_count( size_t vertCount ) { m_data->resize_with_exponential_growth( vertCount * m_primitiveSize ); }

    void reserve_vertices( size_t vertCount ) { m_data->reserve( vertCount * m_primitiveSize ); }

    void add_faces( std::size_t numFacesToAdd ) { m_faces->resize( m_faces->size() + numFacesToAdd ); }

    void set_face_count( std::size_t faceCount ) { m_faces->resize( faceCount ); }
};

///////////////////////
// CONST VERTEX CHANNEL GENERAL ACCESSOR
///////////////////////

// The const accessor just uses all the methods in the base accessor class
class const_trimesh3_vertex_channel_general_accessor
    : public trimesh3_vertex_channel_general_accessor_base<const std::vector<vector3>,
                                                           const frantic::graphics::raw_byte_buffer> {
  protected:
    const_trimesh3_vertex_channel_general_accessor( const std::vector<vector3>* faces,
                                                    const frantic::graphics::raw_byte_buffer* data, std::size_t arity,
                                                    data_type_t dataType, bool hasCustomFaces )
        : trimesh3_vertex_channel_general_accessor_base<const std::vector<vector3>,
                                                        const frantic::graphics::raw_byte_buffer>(
              faces, data, arity, dataType, hasCustomFaces ) {}

    friend class trimesh3_vertex_channel;

  public:
    const_trimesh3_vertex_channel_general_accessor() {}

    const_trimesh3_vertex_channel_general_accessor( const trimesh3_vertex_channel_general_accessor& rhs )
        : trimesh3_vertex_channel_general_accessor_base<const std::vector<vector3>,
                                                        const frantic::graphics::raw_byte_buffer>(
              rhs.m_faces, rhs.m_data, rhs.m_arity, rhs.m_dataType, rhs.m_hasCustomFaces ) {}
};

///////////////////////
// VERTEX CHANNEL ACCESSOR
///////////////////////

// Forward declaration of the const version of this class, so it has access to the members
template <class DataType>
class const_trimesh3_vertex_channel_accessor;

template <class DataType>
class trimesh3_vertex_channel_accessor : public trimesh3_vertex_channel_general_accessor {
  protected:
    trimesh3_vertex_channel_accessor( std::vector<vector3>* faces, frantic::graphics::raw_byte_buffer* data,
                                      std::size_t arity, data_type_t dataType, bool hasCustomFaces )
        : trimesh3_vertex_channel_general_accessor( faces, data, arity, dataType, hasCustomFaces ) {}

    friend class const_trimesh3_vertex_channel_accessor<DataType>;
    friend class trimesh3_vertex_channel;

  public:
    trimesh3_vertex_channel_accessor() {}

    ///////////
    // Access to the data
    ///////////

    const DataType& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }

    DataType& operator[]( std::size_t i ) {
        return *reinterpret_cast<DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }

    void add_vertex( const DataType& data ) {
        *reinterpret_cast<DataType*>( m_data->add_element( sizeof( DataType ) ) ) = data;
    }

    DataType get_barycentric( int faceIndex, const frantic::graphics::vector3f& barycentricCoords ) const {
        vector3 face = ( *m_faces )[faceIndex];
        return frantic::channels::channel_data_type_traits<DataType>::barycentric_combine(
            barycentricCoords, *reinterpret_cast<const DataType*>( m_data->ptr_at( face.x * sizeof( DataType ) ) ),
            *reinterpret_cast<const DataType*>( m_data->ptr_at( face.y * sizeof( DataType ) ) ),
            *reinterpret_cast<const DataType*>( m_data->ptr_at( face.z * sizeof( DataType ) ) ) );
    }
};

///////////////////////
// CONST VERTEX CHANNEL ACCESSOR
///////////////////////

// The const accessor just uses all the methods in the base accessor class
template <class DataType>
class const_trimesh3_vertex_channel_accessor : public const_trimesh3_vertex_channel_general_accessor {

    const_trimesh3_vertex_channel_accessor( const std::vector<frantic::graphics::vector3>* faces,
                                            const frantic::graphics::raw_byte_buffer* data, std::size_t arity,
                                            channels::data_type_t dataType, bool hasCustomFaces )
        : const_trimesh3_vertex_channel_general_accessor( faces, data, arity, dataType, hasCustomFaces ) {}

    friend class trimesh3_vertex_channel;

  public:
    const_trimesh3_vertex_channel_accessor() {}

    const DataType& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }

    DataType get_barycentric( int faceIndex, const frantic::graphics::vector3f& barycentricCoords ) const {
        vector3 face = ( *m_faces )[faceIndex];
        return frantic::channels::channel_data_type_traits<DataType>::barycentric_combine(
            barycentricCoords, *reinterpret_cast<const DataType*>( m_data->ptr_at( face.x * sizeof( DataType ) ) ),
            *reinterpret_cast<const DataType*>( m_data->ptr_at( face.y * sizeof( DataType ) ) ),
            *reinterpret_cast<const DataType*>( m_data->ptr_at( face.z * sizeof( DataType ) ) ) );
    }
};

///////////////////////
// VERTEX CHANNEL CVT ACCESSOR
///////////////////////

template <class T>
class trimesh3_vertex_channel_cvt_accessor
    : protected trimesh3_vertex_channel_general_accessor_base<std::vector<vector3>,
                                                              frantic::graphics::raw_byte_buffer> {
  protected:
    channels::channel_type_convertor_function_t m_convertGet;
    channels::channel_type_convertor_function_t m_convertSet;

    trimesh3_vertex_channel_cvt_accessor( std::vector<vector3>* faces, frantic::graphics::raw_byte_buffer* data,
                                          std::size_t sourceDataArity, data_type_t sourceDataType, bool hasCustomFaces,
                                          const frantic::tstring& channelNameForErrorMessage )
        : trimesh3_vertex_channel_general_accessor_base<std::vector<vector3>, frantic::graphics::raw_byte_buffer>(
              faces, data, sourceDataArity, sourceDataType, hasCustomFaces ) {
        data_type_t targetDataType = frantic::channels::channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = frantic::channels::channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity )
            throw std::runtime_error(
                std::string() + "trimesh3_vertex_channel_cvt_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );

        m_convertGet =
            get_channel_type_convertor_function( sourceDataType, targetDataType, channelNameForErrorMessage );
        m_convertSet =
            get_channel_type_convertor_function( targetDataType, sourceDataType, channelNameForErrorMessage );
    }

    char* data( std::size_t i ) { return m_data->ptr_at( i * m_primitiveSize ); }

    const char* data( std::size_t i ) const { return m_data->ptr_at( i * m_primitiveSize ); }

    friend class trimesh3_vertex_channel;

  public:
    typedef T value_type;

    std::size_t size() const {
        return trimesh3_vertex_channel_general_accessor_base<std::vector<vector3>,
                                                             frantic::graphics::raw_byte_buffer>::size();
    }

    bool has_custom_faces( void ) const {
        return trimesh3_vertex_channel_general_accessor_base<std::vector<vector3>,
                                                             frantic::graphics::raw_byte_buffer>::has_custom_faces();
    }

    std::size_t face_count( void ) const {
        return trimesh3_vertex_channel_general_accessor_base<std::vector<vector3>,
                                                             frantic::graphics::raw_byte_buffer>::face_count();
    }

    vector3& face( std::size_t faceNumber ) { return ( *m_faces )[faceNumber]; }

    const T get( std::size_t i ) const {
        T result;

        m_convertGet( (char*)&result, data( i ), frantic::channels::channel_data_type_traits<T>::arity() );

        return result;
    }

    void set( std::size_t i, const T& val ) {
        m_convertSet( data( i ), (char*)&val, frantic::channels::channel_data_type_traits<T>::arity() );
    }

    const T get_barycentric( std::size_t faceNumber, const frantic::graphics::vector3f& barycentricCoords ) {
        const vector3& vertexNumbers( face( faceNumber ) );
        const T a = get( vertexNumbers.x );
        const T b = get( vertexNumbers.y );
        const T c = get( vertexNumbers.z );

        return frantic::channels::channel_data_type_traits<T>::barycentric_combine( barycentricCoords, a, b, c );
    }

    void set_vertex_count( size_t vertCount ) { m_data->resize_with_exponential_growth( vertCount * m_primitiveSize ); }

    void reserve_vertices( size_t vertCount ) { m_data->reserve( vertCount ); }
};

///////////////////////
// CONST VERTEX CHANNEL CVT ACCESSOR
///////////////////////

template <class T>
class const_trimesh3_vertex_channel_cvt_accessor
    : protected trimesh3_vertex_channel_general_accessor_base<const std::vector<vector3>,
                                                              const frantic::graphics::raw_byte_buffer> {
  protected:
    channels::channel_type_convertor_function_t m_convertGet;

    const_trimesh3_vertex_channel_cvt_accessor( const std::vector<vector3>* faces,
                                                const frantic::graphics::raw_byte_buffer* data,
                                                std::size_t sourceDataArity, data_type_t sourceDataType,
                                                bool hasCustomFaces,
                                                const frantic::tstring& channelNameForErrorMessage )
        : trimesh3_vertex_channel_general_accessor_base<const std::vector<vector3>,
                                                        const frantic::graphics::raw_byte_buffer>(
              faces, data, sourceDataArity, sourceDataType, hasCustomFaces ) {
        data_type_t targetDataType = frantic::channels::channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = frantic::channels::channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity )
            throw std::runtime_error(
                std::string() + "const_trimesh3_vertex_channel_cvt_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );

        m_convertGet =
            get_channel_type_convertor_function( sourceDataType, targetDataType, channelNameForErrorMessage );
    }

    friend class trimesh3_vertex_channel;

  public:
    typedef T value_type;

    std::size_t size() const {
        return trimesh3_vertex_channel_general_accessor_base<const std::vector<vector3>,
                                                             const frantic::graphics::raw_byte_buffer>::size();
    }

    bool has_custom_faces( void ) const {
        return trimesh3_vertex_channel_general_accessor_base<
            const std::vector<vector3>, const frantic::graphics::raw_byte_buffer>::has_custom_faces();
    }

    std::size_t face_count( void ) const {
        return trimesh3_vertex_channel_general_accessor_base<const std::vector<vector3>,
                                                             const frantic::graphics::raw_byte_buffer>::face_count();
    }

    const vector3& face( std::size_t faceNumber ) const {
        return trimesh3_vertex_channel_general_accessor_base<
            const std::vector<vector3>, const frantic::graphics::raw_byte_buffer>::face( faceNumber );
    }

    const T get( std::size_t i ) const {
        T result;

        m_convertGet( (char*)&result, data( i ), frantic::channels::channel_data_type_traits<T>::arity() );

        return result;
    }

    const T operator[]( std::size_t i ) const { return get( i ); }

    const T get_barycentric( std::size_t faceNumber, const frantic::graphics::vector3f& barycentricCoords ) {
        const vector3& vertexNumbers( face( faceNumber ) );
        const T a = get( vertexNumbers.x );
        const T b = get( vertexNumbers.y );
        const T c = get( vertexNumbers.z );

        return frantic::channels::channel_data_type_traits<T>::barycentric_combine( barycentricCoords, a, b, c );
    }
};

class trimesh3_vertex_channel {
    frantic::tstring m_name; // string name of the channel
    std::size_t m_arity;     // # of data values for this channel per particle
    data_type_t m_dataType;  // primitive data type

    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    bool m_hasCustomFaces;
    std::vector<vector3> m_faces;
    frantic::graphics::raw_byte_buffer m_data;

    trimesh3_vertex_channel( const frantic::tstring& name, std::size_t arity, data_type_t dataType )
        : m_name( name )
        , m_arity( arity )
        , m_dataType( dataType )
        , m_hasCustomFaces( false ) {
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    void set( const frantic::tstring& name, std::size_t arity, data_type_t dataType, bool hasCustomFaces ) {
        m_name = name;
        m_arity = arity;
        m_dataType = dataType;
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
        m_hasCustomFaces = hasCustomFaces;
    }

    friend class trimesh3;

  public:
    trimesh3_vertex_channel& operator=( const trimesh3_vertex_channel& rhs ) {
        m_name = rhs.m_name;
        m_arity = rhs.m_arity;
        m_dataType = rhs.m_dataType;
        m_primitiveSize = rhs.m_primitiveSize;
        m_faces = rhs.m_faces;
        m_data = rhs.m_data;
        m_hasCustomFaces = rhs.m_hasCustomFaces;

        return *this;
    }

    std::size_t size() const { return m_data.size() / m_primitiveSize; }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    /**
     * This returns the size of the primitive of this channel.  That's equivalent to arity() * (data type size).
     */
    std::size_t primitive_size() const { return m_primitiveSize; }

    template <typename T>
    bool supports_type() const {
        return channels::channel_data_type_traits<T>::arity() == m_arity &&
               channels::channel_data_type_traits<T>::dataType == m_dataType;
    }

    ///////////
    // Access to the faces
    ///////////

    bool has_custom_faces() const { return m_hasCustomFaces; }

    std::size_t face_count() const { return m_faces.size(); }

    const vector3& face( std::size_t f ) const { return m_faces[f]; }

    vector3& face( std::size_t f ) { return m_faces[f]; }

    ///////////
    // Access to the data
    ///////////

    trimesh3_vertex_channel_general_accessor get_general_accessor( std::vector<vector3>* primaryFaces ) {
        if( !m_hasCustomFaces ) {
            return trimesh3_vertex_channel_general_accessor( primaryFaces, &m_data, m_arity, m_dataType, false );
        } else if( m_faces.size() == primaryFaces->size() ) {
            return trimesh3_vertex_channel_general_accessor( &m_faces, &m_data, m_arity, m_dataType, true );
        } else {
            throw std::runtime_error(
                "trimesh3_vertex_channel.get_general_accessor: The face count of vertex channel \"" +
                frantic::strings::to_string( m_name ) +
                "\" was not kept synchronized with that of the primary mesh.  The mesh's face count is " +
                boost::lexical_cast<std::string>( primaryFaces->size() ) + ", but the channel's face count is " +
                boost::lexical_cast<std::string>( m_faces.size() ) + "." );
        }
    }

    const_trimesh3_vertex_channel_general_accessor
    get_general_accessor( const std::vector<vector3>* primaryFaces ) const {
        if( !m_hasCustomFaces ) {
            return const_trimesh3_vertex_channel_general_accessor( primaryFaces, &m_data, m_arity, m_dataType, false );
        } else if( m_faces.size() == primaryFaces->size() ) {
            return const_trimesh3_vertex_channel_general_accessor( &m_faces, &m_data, m_arity, m_dataType, true );
        } else {
            throw std::runtime_error(
                "trimesh3_vertex_channel.get_general_accessor: The face count of vertex channel \"" +
                frantic::strings::to_string( m_name ) +
                "\" was not kept synchronized with that of the primary mesh.  The mesh's face count is " +
                boost::lexical_cast<std::string>( primaryFaces->size() ) + ", but the channel's face count is " +
                boost::lexical_cast<std::string>( m_faces.size() ) + "." );
        }
    }

    // This returns an accessor to this channel, viewed as the requested type.  When the trimesh3 creates an accessor
    // for someone, it needs to pass a pointer to its faces array.  This is done so that the accessor can assume the
    // faces array is valid for performance reasons.
    template <class T>
    trimesh3_vertex_channel_accessor<T> get_accessor( std::vector<vector3>* primaryFaces ) {
        if( channels::channel_data_type_traits<T>::arity() == m_arity &&
            channels::channel_data_type_traits<T>::data_type() == m_dataType ) {
            if( !m_hasCustomFaces ) {
                return trimesh3_vertex_channel_accessor<T>( primaryFaces, &m_data, m_arity, m_dataType, false );
            } else if( m_faces.size() == primaryFaces->size() ) {
                return trimesh3_vertex_channel_accessor<T>( &m_faces, &m_data, m_arity, m_dataType, true );
            } else {
                throw std::runtime_error(
                    "trimesh3_vertex_channel.get_accessor: The face count of vertex channel \"" +
                    frantic::strings::to_string( m_name ) +
                    "\" was not kept synchronized with that of the primary mesh.  The mesh's face count is " +
                    boost::lexical_cast<std::string>( primaryFaces->size() ) + ", but the channel's face count is " +
                    boost::lexical_cast<std::string>( m_faces.size() ) + "." );
            }
        } else {
            throw std::runtime_error( "trimesh3_vertex_channel.get_accessor: Could not access channel \"" +
                                      frantic::strings::to_string( m_name ) +
                                      "\", because the requested type is incompatible." );
        }
    }

    // This returns an accessor to this channel, viewed as the requested type.  When the trimesh3 creates an accessor
    // for someone, it needs to pass a pointer to its faces array.  This is done so that the accessor can assume the
    // faces array is valid for performance reasons.
    template <class T>
    const_trimesh3_vertex_channel_accessor<T> get_accessor( const std::vector<vector3>* primaryFaces ) const {
        if( channels::channel_data_type_traits<T>::arity() == m_arity &&
            channels::channel_data_type_traits<T>::data_type() == m_dataType ) {
            if( !m_hasCustomFaces ) {
                return const_trimesh3_vertex_channel_accessor<T>( primaryFaces, &m_data, m_arity, m_dataType, false );
            } else if( m_faces.size() == primaryFaces->size() ) {
                return const_trimesh3_vertex_channel_accessor<T>( &m_faces, &m_data, m_arity, m_dataType, true );
            } else {
                throw std::runtime_error(
                    "trimesh3_vertex_channel.get_accessor: The face count of vertex channel \"" +
                    frantic::strings::to_string( m_name ) +
                    "\" was not kept synchronized with that of the primary mesh.  The mesh's face count is " +
                    boost::lexical_cast<std::string>( primaryFaces->size() ) + ", but the channel's face count is " +
                    boost::lexical_cast<std::string>( m_faces.size() ) + "." );
            }
        } else {
            throw std::runtime_error(
                "trimesh3_vertex_channel.get_accessor: Could not access channel \"" +
                frantic::strings::to_string( m_name ) + "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                ", is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }

    template <class T>
    trimesh3_vertex_channel_cvt_accessor<T> get_cvt_accessor( std::vector<vector3>* primaryFaces ) {
        if( channels::channel_data_type_traits<T>::arity() == m_arity ) {
            if( !m_hasCustomFaces ) {
                return trimesh3_vertex_channel_cvt_accessor<T>( primaryFaces, &m_data, m_arity, m_dataType, false,
                                                                m_name );
            } else if( m_faces.size() == primaryFaces->size() ) {
                return trimesh3_vertex_channel_cvt_accessor<T>( &m_faces, &m_data, m_arity, m_dataType, true, m_name );
            } else {
                throw std::runtime_error(
                    "trimesh3_vertex_channel.get_cvt_accessor: The face count of vertex channel \"" +
                    frantic::strings::to_string( m_name ) +
                    "\" was not kept synchronized with that of the primary mesh.  The mesh's face count is " +
                    boost::lexical_cast<std::string>( primaryFaces->size() ) + ", but the channel's face count is " +
                    boost::lexical_cast<std::string>( m_faces.size() ) + "." );
            }
        } else {
            throw std::runtime_error(
                "trimesh3_vertex_channel.get_cvt_accessor: Could not access channel \"" +
                frantic::strings::to_string( m_name ) + "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                ", is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }

    // This returns an accessor to this channel, viewed as the requested type.  When the trimesh3 creates an accessor
    // for someone, it needs to pass a pointer to its faces array.  This is done so that the accessor can assume the
    // faces array is valid for performance reasons.
    template <class T>
    const_trimesh3_vertex_channel_cvt_accessor<T> get_cvt_accessor( const std::vector<vector3>* primaryFaces ) const {
        if( channels::channel_data_type_traits<T>::arity() == m_arity ) {
            if( !m_hasCustomFaces ) {
                return const_trimesh3_vertex_channel_cvt_accessor<T>( primaryFaces, &m_data, m_arity, m_dataType, false,
                                                                      m_name );
            } else if( m_faces.size() == primaryFaces->size() ) {
                return const_trimesh3_vertex_channel_cvt_accessor<T>( &m_faces, &m_data, m_arity, m_dataType, true,
                                                                      m_name );
            } else {
                throw std::runtime_error(
                    "trimesh3_vertex_channel.get_cvt_accessor: The face count of vertex channel \"" +
                    frantic::strings::to_string( m_name ) +
                    "\" was not kept synchronized with that of the primary mesh.  The mesh's face count is " +
                    boost::lexical_cast<std::string>( primaryFaces->size() ) + ", but the channel's face count is " +
                    boost::lexical_cast<std::string>( m_faces.size() ) + "." );
            }
        } else {
            throw std::runtime_error(
                "trimesh3_vertex_channel.get_cvt_accessor: Could not access channel \"" +
                frantic::strings::to_string( m_name ) + "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                ", is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }
};

///////////////////////
// FACE CHANNEL GENERAL ACCESSOR BASE
///////////////////////

// This class is used to provide faster access within a named channel.
template <class DataVector>
class trimesh3_face_channel_general_accessor_base {
  protected:
    DataVector* m_data;

    std::size_t m_arity;         // # of data values for this channel per particle
    data_type_t m_dataType;      // primitive data type
    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    // Function which knows how to do weighted sum combinations of this data type
    frantic::channels::channel_weighted_sum_combine_function_t m_weightedSumFunction;

    trimesh3_face_channel_general_accessor_base( DataVector* data, std::size_t arity, data_type_t dataType )
        : m_data( data )
        , m_arity( arity )
        , m_dataType( dataType ) {
        m_weightedSumFunction = channel_weighted_sum_combine_function( dataType );
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    friend class trimesh3;

  public:
    trimesh3_face_channel_general_accessor_base()
        : m_data( 0 )
        , m_arity( 0 )
        , m_dataType( channels::data_type_int8 )
        , m_primitiveSize( 0 )
        , m_weightedSumFunction( 0 ) {}

    std::size_t size() const { return m_data ? ( m_data->size() / m_primitiveSize ) : 0; }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    bool valid() const { return m_data != 0; }

    operator bool() const { return m_data != 0; }

    frantic::channels::channel_weighted_sum_combine_function_t get_weighted_sum_combine_function() const {
        return m_weightedSumFunction;
    }

    /**
     * This returns the size of the primitive of this channel.  That's equivalent to arity() * (data type size).
     */
    std::size_t primitive_size() const { return m_primitiveSize; }

    void dump( std::ostream& out ) const {
        if( m_data ) {
            out << "m_data size: " << m_data->size() << "\n";
        } else {
            out << "m_data is null\n";
        }
        out << "m_arity: " << m_arity << "\n";
        out << "m_dataType: " << m_dataType << "\n";
        out << "m_dataType str: "
            << frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) << "\n";
        out << "m_primitiveSize: " << m_primitiveSize << "\n";
        out << "size(): " << size() << "\n";
    }

    ///////////
    // Access to the data
    ///////////

    void print( std::ostream& out, std::size_t i ) const {
        channels::channel_data_type_print( out, m_arity, m_dataType, data( i ) );
    }

    const char* data( std::size_t i ) const { return m_data->ptr_at( i * m_primitiveSize ); }
};

///////////////////////
// FACE CHANNEL GENERAL ACCESSOR
///////////////////////

// Forward declaration of the const version of this class, so it has access to the members
class const_trimesh3_face_channel_general_accessor;

// The non-const accessor adds additional methods for writing to the channel
class trimesh3_face_channel_general_accessor
    : public trimesh3_face_channel_general_accessor_base<frantic::graphics::raw_byte_buffer> {
  protected:
    trimesh3_face_channel_general_accessor( frantic::graphics::raw_byte_buffer* data, std::size_t arity,
                                            data_type_t dataType )
        : trimesh3_face_channel_general_accessor_base<frantic::graphics::raw_byte_buffer>( data, arity, dataType ) {}

    friend class const_trimesh3_face_channel_general_accessor;
    friend class trimesh3_face_channel;

  public:
    trimesh3_face_channel_general_accessor() {}

    ///////////
    // Access to the data
    ///////////

    char* data( std::size_t i ) { return m_data->ptr_at( i * m_primitiveSize ); }

    char* add_face() { return m_data->add_element( m_primitiveSize ); }

    void add_faces( size_t numFacesToAdd ) {
        m_data->resize_with_exponential_growth( m_data->size() + numFacesToAdd * m_primitiveSize );
    }

    void set_face_count( size_t faceCount ) { m_data->resize_with_exponential_growth( faceCount * m_primitiveSize ); }

    void reserve_faces( size_t numFacesToReserve ) { m_data->reserve( numFacesToReserve * m_primitiveSize ); }
};

///////////////////////
// CONST FACE CHANNEL GENERAL ACCESSOR
///////////////////////

// The const accessor just uses all the methods in the base accessor class
class const_trimesh3_face_channel_general_accessor
    : public trimesh3_face_channel_general_accessor_base<const frantic::graphics::raw_byte_buffer> {
  protected:
    const_trimesh3_face_channel_general_accessor( const frantic::graphics::raw_byte_buffer* data, std::size_t arity,
                                                  data_type_t dataType )
        : trimesh3_face_channel_general_accessor_base<const frantic::graphics::raw_byte_buffer>( data, arity,
                                                                                                 dataType ) {}

    friend class trimesh3_face_channel;

  public:
    const_trimesh3_face_channel_general_accessor() {}

    const_trimesh3_face_channel_general_accessor( const trimesh3_face_channel_general_accessor& rhs )
        : trimesh3_face_channel_general_accessor_base<const frantic::graphics::raw_byte_buffer>(
              rhs.m_data, rhs.m_arity, rhs.m_dataType ) {}
};

///////////////////////
// FACE CHANNEL ACCESSOR
///////////////////////

// Forward declaration of the const version of this class, so it has access to the members
template <class DataType>
class const_trimesh3_face_channel_accessor;

template <class DataType>
class trimesh3_face_channel_accessor : public trimesh3_face_channel_general_accessor {
  protected:
    trimesh3_face_channel_accessor( frantic::graphics::raw_byte_buffer* data, std::size_t arity, data_type_t dataType )
        : trimesh3_face_channel_general_accessor( data, arity, dataType ) {}

    friend class const_trimesh3_face_channel_accessor<DataType>;
    friend class trimesh3_face_channel;

  public:
    trimesh3_face_channel_accessor() {}

    ///////////
    // Access to the data
    ///////////

    const DataType& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }

    DataType& operator[]( std::size_t i ) {
        return *reinterpret_cast<DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }
};

///////////////////////
// CONST FACE CHANNEL ACCESSOR
///////////////////////

// The const accessor just uses all the methods in the base accessor class
template <class DataType>
class const_trimesh3_face_channel_accessor : public const_trimesh3_face_channel_general_accessor {

    const_trimesh3_face_channel_accessor( const frantic::graphics::raw_byte_buffer* data, std::size_t arity,
                                          channels::data_type_t dataType )
        : const_trimesh3_face_channel_general_accessor( data, arity, dataType ) {}

    friend class trimesh3_face_channel;

  public:
    const_trimesh3_face_channel_accessor() {}

    const_trimesh3_face_channel_accessor( const trimesh3_face_channel_accessor<DataType>& rhs )
        : const_trimesh3_face_channel_general_accessor( rhs.m_data, rhs.arity(), rhs.data_type() ) {}

    const DataType& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const DataType*>( m_data->ptr_at( i * sizeof( DataType ) ) );
    }
};

///////////////////////
// FACE CHANNEL CVT ACCESSOR
///////////////////////

template <class T>
class trimesh3_face_channel_cvt_accessor
    : protected trimesh3_face_channel_general_accessor_base<frantic::graphics::raw_byte_buffer> {
  protected:
    channels::channel_type_convertor_function_t m_convertGet;
    channels::channel_type_convertor_function_t m_convertSet;

    trimesh3_face_channel_cvt_accessor( frantic::graphics::raw_byte_buffer* data, std::size_t sourceDataArity,
                                        data_type_t sourceDataType, const frantic::tstring& channelNameForErrorMessage )
        : trimesh3_face_channel_general_accessor_base<frantic::graphics::raw_byte_buffer>( data, sourceDataArity,
                                                                                           sourceDataType ) {
        data_type_t targetDataType = channels::channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = channels::channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity )
            throw std::runtime_error(
                std::string() + "trimesh3_face_channel_cvt_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );

        m_convertGet =
            get_channel_type_convertor_function( sourceDataType, targetDataType, channelNameForErrorMessage );
        m_convertSet =
            get_channel_type_convertor_function( targetDataType, sourceDataType, channelNameForErrorMessage );
    }

    char* data( std::size_t i ) { return m_data->ptr_at( i * m_primitiveSize ); }

    const char* data( std::size_t i ) const { return m_data->ptr_at( i * m_primitiveSize ); }

    friend class trimesh3_face_channel;

  public:
    typedef T value_type;

    std::size_t size() const {
        return trimesh3_face_channel_general_accessor_base<frantic::graphics::raw_byte_buffer>::size();
    }

    const T get( std::size_t i ) const {
        T result;

        m_convertGet( (char*)&result, data( i ), channels::channel_data_type_traits<T>::arity() );

        return result;
    }

    void set( std::size_t i, const T& val ) {
        m_convertSet( data( i ), (char*)&val, channels::channel_data_type_traits<T>::arity() );
    }
};

///////////////////////
// CONST FACE CHANNEL CVT ACCESSOR
///////////////////////

template <class T>
class const_trimesh3_face_channel_cvt_accessor
    : protected trimesh3_face_channel_general_accessor_base<const frantic::graphics::raw_byte_buffer> {
  protected:
    channels::channel_type_convertor_function_t m_convertGet;

    const_trimesh3_face_channel_cvt_accessor( const frantic::graphics::raw_byte_buffer* data,
                                              std::size_t sourceDataArity, data_type_t sourceDataType,
                                              const frantic::tstring& channelNameForErrorMessage )
        : trimesh3_face_channel_general_accessor_base<const frantic::graphics::raw_byte_buffer>( data, sourceDataArity,
                                                                                                 sourceDataType ) {
        data_type_t targetDataType = channels::channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = channels::channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity )
            throw std::runtime_error(
                std::string() + "const_trimesh3_face_channel_cvt_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );

        m_convertGet =
            get_channel_type_convertor_function( sourceDataType, targetDataType, channelNameForErrorMessage );
    }

    friend class trimesh3_face_channel;

  public:
    typedef T value_type;

    std::size_t size() const {
        return trimesh3_face_channel_general_accessor_base<const frantic::graphics::raw_byte_buffer>::size();
    }

    const T get( std::size_t i ) const {
        T result;

        m_convertGet( (char*)&result, data( i ), channels::channel_data_type_traits<T>::arity() );

        return result;
    }

    const T operator[]( std::size_t i ) const { return get( i ); }
};

class trimesh3_face_channel {
    frantic::tstring m_name; // string name of the channel
    std::size_t m_arity;     // # of data values for this channel per particle
    data_type_t m_dataType;  // primitive data type

    std::size_t m_primitiveSize; // total channel size = arity * the size of m_dataType

    frantic::graphics::raw_byte_buffer m_data;

    trimesh3_face_channel( const frantic::tstring& name, std::size_t arity, data_type_t dataType )
        : m_name( name )
        , m_arity( arity )
        , m_dataType( dataType ) {
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    void set( const frantic::tstring& name, std::size_t arity, data_type_t dataType ) {
        m_name = name;
        m_arity = arity;
        m_dataType = dataType;
        m_primitiveSize = channels::sizeof_channel_data_type( m_dataType ) * m_arity;
    }

    friend class trimesh3;

  public:
    trimesh3_face_channel& operator=( const trimesh3_face_channel& rhs ) {
        m_name = rhs.m_name;
        m_arity = rhs.m_arity;
        m_dataType = rhs.m_dataType;
        m_primitiveSize = rhs.m_primitiveSize;
        m_data = rhs.m_data;

        return *this;
    }

    std::size_t size() const { return m_data.size() / m_primitiveSize; }

    std::size_t arity() const { return m_arity; }

    data_type_t data_type() const { return m_dataType; }

    /**
     * This returns the size of the primitive of this channel.  That's equivalent to arity() * (data type size).
     */
    std::size_t primitive_size() const { return m_primitiveSize; }

    template <typename T>
    bool supports_type() const {
        return channels::channel_data_type_traits<T>::arity() == m_arity &&
               channels::channel_data_type_traits<T>::dataType == m_dataType;
    }

    ///////////
    // Access to the data
    ///////////

    trimesh3_face_channel_general_accessor get_general_accessor() {
        return trimesh3_face_channel_general_accessor( &m_data, m_arity, m_dataType );
    }

    const_trimesh3_face_channel_general_accessor get_general_accessor() const {
        return const_trimesh3_face_channel_general_accessor( &m_data, m_arity, m_dataType );
    }

    // This returns an accessor to this channel, viewed as the requested type.  When the trimesh3 creates an accessor
    // for someone, it needs to pass a pointer to its faces array.  This is done so that the accessor can assume the
    // faces array is valid for performance reasons.
    template <class T>
    trimesh3_face_channel_accessor<T> get_accessor() {
        if( channels::channel_data_type_traits<T>::arity() == m_arity &&
            channels::channel_data_type_traits<T>::data_type() == m_dataType ) {
            return trimesh3_face_channel_accessor<T>( &m_data, m_arity, m_dataType );
        } else {
            throw std::runtime_error( "trimesh3_face_channel.get_accessor: Could not access face channel \"" +
                                      frantic::strings::to_string( m_name ) +
                                      "\", because the requested type is incompatible." );
        }
    }

    // This returns an accessor to this channel, viewed as the requested type.  When the trimesh3 creates an accessor
    // for someone, it needs to pass a pointer to its faces array.  This is done so that the accessor can assume the
    // faces array is valid for performance reasons.
    template <class T>
    const_trimesh3_face_channel_accessor<T> get_accessor() const {
        if( channels::channel_data_type_traits<T>::arity() == m_arity &&
            channels::channel_data_type_traits<T>::data_type() == m_dataType ) {
            return const_trimesh3_face_channel_accessor<T>( &m_data, m_arity, m_dataType );
        } else {
            throw std::runtime_error(
                "trimesh3_face_channel.get_accessor: Could not access face channel \"" +
                frantic::strings::to_string( m_name ) + "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                ", is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }

    // This returns an accessor to this channel, viewed as the requested type.  When the trimesh3 creates an accessor
    // for someone, it needs to pass a pointer to its faces array.  This is done so that the accessor can assume the
    // faces array is valid for performance reasons.
    template <class T>
    trimesh3_face_channel_cvt_accessor<T> get_cvt_accessor() {
        if( channels::channel_data_type_traits<T>::arity() == m_arity ) {
            return trimesh3_face_channel_cvt_accessor<T>( &m_data, m_arity, m_dataType, m_name );
        } else {
            throw std::runtime_error(
                "trimesh3_face_channel.get_cvt_accessor: Could not access face channel \"" +
                frantic::strings::to_string( m_name ) + "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                ", is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }

    // This returns an accessor to this channel, viewed as the requested type.  When the trimesh3 creates an accessor
    // for someone, it needs to pass a pointer to its faces array.  This is done so that the accessor can assume the
    // faces array is valid for performance reasons.
    template <class T>
    const_trimesh3_face_channel_cvt_accessor<T> get_cvt_accessor() const {
        if( channels::channel_data_type_traits<T>::arity() == m_arity ) {
            return const_trimesh3_face_channel_cvt_accessor<T>( &m_data, m_arity, m_dataType, m_name );
        } else {
            throw std::runtime_error(
                "trimesh3_face_channel.get_cvt_accessor: Could not access face channel \"" +
                frantic::strings::to_string( m_name ) + "\", because the requested type, " +
                frantic::strings::to_string( channels::channel_data_type_traits<T>::type_str() ) +
                ", is incompatible with the channel type, " +
                frantic::strings::to_string( channels::channel_data_type_str( m_arity, m_dataType ) ) + "." );
        }
    }
};

} // namespace geometry
} // namespace frantic
