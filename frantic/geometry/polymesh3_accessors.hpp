// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/named_channel_data.hpp>

namespace frantic {
namespace geometry {

typedef int* polymesh3_face_iterator;
typedef const int* polymesh3_const_face_iterator;
typedef std::pair<polymesh3_face_iterator, polymesh3_face_iterator> polymesh3_face_range;
typedef std::pair<polymesh3_const_face_iterator, polymesh3_const_face_iterator> polymesh3_const_face_range;

class polymesh3_vertex_accessor_base {
    std::size_t m_numFaces;
    int* m_pCustomFaces;
    const int* m_pCustomFaceOffsets;

  protected:
    /**
     * This constructor is only available to base classes.
     * @param numFaces The number of custom faces this channel has. If 0, this channel has no custom faces
     *                  and is mapped to each vertex of the mesh. If non-zero, this is the length of the
     *                  pFaceOffset array.
     * @param pCustomFaces An array of integers denoting the indices into the vertex channel of each vertex of
     *                      each polygon in the polymesh3.
     * @param pFaceOffsets An array of integers (one per polygon in the polymesh3) which contains the index
     *                     of the end of the i'th polygon.
     */
    polymesh3_vertex_accessor_base( std::size_t numFaces = 0, int* pCustomFaces = NULL, const int* pFaceOffsets = NULL )
        : m_numFaces( numFaces )
        , m_pCustomFaces( pCustomFaces )
        , m_pCustomFaceOffsets( pFaceOffsets ) {}

  public:
    typedef polymesh3_face_iterator face_iterator;
    typedef polymesh3_const_face_iterator const_face_iterator;
    typedef polymesh3_face_range face_range;
    typedef polymesh3_const_face_range const_face_range;

  public:
    //******* Custom face accessing *******
    /**
     * Returns true if face data was supplied by the polymesh3 when creating this accessor. If true,
     * The vertex channel is mapped per-face by an array of vertex indices that are only valid in
     * this vertex channel (ie. Not related to the polymesh3 'verts' channel in any way). If false,
     * this channel must correspong 1-to-1 with the polymesh3's 'verts' channel.
     * @return true iff there are custom faces mapping this vertex channel.
     */
    bool has_custom_faces() const { return ( m_numFaces != 0 ); }

    /**
     * Returns the number of custom faces mapping this vertex channel. If non-zero, this must match the number
     * of faces in the polymesh3.
     * @return The number of custom faces mapping this vertex channel.
     */
    std::size_t face_count() const { return m_numFaces; }

    /**
     * Returns an iterator pair for the face with the supplied index. These can be used like STL iterators
     * to iterate over the vertex indices in the i'th face.
     * @note Usage example:
     *   polymesh3_vertex_accessor::face_range r = acc.get_face(i);
     *   for( polymesh3_vertex_accessor::face_iterator it = r.first; it != r.second; ++it )
     *		std::cout << "Vertex: " << (*it) << std::endl;
     *
     *   std::size_t polyCount = (r.second - r.first);
     *   for( std::size_t i = 0; i < polyCount; ++i )
     *		std::cout << "Vertex " << i << ": " << r.first[i] << std::endl;
     *
     * @param i The index of the face to retrieve
     * @result An iterator pair denoting the indices associated with the i'th face.
     */
    face_range get_face( std::size_t i ) {
        if( i == 0 )
            return face_range( m_pCustomFaces, m_pCustomFaces + m_pCustomFaceOffsets[0] );
        else
            return face_range( m_pCustomFaces + m_pCustomFaceOffsets[i - 1], m_pCustomFaces + m_pCustomFaceOffsets[i] );
    }

    /**
     * @overload
     */
    const_face_range get_face( std::size_t i ) const {
        if( i == 0 )
            return const_face_range( m_pCustomFaces, m_pCustomFaces + m_pCustomFaceOffsets[0] );
        else
            return const_face_range( m_pCustomFaces + m_pCustomFaceOffsets[i - 1],
                                     m_pCustomFaces + m_pCustomFaceOffsets[i] );
    }

    std::size_t get_face_degree( std::size_t i ) const {
        const_face_range r = get_face( i );
        return r.second - r.first;
    }
};

class polymesh3_const_vertex_accessor_base {
    std::size_t m_numFaces;
    const int* m_pCustomFaces;
    const int* m_pCustomFaceOffsets;

  protected:
    /**
     * This constructor is only available to base classes.
     * @param numFaces The number of custom faces this channel has. If 0, this channel has no custom faces
     *                  and is mapped to each vertex of the mesh. If non-zero, this is the length of the
     *                  pFaceOffset array.
     * @param pCustomFaces An array of integers denoting the indices into the vertex channel of each vertex of
     *                      each polygon in the polymesh3.
     * @param pFaceOffsets An array of integers (one per polygon in the polymesh3) which contains the index
     *                     of the end of the i'th polygon.
     */
    polymesh3_const_vertex_accessor_base( std::size_t numFaces = 0, const int* pCustomFaces = NULL,
                                          const int* pFaceOffsets = NULL )
        : m_numFaces( numFaces )
        , m_pCustomFaces( pCustomFaces )
        , m_pCustomFaceOffsets( pFaceOffsets ) {}

  public:
    typedef polymesh3_const_face_iterator const_face_iterator;
    typedef polymesh3_const_face_range const_face_range;

  public:
    //******* Custom face accessing *******
    /**
     * Returns true if face data was supplied by the polymesh3 when creating this accessor. If true,
     * The vertex channel is mapped per-face by an array of vertex indices that are only valid in
     * this vertex channel (ie. Not related to the polymesh3 'verts' channel in any way). If false,
     * this channel must correspong 1-to-1 with the polymesh3's 'verts' channel.
     * @return true iff there are custom faces mapping this vertex channel.
     */
    bool has_custom_faces() const { return ( m_numFaces != 0 ); }

    /**
     * Returns the number of custom faces mapping this vertex channel. If non-zero, this must match the number
     * of faces in the polymesh3.
     * @return The number of custom faces mapping this vertex channel.
     */
    std::size_t face_count() const { return m_numFaces; }

    /**
     * Returns an iterator pair for the face with the supplied index. These can be used like STL iterators
     * to iterate over the vertex indices in the i'th face.
     * @note Usage example:
     *   polymesh3_const_vertex_accessor::const_face_range r = acc.get_face(i);
     *   for( polymesh3_const_vertex_accessor::const_face_iterator it = r.first; it != r.second; ++it )
     *		std::cout << "Vertex: " << (*it) << std::endl;
     *
     *   std::size_t polyCount = (r.second - r.first);
     *   for( std::size_t i = 0; i < polyCount; ++i )
     *		std::cout << "Vertex " << i << ": " << r.first[i] << std::endl;
     *
     * @param i The index of the face to retrieve
     * @result An iterator pair denoting the indices associated with the i'th face.
     */
    const_face_range get_face( std::size_t i ) const {
        if( i == 0 )
            return const_face_range( m_pCustomFaces, m_pCustomFaces + m_pCustomFaceOffsets[0] );
        else
            return const_face_range( m_pCustomFaces + m_pCustomFaceOffsets[i - 1],
                                     m_pCustomFaces + m_pCustomFaceOffsets[i] );
    }

    std::size_t get_face_degree( std::size_t i ) const {
        const_face_range r = get_face( i );
        return r.second - r.first;
    }
};

/**
 * This class provides a mechanism for accessing the vertices of polymesh3 vertex channel, which
 * may have custom faces.
 * @tparam T The expected type of the underlying vertex channel.
 */
template <typename T>
class polymesh3_vertex_accessor : public polymesh3_vertex_accessor_base {
    friend class polymesh3;

    T* m_pVertexData;
    std::size_t m_numVerts;

    /**
     * This constructor is private, because it is expected that the user will only retrieve an accessor
     * from the owning polymesh3.
     * @param numVerts The number of verts in the vertex array
     * @param pData A pointer to the array of vertex data
     * @param numFaces The number of custom faces this channel has. If 0, this channel has no custom faces
     *                  and is mapped to each vertex of the mesh. If non-zero, this is the length of the
     *                  pFaceOffset array.
     * @param pCustomFaces An array of integers denoting the indices into the vertex channel of each vertex of
     *                      each polygon in the polymesh3.
     * @param pFaceOffsets An array of integers (one per polygon in the polymesh3) which contains the index
     *                     of the end of the i'th polygon.
     */
    polymesh3_vertex_accessor( std::size_t numVerts, char* pData, std::size_t numFaces = 0, int* pCustomFaces = NULL,
                               const int* pFaceOffsets = NULL )
        : polymesh3_vertex_accessor_base( numFaces, pCustomFaces, pFaceOffsets )
        , m_pVertexData( (T*)pData )
        , m_numVerts( numVerts ) {}

  public:
    /**
     * The default constructor puts this accessor into a state where is_valid() returns false.
     * Attempting to access any data from an accessor intialized like this will cause undefined
     * behaviour.
     */
    polymesh3_vertex_accessor()
        : polymesh3_vertex_accessor_base( 0, NULL, NULL )
        , m_pVertexData( NULL )
        , m_numVerts( 0 ) {}

    /**
     * @return true if the vertex data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pVertexData != 0; }

    /**
     * @return The number of vertices in this vertex channel.
     */
    std::size_t vertex_count() const { return m_numVerts; }

    /**
     * @return The i'th vertex in this vertex channel
     */
    T& get_vertex( std::size_t i ) { return m_pVertexData[i]; }

    /**
     * @overload
     */
    const T& get_vertex( std::size_t i ) const { return m_pVertexData[i]; }
};

/**
 * This class provides a mechanism for accessing the vertices of polymesh3 vertex channel, which
 * may have custom faces.
 * @tparam T The expected type of the underlying vertex channel.
 */
template <typename T>
class polymesh3_const_vertex_accessor : public polymesh3_const_vertex_accessor_base {
    friend class polymesh3;

    const T* m_pVertexData;
    std::size_t m_numVerts;

    /**
     * This constructor is private, because it is expected that the user will only retrieve an accessor
     * from the owning polymesh3.
     * @param numVerts The number of verts in the vertex array
     * @param pData A pointer to the array of vertex data
     * @param numFaces The number of custom faces this channel has. If 0, this channel has no custom faces
     *                  and is mapped to each vertex of the mesh. If non-zero, this is the length of the
     *                  pFaceOffset array.
     * @param pCustomFaces An array of integers denoting the indices into the vertex channel of each vertex of
     *                      each polygon in the polymesh3.
     * @param pFaceOffsets An array of integers (one per polygon in the polymesh3) which contains the index
     *                     of the end of the i'th polygon.
     */
    polymesh3_const_vertex_accessor( std::size_t numVerts, const char* pData, std::size_t numFaces = 0,
                                     const int* pCustomFaces = NULL, const int* pFaceOffsets = NULL )
        : polymesh3_const_vertex_accessor_base( numFaces, pCustomFaces, pFaceOffsets )
        , m_pVertexData( (T*)pData )
        , m_numVerts( numVerts ) {}

  public:
    /**
     * The default constructor puts this accessor into a state where is_valid() returns false.
     * Attempting to access any data from an accessor intialized like this will cause undefined
     * behaviour.
     */
    polymesh3_const_vertex_accessor()
        : polymesh3_const_vertex_accessor_base( 0, NULL, NULL )
        , m_pVertexData( NULL )
        , m_numVerts( 0 ) {}

    /**
     * @return true if the vertex data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pVertexData != 0; }

    /**
     * @return The number of vertices in this vertex channel.
     */
    std::size_t vertex_count() const { return m_numVerts; }

    /**
     * @return The i'th vertex in this vertex channel
     */
    const T& get_vertex( std::size_t i ) const { return m_pVertexData[i]; }
};

template <>
class polymesh3_vertex_accessor<void> : public polymesh3_vertex_accessor_base {
    friend class polymesh3;

    char* m_pVertexData;
    std::size_t m_numVerts;

    frantic::channels::data_type_t m_type;
    std::size_t m_arity;
    std::size_t m_elementSize;

    /**
     * This constructor is private, because it is expected that the user will only retrieve an accessor
     * from the owning polymesh3.
     * @param numVerts The number of verts in the vertex array
     * @param type The type of the data stored in this array
     * @param arity The arity of data stored per-element in this array
     * @param pData A pointer to the array of vertex data
     * @param numFaces The number of custom faces this channel has. If 0, this channel has no custom faces
     *                  and is mapped to each vertex of the mesh. If non-zero, this is the length of the
     *                  pFaceOffset array.
     * @param pCustomFaces An array of integers denoting the indices into the vertex channel of each vertex of
     *                      each polygon in the polymesh3.
     * @param pFaceOffset An array of integers (one per polygon in the polymesh3) which contains the index
     *                     of the end of the i'th polygon.
     */
    polymesh3_vertex_accessor( std::size_t numVerts, frantic::channels::data_type_t type, std::size_t arity,
                               char* pData, std::size_t numFaces = 0, int* pCustomFaces = NULL,
                               const int* pFaceOffsets = NULL )
        : polymesh3_vertex_accessor_base( numFaces, pCustomFaces, pFaceOffsets )
        , m_pVertexData( pData )
        , m_numVerts( numVerts )
        , m_type( type )
        , m_arity( arity ) {
        m_elementSize = arity * frantic::channels::sizeof_channel_data_type( type );
    }

  public:
    /**
     * The default constructor puts this accessor into a state where is_valid() returns false.
     * Attempting to access any data from an accessor intialized like this will cause undefined
     * behaviour.
     */
    polymesh3_vertex_accessor()
        : polymesh3_vertex_accessor_base( 0, NULL, NULL )
        , m_pVertexData( NULL )
        , m_numVerts( 0 )
        , m_type( frantic::channels::data_type_invalid )
        , m_arity( 0 ) {}

    /**
     * @return true if the vertex data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pVertexData != 0; }

    /**
     * @return The type of the data stored in this channel
     */
    frantic::channels::data_type_t get_type() const { return m_type; }

    /**
     * @return The arity of the data stored per-element in this channel
     */
    std::size_t get_arity() const { return m_arity; }

    /**
     * @return The number of vertices in this vertex channel.
     */
    std::size_t vertex_count() const { return m_numVerts; }

    /**
     * @return Pointer to the i'th vertex in this vertex channel
     */
    char* get_vertex( std::size_t i ) { return m_pVertexData + ( m_elementSize * i ); }

    /**
     * @overload
     */
    const char* get_vertex( std::size_t i ) const { return m_pVertexData + ( m_elementSize * i ); }
};

template <>
class polymesh3_const_vertex_accessor<void> : public polymesh3_const_vertex_accessor_base {
    friend class polymesh3;

    const char* m_pVertexData;
    std::size_t m_numVerts;

    frantic::channels::data_type_t m_type;
    std::size_t m_arity;
    std::size_t m_elementSize;

    /**
     * This constructor is private, because it is expected that the user will only retrieve an accessor
     * from the owning polymesh3.
     * @param numVerts The number of verts in the vertex array
     * @param type The type of the data stored in this array
     * @param arity The arity of data stored per-element in this array
     * @param pData A pointer to the array of vertex data
     * @param numFaces The number of custom faces this channel has. If 0, this channel has no custom faces
     *                  and is mapped to each vertex of the mesh. If non-zero, this is the length of the
     *                  pFaceOffset array.
     * @param pCustomFaces An array of integers denoting the indices into the vertex channel of each vertex of
     *                      each polygon in the polymesh3.
     * @param pFaceOffset An array of integers (one per polygon in the polymesh3) which contains the index
     *                     of the end of the i'th polygon.
     */
    polymesh3_const_vertex_accessor( std::size_t numVerts, frantic::channels::data_type_t type, std::size_t arity,
                                     const char* pData, std::size_t numFaces = 0, const int* pCustomFaces = NULL,
                                     const int* pFaceOffsets = NULL )
        : polymesh3_const_vertex_accessor_base( numFaces, pCustomFaces, pFaceOffsets )
        , m_pVertexData( pData )
        , m_numVerts( numVerts )
        , m_type( type )
        , m_arity( arity ) {
        m_elementSize = arity * frantic::channels::sizeof_channel_data_type( type );
    }

  public:
    /**
     * The default constructor puts this accessor into a state where is_valid() returns false.
     * Attempting to access any data from an accessor intialized like this will cause undefined
     * behaviour.
     */
    polymesh3_const_vertex_accessor()
        : polymesh3_const_vertex_accessor_base( 0, NULL, NULL )
        , m_pVertexData( NULL )
        , m_numVerts( 0 )
        , m_type( frantic::channels::data_type_invalid )
        , m_arity( 0 ) {}

    /**
     * @return true if the vertex data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pVertexData != 0; }

    /**
     * @return The type of the data stored in this channel
     */
    frantic::channels::data_type_t get_type() const { return m_type; }

    /**
     * @return The arity of the data stored per-element in this channel
     */
    std::size_t get_arity() const { return m_arity; }

    /**
     * @return The number of vertices in this vertex channel.
     */
    std::size_t vertex_count() const { return m_numVerts; }

    /**
     * @return Pointer to the i'th vertex in this vertex channel
     */
    const char* get_vertex( std::size_t i ) const { return m_pVertexData + ( m_elementSize * i ); }
};

template <class T>
class polymesh3_cvt_vertex_accessor : public polymesh3_vertex_accessor_base {
    friend class polymesh3;

  protected:
    char* m_pData;
    std::size_t m_primitiveSize;
    std::size_t m_numVerts;

    frantic::channels::channel_type_convertor_function_t m_convertGet;
    frantic::channels::channel_type_convertor_function_t m_convertSet;

    const char* data( std::size_t i ) const { return m_pData + i * m_primitiveSize; }

    char* data( std::size_t i ) { return m_pData + i * m_primitiveSize; }

    polymesh3_cvt_vertex_accessor( std::size_t numVerts, frantic::channels::data_type_t sourceDataType,
                                   std::size_t sourceDataArity, char* pData,
                                   const frantic::tstring& channelNameForErrorMessage, std::size_t numFaces = 0,
                                   int* pCustomFaces = NULL, const int* pFaceOffset = NULL )
        : polymesh3_vertex_accessor_base( numFaces, pCustomFaces, pFaceOffset )
        , m_pData( pData )
        , m_numVerts( numVerts ) {
        m_primitiveSize = frantic::channels::sizeof_channel_data_type( sourceDataType ) * sourceDataArity;

        frantic::channels::data_type_t targetDataType = frantic::channels::channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = frantic::channels::channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity ) {
            throw std::runtime_error(
                std::string() + "polymesh3_cvt_vertex_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );
        }

        m_convertGet = frantic::channels::get_channel_type_convertor_function( sourceDataType, targetDataType,
                                                                               channelNameForErrorMessage );
        m_convertSet = frantic::channels::get_channel_type_convertor_function( targetDataType, sourceDataType,
                                                                               channelNameForErrorMessage );
    }

  public:
    polymesh3_cvt_vertex_accessor()
        : m_numVerts( 0 )
        , m_pData( 0 )
        , m_primitiveSize( 0 ) {}

    /**
     * @return true if the vertex data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pData != 0; }

    std::size_t vertex_count() const { return m_numVerts; }

    const T get_vertex( std::size_t i ) const {
        T result;
        m_convertGet( (char*)&result, data( i ), frantic::channels::channel_data_type_traits<T>::arity() );
        return result;
    }

    void set_vertex( std::size_t i, const T& val ) {
        m_convertSet( data( i ), (char*)&val, frantic::channels::channel_data_type_traits<T>::arity() );
    }
};

template <class T>
class polymesh3_const_cvt_vertex_accessor : public polymesh3_const_vertex_accessor_base {
    friend class polymesh3;

  protected:
    const char* m_pData;
    std::size_t m_primitiveSize;
    std::size_t m_numVerts;

    frantic::channels::channel_type_convertor_function_t m_convertGet;

    const char* data( std::size_t i ) const { return m_pData + i * m_primitiveSize; }

    polymesh3_const_cvt_vertex_accessor( std::size_t numVerts, frantic::channels::data_type_t sourceDataType,
                                         std::size_t sourceDataArity, const char* pData,
                                         const frantic::tstring& channelNameForErrorMessage, std::size_t numFaces = 0,
                                         const int* pCustomFaces = NULL, const int* pFaceOffset = NULL )
        : polymesh3_const_vertex_accessor_base( numFaces, pCustomFaces, pFaceOffset )
        , m_pData( pData )
        , m_numVerts( numVerts ) {
        m_primitiveSize = frantic::channels::sizeof_channel_data_type( sourceDataType ) * sourceDataArity;

        frantic::channels::data_type_t targetDataType = frantic::channels::channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = frantic::channels::channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity ) {
            throw std::runtime_error(
                std::string() + "polymesh3_cvt_vertex_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );
        }

        m_convertGet = frantic::channels::get_channel_type_convertor_function( sourceDataType, targetDataType,
                                                                               channelNameForErrorMessage );
    }

  public:
    polymesh3_const_cvt_vertex_accessor()
        : m_pData( 0 )
        , m_primitiveSize( 0 )
        , m_numVerts( 0 ) {}

    /**
     * @return true if the vertex data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pData != 0; }

    std::size_t vertex_count() const { return m_numVerts; }

    const T get_vertex( std::size_t i ) const {
        T result;
        m_convertGet( (char*)&result, data( i ), frantic::channels::channel_data_type_traits<T>::arity() );
        return result;
    }
};

/**
 * This class provides a mechanism for accessing the data of polymesh3 face channel.
 * @tparam T The expected type of the underlying face channel.
 */
template <typename T>
class polymesh3_face_accessor {
    friend class polymesh3;

    T* m_pData;
    std::size_t m_numFaces;

    /**
     * This constructor is private, because it is expected that the user will only retrieve an accessor
     * from the owning polymesh3.
     * @param numFaces The number of faces in the data array
     * @param pData A pointer to the array of face data
     */
    polymesh3_face_accessor( std::size_t numFaces, char* pData )
        : m_pData( (T*)pData )
        , m_numFaces( numFaces ) {}

  public:
    /**
     * The default constructor puts this accessor into a state where is_valid() returns false.
     * Attempting to access any data from an accessor intialized like this will cause undefined
     * behaviour.
     */
    polymesh3_face_accessor()
        : m_pData( NULL )
        , m_numFaces( 0 ) {}

    /**
     * @return true if the face data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pData != 0; }

    /**
     * @return The number of faces in this face channel.
     */
    std::size_t face_count() const { return m_numFaces; }

    /**
     * @return The i'th face in this face channel
     */
    T& get_face( std::size_t i ) { return m_pData[i]; }

    /**
     * @overload
     */
    const T& get_face( std::size_t i ) const { return m_pData[i]; }
};

/**
 * This class provides a mechanism for accessing the data of polymesh3 face channel.
 * @tparam T The expected type of the underlying face channel.
 */
template <typename T>
class polymesh3_const_face_accessor {
    friend class polymesh3;

    const T* m_pData;
    std::size_t m_numFaces;

    /**
     * This constructor is private, because it is expected that the user will only retrieve an accessor
     * from the owning polymesh3.
     * @param numFaces The number of faces in the data array
     * @param pData A pointer to the array of face data
     */
    polymesh3_const_face_accessor( std::size_t numFaces, const char* pData )
        : m_pData( (T*)pData )
        , m_numFaces( numFaces ) {}

  public:
    /**
     * The default constructor puts this accessor into a state where is_valid() returns false.
     * Attempting to access any data from an accessor intialized like this will cause undefined
     * behaviour.
     */
    polymesh3_const_face_accessor()
        : m_pData( NULL )
        , m_numFaces( 0 ) {}

    /**
     * @return true if the face data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pData != 0; }

    /**
     * @return The number of faces in this face channel.
     */
    std::size_t face_count() const { return m_numFaces; }

    /**
     * @return The i'th face in this face channel
     */
    const T& get_face( std::size_t i ) const { return m_pData[i]; }
};

template <>
class polymesh3_face_accessor<void> {
    friend class polymesh3;

    char* m_pData;
    std::size_t m_numFaces;

    frantic::channels::data_type_t m_type;
    std::size_t m_arity;
    std::size_t m_elementSize;

    /**
     * This constructor is private, because it is expected that the user will only retrieve an accessor
     * from the owning polymesh3.
     * @param numFaces The number of faces in the data array
     * @param type The type of the data stored in this array
     * @param arity The arity of data stored per-element in this array
     * @param pData A pointer to the array of face data
     */
    polymesh3_face_accessor( std::size_t numFaces, frantic::channels::data_type_t type, std::size_t arity, char* pData )
        : m_pData( pData )
        , m_numFaces( numFaces )
        , m_type( type )
        , m_arity( arity ) {
        m_elementSize = arity * frantic::channels::sizeof_channel_data_type( type );
    }

  public:
    /**
     * The default constructor puts this accessor into a state where is_valid() returns false.
     * Attempting to access any data from an accessor intialized like this will cause undefined
     * behaviour.
     */
    polymesh3_face_accessor()
        : m_pData( NULL )
        , m_numFaces( 0 )
        , m_type( frantic::channels::data_type_invalid )
        , m_arity( 0 ) {}

    /**
     * @return true if the face data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pData != 0; }

    /**
     * @return The type of the data stored in this channel
     */
    frantic::channels::data_type_t get_type() const { return m_type; }

    /**
     * @return The arity of the data stored per-element in this channel
     */
    std::size_t get_arity() const { return m_arity; }

    /**
     * @return The number of faces in this face channel.
     */
    std::size_t face_count() const { return m_numFaces; }

    /**
     * @return Pointer to the i'th face in this face channel
     */
    char* get_face( std::size_t i ) { return m_pData + ( m_elementSize * i ); }

    /**
     * @overload
     */
    const char* get_face( std::size_t i ) const { return m_pData + ( m_elementSize * i ); }
};

template <>
class polymesh3_const_face_accessor<void> {
    friend class polymesh3;

    const char* m_pData;
    std::size_t m_numFaces;

    frantic::channels::data_type_t m_type;
    std::size_t m_arity;
    std::size_t m_elementSize;

    /**
     * This constructor is private, because it is expected that the user will only retrieve an accessor
     * from the owning polymesh3.
     * @param numFaces The number of faces in the data array
     * @param type The type of the data stored in this array
     * @param arity The arity of data stored per-element in this array
     * @param pData A pointer to the array of face data
     */
    polymesh3_const_face_accessor( std::size_t numFaces, frantic::channels::data_type_t type, std::size_t arity,
                                   const char* pData )
        : m_pData( pData )
        , m_numFaces( numFaces )
        , m_type( type )
        , m_arity( arity ) {
        m_elementSize = arity * frantic::channels::sizeof_channel_data_type( type );
    }

  public:
    /**
     * The default constructor puts this accessor into a state where is_valid() returns false.
     * Attempting to access any data from an accessor intialized like this will cause undefined
     * behaviour.
     */
    polymesh3_const_face_accessor()
        : m_pData( NULL )
        , m_numFaces( 0 )
        , m_type( frantic::channels::data_type_invalid )
        , m_arity( 0 ) {}

    /**
     * @return true if the face data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pData != 0; }

    /**
     * @return The type of the data stored in this channel
     */
    frantic::channels::data_type_t get_type() const { return m_type; }

    /**
     * @return The arity of the data stored per-element in this channel
     */
    std::size_t get_arity() const { return m_arity; }

    /**
     * @return The number of faces in this face channel.
     */
    std::size_t face_count() const { return m_numFaces; }

    /**
     * @return Pointer to the i'th face in this face channel
     */
    const char* get_face( std::size_t i ) const { return m_pData + ( m_elementSize * i ); }
};

template <class T>
class polymesh3_cvt_face_accessor {
    friend class polymesh3;

  protected:
    char* m_pData;
    std::size_t m_primitiveSize;
    frantic::channels::channel_type_convertor_function_t m_convertGet;
    frantic::channels::channel_type_convertor_function_t m_convertSet;
    std::size_t m_numFaces;

    const char* data( std::size_t i ) const { return m_pData + i * m_primitiveSize; }

    char* data( std::size_t i ) { return m_pData + i * m_primitiveSize; }

    polymesh3_cvt_face_accessor( std::size_t numFaces, frantic::channels::data_type_t sourceDataType,
                                 std::size_t sourceDataArity, char* pData,
                                 const frantic::tstring& channelNameForErrorMessage )
        : m_pData( pData )
        , m_numFaces( numFaces ) {
        m_primitiveSize = frantic::channels::sizeof_channel_data_type( sourceDataType ) * sourceDataArity;

        frantic::channels::data_type_t targetDataType = frantic::channels::channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = frantic::channels::channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity ) {
            throw std::runtime_error(
                std::string() + "polymesh3_cvt_face_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );
        }

        m_convertGet = frantic::channels::get_channel_type_convertor_function( sourceDataType, targetDataType,
                                                                               channelNameForErrorMessage );
        m_convertSet = frantic::channels::get_channel_type_convertor_function( targetDataType, sourceDataType,
                                                                               channelNameForErrorMessage );
    }

  public:
    polymesh3_cvt_face_accessor()
        : m_pData( 0 )
        , m_primitiveSize( 0 )
        , m_numFaces( 0 ) {}

    /**
     * @return true if the face data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pData != 0; }

    /**
     * @return The number of faces in this face channel.
     */
    std::size_t face_count() const { return m_numFaces; }

    const T get_face( std::size_t i ) const {
        T result;
        m_convertGet( (char*)&result, data( i ), frantic::channels::channel_data_type_traits<T>::arity() );
        return result;
    }

    void set_face( std::size_t i, const T& val ) {
        m_convertSet( data( i ), (char*)&val, frantic::channels::channel_data_type_traits<T>::arity() );
    }
};

template <class T>
class polymesh3_const_cvt_face_accessor {
    friend class polymesh3;

  protected:
    const char* m_pData;
    std::size_t m_primitiveSize;
    frantic::channels::channel_type_convertor_function_t m_convertGet;
    std::size_t m_numFaces;

    const char* data( std::size_t i ) const { return m_pData + i * m_primitiveSize; }

    polymesh3_const_cvt_face_accessor( std::size_t numFaces, frantic::channels::data_type_t sourceDataType,
                                       std::size_t sourceDataArity, const char* pData,
                                       const frantic::tstring& channelNameForErrorMessage )
        : m_pData( pData )
        , m_numFaces( numFaces ) {
        m_primitiveSize = frantic::channels::sizeof_channel_data_type( sourceDataType ) * sourceDataArity;

        frantic::channels::data_type_t targetDataType = frantic::channels::channel_data_type_traits<T>::data_type();
        std::size_t targetDataArity = frantic::channels::channel_data_type_traits<T>::arity();

        if( sourceDataArity != targetDataArity ) {
            throw std::runtime_error(
                std::string() + "polymesh3_cvt_face_accessor() - Cannot convert channel: \"" +
                frantic::strings::to_string( channelNameForErrorMessage ) + "\" from type: " +
                frantic::strings::to_string( channel_data_type_str( sourceDataArity, sourceDataType ) ) + " to " +
                frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to differing arity" );
        }

        m_convertGet = frantic::channels::get_channel_type_convertor_function( sourceDataType, targetDataType,
                                                                               channelNameForErrorMessage );
    }

  public:
    polymesh3_const_cvt_face_accessor()
        : m_pData( 0 )
        , m_primitiveSize( 0 )
        , m_numFaces( 0 ) {}

    /**
     * @return true if the face data is valid.  If false, data access will cause undefined
     * behavior.
     */
    bool is_valid() const { return m_pData != 0; }

    /**
     * @return The number of faces in this face channel.
     */
    std::size_t face_count() const { return m_numFaces; }

    const T get_face( std::size_t i ) const {
        T result;
        m_convertGet( (char*)&result, data( i ), frantic::channels::channel_data_type_traits<T>::arity() );
        return result;
    }
};
} // namespace geometry
} // namespace frantic
