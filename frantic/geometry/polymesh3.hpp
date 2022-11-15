// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <map>
#include <vector>

#include <boost/intrusive_ptr.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/channels/named_channel_data.hpp>
#include <frantic/geometry/polymesh3_accessors.hpp>
#include <frantic/geometry/xmesh_metadata.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/graphics/vector3f.hpp>

#include <frantic/logging/progress_logger.hpp>

namespace frantic {
namespace geometry {

// forward declaration
class mesh_interface;
class polymesh3;
class xmesh_reader;

/**
 *  Holds a shared reference to polymesh3 channel data.
 *
 *  This class uses copy-on-write.  If you invoke get_writable(), and there is
 * more than one reference to the shared data, then we make a unique copy of
 * the shared data before returning it.
 *
 *  This is only intended to be used by polymesh3, and users should rarely
 * interact with it.
 */
class polymesh3_channel_data {
  public:
    polymesh3_channel_data();

    // Steals data.  Caller's data is empty after this call.
    polymesh3_channel_data( frantic::graphics::raw_byte_buffer& data, frantic::channels::data_type_t type,
                            std::size_t arity );

    // The number of primitives per-element
    std::size_t arity() const;

    // The total size of an element
    std::size_t element_size() const;

    // The primitive type of the element
    frantic::channels::data_type_t type() const;

    const frantic::graphics::raw_byte_buffer& get() const;
    frantic::graphics::raw_byte_buffer& get_writable();

    bool is_valid() const;

    bool is_shared() const;

  private:
    void assert_data() const;

    boost::shared_ptr<frantic::graphics::raw_byte_buffer> m_pData;

    std::size_t m_arity;
    std::size_t m_elementSize;
    frantic::channels::data_type_t m_type;
};

/**
 *  Holds a shared reference to polymesh3 face index data.
 *
 *  This class uses copy-on-write.  If you invoke get_writable(), and there is
 * more than one reference to the shared data, then we make a unique copy of
 * the shared data before returning it.
 *
 *  This is only intended to be used by polymesh3, and users should rarely
 * interact with it.
 */
class polymesh3_channel_faces {
  public:
    polymesh3_channel_faces();
    // Steals faces.  Caller's faces is empty after this call.
    polymesh3_channel_faces( std::vector<int>& faces );

    const std::vector<int>& get() const;
    std::vector<int>& get_writable();

    bool is_valid() const;

    bool is_shared() const;

  private:
    void assert_data() const;

    boost::shared_ptr<std::vector<int>> m_pFaces;
};

/**
 * This struct holds the information required for a single channel in a polymesh3 object.
 * It is only intended to be used by a polymesh3, and users should rarely interact with it.
 */
class polymesh3_channel {
  public:
    // True if this is a a vertex channel, false if it is a face channel.
    bool is_vertex_channel() const;

    // The number of primitives per-element
    std::size_t arity() const;

    // The total size of an element
    std::size_t element_size() const;

    // The primitive type of the element
    frantic::channels::data_type_t type() const;

    // Storage for channel data
    const frantic::graphics::raw_byte_buffer& get_data() const;
    frantic::graphics::raw_byte_buffer& get_writable_data();

    // Only used by vertex channels, and iff they have custom faces.
    const std::vector<int>& get_faces() const;
    std::vector<int>& get_writable_faces();

    // Return the number of elements in the channel
    std::size_t size() const { return data.get().size() / element_size(); }

  private:
    friend class polymesh3;

    polymesh3_channel_data data;
    polymesh3_channel_faces faces;

    bool m_isVertexChannel;
};

/**
 * This class stores a polygonal mesh object, and all the auxillary data associated with it.
 * It heavily utilizes the named channel systems from the frantic::channels namespace.
 * This implementation should always be created on the heap and used with a boost::intrusive_ptr
 * in order to achieve optimal performance, while also easily sharing this data structure with
 * multiple parent objects.
 *
 * @note Use the typedef polymesh3_ptr to correctly and easily store a polymesh3 object. This will
 *        automatically handle the reference counting.
 *
 * @note Usage Example:
 *			polymesh3_ptr pMesh = new polymesh3;
 */
class polymesh3 {
  public:
    typedef boost::intrusive_ptr<polymesh3> ptr_type;
    typedef boost::intrusive_ptr<const polymesh3> const_ptr_type;
    typedef std::map<frantic::tstring, polymesh3_channel>::const_iterator iterator;

  private:
    // This is the number of objects referencing the polymesh3. When it goes to 0, this object will be
    // deleted.
    mutable int m_refCount;

    // This is shared by all vertex channels with custom faces. This array stores the end
    // index for the I'th polygon in element I. You can figure out the range of indices associated
    // with polygon I be reading element I-1 and I, which describe the start and end of the polygon
    // respectively.
    const std::vector<int>* m_pFaceEndOffsets;
    // Hold a reference to the faceEndOffsets so it doesn't get deallocated
    // on us.
    polymesh3_channel_faces m_faceEndOffsetsHoldRef;

    // Reserved channels are:
    //  verts <-- Stores the geometry vertices as float32[3], and the geometry faces as the custom faces.
    std::map<frantic::tstring, polymesh3_channel> m_channels;

    // Store an accessor into the geometry channel in order to quickly access information
    // such as the number of vertices.
    frantic::geometry::polymesh3_channel* m_pGeomChannel;

  private:
    polymesh3( polymesh3_channel_data vertData, polymesh3_channel_faces geomPolyIndices,
               polymesh3_channel_faces geomPolyEndIndices );

    polymesh3( frantic::graphics::raw_byte_buffer& vertData, std::vector<int>& geomPolyIndices,
               std::vector<int>& geomPolyEndIndices );

    friend void intrusive_ptr_add_ref( const polymesh3* pPolyMesh );
    friend void intrusive_ptr_release( const polymesh3* pPolyMesh );
    friend class polymesh3_builder;
    friend class const_shared_polymesh3_builder;

    // TODO: Do these need to be friends?
    friend ptr_type load_xmesh_polymesh_file( const frantic::geometry::xmesh_reader& reader,
                                              const frantic::channels::channel_propagation_policy& cpp,
                                              bool loadFaces );
    friend void write_xmesh_polymesh_file( const frantic::tstring& path, const ptr_type& polymesh,
                                           const frantic::geometry::xmesh_metadata* metadata );

    void acquire() const;
    void release() const;

    void init( polymesh3_channel_data vertData, polymesh3_channel_faces geomPolyIndices,
               polymesh3_channel_faces geomPolyEndIndices );

    void add_vertex_channel( const frantic::tstring& channel, polymesh3_channel_data& inVertexBuffer,
                             polymesh3_channel_faces* pInFaceBuffer = (polymesh3_channel_faces*)NULL );

    void add_face_channel( const frantic::tstring& channel, polymesh3_channel_data& inFaceBuffer );

    polymesh3_channel& get_vertex_channel( const frantic::tstring& channel );
    const polymesh3_channel& get_const_vertex_channel( const frantic::tstring& channel ) const;

    template <typename T>
    polymesh3_channel& get_vertex_channel( const frantic::tstring& channel ) {
        polymesh3_channel& ch = get_vertex_channel( channel );
        if( ch.type() != frantic::channels::channel_data_type_traits<T>::data_type() ||
            ch.arity() != frantic::channels::channel_data_type_traits<T>::arity() ) {
            throw std::runtime_error(
                "polymesh3::get_vertex_channel() The channel named: \"" + frantic::strings::to_string( channel ) +
                "\" did not have the type expected. Found: " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str( ch.arity(), ch.type() ) ) +
                ", expected: " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) );
        }
        return ch;
    }

    template <typename T>
    const polymesh3_channel& get_const_vertex_channel( const frantic::tstring& channel ) const {
        const polymesh3_channel& ch = get_const_vertex_channel( channel );
        if( ch.type() != frantic::channels::channel_data_type_traits<T>::data_type() ||
            ch.arity() != frantic::channels::channel_data_type_traits<T>::arity() ) {
            throw std::runtime_error(
                "polymesh3::get_const_vertex_channel() The channel named: \"" + frantic::strings::to_string( channel ) +
                "\" did not have the type expected. Found: " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str( ch.arity(), ch.type() ) ) +
                ", expected: " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) );
        }
        return ch;
    }

    polymesh3_channel& get_face_channel( const frantic::tstring& channel );
    const polymesh3_channel& get_const_face_channel( const frantic::tstring& channel ) const;

    template <typename T>
    polymesh3_channel& get_face_channel( const frantic::tstring& channel ) {
        polymesh3_channel& ch = get_face_channel( channel );
        if( ch.type() != frantic::channels::channel_data_type_traits<T>::data_type() ||
            ch.arity() != frantic::channels::channel_data_type_traits<T>::arity() ) {
            throw std::runtime_error(
                "polymesh3::get_face_channel() The channel named: \"" + frantic::strings::to_string( channel ) +
                "\" did not have the type expected. Found: " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str( ch.arity(), ch.type() ) ) +
                ", expected: " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) );
        }
        return ch;
    }

    template <typename T>
    const polymesh3_channel& get_const_face_channel( const frantic::tstring& channel ) const {
        const polymesh3_channel& ch = get_const_face_channel( channel );
        if( ch.type() != frantic::channels::channel_data_type_traits<T>::data_type() ||
            ch.arity() != frantic::channels::channel_data_type_traits<T>::arity() ) {
            throw std::runtime_error(
                "polymesh3::get_const_face_channel() The channel named: \"" + frantic::strings::to_string( channel ) +
                "\" did not have the type expected. Found: " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str( ch.arity(), ch.type() ) ) +
                ", expected: " +
                frantic::strings::to_string( frantic::channels::channel_data_type_traits<T>::type_str() ) );
        }
        return ch;
    }

  public:
    virtual ~polymesh3() {}

    /**
     * @return The number of float32[3] geometric vertices in this polymesh3. A cube with 6 quads will return 8.
     */
    std::size_t vertex_count() const { return m_pGeomChannel->size(); }

    /**
     * @return the number of polygons stored in this polymesh3. A cube with 6 quads will return 6.
     */
    std::size_t face_count() const { return m_pFaceEndOffsets->size(); }

    /**
     * @return the number of vertices in each face, summed over all faces in the mesh.
     */
    std::size_t face_vertex_count() const;

    /**
     * @return the i'th geometric vertex in this polymesh3. This is a convenience function.
     */
    frantic::graphics::vector3f get_vertex( std::size_t i ) const {
        using frantic::graphics::vector3f;
        return reinterpret_cast<const vector3f*>( m_pGeomChannel->get_data().ptr_at( 0 ) )[i];
    }

    /**
     * @return An iterator over the list of channels in this polymesh3. It does not allow changes to the channels,
     *          only inspection of various channel properties.
     */
    iterator begin() const { return m_channels.begin(); }

    /**
     * @return An iterator one past the end of the list of channels in this polymesh3
     */
    iterator end() const { return m_channels.end(); }

    /**
     * @return the number of channels in this polymesh
     */
    size_t channel_count() const { return m_channels.size(); }

    /**
     * @param outNames output vector where all vertex channel names will be stored.  The reserved channel 'verts' will
     * be omitted.
     */
    void get_vertex_channel_names( std::vector<frantic::tstring>& outNames ) const;

    /**
     * @param outNames output vector where all face channel names will be stored.
     */
    void get_face_channel_names( std::vector<frantic::tstring>& outNames ) const;

    /**
     * @param outOffsets output vector where all face end points will be stored
     */
    void get_face_end_offsets( std::vector<int>& outOffsets ) const;

    /**
     * Return true iff a channel with the given name exists
     * @return true if the specified channel exists, false otherwise
     */
    bool has_channel( const frantic::tstring& name ) const { return ( m_channels.find( name ) != m_channels.end() ); }

    /**
     * Return true iff a vertex channel with the given name exists
     * @return true if the specified channel exists and it is a vertex channel, false otherwise
     */
    bool has_vertex_channel( const frantic::tstring& name ) const;

    /**
     * Return true iff a face channel with the given name exists
     * @return true if the specified channel exists and it is a face channel, false otherwise
     */
    bool has_face_channel( const frantic::tstring& name ) const;

    /**
     * This function retrieves the information about the channel with name 'name'
     * @param name The name of the channel to query about
     * @return The data structure defining the channel with name 'name'
     */
    const polymesh3_channel& get_channel_info( const frantic::tstring& name ) const {
        std::map<frantic::tstring, polymesh3_channel>::const_iterator it = m_channels.find( name );
        if( it == m_channels.end() )
            throw std::runtime_error( "polymesh3::get_channel_info() There is no channel named: \"" +
                                      frantic::strings::to_string( name ) + "\"" );
        return it->second;
    }

    /**
     * This function will return true iff every polygon in this mesh is a triangle.
     * @return True if this mesh only contains triangles. False otherwise.
     */
    bool is_triangle_mesh() const;

    /**
     * This function will add a new vertex channel to the polymesh. It will throw an exception if custom faces are
     * provided but the array length doesn't match the 'verts:Face' array length, or if no custom faces are provided
     * but the number of channel vertices doesn't match the 'verts' channel.
     * @param channel The name of the new channel to add. Must not already exist.
     * @param type The per-data item type stored in this vertex channel.
     * @param arity The number of primitives per vertex element.
     * @param inVertexBuffer A buffer storing the vertex data for this channel. This buffer is stolen from the supplied
     *                        raw_byte_buffer. (ie. after this call, inVertexBuffer.size() == 0)
     * @param pInFaceBuffer An optional buffer of int32[1] that represent the custom faces to map this vertex channel.
     *                       pInFaceBuffer->size()/4 must match face_vertex_count(), or be NULL. This buffer is stolen
     *                       from the supplied raw_byte_buffer. (ie. after this call, pInFaceBuffer->size() == 0)
     */
    void add_vertex_channel( const frantic::tstring& channel, frantic::channels::data_type_t type, std::size_t arity,
                             frantic::graphics::raw_byte_buffer& inVertexBuffer,
                             std::vector<int>* pInFaceBuffer = NULL );
    template <typename T, int N>
    void add_vertex_channel( const frantic::tstring& channel, const T ( &inVertexData )[N],
                             std::vector<int>* pInFaceBuffer = NULL ) {
        graphics::raw_byte_buffer buffer( inVertexData, sizeof( T[N] ) );
        add_vertex_channel( channel, channels::channel_data_type_traits<T>::data_type(),
                            channels::channel_data_type_traits<T>::arity(), buffer, pInFaceBuffer );
    }

    /**
     * This function will add a new face channel to the polymesh. It will throw an exception if the number of entries
     * in 'inFaceBuffer' does not match the number of polygons in the mesh. The supplied buffer is stolen.
     * @param channel The name of the new channel to add. Must not already exist.
     * @param type The per-data item type stored in this channel.
     * @param arity The number of primitives per element.
     * @param inFaceBuffer A buffer storing the face data for this channel. This buffer is stolen from the supplied
     *                        raw_byte_buffer. (ie. after this call, 'inFaceBuffer'.size() == 0)
     */
    void add_face_channel( const frantic::tstring& channel, frantic::channels::data_type_t type, std::size_t arity,
                           frantic::graphics::raw_byte_buffer& inFaceBuffer );
    template <typename T, int N>
    void add_face_channel( const frantic::tstring& channel, const T ( &inFaceData )[N] ) {
        graphics::raw_byte_buffer buffer( inFaceData, sizeof( T[N] ) );
        add_face_channel( channel, channels::channel_data_type_traits<T>::data_type(),
                          channels::channel_data_type_traits<T>::arity(), buffer );
    }

    /**
     * This function will add a new vertex channel to the polymesh, leaving all elements uninitialized. The new
     * channel will NOT have custom faces, therefore there will be one vertex for each vertex in the 'verts' channel.
     * @param channel The name of the new channel
     * @param type The type of the new channel
     * @param arity The number of primitive data elements per vertex element
     */
    void add_empty_vertex_channel( const frantic::tstring& channel, frantic::channels::data_type_t type,
                                   std::size_t arity );

    /**
     * This function will add a new vertex channel to the polymesh, leaving all elements uninitialized. The new
     * channel WILL have custom faces, which will be all initialized to 0. The vertex channel will be allocated
     * for 'numVertices' elements, but they content of each will not be defined. The custom face channel will
     * be allocated to the same size as the 'verts:Faces' channel.
     * @param channel The name of the new channel
     * @param type The type of the new channel
     * @param arity The number of primitive data elements per vertex element
     * @param numVertices The number of vertex elements for this channel
     */
    void add_empty_vertex_channel( const frantic::tstring& channel, frantic::channels::data_type_t type,
                                   std::size_t arity, std::size_t numVertices );

    /**
     * This function will add a new face channel to the polymesh, leaving all elements uninitialized.
     * @param channel The name of the new channel
     * @param type The type of the new channel
     * @param arity The number of primitive data elements per face element
     */
    void add_empty_face_channel( const frantic::tstring& channel, frantic::channels::data_type_t type,
                                 std::size_t arity );

    /**
     * This function will remove an existing vertex channel from the polymesh.
     * @param channel The name of the channel to remove
     */
    void erase_vertex_channel( const frantic::tstring& channel );

    /**
     * This function will remove an existing face channel from the polymesh.
     * @param channel The name of the channel to remove
     */
    void erase_face_channel( const frantic::tstring& channel );

    /**
     * This function will copy an existing channel into a different channel, creating a channel
     * with the name 'destChannel' if one does not already exist.
     * @param destChannel The name of the destination channel. Will be created if it does not already exist.
     * @param srcChannel The name of the source channel.
     */
    void copy_channel( const frantic::tstring& destChannel, const frantic::tstring& srcChannel );

    // Builds the vertex normals based on incident angle.
    void build_vertex_normals( const frantic::tstring& vertexChannelName = _T("Normal"), bool normalizeNormals = true );

    /**
     *  Build a channel which holds the face normals.
     *
     *  The resulting normals are usually unit vectors.  However, if a face has
     * coincident vertices, then its normal will have magnitude 0.
     *
     * @param faceChannelName the name of a float32[3] face channel to create
     *		and fill with the mesh's face normals.
     */
    void build_face_normals( const frantic::tstring& faceChannelName = _T("FaceNormal") );

    /**
     * This function will create and return an accessor object capable of reading the given vertex channel
     * and interpreting its data.
     * @tparam T The type of data expected to be in the channel
     * @param  channel The name of the vertex channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    template <typename T>
    polymesh3_vertex_accessor<T> get_vertex_accessor( const frantic::tstring& channel ) {
        polymesh3_channel& ch = get_vertex_channel<T>( channel );

        if( ch.get_faces().size() == 0 ) {
            return polymesh3_vertex_accessor<T>( ch.size(), ch.get_writable_data().ptr_at( 0 ) );
        } else {
            return polymesh3_vertex_accessor<T>( ch.size(), ch.get_writable_data().ptr_at( 0 ),
                                                 m_pFaceEndOffsets->size(), &ch.get_writable_faces()[0],
                                                 &( *m_pFaceEndOffsets )[0] );
        }
    }

    /**
     * Create and return an accessor object capable of reading the given vertex channel
     * and interpreting its data.
     * @tparam T The type of data expected to be in the channel
     * @param  channel The name of the vertex channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    template <typename T>
    polymesh3_const_vertex_accessor<T> get_const_vertex_accessor( const frantic::tstring& channel ) const {
        const polymesh3_channel& ch = get_const_vertex_channel<T>( channel );

        if( ch.get_faces().size() == 0 ) {
            return polymesh3_const_vertex_accessor<T>( ch.size(), ch.get_data().ptr_at( 0 ) );
        } else {
            return polymesh3_const_vertex_accessor<T>( ch.size(), ch.get_data().ptr_at( 0 ), m_pFaceEndOffsets->size(),
                                                       &ch.get_faces()[0], &( *m_pFaceEndOffsets )[0] );
        }
    }

    /**
     * This function will create and return an accessor object capable of reading the given vertex channel
     * and interpreting its data, but without knowing the expected type at compile time.
     * @param  channel The name of the vertex channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    polymesh3_vertex_accessor<void> get_vertex_accessor( const frantic::tstring& channel ) {
        polymesh3_channel& ch = get_vertex_channel( channel );

        if( ch.get_faces().size() == 0 ) {
            return polymesh3_vertex_accessor<void>( ch.size(), ch.type(), ch.arity(),
                                                    ch.get_writable_data().ptr_at( 0 ) );
        } else {
            return polymesh3_vertex_accessor<void>( ch.size(), ch.type(), ch.arity(),
                                                    ch.get_writable_data().ptr_at( 0 ), m_pFaceEndOffsets->size(),
                                                    &ch.get_writable_faces()[0], &( *m_pFaceEndOffsets )[0] );
        }
    }

    /**
     * This function will create and return an accessor object capable of reading the given vertex channel
     * and interpreting its data, but without knowing the expected type at compile time.
     * @param  channel The name of the vertex channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    polymesh3_const_vertex_accessor<void> get_const_vertex_accessor( const frantic::tstring& channel ) const {
        const polymesh3_channel& ch = get_const_vertex_channel( channel );

        if( ch.get_faces().size() == 0 ) {
            return polymesh3_const_vertex_accessor<void>( ch.size(), ch.type(), ch.arity(), ch.get_data().ptr_at( 0 ) );
        } else {
            return polymesh3_const_vertex_accessor<void>( ch.size(), ch.type(), ch.arity(), ch.get_data().ptr_at( 0 ),
                                                          m_pFaceEndOffsets->size(), &ch.get_faces()[0],
                                                          &( *m_pFaceEndOffsets )[0] );
        }
    }

    /**
     * This function will create and return an accessor object capable of reading the given vertex channel
     * and interpreting its data, automatically converting it to the requested type.
     * @tparam T The type to convert the channel's data to.
     * @param  channel The name of the vertex channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    template <typename T>
    polymesh3_cvt_vertex_accessor<T> get_cvt_vertex_accessor( const frantic::tstring& channel ) {
        polymesh3_channel& ch = get_vertex_channel( channel );

        if( ch.get_faces().size() == 0 ) {
            return polymesh3_cvt_vertex_accessor<T>( ch.size(), ch.type(), ch.arity(),
                                                     ch.get_writable_data().ptr_at( 0 ), channel );
        } else {
            return polymesh3_cvt_vertex_accessor<T>(
                ch.size(), ch.type(), ch.arity(), ch.get_writable_data().ptr_at( 0 ), channel,
                m_pFaceEndOffsets->size(), &ch.get_writable_faces()[0], &( *m_pFaceEndOffsets )[0] );
        }
    }

    /**
     * This function will create and return an accessor object capable of reading the given vertex channel
     * and interpreting its data, automatically converting it to the requested type.
     * @param  channel The name of the vertex channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    template <typename T>
    polymesh3_const_cvt_vertex_accessor<T> get_const_cvt_vertex_accessor( const frantic::tstring& channel ) const {
        const polymesh3_channel& ch = get_const_vertex_channel( channel );

        if( ch.get_faces().size() == 0 ) {
            return polymesh3_const_cvt_vertex_accessor<T>( ch.size(), ch.type(), ch.arity(), ch.get_data().ptr_at( 0 ),
                                                           channel );
        } else {
            return polymesh3_const_cvt_vertex_accessor<T>( ch.size(), ch.type(), ch.arity(), ch.get_data().ptr_at( 0 ),
                                                           channel, m_pFaceEndOffsets->size(), &ch.get_faces()[0],
                                                           &( *m_pFaceEndOffsets )[0] );
        }
    }

    /**
     * This function will create and return an accessor object capable of reading the given face channel
     * and interpreting its data.
     * @tparam T The type of data expected to be in the channel
     * @param  channel The name of the face channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    template <typename T>
    polymesh3_face_accessor<T> get_face_accessor( const frantic::tstring& channel ) {
        polymesh3_channel& ch = get_face_channel<T>( channel );

        return polymesh3_face_accessor<T>( ch.size(), ch.get_writable_data().ptr_at( 0 ) );
    }

    /**
     * This function will create and return an accessor object capable of reading the given face channel
     * and interpreting its data.
     * @tparam T The type of data expected to be in the channel
     * @param  channel The name of the face channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    template <typename T>
    polymesh3_const_face_accessor<T> get_const_face_accessor( const frantic::tstring& channel ) const {
        const polymesh3_channel& ch = get_const_face_channel<T>( channel );

        return polymesh3_const_face_accessor<T>( ch.size(), ch.get_data().ptr_at( 0 ) );
    }

    /**
     * This function will create and return an accessor object capable of reading the given face channel
     * and interpreting its data, but without knowing the expected type at compile time.
     * @param  channel The name of the face channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    polymesh3_face_accessor<void> get_face_accessor( const frantic::tstring& channel ) {
        polymesh3_channel& ch = get_face_channel( channel );

        return polymesh3_face_accessor<void>( ch.size(), ch.type(), ch.arity(), ch.get_writable_data().ptr_at( 0 ) );
    }

    /**
     * This function will create and return an accessor object capable of reading the given face channel
     * and interpreting its data, but without knowing the expected type at compile time.
     * @param  channel The name of the face channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    polymesh3_const_face_accessor<void> get_const_face_accessor( const frantic::tstring& channel ) const {
        const polymesh3_channel& ch = get_const_face_channel( channel );

        return polymesh3_const_face_accessor<void>( ch.size(), ch.type(), ch.arity(), ch.get_data().ptr_at( 0 ) );
    }

    /**
     * This function will create and return an accessor object capable of reading the given face channel
     * and interpreting its data, automatically converting it to the requested type.
     * @tparam T The type to convert the channel's data to.
     * @param  channel The name of the face channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    template <typename T>
    polymesh3_cvt_face_accessor<T> get_cvt_face_accessor( const frantic::tstring& channel ) {
        polymesh3_channel& ch = get_face_channel( channel );

        return polymesh3_cvt_face_accessor<T>( ch.size(), ch.type(), ch.arity(), ch.get_writable_data().ptr_at( 0 ),
                                               channel );
    }

    /**
     * This function will create and return an accessor object capable of reading the given face channel
     * and interpreting its data, automatically converting it to the requested type.
     * @tparam T The type to convert the channel's data to.
     * @param  channel The name of the face channel to create an accessor for
     * @return An accessor capable of interpreting and retreiving the data in the given channel.
     */
    template <typename T>
    polymesh3_const_cvt_face_accessor<T> get_const_cvt_face_accessor( const frantic::tstring& channel ) const {
        const polymesh3_channel& ch = get_const_face_channel( channel );

        return polymesh3_const_cvt_face_accessor<T>( ch.size(), ch.type(), ch.arity(), ch.get_data().ptr_at( 0 ),
                                                     channel );
    }
};

/**
 * This typedef should be used in order to store an instance of a polymesh3. Any other method of storing a polymesh3
 * is highly discouraged.
 */
typedef polymesh3::ptr_type polymesh3_ptr;

typedef polymesh3::const_ptr_type const_polymesh3_ptr;

/**
 *  Apply a multiplicative scale to the specified channel.
 *
 * @param[in,out] mesh The mesh that contains the channel to scale.
 * @param channelName The channel to scale.
 * @param scale The multiplicative scale to apply.
 */
void scale_channel( frantic::geometry::polymesh3_ptr mesh, const frantic::tstring& channelName, double scale );

/**
 *  Transform a polymesh3 in place, using the provided transform matrix
 * and time-derivative transform matrix.
 *
 *  This affects the vertex positions, as well as the Velocity and Normal
 * channels if they exist.
 *
 *  A possible future extension would be to add a xformTimeSecondDerivative
 * parameter, and use it in concert with the other transforms to affect the
 * Acceleration channel.  Currently this function doesn't modify the
 * Acceleration channel at all.
 *
 * @param[in,out] mesh The mesh to transform.
 * @param xform The transform matrix to apply to the polymesh3.
 * @param xformTimeDerivative The time derivative of xform, which affects
 *		the Velocity channel.
 */
void transform( frantic::geometry::polymesh3_ptr mesh, const frantic::graphics::transform4f& xform,
                const frantic::graphics::transform4f& xformTimeDerivative = frantic::graphics::transform4f::zero() );

/**
 *  Scale the polymesh3 in place, using the provided scale factor.
 * This affects the vertex positions, the Normal channel, and the Velocity
 * channel.
 *
 * @param[in,out] mesh the mesh to scale.
 * @param scale the scaling factor applied to the mesh.
 */
void scale( frantic::geometry::polymesh3_ptr mesh, float scale );

/**
 *  Reverse the face winding of all faces in the polymesh.
 *
 * @param[in,out] mesh the mesh in which to reverse the face winding.
 */
void reverse_face_winding( frantic::geometry::polymesh3_ptr mesh );

/**
 *  Perform linear interpolation beween two polymeshes, including all
 * vertex and face channels that are defined in both meshes.  If a channel
 * is defined in one mesh, but not in the other, then the channel is not
 * included in the final mesh.
 *
 *  Both meshes must have the same topology (vertex count and faces must be
 * identical).
 *
 * @param mesh1 the first mesh operand.
 * @param mesh2 the second mesh operand.
 * @param alpha the fractional distance between the two samples at which to
 *		interpolate the final sample.
 * @return the interpolated mesh.
 */
frantic::geometry::polymesh3_ptr linear_interpolate( const frantic::geometry::const_polymesh3_ptr mesh1,
                                                     const frantic::geometry::const_polymesh3_ptr mesh2, float alpha );

/**
 *  Merge all of the faces and vertices of the meshes into a single mesh.
 *
 *  Channels which are newly created due to the combination have their
 * values defaulted to zero.
 *
 * @param meshes the meshes to combine.
 * @return the combined mesh.
 */
frantic::geometry::polymesh3_ptr combine( const std::vector<polymesh3_ptr>& meshes );

/**
 * Merge all of the faces and vertices of the meshes into a single mesh.  Attempts to combine vertices that are in the
 * same position
 *
 *  Channels which are newly created due to the combination have their
 * values defaulted to zero.
 *  This can produce degenerate faces if multiple corners of a polygon are
 * welded together, use fix_degenerate_faces() in geometry/mesh_interface_utils.hpp
 * to remedy this.
 *
 * @param meshes the meshes to weld.
 * @param logger progress and cancellation logger
 * @param errorTolerance the error tolerance when determining whether two vertices are the same
 * @return the welded mesh.
 */
frantic::geometry::polymesh3_ptr weld( const std::vector<polymesh3_ptr>& meshes, float errorTolerance = 1e-8f );
frantic::geometry::polymesh3_ptr weld( const std::vector<polymesh3_ptr>& meshes,
                                       ::frantic::logging::progress_logger& logger, float errorTolerance = 1e-8f );

/**
 *  Create a new mesh, keeping only the faces indicated in the keepFaces array.
 *
 * @param mesh the mesh to cull.
 * @param keepFaces a boolean for each face in the mesh, indicating which faces
 *    should be copied to the output mesh.
 * @return the culled mesh.
 */
frantic::geometry::polymesh3_ptr cull_faces( frantic::geometry::polymesh3_ptr mesh, const std::vector<bool>& keepFaces,
                                             frantic::logging::progress_logger& logger );
frantic::geometry::polymesh3_ptr cull_faces( frantic::geometry::polymesh3_ptr mesh,
                                             const std::vector<bool>& keepFaces );

/**
 * Create a new mesh, keeping only the faces indicated in the keepFaces array.
 *
 * @param mesh the mesh to cull.
 * @param keepFacesFunctor a boolean predicate for each face in the mesh, indicating
 *    which faces should be copied to the output mesh.
 * @return the culled mesh.
 */
template <class Functor>
frantic::geometry::polymesh3_ptr cull_faces( frantic::geometry::polymesh3_ptr mesh, Functor keepFacesFunctor,
                                             frantic::logging::progress_logger& logger ) {
    std::vector<bool> keepFacesList( mesh->face_count() );

    for( size_t i = 0; i < keepFacesList.size(); ++i ) {
        keepFacesList[i] = keepFacesFunctor( i );
    }

    return cull_faces( mesh, keepFacesList, logger );
}

template <class Functor>
frantic::geometry::polymesh3_ptr cull_faces( frantic::geometry::polymesh3_ptr mesh, Functor keepFacesFunctor ) {
    frantic::logging::null_progress_logger nullLogger;
    return cull_faces( mesh, keepFacesFunctor, nullLogger );
}

/**
 *  Creates a new mesh by removing any faces whose center it not within the given bounding box
 *
 * @param mesh source mesh to cull.
 * @param bbox bounding box region to use for culling.
 * @param bboxTransform transform of the bounding box region
 * @return the culled mesh.
 */
frantic::geometry::polymesh3_ptr cull_geometry( polymesh3_ptr mesh, const frantic::graphics::boundbox3f& bbox,
                                                const frantic::graphics::transform4f& bboxTransform );

/**
 * Check if the given mesh should be culled based on the bounding box
 *
 * @param v1 first vertex of triangle
 * @param v2 second vertex of triangle
 * @param v3 third vertex of triangle
 * @param bbox bounding box region to use for culling
 * @param bboxInverseTransform inverse transform of the bounding box region
 * @param outIsOnBoundary set to true if triangle is on the region boundary
 * @return false if mesh is inside, true if mesh is outside
 */
bool check_face_region( const frantic::graphics::vector3f& v1, const frantic::graphics::vector3f& v2,
                        const frantic::graphics::vector3f& v3, const frantic::graphics::boundbox3f& bbox,
                        const frantic::graphics::transform4f& bboxInverseTransform, bool& outIsOnBoundary );

/**
 *  Create a polymesh3 from a copy of the data in the meshInterface.
 *
 * @param meshInterface the mesh data to copy.
 * @return a new mesh containing a copy of the data in the meshInterface.
 */
frantic::geometry::polymesh3_ptr create_polymesh3( const frantic::geometry::mesh_interface* meshInterface );

/**
 *  Return an apparent copy of a mesh, but without any vertex channel custom
 * face indices.
 *
 *  Vertices that use custom indices are split, such that the output mesh
 * has the same vertex channel values in its face corners as the input mesh,
 * but does not have any custom indices.
 *
 * @note In general, the output mesh will have a different vertex count and
 *     topology than the input mesh.
 *
 * @note As a side effect, this function will remove any vertices that are not
 *     referenced by a face.
 *
 * @param mesh the mesh to operate on.
 * @param cpp the channels to propagate to the returned mesh.
 * @return a new mesh without vertex channel custom face indices.
 */
frantic::geometry::polymesh3_ptr explode_custom_faces(
    const frantic::geometry::mesh_interface* mesh,
    const frantic::channels::channel_propagation_policy& cpp = frantic::channels::channel_propagation_policy() );

/**
 *  Return true if the two meshes have consistent topology, and false
 * otherwise.
 *
 * @param mesh1 the first mesh operand.
 * @param mesh2 the second mesh operand.
 * @return true if the two meshes have consistent topology, and false
 *		otherwise.
 */
bool is_consistent_topology( const frantic::geometry::polymesh3_ptr mesh1,
                             const frantic::geometry::polymesh3_ptr mesh2 );

/**
 * Return true if the two meshes have equal vertices, faces, and channels,
 * and false otherwise.
 *
 * @param mesh1 the first mesh operand.
 * @param mesh2 the second mesh operand.
 * @return true if the two meshes have equal vertices, faces, and channels,
 *		and false otherwise.
 */
bool is_equal( const frantic::geometry::const_polymesh3_ptr& mesh1,
               const frantic::geometry::const_polymesh3_ptr& mesh2 );

/**
 *  Create a bounding box that contains all of the mesh's vertex positions.
 *
 * @param mesh the mesh to build a boundbox for.
 * @return a bounding box that contains all of the mesh's vertex positions.
 */
frantic::graphics::boundbox3f compute_boundbox( const_polymesh3_ptr mesh );

/**
 * Some polygons are removed by Maya when creating a mesh.  This function
 * returns true if the face will be added by Maya, and false if it will be
 * removed.
 *
 * This function was developed through experimentation.  I'm not sure if it's
 * accurate, but I brute-force tested it with polygons of degree [3..8] in
 * Maya 2014.
 *
 * @todo optimize this?
 *
 * @param faceRange a polymesh face to test.
 * @return true if the face can be added to a Maya mesh, and false otherwise.
 */
bool is_valid_maya_face( const frantic::geometry::polymesh3_const_face_range& faceRange );

/**
 * Ignore this. This is required for boost::intrusive_ptr<polymesh3> to work.
 */
inline void intrusive_ptr_add_ref( const polymesh3* pPolyMesh ) { pPolyMesh->acquire(); }

/**
 * Ignore this. This is required for boost::intrusive_ptr<polymesh3> to work.
 */
inline void intrusive_ptr_release( const polymesh3* pPolyMesh ) { pPolyMesh->release(); }

} // namespace geometry
} // namespace frantic
