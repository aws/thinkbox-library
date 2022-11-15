// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 *	Abstract base classes for general mesh access.  The intent is to provide a consistent interface for
 *	access to channel-wise data in an underlying mesh while abstracting the particular mesh implementation.
 */
#pragma once

#include <map>
#include <vector>

#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>

// The 3ds Max SDK defines this absoultely reprehensible macro, which interacts poorly with boost::ptr_map.
#ifdef base_type
#undef base_type
#endif

#include <boost/predef.h>

#include <boost/aligned_storage.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace geometry {

class mesh_channel;
class vertex_adjacency_interface;
class face_adjacency_interface;

namespace detail {
BOOST_STATIC_CONSTANT( std::size_t, iterator_storage_size = sizeof( void* ) * 2 );
}

struct vertex_iterator {
    boost::aligned_storage<detail::iterator_storage_size>
        m_data; // Some storage for the hidden implementation of the iterator.
};

struct face_iterator {
    boost::aligned_storage<detail::iterator_storage_size>
        m_data; // Some storage for the hidden implementation of the iterator.
};

/**
 * This abstract class defines a geometric mesh consisting of a series of vertices, and a number of faces that are
 * defined by connecting a
 * co-planar group of vertices. A face can be constituted of 3 or more vertices (ie. not limited to triangles). Data is
 * associated with
 * each vertex, each face, and each vertex within a face.
 */
class mesh_interface {
  public:
    /**
     * This internal class is used to store a number mesh_channels constituting the data associated with a
     * type of geometry. It has a similar interface to frantic::channels::channel_map and std::map.
     */
    class mesh_channel_map {
        typedef boost::ptr_map<frantic::tstring, mesh_channel> impl_type;

        impl_type m_impl;

      private:
        friend class mesh_interface;

        inline void append_channel( std::unique_ptr<mesh_channel> ch );
        inline void append_channel( mesh_channel* ch );
        inline mesh_channel* remove_channel( const frantic::tstring& channelName );

        void reset() { m_impl.clear(); }

      public:
        typedef impl_type::const_iterator const_iterator;

        const_iterator begin() const { return m_impl.begin(); }

        const_iterator end() const { return m_impl.end(); }

        bool has_channel( const frantic::tstring& channelName ) const {
            return ( m_impl.find( channelName ) != m_impl.end() );
        }

        const mesh_channel* get_channel( const frantic::tstring& channelName ) const {
            impl_type::const_iterator it = m_impl.find( channelName );
            if( it != m_impl.end() )
                return it->second;
            return NULL;
        }
    };

    typedef boost::shared_ptr<mesh_interface> ptr_type;

/// Sentinel index to indicate a 'hole' when accessing faces through the iterators
#if BOOST_COMP_CLANG || __GNUC__
    enum { HOLE_INDEX = static_cast<std::size_t>( -1 ) };
#else
    static const std::size_t HOLE_INDEX = static_cast<std::size_t>( -1 );
#endif

  private:
    mesh_channel_map m_vertexChannels, m_faceChannels;

  protected:
    /**
     * Subclasses register their mesh_channel instances here so that clients may access, query and modify the meshes
     * channels
     * @param ch The vertex channel to expose.
     */
    template <class T>
    void append_vertex_channel( std::unique_ptr<T> ch ) {
        m_vertexChannels.append_channel( std::unique_ptr<mesh_channel>( ch.release() ) );
    }

    /**
     * Subclasses register their mesh_channel instances here so that clients may access, query and modify the meshes
     * channels
     * @param ch The vertex channel to expose.
     */
    template <class T>
    void append_face_channel( std::unique_ptr<T> ch ) {
        m_faceChannels.append_channel( std::unique_ptr<mesh_channel>( ch.release() ) );
    }

    /**
     * Delegate method for the subclass to create a vertex channel of the requested name and type.
     */
    virtual mesh_channel* create_vertex_channel( const frantic::tstring& /*channelName*/,
                                                 frantic::channels::data_type_t /*dataType*/, size_t /*arity*/ ) {
        throw std::runtime_error( "mesh_interface::create_vertex_channel -- This mesh_interface implementation does "
                                  "not support this function." );
    }

    /**
     * Delegate method for the subclass to create a vertex channel of the requested name and type, with custom faces.
     */
    virtual mesh_channel* create_vertex_channel_custom_faces( const frantic::tstring& /*channelName*/,
                                                              frantic::channels::data_type_t /*dataType*/,
                                                              size_t /*arity*/, size_t /*vertexCount*/ ) {
        throw std::runtime_error( "mesh_interface::create_vertex_channel_custom_faces -- This mesh_interface "
                                  "implementation does not support this function." );
    }

    /**
     * Delegate to destroy a previously created channel.
     */
    virtual void destroy_vertex_channel( mesh_channel* /*channel*/ ) {
        throw std::runtime_error(
            "mesh_interface::destroy_vertex_channel -- This mesh_interface implementation does not "
            "support this function." );
    }

    /**
     * Delegate method for the subclass to create a face channel of the requested name and type
     */
    virtual mesh_channel* create_face_channel( const frantic::tstring& /*channelName*/,
                                               frantic::channels::data_type_t /*dataType*/, size_t /*arity*/ ) {
        throw std::runtime_error( "mesh_interface::create_face_channel -- This mesh_interface implementation does not "
                                  "support this function." );
    }

    /**
     * Delegate to destroy a previously created channel.
     */
    virtual void destroy_face_channel( mesh_channel* /*channel*/ ) {
        throw std::runtime_error( "mesh_interface::destroy_face_channel -- This mesh_interface implementation does not "
                                  "support this function." );
    }

    /**
     * Clears all channels.
     */
    void reset() {
        m_vertexChannels.reset();
        m_faceChannels.reset();
    }

  public:
    virtual ~mesh_interface() {}

    mesh_channel_map& get_vertex_channels() { return m_vertexChannels; }
    mesh_channel_map& get_face_channels() { return m_faceChannels; }

    const mesh_channel_map& get_vertex_channels() const { return m_vertexChannels; }
    const mesh_channel_map& get_face_channels() const { return m_faceChannels; }

    /**
     * @return True if the mesh_interface is bound to a legitimate mesh object underneath.
     */
    virtual bool is_valid() const = 0;

    /**
     * When a channel is not accessible via get_vertex_channels().get_channel() or get_face_channels().get_channel()
     * this method can be used to have the object try and populate that channel if it is possible.
     * @param channelName The name of the channel we want to populate
     * @param vertexChannel If true, we are expecting a per-vertex channel. If false we are expecting a per-face
     * channel.
     * @param forOutput If true, we are expecting to write to this channel. If false we will only read from it. Some
     * channels cannot be written to.
     * @param throwOnError If true and the channel is unknown (or readonly and 'forOutput' was true) then an exception
     * is throw. Otherwise it returns false.
     * @return True if the channel requested was able to populated, false otherwise.
     */
    virtual bool request_channel( const frantic::tstring& channelName, bool vertexChannel, bool forOutput,
                                  bool throwOnError = true ) = 0;

    /**
     * Insert a vertex channel with the given name. Type and size are determined from the
     * templated data type.
     *
     * @param name The name to give this channel
     * @param vertexCount the number of custom vertex values to allocate
     */
    template <class DataType>
    mesh_channel* add_vertex_channel( const frantic::tstring& name ) {
        return add_vertex_channel( name, frantic::channels::channel_data_type_traits<DataType>::data_type(),
                                   frantic::channels::channel_data_type_traits<DataType>::arity() );
    }

    /**
     * Insert a vertex channel with the given name. Type and size are determined from the
     * templated data type.
     *
     * @param name The name to give this channel
     * @param vertexCount the number of custom vertex values to allocate
     */
    template <class DataType>
    mesh_channel* add_vertex_channel_custom_faces( const frantic::tstring& name, std::size_t vertexCount ) {
        return add_vertex_channel_custom_faces(
            name, frantic::channels::channel_data_type_traits<DataType>::data_type(),
            frantic::channels::channel_data_type_traits<DataType>::arity(), vertexCount );
    }

    /**
     * Inserts a vertex channel with the given name and data type.
     *
     * @param channelName Name of the channel to be inserted
     * @param dataType The type of data to store
     * @param arity The multiplicity of the requested data type to store
     */
    virtual mesh_channel* add_vertex_channel( const frantic::tstring& channelName,
                                              frantic::channels::data_type_t dataType, size_t arity ) {
        mesh_channel* channel = create_vertex_channel( channelName, dataType, arity );
        m_vertexChannels.append_channel( channel );
        return channel;
    }

    /**
     * Insert a vertex channel with the given name, type and size, with a custom face index mapping.
     *
     * @param channelName Name of the channel to be inserted
     * @param dataType The type of data to store
     * @param arity The multiplicity of the requested data type to store
     * @param vertexCount the number of custom vertex values to allocate
     */
    virtual mesh_channel* add_vertex_channel_custom_faces( const frantic::tstring& channelName,
                                                           frantic::channels::data_type_t dataType, size_t arity,
                                                           std::size_t vertexCount ) {
        mesh_channel* channel = create_vertex_channel_custom_faces( channelName, dataType, arity, vertexCount );
        m_vertexChannels.append_channel( channel );
        return channel;
    }

    /**
     * Insert a face channel with the given name. Type and size are determined from the
     * templated data type.
     *
     * @param channelName Name of the channel to be inserted
     */
    template <class DataType>
    mesh_channel* add_face_channel( const frantic::tstring& channelName ) {
        return add_face_channel( channelName, frantic::channels::channel_data_type_traits<DataType>::data_type(),
                                 frantic::channels::channel_data_type_traits<DataType>::arity() );
    }

    /**
     * Insert a face channel with the given name, type and size
     *
     * @param channelName Name of the channel to be inserted
     * @param dataType The type of data to store
     * @param arity The multiplicity of the requested data type to store
     */
    virtual mesh_channel* add_face_channel( const frantic::tstring& channelName,
                                            frantic::channels::data_type_t dataType, size_t arity ) {
        mesh_channel* channel = create_face_channel( channelName, dataType, arity );
        m_faceChannels.append_channel( channel );
        return channel;
    }

    /**
     * Remove a vertex channel with the given name
     *
     * @param channelName Name of the vertex channel to be removed
     */
    virtual void erase_vertex_channel( const frantic::tstring& channelName ) {
        mesh_channel* underlyingChannel = m_vertexChannels.remove_channel( channelName );
        if( underlyingChannel != NULL ) {
            destroy_vertex_channel( underlyingChannel );
        } else {
            throw std::runtime_error( "mesh_interface::remove_vertex_channel -- No such channel : " +
                                      frantic::strings::to_string( channelName ) );
        }
    }

    /**
     * Remove a face channel with the given name
     *
     * @param channelName Name of the face channel to be removed
     */
    virtual void erase_face_channel( const frantic::tstring& channelName ) {
        mesh_channel* underlyingChannel = m_faceChannels.remove_channel( channelName );
        if( underlyingChannel != NULL ) {
            destroy_face_channel( underlyingChannel );
        } else {
            throw std::runtime_error( "mesh_interface::remove_face_channel -- No such channel : " +
                                      frantic::strings::to_string( channelName ) );
        }
    }

    /**
     * Check if the given vertex channel exists
     *
     * @param channelName Name of the vertex channel to query
     */
    bool has_vertex_channel( const frantic::tstring& channelName ) const {
        return m_vertexChannels.has_channel( channelName );
    }

    /**
     * Return the given vertex channel, or NULL if it doesn't exist
     *
     * @param channelName Name of the vertex channel to query
     */
    const mesh_channel* get_vertex_channel( const frantic::tstring& channelName ) const {
        return m_vertexChannels.get_channel( channelName );
    }

    /**
     * Check if the given face channel exists
     *
     * @param channelName Name of the vertex channel to query
     */
    bool has_face_channel( const frantic::tstring& channelName ) const {
        return m_faceChannels.has_channel( channelName );
    }

    /**
     * Return the given face channel, or NULL if it doesn't exist
     *
     * @param channelName Name of the face channel to query
     */
    const mesh_channel* get_face_channel( const frantic::tstring& channelName ) const {
        return m_faceChannels.get_channel( channelName );
    }

    /**
     * @return The number of vertices in the mesh.
     */
    virtual std::size_t get_num_verts() const = 0;

    /**
     * Gets a single vertex
     * @param index Which vertex
     * @param outValues Result is stored here BTW, that is the syntax for a reference to a 3 float array.
     */
    virtual void get_vert( std::size_t index, float ( &outValues )[3] ) const = 0;

    /**
     * Convenience to get a vertex as a vector3f
     * @param index Vertex to get
     * @return the location of the vertex
     */
    frantic::graphics::vector3f get_vert( std::size_t index ) const {
        float v[3];
        get_vert( index, v );
        return frantic::graphics::vector3f( v[0], v[1], v[2] );
    }

    virtual bool is_read_only() const { return true; }

    /**
     * Set the geometry of this vertex.
     * This is meant to be an optional method to implement if the mesh is modifyable
     *
     * @param index the vertex to modify
     * @param v a 3-tuple of floats to set the vertex position
     */
    virtual void set_vert( std::size_t /*index*/, const float* /*v*/ ) {
        throw std::runtime_error( "mesh_interface::set_vert: Error, mesh is read-only" );
    }

    void set_vert( std::size_t index, const frantic::graphics::vector3f& v ) { set_vert( index, &v[0] ); }

    /**
     * @return The number of faces in the mesh.
     */
    virtual std::size_t get_num_faces() const = 0;

    /**
     * @return The number of vertices associated with the specified face. ie. 3 for triangles, 4 for quads, etc.
     */
    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const = 0;

    /**
     * @param faceIndex Which face
     * @param fvertIndex Which vertex in the face. In [ 0, get_num_face_verts(faceIndex) )
     * @return The vertex index that is at the specified slot in the face.
     */
    virtual std::size_t get_face_vert_index( std::size_t faceIndex, std::size_t fvertIndex ) const = 0;

    frantic::graphics::vector3f get_face_vert( std::size_t faceIndex, std::size_t fvertIndex ) const {
        return get_vert( get_face_vert_index( faceIndex, fvertIndex ) );
    }

    /**
     * Gets the vertex indices for all verts in the face.
     * @param faceIndex Which face
     * @param outValues Vertex indices for this face are stored here. Must be at least std::size_t[
     * get_num_face_verts(faceIndex) ].
     */
    virtual void get_face_vert_indices( std::size_t faceIndex, std::size_t outValues[] ) const = 0;

    /**
     * Gets the vertex positions for all verts in the face.
     * @param faceIndex Which face
     * @param outValues Result is stored here. This array must be at least float[ get_num_face_verts(faceIndex) ][3].
     */
    virtual void get_face_verts( std::size_t faceIndex, float outValues[][3] ) const = 0;

    /**
     * A mesh "element" is a connected region of polygons that may or may not be closed. If there are multiple elements
     * then that indicates
     * multiple unconnected regions of connected polygons.
     * @note We don't provide an obvious iteration method for elements because a mesh is stored as a collection of faces
     * which is a lower-level
     *       abstraction than elements. The simplest way to work with elements is to query the FaceElement channel which
     * gives an element index
     *       for each face.
     * @return The number of disjoint collections of connected plygons in the mesh.
     */
    virtual std::size_t get_num_elements() const = 0;

    /**
     *  A mesh "element" is a connected region of polygons that may or may not be closed. This method queries the
     * element index that a face is a part of.
     * @note In order to use this channel, you must first call request_channel( "FaceElementIndex", false, false, false
     * )
     * @param faceIndex The face to query
     * @return The index of the element this face is a part of.
     */
    virtual std::size_t get_face_element_index( std::size_t faceIndex ) const = 0;

    /**
     * Initializes any data structures required for iterating over the neighbors of vertices and faces. This method
     * should be called once before any any iterators are used. By default no adjacency information is stored until it
     * is requested.
     */
    virtual void init_adjacency() = 0;

    virtual bool has_adjacency() const = 0;

    /**
     * Initializes an iterator that visits all the edges connected to a specific vertex.
     * \note Iterators are implemented as obscured objects so that no heap allocation is required at runtime.
     * \return True if the vertex is associated with any edge, false if the vertex is disconnected and has no
     * edges/faces using it.
     */
    virtual bool init_vertex_iterator( frantic::geometry::vertex_iterator& vIt, std::size_t vertexIndex ) const = 0;

    /**
     * Moves an iterator to the next edge in the 1-ring of the initial vertex.
     * \return True if the iterator was able to be advanced, or false if the iteration is complete.
     */
    virtual bool advance_vertex_iterator( frantic::geometry::vertex_iterator& vIt ) const = 0;

    /**
     * Retrieves the vertex on the other end of the currently visited edge.
     * \return The vertex id.
     */
    virtual std::size_t get_edge_endpoint( frantic::geometry::vertex_iterator& vIt ) const = 0;

    /**
     * Retrieves the face to the left of the currently visited edge.
     * \return The face id.
     */
    virtual std::size_t get_edge_left_face( frantic::geometry::vertex_iterator& vIt ) const = 0;

    /**
     * Retrieves the face to the right of the currently visited edge.
     * \return The face id.
     */
    virtual std::size_t get_edge_right_face( frantic::geometry::vertex_iterator& vIt ) const = 0;

    /**
     * Determines if the current edge being iterated on is marked as "visible".
     * \return True if the current edge is visible, false if it is invisible.
     */
    virtual bool is_edge_visible( frantic::geometry::vertex_iterator& vIt ) const = 0;

    /**
     * Determines if the current edge being iterated on is a boundary edge.
     * \note This is equivalent to: ( get_edge_left_face() == -1 || get_edge_right_face() == -1 )
     * \return True if the current edge is a boundary edge.
     */
    virtual bool is_edge_boundary( frantic::geometry::vertex_iterator& vIt ) const = 0;

    /**
     * Initializes an iterator that visits all the faces neighboring a specific face.
     * \note Iterators are implemented as obscured objects so that no heap allocation is required at runtime.
     */
    virtual void init_face_iterator( frantic::geometry::face_iterator& fIt, std::size_t faceIndex ) const = 0;

    /**
     * Moves an iterator to the next edge in the face so that its neighbor can be queried. Once this method returns
     * false, you must not call get_face_neighbor() with the face_iterator object. \return True if the iterator was able
     * to be advanced, or false if the iteration is complete.
     */
    virtual bool advance_face_iterator( frantic::geometry::face_iterator& fIt ) const = 0;

    /**
     * Retrieves the id of the other face connected to the current edge of the original face. This can be -1 for
     * boundary edges (ie. there is no other face using this edge).
     * \return The id of the neighbor face, or -1 if there is no neighbor.
     */
    virtual std::size_t get_face_neighbor( frantic::geometry::face_iterator& fIt ) const = 0;

    /**
     * Retrieves the id of the vertex at the start of the current edge.
     * @return The id of the vertex
     */
    virtual std::size_t get_face_prev_vertex( frantic::geometry::face_iterator& fIt ) const = 0;

    /**
     * Retrieves the id of the vertex at the end of the current edge.
     * @return The id of the vertex
     */
    virtual std::size_t get_face_next_vertex( frantic::geometry::face_iterator& fIt ) const = 0;
};

/**
 * This abstract interface describes the manner in which one would access channel information from a mesh_interface.
 * There are three principle types of channels we work with: per-vertex, per-face, and per-vertex with a custom face
 * mapping. Data is stored with one value for each vertex (ex. position, selection, etc.), or one value per-face (ex.
 * material index, area, etc.), or a distinct value for each vertex in a custom face (ex. normals, texture coordinates,
 * etc.). In 3dsMax this is how map channels are implemented.
 */
class mesh_channel_interface {
  public:
    virtual ~mesh_channel_interface() {}

    /**
     * Copies the value of the specified vertex (or face depending on the channel) into 'outValue'.
     * @param index The index into this channel's data.
     * @param outValue Pointer to a variable to assign this value to. It is expected to be of the correct type, and some
     *                 other mechanism is required to determine what this type is.
     */
    virtual void get_value( std::size_t index, void* outValue ) const = 0;

    /**
     * Assigns the indexed face's or vertex's value.
     * @param index The index of the face or vertex to assign to.
     * @param value A pointer to new value. It is expected to be of the correct type, and some other mechanism
     *              is required to determine what this type is.
     */
    virtual void set_value( std::size_t index, const void* value ) const = 0;

    /**
     * Some meshes use separate indexes per-vertex of each face into a list of values for that channel. It simplifies
     * some situations where some vertices have the same value for all faces it is a part of, but other vertices have
     * discontinuities across face boundaries (ex. texture coords across the edges of a cube vs. across the vertices in
     * the center of a side). This method provides for this extra layer of indirection, translating a face&face-vertex
     * index pair into an index that can be passed to get_value().
     *
     * @param faceIndex The index of the face to inspect
     * @param fvertIndex The index of the vertex within the face to inspect. This ranges from 0 to the number of
     * vertices in a face (ie. 3 for a triangle mesh).
     */
    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const = 0;

    /**
     * Allows assignment of face vertex values. This should throw an exception if the channel in question does not
     * support custom faces.
     *
     * @param faceIndex The index of the face
     * @param fvertIndex The index of the vertex on the face.
     * @param value The vertex index to assign to this face vertex
     */
    virtual void set_fv_index( std::size_t /*faceIndex*/, std::size_t /*fvertIndex*/, std::size_t /*value*/ ) const {
        throw std::runtime_error(
            "mesh_interface::set_fv_index -- This mesh_interface implementation does not support this function." );
    }

    /**
     * Gets the number of vertices contained by the specified face.
     * @param faceIndex Which face to query.
     * @return The number of vertices included in the specified face. ie. 3 for triangles, 4 for quads, etc.
     */
    virtual std::size_t get_num_face_verts( std::size_t faceIndex ) const = 0;

    /**
     * This helper function is useful when you know (at compile time) which type the channel stores.
     */
    template <class T>
    inline T get_value( std::size_t index ) const {
        T result;
        this->get_value( index, &result );
        return result;
    }

    /**
     * This helper function will do a linear combination of this channel's values at the vertices of a face.
     * @param faceIndex which face
     * @param numVerts The number of vertices in the requested face. This is a bit weird, and should probably be
     * accessible via 'this'.
     * @param weights An array ofof 'numVerts' weights.
     */
    template <class T>
    inline T get_barycentric_value( std::size_t faceIndex, std::size_t numVerts, float weights[] ) const {
        T result;
        for( std::size_t i = 0; i < numVerts; ++i )
            result += weights[i] * this->get_value<T>( this->get_fv_index( faceIndex, i ) );
        return result;
    }
};

/**
 * This abstract class decorates the mesh_channel_interface with meta-data information about the channel.
 */
class mesh_channel : public mesh_channel_interface {
  public:
    enum channel_type {
        // Contains one data value per mesh vertex, exactly like the vertex positions.
        vertex,
        // Contains one data value per face.
        face,
        // Like 'vertex', but with a custom array of face indices that are called "custom faces."
        face_vertex,
        // If the mesh contains a "FaceElement" channel, which indicates which element each face
        //    belongs to, this channel contains one data value per mesh element as indexed by that channel.
        // See FranticMaxLibrary for example usage.
        element
    };

    inline static bool is_stored_at_vertex( channel_type channelType ) {
        return channelType == vertex || channelType == face_vertex;
    }

    struct transform_type {
        enum enum_t {
            none,
            point,
            vector,
            normal,
            scale // Different ways a channel's value is affected by transforms. Only applies to float16/32/64[3]
                  // channels, except 'scale' which applies to float16/32/64[1] channels.
        };
    };

  private:
    frantic::tstring m_name;
    channel_type m_objectType;
    frantic::channels::data_type_t m_dataType;
    std::size_t m_dataArity, m_numData, m_numFaces, m_elementSize;
    bool m_writeable;
    transform_type::enum_t m_transformType;

  protected:
    mesh_channel( const frantic::tstring& name, channel_type objectType, frantic::channels::data_type_t dataType,
                  std::size_t dataArity, std::size_t numData, std::size_t numFaces, bool readOnly = false )
        : m_name( name )
        , m_transformType( transform_type::none ) {
        m_objectType = objectType;
        m_dataType = dataType;
        m_dataArity = dataArity;
        m_numData = numData;
        m_numFaces = numFaces;
        m_writeable = !readOnly;
        m_elementSize = dataArity * frantic::channels::sizeof_channel_data_type( dataType );
    }

    void set_transform_type( transform_type::enum_t transformType ) {
        // We allow anything to be set to none, but only float16/32/64[3] for everything except float16/32/64[1] for
        // scale.
        if( transformType == transform_type::none ||
            ( frantic::channels::is_channel_data_type_float( m_dataType ) &&
              ( ( m_dataArity == 1 && transformType == transform_type::scale ) ||
                ( m_dataArity == 3 && transformType != transform_type::scale ) ) ) )
            m_transformType = transformType;
    }

    void set_num_elements( std::size_t numElements ) { m_numData = numElements; }

    void set_num_faces( std::size_t numFaces ) { m_numFaces = numFaces; }

  public:
    virtual ~mesh_channel() {}

    /**
     * @return The name of this channel.
     */
    inline const frantic::tstring& get_name() const { return m_name; }

    /**
     * @return An enum representing the geometry this channel is associated with.
     */
    inline channel_type get_channel_type() const { return m_objectType; }

    /**
     * @return The number of values stored in this channel.
     */
    inline std::size_t get_num_elements() const { return m_numData; }

    /**
     * @return The number of faces in the mesh associated with this channel
     */
    inline std::size_t get_num_faces() const { return m_numFaces; }

    /**
     * @return The data type that contains the value at each element of this channel. ex. A float[3] per vertex would
     * return data_type_float32 for get_data_type().
     */
    inline frantic::channels::data_type_t get_data_type() const { return m_dataType; }

    /**
     * @return The number of values that are stored per element of this channel. ex. A float[3] per vertex would return
     * 3 for get_data_arity().
     */
    inline std::size_t get_data_arity() const { return m_dataArity; }

    /**
     * @return The size of the data associated with one vertex or face. sizeof(get_data_type()) * get_data_arity().
     */
    inline std::size_t get_element_size() const { return m_elementSize; }

    /**
     * @return True if the channel is able to be written to. If false, any calls to this->set_value() are discarded.
     */
    inline bool is_writeable() const { return m_writeable; }

    /**
     * @return The type of transformation to apply to these channels values when transforming between spaces.
     */
    inline transform_type::enum_t get_transform_type() const { return m_transformType; }

  public:
    // See mesh_channel_interface::get_value()
    virtual void get_value( std::size_t index, void* outValue ) const = 0;

    // See mesh_channel_interface::set_value()
    virtual void set_value( std::size_t index, const void* value ) const = 0;

    // See mesh_channel_interface::get_fv_index()
    virtual std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const = 0;
};

/*
 * A converting accessor for mesh_channel
 */
template <class T>
class mesh_channel_cvt {

  protected:
    const frantic::geometry::mesh_channel* m_channel;
    frantic::channels::channel_type_convertor_function_t m_convertGet;
    frantic::channels::channel_type_convertor_function_t m_convertSet;

    void init() {
        if( m_channel == NULL )
            throw std::runtime_error( "mesh_channel_cvt() - was initialized with a NULL channel" );

        frantic::channels::data_type_t targetDataType = frantic::channels::channel_data_type_traits<T>::data_type();
        const size_t targetDataArity = frantic::channels::channel_data_type_traits<T>::arity();

        // TODO: we might consider revisiting mapping to a higher arity type when we make separate
        //   const/non-const accessors. For now, it is an error.
        if( targetDataArity > m_channel->get_data_arity() ) {
            throw std::runtime_error(
                "mesh_channel_cvt() - Cannot convert channel: \"" +
                frantic::strings::to_string( m_channel->get_name() ) + "\" from type: " +
                frantic::strings::to_string(
                    channel_data_type_str( m_channel->get_data_arity(), m_channel->get_data_type() ) ) +
                " to " + frantic::strings::to_string( channel_data_type_str( targetDataArity, targetDataType ) ) +
                " due to increasing arity" );
        }

        m_convertGet = frantic::channels::get_channel_type_convertor_function( m_channel->get_data_type(),
                                                                               targetDataType, m_channel->get_name() );
        m_convertSet = frantic::channels::get_channel_type_convertor_function(
            targetDataType, m_channel->get_data_type(), m_channel->get_name() );
    }

  public:
    mesh_channel_cvt()
        : m_channel( 0 ) {}

    mesh_channel_cvt( const frantic::geometry::mesh_channel* channel )
        : m_channel( channel ) {
        init();
    }

    mesh_channel_cvt( const mesh_channel_cvt& other )
        : m_channel( other.m_channel ) {
        init();
    }

    mesh_channel_cvt( const frantic::geometry::mesh_interface::mesh_channel_map& meshChannelMap,
                      const frantic::tstring& channelName ) {
        m_channel = meshChannelMap.get_channel( channelName );
        init();
    }

    bool is_valid() const { return m_channel != 0; }
    const frantic::tstring& get_name() const { return m_channel->get_name(); }
    frantic::geometry::mesh_channel::channel_type get_channel_type() const { return m_channel->get_channel_type(); }
    std::size_t get_num_elements() const { return m_channel->get_num_elements(); }
    std::size_t get_num_faces() const { return m_channel->get_num_faces(); }
    frantic::channels::data_type_t get_data_type() const { return m_channel->get_data_type(); }
    std::size_t get_data_arity() const { return m_channel->get_data_arity(); }
    std::size_t get_element_size() const { return m_channel->get_element_size(); }
    bool is_writeable() const { return m_channel->is_writeable(); }
    frantic::geometry::mesh_channel::transform_type::enum_t get_transform_type() const {
        return m_channel->get_transform_type();
    }

    // @return T The converted value at index
    const T get_value( std::size_t index ) const {
        T result;
        char* temp = (char*)alloca( m_channel->get_element_size() );
        std::memset( temp, 0, m_channel->get_element_size() );
        m_channel->get_value( index, temp );
        m_convertGet( (char*)&result, temp,
                      std::min( frantic::channels::channel_data_type_traits<T>::arity(), get_data_arity() ) );
        return result;
    }

    const T get_fv_value( size_t face, size_t fvertIndex ) const {
        return get_value( get_fv_index( face, fvertIndex ) );
    }

    void set_value( std::size_t index, const T& value ) {
        char* temp = (char*)alloca( m_channel->get_element_size() );
        std::memset( temp, 0, m_channel->get_element_size() );
        m_convertSet( temp, (const char*)&value,
                      std::min( frantic::channels::channel_data_type_traits<T>::arity(), get_data_arity() ) );
        m_channel->set_value( index, temp );
    }

    void set_fv_value( size_t face, size_t fvertIndex, const T& value ) {
        set_value( get_fv_index( face, fvertIndex ), value );
    }

    std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
        return m_channel->get_fv_index( faceIndex, fvertIndex );
    }

    void set_fv_index( size_t faceIndex, size_t fvertIndex, size_t value ) const {
        m_channel->set_fv_index( faceIndex, fvertIndex, value );
    }
};

void mesh_interface::mesh_channel_map::append_channel( std::unique_ptr<mesh_channel> ch ) {
    // Cannot do this, since copy of ch (and hence invalidation) was happening before get_name() ... Then KABOOM!
    // m_impl.insert( ch->get_name(), ch );

    const frantic::tstring& name = ch->get_name();

    m_impl.insert( name, std::move( ch ) );
}

void mesh_interface::mesh_channel_map::append_channel( mesh_channel* ch ) {
    // This is an exceptionally confusing part of the ptr_map interface where the non-unique_ptr
    // version takes a non-const reference to the key

    frantic::tstring name = ch->get_name();
    m_impl.insert( name, ch );
}

mesh_channel* mesh_interface::mesh_channel_map::remove_channel( const frantic::tstring& channelName ) {
    impl_type::iterator it = m_impl.find( channelName );

    if( it != m_impl.end() ) {
        return m_impl.release( it ).release();
    } else {
        throw std::runtime_error( "mesh_interface::mesh_channel_map::remove_channel -- No such channel: " +
                                  frantic::strings::to_string( channelName ) );
        // return NULL;
    }
}

typedef mesh_interface::ptr_type mesh_interface_ptr;
} // namespace geometry
} // namespace frantic
