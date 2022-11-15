// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/mesh_interface_utils.hpp>

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/geometry/connected_components.hpp>
#include <frantic/geometry/delegated_const_mesh_interface.hpp>
#include <frantic/geometry/polygon_utils.hpp>
#include <frantic/geometry/polymesh3.hpp>
#include <frantic/geometry/polymesh3_builder.hpp>
#include <frantic/geometry/polymesh3_interface.hpp>
#include <frantic/misc/utility.hpp>
#include <frantic/particles/particle_array.hpp>

#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/pending/disjoint_sets.hpp>

namespace frantic {
namespace geometry {
bool is_equal( const frantic::geometry::mesh_interface* mesh, const frantic::geometry::mesh_interface* otherMesh ) {
    if( !mesh || !otherMesh )
        throw std::runtime_error( "is_equal: was passed a null mesh" );

    // check the geom verts
    if( mesh->get_num_verts() != otherMesh->get_num_verts() )
        return false;
    for( std::size_t i = 0, ie = mesh->get_num_verts(); i != ie; ++i ) {
        // get the verts at index i
        float curVert[3];
        mesh->get_vert( i, curVert );
        float otherCurVert[3];
        otherMesh->get_vert( i, otherCurVert );
        // check the vert values
        if( curVert[0] != otherCurVert[0] || curVert[1] != otherCurVert[1] || curVert[2] != otherCurVert[2] )
            return false;
    }

    // check the geom faces
    if( mesh->get_num_faces() != otherMesh->get_num_faces() )
        return false;
    for( std::size_t curFace = 0, numFaces = mesh->get_num_faces(); curFace < numFaces; ++curFace ) {
        if( mesh->get_num_face_verts( curFace ) != otherMesh->get_num_face_verts( curFace ) )
            return false;
        for( std::size_t curFaceVert = 0, numFaceVerts = mesh->get_num_face_verts( curFace );
             curFaceVert < numFaceVerts; ++curFaceVert ) {
            if( mesh->get_face_vert_index( curFace, curFaceVert ) !=
                otherMesh->get_face_vert_index( curFace, curFaceVert ) )
                return false;
        }
    }

    // check the vert channels
    const frantic::geometry::mesh_interface::mesh_channel_map& vertexChannels = mesh->get_vertex_channels();
    frantic::geometry::mesh_interface::mesh_channel_map::const_iterator it = vertexChannels.begin();
    frantic::geometry::mesh_interface::mesh_channel_map::const_iterator itEnd = vertexChannels.end();
    // make sure they have the same channels and number of channels
    std::size_t numVertexChannels = 0;
    for( ; it != itEnd; ++it ) {
        ++numVertexChannels;
        if( !otherMesh->has_vertex_channel( it->second->get_name() ) )
            return false;
    }
    const frantic::geometry::mesh_interface::mesh_channel_map& otherVertexChannels = otherMesh->get_vertex_channels();
    it = otherVertexChannels.begin();
    itEnd = otherVertexChannels.end();

    std::size_t numOtherVertexChannels = std::distance( it, itEnd );
    if( numVertexChannels != numOtherVertexChannels )
        return false;
    // make sure they have the same channel data
    it = vertexChannels.begin();
    itEnd = vertexChannels.end();
    for( ; it != itEnd; ++it ) {
        const frantic::geometry::mesh_channel* channel = it->second;
        const frantic::geometry::mesh_channel* otherChannel = otherVertexChannels.get_channel( channel->get_name() );
        // check the channels elements
        if( channel->get_num_elements() != otherChannel->get_num_elements() )
            return false;
        if( channel->get_data_arity() != otherChannel->get_data_arity() )
            return false;
        if( channel->get_data_type() != otherChannel->get_data_type() )
            return false;
        if( channel->get_element_size() != otherChannel->get_element_size() )
            return false;
        std::vector<char> value( channel->get_element_size() );
        std::vector<char> otherValue( otherChannel->get_element_size() );
        for( std::size_t elementIndex = 0, numElements = channel->get_num_elements(); elementIndex < numElements;
             ++elementIndex ) {
            channel->get_value( elementIndex, &value[0] );
            otherChannel->get_value( elementIndex, &otherValue[0] );
            for( std::size_t i = 0, ie = value.size(); i < ie; ++i ) {
                if( value[i] != otherValue[i] )
                    return false;
            }
        }
        // check the channels faces
        if( channel->get_num_faces() != otherChannel->get_num_faces() )
            return false;
        for( std::size_t curFace = 0, numFaces = channel->get_num_faces(); curFace < numFaces; ++curFace ) {
            if( channel->get_num_face_verts( curFace ) != otherChannel->get_num_face_verts( curFace ) )
                return false;
            for( std::size_t curFaceVert = 0, numFaceVerts = channel->get_num_face_verts( curFace );
                 curFaceVert < numFaceVerts; ++curFaceVert ) {
                if( channel->get_fv_index( curFace, curFaceVert ) !=
                    otherChannel->get_fv_index( curFace, curFaceVert ) )
                    return false;
            }
        }
    }

    // check the face channels
    const frantic::geometry::mesh_interface::mesh_channel_map& faceChannels = mesh->get_face_channels();
    it = faceChannels.begin();
    itEnd = faceChannels.end();
    // make sure they have the same channels and number of channels
    std::size_t numFaceChannels = 0;
    for( ; it != itEnd; ++it ) {
        ++numFaceChannels;
        if( !otherMesh->has_face_channel( it->second->get_name() ) )
            return false;
    }
    const frantic::geometry::mesh_interface::mesh_channel_map& otherFaceChannels = otherMesh->get_face_channels();
    it = otherFaceChannels.begin();
    itEnd = otherFaceChannels.end();
    std::size_t numOtherFaceChannels = 0;
    for( ; it != itEnd; ++it ) {
        ++numOtherFaceChannels;
    }
    if( numFaceChannels != numOtherFaceChannels )
        return false;
    // make sure they have the same channel data
    it = faceChannels.begin();
    itEnd = faceChannels.end();
    for( ; it != itEnd; ++it ) {
        const frantic::geometry::mesh_channel* channel = it->second;
        const frantic::geometry::mesh_channel* otherChannel = otherFaceChannels.get_channel( channel->get_name() );
        // check the channels elements
        if( channel->get_num_elements() != otherChannel->get_num_elements() )
            return false;
        if( channel->get_data_arity() != otherChannel->get_data_arity() )
            return false;
        if( channel->get_data_type() != otherChannel->get_data_type() )
            return false;
        if( channel->get_element_size() != otherChannel->get_element_size() )
            return false;
        std::vector<char> value( channel->get_element_size() );
        std::vector<char> otherValue( otherChannel->get_element_size() );
        for( std::size_t elementIndex = 0, numElements = channel->get_num_elements(); elementIndex < numElements;
             ++elementIndex ) {
            channel->get_value( elementIndex, &value[0] );
            otherChannel->get_value( elementIndex, &otherValue[0] );
            for( std::size_t i = 0, ie = value.size(); i < ie; ++i ) {
                if( value[i] != otherValue[i] )
                    return false;
            }
        }
    }

    // we didn't find any differences so we return true
    return true;
}

frantic::graphics::boundbox3f compute_boundbox( const frantic::geometry::mesh_interface* mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "compute_boundbox Error: mesh is NULL" );
    }

    frantic::graphics::boundbox3f result;

    for( std::size_t i = 0, ie = mesh->get_num_verts(); i != ie; ++i ) {
        const frantic::graphics::vector3f p = mesh->get_vert( i );
        result += p;
    }

    return result;
}

std::size_t get_face_vertex_count( const frantic::geometry::mesh_interface* mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "get_face_vertex_count Error: mesh is NULL" );
    }
    if( !mesh->is_valid() ) {
        throw std::runtime_error( "get_face_vertex_count Error: mesh is not valid" );
    }

    std::size_t sum = 0;

    for( std::size_t i = 0, ie = mesh->get_num_faces(); i < ie; ++i ) {
        sum += mesh->get_num_face_verts( i );
    }

    return sum;
}

bool is_closed_manifold( const frantic::geometry::mesh_interface* mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "is_closed_manifold Error: mesh is NULL" );
    }

    typedef std::pair<std::size_t, std::size_t> edge_t;

    typedef std::map<edge_t, std::size_t> map_t;
    map_t edgeIncidenceCount;

    // count the number faces incident on each edge
    for( std::size_t faceIndex = 0, faceIndexEnd = mesh->get_num_faces(); faceIndex < faceIndexEnd; ++faceIndex ) {
        const std::size_t cornerCount = mesh->get_num_face_verts( faceIndex );
        for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
            const std::size_t nextCornerIndex = ( cornerIndex + 1 ) % cornerCount;

            const std::size_t cornerVertexIndex = mesh->get_face_vert_index( faceIndex, cornerIndex );
            const std::size_t nextCornerVertexIndex = mesh->get_face_vert_index( faceIndex, nextCornerIndex );

            const edge_t edge = make_sorted_pair( cornerVertexIndex, nextCornerVertexIndex );

            map_t::iterator i = edgeIncidenceCount.lower_bound( edge );

            if( i != edgeIncidenceCount.end() && i->first == edge ) {
                ++( i->second );
            } else {
                edgeIncidenceCount.insert( i, map_t::value_type( edge, 1 ) );
            }
        }
    }

    // make sure all edges are incident on exactly two faces
    for( map_t::iterator i = edgeIncidenceCount.begin(), ie = edgeIncidenceCount.end(); i != ie; ++i ) {
        if( i->second != 2 ) {
            return false;
        }
    }

    return true;
}

bool has_degenerate_faces( const mesh_interface_ptr& mesh ) {
    if( !mesh ) {
        throw std::runtime_error( "has_degenerate_faces Error: mesh is NULL" );
    }

    std::vector<std::size_t> face;

    for( std::size_t faceIndex = 0, faceCount = mesh->get_num_faces(); faceIndex < faceCount; ++faceIndex ) {
        const std::size_t fvCount = mesh->get_num_face_verts( faceIndex );

        face.resize( fvCount );
        mesh->get_face_vert_indices( faceIndex, &face[0] );

        std::sort( face.begin(), face.end() );
        if( std::adjacent_find( face.begin(), face.end() ) != face.end() )
            return true;
    }

    return false;
}

// A set of 0 or more faces which replace a degenerate face in a mesh. The values in these faces are indices into the
// list of face-vertex indices of the original face, they are not face-vertex indicies themselves.
typedef std::vector<std::vector<std::size_t>> face_replacement;

static void reconstruct_vertex_channel( const mesh_channel* channel,
                                        const std::vector<std::pair<std::size_t, face_replacement>>& changeset,
                                        std::size_t newFaceVertCount, polymesh3_ptr& targetMesh ) {
    const std::size_t elementSize = channel->get_element_size();
    const std::size_t elementCount = channel->get_num_elements();
    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( elementSize * elementCount );

    // Channel data is unaffected, copy it directly
    for( std::size_t elementIndex = 0; elementIndex < elementCount; ++elementIndex ) {
        channel->get_value( elementIndex, buffer.ptr_at( elementSize * elementIndex ) );
    }

    if( channel->get_channel_type() == mesh_channel::face_vertex ) {
        const std::size_t originalFaceCount = channel->get_num_faces();
        std::vector<int> faceBuffer;
        faceBuffer.reserve( newFaceVertCount );

        std::size_t originalFaceIndex = 0;
        for( std::vector<std::pair<std::size_t, face_replacement>>::const_iterator change = changeset.begin(),
                                                                                   changeEnd = changeset.end();
             change < changeEnd; ++change ) {
            std::size_t nextReplacement = change->first;

            // Directly copy faces until the next modified face
            for( std::size_t faceIndex = originalFaceIndex; faceIndex < nextReplacement; ++faceIndex ) {
                for( std::size_t fvIndex = 0, fvCount = channel->get_num_face_verts( faceIndex ); fvIndex < fvCount;
                     ++fvIndex ) {
                    faceBuffer.push_back( static_cast<int>( channel->get_fv_index( faceIndex, fvIndex ) ) );
                }
            }
            // Mimic the new topology for the modified face
            for( std::size_t faceIndex = 0, faceCount = change->second.size(); faceIndex < faceCount; ++faceIndex ) {
                for( std::size_t fvIndex = 0, fvCount = change->second[faceIndex].size(); fvIndex < fvCount;
                     ++fvIndex ) {
                    faceBuffer.push_back( static_cast<int>(
                        channel->get_fv_index( change->first, change->second[faceIndex][fvIndex] ) ) );
                }
            }
            originalFaceIndex = change->first + 1;
        }
        // Directly copy faces until the next modified face
        for( std::size_t faceIndex = originalFaceIndex; faceIndex < originalFaceCount; ++faceIndex ) {
            for( std::size_t fvIndex = 0, fvCount = channel->get_num_face_verts( faceIndex ); fvIndex < fvCount;
                 ++fvIndex ) {
                faceBuffer.push_back( static_cast<int>( channel->get_fv_index( faceIndex, fvIndex ) ) );
            }
        }

        targetMesh->add_vertex_channel( channel->get_name(), channel->get_data_type(), channel->get_data_arity(),
                                        buffer, &faceBuffer );
    } else {
        targetMesh->add_vertex_channel( channel->get_name(), channel->get_data_type(), channel->get_data_arity(),
                                        buffer );
    }
}

static void reconstruct_face_channel( const mesh_channel* channel,
                                      const std::vector<std::pair<std::size_t, face_replacement>>& changeset,
                                      std::size_t newFaceCount, polymesh3_ptr& targetMesh ) {
    const std::size_t elementSize = channel->get_element_size();
    const std::size_t elementCount = channel->get_num_elements();
    frantic::graphics::raw_byte_buffer buffer;
    buffer.resize( elementSize * newFaceCount );

    std::size_t totalFaceIndex = 0;
    std::size_t originalFaceIndex = 0;
    for( std::vector<std::pair<std::size_t, face_replacement>>::const_iterator change = changeset.begin(),
                                                                               changeEnd = changeset.end();
         change < changeEnd; ++change ) {
        std::size_t nextReplacement = change->first;

        // Directly copy elements until the next modified face
        for( std::size_t faceIndex = originalFaceIndex; faceIndex < nextReplacement; ++faceIndex, ++totalFaceIndex ) {
            channel->get_value( faceIndex, buffer.ptr_at( elementSize * totalFaceIndex ) );
        }

        // Copy the element once for each face replacing the old one
        for( std::size_t targetIndex = totalFaceIndex + change->second.size(); totalFaceIndex < targetIndex;
             ++totalFaceIndex ) {
            channel->get_value( change->first, buffer.ptr_at( elementSize * totalFaceIndex ) );
        }
        originalFaceIndex = change->first + 1;
    }
    // Directly copy elements until the next modified face
    for( std::size_t faceIndex = originalFaceIndex; faceIndex < elementCount; ++faceIndex, ++totalFaceIndex ) {
        channel->get_value( faceIndex, buffer.ptr_at( elementSize * totalFaceIndex ) );
    }

    targetMesh->add_face_channel( channel->get_name(), channel->get_data_type(), channel->get_data_arity(), buffer );
}

polymesh3_ptr fix_degenerate_faces( const mesh_interface_ptr& mesh ) {
    using frantic::graphics::raw_byte_buffer;

    if( !mesh ) {
        throw std::runtime_error( "fix_degenerate_faces Error: mesh is NULL" );
    }

    polymesh3_builder builder;

    // Vertex data is unaffected, copy it directly
    for( std::size_t vertexIndex = 0, vertexCount = mesh->get_num_verts(); vertexIndex < vertexCount; ++vertexIndex ) {
        builder.add_vertex( mesh->get_vert( vertexIndex ) );
    }

    // Check and fix faces
    std::size_t newFaceCount = 0;
    std::size_t newFaceVertCount = 0;
    std::vector<std::pair<std::size_t, face_replacement>> changeset;

    std::vector<std::size_t> oldFace;
    std::vector<int> newFace;
    std::vector<std::size_t> replacementTopology;

    for( std::size_t faceIndex = 0, faceCount = mesh->get_num_faces(); faceIndex < faceCount; ++faceIndex ) {
        const std::size_t fvCount = mesh->get_num_face_verts( faceIndex );

        oldFace.resize( fvCount );
        mesh->get_face_vert_indices( faceIndex, &oldFace[0] );
        newFace.clear();
        replacementTopology.clear();

        // Copy the face while checking for duplicate vertices
        std::size_t fvIndex;
        for( fvIndex = 0; fvIndex < fvCount; ++fvIndex ) {
            int faceVert = static_cast<int>( oldFace[fvIndex] );
            std::vector<int>::iterator lastOccurrence = std::find( newFace.begin(), newFace.end(), faceVert );

            if( lastOccurrence == newFace.end() ) {
                newFace.push_back( faceVert );
            } else {
                // When a duplicate is found, populate a replacement topology and break into the replacement
                // construction loop
                replacementTopology.reserve( fvCount );
                for( std::size_t index = 0; index < fvIndex; ++index ) {
                    replacementTopology.push_back( index );
                }
                break;
            }
        }

        if( replacementTopology.size() > 0 ) {
            face_replacement replacement;

            // Note: This assumes that the polygon is convex. If a concave face has a duplicate vertex which closes an
            // internal loop, this will create an extra face for that loop rather than hollowing it out.
            for( ; fvIndex < fvCount; ++fvIndex ) {
                int faceVert = static_cast<int>( oldFace[fvIndex] );
                std::vector<int>::iterator lastOccurrence = std::find( newFace.begin(), newFace.end(), faceVert );

                if( lastOccurrence == newFace.end() ) {
                    newFace.push_back( faceVert );
                } else {
                    // Give the first occurrence of the duplicate vertex to the split face
                    std::vector<int> splitFace( lastOccurrence, newFace.end() );
                    newFace.resize( newFace.size() - splitFace.size() + 1 );

                    // Update the topology as we split the face
                    std::vector<std::size_t> splitFaceTopology(
                        replacementTopology.begin() + ( lastOccurrence - newFace.begin() ), replacementTopology.end() );
                    replacementTopology.resize( replacementTopology.size() - splitFaceTopology.size() );

                    if( splitFace.size() >= 3 ) {
                        builder.add_polygon( splitFace );
                        ++newFaceCount;
                        newFaceVertCount += splitFace.size();
                        replacement.push_back( splitFaceTopology );
                    }
                }

                replacementTopology.push_back( fvIndex );
            }

            if( newFace.size() >= 3 ) {
                builder.add_polygon( newFace );
                ++newFaceCount;
                newFaceVertCount += newFace.size();
                replacement.push_back( replacementTopology );
            }

            changeset.push_back( std::pair<std::size_t, face_replacement>( faceIndex, replacement ) );
        } else if( newFace.size() >= 3 ) {
            builder.add_polygon( newFace );
            ++newFaceCount;
            newFaceVertCount += newFace.size();
        }
    }

    polymesh3_ptr result = builder.finalize();

    // Reconstruct channels with new topology
    const mesh_interface::mesh_channel_map& vertexChannelMap = mesh->get_vertex_channels();
    for( mesh_interface::mesh_channel_map::const_iterator i = vertexChannelMap.begin(), ie = vertexChannelMap.end();
         i != ie; ++i ) {
        reconstruct_vertex_channel( i->second, changeset, newFaceVertCount, result );
    }

    const mesh_interface::mesh_channel_map& faceChannelMap = mesh->get_face_channels();
    for( mesh_interface::mesh_channel_map::const_iterator i = faceChannelMap.begin(), ie = faceChannelMap.end();
         i != ie; ++i ) {
        reconstruct_face_channel( i->second, changeset, newFaceCount, result );
    }

    return result;
}

void compute_vertex_normals( const frantic::geometry::mesh_interface* mesh,
                             std::vector<frantic::graphics::vector3f>& outVertexNormals ) {
    using namespace frantic::graphics;
    if( !mesh ) {
        throw std::runtime_error( "compute_vertex_normals Error: mesh is NULL" );
    }

    const std::size_t vertexCount = mesh->get_num_verts();
    const std::size_t faceCount = mesh->get_num_faces();

    outVertexNormals.clear();
    outVertexNormals.resize( vertexCount );

    std::vector<std::size_t> indexVector;
    std::vector<float> weightVector;

    // coped from frantic::geometry::trimesh3::build_vertex_normals()
    for( size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        const std::size_t cornerCount = mesh->get_num_face_verts( faceIndex );
        if( cornerCount < 3 ) {
            throw std::runtime_error(
                "calculate_normal_per_vertex Error: expected face with three or more vertices, but found face with " +
                boost::lexical_cast<std::string>( cornerCount ) + " vertices." );
        } else if( cornerCount == 3 ) {
            // Get the triangle
            std::size_t face[3];
            mesh->get_face_vert_indices( faceIndex, face );

            const vector3f a = mesh->get_vert( face[0] );
            const vector3f b = mesh->get_vert( face[1] );
            const vector3f c = mesh->get_vert( face[2] );

            // Get the normal, and use the triangle angles as the weights
            const vector3f geoNormal = triangle_normal( a, b, c );
            const vector3f weights = get_triangle_angles( a, b, c );
            outVertexNormals[face[0]] += weights.x * geoNormal;
            outVertexNormals[face[1]] += weights.y * geoNormal;
            outVertexNormals[face[2]] += weights.z * geoNormal;
        } else {
            // For now, using Newell's Method for polygons with more than three vertices, in
            // case they have collinear vertices.
            if( cornerCount > indexVector.size() ) {
                indexVector.resize( cornerCount );
                weightVector.resize( cornerCount );
            }
            mesh->get_face_vert_indices( faceIndex, &indexVector[0] );

            const vector3f first = mesh->get_vert( indexVector[0] );
            const vector3f last = mesh->get_vert( indexVector.back() );

            vector3f normal;
            vector3f previous;
            vector3f current = last;
            vector3f next = first;
            for( std::size_t cornerIndex = 0, cornerIndexEnd = cornerCount - 1; cornerIndex < cornerIndexEnd;
                 ++cornerIndex ) {
                previous = current;
                current = next;
                next = mesh->get_vert( indexVector[cornerIndex + 1] );

                normal.x += ( current.y - next.y ) * ( current.z + next.z );
                normal.y += ( current.z - next.z ) * ( current.x + next.x );
                normal.z += ( current.x - next.x ) * ( current.y + next.y );
                weightVector[cornerIndex] = get_triangle_angle( previous, current, next );
            }

            previous = current;
            current = next;
            next = first;

            normal.x += ( current.y - next.y ) * ( current.z + next.z );
            normal.y += ( current.z - next.z ) * ( current.x + next.x );
            normal.z += ( current.x - next.x ) * ( current.y + next.y );
            weightVector[cornerCount - 1] = get_triangle_angle( previous, current, next );

            const float length = normal.get_magnitude();
            if( length != 0 ) {
                normal /= length;
            }

            for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
                outVertexNormals[indexVector[cornerIndex]] += weightVector[cornerIndex] * normal;
            }
        }
    }

    BOOST_FOREACH( vector3f& v, outVertexNormals ) {
        v.normalize();
    }
}

void compute_face_normals( const frantic::geometry::mesh_interface* mesh,
                           std::vector<frantic::graphics::vector3f>& outFaceNormals ) {
    outFaceNormals.resize( mesh->get_num_faces() );

    for( size_t faceIndex = 0; faceIndex < mesh->get_num_faces(); ++faceIndex ) {
        outFaceNormals[faceIndex] =
            get_polygon_normal3( face_vertex_iterator( mesh, faceIndex, face_vertex_iterator::begin ),
                                 face_vertex_iterator( mesh, faceIndex, face_vertex_iterator::end ) );
    }
}

namespace {

class vertex_normal_mesh_channel : public frantic::geometry::mesh_channel {
  public:
    vertex_normal_mesh_channel( const mesh_interface* mesh )
        : mesh_channel( _T("Normal"), mesh_channel::vertex, frantic::channels::data_type_float32, 3,
                        mesh->get_num_verts(), mesh->get_num_faces(), true )
        , m_mesh( mesh ) {
        if( !m_mesh ) {
            throw std::runtime_error( "vertex_normal_mesh_channel() mesh is NULL" );
        }

        frantic::geometry::compute_vertex_normals( mesh, m_vertexNormals );
    }

    void get_value( std::size_t index, void* outValue ) const {
        assert( index < m_vertexNormals.size() );
        assert( outValue );

        memcpy( outValue, &m_vertexNormals[index], get_element_size() );
    }

    void set_value( std::size_t /*index*/, const void* /*value*/ ) const {
        // pass
    }

    std::size_t get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
        assert( m_mesh );
        return m_mesh->get_face_vert_index( faceIndex, fvertIndex );
    }

    std::size_t get_num_face_verts( std::size_t faceIndex ) const {
        assert( m_mesh );
        return m_mesh->get_num_face_verts( faceIndex );
    }

  private:
    const mesh_interface* m_mesh;
    std::vector<frantic::graphics::vector3f> m_vertexNormals;
};

class create_vertex_normal_channel_mesh_interface_impl : public delegated_const_mesh_interface {
  public:
    create_vertex_normal_channel_mesh_interface_impl( const mesh_interface* meshInterface )
        : delegated_const_mesh_interface( meshInterface ) {
        if( !meshInterface ) {
            throw std::runtime_error( "create_vertex_normal_channel_mesh_interface() mesh is NULL" );
        }

        // Add all delegate channels except the Normal channel to our mesh
        frantic::channels::channel_propagation_policy cpp;
        cpp.add_channel( _T("Normal") );

        append_delegate_channels( cpp );

        std::unique_ptr<vertex_normal_mesh_channel> normalChannel( new vertex_normal_mesh_channel( meshInterface ) );
        append_vertex_channel( std::move( normalChannel ) );
    }
};

} // anonymous namespace

std::unique_ptr<mesh_interface> create_vertex_normal_channel_mesh_interface( const mesh_interface* meshInterface ) {
    if( !meshInterface ) {
        throw std::runtime_error( "create_vertex_normal_channel_mesh_interface Error: mesh is NULL" );
    }

    std::unique_ptr<mesh_interface> result( new create_vertex_normal_channel_mesh_interface_impl( meshInterface ) );
    return result;
}

namespace {

class face_normal_mesh_channel : public frantic::geometry::mesh_channel {
  public:
    face_normal_mesh_channel( const mesh_interface* mesh )
        : mesh_channel( _T("Normal"), mesh_channel::face, frantic::channels::data_type_float32, 3,
                        mesh->get_num_faces(), 0, true )
        , m_mesh( mesh ) {
        if( !m_mesh ) {
            throw std::runtime_error( "face_normal_mesh_channel() mesh is NULL" );
        }

        frantic::geometry::compute_face_normals( mesh, m_faceNormals );
    }

    void get_value( std::size_t index, void* outValue ) const {
        assert( index < m_faceNormals.size() );
        assert( outValue );

        memcpy( outValue, &m_faceNormals[index], get_element_size() );
    }

    void set_value( std::size_t /*index*/, const void* /*value*/ ) const {
        throw std::runtime_error( "face_normal_mesh_channel::set_value: Error, this channel is read-only." );
    }

    // Technically these methods don't make any sense on a face channel.
    std::size_t get_fv_index( std::size_t /*faceIndex*/, std::size_t /*fvertIndex*/ ) const {
        throw std::runtime_error(
            "face_normal_mesh_channel::set_value: Error, this channel is a face-element channel." );
    }

    std::size_t get_num_face_verts( std::size_t /*faceIndex*/ ) const {
        throw std::runtime_error(
            "face_normal_mesh_channel::set_value: Error, this channel is a face-element channel." );
    }

  private:
    const mesh_interface* m_mesh;
    std::vector<frantic::graphics::vector3f> m_faceNormals;
};

class create_face_normal_channel_mesh_interface_impl : public delegated_const_mesh_interface {
  public:
    create_face_normal_channel_mesh_interface_impl( const mesh_interface* meshInterface )
        : delegated_const_mesh_interface( meshInterface ) {
        if( !meshInterface ) {
            throw std::runtime_error( "create_face_normal_channel_mesh_interface_impl() mesh is NULL" );
        }

        // Add all delegate channels except the Normal channel to our mesh
        frantic::channels::channel_propagation_policy cpp;
        cpp.add_channel( _T("Normal") );

        append_delegate_channels( cpp );

        std::unique_ptr<face_normal_mesh_channel> normalChannel( new face_normal_mesh_channel( meshInterface ) );
        append_face_channel( std::move( normalChannel ) );
    }
};

} // anonymous namespace

std::unique_ptr<frantic::geometry::mesh_interface>
create_face_normal_channel_mesh_interface( const frantic::geometry::mesh_interface* meshInterface ) {
    if( !meshInterface ) {
        throw std::runtime_error( "create_face_normal_channel_mesh_interface Error: mesh is NULL" );
    }

    std::unique_ptr<mesh_interface> result( new create_face_normal_channel_mesh_interface_impl( meshInterface ) );
    return result;
}

namespace {

#if defined( _MSC_VER )
#pragma warning( push )
// Warning: 'argument' : conversion from 'double' to 'float', possible loss of data
// Because the `half` class does not have a constructor from `double`.
#pragma warning( disable : 4244 )
#endif

template <class FloatType>
void compute_face_areas( const frantic::geometry::mesh_interface* mesh, std::vector<FloatType>& outFaceAreas ) {
    outFaceAreas.resize( mesh->get_num_faces() );

    for( size_t faceIndex = 0; faceIndex < mesh->get_num_faces(); ++faceIndex ) {
        face_vertex_range faceRange = make_face_vertex_range( mesh, faceIndex );
        outFaceAreas[faceIndex] = (FloatType)polygon3_area( faceRange.begin(), faceRange.end() );
    }
}

#if defined( _MSC_VER )
#pragma warning( pop )
#endif

template <class FloatType>
class face_area_mesh_channel : public frantic::geometry::mesh_channel {
  public:
    face_area_mesh_channel( const mesh_interface* mesh, const frantic::tstring& channelName )
        : mesh_channel( channelName, frantic::geometry::mesh_channel::face,
                        frantic::channels::channel_data_type_traits<FloatType>::data_type(), 1, mesh->get_num_faces(),
                        true )
        , m_faceAreas( mesh->get_num_faces() ) {
        if( !mesh ) {
            throw std::runtime_error( "create_face_normal_channel_mesh_interface_impl() mesh is NULL" );
        }

        compute_face_areas( mesh, m_faceAreas );
    }

    void get_value( std::size_t index, void* outValue ) const {
        assert( index < m_faceAreas.size() );
        assert( outValue );

        memcpy( outValue, &m_faceAreas[index], get_element_size() );
    }

    void set_value( std::size_t /*index*/, const void* /*value*/ ) const {
        throw std::runtime_error( "face_area_mesh_channel::set_value: Error, this channel is read-only." );
    }

    // Technically these methods don't make any sense on a face channel.
    std::size_t get_fv_index( std::size_t /*faceIndex*/, std::size_t /*fvertIndex*/ ) const {
        throw std::runtime_error( "face_area_mesh_channel::set_value: Error, this channel is a face-element channel." );
    }

    std::size_t get_num_face_verts( std::size_t /*faceIndex*/ ) const {
        throw std::runtime_error( "face_area_mesh_channel::set_value: Error, this channel is a face-element channel." );
    }

  private:
    const mesh_interface* m_mesh;
    std::vector<FloatType> m_faceAreas;
};

template <class FloatType>
class create_face_area_channel_mesh_interface_impl : public delegated_const_mesh_interface {
  public:
    create_face_area_channel_mesh_interface_impl( const mesh_interface* mesh, const frantic::tstring& channelName )
        : delegated_const_mesh_interface( mesh ) {
        if( !mesh ) {
            throw std::runtime_error( "create_face_normal_channel_mesh_interface_impl() mesh is NULL" );
        }

        // Add all delegate channels except the Normal channel to our mesh
        frantic::channels::channel_propagation_policy cpp;
        cpp.add_channel( channelName );

        append_delegate_channels( cpp );

        std::unique_ptr<face_area_mesh_channel<FloatType>> areaChannel(
            new face_area_mesh_channel<FloatType>( mesh, channelName ) );
        append_face_channel( std::move( areaChannel ) );
    }
};

template <class FloatType>
std::unique_ptr<mesh_interface>
create_face_area_channel_mesh_interface_template( const mesh_interface* mesh, const frantic::tstring& channelName ) {
    return std::unique_ptr<mesh_interface>(
        new create_face_area_channel_mesh_interface_impl<FloatType>( mesh, channelName ) );
}

} // namespace

std::unique_ptr<mesh_interface> create_face_area_channel_mesh_interface( const mesh_interface* mesh,
                                                                         frantic::channels::data_type_t dataType,
                                                                         const frantic::tstring& channelName ) {
    using namespace frantic::channels;

    switch( dataType ) {
    case data_type_float16:
        return create_face_area_channel_mesh_interface_template<half>( mesh, channelName );
    case data_type_float32:
        return create_face_area_channel_mesh_interface_template<float>( mesh, channelName );
    case data_type_float64:
        return create_face_area_channel_mesh_interface_template<double>( mesh, channelName );
    default:
        throw std::runtime_error( "frantic::geometry::create_face_area_channel_mesh_interface: Error, could not create "
                                  "face area channel of data type \"" +
                                  frantic::strings::to_string( channel_data_type_str( dataType ) ) +
                                  "\", must be a floating-point type." );
    }
}

namespace {

class create_modified_geometry_mesh_interface_impl : public delegated_const_mesh_interface {
  public:
    create_modified_geometry_mesh_interface_impl( const mesh_interface* meshInterface )
        : delegated_const_mesh_interface( meshInterface )
        , m_vertices( meshInterface->get_num_verts() ) {
        for( size_t vertexId = 0; vertexId < meshInterface->get_num_verts(); ++vertexId ) {
            m_vertices[vertexId] = meshInterface->get_vert( vertexId );
        }
    }

    virtual void get_vert( std::size_t index, float ( &outValues )[3] ) const {
        const frantic::graphics::vector3f& v = m_vertices.at( index );
        for( size_t i = 0; i < 3; ++i ) {
            outValues[i] = v[i];
        }
    }

    virtual void set_vert( std::size_t index, const float* v ) {
        m_vertices.at( index ) = frantic::graphics::vector3f( v[0], v[1], v[2] );
    }

    virtual bool is_read_only() const { return false; }

  private:
    // TODO: eventually the geometry channel should not be limited to only single-precision floats
    std::vector<frantic::graphics::vector3f> m_vertices;
};
} // namespace

std::unique_ptr<frantic::geometry::mesh_interface>
create_modified_geometry_mesh_interface( const frantic::geometry::mesh_interface* meshInterface ) {
    if( !meshInterface ) {
        throw std::runtime_error( "create_modified_geometry_mesh_interface: Error: mesh is NULL" );
    }

    std::unique_ptr<mesh_interface> result( new create_modified_geometry_mesh_interface_impl( meshInterface ) );
    return result;
}

namespace {

class particle_array_mesh_channel : public frantic::geometry::mesh_channel, boost::noncopyable {
  public:
    particle_array_mesh_channel( const frantic::particles::particle_array& particles,
                                 const frantic::channels::channel& ch )
        : frantic::geometry::mesh_channel( ch.name(), frantic::geometry::mesh_channel::vertex, ch.data_type(),
                                           ch.arity(), particles.particle_count(), 0, true )
        , m_particles( particles )
        , m_acc( particles.get_channel_map().get_general_accessor( ch.name() ) )

    {}

    void get_value( std::size_t index, void* outValue ) const {
        assert( index < m_particles.particle_count() );
        assert( outValue );

        memcpy( outValue, m_acc.get_channel_data_pointer( m_particles[index] ), get_element_size() );
    }

    void set_value( std::size_t /*index*/, const void* /*value*/ ) const {
        // pass
    }

    std::size_t get_fv_index( std::size_t /*faceIndex*/, std::size_t /*fvertIndex*/ ) const {
        throw std::runtime_error( "particle_array_mesh_channel::get_fv_index Error: not implemented" );
    }

    std::size_t get_num_face_verts( std::size_t /*faceIndex*/ ) const {
        throw std::runtime_error( "particle_array_mesh_channel::get_num_face_verts Error: not implemented" );
    }

  private:
    const frantic::particles::particle_array& m_particles;
    frantic::channels::channel_general_accessor m_acc;
};

class particle_array_mesh_interface : public frantic::geometry::mesh_interface {
  public:
    particle_array_mesh_interface( BOOST_RV_REF( frantic::particles::particle_array ) particles )
        : m_particles( particles ) {
        const frantic::channels::channel_map& channelMap = m_particles.get_channel_map();

        if( !channelMap.has_channel( _T("Position") ) ) {
            throw std::runtime_error( "particle_array_mesh_interface Error: "
                                      "particles do not have a \"Position\" channel." );
        }
        m_positionAcc = channelMap.get_accessor<frantic::graphics::vector3f>( _T("Position") );

        for( std::size_t i = 0, ie = channelMap.channel_count(); i < ie; ++i ) {
            const frantic::channels::channel& ch = channelMap[i];
            if( ch.name() == _T("Position") ) {
                continue;
            }
            std::unique_ptr<particle_array_mesh_channel> meshChannel(
                new particle_array_mesh_channel( m_particles, ch ) );
            append_vertex_channel( std::move( meshChannel ) );
        }
    }

    bool is_valid() const { return true; }

    bool request_channel( const frantic::tstring& channelName, bool vertexChannel, bool forOutput, bool throwOnError ) {
        if( forOutput ) {
            if( throwOnError ) {
                throw std::runtime_error( "particle_array_mesh_interface::request_channel() "
                                          "Cannot provide writable channel: "
                                          "\"" +
                                          frantic::strings::to_string( channelName ) + "\"" );
            }
            return false;
        }

        const mesh_channel* ch = 0;

        if( vertexChannel ) {
            ch = this->get_vertex_channels().get_channel( channelName );
        } else {
            ch = this->get_face_channels().get_channel( channelName );
        }

        if( ch ) {
            return true;
        } else {
            if( throwOnError ) {
                throw std::runtime_error( "particle_array_mesh_interface::request_channel() "
                                          "Cannot add channel: "
                                          "\"" +
                                          frantic::strings::to_string( channelName ) + "\"" );
            }
            return false;
        }
    }

    std::size_t get_num_verts() const { return m_particles.particle_count(); }

    void get_vert( std::size_t index, float ( &outValues )[3] ) const {
        assert( index < m_particles.particle_count() );
        assert( m_positionAcc.is_valid() );

        const frantic::graphics::vector3f position = m_positionAcc( m_particles[index] );
        memcpy( outValues, &position[0], 3 * sizeof( float ) );
    }

    std::size_t get_num_faces() const { return 0; }

    std::size_t get_num_face_verts( std::size_t /*faceIndex*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_num_face_verts Error: not implemented" );
    }

    std::size_t get_face_vert_index( std::size_t /*faceIndex*/, std::size_t /*fvertIndex*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_face_vert_index Error: not implemented" );
    }

    void get_face_vert_indices( std::size_t /*faceIndex*/, std::size_t /*outValues*/[] ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_face_vert_indices Error: not implemented" );
    }

    void get_face_verts( std::size_t /*faceIndex*/, float /*outValues*/[][3] ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_face_verts Error: not implemented" );
    }

    std::size_t get_num_elements() const { return 0; }

    std::size_t get_face_element_index( std::size_t /*faceIndex*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_face_element_index Error: not implemented" );
    }

    void init_adjacency() {
        throw std::runtime_error( "particle_array_mesh_interface::init_adjacency Error: not implemented" );
    }

    bool has_adjacency() const { return false; }

    bool init_vertex_iterator( frantic::geometry::vertex_iterator& /*vIt*/, std::size_t /*vertexIndex*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::init_vertex_iterator Error: not implemented" );
    }

    bool advance_vertex_iterator( frantic::geometry::vertex_iterator& /*vIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::advance_vertex_iterator Error: not implemented" );
    }

    std::size_t get_edge_endpoint( frantic::geometry::vertex_iterator& /*vIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_edge_endpoint Error: not implemented" );
    }

    std::size_t get_edge_left_face( frantic::geometry::vertex_iterator& /*vIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_edge_left_face Error: not implemented" );
    }

    std::size_t get_edge_right_face( frantic::geometry::vertex_iterator& /*vIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_edge_right_face Error: not implemented" );
    }

    bool is_edge_visible( frantic::geometry::vertex_iterator& /*vIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::is_edge_visible Error: not implemented" );
    }

    bool is_edge_boundary( frantic::geometry::vertex_iterator& /*vIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::is_edge_boundary Error: not implemented" );
    }

    void init_face_iterator( frantic::geometry::face_iterator& /*fIt*/, std::size_t /*faceIndex*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::init_face_iterator Error: not implemented" );
    }

    bool advance_face_iterator( frantic::geometry::face_iterator& /*fIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::advance_face_iterator Error: not implemented" );
    }

    std::size_t get_face_neighbor( frantic::geometry::face_iterator& /*fIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_face_neighbor Error: not implemented" );
    }

    std::size_t get_face_next_vertex( frantic::geometry::face_iterator& /*fIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_face_next_vertex Error: not implemented" );
    }

    std::size_t get_face_prev_vertex( frantic::geometry::face_iterator& /*fIt*/ ) const {
        throw std::runtime_error( "particle_array_mesh_interface::get_face_previous_vertex Error: not implemented" );
    }

  private:
    frantic::particles::particle_array m_particles;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_positionAcc;
};

} // anonymous namespace

std::unique_ptr<frantic::geometry::mesh_interface>
create_particle_array_mesh_interface( BOOST_RV_REF( frantic::particles::particle_array ) particles ) {
    std::unique_ptr<mesh_interface> result( new particle_array_mesh_interface( boost::move( particles ) ) );
    return result;
}

std::unique_ptr<mesh_interface> make_submesh( const mesh_interface* principleMesh, const std::vector<size_t>& faces ) {
    std::vector<int> vertexRemap( principleMesh->get_num_verts(), -1 );

    polymesh3_builder builder;

    std::vector<int> currentFaceIndices( 3 );

    size_t currentVertexCount = 0;

    for( size_t i = 0; i < faces.size(); ++i ) {
        const size_t currentFace = faces[i];
        const size_t faceSize = principleMesh->get_num_face_verts( currentFace );

        if( faceSize > currentFaceIndices.size() ) {
            currentFaceIndices.resize( faceSize );
        }

        for( size_t j = 0; j < faceSize; ++j ) {
            const size_t vertexIndex = principleMesh->get_face_vert_index( currentFace, j );

            // the first time we encounter a given vertex
            if( vertexRemap[vertexIndex] == -1 ) {

                // add it to the builder, and set its re-mapping index
                builder.add_vertex( principleMesh->get_vert( vertexIndex ) );

                vertexRemap[vertexIndex] = int( currentVertexCount );
                ++currentVertexCount;
            }

            currentFaceIndices[j] = vertexRemap[vertexIndex];
        }

        builder.add_polygon( &currentFaceIndices[0], faceSize );
    }

    return std::unique_ptr<mesh_interface>( polymesh3_interface::create_instance( builder.finalize() ).release() );
}

mesh_interface_geom_channel::mesh_interface_geom_channel( mesh_interface* mesh )
    : mesh_channel( _T( "verts" ), mesh_channel::element, frantic::channels::data_type_float32, 3,
                    mesh->get_num_verts(), 0, mesh->is_read_only() )
    , m_mesh( mesh ) {
    set_transform_type( mesh_channel::transform_type::point );
}

void mesh_interface_geom_channel::get_value( std::size_t index, void* outValue ) const {
    float* v = static_cast<float*>( outValue );
    frantic::graphics::vector3f vertex = m_mesh->get_vert( index );
    v[0] = vertex[0];
    v[1] = vertex[1];
    v[2] = vertex[2];
}

void mesh_interface_geom_channel::set_value( std::size_t index, const void* value ) const {
    mesh_interface* nonConst = const_cast<mesh_interface*>( m_mesh );
    nonConst->set_vert( index, static_cast<const float*>( value ) );
}

std::size_t mesh_interface_geom_channel::get_fv_index( std::size_t faceIndex, std::size_t fvertIndex ) const {
    return m_mesh->get_face_vert_index( faceIndex, fvertIndex );
}

std::size_t mesh_interface_geom_channel::get_num_face_verts( size_t faceIndex ) const {
    return m_mesh->get_num_face_verts( faceIndex );
}

bool is_triangle_mesh( const mesh_interface* mesh ) {
    const size_t numFaces = mesh->get_num_faces();

    for( size_t i = 0; i < numFaces; ++i ) {
        if( mesh->get_num_face_verts( i ) != 3 ) {
            return false;
        }
    }
    return true;
}

void assert_valid_indices( const mesh_interface* mesh, const frantic::tstring& nameForErrorMessage ) {
    typedef mesh_interface::mesh_channel_map channel_map_t;

    using frantic::strings::to_string;

    const std::string meshName = to_string( nameForErrorMessage );

    if( !mesh ) {
        throw std::runtime_error( "assert_valid_indices Error: mesh is NULL" );
    }

    const std::size_t faceCount = mesh->get_num_faces();
    const std::size_t vertexCount = mesh->get_num_verts();

    std::vector<std::size_t> indices;
    for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
        const std::size_t cornerCount = mesh->get_num_face_verts( faceIndex );
        if( cornerCount > 0 ) {
            indices.resize( cornerCount );

            mesh->get_face_vert_indices( faceIndex, &indices[0] );

            for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
                const std::size_t vertexIndex = indices[cornerIndex];
                if( vertexIndex >= vertexCount ) {
                    throw exception_stream() << "assert_valid_indices Error: "
                                             << "vertex index out of bounds "
                                             << "(" << vertexIndex << " >= " << vertexCount << ") "
                                             << "in " << meshName;
                }
            }
        }
    }

    const channel_map_t& vertexChannels = mesh->get_vertex_channels();
    for( channel_map_t::const_iterator i = vertexChannels.begin(), ie = vertexChannels.end(); i != ie; ++i ) {
        const mesh_channel& ch = *i->second;
        const mesh_channel::channel_type channelType = ch.get_channel_type();
        if( channelType == mesh_channel::vertex ) {
            const std::size_t elementCount = ch.get_num_elements();
            if( elementCount != vertexCount ) {
                throw exception_stream() << "assert_valid_indices Error: "
                                         << "mismatch between number of data elements and number of mesh vertices "
                                         << "(" << elementCount << "vs. " << vertexCount << ") "
                                         << "in vertex channel " << to_string( ch.get_name() ) << " "
                                         << "of " << meshName;
            }
        } else if( channelType == mesh_channel::face_vertex ) {
            if( ch.get_num_faces() != faceCount ) {
                throw exception_stream() << "assert_valid_indices Error: "
                                         << "mismatch between number of channel faces and number of mesh faces "
                                         << "(" << ch.get_num_faces() << " vs. " << faceCount << ") "
                                         << "in vertex channel " << to_string( ch.get_name() ) << ") "
                                         << "of " << meshName;
            }
            const std::size_t elementCount = ch.get_num_elements();
            for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
                const std::size_t cornerCount = ch.get_num_face_verts( faceIndex );
                for( std::size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex ) {
                    const std::size_t elementIndex = ch.get_fv_index( faceIndex, cornerIndex );
                    if( elementIndex >= elementCount ) {
                        throw exception_stream() << "assert_valid_indices Error: "
                                                 << "data index is out of bounds "
                                                 << "(" << elementIndex << " >= " << elementCount << ") "
                                                 << "in vertex channel " << to_string( ch.get_name() ) << " "
                                                 << "of " << meshName;
                    }
                }
            }
        } else {
            throw exception_stream() << "assert_valid_indices Error: "
                                     << "unexpected channel type (" << channelType << ") for vertex channel";
        }
    }

    const channel_map_t& faceChannels = mesh->get_face_channels();
    for( channel_map_t::const_iterator i = faceChannels.begin(), ie = faceChannels.end(); i != ie; ++i ) {
        const mesh_channel& ch = *i->second;
        const mesh_channel::channel_type channelType = ch.get_channel_type();
        if( channelType == mesh_channel::face ) {
            const std::size_t elementCount = ch.get_num_elements();
            if( elementCount != faceCount ) {
                throw exception_stream() << "assert_valid_indices Error: "
                                         << "mismatch between number of data elements and number of mesh faces "
                                         << "(" << elementCount << " vs. " << faceCount << ")"
                                         << "in face channel " << to_string( ch.get_name() ) << " "
                                         << "of " << meshName;
            }
        } else {
            throw exception_stream() << "assert_valid_indices Error: "
                                     << "unexpected channel type (" << channelType << ") for face channel";
        }
    }
}

void get_connected_components_segmentation( const mesh_interface* iMesh, range_segmentation& outSegmentation ) {
    std::vector<size_t> segmentLabels;
    get_face_connected_components( iMesh, segmentLabels );
    outSegmentation.assign( segmentLabels.begin(), segmentLabels.end() );
}

void get_channel_segmentation( const mesh_interface* mesh, const tstring& customChannelName,
                               range_segmentation& outSegmentation ) {

    if( !mesh->get_vertex_channels().has_channel( customChannelName ) ) {
        throw std::runtime_error( "get_channel_segmentation : Error, mesh does not have channel" +
                                  strings::to_string( customChannelName ) );
    }

    const mesh_channel* customChannel = mesh->get_vertex_channels().get_channel( customChannelName );

    get_channel_segmentation( mesh, customChannel, outSegmentation );
}

void get_channel_segmentation( const mesh_interface* mesh, const mesh_channel* customChannel,
                               range_segmentation& outSegmentation ) {

    if( customChannel->get_channel_type() != mesh_channel::vertex &&
        customChannel->get_channel_type() != mesh_channel::face_vertex ) {
        throw std::runtime_error( "get_channel_segmentation : Error, mesh channel" +
                                  strings::to_string( customChannel->get_name() ) + " was not a vertex channel." );
    }

    std::vector<size_t> parents( mesh->get_num_faces() );
    std::vector<size_t> ranks( mesh->get_num_faces() );

    boost::disjoint_sets<std::size_t*, std::size_t*> chartSets( &ranks[0], &parents[0] );

    for( size_t i = 0; i < mesh->get_num_faces(); ++i ) {
        chartSets.make_set( i );
    }

    std::vector<size_t> faceVerts( 3 );

    for( size_t i = 0; i < mesh->get_num_faces(); ++i ) {
        const size_t faceSize = mesh->get_num_face_verts( i );

        if( faceVerts.size() < faceSize ) {
            faceVerts.resize( faceSize );
        }

        for( size_t j = 0; j < faceSize; ++j ) {
            faceVerts[j] = customChannel->get_fv_index( i, j );
        }

        std::sort( &faceVerts[0], &faceVerts[faceSize] );

        face_iterator fIt;
        mesh->init_face_iterator( fIt, i );

        do {
            const size_t neighbourFace = mesh->get_face_neighbor( fIt );

            // only process this face if its not a boundary, its index is lower than the current face, and they aren't
            // already known to be in the same chart
            if( neighbourFace != mesh_interface::HOLE_INDEX && neighbourFace < i &&
                chartSets.find_set( i ) != chartSets.find_set( neighbourFace ) ) {
                const size_t oppositeFaceSize = mesh->get_num_face_verts( neighbourFace );

                for( size_t j = 0; j < oppositeFaceSize; ++j ) {
                    if( std::binary_search( &faceVerts[0], &faceVerts[faceSize],
                                            customChannel->get_fv_index( neighbourFace, j ) ) ) {
                        chartSets.union_set( i, neighbourFace );
                        break;
                    }
                }
            }
        } while( mesh->advance_face_iterator( fIt ) );
    }

    std::vector<size_t> components( mesh->get_num_faces() );

    for( size_t i = 0; i < mesh->get_num_faces(); ++i ) {
        components[i] = chartSets.find_set( i );
    }

    outSegmentation.assign( components );
}

} // namespace geometry
} // namespace frantic
