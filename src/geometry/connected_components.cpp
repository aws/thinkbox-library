// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/connected_components.hpp>

#include <frantic/logging/logging_level.hpp>
#include <frantic/misc/algorithm.hpp>

#include <boost/pending/disjoint_sets.hpp>

#include <iterator>

namespace frantic {
namespace geometry {

void get_face_connected_components( const frantic::geometry::mesh_interface_ptr mesh,
                                    std::vector<size_t>& outFaceLabels ) {
    get_face_connected_components( mesh.get(), outFaceLabels );
}

void get_face_connected_components( const frantic::geometry::mesh_interface* mesh,
                                    std::vector<size_t>& outFaceLabels ) {
    const size_t numFaces = mesh->get_num_faces();

    if( !mesh->has_adjacency() ) {
        throw std::runtime_error(
            "get_face_connected_components: Error, mesh must have connectivity information to get "
            "face connected components" );
    }

    outFaceLabels.resize( numFaces );

    std::vector<size_t> ranks( numFaces );
    std::vector<size_t> parents( numFaces );
    boost::disjoint_sets<size_t*, size_t*> disjointSets( &ranks[0], &parents[0] );

    for( size_t i = 0; i < numFaces; ++i ) {
        disjointSets.make_set( i );
    }

    for( size_t i = 0; i < numFaces; ++i ) {
        frantic::geometry::face_iterator fIt;
        mesh->init_face_iterator( fIt, i );

        do {
            size_t oppositeFace = mesh->get_face_neighbor( fIt );
            if( oppositeFace != frantic::geometry::mesh_interface::HOLE_INDEX ) {
                disjointSets.union_set( i, oppositeFace );
            }
        } while( mesh->advance_face_iterator( fIt ) );
    }

    for( size_t i = 0; i < numFaces; ++i ) {
        outFaceLabels[i] = disjointSets.find_set( i );
    }
}

void group_face_components( const std::vector<size_t>& outFaceLabels, std::vector<size_t>& outFaceLists,
                            std::vector<size_t>& outListOffsets ) {

    outFaceLists.clear();
    outFaceLists.reserve( outFaceLabels.size() );
    outListOffsets.clear();

    frantic::group_components( outFaceLabels.begin(), outFaceLabels.end(), std::back_inserter( outFaceLists ),
                               std::back_inserter( outListOffsets ) );
}

void separate_segmentation_components( const dcel& edgeStructure, std::vector<size_t>& faceLabels ) {
    const size_t numFaces = faceLabels.size();

    std::vector<size_t> ranks( numFaces );
    std::vector<size_t> parents( numFaces );
    boost::disjoint_sets<size_t*, size_t*> disjointSets( &ranks[0], &parents[0] );

    for( size_t faceId = 0; faceId < numFaces; ++faceId ) {
        disjointSets.make_set( faceId );
    }

    for( size_t faceId = 0; faceId < numFaces; ++faceId ) {
        dcel::const_halfedge_handle edgeHandle = edgeStructure.get_face_halfedge( faceId );

        do {
            size_t oppositeFace = edgeHandle.opposite_face();
            if( oppositeFace < numFaces && faceLabels[faceId] == faceLabels[oppositeFace] ) {
                disjointSets.union_set( faceId, oppositeFace );
            }
            edgeHandle = edgeHandle.face_next();
        } while( edgeHandle != edgeStructure.get_face_halfedge( faceId ) );
    }

    for( size_t faceId = 0; faceId < numFaces; ++faceId ) {
        faceLabels[faceId] = disjointSets.find_set( faceId );
    }
}

void separate_segmentation_components( const mesh_interface* mesh, const range_segmentation& initialSegmentation,
                                       range_segmentation& outSegmentation ) {
    using namespace graphics;
    using namespace geometry;

    if( !mesh->has_adjacency() ) {
        throw std::runtime_error(
            "clean_segmentation: Error, mesh must have adjacency information to use this method." );
    }

    const size_t numFaces = mesh->get_num_faces();

    std::vector<size_t> ranks( numFaces );
    std::vector<size_t> parents( numFaces );
    boost::disjoint_sets<size_t*, size_t*> disjointSets( &ranks[0], &parents[0] );

    for( size_t faceId = 0; faceId < numFaces; ++faceId ) {
        disjointSets.make_set( faceId );
    }

    for( size_t faceId = 0; faceId < numFaces; ++faceId ) {
        face_iterator fIt;
        mesh->init_face_iterator( fIt, faceId );

        do {
            size_t oppositeFace = mesh->get_face_neighbor( fIt );
            if( oppositeFace != mesh_interface::HOLE_INDEX &&
                initialSegmentation.get_face_subset( faceId ) == initialSegmentation.get_face_subset( oppositeFace ) ) {
                disjointSets.union_set( faceId, oppositeFace );
            }
        } while( mesh->advance_face_iterator( fIt ) );
    }

    std::vector<size_t> faceLabels( numFaces );

    for( size_t faceId = 0; faceId < numFaces; ++faceId ) {
        faceLabels[faceId] = disjointSets.find_set( faceId );
    }

    outSegmentation.assign( faceLabels.begin(), faceLabels.end() );
}

void separate_segmentation_components( const mesh_interface* mesh, range_segmentation& inoutSegmentation ) {
    separate_segmentation_components( mesh, inoutSegmentation, inoutSegmentation );
}

} // namespace geometry
} // namespace frantic
