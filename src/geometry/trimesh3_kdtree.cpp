// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>

//#define DEBUG_TRIMESH3_KDTREE

using namespace std;

namespace frantic {
namespace geometry {

frantic::diagnostics::profiling_section p_sorting( _T("Sorting") ), p_triangleBox( _T("Triangle Box Intersection") ),
    p_total( _T("Total") );

trimesh3_kdtree::trimesh3_kdtree( trimesh3& mesh )
    : m_mesh( &mesh ) {
    if( m_mesh == 0 )
        throw runtime_error(
            "trimesh3_kdtree constructor: The mesh passed to the constructor was null, a condition which "
            "should never occur." );

    // Make sure the mesh has the required helper data structures computed.
    m_mesh->build_face_planes( false );
    m_mesh->build_face_bound_boxes( false );
    m_mesh->build_barycentric_helpers( false );

    m_rootBounds = m_mesh->compute_bound_box();
    m_rootBounds.expand_fractional( 0.00001f );

    vector<int> allIntersectableFaces;
    allIntersectableFaces.reserve( m_mesh->face_count() );

    for( size_t f = 0; f < m_mesh->face_count(); ++f ) {
        vector3 face = m_mesh->get_face( f );
        if( face.x != face.y && face.x != face.z && face.y != face.z )
            allIntersectableFaces.push_back( (int)f );
    }

    p_sorting.reset();
    p_triangleBox.reset();
    p_total.reset();

    p_total.enter();
    m_rootNode = new trimesh3_kdtree_node();
    // construct_kdtree_n_log2_n( allIntersectableFaces, m_rootBounds, 0, *m_rootNode );
    build_kdtree_greedy_SAH_nlogn( *m_rootNode, *m_mesh, m_rootBounds );
    p_total.exit();

    //		cout << p_sorting << endl;
    //		cout << p_triangleBox << endl;
    //		cout << p_total << endl;
}

//////////////////////////////////////////////////////////
// O(n log^2 n) algorithm for SAH tree generation
//////////////////////////////////////////////////////////

// A structure holding a split candidate.  It will sort in the correct order for the sweep algorithm of finding the best
// split.
struct split_candidate {
    // The position at which the split candidate sits
    float splitValue;
    // Whether the split has the new triangle to the right (+1, "start event"), left (-1, "end event"), or the
    // triangle lies entirely on the plane of the split (0, "planar event")
    int splitType;
    // Which face generated the event
    int splitFaceIndex;

    split_candidate( float split, int type, int faceIndex ) {
        splitValue = split;
        splitType = type;
        splitFaceIndex = faceIndex;
    }

    bool operator<( const split_candidate& rhs ) const {
        return splitValue < rhs.splitValue || ( splitValue == rhs.splitValue && splitType < rhs.splitType );
    }
};

ostream& operator<<( ostream& out, const split_candidate& sc ) {
    out << "split: " << sc.splitValue << ", type: " << sc.splitType << ", faceindex: " << sc.splitFaceIndex;
    return out;
}

// This gets all the potential splits based on the primitives
int get_split_candidates( const trimesh3& mesh, const vector<int>& faceIndices, const boundbox3f& box,
                          vector<split_candidate>* splitCandidatesArray ) {
    // Allocate enough memory for our split candidates.
    unsigned expectedCandidateCount = 2 * (unsigned)faceIndices.size();
    splitCandidatesArray[0].reserve( expectedCandidateCount );
    splitCandidatesArray[1].reserve( expectedCandidateCount );
    splitCandidatesArray[2].reserve( expectedCandidateCount );

    boundbox3f faceBounds;
    int totalFaceCount = 0;
    p_triangleBox.enter();
    for( unsigned i = 0; i < faceIndices.size(); ++i ) {
        int faceIndex = faceIndices[i];
        vector3 face = mesh.get_face( faceIndex );
        vector3f vert0 = mesh.get_vertex( face.x ), vert1 = mesh.get_vertex( face.y ),
                 vert2 = mesh.get_vertex( face.z );
        if( box.intersect_with_triangle( vert0, vert1, vert2, faceBounds ) ) {
            ++totalFaceCount;
            for( int splitAxis = 0; splitAxis < 3; ++splitAxis ) {
                float splitCandidateA = faceBounds.minimum()[splitAxis],
                      splitCandidateB = faceBounds.maximum()[splitAxis];
                if( splitCandidateA == splitCandidateB ) {
                    splitCandidatesArray[splitAxis].push_back( split_candidate( splitCandidateA, 0, faceIndex ) );
                } else {
                    splitCandidatesArray[splitAxis].push_back( split_candidate( splitCandidateA, +1, faceIndex ) );
                    splitCandidatesArray[splitAxis].push_back( split_candidate( splitCandidateB, -1, faceIndex ) );
                }
            }
        }
    }
    p_triangleBox.exit();

    p_sorting.enter();
    for( int splitAxis = 0; splitAxis < 3; ++splitAxis ) {
        sort( splitCandidatesArray[splitAxis].begin(), splitCandidatesArray[splitAxis].end() );
    }
    p_sorting.exit();
    return totalFaceCount;
}

// Recursively generates the kd-tree from the mesh.  For efficiency, this destructively swaps out the passed
// in faceIndices array if it is needed.
void trimesh3_kdtree::construct_kdtree_n_log2_n( vector<int>& faceIndices, const boundbox3f& bounds, int depth,
                                                 trimesh3_kdtree_node& outNode ) const {
    // An arbitrarily small number of faces provides the cutoff point for subdivision
    // int leafSizeCutoff = 2;
    // if( m_mesh->face_count() > 100000 )
    //	leafSizeCutoff = 10;

    if( (int)faceIndices.size() <= 8 ) {
        outNode.initialize( faceIndices );
        return;
    }

    //		if(depth == 4 && faceIndices.size() == 3 && faceIndices[0] == 4 && faceIndices[1] == 0 && faceIndices[2] ==
    //3) 			depth = 4;

    // The 3 arrays of split candidates for X, Y, and Z.
    // The face count which actually intersected the bound box could conceivably be less than faceIndices.count(), so we
    // get that value here as well.
    vector<split_candidate> splitCandidates[3];
    int totalFaceCount = get_split_candidates( *m_mesh, faceIndices, bounds, splitCandidates );

    /*
    if( depth > 100 ) {
      for( int axis = 0; axis < 3; ++axis ) {
        cout << "Axis: " << axis << endl;
        for( unsigned i = 0; i < splitCandidates[axis].size(); ++i ) {
          cout << splitCandidates[axis][i] << endl;
        }
      }
      exit(1);
    }
    */

    // The best split value found.
    float bestSplitValue = 0;
    // The axis of the best split value.
    int bestAxis = -1;
    // The cost of the best split value.
    float bestCost = ( numeric_limits<float>::max )();
    // Which side of the tree to put the perp faces on.
    int bestPerpSide = 0;

    // Now do the sweep through the split candidates to find the best sweep.
    for( int testSplitAxis = 0; testSplitAxis < 3; ++testSplitAxis ) {
#ifdef DEBUG_TRIMESH3_KDTREE
        if( depth > 100 )
            cout << "Axis: " << testSplitAxis << endl;
#endif // DEBUG_TRIMESH3_KDTREE
       //  Only do this if the dimensions of the bounding box along this axis are positive.
        if( bounds.maximum()[testSplitAxis] - bounds.minimum()[testSplitAxis] > 0 ) {
            vector<split_candidate>& candidates = splitCandidates[testSplitAxis];
            // Start with all the faces on the right
            int leftFaceCount = 0, perpFaceCount = 0, rightFaceCount = totalFaceCount;
            boundbox3f leftBounds = bounds, rightBounds = bounds;
            unsigned i = 0, candidatesSize = (unsigned)candidates.size();
            while( i < candidatesSize ) {
                float splitValue = candidates[i].splitValue;
                int endEventCount = 0, perpEventCount = 0, startEventCount = 0;
                // Go through all the events coincident with this position.  This leaves i pointing at the next split
                // candidate, or past the end.
                while( i < candidatesSize && candidates[i].splitValue == splitValue && candidates[i].splitType == -1 ) {
                    ++endEventCount;
                    ++i;
                }
                while( i < candidatesSize && candidates[i].splitValue == splitValue && candidates[i].splitType == 0 ) {
                    ++perpEventCount;
                    ++i;
                }
                while( i < candidatesSize && candidates[i].splitValue == splitValue && candidates[i].splitType == +1 ) {
                    ++startEventCount;
                    ++i;
                }
                // Adjust all the face counts
                perpFaceCount = perpEventCount;
                rightFaceCount -= perpEventCount;
                rightFaceCount -= endEventCount;
                // If a split is at the edge of the box, and there are no faces that could go in this empty box, then
                // don't consider its cost
                if( !( ( splitValue <= bounds.minimum()[testSplitAxis] && ( leftFaceCount + perpFaceCount == 0 ) ) ||
                       ( splitValue >= bounds.maximum()[testSplitAxis] &&
                         ( rightFaceCount + perpFaceCount == 0 ) ) ) ) {
                    // Adjust the bounding boxes for this split
                    leftBounds.maximum()[testSplitAxis] = splitValue;
                    rightBounds.minimum()[testSplitAxis] = splitValue;
                    // Compute the surface area heuristic for this split, and check whether it's better than our current
                    // best one.
                    int perpSide = 0;
                    float cost = greedy_sah_cost( bounds.get_surface_area(), leftBounds.get_surface_area(),
                                                  rightBounds.get_surface_area(), leftFaceCount, perpFaceCount,
                                                  rightFaceCount, perpSide );
                    if( cost < bestCost ) {
                        bestSplitValue = splitValue;
                        bestAxis = testSplitAxis;
                        bestCost = cost;
                        bestPerpSide = perpSide;
                    }

#ifdef DEBUG_TRIMESH3_KDTREE
                    if( depth > 100 )
                        cout << "split: " << splitValue << ", cost: " << cost << ", left: " << leftFaceCount
                             << ", perp: " << perpFaceCount << ", right: " << rightFaceCount << ", side: " << perpSide
                             << ", sal: " << leftBounds.get_surface_area()
                             << ", sar: " << rightBounds.get_surface_area() << "sa: " << bounds.get_surface_area()
                             << endl;
#endif // DEBUG_TRIMESH3_KDTREE
                }
                leftFaceCount += startEventCount;
                leftFaceCount += perpEventCount;
                perpFaceCount = 0;
            }
        }
    }

#ifdef DEBUG_TRIMESH3_KDTREE
    if( depth > 100 ) {
        cout << "best split: " << bestSplitValue << endl;
        cout << "cost: " << bestCost << endl;
        cout << "perp side: " << bestPerpSide << endl;
        cout << "axis: " << bestAxis << endl;
        cout << "face indices size: " << faceIndices.size() << endl;
        // exit(1);
    }
#endif // DEBUG_TRIMESH3_KDTREE

    //		cout << "bestCost " << bestCost << " costThreshold " << faceIndices.size() << endl;
    // If we didn't find an acceptable split, (i.e. just running through all the faces costs less), then return a leaf
    // node. The 1.01f is in here because floating point error was causing the cost to be less than the size, even when
    // it should have been equal to it.
    if( bestCost * 1.01f >= faceIndices.size() ) {
        outNode.initialize( faceIndices );
        return;
    }

    // Compute the bounding boxes of both sides
    vector3f newLeftMaximum = bounds.maximum();
    newLeftMaximum[bestAxis] = bestSplitValue;
    boundbox3f boundsLeft( bounds.minimum(), newLeftMaximum );

    vector3f newRightMinimum = bounds.minimum();
    newRightMinimum[bestAxis] = bestSplitValue;
    boundbox3f boundsRight( newRightMinimum, bounds.maximum() );

    vector<int> leftFaceIndices, rightFaceIndices;
    vector<split_candidate>& candidates = splitCandidates[bestAxis];
    unsigned candidatesSize = (unsigned)candidates.size();
    // Now go through the split candidate array, and collect the face indices
    unsigned i = 0;
    while( i < candidatesSize && candidates[i].splitValue < bestSplitValue ) {
        // Add all faces that had a start event or perp event in this range to the left face list
        if( candidates[i].splitType >= 0 )
            leftFaceIndices.push_back( candidates[i].splitFaceIndex );
        ++i;
    }
    while( i < candidatesSize && candidates[i].splitValue == bestSplitValue ) {
        // Add all the faces that had a perp event to the left or right event as determined by the SAHs
        if( candidates[i].splitType == 0 ) {
            if( bestPerpSide < 0 ) {
                leftFaceIndices.push_back( candidates[i].splitFaceIndex );
            } else {
                rightFaceIndices.push_back( candidates[i].splitFaceIndex );
            }
        }
        ++i;
    }
    while( i < candidatesSize ) {
        // Add all the faces that had an end event or perp event in this range to the right face list
        if( candidates[i].splitType <= 0 )
            rightFaceIndices.push_back( candidates[i].splitFaceIndex );
        ++i;
    }

    // Now free the memory for the candidate splits.
    for( int axis = 0; axis < 3; ++axis ) {
        vector<split_candidate> freeArray;
        splitCandidates[axis].swap( freeArray );
    }

    //		cout << "Splitting interior node at depth " << depth << " with bounds " << bounds << ", axis " << bestAxis <<
    //", split " << bestSplitValue << ", perpSide " << bestPerpSide << endl;

    // These recursive calls to construct_kdtree may destructively swap out the leftFaceIndices and rightFaceIndices
    // contents, so those vectors aren't valid after this call.
    //		cout << "creating left child at depth " << depth+1 << " with bounds " << boundsLeft << ", " <<
    // leftFaceIndices.size() << " nodes " << endl;
    trimesh3_kdtree_node* children = new trimesh3_kdtree_node[2];
    construct_kdtree_n_log2_n( leftFaceIndices, boundsLeft, depth + 1, children[0] );
    //		cout << "creating right child at depth " << depth+1 << " with bounds " << boundsRight << ", " <<
    // rightFaceIndices.size() << " nodes " << endl;
    construct_kdtree_n_log2_n( rightFaceIndices, boundsRight, depth + 1, children[1] );

    outNode.initialize( children, bestAxis, bestSplitValue );
}

} // namespace geometry
} // namespace frantic
