// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// This class provides a kdtree acceleration structure around child primitives that expose a ray intersection interface
// as well as an axis-aligned bounding box.  An example of that is the trimesh3_kdtree, and this class itself.  Nesting
// primitive_kdtree objects isn't a good idea, though.

// PrimitiveType requirements:
//
// Must have the following functions
//
// bool intersect_ray( const ray3f& ray, float segmentLength, raytrace_intersection& outIntersection ) const;
// void intersect_ray_all( const ray3f& ray, float segmentLength, std::vector<raytrace_intersection>& outIntersections )
// const; bool intersects_ray_segment( const vector3f& start, const vector3f& end ) const; bool find_nearest_point(
// const vector3f& point, float maxDistance, nearest_point_search_result& outNearestPoint ) const;
//
// boundbox3f get_bounds() const;

//#define DEBUG_PRIMITIVE_KDTREE

#pragma once

#include <frantic/diagnostics/profiling_section.hpp>
#include <frantic/geometry/primitive_kdtree_node.hpp>

namespace frantic {
namespace geometry {

template <class PrimitiveType>
class primitive_kdtree {
    std::vector<boost::shared_ptr<PrimitiveType>> m_primitives;
    std::vector<boundbox3f> m_primitiveBounds;
    boundbox3f m_rootBounds;
    primitive_kdtree_node<PrimitiveType>* m_rootNode;

  public:
    primitive_kdtree( const std::vector<boost::shared_ptr<PrimitiveType>>& primitives )
        : m_primitives( primitives ) {
        // Get all the bounding boxes from the primitives, and build up the bounding box of all the geometry
        m_primitiveBounds.reserve( primitives.size() );
        for( unsigned i = 0; i < primitives.size(); ++i ) {
            m_primitiveBounds.push_back( primitives[i]->get_bounds() );
            m_rootBounds += m_primitiveBounds.back();
        }

        // Initialize an array of all the primitives.
        std::vector<int> allIntersectablePrimitives;
        allIntersectablePrimitives.reserve( m_primitives.size() );
        for( unsigned i = 0; i < m_primitives.size(); ++i ) {
            allIntersectablePrimitives.push_back( i );
        }

        // Build the kd-tree
        m_rootNode = new primitive_kdtree_node<PrimitiveType>();
        construct_kdtree_n_log2_n( allIntersectablePrimitives, m_rootBounds, 0, *m_rootNode );

        std::cout << "Created primitive kd-tree: " << std::endl;
        std::cout << " node count: " << get_node_count() << std::endl;
        std::cout << " max depth: " << get_maximum_depth() << std::endl;
        std::cout << " largest leaf: " << get_largest_leaf_size() << std::endl;

        //		dump_tree(std::cout);
        // exit(1);
    }

    ~primitive_kdtree() {
        if( m_rootNode != 0 ) {
            delete m_rootNode;
            m_rootNode = 0;
        }
    }

    // Returns the total number of nodes in the tree.
    int get_node_count() const { return m_rootNode->get_node_count(); }

    // Returns the depth of the deepest leaf node.
    int get_maximum_depth() const { return m_rootNode->get_maximum_depth(); }

    // Returns the largest count of primitives in a leaf node.
    int get_largest_leaf_size() const { return m_rootNode->get_largest_leaf_size(); }

    void dump_tree( std::ostream& out ) const {
        out << "----- Dumping Primitive KD-Tree\n";
        out << " boundbox: " << m_rootBounds << "\n";
        m_rootNode->dump_tree( out, 0, m_primitives );
        out << "----- Finished dumping Primitive KD-Tree" << std::endl;
    }

    ///// Queries/////

    bool is_point_in_volume( const vector3f& point ) const {
        return m_rootNode->is_point_in_volume( point, m_primitives );
    }

    bool find_nearest_point( const vector3f& point, float maxDistance,
                             nearest_point_search_result& outNearestPoint ) const {
        boundbox3f bounds( m_rootBounds );

        return m_rootNode->find_nearest_point( point, maxDistance, bounds, m_primitives, outNearestPoint );
    }

    // This returns the first intersection between the ray and one of the primitives, or false if none is found.
    bool intersect_ray( const ray3f& ray, double tMin, double tMax, raytrace_intersection& outIntersection ) const {

        return m_rootNode->intersect_ray( ray, tMin, tMax, m_primitives, outIntersection );
    }

    // This adds the intersections found to the array, and maintains its sorted order
    void intersect_ray_all( const ray3f& ray, float segmentLength,
                            std::vector<raytrace_intersection>& outIntersections ) const {
        m_rootNode->intersect_ray_all( ray, 0, segmentLength, m_primitives, outIntersections );
        std::sort( outIntersections.begin(), outIntersections.end() );
    }

    // This returns whether or not the ray intersects with the mesh within the given distance.
    bool intersects_ray_segment( const vector3f& start, const vector3f& end ) const {
        ray3f ray( start, end - start );
        return m_rootNode->intersects_ray_segment( ray, 0, 1, m_primitives );
    }

    // This returns whether or not the ray intersects with the mesh within the given distance.
    bool intersects_ray_segment( const ray3f& ray, double tMin, double tMax ) const {
        return m_rootNode->intersects_ray_segment( ray, tMin, tMax, m_primitives );
    }

    //////////////////////////////////////////////////////////
    // O(n log^2 n) algorithm for SAH tree generation
    //////////////////////////////////////////////////////////

    // A structure holding a split candidate.  It will sort in the correct order for the sweep algorithm of finding the
    // best split.
    struct split_candidate {
        // The position at which the split candidate sits
        float splitValue;
        // Whether the split has the new triangle to the right (+1, "start event"), left (-1, "end event"), or the
        // triangle lies entirely on the plane of the split (0, "planar event")
        int splitType;
        // Which face generated the event
        int splitPrimitiveIndex;

        split_candidate( float split, int type, int primitiveIndex ) {
            splitValue = split;
            splitType = type;
            splitPrimitiveIndex = primitiveIndex;
        }

        bool operator<( const split_candidate& rhs ) const {
            return splitValue < rhs.splitValue || ( splitValue == rhs.splitValue && splitType < rhs.splitType );
        }

        void print( std::ostream& out ) const {
            out << "split: " << splitValue << ", type: " << splitType << ", primitiveindex: " << splitPrimitiveIndex
                << std::endl;
        }
    };

    // This gets all the potential splits based on the primitives
    static int get_split_candidates( const std::vector<boundbox3f>& primitiveBoundsArray,
                                     const std::vector<int>& primitiveIndices, const boundbox3f& box,
                                     std::vector<split_candidate>* splitCandidatesArray ) {
        // Allocate enough memory for our split candidates.
        unsigned expectedCandidateCount = 2 * (unsigned)primitiveIndices.size();
        splitCandidatesArray[0].reserve( expectedCandidateCount );
        splitCandidatesArray[1].reserve( expectedCandidateCount );
        splitCandidatesArray[2].reserve( expectedCandidateCount );

        boundbox3f primitiveBounds;
        int totalPrimitiveCount = 0;
        for( unsigned i = 0; i < primitiveIndices.size(); ++i ) {
            int primitiveIndex = primitiveIndices[i];
            primitiveBounds = primitiveBoundsArray[primitiveIndex];
            primitiveBounds.intersect_with( box );
            if( !primitiveBounds.is_empty() ) {
                ++totalPrimitiveCount;
                for( int splitAxis = 0; splitAxis < 3; ++splitAxis ) {
                    float splitCandidateA = primitiveBounds.minimum()[splitAxis],
                          splitCandidateB = primitiveBounds.maximum()[splitAxis];
                    if( splitCandidateA == splitCandidateB ) {
                        splitCandidatesArray[splitAxis].push_back(
                            split_candidate( splitCandidateA, 0, primitiveIndex ) );
                    } else {
                        splitCandidatesArray[splitAxis].push_back(
                            split_candidate( splitCandidateA, +1, primitiveIndex ) );
                        splitCandidatesArray[splitAxis].push_back(
                            split_candidate( splitCandidateB, -1, primitiveIndex ) );
                    }
                }
            }
        }

        for( int splitAxis = 0; splitAxis < 3; ++splitAxis ) {
            std::sort( splitCandidatesArray[splitAxis].begin(), splitCandidatesArray[splitAxis].end() );
        }
        return totalPrimitiveCount;
    }

    // Recursively generates the kd-tree from the mesh.  For efficiency, this destructively swaps out the passed
    // in primitiveIndices array if it is needed.
    void construct_kdtree_n_log2_n( std::vector<int>& primitiveIndices, const boundbox3f& bounds, int depth,
                                    primitive_kdtree_node<PrimitiveType>& outNode ) const {
        // An arbitrarily small number of faces provides the cutoff point for subdivision
        if( primitiveIndices.size() <= 1 ) {
            outNode.initialize( primitiveIndices );
            return;
        }

        // The 3 arrays of split candidates for X, Y, and Z.
        // The face count which actually intersected the bound box could conceivably be less than
        // primitiveIndices.count(), so we get that value here as well.
        std::vector<split_candidate> splitCandidates[3];
        int totalFaceCount = get_split_candidates( m_primitiveBounds, primitiveIndices, bounds, splitCandidates );

        /*
        if( depth == 10 ) {
          for( int axis = 0; axis < 3; ++axis ) {
            cout << "Axis: " << axis << endl;
            for( unsigned i = 0; i < splitCandidates[axis].size(); ++i ) {
              splitCandidates[axis][i].print(cout);
            }
          }
          exit(1);
        }
        //*/

        // The best split value found.
        float bestSplitValue = 0;
        // The axis of the best split value.
        int bestAxis = -1;
        // The cost of the best split value.
        float bestCost = ( std::numeric_limits<float>::max )();
        // Which side of the tree to put the perp faces on.
        int bestPerpSide = 0;

        // Now do the sweep through the split candidates to find the best sweep.
        for( int testSplitAxis = 0; testSplitAxis < 3; ++testSplitAxis ) {
#ifdef DEBUG_PRIMITIVE_KDTREE
            if( depth == 100 )
                cout << "Axis: " << testSplitAxis << endl;
#endif // DEBUG_PRIMITIVE_KDTREE
       //  Only do this if the dimensions of the bounding box along this axis are positive.
            if( bounds.maximum()[testSplitAxis] - bounds.minimum()[testSplitAxis] > 0 ) {
                std::vector<split_candidate>& candidates = splitCandidates[testSplitAxis];
                // Start with all the faces on the right
                int leftFaceCount = 0, perpFaceCount = 0, rightFaceCount = totalFaceCount;
                boundbox3f leftBounds = bounds, rightBounds = bounds;
                unsigned i = 0, candidatesSize = (unsigned)candidates.size();
                while( i < candidatesSize ) {
                    float splitValue = candidates[i].splitValue;
                    int endEventCount = 0, perpEventCount = 0, startEventCount = 0;
                    // Go through all the events coincident with this position.  This leaves i pointing at the next
                    // split candidate, or past the end.
                    while( i < candidatesSize && candidates[i].splitValue == splitValue &&
                           candidates[i].splitType == -1 ) {
                        ++endEventCount;
                        ++i;
                    }
                    while( i < candidatesSize && candidates[i].splitValue == splitValue &&
                           candidates[i].splitType == 0 ) {
                        ++perpEventCount;
                        ++i;
                    }
                    while( i < candidatesSize && candidates[i].splitValue == splitValue &&
                           candidates[i].splitType == +1 ) {
                        ++startEventCount;
                        ++i;
                    }
                    // Adjust all the face counts
                    perpFaceCount = perpEventCount;
                    rightFaceCount -= perpEventCount;
                    rightFaceCount -= endEventCount;

                    // If a split is at the edge of the box, and there are no faces that could go in this empty box,
                    // then don't consider its cost
                    // TODO:  This shouldn't be using a hard coded epsilon value for the float comparison
                    if( splitValue - bounds.minimum()[testSplitAxis] > 0.0001f &&
                        splitValue - bounds.maximum()[testSplitAxis] < -0.0001f ) {
                        // Adjust the bounding boxes for this split
                        leftBounds.maximum()[testSplitAxis] = splitValue;
                        rightBounds.minimum()[testSplitAxis] = splitValue;
                        // Compute the surface area heuristic for this split, and check whether it's better than our
                        // current best one.
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

#ifdef DEBUG_PRIMITIVE_KDTREE
                        if( depth == 100 )
                            cout << "split: " << splitValue << ", cost: " << cost << ", left: " << leftFaceCount
                                 << ", perp: " << perpFaceCount << ", right: " << rightFaceCount
                                 << ", side: " << perpSide << endl;
#endif // DEBUG_PRIMITIVE_KDTREE
                    }
                    leftFaceCount += startEventCount;
                    leftFaceCount += perpEventCount;
                    perpFaceCount = 0;
                }
            }
        }

#ifdef DEBUG_PRIMITIVE_KDTREE
        if( depth == 100 ) {
            cout << "best split: " << bestSplitValue << endl;
            cout << "cost: " << bestCost << endl;
            cout << "perp side: " << bestPerpSide << endl;
            cout << "axis: " << bestAxis << endl;
            cout << "face indices size: " << primitiveIndices.size() << endl;
            exit( 1 );
        }
#endif // DEBUG_PRIMITIVE_KDTREE

        //		cout << "bestCost " << bestCost << " costThreshold " << primitiveIndices.size() << endl;
        // If we didn't find an acceptable split, (i.e. just running through all the faces costs less), then return a
        // leaf node. The 1.0001f is in here because floating point error was causing the cost to be less than the size,
        // even when it should have been equal to it.
        if( bestCost * 1.0001f >= primitiveIndices.size() ) {
            outNode.initialize( primitiveIndices );
            return;
        }

        // Compute the bounding boxes of both sides
        vector3f newLeftMaximum = bounds.maximum();
        newLeftMaximum[bestAxis] = bestSplitValue;
        boundbox3f boundsLeft( bounds.minimum(), newLeftMaximum );

        vector3f newRightMinimum = bounds.minimum();
        newRightMinimum[bestAxis] = bestSplitValue;
        boundbox3f boundsRight( newRightMinimum, bounds.maximum() );

        std::vector<int> leftPrimitiveIndices, rightPrimitiveIndices;
        std::vector<split_candidate>& candidates = splitCandidates[bestAxis];
        unsigned candidatesSize = (unsigned)candidates.size();
        // Now go through the split candidate array, and collect the face indices
        unsigned i = 0;
        while( i < candidatesSize && candidates[i].splitValue < bestSplitValue ) {
            // Add all faces that had a start event or perp event in this range to the left face list
            if( candidates[i].splitType >= 0 )
                leftPrimitiveIndices.push_back( candidates[i].splitPrimitiveIndex );
            ++i;
        }
        while( i < candidatesSize && candidates[i].splitValue == bestSplitValue ) {
            // Add all the faces that had a perp event to the left or right event as determined by the SAHs
            if( candidates[i].splitType == 0 ) {
                if( bestPerpSide < 0 ) {
                    leftPrimitiveIndices.push_back( candidates[i].splitPrimitiveIndex );
                } else {
                    rightPrimitiveIndices.push_back( candidates[i].splitPrimitiveIndex );
                }
            }
            ++i;
        }
        while( i < candidatesSize ) {
            // Add all the faces that had an end event or perp event in this range to the right face list
            if( candidates[i].splitType <= 0 )
                rightPrimitiveIndices.push_back( candidates[i].splitPrimitiveIndex );
            ++i;
        }

        // Now free the memory for the candidate splits.
        for( int axis = 0; axis < 3; ++axis ) {
            std::vector<split_candidate> freeArray;
            splitCandidates[axis].swap( freeArray );
        }

        //		cout << "Splitting interior node at depth " << depth << " with bounds " << bounds << ", axis "
        //<< bestAxis <<
        //", split " << bestSplitValue << ", perpSide " << bestPerpSide << endl;

        // These recursive calls to construct_kdtree may destructively swap out the leftPrimitiveIndices and
        // rightPrimitiveIndices contents, so those vectors aren't valid after this call.
        //		cout << "creating left child at depth " << depth+1 << " with bounds " << boundsLeft << ", " <<
        // leftPrimitiveIndices.size() << " nodes " << endl;
        primitive_kdtree_node<PrimitiveType>* children = new primitive_kdtree_node<PrimitiveType>[2];
        construct_kdtree_n_log2_n( leftPrimitiveIndices, boundsLeft, depth + 1, children[0] );
        //		cout << "creating right child at depth " << depth+1 << " with bounds " << boundsRight << ", " <<
        // rightPrimitiveIndices.size() << " nodes " << endl;
        construct_kdtree_n_log2_n( rightPrimitiveIndices, boundsRight, depth + 1, children[1] );

        outNode.initialize( children, bestAxis, bestSplitValue );
    }
};

} // namespace geometry
} // namespace frantic
