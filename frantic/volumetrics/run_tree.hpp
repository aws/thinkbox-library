// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/shared_array.hpp>
#include <cstdlib>

namespace frantic {
namespace volumetrics {

class run_tree {

    // Node state enum.  A node has 3 states.
    enum {
        CHILDREN = 0x02,  // the node has children and is neither defined or undefined
        UNDEFINED = 0x00, // the node is part of an undefined run
        DEFINED = 0x01,   // the node is part of a defined run
    };

    boost::shared_array<unsigned char> m_nodes;
    int m_extentSize;    // the actual extent size encoded
    int m_maxExtentSize; // the next biggest power of 2 extent size

    /**
     * Split a run interval into a left and right interval around a split index.
     *
     * @param run			The run to be split.
     * @param splitIndex	The index to be split around.
     * @param leftRun		The left run return value.
     * @param rightRun		The right run return value.
     *
     */
    inline void split_run( const std::pair<int, int>& run, int splitIndex, std::pair<int, int>& leftRun,
                           std::pair<int, int>& rightRun ) {
        leftRun.first = std::min( run.first, splitIndex );
        leftRun.second = std::min( run.second, splitIndex - 1 );
        rightRun.first = std::max( splitIndex, run.first ) - splitIndex;
        rightRun.second = std::max( run.second, splitIndex - 1 ) - splitIndex;
    }

    /**
     * Recursively inserts a run into the run tree.  The run is just an interval along the
     * current extentSize. If the run is the same as the extent size the node should cover,
     * the node is updated as defined and recursion stops.  If the run is empty, recursion stops.
     * Otherwise, the run is split and recursive calls are made for the left and right halves.
     * If an insertion causes all nodes below a given node to become defined, the node definition
     * propagates back upwards.
     *
     * @param run	A standard pair containing the start and end index in the current extentsize
     *				of the run interval.
     * @param extentSize	The current extent size, depending on the level of recursion.
     * @param nodeIndex		The index of the current node, also depending on the level of recursion.
     * @param extentData	Pointers to the start of the extent in each data channel.
     * @param initData		Initialization data for each channel.  Each sub-vector of chars represents
     *						the data to initialize with, and are of the appropriate length for the
     *datasize.
     */
    inline void insert_run( std::pair<int, int> run, int extentSize, int nodeIndex, int dataOffset,
                            std::vector<char*>& extentData, const std::vector<std::vector<char>>& initData ) {

        if( m_nodes[nodeIndex] == DEFINED )
            return;

        // If the run size equals the run extent, set the node to defined and initialize the data if
        // the node doesn't have children.  If it does, I have to keep going because some of it
        // is initialized already.
        if( run.second - run.first == extentSize - 1 && m_nodes[nodeIndex] != CHILDREN ) {

            for( size_t j = 0; j < extentData.size(); ++j ) {
                char* p = extentData[j] + (int)initData[j].size() * dataOffset;
                for( int i = 0; i < extentSize; ++i ) {
                    memcpy( p, &( initData[j][0] ), initData[j].size() );
                    p += initData[j].size();
                }
            }

            m_nodes[nodeIndex] = DEFINED;
            return;
        }

        // if the run is empty, just return
        if( run.second - run.first < 0 )
            return;

        // otherwise, you've got a run that needs to be split and sent to the children
        std::pair<int, int> leftRun, rightRun;
        int newExtentSize = extentSize / 2;
        split_run( run, newExtentSize, leftRun, rightRun );

        int leftChild = 2 * nodeIndex + 1, rightChild = 2 * nodeIndex + 2;

        // If the node didn't have children previously, initialize them to undefined and recurse
        if( m_nodes[nodeIndex] == UNDEFINED ) {
            m_nodes[leftChild] = UNDEFINED;
            m_nodes[rightChild] = UNDEFINED;
        }
        m_nodes[nodeIndex] = CHILDREN;

        insert_run( leftRun, newExtentSize, leftChild, dataOffset, extentData, initData );

        // You don't have to touch the extent data pointer when you go left, but when you go right
        // you have to bump it over.
        dataOffset += newExtentSize;
        insert_run( rightRun, newExtentSize, rightChild, dataOffset, extentData, initData );

        // if recursive insertions come back with 2 defined children, this node needs to be updated
        // to defined.
        if( m_nodes[leftChild] == DEFINED && m_nodes[rightChild] == DEFINED ) {
            m_nodes[nodeIndex] = DEFINED;
        }

        return;
    }

    /**
     *	Recursively traverses tree nodes until it finds a defined/undefined node containing the
     *  index and returns appropriately.
     *
     *	@param	index	The index to determine whether it is in a defined run.
     */
    inline bool is_defined( int index, int nodeIndex, int splitIndex ) {

        switch( m_nodes[nodeIndex] ) {
        case DEFINED:
            return true;
            break;
        case UNDEFINED:
            return false;
            break;
        case CHILDREN:
            if( index < splitIndex )
                return is_defined( index, nodeIndex * 2 + 1, splitIndex / 2 );
            else
                return is_defined( index - splitIndex, nodeIndex * 2 + 2, splitIndex / 2 );
            break;
        default:
            throw std::runtime_error( "run_tree::is_defined() - Unknown node state (" +
                                      boost::lexical_cast<std::string>( m_nodes[nodeIndex] ) + ") for node " +
                                      boost::lexical_cast<std::string>( nodeIndex ) + "." );
        }
    }

    /**
     *	Dumps the node state at a given node index to a out file stream.
     *
     *	@param	nodeIndex	the node to dump
     *	@param	ofs			the stream to dump to
     */
    void dump_node( int nodeIndex, std::ostream& ofs ) {
        ofs << std::endl << "node " << nodeIndex << ":  " << std::endl;
        switch( m_nodes[nodeIndex] ) {
        case CHILDREN:
            ofs << "CHILDREN" << std::endl;
            break;
        case UNDEFINED:
            ofs << "NO CHILDREN - UNDEFINED" << std::endl;
            break;
        case DEFINED:
            ofs << "NO CHILDREN - DEFINED" << std::endl;
            break;
        default:
            ofs << "UNKNOWN NODE STATE" << std::endl;
            break;
        }
    }

  public:
    /**
     * Default constructor makes a tree covering extent of size 0.
     */
    run_tree() {
        m_nodes.reset();
        m_extentSize = 0;
        m_maxExtentSize = 0;
    }

    /**
     * Constructor takes the size of the extent to be run length encoded and allocates
     * enough storage for nodes up to the nearest power of 2.
     *
     *	@param  extentSize  the size of the extent to be run length encoded
     */
    run_tree( int extentSize ) {
        if( extentSize < 0 )
            throw std::runtime_error( "run_tree - Extent size to encode must be non negative." );

        m_extentSize = extentSize;

        // Check if the given extent size is a power of 2.  If it isn't
        // we need to pad it so that it is.
        int i = 1, nextPowerOf2 = 0, bitsFlagged = 0;
        for( int c = 0; c < sizeof( int ) * 8 - 1; c++ ) {
            if( m_extentSize & i ) {
                nextPowerOf2 = i << 1;
                bitsFlagged++;
            }
            i = i << 1;
        }

        if( bitsFlagged > 1 )
            m_maxExtentSize = nextPowerOf2; // not a power of 2, need the next biggest power of 2
        else
            m_maxExtentSize = m_extentSize;

        m_nodes.reset( new unsigned char[2 * m_maxExtentSize - 1] );
        m_nodes[0] = UNDEFINED;
    }

    /**
     *	Checks if the given index is within a defined run.
     *
     *	@param	index	The index to determine whether it is in a defined run.
     */
    inline bool is_defined( int index ) {

        if( index >= m_extentSize )
            throw std::runtime_error( "run_tree::is_defined() - Requested index " +
                                      boost::lexical_cast<std::string>( index ) + " is larger than the extent size (" +
                                      boost::lexical_cast<std::string>( m_extentSize ) + ")." );
        return is_defined( index, 0, m_maxExtentSize / 2 );
    }

    /**
     * Checks if the tree is empty
     *
     * @return bool		True if the tree has no defined runs, false otherwise.
     */
    inline bool is_empty() { return ( m_extentSize == 0 || m_nodes[0] == UNDEFINED ); }

    /**
     * Inserts a run into the tree.
     *
     * @param run	A standard pair of start/end indices of the run.
     */
    inline void insert_run( std::pair<int, int> run, std::vector<char*>& extentData,
                            const std::vector<std::vector<char>>& initData ) {

        if( run.second < run.first )
            throw std::runtime_error( "runtree::insert_run() - Attempted to insert a run with a starting index greater "
                                      "than the ending index: [" +
                                      boost::lexical_cast<std::string>( run.first ) + "," +
                                      boost::lexical_cast<std::string>( run.second ) + "]" );
        if( run.second - run.first > m_extentSize )
            throw std::runtime_error( "runtree::insert_run() - Attempted to insert run [" +
                                      boost::lexical_cast<std::string>( run.first ) + "," +
                                      boost::lexical_cast<std::string>( run.second ) + "] with length " +
                                      boost::lexical_cast<std::string>( run.second - run.first ) +
                                      ", which is greater than the tree's covered extent: " +
                                      boost::lexical_cast<std::string>( m_extentSize ) );
        // pad the run if it hits the end of the extent
        if( run.second == m_extentSize - 1 )
            run.second = m_maxExtentSize - 1;

        insert_run( run, m_maxExtentSize, 0, 0, extentData, initData );
    }

    /**
     * A non-recursive in-order tree traversal to output the run length encoded extent that the tree
     * covers.  The run length encoding here takes the form simply of a set of ordered indices,
     * indicating the starts of runs of defined/undefined extents of data.  The first index
     * indicates the start of the first defined run, and so the following indices indicate
     *
     * @param runs	A vector of indices that will return the run length encoded information.
     */
    inline void traverse_runs( std::vector<std::pair<int, int>>& runs ) {
        runs.clear();

        if( !m_nodes || m_nodes[0] == UNDEFINED )
            return;

        if( m_nodes[0] == DEFINED ) {
            runs.push_back( std::pair<int, int>( 0, m_extentSize - 1 ) );
            return;
        }

        // some state information required to keep track of the current node index and its run size
        unsigned int nodeIndex = 0; // Current node index

        int currentRunIndex = 0, currentExtentSize = m_maxExtentSize;
        bool runDefined = false, done = false;

        // go all the way down the left side as far as you can
        while( m_nodes[nodeIndex] == CHILDREN ) {
            nodeIndex = nodeIndex * 2 + 1;
            currentExtentSize /= 2;
        }

        // main traversal loop
        while( !done ) {
            // go to the right if you can
            if( m_nodes[nodeIndex] == CHILDREN ) {
                nodeIndex = nodeIndex * 2 + 2;

                // then go to the left as far as you can
                while( m_nodes[nodeIndex] == CHILDREN ) {
                    nodeIndex = nodeIndex * 2 + 1;
                    currentExtentSize /= 2;
                }
            } else {

                // run start
                if( m_nodes[nodeIndex] == DEFINED && !runDefined ) {
                    runDefined = !runDefined;
                    runs.push_back( std::pair<int, int>( currentRunIndex, -1 ) );
                }
                // run end
                else if( m_nodes[nodeIndex] == UNDEFINED && runDefined ) {
                    runDefined = !runDefined;
                    runs.back().second = currentRunIndex - 1;
                }

                currentRunIndex += currentExtentSize;

                // while i'm coming from a right node, keep going back up
                while( ( nodeIndex & 1 ) == 0 ) {
                    currentExtentSize *= 2;
                    nodeIndex = ( nodeIndex - 1 ) / 2;
                    if( nodeIndex == 0 )
                        done = true;
                }

                // i'm no longer going up to a parent from the right, so move up one more
                // from the left, so i can start going down the right again
                nodeIndex = ( nodeIndex - 1 ) / 2;
            }
        }

        // if we have an odd number of indices, that means we need to close the last defined run
        if( runs.back().second == -1 )
            runs.back().second = m_extentSize - 1;
    }
};
} // namespace volumetrics
} // namespace frantic
