// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/named_channel_data.hpp>

namespace frantic {
namespace channels {

namespace detail {
/**
 * The result of channel_op_node::compile(), stores a function ptr and some arbitrary
 * data allocated using malloc. The codeData member will be passed in to calls to
 * codePtr, as the third argument. The first and second arguments to codePtr are the
 * input channels, and temporary stack respectively.
 */
struct code_segment {
    void ( *codePtr )( char*, char*, void* );
    void* codeData;
};

/**
 * A helper function that automatically wraps a char* into a std::string
 * @param type The data_type_t to get a string representation of.
 * @return The std::string representation of the supplied type.
 */
inline std::string data_type_string( data_type_t type ) {
    return frantic::strings::to_string( frantic::channels::channel_data_type_str( type ) );
}
/**
 * @overload
 */
inline std::string data_type_string( int arity, data_type_t type ) {
    return frantic::strings::to_string( frantic::channels::channel_data_type_str( arity, type ) );
}
} // namespace detail

/**
 * An exception class wrapping any error occuring during the compilation of an AST node. It
 * stores the exception message, as well as the node responsible.
 */
class channel_compiler_error : public std::runtime_error {
    int m_nodeIndex;

  public:
    channel_compiler_error( int nodeIndex, const std::string& message )
        : runtime_error( message )
        , m_nodeIndex( nodeIndex ) {}

    int which_node() const throw() { return m_nodeIndex; }
};

/**
 * This object stores the type, arity and stack offset of temporaries used during execution
 * of a channel_operation_compiler::eval() call.
 */
struct temporary_result {
    int offset;
    int arity;
    data_type_t type;

    temporary_result() {}

    temporary_result( const temporary_result& rhs ) {
        offset = rhs.offset;
        arity = rhs.arity;
        type = rhs.type;
    }

    temporary_result( int _offset, int _arity, data_type_t _type )
        : offset( _offset )
        , arity( _arity )
        , type( _type ) {}
};

template <class T>
struct deleter {
    static void impl( void* ptr ) { delete(T*)ptr; }
};

/**
 * This class is responsible for holding the state of input and output channels, the temporary
 * results of each node, and the individual code segments produced by the AST.
 */
class channel_operation_compiler {
    // Stores output code segments, in the order that the should be evaluated.
    std::vector<detail::code_segment> codeSegs;

    // a set of objects and deleter functions to be deleted when the compiler is deleted.
    std::vector<std::pair<void*, void ( * )( void* )>> scopedObjs;

    // Stores the offset of the ouput of each node after it runs. This is used to prevent
    // recursive compilation by checking if this node has already produced output. This is
    // also used to retrieve the results of child AST nodes, or to determine if they have
    // not already been compiled.
    std::map<int, temporary_result> tempOffsets;

    // List of all possible input channels. New output channels should append themselves
    // to this channel map so that external consumers are aware of the results of this operation.
    frantic::channels::channel_map nativePcm;

    // List of all current input channels. Can be modified by appending extra channels if they
    // are required as inputs. Channels should only be appended if they exist in the nativePcm.
    frantic::channels::channel_map internalPcm;

    // Stores offset of next temporary in the stack. After compilation is complete, this is
    // the overall size of the stack required by this operation.
    int nextOutputOffset;

    // a flag set to avoid having to check if the index channel exists everytime eval is called
    // when it is set to true the compiler will assume it exists
    bool hasIndexChannel;

    // an accessor to the index channel
    // this avoids getting this accessor every time eval is run
    // this is set when hasIndexChannel is set to true
    frantic::channels::channel_cvt_accessor<int> indexAcc;

  public:
    /**
     * This default constructor is provided for container support only, do not use in your code.
     */
    channel_operation_compiler() {}

    channel_operation_compiler( const frantic::channels::channel_map& pcm,
                                const frantic::channels::channel_map& nativePcm );
    ~channel_operation_compiler();

    /**
     * Resets the compiler to the same state as if a new object was created.
     * @param pcm Same as in channel_operation_compiler( pcm, nativePcm )
     * @param nativePcm Same as in channel_operation_compiler( pcm, nativePcm )
     */
    void reset( const frantic::channels::channel_map& pcm, const frantic::channels::channel_map& nativePcm );

    frantic::channels::channel_map& get_channel_map();
    frantic::channels::channel_map& get_native_channel_map();

    const frantic::channels::channel_map& get_channel_map() const;
    const frantic::channels::channel_map& get_native_channel_map() const;

    /**
     * Sets the has index channel flag
     * @param hasChannel true if the channel exists
     */
    void set_has_index_channel( bool hasChannel );

    /**
     * Prints a debug representation of the stack to a stream.
     * @param stream The ostream to write to.
     */
    // void debug_print( std::ostream& stream ) const; //Disabled due to Unicode switching

    /**
     * Returns the type, arity and stack offset of the data produced by the node.
     *
     * @param nodeId The id of the node producing the data.
     * @return NULL if the node has not been processed yet, or a pointer to the structure
     *         describing the output of the specified AST node.
     */
    temporary_result* get_node_results( int nodeId );

    /**
     * Copies the result information from the node with id sourceNodeId to the node with id nodeId.
     * The end result is that the node w/ id nodeId has the same output info as the node w/ id sourceNodeId.
     * @note This is primarily used by pass through nodes that don't acutally generate their own results.
     *       It might be advisable to re-evaluate that architecture.
     *
     * @param destNodeId The node which will obtain the new copied output.
     * @param sourceNodeId The node which provides the output to be copied.
     */
    void copy_node_results( int destNodeId, int sourceNodeId );

    /**
     * Allocates stack space of a given data type and arity for the specified node. Returns
     * this information in the same format as get_node_results()
     *
     * @param nodeId The id of the node producing the data.
     * @param type The data type (float32, int32) of the node's output.
     * @param arity The arity of the node's output.
     * @result A pointer to a structure describing the results of the specified node.
     */
    temporary_result* allocate_temporary( int nodeId, data_type_t type, int arity );

    /**
     * Adds a new code segment that will be executed after all previous code segments.
     * The member 'segment.codeData' must be NULL or allocated with malloc().
     *
     * @param segment The code segment to append to the internal list.
     */
    void append_code_segment( const detail::code_segment& segment );

    /**
     * This function will evaluate the code segments produced by compiling an AST, using the
     * supplied pointer as input values laid out according to internalPcm. The output of the
     * channel operation will be stored in one of the particles channels.
     *
     * @param particle A pointer to a particle with memory layout as described by internalPcm.
     */
    void eval( char* particle, std::size_t index ) const;

    /**
     * This function will register an object allocated with operator new() tsuch that
     * the object will be deleted when this compiler object is destroyed.
     */
    template <class T>
    void register_scoped_object( T* obj ) {
        scopedObjs.push_back( std::pair<void*, void ( * )( void* )>( obj, &deleter<T>::impl ) );
    }
};

/**
 * Abstract base class for the AST nodes of the channel operation compiler. There should be
 * one node per operation, such as Addition, Type Casting, etc. Shared information among all
 * AST nodes goes here, for example the node ID (used for referencing the output of other nodes).
 */
class channel_op_node {
  protected:
    // The reference ID of this node.
    int m_nodeId;

  public:
    channel_op_node( int id )
        : m_nodeId( id ) {}

    virtual ~channel_op_node() {}

    /**
     * This is the workhorse of the AST. Calling this function will generate the code segment that this
     * node represents, recursively evaluating its input nodes. The results of this node's compilation
     * are passed into 'inoutCompData'.
     *
     * @param expressionTree A list of all other nodes in the overall expression. The index into this list
     *                       is also the node's ID.
     * @param inoutCompData Used for allocating space for results, retrieving the results of child nodes,
     *                      and storing the code segment output from this compilation step.
     */
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData ) = 0;
};

} // namespace channels
} // namespace frantic
