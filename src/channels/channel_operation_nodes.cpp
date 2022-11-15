// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

// TODO: This file is HUGE, at > 2700 lines, should be split up into multiple files, probably its own directory

#include <boost/mpl/assert.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/thread/tss.hpp>

#include <frantic/channels/channel_operation_nodes.hpp>
#include <frantic/graphics/quat4f.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/math/perlin_noise_generator.hpp>

#include <limits>

using frantic::graphics::transform4f;
using frantic::graphics::vector3f;

namespace std {
float sqrt( half _X ) { return sqrtf( _X ); }
} // namespace std

namespace {
using namespace frantic::channels;

typedef void ( *bind_template_fn_return_type )( char*, char*, void* );

/**
 * This template function will bind a binary template function T::eval to
 * T::eval<T1,T2> where T1 is the type corresponding to type1 and T2 is the type
 * corresponding to type2.
 */
template <class T>
inline bind_template_fn_return_type bind_template_fn( data_type_t type1, data_type_t type2 ) {
    switch( type1 ) {
    case data_type_int16:
        switch( type2 ) {
        case data_type_int8:
            return &T::template eval<boost::int16_t, boost::int8_t>;
        case data_type_int16:
            return &T::template eval<boost::int16_t, boost::int16_t>;
        case data_type_int32:
            return &T::template eval<boost::int16_t, boost::int32_t>;
        case data_type_int64:
            return &T::template eval<boost::int16_t, boost::int64_t>;
        default:
            return NULL;
        }
    case data_type_int32:
        switch( type2 ) {
        case data_type_int8:
            return &T::template eval<boost::int32_t, boost::int8_t>;
        case data_type_int16:
            return &T::template eval<boost::int32_t, boost::int16_t>;
        case data_type_int32:
            return &T::template eval<boost::int32_t, boost::int32_t>;
        case data_type_int64:
            return &T::template eval<boost::int32_t, boost::int64_t>;
        default:
            return NULL;
        }
    case data_type_float16:
        switch( type2 ) {
        case data_type_float16:
            return &T::template eval<half, half>;
        case data_type_float32:
            return &T::template eval<half, float>;
        case data_type_float64:
            return &T::template eval<half, double>;
        default:
            return NULL;
        }
    case data_type_float32:
        switch( type2 ) {
        case data_type_float16:
            return &T::template eval<float, half>;
        case data_type_float32:
            return &T::template eval<float, float>;
        case data_type_float64:
            return &T::template eval<float, double>;
        default:
            return NULL;
        }
    case data_type_float64:
        switch( type2 ) {
        case data_type_float16:
            return &T::template eval<double, half>;
        case data_type_float32:
            return &T::template eval<double, float>;
        case data_type_float64:
            return &T::template eval<double, double>;
        default:
            return NULL;
        }
    default:
        return NULL;
    }
}

template <class T, class TypeList, bool NotDone = ( boost::mpl::size<TypeList>::value > 0 )>
struct bind_template_fn_impl {};

template <class T, class TypeList>
struct bind_template_fn_impl<T, TypeList, true> {
    typedef void ( *result_type )( char*, char*, void* );

    static result_type eval( data_type_t type ) {
        typedef typename boost::mpl::front<TypeList>::type head_type;
        typedef typename boost::mpl::pop_front<TypeList>::type tail_types;

        if( frantic::channels::channel_data_type_traits<head_type>::data_type() == type )
            return &T::template eval<head_type>;

        return bind_template_fn_impl<T, tail_types>::eval( type );
    }
};

template <class T, class TypeList>
struct bind_template_fn_impl<T, TypeList, false> {
    static void ( *eval( data_type_t /*type*/ ) )( char*, char*, void* ) { return NULL; }
};

/**
 * @fn bind_template_fn
 * This function will choose a specialization of the static template function T::eval()
 * based on the runtime supplied data_type_t. The available specializations are specified
 * by the boost::mpl::vector of types TypeList.
 *
 * @tparam T The container class that has a static template function eval()
 * @tparam TypeList The list of acceptable types to bind T::eval() to
 * @param type The runtime determined type for which to retrieve a specific specialization
 *             of T::eval()
 * @return A function pointer with signature void(char*,char*,void*) that is a specific
 *         template specialization of the template function T::eval()
 */
template <class T, class TypeList>
inline bind_template_fn_return_type bind_template_fn( data_type_t type ) {
    return bind_template_fn_impl<T, TypeList>::eval( type );
}

/**
 * This function will check if a node has already been compiled, and if so retrieve its results
 * description. If it has not been compiled, it will do so and return its results description.
 *
 * @param id The id of the node to compile/get results from
 * @param expressionTree The collection of AST nodes that make up the expression
 * @param inoutCompData The compiler where the results are allocated and code segments collected.
 * @return The descriptor for the output of this the node with the specified ID
 */
inline temporary_result* get_or_create_node_result( int id, const std::vector<channel_op_node*>& expressionTree,
                                                    channel_operation_compiler& inoutCompData ) {
    temporary_result* result = inoutCompData.get_node_results( id );
    if( !result ) {
        expressionTree[id]->compile( expressionTree, inoutCompData );
        result = inoutCompData.get_node_results( id );
    }
    return result;
}

template <int NumInputs>
struct code_gen {
    int m_srcIndex[NumInputs];
    int m_destIndex;
    int m_arity;
    void init( int arity, int destIndex, int srcIndices[] ) {
        m_arity = arity;
        m_destIndex = destIndex;
        memcpy( m_srcIndex, srcIndices, sizeof( int ) * NumInputs );
    }
    void init( int arity, int destIndex, int src1Index ) {
        boost::mpl::assert<NumInputs == 1>();
        m_arity = arity;
        m_destIndex = destIndex;
        m_srcIndex[0] = src1Index;
    }
    void init( int arity, int destIndex, int src1Index, int src2Index ) {
        boost::mpl::assert<NumInputs == 2>();
        m_arity = arity;
        m_destIndex = destIndex;
        m_srcIndex[0] = src1Index;
        m_srcIndex[1] = src2Index;
    }
    void init( int arity, int destIndex, int src1Index, int src2Index, int src3Index ) {
        boost::mpl::assert<NumInputs == 3>();
        m_arity = arity;
        m_destIndex = destIndex;
        m_srcIndex[0] = src1Index;
        m_srcIndex[1] = src2Index;
        m_srcIndex[2] = src3Index;
    }
};

/**
 * This constraint terminates a list of constraints
 */
struct EmptyConstraint {
    static void check( int /*nodeId*/, temporary_result** /*pTemps*/ ) {}
};

/**
 * This constraint ensures the RHS input type matches the LHS input type
 */
template <int LHS, int RHS, class NextConstraint = EmptyConstraint>
struct TypeConstraintMatch {
    static void check( int nodeId, temporary_result** pTemps ) {
        if( pTemps[LHS]->type != pTemps[RHS]->type )
            throw channel_compiler_error( nodeId, "Unexpected Input " + boost::lexical_cast<std::string>( RHS + 1 ) +
                                                      " type: " + detail::data_type_string( pTemps[RHS]->type ) +
                                                      ", expected: " + detail::data_type_string( pTemps[LHS]->type ) );
        NextConstraint::check( nodeId, pTemps );
    }
};

/**
 * This constraint ensures the LHS input type matches the supplied constant
 */
template <int LHS, data_type_t Type, class NextConstraint = EmptyConstraint>
struct TypeConstraintEq {
    static void check( int nodeId, temporary_result** pTemps ) {
        if( pTemps[LHS]->type != Type )
            throw channel_compiler_error( nodeId, "Unexpected Input " + boost::lexical_cast<std::string>( LHS + 1 ) +
                                                      " type: " + detail::data_type_string( pTemps[LHS]->type ) +
                                                      ", expected: " + detail::data_type_string( Type ) );
        NextConstraint::check( nodeId, pTemps );
    }
};

/**
 * This constraint ensures the RHS input arity matches the LHS input arity
 */
template <int LHS, int RHS, class NextConstraint = EmptyConstraint>
struct ArityConstraintMatch {
    static void check( int nodeId, temporary_result** pTemps ) {
        if( pTemps[LHS]->arity != pTemps[RHS]->arity )
            throw channel_compiler_error( nodeId,
                                          "Unexpected Input " + boost::lexical_cast<std::string>( RHS + 1 ) +
                                              " arity: " + boost::lexical_cast<std::string>( pTemps[RHS]->arity ) +
                                              ", expected: " + boost::lexical_cast<std::string>( pTemps[LHS]->arity ) );
        NextConstraint::check( nodeId, pTemps );
    }
};

/**
 * This constraint ensures the LHS input arity matches the supplied constant
 */
template <int LHS, int Arity, class NextConstraint = EmptyConstraint>
struct ArityConstraintEq {
    static void check( int nodeId, temporary_result** pTemps ) {
        if( pTemps[LHS]->arity != Arity )
            throw channel_compiler_error( nodeId,
                                          "Unexpected Input " + boost::lexical_cast<std::string>( LHS + 1 ) +
                                              " arity: " + boost::lexical_cast<std::string>( pTemps[LHS]->arity ) +
                                              ", expected: " + boost::lexical_cast<std::string>( Arity ) );
        NextConstraint::check( nodeId, pTemps );
    }
};

/**
 * This constraint ensures the LHS input arity is not greater than the supplied constant
 */
template <int LHS, int Arity, class NextConstraint = EmptyConstraint>
struct ArityConstraintLessEq {
    static void check( int nodeId, temporary_result** pTemps ) {
        if( pTemps[LHS]->arity > Arity )
            throw channel_compiler_error(
                nodeId, "Unexpected Input " + boost::lexical_cast<std::string>( LHS + 1 ) +
                            " arity: " + boost::lexical_cast<std::string>( pTemps[LHS]->arity ) +
                            ", expected an arity less than or equal to: " + boost::lexical_cast<std::string>( Arity ) );
        NextConstraint::check( nodeId, pTemps );
    }
};

/**
 * This constraint ensures the LHS input arity is not less than the supplied constant
 */
template <int LHS, int Arity, class NextConstraint = EmptyConstraint>
struct ArityConstraintGreaterEq {
    static void check( int nodeId, temporary_result** pTemps ) {
        if( pTemps[LHS]->arity < Arity )
            throw channel_compiler_error(
                nodeId,
                "Unexpected Input " + boost::lexical_cast<std::string>( LHS + 1 ) +
                    " arity: " + boost::lexical_cast<std::string>( pTemps[LHS]->arity ) +
                    ", expected an arity greater than or equal to: " + boost::lexical_cast<std::string>( Arity ) );
        NextConstraint::check( nodeId, pTemps );
    }
};

/**
 * @class ArityConstraintMatchOrEq
 * This constraint ensures the RHS input arity matches the LHS input arity, or the supplied
 * constant.
 *
 * @tparam LHS Index of the input to check arity against
 * @tparam RHS Index of the input to check arity of. It must match input LHS's arity or be Arity
 * @tparam Arity A constant arity that input RHS must have if it doesn't match input LHS's arity
 * @tparam NextConstraint The next constraint to apply after this one.
 */
template <int LHS, int RHS, int Arity, class NextConstraint = EmptyConstraint>
struct ArityConstraintMatchOrEq {
    static void check( int nodeId, temporary_result** pTemps ) {
        if( pTemps[RHS]->arity != pTemps[LHS]->arity && pTemps[RHS]->arity != Arity )
            throw channel_compiler_error( nodeId,
                                          "Unexpected Input " + boost::lexical_cast<std::string>( LHS + 1 ) +
                                              " arity: " + boost::lexical_cast<std::string>( pTemps[RHS]->arity ) +
                                              ", expected: " + boost::lexical_cast<std::string>( pTemps[LHS]->arity ) +
                                              " or: " + boost::lexical_cast<std::string>( Arity ) );
        NextConstraint::check( nodeId, pTemps );
    }
};

/**
 * @fn generate_code_segment
 * This function encapsulates the general case of generating code segments. The caller supplies various
 * template traits class that define the behaviour of this function.
 *
 * @tparam Op This class must have several members.
 *             Op::eval must be a static function with signature void(char*,char*,void*)
 *             Op::init must be a member function with signature void(int, int, int[]) which initializes the operator's
 * malloc'd storage.
 *             Op::NumSources must be the number of inputs it requires
 *             Op::InputTypes must be a boost::mpl::sequence of types that are acceptable inputs
 *             Op::Constraints must be a linked-list of constraint objects (ex. TypeConstraintEq)
 *             Op::OutputArity must be an integer arity for the output, or 0 to take the first input's arity
 *             Op::OutputType must be a data_type_t for the output type, or data_type_invalid to use the first input's
 * type
 * @param nodeId The id of the node being compiled.
 * @param srcIndex The ids of the nodes connected to the current node.
 * @param expressionTree The list of nodes in the expression.
 * @param inoutCompData The compiler generating stack results and storing code segements.
 */
template <class Op>
detail::code_segment generate_code_segment( int nodeId, int srcIndex[],
                                            const std::vector<channel_op_node*>& expressionTree,
                                            channel_operation_compiler& inoutCompData ) {
    temporary_result* r = inoutCompData.get_node_results( nodeId );
    if( r != NULL )
        throw channel_compiler_error( nodeId, "Graph cycle detected" );

    temporary_result* src[Op::NumSources];
    for( int i = 0; i < Op::NumSources; ++i ) {
        if( srcIndex[i] < 0 || srcIndex[i] >= (int)expressionTree.size() )
            throw channel_compiler_error( nodeId, "Invalid connection id for connection #" +
                                                      boost::lexical_cast<std::string>( i + 1 ) );
        src[i] = get_or_create_node_result( srcIndex[i], expressionTree, inoutCompData );
    }
    Op::Constraints::check( nodeId, src );

    int outArity = ( Op::OutputArity > 0 ) ? Op::OutputArity : src[0]->arity;
    data_type_t outType = ( Op::OutputType != data_type_invalid ) ? Op::OutputType : src[0]->type;

    r = inoutCompData.allocate_temporary( nodeId, outType, outArity );

    int srcOffset[Op::NumSources];
    for( int i = 0; i < Op::NumSources; ++i )
        srcOffset[i] = src[i]->offset;

    detail::code_segment theSeg;
    theSeg.codePtr = bind_template_fn<Op, typename Op::InputTypes>( src[0]->type );
    theSeg.codeData = (Op*)malloc( sizeof( Op ) );
    reinterpret_cast<Op*>( theSeg.codeData )->init( src[0]->arity, r->offset, srcOffset );

    return theSeg;
}
} // namespace

namespace frantic {
namespace channels {

disabled_op_node::disabled_op_node( int id, int passThroughIndex )
    : channel_op_node( id )
    , m_passThroughIndex( passThroughIndex ) {}

disabled_op_node::~disabled_op_node() {}

void disabled_op_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                channel_operation_compiler& inoutCompData ) {
    if( m_passThroughIndex < 0 )
        throw channel_compiler_error( m_nodeId, "There is no input to pass through" );
    temporary_result* src1 = inoutCompData.get_node_results( m_passThroughIndex );
    if( !src1 )
        expressionTree[m_passThroughIndex]->compile( expressionTree, inoutCompData );
    inoutCompData.copy_node_results( m_nodeId, m_passThroughIndex );
}

output_channel_op_node::output_channel_op_node( int id, int inputIndex, const std::string& outChannel, int outArity,
                                                data_type_t outType )
    : channel_op_node( id )
    , m_outputChannel( outChannel )
    , m_outputType( outType )
    , m_outputArity( outArity )
    , m_inputIndex( inputIndex ) {}

output_channel_op_node::~output_channel_op_node() {}

const std::string& output_channel_op_node::get_channel_name() const { return m_outputChannel; }

struct output_node_channel_impl {
    frantic::channels::channel_general_accessor m_destAcc;
    frantic::channels::channel_general_accessor m_srcAcc;

    void init( const frantic::channels::channel_general_accessor& destAcc,
               const frantic::channels::channel_general_accessor& srcAcc ) {
        m_destAcc = destAcc;
        m_srcAcc = srcAcc;
    }

    template <class TDest, class TSrc>
    static void eval( char* particle, char*, void* pVoidData ) {
#if defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif
        output_node_channel_impl* thisData = reinterpret_cast<output_node_channel_impl*>( pVoidData );
        TSrc* pSrc = reinterpret_cast<TSrc*>( thisData->m_srcAcc.get_channel_data_pointer( particle ) );
        TDest* pDest = reinterpret_cast<TDest*>( thisData->m_destAcc.get_channel_data_pointer( particle ) );
        for( std::size_t i = 0; i < thisData->m_srcAcc.arity(); ++i )
            pDest[i] = TDest( pSrc[i] );
#if defined( _MSC_VER )
#pragma warning( pop )
#endif
    }
};

struct output_node_stack_impl {
    frantic::channels::channel_general_accessor m_destAcc;
    int m_srcOffset;

    void init( const frantic::channels::channel_general_accessor& destAcc, int srcOffset ) {
        m_destAcc = destAcc;
        m_srcOffset = srcOffset;
    }

    template <class TDest, class TSrc>
    static void eval( char* particle, char* stack, void* pVoidData ) {
#if defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif
        output_node_stack_impl* thisData = reinterpret_cast<output_node_stack_impl*>( pVoidData );
        TSrc* pSrc = reinterpret_cast<TSrc*>( stack + thisData->m_srcOffset );
        TDest* pDest = reinterpret_cast<TDest*>( thisData->m_destAcc.get_channel_data_pointer( particle ) );
        for( std::size_t i = 0; i < thisData->m_destAcc.arity(); ++i )
            pDest[i] = TDest( pSrc[i] );
#if defined( _MSC_VER )
#pragma warning( pop )
#endif
    }
};

void output_channel_op_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                      channel_operation_compiler& inoutCompData ) {
    using frantic::channels::is_channel_data_type_float;
    using frantic::channels::is_channel_data_type_int;
    using frantic::channels::is_channel_data_type_pod;

    if( m_inputIndex < 0 || m_inputIndex >= (int)expressionTree.size() )
        throw channel_compiler_error( m_nodeId, "Invalid connection id for connection #" +
                                                    boost::lexical_cast<std::string>( m_inputIndex + 1 ) );

    frantic::tstring outputChannel = frantic::strings::to_tstring( m_outputChannel );

    if( !inoutCompData.get_channel_map().has_channel( outputChannel ) )
        inoutCompData.get_channel_map().append_channel( outputChannel, m_outputArity, m_outputType );
    if( !inoutCompData.get_native_channel_map().has_channel( outputChannel ) )
        inoutCompData.get_native_channel_map().append_channel( outputChannel, m_outputArity, m_outputType );
    frantic::channels::channel_general_accessor outAcc =
        inoutCompData.get_channel_map().get_general_accessor( outputChannel );

    if( input_channel_op_node* inputNode = dynamic_cast<input_channel_op_node*>( expressionTree[m_inputIndex] ) ) {
        typedef output_node_channel_impl Op;

        frantic::tstring inputChannel = frantic::strings::to_tstring( inputNode->get_channel() );

        // the index channel is a special case since it will be filled when eval is called on the compiler
        if( inputChannel == _T("Index") ) {
            if( !inoutCompData.get_channel_map().has_channel( _T("Index") ) )
                inoutCompData.get_channel_map().append_channel( _T("Index"), 1, data_type_int32 );
            inoutCompData.set_has_index_channel( 1 );
        }

        if( !inoutCompData.get_channel_map().has_channel( inputChannel ) ) {
            if( !inoutCompData.get_native_channel_map().has_channel( inputChannel ) )
                throw channel_compiler_error( m_nodeId, "The channel: \"" + inputNode->get_channel() +
                                                            "\" was not available in this stream" );
            const frantic::channels::channel& ch = inoutCompData.get_native_channel_map()[inputChannel];
            inoutCompData.get_channel_map().append_channel( inputChannel, ch.arity(), ch.data_type() );
        }

        frantic::channels::channel_general_accessor inAcc =
            inoutCompData.get_channel_map().get_general_accessor( inputChannel );

        if( outAcc.arity() != inAcc.arity() )
            throw channel_compiler_error(
                m_nodeId, "Unexpected Input1 arity: " + boost::lexical_cast<std::string>( inAcc.arity() ) +
                              ", expected: " + boost::lexical_cast<std::string>( outAcc.arity() ) );

        if( is_channel_data_type_float( outAcc.data_type() ) && !is_channel_data_type_float( inAcc.data_type() ) )
            throw channel_compiler_error( m_nodeId,
                                          "Unexpected Input1 type: " + detail::data_type_string( inAcc.data_type() ) +
                                              ", expected a floating point type" );
        else if( is_channel_data_type_int( outAcc.data_type() ) && !is_channel_data_type_int( inAcc.data_type() ) )
            throw channel_compiler_error( m_nodeId,
                                          "Unexpected Input1 type: " + detail::data_type_string( inAcc.data_type() ) +
                                              ", expected an integer type" );

        detail::code_segment theSeg;
        theSeg.codePtr = bind_template_fn<Op>( outAcc.data_type(), inAcc.data_type() );
        theSeg.codeData = malloc( sizeof( Op ) );
        reinterpret_cast<Op*>( theSeg.codeData )->init( outAcc, inAcc );

        inoutCompData.append_code_segment( theSeg );
    } else {
        typedef output_node_stack_impl Op;

        temporary_result* src1 = inoutCompData.get_node_results( m_inputIndex );
        if( !src1 ) {
            expressionTree[m_inputIndex]->compile( expressionTree, inoutCompData );
            src1 = inoutCompData.get_node_results( m_inputIndex );
        }

        if( src1->arity != (int)outAcc.arity() )
            throw channel_compiler_error(
                m_nodeId, "Unexpected Input1 arity: " + boost::lexical_cast<std::string>( src1->arity ) +
                              ", expected: " + boost::lexical_cast<std::string>( outAcc.arity() ) );

        if( is_channel_data_type_float( outAcc.data_type() ) && src1->type != data_type_float32 )
            throw channel_compiler_error(
                m_nodeId, "Unexpected Input1 type: " + detail::data_type_string( src1->type ) + ", expected: float32" );
        else if( is_channel_data_type_int( outAcc.data_type() ) && src1->type != data_type_int32 )
            throw channel_compiler_error(
                m_nodeId, "Unexpected Input1 type: " + detail::data_type_string( src1->type ) + ", expected: int32" );

        detail::code_segment theSeg;
        theSeg.codePtr = bind_template_fn<Op>( outAcc.data_type(), src1->type );
        theSeg.codeData = malloc( sizeof( Op ) );
        reinterpret_cast<Op*>( theSeg.codeData )->init( outAcc, src1->offset );

        inoutCompData.append_code_segment( theSeg );
    }
}

input_channel_op_node::input_channel_op_node( int id, const std::string& inputChannel )
    : channel_op_node( id )
    , m_inputChannel( inputChannel ) {}

input_channel_op_node::~input_channel_op_node() {}

struct input_channel_node_impl {
    frantic::channels::channel_general_accessor m_srcAcc;
    int m_destOffset;
    void init( int destOffset, const frantic::channels::channel_general_accessor& srcAcc ) {
        m_destOffset = destOffset;
        m_srcAcc = srcAcc;
    }
    template <class TDest, class TSrc>
    static void eval( char* particle, char* stack, void* data ) {
#if defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif
        input_channel_node_impl* thisData = static_cast<input_channel_node_impl*>( data );
        TDest* pDest = reinterpret_cast<TDest*>( stack + thisData->m_destOffset );
        TSrc* pSrc = reinterpret_cast<TSrc*>( thisData->m_srcAcc.get_channel_data_pointer( particle ) );
        for( std::size_t i = 0; i < thisData->m_srcAcc.arity(); ++i )
            pDest[i] = TDest( pSrc[i] );
#if defined( _MSC_VER )
#pragma warning( pop )
#endif
    }
};

void input_channel_op_node::compile( const std::vector<channel_op_node*>& /*expressionTree*/,
                                     channel_operation_compiler& inoutCompData ) {
    typedef input_channel_node_impl Op;

    // the index channel is a special case since it will be filled when eval is called on the compiler
    if( m_inputChannel == "Index" ) {
        if( !inoutCompData.get_channel_map().has_channel( _T("Index") ) )
            inoutCompData.get_channel_map().append_channel( _T("Index"), 1, data_type_int32 );
        inoutCompData.set_has_index_channel( 1 );
    }

    frantic::tstring inputChannel = frantic::strings::to_tstring( m_inputChannel );

    if( !inoutCompData.get_channel_map().has_channel( inputChannel ) ) {
        if( !inoutCompData.get_native_channel_map().has_channel( inputChannel ) )
            throw channel_compiler_error( m_nodeId,
                                          "The channel: \"" + m_inputChannel + "\" was not available in this stream" );
        frantic::channels::channel_general_accessor nativeAcc =
            inoutCompData.get_native_channel_map().get_general_accessor( inputChannel );
        inoutCompData.get_channel_map().append_channel( inputChannel, nativeAcc.arity(), nativeAcc.data_type() );
    }

    frantic::channels::channel_general_accessor inAcc =
        inoutCompData.get_channel_map().get_general_accessor( inputChannel );

    data_type_t outType;
    if( frantic::channels::is_channel_data_type_float( inAcc.data_type() ) )
        outType = data_type_float32;
    else if( frantic::channels::is_channel_data_type_int( inAcc.data_type() ) )
        outType = data_type_int32;
    else
        throw channel_compiler_error( m_nodeId,
                                      "Unexpected Input1 type: " + detail::data_type_string( inAcc.data_type() ) );

    temporary_result* r = inoutCompData.get_node_results( m_nodeId );
    if( r != NULL )
        throw channel_compiler_error( m_nodeId, "Graph cycle detected" );
    r = inoutCompData.allocate_temporary( m_nodeId, outType, (int)inAcc.arity() );

    detail::code_segment theSeg;
    theSeg.codePtr = bind_template_fn<Op>( outType, inAcc.data_type() );
    theSeg.codeData = malloc( sizeof( Op ) );
    reinterpret_cast<Op*>( theSeg.codeData )->init( r->offset, inAcc );

    inoutCompData.append_code_segment( theSeg );
}

template <class DestType, unsigned OutArity>
convert_node<DestType, OutArity>::convert_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

template <class DestType, unsigned OutArity>
convert_node<DestType, OutArity>::~convert_node() {}

struct convert_node_impl {
    int m_destOffset;
    int m_srcOffset;
    int m_arity;

    void init( int arity, int destOffset, int srcOffset ) {
        m_destOffset = destOffset;
        m_srcOffset = srcOffset;
        m_arity = arity;
    }

    template <class TDest, class TSrc>
    static void eval( char*, char* stack, void* pVoidData ) {
#if defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif
        convert_node_impl* pData = reinterpret_cast<convert_node_impl*>( pVoidData );
        TDest* pDest = reinterpret_cast<TDest*>( stack + pData->m_destOffset );
        TSrc* pSrc = reinterpret_cast<TSrc*>( stack + pData->m_srcOffset );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = TDest( pSrc[i] );
#if defined( _MSC_VER )
#pragma warning( pop )
#endif
    }
};

template <class DestType, unsigned OutArity>
void convert_node<DestType, OutArity>::compile( const std::vector<channel_op_node*>& expressionTree,
                                                channel_operation_compiler& inoutCompData ) {
    // Prevent this node from being compiled twice.
    temporary_result* r = inoutCompData.get_node_results( m_nodeId );
    if( r != NULL )
        throw channel_compiler_error( m_nodeId, "Graph cycle detected" );

    const int NumSources = 1;

    temporary_result* src[NumSources];
    for( int i = 0; i < NumSources; ++i ) {
        if( m_srcIndex[i] < 0 || m_srcIndex[i] >= (int)expressionTree.size() )
            throw channel_compiler_error( m_nodeId, "Invalid connection id for connection #" +
                                                        boost::lexical_cast<std::string>( i + 1 ) );
        src[i] = get_or_create_node_result( m_srcIndex[i], expressionTree, inoutCompData );
    }

    data_type_t outType = channel_data_type_traits<DestType>::data_type();

#if defined( _MSC_VER )
#pragma warning( push, 3 )
#pragma warning( disable : 4127 )
#endif
    if( OutArity > 0 && src[0]->arity != OutArity )
        throw channel_compiler_error( m_nodeId,
                                      "Unexpected Input1 arity: " + boost::lexical_cast<std::string>( src[0]->arity ) +
                                          ", expected: " + boost::lexical_cast<std::string>( OutArity ) );
#if defined( _MSC_VER )
#pragma warning( pop )
#endif

    // Don't do anything if this is an identity conversion.
    if( outType == src[0]->type ) {
        inoutCompData.copy_node_results( m_nodeId, m_srcIndex[0] );
        return;
    }

    detail::code_segment theSeg;
    theSeg.codePtr = NULL;

    switch( outType ) {
    case data_type_float32:
        if( src[0]->type == data_type_int32 )
            theSeg.codePtr = &convert_node_impl::eval<float, int>;
        break;
    case data_type_int32:
        if( src[0]->type == data_type_float32 )
            theSeg.codePtr = &convert_node_impl::eval<int, float>;
        break;
    default:
        throw channel_compiler_error( m_nodeId,
                                      "Unexpected conversion target type: " + detail::data_type_string( outType ) );
    }

    if( !theSeg.codePtr )
        throw channel_compiler_error( m_nodeId, "Unexpected conversion source type: " +
                                                    detail::data_type_string( src[0]->type ) );

    r = inoutCompData.allocate_temporary( m_nodeId, outType, OutArity );

    theSeg.codeData = malloc( sizeof( convert_node_impl ) );
    reinterpret_cast<convert_node_impl*>( theSeg.codeData )->init( src[0]->arity, r->offset, src[0]->offset );

    inoutCompData.append_code_segment( theSeg );
}

// Instantiate conversions for floats and integer scalars.
template class convert_node<float, 1>;
template class convert_node<int, 1>;

negate_node::negate_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

negate_node::~negate_node() {}

struct negate_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_invalid;
    typedef EmptyConstraint Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        negate_node_impl* pData = reinterpret_cast<negate_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = -pSrc[i];
    }
};

void negate_node::compile( const std::vector<channel_op_node*>& expressionTree,
                           channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<negate_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

component_sum_node::component_sum_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

component_sum_node::~component_sum_node() {}

struct component_sum_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_invalid;
    typedef EmptyConstraint Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        component_sum_node_impl* pData = reinterpret_cast<component_sum_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        pDest[0] = pSrc[0];
        for( int i = 1; i < pData->m_arity; ++i )
            pDest[0] += pSrc[i];
    }
};

void component_sum_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                  channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<component_sum_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

transform_by_quat_node::transform_by_quat_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

transform_by_quat_node::~transform_by_quat_node() {}

struct transform_by_quat_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 3;
    static const data_type_t OutputType = data_type_float32;
    typedef ArityConstraintEq<0, 3, ArityConstraintEq<1, 4>> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        using frantic::graphics::quat4f;
        using frantic::graphics::transform4f;
        using frantic::graphics::vector3f;

        transform_by_quat_node_impl* pData = reinterpret_cast<transform_by_quat_node_impl*>( pVoidData );
        vector3f* pDest = reinterpret_cast<vector3f*>( stack + pData->m_destIndex );
        vector3f* pSrc1 = reinterpret_cast<vector3f*>( stack + pData->m_srcIndex[0] );
        quat4f* pSrc2 = reinterpret_cast<quat4f*>( stack + pData->m_srcIndex[1] );

        transform4f tm;
        pSrc2->as_transform4f( tm );

        *pDest = tm * ( *pSrc1 );
    }
};

struct transform_quat_by_quat_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 4;
    static const data_type_t OutputType = data_type_float32;
    typedef ArityConstraintEq<0, 4, ArityConstraintEq<1, 4>> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        using frantic::graphics::quat4f;
        using frantic::graphics::transform4f;
        using frantic::graphics::vector3f;

        transform_by_quat_node_impl* pData = reinterpret_cast<transform_by_quat_node_impl*>( pVoidData );
        quat4f* pDest = reinterpret_cast<quat4f*>( stack + pData->m_destIndex );
        quat4f* pSrc1 = reinterpret_cast<quat4f*>( stack + pData->m_srcIndex[0] );
        quat4f* pSrc2 = reinterpret_cast<quat4f*>( stack + pData->m_srcIndex[1] );

        *pDest = ( *pSrc1 ) * ( *pSrc2 );
    }
};

void transform_by_quat_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                      channel_operation_compiler& inoutCompData ) {
    // We must check for graph cycles before attempting to evaluate any inputs.
    temporary_result* r = inoutCompData.get_node_results( m_nodeId );
    if( r != NULL )
        throw channel_compiler_error( m_nodeId, "Graph cycle detected" );

    // We determine which implementation to use based on the arity of the first input.
    if( m_srcIndex[0] < 0 || (size_t)m_srcIndex[0] >= expressionTree.size() )
        throw channel_compiler_error( m_nodeId, "Connection #1 was invalid. Got: " +
                                                    boost::lexical_cast<std::string>( m_srcIndex[0] ) );

    temporary_result* src = get_or_create_node_result( m_srcIndex[0], expressionTree, inoutCompData );

    if( src->arity == 3 ) {
        inoutCompData.append_code_segment(
            generate_code_segment<transform_by_quat_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
    } else if( src->arity == 4 ) {
        inoutCompData.append_code_segment( generate_code_segment<transform_quat_by_quat_node_impl>(
            m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
    } else
        throw channel_compiler_error( m_nodeId,
                                      "Unexpected input1 arity: " + boost::lexical_cast<std::string>( src->arity ) +
                                          ", expected 3 (Vector) or 4 (Quaternion)" );
}

magnitude_node::magnitude_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

magnitude_node::~magnitude_node() {}

struct magnitude_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        magnitude_node_impl* pData = reinterpret_cast<magnitude_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        pDest[0] = pSrc[0] * pSrc[0];
        for( int i = 1; i < pData->m_arity; ++i )
            pDest[0] += pSrc[i] * pSrc[i];
        pDest[0] = std::sqrt( pDest[0] );
    }
};

void magnitude_node::compile( const std::vector<channel_op_node*>& expressionTree,
                              channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<magnitude_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

normalize_node::normalize_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

normalize_node::~normalize_node() {}

struct normalize_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        normalize_node_impl* pData = reinterpret_cast<normalize_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T len = pSrc[0] * pSrc[0];
        for( int i = 1; i < pData->m_arity; ++i )
            len += pSrc[i] * pSrc[i];
        len = std::sqrt( len );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = pSrc[i] / len;
    }
};

void normalize_node::compile( const std::vector<channel_op_node*>& expressionTree,
                              channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<normalize_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

abs_node::abs_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

abs_node::~abs_node() {}

struct abs_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_invalid;
    typedef EmptyConstraint Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        abs_node_impl* pData = reinterpret_cast<abs_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::abs( pSrc[i] );
    }
};

void abs_node::compile( const std::vector<channel_op_node*>& expressionTree,
                        channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<abs_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

floor_node::floor_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

floor_node::~floor_node() {}

struct floor_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        floor_node_impl* pData = reinterpret_cast<floor_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::floor( pSrc[i] );
    }
};

void floor_node::compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<floor_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

ceil_node::ceil_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

ceil_node::~ceil_node() {}

struct ceil_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        ceil_node_impl* pData = reinterpret_cast<ceil_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::ceil( pSrc[i] );
    }
};

void ceil_node::compile( const std::vector<channel_op_node*>& expressionTree,
                         channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<ceil_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

sin_trig_node::sin_trig_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

sin_trig_node::~sin_trig_node() {}

struct sin_trig_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        sin_trig_node_impl* pData = reinterpret_cast<sin_trig_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::sin( pSrc[i] );
    }
};

void sin_trig_node::compile( const std::vector<channel_op_node*>& expressionTree,
                             channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<sin_trig_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

cos_trig_node::cos_trig_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

cos_trig_node::~cos_trig_node() {}

struct cos_trig_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        cos_trig_node_impl* pData = reinterpret_cast<cos_trig_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::cos( pSrc[i] );
    }
};

void cos_trig_node::compile( const std::vector<channel_op_node*>& expressionTree,
                             channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<cos_trig_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

tan_trig_node::tan_trig_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

tan_trig_node::~tan_trig_node() {}

struct tan_trig_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        tan_trig_node_impl* pData = reinterpret_cast<tan_trig_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::tan( pSrc[i] );
    }
};

void tan_trig_node::compile( const std::vector<channel_op_node*>& expressionTree,
                             channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<tan_trig_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

asin_trig_node::asin_trig_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

asin_trig_node::~asin_trig_node() {}

struct asin_trig_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        asin_trig_node_impl* pData = reinterpret_cast<asin_trig_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::asin( pSrc[i] );
    }
};

void asin_trig_node::compile( const std::vector<channel_op_node*>& expressionTree,
                              channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<asin_trig_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

acos_trig_node::acos_trig_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

acos_trig_node::~acos_trig_node() {}

struct acos_trig_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        acos_trig_node_impl* pData = reinterpret_cast<acos_trig_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::acos( pSrc[i] );
    }
};

void acos_trig_node::compile( const std::vector<channel_op_node*>& expressionTree,
                              channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<acos_trig_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

atan_trig_node::atan_trig_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

atan_trig_node::~atan_trig_node() {}

struct atan_trig_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        atan_trig_node_impl* pData = reinterpret_cast<atan_trig_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::atan( pSrc[i] );
    }
};

void atan_trig_node::compile( const std::vector<channel_op_node*>& expressionTree,
                              channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<atan_trig_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

atan2_trig_node::atan2_trig_node( int id, int srcIndex1, int srcIndex2 )
    : operation_node<2>( id, srcIndex1, srcIndex2 ) {}

atan2_trig_node::~atan2_trig_node() {}

struct atan2_trig_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32, TypeConstraintEq<1, data_type_float32>> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        atan2_trig_node_impl* pData = reinterpret_cast<atan2_trig_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = atan2( pSrc1[i], pSrc2[i] );
    }
};

void atan2_trig_node::compile( const std::vector<channel_op_node*>& expressionTree,
                               channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<atan2_trig_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

transform_point_node::transform_point_node( int id, int src1Index, const transform4f& tm )
    : operation_node<1>( id, src1Index )
    , m_transform( tm ) {}

transform_point_node::~transform_point_node() {}

struct transform_point_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 3;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32, ArityConstraintEq<0, 3>> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    transform4f m_transform;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        transform_point_node_impl* pData = reinterpret_cast<transform_point_node_impl*>( pVoidData );
        vector3f* pDest = reinterpret_cast<vector3f*>( stack + pData->m_destIndex );
        vector3f* pSrc = reinterpret_cast<vector3f*>( stack + pData->m_srcIndex[0] );
        *pDest = pData->m_transform * ( *pSrc );
    }
};

void transform_point_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                    channel_operation_compiler& inoutCompData ) {
    typedef transform_point_node_impl Op;

    detail::code_segment theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );
    reinterpret_cast<Op*>( theSeg.codeData )->m_transform = m_transform;

    inoutCompData.append_code_segment( theSeg );
}

transform_vector_node::transform_vector_node( int id, int src1Index, const transform4f& tm )
    : operation_node<1>( id, src1Index )
    , m_transform( tm ) {}

transform_vector_node::~transform_vector_node() {}

struct transform_vector_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 3;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32, ArityConstraintEq<0, 3>> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    transform4f m_transform;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        transform_vector_node_impl* pData = reinterpret_cast<transform_vector_node_impl*>( pVoidData );
        vector3f* pDest = reinterpret_cast<vector3f*>( stack + pData->m_destIndex );
        vector3f* pSrc = reinterpret_cast<vector3f*>( stack + pData->m_srcIndex[0] );
        *pDest = pData->m_transform.transform_no_translation( *pSrc );
    }
};

void transform_vector_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                     channel_operation_compiler& inoutCompData ) {
    typedef transform_vector_node_impl Op;

    detail::code_segment theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );
    reinterpret_cast<Op*>( theSeg.codeData )->m_transform = m_transform;

    inoutCompData.append_code_segment( theSeg );
}

transform_normal_node::transform_normal_node( int id, int src1Index, const transform4f& tm )
    : operation_node<1>( id, src1Index )
    , m_transform( tm ) {}

transform_normal_node::~transform_normal_node() {}

struct transform_normal_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 3;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32, ArityConstraintEq<0, 3>> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    transform4f m_transform;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        transform_normal_node_impl* pData = reinterpret_cast<transform_normal_node_impl*>( pVoidData );
        vector3f* pDest = reinterpret_cast<vector3f*>( stack + pData->m_destIndex );
        vector3f* pSrc = reinterpret_cast<vector3f*>( stack + pData->m_srcIndex[0] );
        *pDest = pData->m_transform.transpose_transform_no_translation( *pSrc );
    }
};

void transform_normal_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                     channel_operation_compiler& inoutCompData ) {
    typedef transform_normal_node_impl Op;

    detail::code_segment theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );
    reinterpret_cast<Op*>( theSeg.codeData )->m_transform = m_transform.to_inverse();

    inoutCompData.append_code_segment( theSeg );
}

static boost::thread_specific_ptr<frantic::math::perlin_noise_generator> g_pNoiseGen;

noise_node::noise_node( int id, int src1Index, int src2Index, int src3Index )
    : operation_node<3>( id, src1Index, src2Index, src3Index ) {}

noise_node::~noise_node() {}

struct noise_node_impl : public code_gen<3> {
    static const int NumSources = 3;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<
        0, data_type_float32,
        ArityConstraintLessEq<
            0, 3,
            TypeConstraintEq<1, data_type_int32,
                             ArityConstraintEq<1, 1, TypeConstraintEq<2, data_type_float32, ArityConstraintEq<2, 1>>>>>>
        Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        noise_node_impl* pData = reinterpret_cast<noise_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        int* pSrc2 = reinterpret_cast<int*>( stack + pData->m_srcIndex[1] );
        float* pSrc3 = reinterpret_cast<float*>( stack + pData->m_srcIndex[2] );

        if( !g_pNoiseGen.get() )
            g_pNoiseGen.reset( new frantic::math::perlin_noise_generator );
        g_pNoiseGen->set_octaves( *pSrc2 );
        g_pNoiseGen->set_persistence( *pSrc3 );

        switch( pData->m_arity ) {
        case 1:
            pDest[0] = g_pNoiseGen->noise( pSrc1[0] );
            break;
        case 2:
            pDest[0] = g_pNoiseGen->noise( pSrc1[0], pSrc1[1] );
            break;
        case 3:
            pDest[0] = g_pNoiseGen->noise( pSrc1[0], pSrc1[1], pSrc1[2] );
            break;
        default:
#ifdef _WIN32
            __assume( 0 );
#else
            break;
#endif
        };
    }
};

void noise_node::compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<noise_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

dnoise_node::dnoise_node( int id, int src1Index, int src2Index, int src3Index )
    : operation_node<3>( id, src1Index, src2Index, src3Index ) {}

dnoise_node::~dnoise_node() {}

struct dnoise_node_impl : public code_gen<3> {
    static const int NumSources = 3;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<
        0, data_type_float32,
        ArityConstraintLessEq<
            0, 3,
            TypeConstraintEq<1, data_type_int32,
                             ArityConstraintEq<1, 1, TypeConstraintEq<2, data_type_float32, ArityConstraintEq<2, 1>>>>>>
        Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        dnoise_node_impl* pData = reinterpret_cast<dnoise_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        int* pSrc2 = reinterpret_cast<int*>( stack + pData->m_srcIndex[1] );
        float* pSrc3 = reinterpret_cast<float*>( stack + pData->m_srcIndex[2] );

        if( !g_pNoiseGen.get() )
            g_pNoiseGen.reset( new frantic::math::perlin_noise_generator );
        g_pNoiseGen->set_octaves( *pSrc2 );
        g_pNoiseGen->set_persistence( *pSrc3 );

        switch( pData->m_arity ) {
        case 1:
            g_pNoiseGen->dnoise( pSrc1[0], pDest );
            break;
        case 2:
            g_pNoiseGen->dnoise( pSrc1[0], pSrc1[1], pDest, &( pDest[1] ) );
            break;
        case 3:
            g_pNoiseGen->dnoise( pSrc1[0], pSrc1[1], pSrc1[2], pDest, &( pDest[1] ), &( pDest[2] ) );
            break;
        default:
#ifdef _WIN32
            __assume( 0 );
#else
            break;
#endif
        }
    }
};

void dnoise_node::compile( const std::vector<channel_op_node*>& expressionTree,
                           channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<dnoise_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}
less_comparison_node::less_comparison_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

less_comparison_node::~less_comparison_node() {}

struct less_comparison_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef TypeConstraintMatch<0, 1, ArityConstraintEq<0, 1, ArityConstraintEq<1, 1>>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        less_comparison_node_impl* pData = reinterpret_cast<less_comparison_node_impl*>( pVoidData );
        int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        pDest[0] = ( pSrc1[0] < pSrc2[0] ) ? 1 : 0;
    }
};

void less_comparison_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                    channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<less_comparison_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

less_or_equal_comparison_node::less_or_equal_comparison_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

less_or_equal_comparison_node::~less_or_equal_comparison_node() {}

struct less_or_equal_comparison_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef TypeConstraintMatch<0, 1, ArityConstraintEq<0, 1, ArityConstraintEq<1, 1>>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        less_or_equal_comparison_node_impl* pData = reinterpret_cast<less_or_equal_comparison_node_impl*>( pVoidData );
        int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        pDest[0] = ( pSrc1[0] <= pSrc2[0] ) ? 1 : 0;
    }
};

void less_or_equal_comparison_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                             channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment( generate_code_segment<less_or_equal_comparison_node_impl>(
        m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

greater_comparison_node::greater_comparison_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

greater_comparison_node::~greater_comparison_node() {}

struct greater_comparison_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef TypeConstraintMatch<0, 1, ArityConstraintEq<0, 1, ArityConstraintEq<1, 1>>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        greater_comparison_node_impl* pData = reinterpret_cast<greater_comparison_node_impl*>( pVoidData );
        int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        pDest[0] = ( pSrc1[0] > pSrc2[0] ) ? 1 : 0;
    }
};

void greater_comparison_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                       channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<greater_comparison_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

greater_or_equal_comparison_node::greater_or_equal_comparison_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

greater_or_equal_comparison_node::~greater_or_equal_comparison_node() {}

struct greater_or_equal_comparison_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef TypeConstraintMatch<0, 1, ArityConstraintEq<0, 1, ArityConstraintEq<1, 1>>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        greater_or_equal_comparison_node_impl* pData =
            reinterpret_cast<greater_or_equal_comparison_node_impl*>( pVoidData );
        int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        pDest[0] = ( pSrc1[0] >= pSrc2[0] ) ? 1 : 0;
    }
};

void greater_or_equal_comparison_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                                channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment( generate_code_segment<greater_or_equal_comparison_node_impl>(
        m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

equal_comparison_node::equal_comparison_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

equal_comparison_node::~equal_comparison_node() {}

struct equal_comparison_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef TypeConstraintMatch<0, 1, ArityConstraintEq<0, 1, ArityConstraintEq<1, 1>>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        equal_comparison_node_impl* pData = reinterpret_cast<equal_comparison_node_impl*>( pVoidData );
        int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        pDest[0] = ( pSrc1[0] == pSrc2[0] ) ? 1 : 0;
    }
};

void equal_comparison_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                     channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<equal_comparison_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

logical_and_node::logical_and_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

logical_and_node::~logical_and_node() {}

struct logical_and_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef TypeConstraintEq<0, data_type_int32,
                             ArityConstraintEq<0, 1, TypeConstraintEq<1, data_type_int32, ArityConstraintEq<1, 1>>>>
        Constraints;
    typedef boost::mpl::vector<int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        logical_and_node_impl* pData = reinterpret_cast<logical_and_node_impl*>( pVoidData );
        int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
        int* pSrc1 = reinterpret_cast<int*>( stack + pData->m_srcIndex[0] );
        int* pSrc2 = reinterpret_cast<int*>( stack + pData->m_srcIndex[1] );
        pDest[0] = ( pSrc1[0] && pSrc2[0] ) ? 1 : 0;
    }
};

void logical_and_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<logical_and_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

logical_or_node::logical_or_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

logical_or_node::~logical_or_node() {}

struct logical_or_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef TypeConstraintEq<0, data_type_int32,
                             ArityConstraintEq<0, 1, TypeConstraintEq<1, data_type_int32, ArityConstraintEq<1, 1>>>>
        Constraints;
    typedef boost::mpl::vector<int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        logical_or_node_impl* pData = reinterpret_cast<logical_or_node_impl*>( pVoidData );
        int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
        int* pSrc1 = reinterpret_cast<int*>( stack + pData->m_srcIndex[0] );
        int* pSrc2 = reinterpret_cast<int*>( stack + pData->m_srcIndex[1] );
        pDest[0] = ( pSrc1[0] || pSrc2[0] ) ? 1 : 0;
    }
};

void logical_or_node::compile( const std::vector<channel_op_node*>& expressionTree,
                               channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<logical_or_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

logical_not_node::logical_not_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

logical_not_node::~logical_not_node() {}

struct logical_not_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef TypeConstraintEq<0, data_type_int32, ArityConstraintEq<0, 1>> Constraints;
    typedef boost::mpl::vector<int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        logical_not_node_impl* pData = reinterpret_cast<logical_not_node_impl*>( pVoidData );
        int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
        int* pSrc1 = reinterpret_cast<int*>( stack + pData->m_srcIndex[0] );
        pDest[0] = ( pSrc1[0] ) ? 0 : 1;
    }
};

void logical_not_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<logical_not_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

addition_node::addition_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

addition_node::~addition_node() {}

struct addition_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintMatch<0, 1, ArityConstraintMatch<0, 1>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        addition_node_impl* pData = reinterpret_cast<addition_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = pSrc1[i] + pSrc2[i];
    }
};

void addition_node::compile( const std::vector<channel_op_node*>& expressionTree,
                             channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<addition_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

subtraction_node::subtraction_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

subtraction_node::~subtraction_node() {}

struct subtraction_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintMatch<0, 1, ArityConstraintMatch<0, 1>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        subtraction_node_impl* pData = reinterpret_cast<subtraction_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = pSrc1[i] - pSrc2[i];
    }
};

void subtraction_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<subtraction_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

division_node::division_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

division_node::~division_node() {}

struct division_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintMatch<0, 1, ArityConstraintEq<1, 1>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        division_node_impl* pData = reinterpret_cast<division_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );

        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = pSrc1[i] / pSrc2[0];
    }
};

/**
 * A template specialization for integers to prevent divide by 0.
 */
template <>
void division_node_impl::eval<int>( char*, char* stack, void* pVoidData ) {
    division_node_impl* pData = reinterpret_cast<division_node_impl*>( pVoidData );
    int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
    int* pSrc1 = reinterpret_cast<int*>( stack + pData->m_srcIndex[0] );
    int* pSrc2 = reinterpret_cast<int*>( stack + pData->m_srcIndex[1] );

    for( int i = 0; i < pData->m_arity; ++i )
        pDest[i] = ( pSrc2[0] != 0 ) ? ( pSrc1[i] / pSrc2[0] ) : ( pSrc1[i] >= 0 ) ? INT_MAX : INT_MIN;
}

void division_node::compile( const std::vector<channel_op_node*>& expressionTree,
                             channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<division_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

modulo_node::modulo_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

modulo_node::~modulo_node() {}

struct modulo_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintMatch<0, 1, ArityConstraintEq<1, 1>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData );
};

template <>
void modulo_node_impl::eval<float>( char*, char* stack, void* pVoidData ) {
    modulo_node_impl* pData = reinterpret_cast<modulo_node_impl*>( pVoidData );
    float* pDest = reinterpret_cast<float*>( stack + pData->m_destIndex );
    float* pSrc1 = reinterpret_cast<float*>( stack + pData->m_srcIndex[0] );
    float* pSrc2 = reinterpret_cast<float*>( stack + pData->m_srcIndex[1] );

    for( int i = 0; i < pData->m_arity; ++i )
        pDest[i] = std::fmod( pSrc1[i], pSrc2[0] );
}

template <>
void modulo_node_impl::eval<int>( char*, char* stack, void* pVoidData ) {
    modulo_node_impl* pData = reinterpret_cast<modulo_node_impl*>( pVoidData );
    int* pDest = reinterpret_cast<int*>( stack + pData->m_destIndex );
    int* pSrc1 = reinterpret_cast<int*>( stack + pData->m_srcIndex[0] );
    int* pSrc2 = reinterpret_cast<int*>( stack + pData->m_srcIndex[1] );

    for( int i = 0; i < pData->m_arity; ++i )
        pDest[i] = ( pSrc2[0] != 0 ) ? ( pSrc1[i] % pSrc2[0] ) : 0;
}

void modulo_node::compile( const std::vector<channel_op_node*>& expressionTree,
                           channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<modulo_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

sqrt_node::sqrt_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

sqrt_node::~sqrt_node() {}

struct sqrt_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        sqrt_node_impl* pData = reinterpret_cast<sqrt_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i ) {
            pDest[i] = std::sqrt( pSrc[i] );
        }
    }
};

void sqrt_node::compile( const std::vector<channel_op_node*>& expressionTree,
                         channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<sqrt_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

log_node::log_node( int id, int inIndex )
    : operation_node<1>( id, inIndex ) {}

log_node::~log_node() {}

struct log_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        log_node_impl* pData = reinterpret_cast<log_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        for( int i = 0; i < pData->m_arity; ++i ) {
            pDest[i] = log( pSrc[i] );
        }
    }
};

void log_node::compile( const std::vector<channel_op_node*>& expressionTree,
                        channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<log_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

power_node::power_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

power_node::~power_node() {}

struct power_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32, TypeConstraintEq<1, data_type_float32, ArityConstraintEq<1, 1>>>
        Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        power_node_impl* pData = reinterpret_cast<power_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );

        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = std::pow( pSrc1[i], pSrc2[0] );
    }
};

void power_node::compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<power_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

vector_cross_node::vector_cross_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

vector_cross_node::~vector_cross_node() {}

struct vector_cross_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 3;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32,
                             ArityConstraintEq<0, 3, TypeConstraintEq<1, data_type_float32, ArityConstraintEq<1, 3>>>>
        Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        vector_cross_node_impl* pData = reinterpret_cast<vector_cross_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        pDest[0] = pSrc1[1] * pSrc2[2] - pSrc1[2] * pSrc2[1];
        pDest[1] = pSrc1[2] * pSrc2[0] - pSrc1[0] * pSrc2[2];
        pDest[2] = pSrc1[0] * pSrc2[1] - pSrc1[1] * pSrc2[0];
    }
};

void vector_cross_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                 channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<vector_cross_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

vector_dot_node::vector_dot_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

vector_dot_node::~vector_dot_node() {}

struct vector_dot_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintMatch<0, 1, ArityConstraintMatch<0, 1>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        vector_dot_node_impl* pData = reinterpret_cast<vector_dot_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        pDest[0] = pSrc1[0] * pSrc2[0];
        for( int i = 1; i < pData->m_arity; ++i )
            pDest[0] += pSrc1[i] * pSrc2[i];
    }
};

void vector_dot_node::compile( const std::vector<channel_op_node*>& expressionTree,
                               channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<vector_dot_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

vector_to_scalar_node::vector_to_scalar_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

vector_to_scalar_node::~vector_to_scalar_node() {}

// NOTE: This is one-based because the decision was made to match MaxScript which is a stupid piece
//      of software.
struct vector_to_scalar_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintEq<1, data_type_int32, ArityConstraintEq<1, 1>> Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        vector_to_scalar_node_impl* pData = reinterpret_cast<vector_to_scalar_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        int* pSrc2 = reinterpret_cast<int*>( stack + pData->m_srcIndex[1] );

        int index = pSrc2[0];
        if( index <= 0 || index > pData->m_arity )
            // pDest[0] = 0;
            throw std::runtime_error( "ToScalar node: invalid index, must be 1 <= index <= " +
                                      boost::lexical_cast<std::string>( pData->m_arity ) );
        else
            pDest[0] = pSrc1[index - 1];
    }
};

void vector_to_scalar_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                     channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<vector_to_scalar_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

quat_to_vector_node::quat_to_vector_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

quat_to_vector_node::~quat_to_vector_node() {}

struct quat_to_vector_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = 3;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32,
                             ArityConstraintEq<0, 4, TypeConstraintEq<1, data_type_int32, ArityConstraintEq<1, 1>>>>
        Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        quat_to_vector_node_impl* pData = reinterpret_cast<quat_to_vector_node_impl*>( pVoidData );
        frantic::graphics::vector3f* pDest =
            reinterpret_cast<frantic::graphics::vector3f*>( stack + pData->m_destIndex );
        frantic::graphics::quat4f* pSrc1 = reinterpret_cast<frantic::graphics::quat4f*>( stack + pData->m_srcIndex[0] );
        int* pSrc2 = reinterpret_cast<int*>( stack + pData->m_srcIndex[1] );

        int index = pSrc2[0];
        if( index <= 0 || index > pData->m_arity )
            // pDest[0] = 0;
            throw std::runtime_error( "ToScalar node: invalid index, must be 1 <= index <= " +
                                      boost::lexical_cast<std::string>( pData->m_arity ) );

        frantic::graphics::transform4f tm;
        pSrc1->as_transform4f( tm );

        *pDest = tm.get_column( *pSrc2 - 1 );
    }
};

void quat_to_vector_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                   channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<quat_to_vector_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

vectors_to_quat_node::vectors_to_quat_node( int id, int src1Index, int src2Index, int src3Index )
    : operation_node<3>( id, src1Index, src2Index, src3Index ) {}

vectors_to_quat_node::~vectors_to_quat_node() {}

struct vectors_to_quat_node_impl : public code_gen<3> {
    static const int NumSources = 3;
    static const int OutputArity = 4;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<
        0, data_type_float32,
        ArityConstraintEq<
            0, 3,
            TypeConstraintEq<1, data_type_float32,
                             ArityConstraintEq<1, 3, TypeConstraintEq<2, data_type_float32, ArityConstraintEq<2, 3>>>>>>
        Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        using frantic::graphics::quat4f;
        using frantic::graphics::vector3f;

        vectors_to_quat_node_impl* pData = reinterpret_cast<vectors_to_quat_node_impl*>( pVoidData );
        quat4f* pDest = reinterpret_cast<quat4f*>( stack + pData->m_destIndex );
        vector3f* pSrc0 = reinterpret_cast<vector3f*>( stack + pData->m_srcIndex[0] );
        vector3f* pSrc1 = reinterpret_cast<vector3f*>( stack + pData->m_srcIndex[1] );
        vector3f* pSrc2 = reinterpret_cast<vector3f*>( stack + pData->m_srcIndex[2] );

        *pDest = frantic::graphics::quat4f::from_coord_sys( *pSrc0, *pSrc1, *pSrc2 );
    }
};

void vectors_to_quat_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                    channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<vectors_to_quat_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

scalars_to_vector_node::scalars_to_vector_node( int id, int src1Index, int src2Index, int src3Index )
    : operation_node<3>( id, src1Index, src2Index, src3Index ) {}

scalars_to_vector_node::~scalars_to_vector_node() {}

struct scalars_to_vector_node_impl : public code_gen<3> {
    static const int NumSources = 3;
    static const int OutputArity = 3;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintMatch<
        0, 1, TypeConstraintMatch<0, 2, ArityConstraintEq<0, 1, ArityConstraintEq<1, 1, ArityConstraintEq<2, 1>>>>>
        Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        scalars_to_vector_node_impl* pData = reinterpret_cast<scalars_to_vector_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        T* pSrc3 = reinterpret_cast<T*>( stack + pData->m_srcIndex[2] );
        pDest[0] = pSrc1[0];
        pDest[1] = pSrc2[0];
        pDest[2] = pSrc3[0];
    }
};

void scalars_to_vector_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                      channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<scalars_to_vector_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

clamp_node::clamp_node( int id, int src1Index, int src2Index, int src3Index )
    : operation_node<3>( id, src1Index, src2Index, src3Index ) {}

clamp_node::~clamp_node() {}

struct clamp_node_impl : public code_gen<3> {
    static const int NumSources = 3;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintMatch<0, 1, TypeConstraintMatch<0, 2, ArityConstraintEq<1, 1, ArityConstraintEq<2, 1>>>>
        Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        clamp_node_impl* pData = reinterpret_cast<clamp_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        T* pSrc3 = reinterpret_cast<T*>( stack + pData->m_srcIndex[2] );

        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = frantic::math::clamp( pSrc1[i], pSrc2[0], pSrc3[0] );
    }
};

void clamp_node::compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<clamp_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

switch_node::switch_node( int id, int src1Index, int src2Index, int src3Index )
    : operation_node<3>( id, src1Index, src2Index, src3Index ) {}

switch_node::~switch_node() {}

struct switch_node_impl : public code_gen<3> {
    static const int NumSources = 3;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_invalid;
    typedef TypeConstraintMatch<
        0, 1, TypeConstraintEq<2, data_type_int32, ArityConstraintMatch<0, 1, ArityConstraintEq<2, 1>>>>
        Constraints;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        switch_node_impl* pData = reinterpret_cast<switch_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        int* pSrc3 = reinterpret_cast<int*>( stack + pData->m_srcIndex[2] );

        if( pSrc3[0] != 0 )
            memcpy( pDest, pSrc1, sizeof( T ) * pData->m_arity );
        else
            memcpy( pDest, pSrc2, sizeof( T ) * pData->m_arity );
    }
};

void switch_node::compile( const std::vector<channel_op_node*>& expressionTree,
                           channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<switch_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

blend_node::blend_node( int id, int src1Index, int src2Index, int src3Index )
    : operation_node<3>( id, src1Index, src2Index, src3Index ) {}

blend_node::~blend_node() {}

struct blend_node_impl : public code_gen<3> {
    static const int NumSources = 3;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<
        0, data_type_float32,
        TypeConstraintEq<1, data_type_float32,
                         TypeConstraintEq<2, data_type_float32, ArityConstraintMatch<0, 1, ArityConstraintEq<2, 1>>>>>
        Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        blend_node_impl* pData = reinterpret_cast<blend_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        float* pSrc3 = reinterpret_cast<float*>( stack + pData->m_srcIndex[2] );

        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = frantic::math::lerp( pSrc1[i], pSrc2[i], pSrc3[0] );
    }
};

void blend_node::compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData ) {
    inoutCompData.append_code_segment(
        generate_code_segment<blend_node_impl>( m_nodeId, m_srcIndex, expressionTree, inoutCompData ) );
}

bezier_x_to_y_node::bezier_x_to_y_node(
    int id, int src1Index, const std::vector<frantic::math::bezier_curve_point<frantic::graphics2d::vector2f>>& points )
    : operation_node<1>( id, src1Index )
    , m_points( points ) {}

bezier_x_to_y_node::~bezier_x_to_y_node() {}

struct bezier_x_to_y_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 0;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32> Constraints;
    typedef boost::mpl::vector<float> InputTypes;

    std::vector<frantic::math::bezier_curve_point<frantic::graphics2d::vector2f>>* m_points;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        bezier_x_to_y_node_impl* pData = reinterpret_cast<bezier_x_to_y_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );

        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = frantic::math::bezier_curve_x_to_y( *pData->m_points, pSrc[i] );
    }
};

void bezier_x_to_y_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                  channel_operation_compiler& inoutCompData ) {
    typedef bezier_x_to_y_node_impl Op;

    detail::code_segment theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );
    std::vector<frantic::math::bezier_curve_point<frantic::graphics2d::vector2f>>* newPoints =
        new std::vector<frantic::math::bezier_curve_point<frantic::graphics2d::vector2f>>( m_points );
    inoutCompData.register_scoped_object( newPoints );
    reinterpret_cast<Op*>( theSeg.codeData )->m_points = newPoints;

    inoutCompData.append_code_segment( theSeg );
}

geometry_input_node::geometry_input_node( int id )
    : channel_op_node( id ) {}

geometry_input_node::~geometry_input_node() {}

void geometry_input_node::compile( const std::vector<channel_op_node*>&, channel_operation_compiler& ) {}

void geometry_input_node::add_mesh( boost::shared_ptr<frantic::geometry::trimesh3> pMesh,
                                    boost::shared_ptr<frantic::geometry::trimesh3_kdtree> pTree ) {
    m_meshes.push_back( pMesh );
    m_meshTrees.push_back( pTree );
}

// this is what will be given as the output to all geometry queries
struct geometry_query_result {
    graphics::vector3f position;
    float distance;
    graphics::vector3f faceNormal;
    boost::int32_t meshIndex;
    boost::int32_t faceIndex;
    graphics::vector3f barycentricCoords;
    boost::int32_t intersected; // boolean result stating if there was an intersection or not
};

// void geometry_node::add_mesh( boost::shared_ptr< frantic::geometry::trimesh3 > pMesh, boost::shared_ptr<
// frantic::geometry::trimesh3_kdtree > pTree  ){
//	m_meshes.push_back( pMesh );
//	m_trees.push_back( pTree );
//}

ray_intersect_node::ray_intersect_node( int id, int geomSrcIndex, int src1Index, int src2Index )
    : operation_node<2>( id, src1Index, src2Index )
    , geometry_node( geomSrcIndex ) {}

ray_intersect_node::~ray_intersect_node() {}

struct ray_intersect_node_impl : public code_gen<2> {
    static const int NumSources = 2;
    static const int OutputArity = ( sizeof( geometry_query_result ) + 3 ) / 4;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32,
                             TypeConstraintEq<1, data_type_float32, ArityConstraintEq<0, 3, ArityConstraintEq<1, 3>>>>
        Constraints;

    typedef boost::mpl::vector<float> InputTypes;

    std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>* m_meshes;
    std::vector<boost::shared_ptr<frantic::geometry::trimesh3_kdtree>>* m_meshTrees;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        ray_intersect_node_impl* pData = reinterpret_cast<ray_intersect_node_impl*>( pVoidData );

        const graphics::vector3f& inputPos = *reinterpret_cast<graphics::vector3f*>( stack + pData->m_srcIndex[0] );
        const graphics::vector3f& inputDir = *reinterpret_cast<graphics::vector3f*>( stack + pData->m_srcIndex[1] );

        geometry_query_result& result = *reinterpret_cast<geometry_query_result*>( stack + pData->m_destIndex );
        memset( &result, 0, sizeof( channels::geometry_query_result ) );

        graphics::ray3f queryRay( inputPos, inputDir );
        geometry::raytrace_intersection raytraceResult;

        double closestResult = std::numeric_limits<float>::max();
        for( size_t i = 0; i < pData->m_meshTrees->size(); ++i ) {
            if( ( *pData->m_meshTrees )[i]->intersect_ray( queryRay, 0.0f, (float)closestResult, raytraceResult ) &&
                raytraceResult.distance < closestResult ) {
                closestResult = raytraceResult.distance; // update the closest result distance

                result.intersected = 1;
                result.position = raytraceResult.position;
                result.distance = vector3f::dot( queryRay.direction(), raytraceResult.geometricNormal ) <= 0
                                      ? (float)raytraceResult.distance
                                      : -(float)raytraceResult.distance;
                result.faceNormal = raytraceResult.geometricNormal;
                result.meshIndex = (boost::int32_t)i;
                result.faceIndex = raytraceResult.faceIndex;
                result.barycentricCoords = raytraceResult.barycentricCoords;
            }
        }
    }
};

void ray_intersect_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                  channel_operation_compiler& inoutCompData ) {
    typedef ray_intersect_node_impl Op;

    int inputId = get_geometry_input_index();
    if( inputId < 0 || inputId >= (int)expressionTree.size() )
        throw channel_compiler_error( m_nodeId, "Invalid connection id for connection #" +
                                                    boost::lexical_cast<std::string>( inputId + 1 ) );

    geometry_input_node* pGeomInput = dynamic_cast<geometry_input_node*>( expressionTree[get_geometry_input_index()] );
    if( !pGeomInput )
        throw channel_compiler_error( m_nodeId, "Expected Geometry input in connection #" +
                                                    boost::lexical_cast<std::string>( 0 ) );

    detail::code_segment theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );

    Op* pOp = reinterpret_cast<Op*>( theSeg.codeData );
    pOp->m_meshes = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>( pGeomInput->get_meshes() );
    pOp->m_meshTrees =
        new std::vector<boost::shared_ptr<frantic::geometry::trimesh3_kdtree>>( pGeomInput->get_kdtrees() );
    // pOp->m_meshes = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3> >( m_meshes );
    // pOp->m_meshTrees = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3_kdtree> >( m_trees );

    inoutCompData.register_scoped_object( pOp->m_meshes );
    inoutCompData.register_scoped_object( pOp->m_meshTrees );
    inoutCompData.append_code_segment( theSeg );
}

nearest_point_node::nearest_point_node( int id, int geomSrcIndex, int src1Index )
    : operation_node<1>( id, src1Index )
    , geometry_node( geomSrcIndex ) {}

nearest_point_node::~nearest_point_node() {}

struct nearest_point_node_impl : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = ( sizeof( geometry_query_result ) + 3 ) / 4;
    static const data_type_t OutputType = data_type_float32;
    typedef TypeConstraintEq<0, data_type_float32, ArityConstraintEq<0, 3>> Constraints;

    typedef boost::mpl::vector<float> InputTypes;

    std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>* m_meshes;
    std::vector<boost::shared_ptr<frantic::geometry::trimesh3_kdtree>>* m_meshTrees;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        nearest_point_node_impl* pData = reinterpret_cast<nearest_point_node_impl*>( pVoidData );

        const graphics::vector3f& inputPos = *reinterpret_cast<graphics::vector3f*>( stack + pData->m_srcIndex[0] );

        geometry_query_result& result = *reinterpret_cast<geometry_query_result*>( stack + pData->m_destIndex );
        memset( &result, 0, sizeof( channels::geometry_query_result ) );

        geometry::nearest_point_search_result npsResult;

        double closestResult = std::numeric_limits<float>::max();
        for( size_t i = 0; i < pData->m_meshTrees->size(); ++i ) {
            if( ( *pData->m_meshTrees )[i]->find_nearest_point( inputPos, (float)closestResult, npsResult ) &&
                npsResult.distance < closestResult ) {
                closestResult = npsResult.distance; // update the closest result distance

                result.intersected = 1;
                result.position = npsResult.position;
                result.distance = vector3f::dot( ( npsResult.position - inputPos ), npsResult.geometricNormal ) <= 0
                                      ? (float)npsResult.distance
                                      : -(float)npsResult.distance;
                result.faceNormal = npsResult.geometricNormal;
                result.meshIndex = (boost::int32_t)i;
                result.faceIndex = npsResult.faceIndex;
                result.barycentricCoords = npsResult.barycentricCoords;
            }
        }
    }
};

void nearest_point_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                  channel_operation_compiler& inoutCompData ) {
    typedef nearest_point_node_impl Op;

    int inputId = get_geometry_input_index();
    if( inputId < 0 || inputId >= (int)expressionTree.size() )
        throw channel_compiler_error( m_nodeId, "Invalid connection id for connection #" +
                                                    boost::lexical_cast<std::string>( inputId + 1 ) );

    geometry_input_node* pGeomInput = dynamic_cast<geometry_input_node*>( expressionTree[get_geometry_input_index()] );
    if( !pGeomInput )
        throw channel_compiler_error( m_nodeId, "Expected Geometry input in connection #" +
                                                    boost::lexical_cast<std::string>( 0 ) );

    detail::code_segment theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );

    Op* pOp = reinterpret_cast<Op*>( theSeg.codeData );
    pOp->m_meshes = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>( pGeomInput->get_meshes() );
    pOp->m_meshTrees =
        new std::vector<boost::shared_ptr<frantic::geometry::trimesh3_kdtree>>( pGeomInput->get_kdtrees() );
    // pOp->m_meshes = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3> >( m_meshes );
    // Op->m_meshTrees = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3_kdtree> >( m_trees );

    inoutCompData.register_scoped_object( pOp->m_meshes );
    inoutCompData.register_scoped_object( pOp->m_meshTrees );
    inoutCompData.append_code_segment( theSeg );
}

surf_data_value_node::surf_data_value_node( int id, int src1Index, target_result target )
    : operation_node<1>( id, src1Index )
    , m_target( target ) {}

surf_data_value_node::~surf_data_value_node() {}

static const std::string MATERIAL_CHANNEL_NAME = "MaterialID";
static const std::string SMOOTH_GROUP_CHANNEL_NAME = "SmoothingGroup";
static const std::string TEXTURE_COORD_CHANNEL_NAME = "TextureCoord";
static const std::string NORMAL_CHANNEL_NAME = "Normals"; // probably not the correct chnanel name

/*
 I'm not entirely sure what the correct approach here would have been, but attempting to template this on
 multiple types did not work very well, so I just made an impl-class for each possible combination of
 outputs.
*/

// for single integer outputs
struct surf_data_value_node_impl_int1 : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_int32;
    typedef EmptyConstraint Constraints;

    typedef boost::mpl::vector<float> InputTypes;

    std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>* m_meshes;
    surf_data_value_node::target_result m_target;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        surf_data_value_node_impl_int1* pData = reinterpret_cast<surf_data_value_node_impl_int1*>( pVoidData );

        geometry_query_result& prevResult = *reinterpret_cast<geometry_query_result*>( stack + pData->m_srcIndex[0] );

        boost::int32_t& result = *reinterpret_cast<boost::int32_t*>( stack + pData->m_destIndex );

        switch( pData->m_target ) {
        case surf_data_value_node::RESULT_MESH_INDEX:
            result = prevResult.meshIndex;
            break;
        case surf_data_value_node::RESULT_FACE_INDEX:
            result = prevResult.faceIndex;
            break;
        case surf_data_value_node::RESULT_FACE_MAT_ID: {
            const frantic::geometry::trimesh3& targetMesh = *pData->m_meshes->at( prevResult.meshIndex );
            frantic::geometry::const_trimesh3_face_channel_general_accessor materialIDAccess =
                targetMesh.get_face_channel_general_accessor( frantic::strings::to_tstring( MATERIAL_CHANNEL_NAME ) );
            result = *reinterpret_cast<const unsigned short*>( materialIDAccess.data( prevResult.faceIndex ) );
            break;
        } break;
        case surf_data_value_node::RESULT_SMOOTH_GROUP: {
            const frantic::geometry::trimesh3& targetMesh = *pData->m_meshes->at( prevResult.meshIndex );
            frantic::geometry::const_trimesh3_face_channel_general_accessor smoothGroupAccess =
                targetMesh.get_face_channel_general_accessor(
                    frantic::strings::to_tstring( SMOOTH_GROUP_CHANNEL_NAME ) );
            result = *reinterpret_cast<const int*>( smoothGroupAccess.data( prevResult.faceIndex ) );
        } break;
        case surf_data_value_node::RESULT_INTERSECTED:
            result = prevResult.intersected;
            break;
        default:
#ifdef _WIN32
            __assume( 0 );
#else
            break;
#endif
        }
    }
};

// for single float outputs
struct surf_data_value_node_impl_float1 : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 1;
    static const data_type_t OutputType = data_type_float32;
    typedef EmptyConstraint Constraints;

    typedef boost::mpl::vector<float> InputTypes;

    std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>* m_meshes;
    surf_data_value_node::target_result m_target;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        surf_data_value_node_impl_float1* pData = reinterpret_cast<surf_data_value_node_impl_float1*>( pVoidData );

        geometry_query_result& prevResult = *reinterpret_cast<geometry_query_result*>( stack + pData->m_srcIndex[0] );

        float& result = *reinterpret_cast<float*>( stack + pData->m_destIndex );

        result = prevResult.distance;
    }
};

// for vector float outputs
struct surf_data_value_node_impl_float3 : public code_gen<1> {
    static const int NumSources = 1;
    static const int OutputArity = 3;
    static const data_type_t OutputType = data_type_float32;
    typedef EmptyConstraint Constraints;

    typedef boost::mpl::vector<float> InputTypes;

    std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>* m_meshes;
    surf_data_value_node::target_result m_target;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        surf_data_value_node_impl_float3* pData = reinterpret_cast<surf_data_value_node_impl_float3*>( pVoidData );

        geometry_query_result& prevResult = *reinterpret_cast<geometry_query_result*>( stack + pData->m_srcIndex[0] );

        frantic::graphics::vector3f& result =
            *reinterpret_cast<frantic::graphics::vector3f*>( stack + pData->m_destIndex );

        switch( pData->m_target ) {
        case surf_data_value_node::RESULT_POSITION:
            result = prevResult.position;
            break;
        case surf_data_value_node::RESULT_FACE_NORMAL:
            result = prevResult.faceNormal;
            break;
        case surf_data_value_node::RESULT_SMOOTH_NORMAL:
            // not implemented
            break;
        case surf_data_value_node::RESULT_BARY_COORDS:
            result = prevResult.barycentricCoords;
            break;
        case surf_data_value_node::RESULT_TEXTURE_COORD: {
            const frantic::geometry::trimesh3& targetMesh = *pData->m_meshes->at( prevResult.meshIndex );
            frantic::geometry::const_trimesh3_vertex_channel_general_accessor texCoordsAccess =
                targetMesh.get_vertex_channel_general_accessor(
                    frantic::strings::to_tstring( TEXTURE_COORD_CHANNEL_NAME ) );
            const frantic::graphics::vector3 vertexIDs = texCoordsAccess.face( prevResult.faceIndex );
            const frantic::graphics::vector3f& texCoord0 =
                *reinterpret_cast<const frantic::graphics::vector3f*>( texCoordsAccess.data( vertexIDs[0] ) );
            const frantic::graphics::vector3f& texCoord1 =
                *reinterpret_cast<const frantic::graphics::vector3f*>( texCoordsAccess.data( vertexIDs[1] ) );
            const frantic::graphics::vector3f& texCoord2 =
                *reinterpret_cast<const frantic::graphics::vector3f*>( texCoordsAccess.data( vertexIDs[2] ) );

            result = prevResult.barycentricCoords[0] * texCoord0 + prevResult.barycentricCoords[1] * texCoord1 +
                     prevResult.barycentricCoords[2] * texCoord2;
        } break;
        default:
#ifdef _WIN32
            __assume( 0 );
#else
            break;
#endif
        }
    }
};

// a pair of helper routines when checking if a given channel exists in a mesh
static void check_face_channel( int m_nodeId, const frantic::geometry::trimesh3& targetMesh,
                                const std::string& channelName, int idx ) {
    if( !targetMesh.has_face_channel( frantic::strings::to_tstring( channelName ) ) ) {
        throw channel_compiler_error( m_nodeId, "Channel: \"" + channelName + "\" not availiable in mesh #" +
                                                    boost::lexical_cast<std::string>( idx + 1 ) );
    }
}

static void check_vertex_channel( int m_nodeId, const frantic::geometry::trimesh3& targetMesh,
                                  const std::string& channelName, int idx ) {
    if( !targetMesh.has_vertex_channel( frantic::strings::to_tstring( channelName ) ) ) {
        throw channel_compiler_error( m_nodeId, "Channel: \"" + channelName + "\" not availiable in mesh #" +
                                                    boost::lexical_cast<std::string>( idx + 1 ) );
    }
}

void surf_data_value_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                    channel_operation_compiler& inoutCompData ) {
    const int NumSources = 1;
    for( int i = 0; i < NumSources; ++i ) {
        if( m_srcIndex[i] < 0 || m_srcIndex[i] >= (int)expressionTree.size() )
            throw channel_compiler_error( m_nodeId, "Invalid connection id for connection #" +
                                                        boost::lexical_cast<std::string>( i + 1 ) );
    }

    const std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>* pMeshes = NULL;

    // check if the 1st input is one of the geometry query types
    if( geometry_node* pGeomNode = dynamic_cast<geometry_node*>( expressionTree[m_srcIndex[0]] ) ) {
        // pMeshes = &pGeomNode->get_meshes();
        if( geometry_input_node* pGeomInput =
                dynamic_cast<geometry_input_node*>( expressionTree[pGeomNode->get_geometry_input_index()] ) ) {
            pMeshes = &pGeomInput->get_meshes();
        } else
            throw channel_compiler_error(
                m_nodeId,
                "Error, input to surface data node was not a geometry query node (Ex. RayIntersect or NearestPoint)" );
    } else {
        throw channel_compiler_error(
            m_nodeId,
            "Error, input to surface data node was not a geometry query node (Ex. RayIntersect or NearestPoint)" );
    }

    // check if all of the neccessary channels are availiable for the target output value
    for( size_t i = 0; i < pMeshes->size(); ++i ) {
        switch( m_target ) {
        case RESULT_SMOOTH_NORMAL:
            throw channel_compiler_error( m_nodeId, "Not implemented yet" );
            check_face_channel( m_nodeId, *pMeshes->at( i ), SMOOTH_GROUP_CHANNEL_NAME, (int)i );
            check_vertex_channel( m_nodeId, *pMeshes->at( i ), NORMAL_CHANNEL_NAME, (int)i );
            break;
        case RESULT_FACE_MAT_ID:
            check_face_channel( m_nodeId, *pMeshes->at( i ), MATERIAL_CHANNEL_NAME, (int)i );
            break;
        case RESULT_SMOOTH_GROUP:
            check_face_channel( m_nodeId, *pMeshes->at( i ), SMOOTH_GROUP_CHANNEL_NAME, (int)i );
            break;
        case RESULT_TEXTURE_COORD:
            check_vertex_channel( m_nodeId, *pMeshes->at( i ), TEXTURE_COORD_CHANNEL_NAME, (int)i );
            break;
        default:
            break;
        }
    }

    detail::code_segment theSeg;

    // this will insert a new code segment of the appropriate type based on the target value.  Its a reasonably
    // poor solution, but it works
    if( m_target == RESULT_MESH_INDEX || // single integer types
        m_target == RESULT_FACE_INDEX || m_target == RESULT_FACE_MAT_ID || m_target == RESULT_SMOOTH_GROUP ||
        m_target == RESULT_INTERSECTED ) {
        typedef surf_data_value_node_impl_int1 Op;

        theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );

        Op* pOp = reinterpret_cast<Op*>( theSeg.codeData );
        pOp->m_meshes = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>( *pMeshes );
        pOp->m_target = m_target;

        inoutCompData.register_scoped_object( pOp->m_meshes );
    } else if( m_target == RESULT_POSITION || // vec3 types
               m_target == RESULT_FACE_NORMAL || m_target == RESULT_SMOOTH_NORMAL || m_target == RESULT_BARY_COORDS ||
               m_target == RESULT_TEXTURE_COORD ) {
        typedef surf_data_value_node_impl_float3 Op;

        theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );

        Op* pOp = reinterpret_cast<Op*>( theSeg.codeData );
        pOp->m_meshes = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>( *pMeshes );
        pOp->m_target = m_target;

        inoutCompData.register_scoped_object( pOp->m_meshes );
    } else if( m_target == RESULT_SIGNED_DIST ) // single float types
    {
        typedef surf_data_value_node_impl_float1 Op;

        theSeg = generate_code_segment<Op>( m_nodeId, m_srcIndex, expressionTree, inoutCompData );

        Op* pOp = reinterpret_cast<Op*>( theSeg.codeData );
        pOp->m_meshes = new std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>( *pMeshes );
        pOp->m_target = m_target;

        inoutCompData.register_scoped_object( pOp->m_meshes );
    } else {
        throw channel_compiler_error( m_nodeId, "Error, invalid surface data value target: " +
                                                    boost::lexical_cast<std::string>( m_target ) );
    }

    inoutCompData.append_code_segment( theSeg );
}

//****************************************
// Multiplication is a little weird, so we
// generate the code segment ourselves.
//****************************************

multiplication_node::multiplication_node( int id, int lhsIndex, int rhsIndex )
    : operation_node<2>( id, lhsIndex, rhsIndex ) {}

multiplication_node::~multiplication_node() {}

struct multiplication_node_impl : public code_gen<2> {
    static const int NumInputs = 2;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        multiplication_node_impl* pData = reinterpret_cast<multiplication_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = pSrc1[i] * pSrc2[0];
    }
};

struct component_multiplication_node_impl : public code_gen<2> {
    static const int NumInputs = 2;
    typedef boost::mpl::vector<float, int> InputTypes;

    template <class T>
    static void eval( char*, char* stack, void* pVoidData ) {
        component_multiplication_node_impl* pData = reinterpret_cast<component_multiplication_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destIndex );
        T* pSrc1 = reinterpret_cast<T*>( stack + pData->m_srcIndex[0] );
        T* pSrc2 = reinterpret_cast<T*>( stack + pData->m_srcIndex[1] );
        for( int i = 0; i < pData->m_arity; ++i )
            pDest[i] = pSrc1[i] * pSrc2[i];
    }
};

void multiplication_node::compile( const std::vector<channel_op_node*>& expressionTree,
                                   channel_operation_compiler& inoutCompData ) {
    temporary_result* r = inoutCompData.get_node_results( m_nodeId );
    if( r != NULL )
        throw channel_compiler_error( m_nodeId, "Graph cycle detected" );

    const int NumSources = 2;

    temporary_result* src[NumSources];
    for( int i = 0; i < NumSources; ++i ) {
        if( m_srcIndex[i] < 0 || m_srcIndex[i] >= (int)expressionTree.size() )
            throw channel_compiler_error( m_nodeId, "Invalid connection id for connection #" +
                                                        boost::lexical_cast<std::string>( i + 1 ) );
        src[i] = get_or_create_node_result( m_srcIndex[i], expressionTree, inoutCompData );
    }

    if( src[0]->arity < src[1]->arity ) {
        TypeConstraintMatch<0, 1, ArityConstraintMatchOrEq<1, 0, 1>>::check( m_nodeId, src );
        std::swap( src[0], src[1] );
    } else {
        TypeConstraintMatch<0, 1, ArityConstraintMatchOrEq<0, 1, 1>>::check( m_nodeId, src );
    }

    r = inoutCompData.allocate_temporary( m_nodeId, src[0]->type, src[0]->arity );

    if( src[1]->arity == 1 ) {
        typedef multiplication_node_impl Op;

        detail::code_segment theSeg;
        theSeg.codePtr = bind_template_fn<Op, Op::InputTypes>( src[0]->type );
        theSeg.codeData = (Op*)malloc( sizeof( Op ) );
        reinterpret_cast<Op*>( theSeg.codeData )->init( src[0]->arity, r->offset, src[0]->offset, src[1]->offset );

        inoutCompData.append_code_segment( theSeg );
    } else {
        typedef component_multiplication_node_impl Op;

        detail::code_segment theSeg;
        theSeg.codePtr = bind_template_fn<Op, Op::InputTypes>( src[0]->type );
        theSeg.codeData = (Op*)malloc( sizeof( Op ) );
        reinterpret_cast<Op*>( theSeg.codeData )->init( src[0]->arity, r->offset, src[0]->offset, src[1]->offset );

        inoutCompData.append_code_segment( theSeg );
    }
}
} // namespace channels
} // namespace frantic
