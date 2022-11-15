// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/mpl/assert.hpp>
#include <frantic/channels/channel_operation_compiler.hpp>
#include <frantic/math/splines/bezier_spline.hpp>

#include <vector>

#include <frantic/geometry/trimesh3.hpp>
#include <frantic/geometry/trimesh3_kdtree.hpp>

#include <boost/shared_ptr.hpp>

namespace frantic {
namespace channels {

/**
 * This AST node will just forward the output of another node.
 */
class disabled_op_node : public channel_op_node {
    int m_passThroughIndex;

  public:
    disabled_op_node( int id, int passThroughIndex );
    virtual ~disabled_op_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will extract a value from a channel in the channel map of the compiler.
 */
class input_channel_op_node : public channel_op_node {
    std::string m_inputChannel;

  public:
    input_channel_op_node( int id, const std::string& inputChannel );
    virtual ~input_channel_op_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );

    const std::string& get_channel() const { return m_inputChannel; }
};

/**
 * This AST node will insert a value into a channel of the channel map of the compiler.
 */
class output_channel_op_node : public channel_op_node {
    std::string m_outputChannel;
    data_type_t m_outputType;
    int m_outputArity;
    int m_inputIndex;

  public:
    output_channel_op_node( int id, int inputIndex, const std::string& outChannel, int outArity, data_type_t outType );
    virtual ~output_channel_op_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );

    const std::string& get_channel_name() const;
};

/**
 * This AST node will put a constant value into the temporaries stack.
 */
template <class T, int Arity>
class input_value_op_node : public channel_op_node {
    T m_value[Arity];

  public:
    inline input_value_op_node( int id, T value[] );
    inline virtual ~input_value_op_node();
    inline virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                                 channel_operation_compiler& inoutCompData );
};

/**
 * This is a base class for AST nodes that have child nodes that must be compiled
 * before the current node can be compiled.
 */
template <int NumInputs>
class operation_node : public channel_op_node {
  protected:
    int m_srcIndex[NumInputs];

  public:
    operation_node( int nodeId, int src1Index )
        : channel_op_node( nodeId ) {
        boost::mpl::assert<NumInputs == 1>();
        m_srcIndex[0] = src1Index;
    }
    operation_node( int nodeId, int src1Index, int src2Index )
        : channel_op_node( nodeId ) {
        boost::mpl::assert<NumInputs == 2>();
        m_srcIndex[0] = src1Index;
        m_srcIndex[1] = src2Index;
    }
    operation_node( int nodeId, int src1Index, int src2Index, int src3Index )
        : channel_op_node( nodeId ) {
        boost::mpl::assert<NumInputs == 3>();
        m_srcIndex[0] = src1Index;
        m_srcIndex[1] = src2Index;
        m_srcIndex[2] = src3Index;
    }
    virtual ~operation_node() {}
};

/**
 * This AST node will convert a stack temporary into another type. You can specify an arity
 * by setting OutArity to a positive value. If OutArity is 0, any arity is accepted.
 */
template <class DestType, unsigned OutArity = 0>
class convert_node : public operation_node<1> {
  public:
    inline convert_node( int id, int inIndex );
    inline virtual ~convert_node();
    inline virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                                 channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will negate the input value, component-wise.
 */
class negate_node : public operation_node<1> {
  public:
    negate_node( int id, int inIndex );
    virtual ~negate_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will sum the components of the input value.
 */
class component_sum_node : public operation_node<1> {
  public:
    component_sum_node( int id, int inIndex );
    virtual ~component_sum_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will transform a Vector/Quaternion input by a Quaternion input
 */
class transform_by_quat_node : public operation_node<2> {
  public:
    transform_by_quat_node( int id, int src1Index, int scr2Index );
    virtual ~transform_by_quat_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the magnitude (sqrt of the sum of the squares) of the input value.
 */
class magnitude_node : public operation_node<1> {
  public:
    magnitude_node( int id, int inIndex );
    virtual ~magnitude_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will normalize the input value, such that its magnitude will be 1 afterwards.
 */
class normalize_node : public operation_node<1> {
  public:
    normalize_node( int id, int inIndex );
    virtual ~normalize_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the component-wise absolute value of the input.
 */
class abs_node : public operation_node<1> {
  public:
    abs_node( int id, int inIndex );
    virtual ~abs_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the nearest integer less than or equal to than the input.
 */
class floor_node : public operation_node<1> {
  public:
    floor_node( int id, int inIndex );
    virtual ~floor_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the nearest integer greater than or equal to than the input.
 */
class ceil_node : public operation_node<1> {
  public:
    ceil_node( int id, int inIndex );
    virtual ~ceil_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the sine() of the input value (as radians).
 */
class sin_trig_node : public operation_node<1> {
  public:
    sin_trig_node( int id, int inIndex );
    virtual ~sin_trig_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the cosine() of the input value (as radians).
 */
class cos_trig_node : public operation_node<1> {
  public:
    cos_trig_node( int id, int inIndex );
    virtual ~cos_trig_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the tangent() of the input value (as radians).
 */
class tan_trig_node : public operation_node<1> {
  public:
    tan_trig_node( int id, int inIndex );
    virtual ~tan_trig_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the inverse sine() of the input value (as radians).
 */
class asin_trig_node : public operation_node<1> {
  public:
    asin_trig_node( int id, int inIndex );
    virtual ~asin_trig_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the inverse cosine() of the input value (as radians).
 */
class acos_trig_node : public operation_node<1> {
  public:
    acos_trig_node( int id, int inIndex );
    virtual ~acos_trig_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the inverse tangent() of the input value (as radians).
 */
class atan_trig_node : public operation_node<1> {
  public:
    atan_trig_node( int id, int inIndex );
    virtual ~atan_trig_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the inverse tangent() of the input values (as radians). This version of atan() is capable
 * of returning an angle in the range [-pi, pi].
 */
class atan2_trig_node : public operation_node<2> {
  public:
    atan2_trig_node( int id, int srcIndex1, int srcIndex2 );
    virtual ~atan2_trig_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will transform the input point by a 4x4 transformation matrix.
 */
class transform_point_node : public operation_node<1> {
    frantic::graphics::transform4f m_transform;

  public:
    transform_point_node( int id, int src1Index, const frantic::graphics::transform4f& tm );
    virtual ~transform_point_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will transform the input direction by a 4x4 transformation matrix.
 */
class transform_vector_node : public operation_node<1> {
    frantic::graphics::transform4f m_transform;

  public:
    transform_vector_node( int id, int src1Index, const frantic::graphics::transform4f& tm );
    virtual ~transform_vector_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will transform the input normal by a 4x4 transformation matrix.
 */
class transform_normal_node : public operation_node<1> {
    frantic::graphics::transform4f m_transform;

  public:
    transform_normal_node( int id, int src1Index, const frantic::graphics::transform4f& tm );
    virtual ~transform_normal_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will generate a perlin noise value for the given input.
 */
class noise_node : public operation_node<3> {
  public:
    noise_node( int id, int src1Index, int src2Index, int src3Index );
    virtual ~noise_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 *  This AST node gets the derivative(s) of a perlin nose function for a given value
 */

class dnoise_node : public operation_node<3> {
  public:
    dnoise_node( int id, int src1Index, int src2Index, int src3Index );
    virtual ~dnoise_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will less-than compare the two inputs, producing an integer boolean.
 */
class less_comparison_node : public operation_node<2> {
  public:
    less_comparison_node( int id, int lhsIndex, int rhsIndex );
    virtual ~less_comparison_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will less-than-or-equal compare the two inputs, producing an integer boolean.
 */
class less_or_equal_comparison_node : public operation_node<2> {
  public:
    less_or_equal_comparison_node( int id, int lhsIndex, int rhsIndex );
    virtual ~less_or_equal_comparison_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will greater-than compare the two inputs, producing an integer boolean.
 */
class greater_comparison_node : public operation_node<2> {
  public:
    greater_comparison_node( int id, int lhsIndex, int rhsIndex );
    virtual ~greater_comparison_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will greater-than-or-equal compare the two inputs, producing an integer boolean.
 */
class greater_or_equal_comparison_node : public operation_node<2> {
  public:
    greater_or_equal_comparison_node( int id, int lhsIndex, int rhsIndex );
    virtual ~greater_or_equal_comparison_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will equals compare the two inputs, producing an integer boolean.
 */
class equal_comparison_node : public operation_node<2> {
  public:
    equal_comparison_node( int id, int lhsIndex, int rhsIndex );
    virtual ~equal_comparison_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will logical AND the two boolean integer inputs.
 */
class logical_and_node : public operation_node<2> {
  public:
    logical_and_node( int id, int lhsIndex, int rhsIndex );
    virtual ~logical_and_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will logical OR the two boolean integer inputs.
 */
class logical_or_node : public operation_node<2> {
  public:
    logical_or_node( int id, int lhsIndex, int rhsIndex );
    virtual ~logical_or_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will logical NOT the boolean integer input.
 */
class logical_not_node : public operation_node<1> {
  public:
    logical_not_node( int id, int inIndex );
    virtual ~logical_not_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will sum the two inputs.
 */
class addition_node : public operation_node<2> {
  public:
    addition_node( int id, int lhsIndex, int rhsIndex );
    virtual ~addition_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will subtract the two inputs.
 */
class subtraction_node : public operation_node<2> {
  public:
    subtraction_node( int id, int lhsIndex, int rhsIndex );
    virtual ~subtraction_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will multiply the two inputs.
 */
class multiplication_node : public operation_node<2> {
  public:
    multiplication_node( int id, int lhsIndex, int rhsIndex );
    virtual ~multiplication_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will divide the two inputs.
 */
class division_node : public operation_node<2> {
  public:
    division_node( int id, int lhsIndex, int rhsIndex );
    virtual ~division_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the modulus (ie. remainder) of the two inputs.
 */
class modulo_node : public operation_node<2> {
  public:
    modulo_node( int id, int lhsIndex, int rhsIndex );
    virtual ~modulo_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};
/**
 * This AST node will compute the principle square root of x
 */
class sqrt_node : public operation_node<1> {
  public:
    sqrt_node( int id, int inIndex );
    virtual ~sqrt_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the natural logarithm of x
 */
class log_node : public operation_node<1> {
  public:
    log_node( int id, int inIndex );
    virtual ~log_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the exponent x^y, where the first input is raised to the power of the second.
 */
class power_node : public operation_node<2> {
  public:
    power_node( int id, int lhsIndex, int rhsIndex );
    virtual ~power_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the cross product of the input vectos.
 */
class vector_cross_node : public operation_node<2> {
  public:
    vector_cross_node( int id, int lhsIndex, int rhsIndex );
    virtual ~vector_cross_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will compute the cross product of the input vectos.
 */
class vector_dot_node : public operation_node<2> {
  public:
    vector_dot_node( int id, int lhsIndex, int rhsIndex );
    virtual ~vector_dot_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will extract a scalar component of a vector. NOTE: This operation is 1 based.
 */
class vector_to_scalar_node : public operation_node<2> {
  public:
    vector_to_scalar_node( int id, int lhsIndex, int rhsIndex );
    virtual ~vector_to_scalar_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will extrate a vector component of a quaternion. NOTE: This operation is 1 based.
 * NOTE: converting from vectors_to_quat(quat_to_vector(quat,1),quat_to_vector(quat,2),quat_to_vector(quat,3))
 * cannot be guarenteed to produce the initial quaternion.
 */
class quat_to_vector_node : public operation_node<2> {
  public:
    quat_to_vector_node( int id, int lhsIndex, int rhsIndex );
    virtual ~quat_to_vector_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will create a quaternion from three input vectors.
 * NOTE: converting from vectors_to_quat(quat_to_vector(quat,1),quat_to_vector(quat,2),quat_to_vector(quat,3))
 * cannot be guarenteed to produce the initial quaternion.
 */
class vectors_to_quat_node : public operation_node<3> {
  public:
    vectors_to_quat_node( int id, int src1Index, int src2Index, int src3Index );
    virtual ~vectors_to_quat_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will create a vector triple from three input scalars.
 */
class scalars_to_vector_node : public operation_node<3> {
  public:
    scalars_to_vector_node( int id, int src1Index, int src2Index, int src3Index );
    virtual ~scalars_to_vector_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will clamp the first value to be greater than the second and less than the third.
 */
class clamp_node : public operation_node<3> {
  public:
    clamp_node( int id, int src1Index, int src2Index, int src3Index );
    virtual ~clamp_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will pass through input 1 if input 3 is non-zero, otherwise it will
 * pass through input 2.
 */
class switch_node : public operation_node<3> {
  public:
    switch_node( int id, int src1Index, int src2Index, int src3Index );
    virtual ~switch_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will linearly blend between the first and second inputs using the third input
 * as the blend parameter.
 */
class blend_node : public operation_node<3> {
  public:
    blend_node( int id, int src1Index, int src2Index, int src3Index );
    virtual ~blend_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

class bezier_x_to_y_node : public operation_node<1> {
    std::vector<frantic::math::bezier_curve_point<frantic::graphics2d::vector2f>> m_points;

  public:
    bezier_x_to_y_node( int id, int src1Index,
                        const std::vector<frantic::math::bezier_curve_point<frantic::graphics2d::vector2f>>& points );
    virtual ~bezier_x_to_y_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node represents a set of meshes to be queried
 */
class geometry_input_node : public channel_op_node {
  private:
    std::vector<boost::shared_ptr<frantic::geometry::trimesh3>> m_meshes;
    std::vector<boost::shared_ptr<frantic::geometry::trimesh3_kdtree>> m_meshTrees;

    friend class ray_intersect_node;
    friend class nearest_point_node;
    friend class surf_data_value_node;

  public:
    geometry_input_node( int id );
    virtual ~geometry_input_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );

    void add_mesh( boost::shared_ptr<frantic::geometry::trimesh3> pMesh,
                   boost::shared_ptr<frantic::geometry::trimesh3_kdtree> pTree );

    const std::vector<boost::shared_ptr<frantic::geometry::trimesh3>>& get_meshes() const { return m_meshes; }

    const std::vector<boost::shared_ptr<frantic::geometry::trimesh3_kdtree>>& get_kdtrees() const {
        return m_meshTrees;
    }
};

/**
 * This is a base class for nodes that will do geometry queries.
 */
class geometry_node {
  protected:
    // std::vector< boost::shared_ptr< frantic::geometry::trimesh3 > > m_meshes;
    // std::vector< boost::shared_ptr< frantic::geometry::trimesh3_kdtree > > m_trees;

    int m_geomSrcIndex;

  public:
    geometry_node( int geomSrcIndex )
        : m_geomSrcIndex( geomSrcIndex ) {}

    // void add_mesh( boost::shared_ptr< frantic::geometry::trimesh3 > pMesh, boost::shared_ptr<
    // frantic::geometry::trimesh3_kdtree > pTree  );

    int get_geometry_input_index() const { return m_geomSrcIndex; }

    // const std::vector< boost::shared_ptr< frantic::geometry::trimesh3 > >& get_meshes() const {
    //	return m_meshes;
    // }
};

/**
 * This AST node will perfrom a ray-intersection on a set of input geometries
 */
class ray_intersect_node : public operation_node<2>, public geometry_node {
  public:
    ray_intersect_node( int id, int geomSrcIndex, int src1Index, int src2Index );
    virtual ~ray_intersect_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

/**
 * This AST node will find the nearest point on a set of input geometries
 */
class nearest_point_node : public operation_node<1>, public geometry_node {
  public:
    nearest_point_node( int id, int geomSrcIndex, int src1Index );
    virtual ~nearest_point_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

class surf_data_value_node : public operation_node<1> {
  public:
    enum target_result {
        RESULT_INVALID = 0,
        RESULT_POSITION,
        RESULT_SIGNED_DIST,
        RESULT_FACE_NORMAL,
        RESULT_SMOOTH_NORMAL,
        RESULT_MESH_INDEX,
        RESULT_FACE_INDEX,
        RESULT_FACE_MAT_ID,
        RESULT_SMOOTH_GROUP,
        RESULT_BARY_COORDS,
        RESULT_TEXTURE_COORD,
        RESULT_INTERSECTED,

        NUM_TARGET_RESULTS
    };

  private:
    target_result m_target;

  public:
    surf_data_value_node( int id, int src1Index, target_result target );
    virtual ~surf_data_value_node();
    virtual void compile( const std::vector<channel_op_node*>& expressionTree,
                          channel_operation_compiler& inoutCompData );
};

//****************************************
//  Template class member implementations.
//****************************************
template <class T, int Arity>
input_value_op_node<T, Arity>::input_value_op_node( int id, T value[] )
    : channel_op_node( id ) {
    memcpy( m_value, value, sizeof( T[Arity] ) );
}

template <class T, int Arity>
input_value_op_node<T, Arity>::~input_value_op_node() {}

template <class T, int Arity>
struct input_value_op_node_impl {
    T m_value[Arity];
    int m_destOffset;

    void init( int destOffset, T* value ) {
        m_destOffset = destOffset;
        memcpy( m_value, value, sizeof( T ) * Arity );
    }

    static void eval( char*, char* stack, void* pVoidData ) {
        input_value_op_node_impl* pData = reinterpret_cast<input_value_op_node_impl*>( pVoidData );
        T* pDest = reinterpret_cast<T*>( stack + pData->m_destOffset );
        memcpy( pDest, pData->m_value, sizeof( T ) * Arity );
    }
};

template <class T, int Arity>
void input_value_op_node<T, Arity>::compile( const std::vector<channel_op_node*>& /*expressionTree*/,
                                             channel_operation_compiler& inoutCompData ) {
    // Prevent this node from being compiled twice. TODO: Perhaps that is accpetable for inputs?
    temporary_result* r = inoutCompData.get_node_results( m_nodeId );
    if( r != NULL )
        throw channel_compiler_error( m_nodeId, "Graph cycle detected" );

    data_type_t valueType = frantic::channels::channel_data_type_traits<T>::data_type();
    if( valueType != data_type_float32 && valueType != data_type_int32 )
        throw channel_compiler_error( m_nodeId, "Unexpected input type: " + detail::data_type_string( valueType ) );

    r = inoutCompData.allocate_temporary( m_nodeId, valueType, Arity );

    detail::code_segment theSeg;
    theSeg.codePtr = &input_value_op_node_impl<T, Arity>::eval;
    theSeg.codeData = malloc( sizeof( input_value_op_node_impl<T, Arity> ) );
    reinterpret_cast<input_value_op_node_impl<T, Arity>*>( theSeg.codeData )->init( r->offset, m_value );

    inoutCompData.append_code_segment( theSeg );
}

} // namespace channels
} // namespace frantic
