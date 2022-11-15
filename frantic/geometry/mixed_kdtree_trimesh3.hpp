// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <frantic/geometry/mixed_kdtree.hpp>
#include <frantic/geometry/trimesh3.hpp>

namespace frantic {
namespace geometry {

namespace detail {

/**
 *  Return the Velocity channel for the specified vertex as a vector3f.
 * You must ensure that the mesh exists while you are using this object.
 * Should there be a convert accessor for trimesh3 channels, or should I
 * just special-case the float16 velocity channel ?  float32 and float16
 * are the only two I've seen.
 */
class trimesh3_vertex_velocity_accessor {
    // const frantic::geometry::trimesh3 & m_mesh;
    frantic::geometry::const_trimesh3_vertex_channel_general_accessor m_ca;
    frantic::channels::channel_type_convertor_function_t m_convert;
    bool m_hasVelocity;
    bool m_hasConvert;

  public:
    trimesh3_vertex_velocity_accessor( const frantic::geometry::trimesh3& mesh );
    frantic::graphics::vector3f get( const std::size_t vertexNumber ) const;
    bool has_velocity( void ) const { return m_hasVelocity; };
};

class set_channel_from_trimesh3_vertex_channel {
    // convert_and_set_channel m_setter;
    frantic::channels::channel_weighted_sum_combine_function_t m_setFromWeightedSum;
    const_trimesh3_vertex_channel_general_accessor m_inputChannelAccessor;
    std::size_t m_outputOffset;

  public:
    set_channel_from_trimesh3_vertex_channel( channels::channel& outChannel, const frantic::geometry::trimesh3& mesh,
                                              const frantic::tstring& channelName );
    static bool is_valid( channels::channel& outChannel, const trimesh3& mesh, const frantic::tstring& channelName );
    void set( char* outData, const std::size_t faceNumber, const vector3f& barycentricCoord ) const;
};

class set_channel_from_trimesh3_face_channel {
    convert_and_set_channel m_setter;
    const_trimesh3_face_channel_general_accessor m_channelAccessor;

  public:
    set_channel_from_trimesh3_face_channel( channels::channel& outChannel, const trimesh3& mesh,
                                            const frantic::tstring& channelName );
    static bool is_valid( channels::channel& outChannel, const trimesh3& mesh, const frantic::tstring& channelName );
    void set( char* outData, const std::size_t faceNumber ) const;
};

/**
 *  Set channels named "BarycentricCoord", "FaceNumber", "FaceIndex",
 * "GeometricNormal",
 * or with the same name as a vertex or face channel in the mesh.
 */
class set_channel_map_data_from_trimesh3 {
    boost::shared_ptr<frantic::geometry::trimesh3> m_mesh;

    std::vector<set_channel_from_trimesh3_vertex_channel> m_setFromVertexChannel;
    std::vector<set_channel_from_trimesh3_face_channel> m_setFromFaceChannel;

    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_setBarycentricCoord;
    frantic::channels::channel_cvt_accessor<boost::int32_t> m_setFaceIndex;
    frantic::channels::channel_cvt_accessor<boost::int32_t> m_setFaceNumber;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_setGeometricNormal;

    void reset( void );

  public:
    set_channel_map_data_from_trimesh3( boost::shared_ptr<frantic::geometry::trimesh3> mesh );

    void set_channel_map( const frantic::channels::channel_map& channelMap );
    void set_channel_map_data( char* outData, const boost::int32_t faceNumber,
                               const frantic::graphics::vector3f& barycentricCoord,
                               const frantic::graphics::vector3f& geometricNormal ) const;
};

} // namespace detail

class mixed_kdtree_trimesh3_primitive;

class mixed_kdtree_trimesh3 : public mixed_kdtree_primitive_creator, public output_channel_map_listener {
    boost::shared_ptr<frantic::geometry::trimesh3> m_mesh;
    // boost::shared_ptr<trimesh3> m_mesh;

    // todo: probably no need to hold on to this if the kdtree does
    // std::vector<boost::shared_ptr<mixed_kdtree_trimesh3_primitive> > m_primitives;

    // output channel accessors
    detail::set_channel_map_data_from_trimesh3 m_setChannelMapDataFromMesh;
    /*
    std::vector<detail::set_channel_from_trimesh3_vertex_channel> m_setFromVertexChannel;
    std::vector<detail::set_channel_from_trimesh3_face_channel> m_setFromFaceChannel;

    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_setBarycentricCoord;
    frantic::channels::channel_cvt_accessor<boost::int32_t> m_setFaceIndex;
    frantic::channels::channel_cvt_accessor<boost::int32_t> m_setFaceNumber;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_setGeometricNormal;
    */
    // detail::set_channel_from_mem_function1<frantic::graphics::vector3f, mixed_kdtree_trimesh3_primitive, const char
    // *> m_setBarycentricCoord; detail::set_channel_from_mem_function<boost::int32_t, mixed_kdtree_trimesh3_primitive>
    // m_setFaceIndex; detail::set_channel_from_mem_function<boost::int32_t, mixed_kdtree_trimesh3_primitive>
    // m_setFaceNumber; detail::set_channel_from_mem_function1<frantic::graphics::vector3f,
    // mixed_kdtree_trimesh3_primitive, const char *> m_setGeometricNormal;

    // void reset_channel_accessors( void );

    // not implemented
    mixed_kdtree_trimesh3& operator=( const mixed_kdtree_trimesh3& );

    // void initialize_from_mesh( boost::shared_ptr<frantic::geometry::trimesh3> mesh );

    // mixed_kdtree_trimesh3( boost::shared_ptr<frantic::geometry::trimesh3> mesh, bool initFromMesh );

  public:
    mixed_kdtree_trimesh3( boost::shared_ptr<frantic::geometry::trimesh3> mesh );
    ~mixed_kdtree_trimesh3( void );

    // static boost::shared_ptr<mixed_kdtree_trimesh3> from_face_subset( boost::shared_ptr<frantic::geometry::trimesh3>
    // mesh, const std::vector<std::size_t> & faceNumbers );

    std::size_t get_primitive_count( void );
    void get_primitives( mixed_kdtree_primitive_recorder* recorder );
    trimesh3& get_mesh_ref( void );
    void to_channel_map_data( char* data, const mixed_kdtree_trimesh3_primitive* mixed,
                              const mixed_kdtree_point_data* pointData ) const;
    void set_output_channel_map( const frantic::channels::channel_map& channelMap );
};

class mixed_kdtree_trimesh3_constant_velocity_primitive;

class mixed_kdtree_trimesh3_constant_velocity : public mixed_kdtree_primitive_creator,
                                                public output_channel_map_listener {
    boost::shared_ptr<frantic::geometry::trimesh3> m_mesh;

    // todo: probably no need to hold on to this if the kdtree does
    // std::vector<boost::shared_ptr<mixed_kdtree_primitive> > m_primitives;

    boost::shared_ptr<mixed_kdtree_trimesh3> m_staticGeometry;

    double m_t0, m_t1;

    // output channel accessors
    detail::set_channel_map_data_from_trimesh3 m_setChannelMapDataFromMesh;
    /*
    std::vector<detail::set_channel_from_trimesh3_vertex_channel> m_setFromVertexChannel;
    std::vector<detail::set_channel_from_trimesh3_face_channel> m_setFromFaceChannel;

    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_setBarycentricCoord;
    frantic::channels::channel_cvt_accessor<boost::int32_t> m_setFaceIndex;
    frantic::channels::channel_cvt_accessor<boost::int32_t> m_setFaceNumber;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_setGeometricNormal;
    */
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_setDisplacementDuringTimeStep;

    // detail::set_channel_from_mem_function1<frantic::graphics::vector3f,
    // mixed_kdtree_trimesh3_constant_velocity_primitive, const char *> m_setBarycentricCoord;
    // detail::set_channel_from_mem_function<boost::int32_t, mixed_kdtree_trimesh3_constant_velocity_primitive>
    // m_setFaceIndex; detail::set_channel_from_mem_function<boost::int32_t,
    // mixed_kdtree_trimesh3_constant_velocity_primitive> m_setFaceNumber;
    // detail::set_channel_from_mem_function1<frantic::graphics::vector3f,
    // mixed_kdtree_trimesh3_constant_velocity_primitive, const mixed_kdtree_point_data *> m_setGeometricNormal;

    // detail::set_channel_from_mem_function1<vector3f, mixed_kdtree_trimesh3_constant_velocity_primitive, const char *>
    // m_setDisplacementDuringTimeStep;

    // initial position, at the beginning of the time step
    std::vector<frantic::graphics::vector3f> m_vertexInitialPosition;

    // displacement during the time step
    std::vector<frantic::graphics::vector3f> m_vertexDisplacement;

    void init_vertex_position_and_displacement( const frantic::geometry::trimesh3& mesh, const double t0,
                                                const double t1 );
    void reset_channel_accessors( void );

    // not implemented
    mixed_kdtree_trimesh3_constant_velocity& operator=( const mixed_kdtree_trimesh3_constant_velocity& );

  public:
    /**
     *  Insert a mesh with moving vertices.  Vertices are assumed to have
     * constant velocity over the time step.  Only the velocity channel of
     * the mesh is used to form the vertex velocities.
     *
     * @param mesh the mesh whose faces will be inserted into the tree.
     * @param t0 the start of the time step during which the vertices are
     *		moving.  This is time 0 from the perspective of subsequent
     *		ray intersections.
     * @param t1 the end of the timestep during which the vertices are
     *		moving.
     *
     */
    mixed_kdtree_trimesh3_constant_velocity( boost::shared_ptr<trimesh3> mesh, const double t0, const double t1 );

    const trimesh3& get_mesh_ref( void );

    std::size_t get_primitive_count( void );

    void get_primitives( mixed_kdtree_primitive_recorder* recorder );

    vector3f get_vertex_initial_position( const boost::int32_t vertexNumber ) const;

    // todo: maybe change displacement into PositionStart and PositionEnd or
    // something similar
    vector3f get_vertex_displacement( const boost::int32_t vertexNumber ) const;

    // this is nearly shared with the trimesh3 manager
    // it should probably be refactored
    void to_channel_map_data( char* data, const mixed_kdtree_trimesh3_constant_velocity_primitive* mixed,
                              const mixed_kdtree_point_data* pointData ) const;
    void set_output_channel_map( const channels::channel_map& channelMap );
};

class mixed_kdtree_trimesh3_point_data {
  public:
    const frantic::geometry::trimesh3* mesh;
    std::size_t faceNumber;
    frantic::graphics::vector3f barycentricCoords;

    mixed_kdtree_trimesh3_point_data( const frantic::geometry::trimesh3* mesh, const std::size_t faceNumber,
                                      const frantic::graphics::vector3f& barycentricCoords )
        : mesh( mesh )
        , faceNumber( faceNumber )
        , barycentricCoords( barycentricCoords ) {}
};

class mixed_kdtree_trimesh3_constant_velocity_primitive : public mixed_kdtree_primitive {
    mixed_kdtree_trimesh3_constant_velocity& m_helper;

    const frantic::graphics::vector3 m_face;
    const boost::int32_t m_faceNumber;

    // initial vertex positions
    // position at the beginning of the time step
    /*
    frantic::graphics::vector3f m_x1;
    frantic::graphics::vector3f m_x2;
    frantic::graphics::vector3f m_x3;
    */

    // vertex velocity
    // expects [tMin,tMax] = [0,1], so scale accordingly !

    // So velocity here is misleading --
    // this is really the displacement over one timestep.

    // The timestep is [t0,t1] in the constructor.

    // TODO: there's a lot of redundancy here --
    // this storage should be shared between faces.
    // the vertex positions should be shared too....

    // if the t0 is always 0, then the vertices can
    // simply be the mesh vertices.
    // otherwise we can use shared storage for them
    // somewhere (maybe a member of the mixed_kdtree ?)

    // vertex displacement during the time step
    /*
    frantic::graphics::vector3f m_v1;
    frantic::graphics::vector3f m_v2;
    frantic::graphics::vector3f m_v3;
    */

    // not implemented
    mixed_kdtree_trimesh3_constant_velocity_primitive&
    operator=( const mixed_kdtree_trimesh3_constant_velocity_primitive& );
    /**
     *  Test whether a point is inside the given triangle.  Also computes
     * the point's barycentric coordinates if it is inside the triangle.
     *
     * @note: right now it computes the point's barycentric coordindates
     *		always but this will likely change (ie: the
     *      outBarycentricCoord will be valid only when the function
     *      returns true).
     * @note: this needn't be a member function.  It's here for now in
     *		case some precomputed data gets stored in this object.
     *
     * precondition: pt should be coplanar with the triangle
     * otherwise the test is probably not useful for you
     *
     * @param pt the point to look for in the triangle.
     * @param va the triangle's first vertex
     * @param vb the triangle's second vertex
     * @param vc the triangle's third vertex
     * @param[out] outBarycentricCoord if the function returns true,
     *		then this is the barycentric coordinate of pt in the triangle.
     * @return true if pt is inside the triangle, and false otherwise.
     */
    bool is_coplanar_point_in_triangle( const frantic::graphics::vector3f& pt, const frantic::graphics::vector3f& va,
                                        const frantic::graphics::vector3f& vb, const frantic::graphics::vector3f& vc,
                                        frantic::graphics::vector3f& outBarycentricCoord );
    frantic::graphics::vector3f get_x1( void ) const;
    frantic::graphics::vector3f get_x2( void ) const;
    frantic::graphics::vector3f get_x3( void ) const;
    void get_x( frantic::graphics::vector3f& x1, frantic::graphics::vector3f& x2, frantic::graphics::vector3f& x3 );

    frantic::graphics::vector3f get_v1( void ) const;
    frantic::graphics::vector3f get_v2( void ) const;
    frantic::graphics::vector3f get_v3( void ) const;
    void get_v( frantic::graphics::vector3f& v1, frantic::graphics::vector3f& v2, frantic::graphics::vector3f& v3 );

  public:
    // todo: store t0 and t1, and have the object disappear outside of these times?
    // or clamp the motion to within these times ?
    mixed_kdtree_trimesh3_constant_velocity_primitive( mixed_kdtree_trimesh3_constant_velocity& helper,
                                                       const boost::int32_t faceNumber, const double /*t0*/,
                                                       const double /*t1*/ );

    /**
     * assume the ray travels from origin at time 0
     * to origin + direction at time 1
     */
    bool intersect_particle( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                             mixed_kdtree_ray_observer* observer );
    // test for intersection with an instantaneous ray
    // the geometry is offset to the specified time
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                        mixed_kdtree_ray_observer* observer, double time );
    // todo: should the output be restricted to nodebounds ?
    bool find_nearest_point( const frantic::graphics::vector3f& point, boundbox3f& /*nodeBounds*/,
                             mixed_kdtree_distance_observer* observer, double time );
    void get_any_point( vector3f& outPoint, vector3f& outNormal, double time );

    frantic::graphics::boundbox3f get_bounds() const;
    frantic::graphics::boundbox3f intersect_with( const boundbox3f& box ) const;

    std::size_t get_data_size( void ) const;

    boost::int32_t get_face_number( void ) const;
    frantic::graphics::vector3f get_geometric_normal( const mixed_kdtree_point_data* pointData ) const;
    frantic::graphics::vector3f get_barycentric_coord( const char* data ) const;
    frantic::graphics::vector3f get_motion_during_time_step( const char* data ) const;

    void to_raytrace_intersection( frantic::geometry::raytrace_intersection& raytraceIntersection,
                                   const char* data ) const;
    void to_nearest_point_search_result( frantic::geometry::nearest_point_search_result& searchResult,
                                         const char* data ) const;
    void to_channel_map_data( char* outData, const mixed_kdtree_point_data* pointData ) const;
};

class mixed_kdtree_trimesh3_primitive : public mixed_kdtree_primitive {
    mixed_kdtree_trimesh3& m_helper;
    const boost::uint32_t m_faceNumber;
    plane3f m_facePlane;

    char m_barycentric0Axis;
    char m_barycentric1Axis;
    float m_barycentricInverseDeterminant;
    vector3 m_face;
    // vector3f m_vertices[3];

    // float m_pluckerA[6];
    // float m_pluckerB[6];
    // float m_pluckerC[6];
    /*
    int	m_axis1;
    int	m_axis2;
    float	m_nu;
    float	m_nv;
    float	m_nd;

    float	m_bnu;
    float	m_bnb;
    float	m_bd;

    float	m_cnu;
    float	m_cnb;
    float	m_cd;
*/
    // not implemented
    mixed_kdtree_trimesh3_primitive& operator=( mixed_kdtree_trimesh3_primitive& );

    frantic::graphics::vector3f compute_barycentric_coordinates( const vector3f& pt );
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                        mixed_kdtree_ray_observer* observer, bool fixedTime, double time );

  public:
    mixed_kdtree_trimesh3_primitive( mixed_kdtree_trimesh3& helper, const boost::uint32_t faceNumber );
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                        mixed_kdtree_ray_observer* observer, double time );
    bool intersect_particle( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                             mixed_kdtree_ray_observer* observer );
    // todo: should the output be restricted to nodeBounds ?
    bool find_nearest_point( const frantic::graphics::vector3f& point, frantic::graphics::boundbox3f& /*nodeBounds*/,
                             mixed_kdtree_distance_observer* distanceObserver, double time );
    void get_any_point( frantic::graphics::vector3f& outPoint, frantic::graphics::vector3f& outNormal, double time );

    frantic::graphics::boundbox3f get_bounds() const;

    frantic::graphics::boundbox3f intersect_with( const boundbox3f& box ) const;

    // I'm not sure if this makes sense, but I am implementing it for testing purposes
    boost::int32_t get_face_number( void ) const;

    frantic::graphics::vector3f get_barycentric_coord( const char* data ) const;

    frantic::graphics::vector3f get_geometric_normal( const char* ) const;

    std::size_t get_data_size( void ) const;

    void to_raytrace_intersection( frantic::geometry::raytrace_intersection& raytraceIntersection,
                                   const char* data ) const;
    void to_nearest_point_search_result( frantic::geometry::nearest_point_search_result& searchResult,
                                         const char* data ) const;
    void to_channel_map_data( char* outData, const mixed_kdtree_point_data* pointData ) const;
};

bool is_intersecting( const frantic::graphics::boundbox3f& bbox, const frantic::graphics::ray3f& ray, const double t0,
                      const double t1 );

/**
 *  Determine whether the ray intersects the triangle within the available
 * precision.  If an intersection is found, then this function produces
 * a bounding box which contains the intersection.
 *
 * @param ray the ray to check for intersections.
 * @param t0 the first value of the ray parameter.
 * @param t1 the last value of the ray parameter.
 * @param va the first vertex of the triangle to find intersections with.
 * @param vb the second vertex.
 * @param vc the third vertex.
 * @param outBBox if the function returns true, then this boundbox contains
 *		the intersection between the ray and the triangle.
 * @param lastDiagonalLength this is set during recursion; you should simply
 *		use the default value.
 * @return true if the ray intersects the triangle within the available
 *		precision, and false otherwise.
 *
 * Recursive algorithm for detecting ray-volume intersections as described
 * in:
 *  Dammertz and Keller.  "Improving Ray Tracing Precision by Object Space
 * Intersection Computation."
 * Be aware that details in my implementation may be off.
 */
bool intersect_triangle( const frantic::graphics::ray3f& ray, const double t0, const double t1,
                         const frantic::graphics::vector3f& va, const frantic::graphics::vector3f& vb,
                         const frantic::graphics::vector3f& vc, frantic::graphics::boundbox3f& outBBox,
                         const double lastDiagonalLength = std::numeric_limits<double>::max() );

} // namespace geometry
} // namespace frantic
