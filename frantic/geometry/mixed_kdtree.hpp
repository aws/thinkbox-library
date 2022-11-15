// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/geometry/mixed_kdtree_detail.hpp>
#include <frantic/geometry/raytracing.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/graphics/ray3f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/particles/particle_array.hpp>

// Possible alternative to std::vector for internal intersection
// result data.
//#include <boost/pool/pool.hpp>
//#include <boost/pool/pool_alloc.hpp>

#include <array>
#include <fstream>
#include <iomanip>

// This doesn't play nicely with krakatoa/3ds max.
// Presumably this is because krakatoa uses namespace boost --
// I think foreach adds END to the namespace which interferes with 3ds max's
// includes.
//#include <boost/foreach.hpp>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace frantic {
namespace geometry {

class mixed_kdtree_primitive_creator;
class mixed_kdtree_primitive;
class mixed_kdtree_distance_observer;
class mixed_kdtree_ray_observer;
class mixed_kdtree_node;
class mixed_kdtree_point_data;
class mixed_kdtree_ray_intersection;
class mixed_kdtree_nearest_point;
class output_channel_map_listener;

//
// mixed_kdtree_node
//

class mixed_kdtree;

void build_kdtree_greedy_SAH_nlogn( mixed_kdtree_node& node,
                                    const std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& primitives,
                                    const boundbox3f& bounds, std::vector<mixed_kdtree_detail::kdtreeEvent>& events,
                                    std::vector<mixed_kdtree_detail::index_t>& indices,
                                    std::vector<boost::int8_t>& objFlags, const int maxDepth, int depth );
void build_kdtree_greedy_SAH_nlogn( mixed_kdtree_node& node,
                                    std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& primitives,
                                    const frantic::graphics::boundbox3f& bounds );

inline void get_distance_interval_for_ray_splitting_plane_intersection( const float rayOriginAlongAxis,
                                                                        const float rayDirectionAlongAxis,
                                                                        const float split, float& outMinDistance,
                                                                        float& outMaxDistance );

class mixed_kdtree_node {
    friend class mixed_kdtree;
    friend void build_kdtree_greedy_SAH_nlogn( mixed_kdtree_node& node,
                                               const std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& primitives,
                                               const boundbox3f& bounds,
                                               std::vector<mixed_kdtree_detail::kdtreeEvent>& events,
                                               std::vector<mixed_kdtree_detail::index_t>& indices,
                                               std::vector<boost::int8_t>& objFlags, const int maxDepth, int depth );

    union {
        mixed_kdtree_node* m_children;
        std::vector<mixed_kdtree_primitive*>* m_primitivePointers;
    };

    int m_axisAndLeafFlag;
    float m_split;

    mixed_kdtree_node();
    ~mixed_kdtree_node();

    // Constructor for making an interior node
    // NOTE: Only one initialize function should be called once on a node.  No checks are done to make sure things are
    // valid.
    // void initialize( mixed_kdtree_node<mixed_kdtree_primitive>* children, int axis, float split )
    void initialize( mixed_kdtree_node* children, int axis, float split );
    void initialize( std::vector<mixed_kdtree_detail::index_t>& primitiveIndices,
                     const std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& );

    mixed_kdtree_node* left_child();
    mixed_kdtree_node* right_child();

    const mixed_kdtree_node* left_child() const;
    const mixed_kdtree_node* right_child() const;
    int get_node_count() const;
    int get_maximum_depth() const;
    int get_largest_leaf_size() const;
    std::size_t get_populated_leaf_count() const;
    std::size_t get_leaf_count() const;
    void dump_tree( std::ostream& out, int depth ) const;

    // This returns the nearest intersection between the ray and the mesh, or false if none is found.
    template <class ray_observer, class ray_intersection_type>
    inline void traverse_ray( const frantic::graphics::ray3f& ray, const double tMin, double tMax,
                              ray_observer& observer, const ray_intersection_type& rayCollisionTest );

    template <class distance_observer>
    inline void traverse_nearest_points( const frantic::graphics::vector3f& point,
                                         frantic::graphics::boundbox3f& nodeBounds, distance_observer& observer,
                                         double time );
};

class mixed_kdtree {
    mixed_kdtree_node* m_rootNode;
    std::vector<boost::shared_ptr<mixed_kdtree_primitive_creator>> m_objects;
    std::vector<output_channel_map_listener*> m_outputChannelMapListeners;
    // The primitives are owned by the m_objects.
    // Now that I am holding them here too,
    // maybe they should be owned by the kdtree instead..
    std::vector<boost::shared_ptr<mixed_kdtree_primitive>> m_primitives;
    frantic::graphics::boundbox3f m_rootBounds;
    bool m_final;
    std::size_t m_maxDataSize;

    frantic::channels::channel_map m_channelMap;
    std::vector<char> m_defaultChannelMapData;

    detail::common_named_channel_setters m_commonNamedChannelSetters;

    /**
     *  Write the defaults values for all channels to the specified
     * buffer.
     *
     * @note the buffer must have at least the size of the kd-tree's
     *		channel map, which may be accessed by calling get_channel_map().
     *
     * @param[out] data the buffer to write the default channel values to.
     */
    void write_default_channel_map_data( char* data );

    /**
     *  Return an object with setters for commonly used channels,
     * such as set_position().
     *
     * @return an object with setters for commonly-used channels.
     */
    const detail::common_named_channel_setters& get_common_named_channel_setters( void );

    // std::filebuf m_logStream;
    // std::ofstream m_outFile;

    /**
     *  Add the kd-tree's default channels to the channel map.  For now
     * the default channels are Position, Distance and Normal.
     *
     * @param channelMap the channel map to add the kdtree's default
     *		channels to.
     */
    void insert_default_channels( frantic::channels::channel_map& channelMap );
    /**
     *  Construct the kd-tree's default channel map.
     *
     *  Reset the channel map, add the kd-tree's default channels, and
     * finalize the channel map.
     *
     * @param channelMap the channel map which will hold the kdtree's
     *		default channel map.
     */
    void build_default_channel_map( frantic::channels::channel_map& channelMap );

    /**
     *  Initialize the kd-tree's member variables using their default values.
     */
    void initialize( void );

    /**
     *  Add a geometry object to the kd-tree.
     *
     *  This cannot be done after calling finalize().
     */
    void insert( boost::shared_ptr<mixed_kdtree_primitive_creator> object );

    // only tested indirectly
    // I'm not sure if the mixed_kdtree_ray_intersection objects should be
    // available to callers.
    /**
     *  Find all intersections between the specified ray and the kdtree's
     * geometry.
     *
     * @param ray the ray to cast into the kd-tree.  The ray will be followed
     *		from ray.at( tMin ) through ray.at( tMax ).
     * @param tMin the ray parameter for the start of the ray.
     * @param tMax the ray parameter for the end of the ray.
     * @param[out] outIntersections an array of information about
     *		the intersections between the ray and the kdtree's geometry.
     *		This is cleared before adding intersections.
     * @param time the time at which to cast the ray.
     */
    void intersect_ray_all( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                            std::vector<mixed_kdtree_ray_intersection>& outIntersections, double time = 0 ) const;

  public:
    /**
     *  Create an empty kd-tree with a default channel map.
     *
     *  You should insert at least one object and call the
     * mixed_kdtree's finalize() method before performing
     * and searches using the tree.
     */
    mixed_kdtree();
    /**
     *  Create an empty kd-tree with the specified channel map.
     */
    mixed_kdtree( const frantic::channels::channel_map& cm );
    /**
     *  Create and finalize a kd-tree which contains the specified mesh.
     *
     * @param mesh the mesh to insert into the kd-tree.
     */
    mixed_kdtree( boost::shared_ptr<frantic::geometry::trimesh3> mesh );

    ~mixed_kdtree();

    /**
     *  Set the channels to write for searches which result in a
     * char * or a particle_array.
     *
     * @note This resets the default channel map data specified by
     *		set_default_channel_map_data.
     *
     * @param channelMap the channels to output for searches which
     *		result in a char* or a particle_array.
     */
    void set_channel_map( const frantic::channels::channel_map& channelMap );

    /**
     *  Return the channel map which is used by searches which
     * result in a char* or a particle_array.
     */
    const frantic::channels::channel_map& get_channel_map( void ) const;

    /**
     *  Set default values for the output channel map data.
     * This default data is reset by caling set_channel_map, so this function
     * should be called after set_channel_map.
     *
     * @param defaultData the default data to write to output channel map
     *		data.
     */
    void set_default_channel_map_data( const std::vector<char>& defaultData );

    /**
     *  Add a trimesh3 to the kdtree.
     *
     * @param mesh the mesh to add to the kdtree.
     */
    void insert( boost::shared_ptr<frantic::geometry::trimesh3> mesh );

    /**
     *  Add a trimesh3 to the kdtree.  The kdtree will be animated
     * according to its Velocity channel.
     *
     * When performing a particle intersection with such animated
     * geometry, the ray parameter t is measured as a fraction of the
     * time step, 0 = tMin < t < 1 = tMax.  Using tMin, tMax outside
     * of the range [0,1] will probably not behave how you want.
     *
     * Let x and v denote the position and velocity, respectively, of
     * a vertex in the mesh.  At t = 0, the vertex will be located
     * at x + t0 * v, where t0 is the starting time passed to the function.
     * At time t = 1 the vertex will be located at x + t1 * v.
     *
     * @param mesh the mesh to add to the kdtree.
     * @param t0 time offset for the start of the mesh's animation.
     *		This time corresponds to t = 0 during intersect_particle
     *		searches.
     * @param t1 time offset for the end of the mesh's animation.
     *		This time corresponds to t = 1 during intersect_particle
     *		searches.
     */
    void insert_with_constant_velocity( boost::shared_ptr<frantic::geometry::trimesh3> mesh, double t0, double t1 );

    /**
     *  Finalize the tree so that you can perform searches in it.
     * You cannot insert new objects into the tree after finalizing
     * it, and you should not modify objects while they are used in
     * the tree.
     */
    void finalize( void );

    /**
     *  Return a boundbox that contains all of the geometry in the
     * kdtree.  This cannot be called before finalize().
     *
     * @return a boundbox that contains all the geometry in the kdtree.
     */
    boundbox3f get_bounds() const;

    /**
     * @return the total number of nodes in the tree.
     */
    std::size_t get_node_count() const;

    /**
     * @return the maximum depth of a node in the tree.
     */
    std::size_t get_maximum_depth() const;

    /**
     * @return the greatest number of primitives in a node in the tree.
     */
    std::size_t get_largest_leaf_size() const;

    /**
     * @return the total number of leaf nodes in the tree.
     */
    std::size_t get_leaf_count() const;

    /**
     * @return the total number of populated leaf nodes in the tree.
     */
    std::size_t get_populated_leaf_count() const;

    /**
     *  Output a representation of the tree, including all of its nodes,
     * to the specified stream.
     *
     * @param out the stream to write a representation of the tree into.
     */
    void dump_tree( std::ostream& out ) const;

    ///// Queries/////

    /**
     *  Find an intersection between a finite-velocity particle and
     * the geometry in the kd-tree.
     *
     * @note Usually [tMin, tMax] = [0, 1].
     *
     * @param ray a ray which represent the motion of a particle.  The
     *		particle moves from ray.at( tMin ) at time tMin through
     *		ray.at( tMax ) at time tMax.
     * @param tMin the starting time of the particle's motion.
     * @param tMax the ending time of the particle's motion.
     * @param[out] outIntersection if the function returns true, then
     *		this is information about the intersection between
     *		the particle and the geometry.
     * @return true if the particle intersected with the geometry in
     *		the kd-tree, and false otherwise.
     */
    bool intersect_particle( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                             frantic::geometry::raytrace_intersection& outIntersection ) const;

    /**
     *  Find an intersection between a finite-velocity particle and
     * the geometry in the kd-tree.
     *
     * @note Usually [tMin, tMax] = [0, 1].
     *
     * @param ray a ray which represent the motion of a particle.  The
     *		particle moves from ray.at( tMin ) at time tMin through
     *		ray.at( tMax ) at time tMax.
     * @param tMin the starting time of the particle's motion.
     * @param tMax the ending time of the particle's motion.
     * @param[out] outData if the function returns true, then
     *		this is channel-mapped information about the intersection
     *		between the particle and the geometry.  The channel mapped
     *		data follows the form of the kd-tree's channel map which
     *		may be accessed using the get_channel_map() method.  The
     *		buffer passed to outData must have a size large enough
     *		to accomodate the channel map's structure size.
     * @return true if the particle intersected with the geometry in
     *		the kd-tree, and false otherwise.
     */
    bool intersect_particle( const frantic::graphics::ray3f& ray, double tMin, double tMax, char* outData ) const;

    /**
     * Returns the nearest intersection between the ray and one of the
     * primitives, or false if none is found.
     *
     * @note you should use tMin >= 0 and tMax > tMin - I wouldn't trust this
     *		to work correctly otherwise.
     * @todo It may be nice to get rid of the restriction on tMin as a way to
     *		add epsilons.
     *
     * @param ray the ray to cast into the kd-tree.  The ray will be followed
     *		from ray.at( tMin ) through ray.at( tMax ).
     * @param tMin the ray parameter for the start of the ray.
     * @param tMax the ray parameter for the end of the ray.
     * @param[out] outIntersection if the function returns true, then this is
     *		information about the intersection point.
     * @param time the time at which to cast the ray.
     * @return true if the ray intersected with an object and false otherwise.
     */
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                        frantic::geometry::raytrace_intersection& outIntersection, double time = 0 ) const;

    /**
     * Returns the nearest intersection between the ray and one of the
     * primitives, or false if none is found.
     *
     * @note you should use tMin >= 0 and tMax > tMin - I wouldn't trust this
     *		to work correctly otherwise.
     *
     * @param ray the ray to cast into the kd-tree.  The ray will be followed
     *		from ray.at( tMin ) through ray.at( tMax ).
     * @param tMin the ray parameter for the start of the ray.
     * @param tMax the ray parameter for the end of the ray.
     * @param[out] outData if the function returns true, then
     *		this is channel-mapped information about the intersection
     *		between the particle and the geometry.  The channel mapped
     *		data follows the form of the kd-tree's channel map which
     *		may be accessed using the get_channel_map() method.  The
     *		buffer passed to outData must have a size large enough
     *		to accomodate the channel map's structure size.
     * @param time the time at which to cast the ray.
     * @return true if the ray intersected with an object and false otherwise.
     */
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax, char* outData,
                        double time = 0 ) const;

    /**
     * Returns the nearest intersection between the ray and one of the
     * primitives, or false if none is found.
     *
     * @note you should use tMin >= 0 and tMax > tMin - I wouldn't trust this
     *		to work correctly otherwise.
     *
     * @param ray the ray to cast into the kd-tree.  The ray will be followed
     *		from ray.at( tMin ) through ray.at( tMax ).
     * @param tMin the ray parameter for the start of the ray.
     * @param tMax the ray parameter for the end of the ray.
     * @param[out] outData if the function returns true, then
     *		this is channel-mapped information about the intersection
     *		between the particle and the geometry.  The channel mapped
     *		data follows the form of the kd-tree's channel map which
     *		may be accessed using the get_channel_map() method.  The
     *		output array will be resized to match the channel map size
     *		if necessary.
     * @param time the time at which to cast the ray.
     * @return true if the ray intersected with an object and false otherwise.
     */
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax, std::vector<char>& outData,
                        double time = 0 ) const;

    /**
     * Returns true if there is an intersection between the ray and one
     * of the primitives, or false if none is found.
     *
     * @note you should use tMin >= 0 and tMax > tMin - I wouldn't trust this
     *		to work correctly otherwise.
     *
     * @param ray the ray to cast into the kd-tree.  The ray will be followed
     *		from ray.at( tMin ) through ray.at( tMax ).
     * @param tMin the ray parameter for the start of the ray.
     * @param tMax the ray parameter for the end of the ray.
     * @param time the time at which to cast the ray.
     * @return true if the ray intersected with an object and false otherwise.
     */
    bool intersects_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax, double time = 0 ) const;

    /**
     * Returns true if there is an intersection between the ray and one
     * of the primitives, or false if none is found.
     *
     * @note you should use tMin >= 0 and tMax > tMin - I wouldn't trust this
     *		to work correctly otherwise.
     *
     * @param ray the ray to cast into the kd-tree.  The ray will be followed
     *		from ray.at( tMin ) through ray.at( tMax ).
     * @param tMin the ray parameter for the start of the ray.
     * @param tMax the ray parameter for the end of the ray.
     * @param time the time at which to cast the ray.
     * @return true if the ray intersected with an object and false otherwise.
     */
    bool intersects_ray_segment( const frantic::graphics::ray3f& ray, double tMin, double tMax, double time = 0 ) const;

    /**
     * Returns true if there is an intersection between segment and one
     * of the primitives, or false if none is found.
     *
     * @param start the starting point of the segment to search along.
     * @param end the ending point of the segment to search along.
     * @param time the time at which to cast the ray.
     * @return true if the segment intersected with an object and false
     *		otherwise.
     */
    bool intersects_ray_segment( const frantic::graphics::vector3f& start, const frantic::graphics::vector3f& end,
                                 double time = 0 ) const;

    /**
     *  Adds all intersections found between the ray and the geometry, and
     * maintains the intersections in order of increasing distance.
     *
     * @note you should use tMin >= 0 and tMax > tMin - I wouldn't trust this
     *		to work correctly otherwise.
     *
     * @param ray the ray to cast into the kd-tree.  The ray will be followed
     *		from ray.at( tMin ) through ray.at( tMax ).
     * @param tMax the ray parameter for the end of the ray.
     * @param[out] outIntersections an array of information about
     *		the intersections between the ray and the kdtree's geometry.
     *		This is not cleared before adding intersections.  The output
     *		intersections are sorted in order of increasing distance.
     * @param time the time at which to cast the ray.
     */
    void intersect_ray_all( const frantic::graphics::ray3f& ray, double tMax,
                            std::vector<raytrace_intersection>& outIntersections, double time = 0 ) const;

    /**
     *  Find all intersections between the ray and the geometry in the
     * kd-tree.
     *
     * @note you should use tMin >= 0 and tMax > tMin - I wouldn't trust this
     *		to work correctly otherwise.
     *
     * @param ray the ray to cast into the kd-tree.  The ray will be followed
     *		from ray.at( tMin ) through ray.at( tMax ).
     * @param tMax the ray parameter for the end of the ray.
     * @param[out] outIntersections an array of channel mapped
     *		information about the intersections between the ray and the
     *		kdtree's geometry.
     *		The channel mapped data follows the form of the kd-tree's
     *		channel map, which may be accessed using the get_channel_map()
     *		method.  The particle_array must have been created using this
     *		same channel map.
     *		This array not cleared before adding intersections.
     * @param time the time at which to cast the ray.
     */

    void intersect_ray_all( const frantic::graphics::ray3f& ray, double tMax,
                            frantic::particles::particle_array& outIntersections, double time = 0 ) const;

    /**
     *  Find the nearest point in the kd-tree to the specified point,
     * within a specified search radius.
     *
     * @param point the point to search from.
     * @param maxDistance the maximum distance away from the specified point
     *		to search.
     * @param[out] outNearestPoint if the function returns true, then this is
     *		information regarding the nearest point in the kd-tree.
     * @param time the time at which to search.
     * @return true if a point could be found within the specified search
     *		radius, and false otherwise.
     */
    bool find_nearest_point( const frantic::graphics::vector3f& point, float maxDistance,
                             frantic::geometry::nearest_point_search_result& outNearestPoint, double time = 0 ) const;

    /**
     *  Find the nearest point in the kd-tree to the specified point,
     * within a specified search radius.
     *
     * @param point the point to search from.
     * @param maxDistance the maximum distance away from the specified point
     *		to search.
     * @param[out] outData if the function returns true, then this is
     *		information regarding the nearest point in the kd-tree.
     *		The data uses the kd-tree's channel map, which may be accessed
     *		using the get_channel_map() method.  The data array must have
     *		at least the size of the channel map's structure_size().
     * @param time the time at which to search.
     * @return true if a point could be found within the specified search
     *		radius, and false otherwise.
     */
    bool find_nearest_point( const frantic::graphics::vector3f& point, float maxDistance, char* outData,
                             double time = 0 ) const;

    /**
     *  Find the nearest point in the kd-tree to the specified point,
     * within a specified search radius.
     *
     * @param point the point to search from.
     * @param maxDistance the maximum distance away from the specified point
     *		to search.
     * @param[out] outData if the function returns true, then this is
     *		information regarding the nearest point in the kd-tree.
     *		The data uses the kd-tree's channel map, which may be accessed
     *		using the get_channel_map() method.  The array will be resized
     *		to match the channel map if necessary.
     * @param time the time at which to search.
     * @return true if a point could be found within the specified search
     *		radius, and false otherwise.
     */
    bool find_nearest_point( const frantic::graphics::vector3f& point, float maxDistance, std::vector<char>& outData,
                             double time = 0 ) const;

    /**
     *  Find the nearest point in the kd-tree to the specified point,
     * within a specified search radius.
     *
     * @param point the point to search from.
     * @param maxDistance the maximum distance away from the specified point
     *		to search.
     * @param[out] pOutNearestPoint a pointer to the variable which will
     *		hold the nearest point's position if the function returns true.
     *		This pointer can be NULL if you don't want the position.
     * @param[out] pOutDistance a pointer to the variable which will
     *		hold the nearest point's distance if the function returns true.
     *		This pointer can be NULL if you don't want the distance.
     * @param[out] pOutNormal a pointer to the variable which will
     *		hold the nearest point's normal if the function returns true.
     *		This pointer can be NULL if you don't want the normal.
     * @param time the time at which to search.
     * @return true if a point could be found within the specified search
     *		radius, and false otherwise.
     */
    bool find_nearest_point( const vector3f& point, float maxDistance, vector3f* pOutNearestPoint, float* pOutDistance,
                             vector3f* pOutNormal, double time = 0 ) const;

    /**
     *  Collect up to maxPrimitiveCount of the primitives closest to the
     * specified point.
     *
     *  The primitives are returned in order of increasing distance, and the
     * list is cleared before adding primitives.
     *
     * @param point the point from which to search for the nearest primitives.
     * @param maxPrimitiveCount the greatest number of primitives to retrieve.
     * @param[out] outResults information regarding the nearest primitives.
     * @param time the time at which examine the primitives.
     */
    void collect_nearest_primitives( const frantic::graphics::vector3f& point, std::size_t maxPrimitiveCount,
                                     std::vector<frantic::geometry::nearest_point_search_result>& outResults,
                                     double time = 0 ) const;
    /**
     *  Collect all primitives which are closer than maxDistance to the
     * specified point.
     *
     *  The primitives are returned in sorted order, and the list is cleared
     * before adding primitives.
     *
     * @param point the point from which to search for the nearest primitives.
     * @param maxDistance the maximum distance from the specified point to
     *		search for primitives.
     * @param[out] outResults information regarding the nearest primitives.
     * @param time the time at which examine the primitives.
     */
    void collect_primitives_within_range( const frantic::graphics::vector3f& point, float maxDistance,
                                          std::vector<frantic::geometry::nearest_point_search_result>& outResults,
                                          double time = 0 ) const;
    /**
     *  Collect up to maxPrimitiveCount of the primitives closest to the
     * specified point.
     *
     *  The primitives are returned in order of increasing distance, and the
     * list is cleared before adding primitives.
     *
     * @param point the point from which to search for the nearest primitives.
     * @param maxPrimitiveCount the greatest number of primitives to retrieve.
     * @param[out] outResults information regarding the nearest primitives.
     *		The output data follows the form of the kd-tree's
     *		channel map, which may be accessed using the get_channel_map()
     *		method.  The particle_array must have been created using this
     *		same channel map.
     * @param time the time at which examine the primitives.
     */
    void collect_nearest_primitives( const frantic::graphics::vector3f& point, std::size_t maxPrimitiveCount,
                                     frantic::particles::particle_array& outResults, double time = 0 ) const;
    /**
     *  Collect all primitives which are closer than maxDistance to the
     * specified point.
     *
     *  The primitives are returned in sorted order, and the list is cleared
     * before adding primitives.
     *
     * @param point the point from which to search for the nearest primitives.
     * @param maxDistance the maximum distance from the specified point to
     *		search for primitives.
     * @param[out] outResults information regarding the nearest primitives.
     *		The output data follows the form of the kd-tree's
     *		channel map, which may be accessed using the get_channel_map()
     *		method.  The particle_array must have been created using this
     *		same channel map.
     * @param time the time at which examine the primitives.
     */
    void collect_primitives_within_range( const frantic::graphics::vector3f& point, float maxDistance,
                                          frantic::particles::particle_array& outResults, double time = 0 ) const;
    // This collects all triangles indices within the sphere with center "point" and radius "maxDistance"
    /**
     *  Collect the face index number of mesh faces which are within
     * maxDistance of the specified search point.
     *
     * @note this is included only for compatibility with the trimesh3_kdtree.
     *		Its output probably isn't useful for you if the tree contains more
     *		than one mesh, or if it contains objects other than a mesh.
     *
     * @param point the point from which to search for the nearest primitives.
     * @param maxDistance the maximum distance from the specified point to
     *		search for faces.
     * @param[out] outFaceNumbers the face index number of mesh faces within
     *		maxDistance of the search point.
     * @param time the time at which examine the mesh.
     */
    void collect_faces_within_range( const frantic::graphics::vector3f& point, float maxDistance,
                                     std::vector<int>& outFaceNumbers, double time = 0 ) const;
    /**
     *  Collect all of the primitives within  maxDistance of the specified
     * search point.
     *
     * @note I'm not sure if the primitives should be directly available to
     *		callers...
     *
     * @param point the point from which to search for the nearest primitives.
     * @param maxDistance the maximum distance from the specified point to
     *		search for faces.
     * @param[out] outPrimitives the primitives within maxDistance of the search
     *		point.
     * @param time the time at which examine the mesh.
     */
    void collect_primitives_within_range( const frantic::graphics::vector3f& point, float maxDistance,
                                          std::vector<const mixed_kdtree_primitive*>& outPrimitives,
                                          double time = 0 ) const;

    // from volume_collection
    // I'm dropping this for now since its behaviour is covered
    // by find_nearest_point().
    // float get_distance_to_surface( const frantic::graphics::vector3f& point, frantic::graphics::vector3f*
    // pOutNearestPoint = 0, frantic::graphics::vector3f* pOutNormal = 0, double time = 0 ) const;

    /**
     *  Determine whether the specified point is within the geometry's
     * volume.  The inside/outside test is based on the geometry's normal,
     * whereby normals point outward from a volume.
     *
     * @param point the point to examine.
     * @param time the time at which examine the geometry.
     * @return true if the point is inside the geometry's volume, and false
     *		otherwise.
     */
    bool is_point_in_volume( const frantic::graphics::vector3f& point, double time = 0 ) const;

    /**
     *  Return the distance along the specified segment that is within
     * the geometry's volume.
     *
     * @param p1 the start of the segment.
     * @param p2 the end of the segment.
     * @param time the time at which examine the geometry.
     * @return the distance of the segment from p1 through p2 which is within
     *		the volume.
     */
    float get_segment_distance_in_volume( const vector3f& p1, const vector3f& p2, double time = 0 ) const;

    template <class distance_observer>
    inline void volume_traversal( const frantic::graphics::vector3f& point, distance_observer& observer,
                                  const double time = 0 );

    template <class ray_observer>
    inline void ray_traversal( const frantic::graphics::ray3f& ray, const double tMin, const double tMax,
                               ray_observer& observer, const double time = 0 );

    template <class ray_observer>
    inline void particle_traversal( const frantic::graphics::ray3f& ray, const double tMin, const double tMax,
                                    ray_observer& observer );
};

/**
 *  Virtual base class for objects which respond to changes in the
 * output channel map.
 */
class output_channel_map_listener {
  public:
    virtual ~output_channel_map_listener() {}
    virtual void set_output_channel_map( const frantic::channels::channel_map& channelMap ) = 0;
};

/**
 *  Used internally by mixed_kdtree to collect the primitives
 * from each geometry object in the tree.
 */
class mixed_kdtree_primitive_recorder {
    // vector of primitives encountered by the recorder
    // todo: this should maybe be a reference to a vector of some kind of
    // smart pointer, if the kd-tree itself is going to manage the primitives
    std::vector<boost::shared_ptr<mixed_kdtree_primitive>>& m_primitives;
    std::vector<output_channel_map_listener*>& m_outputChannelMapListeners;
    std::size_t m_maxDataSize;
    frantic::graphics::boundbox3f m_bounds;

    // undefined
    mixed_kdtree_primitive_recorder& operator=( mixed_kdtree_primitive_recorder& );

  public:
    mixed_kdtree_primitive_recorder( std::vector<boost::shared_ptr<mixed_kdtree_primitive>>&,
                                     std::vector<output_channel_map_listener*>& );

    /**
     *  Reserve space for the specified number of primitives.
     *
     * @param primitiveCount the number of primitives to reserve space
     *		for.
     */
    void reserve( const std::size_t primitiveCount );
    /**
     *  Add the primitive to the recorder's records.
     *
     * @param primitive the primitive whose shared_ptr and other data will
     *		be recorded.
     */
    void insert( boost::shared_ptr<mixed_kdtree_primitive> primitive );

    /**
     *  Add an output_channel_map_listener to the recorder's records.
     *
     * @param channelMapListener an object which should be notified
     *		when the kd-tree's channel map changes.
     */
    void insert_output_channel_map_listener( output_channel_map_listener* channelMapListener );

    // std::vector<mixed_kdtree_primitive *> & get_primitives_ref( void );
    /**
     *  Get the combined boundbox of all primitives encountered so far.
     */
    frantic::graphics::boundbox3f get_bounds( void ) const;
    /**
     *  Get the maximum data size of all primitives encountered so far.
     */
    std::size_t get_data_size( void ) const;
};

// kdtree objects hold kdtree primitives
// an object is responsible for managing the allocation of primitives
// the kdtree holds a pointer to each primitive, so you must not move them
// around in memory after you have passed them to the primitive recorder
class mixed_kdtree_primitive_creator {
  public:
    virtual ~mixed_kdtree_primitive_creator() {}

    /**
     *  Change the channel map which will be used for output data when
     * reporting intersections.
     *
     * @todo this should probably be separated into another class.
     */
    // virtual void set_output_channel_map( const frantic::channels::channel_map & ) {};
    /**
     *  Return the number of primitives which will be generated by a call
     * to get_primitives.
     *
     *  This isn't really necessary, for now at least I'm just using it to
     * preallocate the primitive vector.
     */
    virtual std::size_t get_primitive_count( void ) = 0;
    /**
     *  Should add all of the object's primitives to the recorder using
     * its insert method.
     *
     * @param recorder the recorder object into which the objects primitives
     *		should be inserted.
     */
    virtual void get_primitives( mixed_kdtree_primitive_recorder* recorder ) = 0;
};

//
// mixed_kdtree_primitive
//

/**
 *  This is an abstract class for primitives which can accept ray traversals
 * and nearest point searches.  Such primitives are stored at the leaf nodes
 * of the kd-tree.
 */
class mixed_kdtree_primitive {
  public:
    virtual ~mixed_kdtree_primitive() {}

    // time is fixed at t
    /**
     *  Attempt to find an intersection between the given ray and the primitive.
     * The ray is instantaneous, and the time is fixed at the given time.
     */
    virtual bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                mixed_kdtree_ray_observer* observer, double time = 0 ) = 0;

    // here the ray's t is the time
    // time is variable in [tMin,tMax]
    /**
     *  Attempt to find an intersection between the primitive and a particle
     * whose motion is described by the ray.  The particle moves from
     * ray.origin() at time 0, and moves by ray.direction() every time step.
     *
     * Here the time is a variable between tMin and tMax which you must solve
     * for if an intersection occured.
     */
    virtual bool intersect_particle( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                     mixed_kdtree_ray_observer* observer ) = 0;

    // cover this behaviour using traversal strategy
    // virtual void intersect_ray_all( const frantic::graphics::ray3f& ray, float segmentLength,
    // std::vector<raytrace_intersection>& outIntersections ) const;
    // virtual bool intersects_ray_segment( const frantic::graphics::vector3f& start, const frantic::graphics::vector3f&
    // end ) const; ?
    /**
     *  Find the nearest point on the primitive to the given point.
     */
    virtual bool find_nearest_point( const frantic::graphics::vector3f& point,
                                     frantic::graphics::boundbox3f& nodeBounds,
                                     mixed_kdtree_distance_observer* observer, double time = 0 ) = 0;

    // added just to support is_point_in_volume
    // todo: find a more appropriate way to implement is_point_in_volume
    // maybe split this into two parts:
    // - One accepts an input point, and returns a region code:
    //   interior/exterior/on surface/unknown
    // - another returns an outPoint + outNormal, perhaps closest to a
    //   query point
    virtual void get_any_point( vector3f& outPoint, vector3f& outNormal, double time = 0 ) = 0;

    /**
     *  Return the boundbox of this primitive.
     */
    virtual frantic::graphics::boundbox3f get_bounds() const = 0;
    /**
     *  Intersect this primitive's boundbox with the given boundbox.  Ideally
     * this should be a "perfect split," which does not include empty space in
     * the returned bounds.
     */
    virtual frantic::graphics::boundbox3f intersect_with( const frantic::graphics::boundbox3f& box ) const = 0;
    /**
     *  Return the number of bytes used by this primitive for its internal
     * representation of intersection points.
     */
    virtual std::size_t get_data_size( void ) const { return 0; }

    /**
     *  Convert this primitive's intersection data into a raytrace_intersection.
     */
    virtual void to_raytrace_intersection( frantic::geometry::raytrace_intersection& raytraceIntersection,
                                           const char* /*mixedData*/ ) const {
        raytraceIntersection.faceIndex = -1;
        raytraceIntersection.primitiveIndex = -1;
        raytraceIntersection.barycentricCoords = vector3f( 1.f / 3 );
    };

    /**
     *  Convert this primitive's intersection data into a
     * nearest_point_search_result.
     */
    virtual void to_nearest_point_search_result( frantic::geometry::nearest_point_search_result& searchResult,
                                                 const char* /*mixedData*/ ) const {
        searchResult.faceIndex = -1;
        searchResult.barycentricCoords = vector3f( 1.f / 3 );
    };

    /**
     *  Convert this primitive's intersection data into channel
     * mapped data.
     *
     */
    virtual void to_channel_map_data( char* outData, const mixed_kdtree_point_data* pointData ) const = 0;
};

/**
 *  Helper functor class to call a primitive's intersect_particle method
 * during a ray traversal.
 */
class do_intersect_particle {
  public:
    bool operator()( mixed_kdtree_primitive* primitive, const frantic::graphics::ray3f& ray, double tMin, double tMax,
                     mixed_kdtree_ray_observer& observer ) const {
        return primitive->intersect_particle( ray, tMin, tMax, &observer );
    }
};

/**
 * Helper functor class to call a primitive's intersect_ray method
 * during a ray traversal.
 */
class do_intersect_ray {
    double m_time;

  public:
    do_intersect_ray( void ) {}
    do_intersect_ray( const double time )
        : m_time( time ) {}
    bool operator()( mixed_kdtree_primitive* primitive, const frantic::graphics::ray3f& ray, double tMin, double tMax,
                     mixed_kdtree_ray_observer& observer ) const {
        return primitive->intersect_ray( ray, tMin, tMax, &observer, m_time );
    }
};

/**
 *  Absract base class for proximity-type searches.
 *
 *  Derived classes are intended to be used in calls to the kdtree nodes'
 * traverse_nearest_points method.
 */
class mixed_kdtree_distance_observer {
  public:
    virtual ~mixed_kdtree_distance_observer() {}
    /**
     *  Return true if the observer is its search.
     */
    virtual bool is_done( void ) = 0;
    /**
     *  Return the maximum range of this search, measured from the starting
     * point.
     */
    virtual double get_max_range( void ) = 0;
    /**
     *  Return true if a point at the specified distance is within the
     * search's range, or false otherwise.
     */
    virtual bool is_in_range( double d ) = 0;
    // virtual bool wants_hit_details( double d ) = 0;
    /**
     *  This is called by primitives to report their nearest point to
     * the search point.  You should store any information you want to
     * keep regarding the found point when this is called.
     */
    virtual void insert( const frantic::graphics::vector3f& pt, const double d, const double time,
                         const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* primitive,
                         const char* data = 0, const std::size_t dataSize = 0 ) = 0;
};

/**
 *  Used to determine the distance between a search point and
 * the nearest object in the kd-tree.
 */
class least_distance_mixed_kdtree_distance_observer : public mixed_kdtree_distance_observer {
    double m_distance;

  public:
    least_distance_mixed_kdtree_distance_observer();
    bool is_done( void );
    double get_max_range( void );
    bool is_in_range( double d );
    // bool wants_hit_details( double d );
    void insert( const frantic::graphics::vector3f& /*pt*/, const double d, const double /*time*/,
                 const frantic::graphics::vector3f& /*normal*/, mixed_kdtree_primitive* /*mixed*/,
                 const char* /*data = 0*/, const std::size_t /*dataSize = 0*/ );
    double get_distance( void );
};

/**
 *  Used to determine the point and corresponding primitive
 * which are nearest to a search point.
 */
class nearest_mixed_kdtree_distance_observer : public mixed_kdtree_distance_observer {
    double m_maxDistance;
    bool m_gotHit;
    mixed_kdtree_nearest_point& m_nearestPoint;

    // not implemented
    nearest_mixed_kdtree_distance_observer& operator=( const nearest_mixed_kdtree_distance_observer& );

  public:
    nearest_mixed_kdtree_distance_observer( mixed_kdtree_nearest_point& outResult,
                                            const double maxDistance = std::numeric_limits<double>::max() );
    bool is_done( void );
    bool got_hit( void );
    double get_max_range( void );
    bool is_in_range( double d );
    // bool wants_hit_details( double d );
    void insert( const frantic::graphics::vector3f& pt, const double d, const double time,
                 const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed, const char* data,
                 const std::size_t dataSize );
};

/**
 *  Find up to a specified number of primitives within a specified search
 * distance.
 */
class within_range_mixed_kdtree_distance_observer : public mixed_kdtree_distance_observer {
    float m_maxDistance;
    std::size_t m_maxCount;
    std::vector<mixed_kdtree_nearest_point>& m_nearestPoints;

    // not implemented
    within_range_mixed_kdtree_distance_observer& operator=( const within_range_mixed_kdtree_distance_observer& );

    bool is_collection_full( void );
    double get_greatest_distance( void );
    void bubble_up( std::size_t i );

  public:
    within_range_mixed_kdtree_distance_observer( std::vector<mixed_kdtree_nearest_point>& outResults,
                                                 const float maxDistance = std::numeric_limits<float>::max(),
                                                 const std::size_t maxCount = std::numeric_limits<std::size_t>::max() );
    bool is_done( void );
    double get_max_range( void );
    bool is_in_range( double d );
    // bool wants_hit_details( double d );
    void insert( const frantic::graphics::vector3f& pt, const double d, const double time,
                 const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed, const char* data = 0,
                 const std::size_t dataSize = 0 );
};

/**
 *  Abstract base class for ray traversal-type searches.
 *
 *  Derived classes are intended to be used in calls to the kd-tree nodes'
 * traverse_ray method.
 */
class mixed_kdtree_ray_observer {
  public:
    virtual ~mixed_kdtree_ray_observer() {}
    /**
     *  Return true if this search is complete.  This is used to
     * end the traversal early.
     */
    virtual bool is_done( void ) = 0;
    /**
     *  Return the maximum ray parameter ( that is the t in ray.at( t ) )
     * of interest in the search.
     */
    virtual double get_max_range( void ) = 0;
    /**
     *  Return true is the specified ray parameter ( t in ray.at( t ) )
     * is within the range of this search.
     */
    virtual bool is_in_range( double t ) = 0;
    // virtual bool wants_hit_details( double t ) = 0;
    // virtual void insert( const mixed_kdtree_ray_intersection & intersection ) = 0;
    /**
     *  This is called by a primitive when it finds an intersection between
     * the ray and itself.  You should store any information you need
     * regarding the intersection when this is called.
     */
    virtual void insert( const frantic::graphics::ray3f& ray, const double t, const double time,
                         const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* primitive, const char* data,
                         const std::size_t dataSize ) = 0; //{};
    virtual void insert( const frantic::graphics::ray3f& ray, const double t, const double time,
                         const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* primitive ) {
        insert( ray, t, time, normal, primitive, 0, 0 );
    };
};

/**
 *  For searches which return only a boolean indicating whether
 * or not a ray touches an object in the scene.
 */
class find_any_intersection_mixed_kdtree_ray_observer : public mixed_kdtree_ray_observer {
    bool m_gotHit;

  public:
    find_any_intersection_mixed_kdtree_ray_observer();
    bool is_done( void );
    double get_max_range( void );
    bool is_in_range( double /*t*/ );
    // bool wants_hit_details( double t );
    bool got_hit( void ) const;
    // void insert( const mixed_kdtree_ray_intersection & intersection );
    void insert( const frantic::graphics::ray3f&, const double t, const double, const frantic::graphics::vector3f&,
                 mixed_kdtree_primitive*, const char*, const std::size_t );
};

/**
 *  A search which finds the closest object along the ray.
 */
class find_nearest_intersection_mixed_kdtree_ray_observer : public mixed_kdtree_ray_observer {
    mixed_kdtree_ray_intersection& m_intersection;
    bool m_gotHit;

    // not implemented
    find_nearest_intersection_mixed_kdtree_ray_observer&
    operator=( const find_nearest_intersection_mixed_kdtree_ray_observer& );

  public:
    find_nearest_intersection_mixed_kdtree_ray_observer( mixed_kdtree_ray_intersection& intersection );
    bool is_done( void );
    double get_max_range( void );
    bool is_in_range( double t );
    // bool wants_hit_details( double t );
    // void insert( const mixed_kdtree_ray_intersection & intersection );
    void insert( const frantic::graphics::ray3f& ray, const double t, const double time,
                 const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed, const char* data,
                 const std::size_t dataSize );
    bool got_hit( void ) const;
};

// for now this just resets the point's time to the ray's t
// todo: this probably needs to be refactored
/*
class find_nearest_intersection_mixed_kdtree_particle_observer : public mixed_kdtree_ray_observer {
  mixed_kdtree_ray_intersection & m_intersection;
  bool m_gotHit;

    // not implemented
  find_nearest_intersection_mixed_kdtree_particle_observer & operator=( const
find_nearest_intersection_mixed_kdtree_particle_observer & );
public:
  find_nearest_intersection_mixed_kdtree_particle_observer( mixed_kdtree_ray_intersection & intersection );
  bool is_done( void );
  double get_max_range( void );
  bool is_in_range( double t );
  //bool wants_hit_details( double t );
  //void insert( const mixed_kdtree_ray_intersection & intersection );
  void insert( const frantic::graphics::ray3f & ray, const double t, const double time, const
frantic::graphics::vector3f & normal, mixed_kdtree_primitive * mixed, const char * data, const std::size_t dataSize );
  bool got_hit( void ) const;
};
*/

/**
 *  A search which finds all intersections between the ray and
 * objects in the kd-tree.
 */
class find_all_intersections_mixed_kdtree_ray_observer : public mixed_kdtree_ray_observer {
    // TODO: probably use something other than a vector for this
    std::vector<mixed_kdtree_ray_intersection>& m_intersections;

    // not implemented
    find_all_intersections_mixed_kdtree_ray_observer&
    operator=( const find_all_intersections_mixed_kdtree_ray_observer& );

  public:
    find_all_intersections_mixed_kdtree_ray_observer( std::vector<mixed_kdtree_ray_intersection>& intersections );
    bool is_done( void );
    bool is_in_range( double /*t*/ );
    double get_max_range( void );
    // bool wants_hit_details( double /*t*/ );
    // void insert( const mixed_kdtree_ray_intersection & intersection );
    void insert( const frantic::graphics::ray3f& ray, const double t, const double time,
                 const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed, const char* data,
                 const std::size_t dataSize );
};

constexpr std::size_t staticDataSize{ 80 };

/**
 *  Holds information about points on objects in the kd-tree.
 */
class mixed_kdtree_point_data {
    frantic::graphics::vector3f m_position;
    double m_time;
    frantic::graphics::vector3f m_normal;
    const mixed_kdtree_primitive* m_primitive;
    // now that this is not reused between queries (such as intersect_ray),
    // this will require a malloc for every query
    // TODO: this fixed buffer is a hack, replace it with something more
    // reasonable
    bool m_useStaticData;
    std::size_t m_primitiveDataSize;
    std::array<char, staticDataSize> m_staticData;
    std::vector<char> m_dynamicData;

    void set_my_data_from( const char* data, const std::size_t dataSize );

  public:
    mixed_kdtree_point_data();
    mixed_kdtree_point_data( const frantic::graphics::vector3f& position, const double time,
                             const frantic::graphics::vector3f& normal, const mixed_kdtree_primitive* mixed );
    mixed_kdtree_point_data( const frantic::graphics::vector3f& position, const double time,
                             const frantic::graphics::vector3f& normal, const mixed_kdtree_primitive* mixed,
                             const char* data, const std::size_t dataSize );

    void reset();
    void reset( const frantic::graphics::vector3f& position, const double time,
                const frantic::graphics::vector3f& normal, const mixed_kdtree_primitive* mixed, const char* data,
                const std::size_t dataSize );
    void reset( const frantic::graphics::vector3f& position, const double time,
                const frantic::graphics::vector3f& normal, const mixed_kdtree_primitive* mixed );

    const mixed_kdtree_primitive* get_primitive( void ) const;
    frantic::graphics::vector3f get_position( void ) const;
    double get_time( void ) const;
    void set_time( double time );
    frantic::graphics::vector3f get_normal( void ) const;
    const char* get_primitive_data( void ) const;
    std::size_t get_primitive_data_size( void ) const;
    void to_raytrace_intersection( frantic::geometry::raytrace_intersection& raytraceIntersection ) const;
    void to_nearest_point_search_result( frantic::geometry::nearest_point_search_result& searchResult ) const;

    template <class channel_map_data_writer>
    void to_channel_map_data( char* data, const channel_map_data_writer& w ) const {
        w.set_position( data, m_position );
        w.set_time( data, m_time );
        w.set_normal( data, m_normal );

        // need a method to set the data from the primitive
        m_primitive->to_channel_map_data( data, this );
    }
};

/**
 *  Holds information from nearest point-type searches.
 */
class mixed_kdtree_nearest_point_data {
  protected:
    double m_distance;

  public:
    mixed_kdtree_nearest_point_data();
    mixed_kdtree_nearest_point_data( const double distance );
    void reset( const double distance );
    double get_distance( void ) const;
    bool operator<( const mixed_kdtree_nearest_point_data& rhs ) const;
    void to_nearest_point_search_result( nearest_point_search_result& searchResult ) const;

    template <class channel_map_data_writer>
    void to_channel_map_data( char* data, const channel_map_data_writer& w ) const {
        w.set_distance( data, m_distance );
    }
};

/**
 *  Holds information from nearest point-type searches.
 *
 * @todo: combine mixed_kdtree_nearest_point_data with this ?
 */
class mixed_kdtree_nearest_point : public mixed_kdtree_nearest_point_data, public mixed_kdtree_point_data {
  public:
    mixed_kdtree_nearest_point();
    mixed_kdtree_nearest_point( const frantic::graphics::vector3f& pt, const double distance, const double time,
                                const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* primitive );
    mixed_kdtree_nearest_point( const frantic::graphics::vector3f& pt, const double distance, const double time,
                                const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* primitive,
                                const char* data, const std::size_t dataSize );
    void reset( const vector3f& pt, const double distance, const double time, const frantic::graphics::vector3f& normal,
                mixed_kdtree_primitive* primitive );
    void reset( const frantic::graphics::vector3f& pt, const double distance, const double time,
                const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed, const char* data,
                const std::size_t dataSize );
    void to_nearest_point_search_result( frantic::geometry::nearest_point_search_result& searchResult ) const;
    template <class channel_map_data_writer>
    void to_channel_map_data( char* data, const channel_map_data_writer& w,
                              const std::vector<char>& defaultChannelData ) const {
        memcpy( data, &defaultChannelData[0], defaultChannelData.size() );
        mixed_kdtree_nearest_point_data::to_channel_map_data( data, w );
        mixed_kdtree_point_data::to_channel_map_data( data, w );
    }
};

/*
// I don't think these are useful any longer ?
// I was using them temporarily to help find a target suitable for
// interior/exterior queries.
class find_nearest_non_coplanar_primitive : public mixed_kdtree_distance_observer {
    bool m_gotPoint;
  frantic::graphics::vector3f m_srcPoint;
    mixed_kdtree_nearest_point m_nearestPoint;
    float m_threshold;
    double m_minimumDistance;

public:
    find_nearest_non_coplanar_primitive( const frantic::graphics::vector3f & srcPoint, float dotThreshold )
        :   m_gotPoint( false ),
            m_srcPoint( srcPoint ),
            m_threshold( dotThreshold ),
            m_minimumDistance( std::numeric_limits<double>::max() )
    {}
    bool is_done() { return false; }
    double get_max_range( void ) { return m_minimumDistance; }
    bool is_in_range( double d ) { return d < m_minimumDistance; }
    //bool wants_hit_details( double d ) { return is_in_range( d ); }
    bool got_hit() { return m_gotPoint; }
    void insert( const frantic::graphics::vector3f & pt, const double d, const double time, const
frantic::graphics::vector3f & normal, mixed_kdtree_primitive * primitive, const char * data = 0, const std::size_t
dataSize = 0 ) {
        frantic::graphics::vector3f srcToDest = ( pt - m_srcPoint ).to_normalized();
        if( is_in_range( d ) && abs( frantic::graphics::vector3f::dot( normal, srcToDest ) ) > m_threshold ) {
            m_nearestPoint.reset( pt, d, time, normal, primitive, data, dataSize );
            m_gotPoint = true;
            m_minimumDistance = d;
        }
    }
    void report( frantic::geometry::nearest_point_search_result & result ) {
        m_nearestPoint.to_nearest_point_search_result( result );
    }
};

class find_any_non_coplanar_primitive : public mixed_kdtree_distance_observer {
    bool m_gotPoint;
  frantic::graphics::vector3f m_srcPoint;
    mixed_kdtree_nearest_point m_nearestPoint;
    float m_threshold;
    double m_minimumDistance;

public:
    find_any_non_coplanar_primitive( const frantic::graphics::vector3f & srcPoint, float dotThreshold )
        :   m_gotPoint( false ),
            m_srcPoint( srcPoint ),
            m_threshold( dotThreshold ),
            m_minimumDistance( std::numeric_limits<double>::max() )
    {}
    bool is_done() { return m_gotPoint; }
    double get_max_range( void ) { return m_minimumDistance; }
    bool is_in_range( double d ) { return d < m_minimumDistance; }
    //bool wants_hit_details( double d ) { return is_in_range( d ); }
    bool got_hit() { return m_gotPoint; }
    void insert( const frantic::graphics::vector3f & , const double d, const double time, const
frantic::graphics::vector3f & , mixed_kdtree_primitive * primitive, const char * data = 0, const std::size_t dataSize =
0 ) {
    vector3f point, normal;
    primitive->get_any_point( point, normal, time );

        frantic::graphics::vector3f srcToDest = ( point - m_srcPoint ).to_normalized();
        if( is_in_range( d ) && abs( frantic::graphics::vector3f::dot( normal, srcToDest ) ) > m_threshold ) {
            m_nearestPoint.reset( point, 0, time, normal, primitive, data, dataSize );
            m_gotPoint = true;
            //m_minimumDistance = d;
        }
    }
    void report( frantic::geometry::nearest_point_search_result & result ) {
        m_nearestPoint.to_nearest_point_search_result( result );
    }
};
*/

/**
 *  Holds information from ray intersections.
 */
class mixed_kdtree_ray_intersection_data {
    double m_t;
    frantic::graphics::ray3f m_ray;
    double m_distance;

  public:
    mixed_kdtree_ray_intersection_data();

    mixed_kdtree_ray_intersection_data( const frantic::graphics::ray3f& ray, const double t );

    void reset( const frantic::graphics::ray3f& ray, const double t );
    frantic::graphics::ray3f get_ray( void ) const;
    double get_t( void ) const;
    double get_distance( void ) const;
    bool operator<( const mixed_kdtree_ray_intersection_data& rhs ) const;
    void to_raytrace_intersection( frantic::geometry::raytrace_intersection& raytraceIntersection ) const;

    template <class channel_map_data_writer>
    void to_channel_map_data( char* data, const channel_map_data_writer& w ) const {
        w.set_ray( data, m_ray );
        w.set_t( data, m_t );
        w.set_distance( data, get_distance() );
    }
};

/**
 *  Holds information from ray traversal-type searches.
 *
 * @todo combine mixed_kdtree_ray_intersection_data into
 *		this ?
 */
class mixed_kdtree_ray_intersection : public mixed_kdtree_ray_intersection_data, public mixed_kdtree_point_data {
  public:
    mixed_kdtree_ray_intersection();
    mixed_kdtree_ray_intersection( const frantic::graphics::ray3f& ray, const double t, const double time,
                                   const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed );
    mixed_kdtree_ray_intersection( const frantic::graphics::ray3f& ray, const double t, const double time,
                                   const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed,
                                   const char* data, const std::size_t dataSize );

    void reset( const mixed_kdtree_ray_intersection& intersection );
    void reset( const frantic::graphics::ray3f& ray, const double t, const double time,
                const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed );
    void reset( const frantic::graphics::ray3f& ray, const double t, const double time,
                const frantic::graphics::vector3f& normal, mixed_kdtree_primitive* mixed, const char* data,
                const std::size_t dataSize );
    void to_raytrace_intersection( frantic::geometry::raytrace_intersection& raytraceIntersection ) const;

    template <class channel_map_data_writer>
    void to_channel_map_data( char* data, const channel_map_data_writer& w,
                              const std::vector<char>& defaultChannelData ) const {
        memcpy( data, &defaultChannelData[0], defaultChannelData.size() );
        mixed_kdtree_ray_intersection_data::to_channel_map_data( data, w );
        mixed_kdtree_point_data::to_channel_map_data( data, w );
    }
};

/**
 *  Determine the ray parameter at the intersection between the given
 * ray and the given axis-aligned splitting plane.  Find the minimum and
 * maximum bounds of this parameter due to rounding error.
 *
 *  The returned [outMinDistance,outMaxDistance] should produce ray points
 * such that ray.at(outMinDistance) and ray.at(outMaxDistance) are on
 * opposite sides of the splitting plane.
 *
 *  I assume that the denominator is finite because of the control paths
 * which reach this in the ray traversal code.
 *
 * @todo I think this is overkill and should be simplified.
 * Maybe just precalculate a pessimistic eps for each ray at the start
 * of traverse_ray ?  as it is now, this relatively slow function is
 * evaluated for every splitting plane encountered during the
 * ray intersection.
 *
 * @param rayOriginAlongAxis the ray's origin on the splitting
 *		plane's axis.
 * @param rayDirectionAlongAxis the ray's direction along the
 *		splitting plane's axis.
 * @param split the location of the splitting plane along the
 *		splitting plane's axis.
 * @param[out] outMinDistance the minimum distance along the ray
 *		for the intersection between the ray and the splitting plane.
 * @param[out] outMaxDistance the maximum distance along the ray
 *		for the intersection between the ray and the splitting plane.
 */
inline void get_distance_interval_for_ray_splitting_plane_intersection( const float rayOriginAlongAxis,
                                                                        const float rayDirectionAlongAxis,
                                                                        const float split, float& outMinDistance,
                                                                        float& outMaxDistance ) {
    const float eps = std::numeric_limits<float>::epsilon();
    const float num = split - rayOriginAlongAxis;
    const float numDelta = eps * ( fabs( split ) + fabs( rayOriginAlongAxis ) );
    const float numMin = num - numDelta;
    const float numMax = num + numDelta;
    const float den = rayDirectionAlongAxis;
    const float denDelta = eps * rayDirectionAlongAxis;
    // const float den0 = den + abs(denDelta);
    // const float den1 = den - abs(denDelta);

    if( den > 0 ) {
        if( numMin > 0 ) {
            // it's positive so move min closer to 0
            outMinDistance = numMin / ( den + denDelta );
        } else {
            outMinDistance = numMin / ( den - denDelta );
        }
        if( numMax > 0 ) {
            // farther
            outMaxDistance = numMax / ( den - denDelta );
        } else {
            outMaxDistance = numMax / ( den + denDelta );
        }
    } else {
        if( numMin > 0 ) {
            // it's negative so move max closer to 0
            outMaxDistance = numMin / ( den + denDelta );
        } else {
            outMaxDistance = numMin / ( den - denDelta );
        }
        if( numMax > 0 ) {
            outMinDistance = numMax / ( den - denDelta );
        } else {
            outMinDistance = numMax / ( den + denDelta );
        }
    }
}

// This returns the nearest intersection between the ray and the mesh, or false if none is found.
template <class ray_observer, class ray_intersection_type>
void mixed_kdtree_node::traverse_ray( const frantic::graphics::ray3f& ray, const double tMin, double tMax,
                                      ray_observer& observer, const ray_intersection_type& rayCollisionTest ) {
    // std::cout << "intersect ray " << tMin << " " << tMax << "\n";
    if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
        // std::cout << "at leaf.\n";

        if( m_primitivePointers != 0 ) {
            tMax = std::min( tMax, observer.get_max_range() );

            for( std::vector<mixed_kdtree_primitive*>::iterator primitive = m_primitivePointers->begin();
                 primitive != m_primitivePointers->end(); ++primitive ) {
                // std::cout << "*** someone at leaf !\n";
                // if( primitive->intersect_ray( ray, tMin, tMax, & observer ) ) {
                // if( ( primitive->*rayCollisionTest )( ray, tMin, tMax, & observer, time ) ) {
                if( rayCollisionTest( *primitive, ray, tMin, tMax, observer ) ) {
                    // std::cout << "hit for ray " << ray << "\n";
                    // observer.insert( intersection );
                    // either get the tMax from the observer's range,
                    // or just have the primitives handle it themselves
                    // std::cout << "found = true\n";
                    tMax = std::min( tMax, observer.get_max_range() );
                    if( observer.is_done() ) {
                        return;
                    }
                }
            }
        } // else {  std::cout << " no one at leaf.\n"; }
    } else {
        // It's an interior node.
        const int axis = m_axisAndLeafFlag & 0x00000003;

        // If the ray is parallel to the plane, the distance computed would be invalid, so we special case it
        const float rayDirectionAlongAxis = ray.direction()[axis];
        if( rayDirectionAlongAxis == 0 ) {
            // When the ray is parallel, we only need to check one of the sides
            // maybe both when it's on ?  this depends on how the triangle-ray
            // intersections behave
            if( ray.origin()[axis] <= m_split ) {
                left_child()->traverse_ray( ray, tMin, tMax, observer, rayCollisionTest );
            } else { // if( ray.origin()[axis] > m_split ) {
                right_child()->traverse_ray( ray, tMin, tMax, observer, rayCollisionTest );
            }
            /*
            else {
              bool found = left_child()->traverse_ray( ray, tMin, tMax, observer );
              return right_child()->traverse_ray( ray, tMin, tMax, observer ) || found;
            }
            */
        } else {
            double distance = ( (double)m_split - ray.origin()[axis] ) / rayDirectionAlongAxis;

            // Traverse the children in the order so the closest child is tried first
            if( rayDirectionAlongAxis > 0 ) {
                // are these cases ok ?
                // I think they are because you would only have problems
                // if you are near your ray's start or end point.
                // Since the splitting planes have epsilons added (below),
                // I think this will only be problematic with the user-
                // specified ray boundaries, so it's a clipping problem
                // that should probably be solved outside of this scope.
                // However I have not tested this.
                if( distance < tMin ) {
                    // The plane is to the left of the ray segment, so only consider the right side
                    right_child()->traverse_ray( ray, tMin, tMax, observer, rayCollisionTest );
                } else if( distance > tMax ) {
                    // The plane is to the right of the ray segment, so only consider the left side
                    left_child()->traverse_ray( ray, tMin, tMax, observer, rayCollisionTest );
                } else {
                    // distance +/- float eps is insufficient, but
                    // this is overkill and may be incorrect.
                    // However it at least catches the test data
                    // that distance +/- eps missed.
                    // This assumes that the distance >= 0,
                    // so minDistance >= 0 and maxDistance >= 0,
                    // but I think this is required elsewhere as well ?
                    /*
                              const float num = m_split - ray.origin()[axis];
                              const float numEps = (m_split + ray.origin()[axis]) * eps;
                              const float den = rayDirectionAlongAxis;
                              const float num0 = num + numEps;
                              const float num1 = num - numEps;
                              const float den1 = den * ( 1 - eps );
                              const float B = num1 / den1;
                              const float D = num0 / den1;
                              const double maxDistance = std::min<double>( std::max(B, D), tMax );
                    */
                    // TODO: find a better way to do this
                    // compensate for error when refining the active ray segment

                    float minDistanceFloat, maxDistanceFloat;
                    get_distance_interval_for_ray_splitting_plane_intersection(
                        ray.origin()[axis], rayDirectionAlongAxis, m_split, minDistanceFloat, maxDistanceFloat );
                    const double minDistance = std::max( tMin, static_cast<double>( minDistanceFloat ) );
                    const double maxDistance = std::min( tMax, static_cast<double>( maxDistanceFloat ) );

                    left_child()->traverse_ray( ray, tMin, maxDistance, observer, rayCollisionTest );
                    if( !observer.is_done() && observer.is_in_range( minDistance ) ) {
                        /*
                        const float den0 = den * ( 1 + eps );
                        const float A = num1 / den0;
                        const float C = num0 / den0;
                        const double minDistance = std::max<double>( std::min(A, C), tMin );
                        */
                        right_child()->traverse_ray( ray, minDistance, tMax, observer, rayCollisionTest );
                    }
                }
            } else {
                if( distance < tMin ) {
                    // The plane is to the left of the ray segment, so only consider the left side
                    left_child()->traverse_ray( ray, tMin, tMax, observer, rayCollisionTest );
                } else if( distance > tMax ) {
                    // The plane is to the right of the ray segment, so only consider the right side
                    right_child()->traverse_ray( ray, tMin, tMax, observer, rayCollisionTest );
                } else {
                    float minDistanceFloat, maxDistanceFloat;
                    get_distance_interval_for_ray_splitting_plane_intersection(
                        ray.origin()[axis], rayDirectionAlongAxis, m_split, minDistanceFloat, maxDistanceFloat );
                    const double minDistance = std::max( tMin, static_cast<double>( minDistanceFloat ) );
                    const double maxDistance = std::min( tMax, static_cast<double>( maxDistanceFloat ) );

                    right_child()->traverse_ray( ray, tMin, maxDistance, observer, rayCollisionTest );

                    if( !observer.is_done() && observer.is_in_range( minDistance ) ) {
                        left_child()->traverse_ray( ray, minDistance, tMax, observer, rayCollisionTest );
                    }
                }
            }
        }
    }
}

template <class distance_observer>
void mixed_kdtree_node::traverse_nearest_points( const frantic::graphics::vector3f& point,
                                                 frantic::graphics::boundbox3f& nodeBounds, distance_observer& observer,
                                                 double time ) {
    if( m_axisAndLeafFlag < 0 ) { // Test if this node is a leaf.
        if( m_primitivePointers ) {
            // nearest_point_search_result testNearestPoint;
            // double bestDistance = std::numeric_limits<double>::infinity();

            for( std::vector<mixed_kdtree_primitive*>::iterator primitive = m_primitivePointers->begin();
                 primitive != m_primitivePointers->end(); ++primitive ) {
                ( *primitive )->find_nearest_point( point, nodeBounds, &observer, time );

                if( observer.is_done() ) {
                    return;
                }
            }
        }
    } else {
        // It's an interior node
        int axis = m_axisAndLeafFlag & 0x00000003;

        // If we're too far away from the bounding box, don't even bother looking.  Use clamp_nothrow to avoid
        // the is_empty check, for efficiency.
        if( vector3f::distance_squared( nodeBounds.clamp_nothrow( point ), point ) >=
            observer.get_max_range() * observer.get_max_range() )
            return;

        // Traverse the children so the closest child is tried first
        if( point[axis] < m_split ) {
            float savedBoundsValue = nodeBounds.maximum()[axis];
            nodeBounds.maximum()[axis] = m_split;
            // First find the nearest point in the left child
            left_child()->traverse_nearest_points( point, nodeBounds, observer,
                                                   time ); // primitives, outNearestPoint );
            nodeBounds.maximum()[axis] = savedBoundsValue;

            // If found, constrict the search radius for the search in the right child
            // if( found )
            // maxDistance = (float)outNearestPoint.distance;

            // todo: logic here can be improved to avoid actually calling
            // the primitives when we're finishing up
            if( !observer.is_done() ) {
                savedBoundsValue = nodeBounds.minimum()[axis];
                nodeBounds.minimum()[axis] = m_split;
                // Then find the nearest point in the right child, using this new constrained radius
                right_child()->traverse_nearest_points( point, nodeBounds, observer,
                                                        time ); // primitives, outNearestPoint ) )
                // if( right_child()->traverse_nearest_points( point, maxDistance, nodeBounds, primitives,
                // outNearestPoint ) ) found = true;
                nodeBounds.minimum()[axis] = savedBoundsValue;
            }
        } else {
            float savedBoundsValue = nodeBounds.minimum()[axis];
            nodeBounds.minimum()[axis] = m_split;
            // First find the nearest point in the right child
            // bool found = right_child()->traverse_nearest_points( point, maxDistance, nodeBounds, primitives,
            // outNearestPoint );
            right_child()->traverse_nearest_points( point, nodeBounds, observer,
                                                    time ); // primitives, outNearestPoint );
            nodeBounds.minimum()[axis] = savedBoundsValue;

            // If found, constrict the search radius for the search in the left child
            // if( found )
            // maxDistance = (float)outNearestPoint.distance;
            if( !observer.is_done() ) {
                savedBoundsValue = nodeBounds.maximum()[axis];
                nodeBounds.maximum()[axis] = m_split;
                // Then find the nearest point in the left child, using this new constrained radius
                // if( left_child()->traverse_nearest_points( point, maxDistance, nodeBounds, primitives,
                // outNearestPoint ) )
                left_child()->traverse_nearest_points( point, nodeBounds, observer,
                                                       time ); // primitives, outNearestPoint ) )
                // found = true;
                nodeBounds.maximum()[axis] = savedBoundsValue;
            }
        }
    }
}

template <class distance_observer>
void mixed_kdtree::volume_traversal( const frantic::graphics::vector3f& point, distance_observer& observer,
                                     const double time ) {
    m_rootNode->traverse_nearest_points( point, m_rootBounds, observer, time );
}

template <class ray_observer>
void mixed_kdtree::ray_traversal( const frantic::graphics::ray3f& ray, const double tMin, const double tMax,
                                  ray_observer& observer, const double time ) {
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_ray( time ) );
}

template <class ray_observer>
void mixed_kdtree::particle_traversal( const frantic::graphics::ray3f& ray, const double tMin, const double tMax,
                                       ray_observer& observer ) {
    m_rootNode->traverse_ray( ray, tMin, tMax, observer, do_intersect_particle() );
}
} // namespace geometry
} // namespace frantic
