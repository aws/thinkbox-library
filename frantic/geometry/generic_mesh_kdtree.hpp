// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/raytracing.hpp>
#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/graphics/ray3f.hpp>

#include <boost/scoped_ptr.hpp>

namespace frantic {
namespace geometry {

// Forward declaration of the type used as a node in generic_mesh_kdtree
struct generic_mesh_kdtree_node;

/**
 * This class provides a kdtree acceleration structure around a generic mesh to accelerate
 * operations such as ray casts and nearest point on the surface. A traits class is required
 * to provide the mesh-specific stuff. The tree only stores indices into the mesh object.
 * @tparam MeshTraits This template parameter supplies typedefs and static functions for manipulating
 *                    the underlying mesh primitive. This includes:
 *
 *                    typedef mesh_type The concrete type used as the source of mesh data
 *                    typedef mesh_type_const_ptr The type used to store a pointer to the mesh_type object. Should
 * manage the lifetime of the data. typedef raytrace_result The type storing information about a raytrace in through the
 * tree. typedef closest_point_result The type storing information about a nearest point search in the tree.
 *
 *                    inline static unsigned get_count( const mesh_type& inst );
 *                             Extracts the number of primitives (faces) in the mesh.
 *                    inline static frantic::graphics::boundbox3f get_bounds();
 *                             Extracts the bounding box of all primitives in the mesh.
 *                    inline static frantic::graphics::boundbox3f get_clipped_bounds( const mesh_type& inst, unsigned
 * index, const frantic::graphics::boundbox3f& bounds ); Extracts the bounds of a single primitive (ie. face) after
 * clipping the primitive to the supplied bounding box. inline static bool intersect_ray( const mesh_type& inst,
 * unsigned index, const frantic::graphics::ray3f& ray, double tMin, double tMax, raytrace_result& outIntersection );
 *                             Interestects the given ray with the specified primitive. Returns true if an intersection
 * occurred in the (tMin,tMax) range. If returning true, then 'outIntersection' is updated with the relevant
 * intersection data.
 */
template <class MeshTraits>
class generic_mesh_kdtree {
    // Typedef to simplify the name of the node type used in the tree.
    typedef generic_mesh_kdtree_node node_type;

    // Simplify the type used for returning raytracing results
    typedef typename MeshTraits::raytrace_result raytrace_result;

    // Simplify the type used for returning nearest point query results
    typedef typename MeshTraits::nearest_point_result nearest_point_result;

    // Store a pointer to the mesh. If this is a boost::shared_ptr, then we don't have to worry about the lifetime of
    // the mesh being different from the tree.
    typename MeshTraits::mesh_type_const_ptr m_mesh;

    // Store a pointer to the root node of the tree.
    boost::scoped_ptr<node_type> m_rootNode;

    // Store the bounding box of the entire tree.
    frantic::graphics::boundbox3f m_rootBounds;

  public:
    /**
     * Default constructor. User must call set_mesh() and finalize() before using the tree.
     */
    generic_mesh_kdtree();

    /**
     * Constructor that initializes the mesh ptr. If final == false, the user must call finalize()
     * before using the tree.
     * @param mesh A ptr to the mesh to create the kdtree for.
     * @param final If true, the tree is built immediately. If false, the user must call finalize()
     *              before searching the tree.
     */
    generic_mesh_kdtree( typename MeshTraits::mesh_type_const_ptr mesh, bool final = true );

    /**
     * Virtual destructor so we can subclass.
     */
    virtual ~generic_mesh_kdtree();

    /**
     * Set the mesh the tree is built on. After calling this, the user must call finalize() to actually
     * build the tree. That doesn't happen automatically so that we can modify the mesh while storing
     * the pointer internally.
     * @param mesh The mesh to build the tree for.
     */
    void set_mesh( typename MeshTraits::mesh_type_const_ptr mesh );

    /**
     * @return A const reference to the mesh data the tree indexes
     */
    const typename MeshTraits::mesh_type& get_mesh() const { return *m_mesh; }

    /**
     * @return A const pointer to the mesh data the tree indexes. This is different from &get_mesh() if using smart
     * pointers.
     */
    typename MeshTraits::mesh_type_const_ptr get_mesh_ptr() const { return m_mesh; }

    /**
     * This function actually builds the kdtree on the mesh set by the constructor or set_mesh(). This must be called
     * once before querying the tree with rays or nearest point searches.
     */
    void finalize();

    ///// Queries/////

    /**
     * This returns the first intersection between the ray and the mesh, or false if none is found.
     * @param ray The ray to search along in the tree
     * @param tMin Any intersection along the ray must be strictly further than 'ray.at(tMin)'. A good default is 0.0f.
     * @param tMax Any intersection along the ray must be strictly closer than 'ray.at(tMax)'. A good default is
     * std::numeric_limits<float>::max().
     * @param[out] outIntersection If an intersection is found, the information about it will be stored here.
     * @returns true if intersection found, false otherwise.
     */
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax, raytrace_result& outIntersection,
                        bool ignoreBackfaces = false ) const;

    /**
     * Find the point on the mesh, closest to the query location within a specified search radius.
     * @param point the point to search from.
     * @param maxDistance the maximum distance away from the specified point to search.
     * @param[out] outNearestPoint If a point is found, the information about it will be stored here.
     * @return true if a point could be found within the specified search radius, and false otherwise.
     */
    bool find_nearest_point( const frantic::graphics::vector3f& point, double maxDistance,
                             nearest_point_result& outNearestPoint, bool ignoreBackfaces = false ) const;

  private:
    bool intersect_ray_recursive( const node_type& curNode, const frantic::graphics::ray3f& ray, double tMin,
                                  double tMax, raytrace_result& outIntersection, bool ignoreBackfaces = false ) const;

    // bool find_nearest_point_recursive( const node_type& curNode, const frantic::graphics::vector3f& p, double tMax,
    // nearest_point_result& outResult ) const;
    bool find_nearest_point_recursive( const node_type& curNode, const frantic::graphics::vector3f& p, double tMax,
                                       frantic::graphics::boundbox3f& currentBounds, nearest_point_result& outResult,
                                       bool ignoreBackfaces = false ) const;
};

} // namespace geometry
} // namespace frantic

// This include is OK being here since its just being used as a convenient way to separate interface and implementation.
// I can't use a .cpp file because this is a template class.
#include <frantic/geometry/generic_mesh_kdtree_impl.hpp>
