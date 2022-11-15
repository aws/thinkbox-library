// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/quickhull_generic.hpp>
#include <frantic/graphics2d/vector2.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#include <frantic/math/exact_predicates.h>

namespace frantic {
namespace geometry {
namespace detail {

class Facet2Comp;
class Facet2Hash;

class QHTraits_2D {
  public:
    typedef int RidgeType;
    typedef std::less<int> RidgeComp;

    typedef Facet<2, RidgeType, RidgeComp> FacetType;
    typedef Facet2Comp FacetComp;
    typedef Facet2Hash FacetHash;

    typedef graphics2d::vector2f VertexType;
    typedef graphics2d::vector2 OutputType;

  private:
    static FacetType& init_facet( FacetType& f, int i0, int i1 ) {
        f.set_vertex( 0, i0 );
        f.set_ridge( 0, i0 );
        f.set_vertex( 1, i1 );
        f.set_ridge( 1, i1 );
        return f;
    }

    static const float distance_to_facet( const VertexType& v0, const VertexType& v1, const VertexType& v2 ) {
        // TODO: maybe add this back if the extra precision is needed
        // return -math::orient2d((float*)&v0, (float*)&v1, (float*)&v2);

        // This is the dot product of the edge normal w/ pt - dot of edge normal w/ the first vertex
        // This is fast because the edge normal is not normalized, meaning the distance is not a real distance
        return ( v1.y - v0.y ) * ( v2.x - v0.x ) + ( v0.x - v1.x ) * ( v2.y - v0.y );
    }

  public:
    static const float distance_to_facet( const std::vector<VertexType>& pts, const FacetType& f, int p );
    static const OutputType primitive_from_facet( const FacetType& f );
    static FacetType facet_from_ridge( const RidgeType& r, const FacetType& f, int p );

    template <class FacetCollection>
    static void build_initial_simplex( const std::vector<VertexType>& pts, FacetCollection& facets );
};

inline const float QHTraits_2D::distance_to_facet( const std::vector<VertexType>& pts, const FacetType& f, int p ) {
    return distance_to_facet( pts[f.get_vertex( 0 )], pts[f.get_vertex( 1 )], pts[p] );
}

inline const QHTraits_2D::OutputType QHTraits_2D::primitive_from_facet( const FacetType& f ) {
    return OutputType( f.get_vertex( 0 ), f.get_vertex( 1 ) );
}

inline QHTraits_2D::FacetType QHTraits_2D::facet_from_ridge( const RidgeType& r, const FacetType& f, int p ) {
    FacetType result;
    if( f.get_vertex( 1 ) == r )
        return init_facet( result, r, p );
    else
        return init_facet( result, p, r );
}

template <class FacetCollection>
inline void QHTraits_2D::build_initial_simplex( const std::vector<VertexType>& pts, FacetCollection& facets ) {
    int i0 = 0, i1 = 1, i2 = 2;
    if( distance_to_facet( pts[i0], pts[i1], pts[i2] ) > 0 )
        std::swap( i1, i2 );

    FacetType result;

    facets.insert( init_facet( result, i0, i1 ) );
    facets.insert( init_facet( result, i1, i2 ) );
    facets.insert( init_facet( result, i2, i0 ) );
}

class Facet2Comp {
  public:
    bool operator()( const QHTraits_2D::FacetType& f1, const QHTraits_2D::FacetType& f2 ) const {
        return f1.get_vertex( 0 ) < f2.get_vertex( 0 ) ||
               ( f1.get_vertex( 0 ) == f2.get_vertex( 0 ) && ( f1.get_vertex( 1 ) < f2.get_vertex( 1 ) ) );
    }
};

class Facet2Hash {
  public:
    size_t operator()( const QHTraits_2D::FacetType& f ) const { return f.get_vertex( 0 ) + 7 * f.get_vertex( 1 ); }
};

} // namespace detail

inline void quickhull( const std::vector<graphics2d::vector2f>& pts, std::vector<graphics2d::vector2>& hull ) {
    detail::quickhull<detail::QHTraits_2D>( pts, hull );
}

} // namespace geometry
} // namespace frantic
