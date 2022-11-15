// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/unordered_set.hpp>
#include <map>

#include <assert.h>
#include <frantic/diagnostics/profiling_manager.hpp>

namespace frantic {
namespace geometry {
namespace detail {

template <int Dim, class RidgeType_, class RidgeComp = std::less<RidgeType_>>
class Facet {
  public:
    typedef RidgeType_ RidgeType;
    static const int DIMENSION = Dim;

  private:
    std::vector<int> outside_pts;
    int far_point;
    float far_dist;

    Facet* links[Dim];
    RidgeType ridges[Dim];
    int verts[Dim];

    unsigned search_id;

  public:
    Facet()
        : far_point( -1 )
        , far_dist( 0 )
        , search_id( 0 ) {}

    static Facet* UNLINKED_VALUE() { return (Facet*)1; }

    void add_outside( float d, int p ) {
        if( d > far_dist ) {
            far_dist = d;
            far_point = p;
        }
        outside_pts.push_back( p );
    }

    bool operator==( const Facet& other ) const {
        bool isEqual = true;

        for( int i = 0; i < Dim; ++i ) {
            isEqual &= verts[i] == other.verts[i];
        }

        return isEqual;
    }

    int num_outside() const { return (int)outside_pts.size(); }
    int get_outside( int i ) const { return outside_pts[i]; }
    int furthest_outside() const { return far_point; }

    void set_search_id( unsigned s ) { search_id = s; }
    unsigned get_search_id() { return search_id; }

    void set_vertex( int i, int v ) { verts[i] = v; }
    int get_vertex( int i ) const { return verts[i]; }

    void set_ridge( int i, const RidgeType& r ) {
        ridges[i] = r;
        links[i] = UNLINKED_VALUE();
    }
    const RidgeType& get_ridge( int i ) const { return ridges[i]; }

    // TODO: Change this linear search to a simple binary tree (std::map) perhaps.
    Facet* get_link_indexed( int i ) const { return links[i]; }
    Facet* get_link( const RidgeType& r ) const {
        RidgeComp comp;
        for( int i = 0; i < Dim; ++i ) {
            if( !comp( r, ridges[i] ) && !comp( ridges[i], r ) )
                return links[i];
        }

        return NULL;
    }

    void set_link_indexed( int i, Facet* f ) { links[i] = f; }
    void set_link( const RidgeType& r, Facet* f ) {
        RidgeComp comp;
        for( int i = 0; i < Dim; ++i ) {
            if( !comp( r, ridges[i] ) && !comp( ridges[i], r ) ) {
                links[i] = f;
                return;
            }
        }

        throw std::runtime_error( "ERROR: Tried to set Facet link by ridge that did not exist." );
    }
};

template <class Traits>
void build_horizon( std::vector<std::pair<const typename Traits::RidgeType*, typename Traits::FacetType*>>& horizon,
                    std::vector<typename Traits::FacetType*>& visible, typename Traits::FacetType& cur,
                    const std::vector<typename Traits::VertexType>& pts, int pt, const unsigned SEARCH_ID ) {
    cur.set_search_id( SEARCH_ID );
    visible.push_back( &cur );

    for( int i = 0; i < Traits::FacetType::DIMENSION; ++i ) {
        typename Traits::FacetType& f = *( cur.get_link_indexed( i ) );

        float d = Traits::distance_to_facet( pts, f, pt );
        if( d > 0 ) {
            if( f.get_search_id() != SEARCH_ID )
                build_horizon<Traits>( horizon, visible, f, pts, pt, SEARCH_ID );
        } else
            horizon.push_back(
                std::pair<const typename Traits::RidgeType*, typename Traits::FacetType*>( &cur.get_ridge( i ), &f ) );
    }
}

#if defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4127 )
#endif

// The Traits class required for the quickhull algorithm is responsible for defining
// various dimension/application specific functions and typedefs
//
// the typedefs needed are:
//	RidgeType;	// The class that represents the boundary of two facets (ie. an edge in 3D)	(D-1 pts)
//	RidgeComp;	// Class that implements the hash_compare std class for use in hash_maps
//	FacetType;	// The class that represents the boundary of two simplices (ie. a triangle in 3D),
//              // a specialization of Facet<> (D pts)
//	FacetComp;	// Class that implements operator() as a less-than predicate for FacetType
//  FacetHash;  // Class that implements operator() as a hash predicate for FacetType
//	VertexType;	// What type of vertices are used
//	OutputType;	// The output type (ie. triangles for 3D convex hull, tetrahedra for 3D delaunay triangulation)
//
// The exact functions to implement are:
//  // Exact or approximation of distance of point to a facet
//	static const float distance_to_facet(const std::vector< VertexType >& pts, FacetType& f, int p);
//  // Extract an OutputType from a Facet
//	static const OutputType primitive_from_facet(FacetType& f);
//  // Create a new Facet connecting to another Facet on the specified ridge
//	static FacetType facet_from_ridge(const RidgeType& r, FacetType&, int p);
//   // Create the first simplex (D+1 pts) to bootstrap the algorithm
//	static void build_initial_simplex<class FacetCollection>(const std::vector< VertexType >& pts,
//                                                           FacetCollection& facets);
template <class Traits>
void quickhull( const std::vector<typename Traits::VertexType>& pts, std::vector<typename Traits::OutputType>& hull ) {
    typedef typename Traits::FacetType Facet;
    typedef boost::unordered_set<Facet, typename Traits::FacetHash> FacetSet;
    typedef std::map<typename Traits::RidgeType, Facet*, typename Traits::RidgeComp> FacetMap;
    typedef boost::unordered_set<Facet*> FacetPointerSet;

    FacetSet facets;

    // A quick and dirty way I sped up the search routine was to have a hash_set of pointers to
    // facets that contained at least one outside point. This hash_set is used for that.
    FacetPointerSet activeFacets;

    if( pts.size() <= Facet::DIMENSION )
        return;

    // Profiling stuff
    /*diagnostics::profiling_manager pfm;
    int initSection = pfm.new_profiling_section("Initialization");
    int loopSection = pfm.new_profiling_section("Loop");
    int loopSearchSection = pfm.new_profiling_section("Loop->Search");
    int loopHorizonSection = pfm.new_profiling_section("Loop->Horizon");
    int loopBuildSection = pfm.new_profiling_section("Loop->Build");
    int loopDestroySection = pfm.new_profiling_section("Loop->Destroy");
    int loopBuildInsertSection = pfm.new_profiling_section("Loop->Build->Insert");
    int loopBuildUpdateSection = pfm.new_profiling_section("Loop->Build->Update");*/

    // pfm.enter_section(initSection);

    // Build the initial simplex;
    Traits::build_initial_simplex( pts, facets );

    // Add each pt to the outside set of exactly one or zero facets.
    for( int i = Facet::DIMENSION + 1; i < (int)pts.size(); ++i ) {
        for( typename FacetSet::iterator it1 = facets.begin(); it1 != facets.end(); ++it1 ) {
            float d = Traits::distance_to_facet( pts, *it1, i );
            if( d > 0 ) {
                const_cast<typename FacetSet::value_type&>( *it1 ).add_outside( d, i );
                break;
            }
        }
    }

    // A quick and dirty way I sped up the search routine was to have a hash_set of pointers to
    // facets that contained at least one outside point.
    for( typename FacetSet::iterator it1 = facets.begin(); it1 != facets.end(); ++it1 ) {
        if( it1->num_outside() > 0 )
            activeFacets.insert( const_cast<typename FacetSet::value_type*>( &*it1 ) );
    }

    // Link all the initial facets together by trying to link each facet to every other ridge.
    for( typename FacetSet::iterator it1 = facets.begin(); it1 != facets.end(); ++it1 ) {
        typename FacetSet::iterator it2 = it1;
        while( ++it2 != facets.end() ) {
            // Check if these two facets share a ridge ... if so link them
            for( int i = 0; i < Facet::DIMENSION; ++i ) {
                const typename Facet::RidgeType& r = it1->get_ridge( i );
                if( it2->get_link( r ) != NULL ) {
                    const_cast<typename FacetSet::value_type&>( *it1 ).set_link(
                        r, const_cast<typename FacetSet::value_type*>( &*it2 ) );
                    const_cast<typename FacetSet::value_type&>( *it2 ).set_link(
                        r, const_cast<typename FacetSet::value_type*>( &*it1 ) );
                    break;
                }
            }
        }
    }

#ifdef _DEBUG
    for( FacetSet::iterator it1 = facets.begin(); it1 != facets.end(); ++it1 ) {
        for( int i = 0; i < Facet::DIMENSION; ++i )
            assert( it1->get_link_indexed( i ) != Facet::UNLINKED_VALUE() );
    }
#endif

    std::vector<std::pair<const typename Facet::RidgeType*, Facet*>> horizon;
    std::vector<Facet*> visible;
    std::vector<Facet*> new_facets;
    FacetMap new_ridges;

    // pfm.exit_section(initSection);
    // pfm.enter_section(loopSection);

    while( 1 ) {
        // pfm.enter_section(loopSearchSection);

        // A quick and dirty way I sped up the search routine was to have a hash_set of pointers to
        // facets that contained at least one outside point. We are done if there are no active facets.
        if( activeFacets.size() == 0 )
            break;

        typename FacetSet::iterator it = facets.find( **activeFacets.begin() );
        // FacetSet::iterator it = facets.begin();
        // while(it != facets.end() && it->num_outside() == 0)
        //	++it;

        // if(it == facets.end())
        //	break;

        int pt = it->furthest_outside();

        static unsigned SEARCH_ID = 0; // 0xFFFFFFFF;
        ++SEARCH_ID;

        // pfm.exit_section(loopSearchSection);
        // pfm.enter_section(loopHorizonSection);

        build_horizon<Traits>( horizon, visible, const_cast<typename FacetSet::value_type&>( *it ), pts, pt,
                               SEARCH_ID );
        assert( horizon.size() >= Facet::DIMENSION );

        // pfm.exit_section(loopHorizonSection);
        // pfm.enter_section(loopBuildSection);
        for( int i = 0; i < (int)horizon.size(); ++i ) {
            // pfm.enter_section(loopBuildInsertSection);

            typename FacetSet::iterator new_f =
                facets.insert( Traits::facet_from_ridge( *horizon[i].first, *horizon[i].second, pt ) ).first;

            // pfm.exit_section(loopBuildInsertSection);
            // pfm.enter_section(loopBuildUpdateSection);

            new_facets.push_back( const_cast<typename FacetSet::value_type*>( &*new_f ) );

            const_cast<typename FacetSet::value_type&>( *new_f ).set_link( *horizon[i].first, horizon[i].second );
            horizon[i].second->set_link( *horizon[i].first, const_cast<typename FacetSet::value_type*>( &*new_f ) );

            for( int j = 0; j < (int)Facet::DIMENSION; ++j ) {
                if( new_f->get_link_indexed( j ) == Facet::UNLINKED_VALUE() ) {
                    const typename Facet::RidgeType& r = new_f->get_ridge( j );
                    typename FacetMap::iterator map_it = new_ridges.find( r );

                    if( map_it != new_ridges.end() ) {
                        // Link the two facets that share a ridge then remove the active ridge from the map
                        const_cast<typename FacetSet::value_type&>( *new_f ).set_link_indexed( j, map_it->second );
                        map_it->second->set_link( r, const_cast<typename FacetSet::value_type*>( &*new_f ) );
                        new_ridges.erase( map_it );
                    } else
                        new_ridges[r] =
                            const_cast<typename FacetSet::value_type*>( &*new_f ); // Add a new active ridge to the map
                }
            }

            // pfm.exit_section(loopBuildUpdateSection);
        }

        assert( new_ridges.size() == 0 );
        new_ridges.clear();
        horizon.clear();

        // pfm.exit_section(loopBuildSection);
        // pfm.enter_section(loopDestroySection);

        // Collect and redistribute the outside pts from the visible set, then destory the facets in the visible set.
        for( int i = 0; i < (int)visible.size(); ++i ) {
            Facet& old_f = *visible[i];
            for( int j = 0; j < old_f.num_outside(); ++j ) {
                int p = old_f.get_outside( j );
                if( p == pt )
                    continue;

                for( int k = 0; k < (int)new_facets.size(); ++k ) {
                    Facet& new_f = *new_facets[k];
                    float d = Traits::distance_to_facet( pts, new_f, p );

                    if( d > 0 ) {
                        new_f.add_outside( d, p );
                        break;
                    }
                }
            }

            activeFacets.erase( &old_f );
            facets.erase( old_f );
        }

        // A quick and dirty way I sped up the search routine was to have a hash_set of pointers to
        // facets that contained at least one outside point.
        for( int k = 0; k < (int)new_facets.size(); ++k ) {
            if( new_facets[k]->num_outside() > 0 )
                activeFacets.insert( new_facets[k] );
        }

        visible.clear();
        new_facets.clear();

        // pfm.exit_section(loopDestroySection);
    } // while(1)

    // pfm.exit_section(loopSection);

    // Convert all the remaining facets into the output type and add to the hull
    for( typename FacetSet::iterator it1 = facets.begin(); it1 != facets.end(); ++it1 )
        hull.push_back( Traits::primitive_from_facet( *it1 ) );

    // std::cout << pfm << "\n";
}

#if defined( _MSC_VER )
#pragma warning( pop )
#endif
} // namespace detail
} // namespace geometry
} // namespace frantic
