// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/particles/particle_array.hpp>
#pragma warning( push )
#pragma warning( disable : 4100 )
#include <boost/random.hpp>
#pragma warning( pop )

#include <algorithm>

namespace frantic {
namespace particles {

/**
 *	A convenience function for checking whether the channels in given channel map are all
 *	doubles, and that the number of channels matches the given dimensionality.  If
 *	not, the function throws an exception.
 *
 *	@param	pcm		The particle channel map to check.
 *	@param	N		The dimensionality to check against.
 */
void float64_particle_check( const frantic::channels::channel_map& pcm, size_t N ) {
    // check to make sure that the particles contain doubles only
    for( size_t i = 0; i < pcm.channel_count(); ++i ) {
        if( pcm[i].data_type() != frantic::channels::data_type_float64 )
            throw std::runtime_error(
                "frantic::particles::float64_particle_check() - Particles are requried to consist only "
                "of float64 channels for the distance calculation." );
    }
    // check to make sure the dimensionality matches the size of the pcm
    if( float( N ) != pcm.structure_size() / 8.f ) {
        std::stringstream ss;
        ss << "frantic::particles::float64_particle_check() - The data dimensionality (" << N
           << ") does not match the particle structure size (" << pcm.structure_size()
           << ").  The particle channel map should contain N double values only, and be of size N*8 (8 bytes per "
              "double)."
           << std::endl;
        throw std::runtime_error( ss.str() );
    }
}

/**
 *	A convenience function for checking whether the channels in given channel map are all
 *	floats, and that the number of channels matches the given dimensionality.  If
 *	not, the function throws an exception.
 *
 *	@param	pcm		The particle channel map to check.
 *	@param	N		The dimensionality to check against.
 */
void float32_particle_check( const frantic::channels::channel_map& pcm, size_t N ) {
    // check to make sure that the particles contain doubles only
    for( size_t i = 0; i < pcm.channel_count(); ++i ) {
        if( pcm[i].data_type() != frantic::channels::data_type_float32 )
            throw std::runtime_error(
                "frantic::particles::float32_particle_check() - Particles are requried to consist only "
                "of float32 channels for the distance calculation." );
    }
    // check to make sure the dimensionality matches the size of the pcm
    if( float( N ) != pcm.structure_size() / 4.f ) {
        std::stringstream ss;
        ss << "frantic::particles::float32_particle_check() - The data dimensionality (" << N
           << ") does not match the particle structure size (" << pcm.structure_size()
           << ").  The particle channel map should contain N float values only, and be of size N*4 (4 bytes per float)."
           << std::endl;
        throw std::runtime_error( ss.str() );
    }
}

/**
 *	This function performs a kmeans clustering of the particle data in the
 *  given particle array.  Note that the particle data in each channel
 *  is assumed to be in generally the same scale, otherwise channels of
 *	a significantly larger scale will skew the clustering in favour
 *  of that channel data.  Particle data in the pa array must consist
 *	of doubles only.
 *
 *	@param	pa			An array of particles to cluster, all of dimension N.
 *	@param	N			The dimensionality of the data.
 *	@param	K			The number of clusters to determine.
 *	@param	outClusters	A return vector of cluster assignments for each particl in pa.
 *	@param	outKMeans	A return particle array of mean cluster vectors.  Values passed in
 *						this array will be used as initial cluster centres.  Note: initializing
 *						two clusters to the same centre is a degenerate case.
 *	@param	tol			Tolerance value for movement of cluster vectors to determine when to stop.
 *	@param	maxIters	Maximum number of iterations before terminating.
 */
void particle_kmeans( const frantic::particles::particle_array& pa, const size_t N, const size_t K,
                      std::vector<int>& outClusters, frantic::particles::particle_array& outKMeans, double tol = 1e-5f,
                      int maxIters = 100 ) {

    frantic::channels::channel_map pcm = pa.get_channel_map();

    // check for non-zero cluster count
    if( K <= 0 )
        throw std::runtime_error(
            "frantic::particles::particle_kmeans() - The number of requested cluster centres must be positive." );

    // check to make sure we have enough particles...
    if( K > pa.size() ) {
        std::stringstream ss;
        ss << "frantic::particles::particle_kmeans() - " << K << " cluster centres have been requested but only "
           << pa.size()
           << " particles have been provided.  There must be at least as many particles provided as the requested "
              "number "
              "of cluster centres."
           << std::endl;
        throw std::runtime_error( ss.str() );
    }

    // check to make sure that the particles contain doubles only
    float64_particle_check( pcm, N );

    // setup return and temp data
    size_t numParticles = pa.size();
    outClusters.resize( numParticles );
    std::vector<int> clusterCount( K );
    frantic::particles::particle_array outKMeansPrev( pcm );
    outKMeansPrev.copy_particles_from( outKMeans );

    // iteratively optimize cluster vector locations by modifying alternately the
    // particle to cluster membership, followed by the cluster mean vectors
    double maxTol, maxVal;
    int iters = 0;
    do {
        // for each particle, determine the distance from it to each cluster mean and
        // assign it to the closest one.
        double *p, *m, minVal, curVal;
        size_t minK = 0;
        for( size_t j = 0; j < numParticles; ++j ) {
            p = (double*)pa[j];

            minVal = std::numeric_limits<double>::max();
            for( size_t k = 0; k < K; ++k ) {
                m = (double*)outKMeans[k];
                curVal = 0;
                for( size_t n = 0; n < N; ++n ) {
                    curVal += pow( p[n] - m[n], 2 );
                }
                if( curVal < minVal ) {
                    minVal = curVal;
                    minK = k;
                }
            }
            outClusters[j] = (int)minK;
        }

        // for each cluster, average the locations of the member particles and reassign
        // the cluster mean to that location.  do a difference of the new mean locations
        // and the previous ones to check against the stopping criterion
        outKMeans.swap( outKMeansPrev );
        std::fill( clusterCount.begin(), clusterCount.end(), 0 );
        memset( outKMeans[0], 0, outKMeans.size_in_memory() );
        for( size_t i = 0; i < numParticles; ++i ) {
            p = (double*)pa[i];
            m = (double*)outKMeans[outClusters[i]];

            ++clusterCount[outClusters[i]];
            for( size_t n = 0; n < N; ++n ) {
                m[n] += p[n];
            }
        }

        double curTol, *mPrev;
        maxVal = -std::numeric_limits<double>::max();
        maxTol = -std::numeric_limits<double>::max();
        for( size_t k = 0; k < K; ++k ) {
            curTol = 0.f;
            m = (double*)outKMeans[k];
            mPrev = (double*)outKMeansPrev[k];
            for( size_t n = 0; n < N; ++n ) {
                if( clusterCount[k] > 0 )
                    m[n] /= clusterCount[k];
                else
                    m[n] = mPrev[n];
                curTol += pow( m[n] - mPrev[n], 2 );
                maxVal = std::max( maxVal, m[n] );
            }
            maxTol = std::max( maxTol, sqrt( curTol ) );
        }
    } while( maxTol / maxVal > tol && ++iters < maxIters );
}

} // namespace particles
}; // namespace frantic
