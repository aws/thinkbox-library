// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace particles {
namespace sph {

/**
 * "Poly 6" kernel function used in SPH
 */
class poly6_kernel {
  private:
    float m_poly6Scale;
    float m_poly6DerivScale;
    float m_compactSupportSquared;

  public:
    poly6_kernel() {}
    poly6_kernel( float compactSupport ) {
        // formula for the poly6 kernel function is m_poly6Scale * ( m_compactSupportSquared - distanceSquared )^3. (see
        // papers on the topic)
        m_poly6Scale = 315.0f / ( 64.0f * (float)M_PI * powf( compactSupport, 9.0f ) );
        m_poly6DerivScale = 945.0f / ( 32.0f * (float)M_PI * powf( compactSupport, 9.0f ) );
        m_compactSupportSquared = compactSupport * compactSupport;
    }

    float kernel( float distanceSquared ) {
        float difference = m_compactSupportSquared - distanceSquared;
        if( difference < 0.0f )
            return 0.0f;
        return m_poly6Scale * difference * difference * difference;
    }

    frantic::graphics::vector3f kernel_gradient( frantic::graphics::vector3f distanceVector ) {
        float difference = m_compactSupportSquared - distanceVector.get_magnitude_squared();
        if( difference < 0.0f )
            return frantic::graphics::vector3f( 0.0f );
        return ( m_poly6DerivScale * difference * difference ) * distanceVector;
    }
};

/**
 * "Spiky" kernel function used in SPH
 */
class spiky_kernel {
  private:
    float m_spikyScale;
    float m_spikyDerivScale;
    float m_compactSupport;

  public:
    spiky_kernel() {}
    spiky_kernel( float compactSupport ) {
        m_spikyScale = 15.0f / ( (float)M_PI * powf( compactSupport, 6.0f ) );
        m_spikyDerivScale = -45.0f / ( (float)M_PI * powf( compactSupport, 6.0f ) );
        m_compactSupport = compactSupport;
    }

    float kernel( float distance ) {
        float difference = m_compactSupport - distance;
        if( difference < 0.0f )
            return 0.0f;
        return m_spikyScale * difference * difference * difference;
    }

    frantic::graphics::vector3f kernel_gradient( frantic::graphics::vector3f distanceVector ) {
        float difference = m_compactSupport - distanceVector.get_magnitude();
        if( difference < 0.0f )
            return frantic::graphics::vector3f( 0.0f );
        return ( m_spikyDerivScale * difference * difference ) * distanceVector;
    }
};

} // namespace sph
} // namespace particles
} // namespace frantic
