// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace frantic {
namespace rendering {
// using frantic::graphics::vector3f;

class shadow_generator {
  public:
    virtual ~shadow_generator() {}

    float compute_attenuation( const frantic::graphics::vector3f& /*start*/,
                               const frantic::graphics::vector3f& /*end*/ ) {
        return 1;
    }
};

// The shadow combiner class has three levels of shadow complexity.
//   Geometry shadows feature quick shadow computation and have a good chance of causing complete
//       attenuation.
//   Fast volumetric shadows feature quick shadow computation. An example of this is a constant density
//       volumetric.
//   Volumetric shadows are everything else.  An example of this is a varying density volumetric with
//       a raymarching shadow calculation.
class shadow_combiner {
    std::vector<boost::shared_ptr<shadow_generator>> m_geometryShadows;
    std::vector<boost::shared_ptr<shadow_generator>> m_fastVolumetricshadows;
    std::vector<boost::shared_ptr<shadow_generator>> m_volumetricShadows;

  public:
    void add_geometry_shadow( const boost::shared_ptr<shadow_generator>& shadow ) {
        m_geometryShadows.push_back( shadow );
    }

    void add_fast_volumetric_shadow( const boost::shared_ptr<shadow_generator>& shadow ) {
        m_fastVolumetricshadows.push_back( shadow );
    }

    void add_volumetric_shadow( const boost::shared_ptr<shadow_generator>& shadow ) {
        m_volumetricShadows.push_back( shadow );
    }

    float compute_attenuation( const frantic::graphics::vector3f& start,
                               const frantic::graphics::vector3f& end ) const {
        float attenuation = 1;
        std::vector<boost::shared_ptr<shadow_generator>>::const_iterator i;
        for( i = m_geometryShadows.begin(); i != m_geometryShadows.end(); ++i ) {
            attenuation *= ( *i )->compute_attenuation( start, end );
            if( attenuation == 0 )
                return 0;
        }
        for( i = m_fastVolumetricshadows.begin(); i != m_fastVolumetricshadows.end(); ++i ) {
            attenuation *= ( *i )->compute_attenuation( start, end );
            if( attenuation == 0 )
                return 0;
        }
        for( i = m_volumetricShadows.begin(); i != m_volumetricShadows.end(); ++i ) {
            attenuation *= ( *i )->compute_attenuation( start, end );
            if( attenuation == 0 )
                return 0;
        }
        return attenuation;
    }
};

} // namespace rendering
} // namespace frantic
