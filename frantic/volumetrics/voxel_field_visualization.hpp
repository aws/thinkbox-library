// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0


#pragma once

#include <vector>

#include <boost/function.hpp>

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/graphics/boundbox3.hpp>
#include <frantic/math/utils.hpp>

#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/volumetrics/levelset/rle_defined_and_adj_iterator.hpp>

// this is totally wrong... perhaps the debug mesh functionality should be moved into fluids ?
// the reason to include it is the debug mesh policy uses face state masks that included throug this header
#include <frantic/fluids/boundary_condition.hpp> // TODO this feels wrong to be including a "fluids" header here...something is wrong
//#include <frantic/fluids/variational_solid_matrix.hpp>

namespace frantic {
namespace volumetrics {
namespace visualization {

/**
 * This color policy class specifies an interface that can be used return a set of vertex colors based on channel data
 */
class debug_mesh_color_policy {
  public:
    debug_mesh_color_policy() {}

    virtual bool
    compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                           frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                           frantic::graphics::color3f* outColors ) = 0;
};

/**
 * This color policy class specifies a coloring suitable for viewing a staggered populated channel
 */
class debug_mesh_populated_color_policy : public debug_mesh_color_policy {
  public:
    debug_mesh_populated_color_policy() {}

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {

        const boost::uint8_t populated =
            *reinterpret_cast<const boost::uint8_t*>( srcChannelAccessor.data( iterator.get_center_data_index() ) );

        frantic::graphics::color3f col( 0.f );

        if( populated == 0 ) {
            col = frantic::graphics::color3f( 1.0f );
        } else if( populated == 1 ) {
            col = frantic::graphics::color3f( 0.f, 0.f, 1.f );
        } else if( populated == 2 ) {
            col = frantic::graphics::color3f( 1.0f, 0.f, 1.f );
        } else if( populated == 3 ) {
            col = frantic::graphics::color3f( 1.0f, 0.f, 0.f );
        }
        for( int v = 0; v < 24; ++v ) {
            outColors[v] = col;
        }

        return true;
    }
};

/**
 * This color policy class specifies a coloring suitable for viewing a staggered populated channel
 */
class debug_mesh_staggered_populated_color_policy : public debug_mesh_color_policy {
  public:
    debug_mesh_staggered_populated_color_policy() {}

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {

        const boost::uint8_t staggeredPopulated =
            *reinterpret_cast<const boost::uint8_t*>( srcChannelAccessor.data( iterator.get_center_data_index() ) );

        int reorderedFaces[6] = { 1, 0, 3, 2, 5, 4 };
        //			int invalidCount = 0;
        frantic::graphics::color3f col( 0.f );

        bool include = false;
        // color the 3 locally stored faces
        for( int face = 0; face < 6; face += 2 ) {

            int v = face * 4;
            const boost::uint8_t facePopulated = staggeredPopulated & ( boost::uint8_t )( 1 << ( face / 2 ) );

            if( facePopulated > 0 ) {
                include = true;
                col = frantic::graphics::color3f( 0.1f, 0.8f, 1.f );
            } else {
                col = frantic::graphics::color3f( 1.f, 0.8f, 0.1f );
            }

            for( size_t i = v; i != (size_t)v + 4; ++i ) {
                outColors[i] = col;
            }
        }

        // color using the adjacent faces
        for( int face = 1; face < 6; face += 2 ) {
            int v = face * 4;
            int adjIndex = iterator.get_adjacent_data_index( reorderedFaces[face] );

            if( adjIndex >= 0 ) {
                const boost::uint8_t staggeredAdjPopulated =
                    *reinterpret_cast<const boost::uint8_t*>( srcChannelAccessor.data( adjIndex ) );
                const boost::uint8_t facePopulated = staggeredAdjPopulated & ( boost::uint8_t )( 1 << ( face / 2 ) );

                if( facePopulated > 0 ) {
                    include = true;
                    col = frantic::graphics::color3f( 0.1f, 0.8f, 1.f );
                } else {
                    col = frantic::graphics::color3f( 1.f, 0.8f, 0.1f );
                }

            } else {
                col = frantic::graphics::color3f( 1.f, 0.8f, 0.1f );
            }

            for( size_t i = v; i != (size_t)v + 4; ++i ) {
                outColors[i] = col;
            }
        }

        return include;
    }
};

/**
 * This color policy class specifies a coloring suitable for viewing a face state channel
 */
class debug_mesh_face_state_color_policy : public debug_mesh_color_policy {
  public:
    debug_mesh_face_state_color_policy() {}

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        const boost::uint32_t faceValues =
            *reinterpret_cast<const boost::uint32_t*>( srcChannelAccessor.data( iterator.get_center_data_index() ) );

        int reorderedFaces[6] = { 1, 0, 3, 2, 5, 4 };
        int invalidCount = 0;
        frantic::graphics::color3f col( 0.f );

        for( int face = 0, v = 0; face < 6; ++face ) {
            boost::uint32_t value = ( faceValues & frantic::fluids::face_masks::get_mask( reorderedFaces[face] ) ) >>
                                    frantic::fluids::face_masks::get_shift( reorderedFaces[face] );
            switch( value ) {
            case frantic::fluids::INVALID:
                col = frantic::graphics::color3f( 0.f, 0.8f, 1.f );
                ++invalidCount;
                break;
            case frantic::fluids::NONE:
                col = frantic::graphics::color3f( 0.f, 0.f, 1.f );
                break;
            case frantic::fluids::DIRICHLET:
                col = frantic::graphics::color3f( 0.f, 1.f, 0.f );
                break;
            case frantic::fluids::NEUMANN:
                col = frantic::graphics::color3f( 1.f, 0.f, 0.f );
                break;
            default:
                throw std::runtime_error(
                    "For voxel " + iterator.get_coord().str() + " the boundary condition value for face " +
                    boost::lexical_cast<std::string>( face ) + " is not one of the set values (0,1,2,4) it is " +
                    boost::lexical_cast<std::string>( value ) );
            }

            for( size_t i = v; i != (size_t)v + 4; ++i ) {
                outColors[i] = col;
            }
            v += 4;
        }

        if( invalidCount != 6 )
            return true;
        else
            return false;
    }
};

/**
 * This color policy class specifies a coloring suitable for viewing a pressure update flag channel
 *
 * The coloring follows that used in debug_mesh_face_state_color_policy:
 *
 *   none - invisible
 *   ghost cell - green
 *   air fluid cell - cyan
 *   fluid cell - blue
 *   other - black
 *
 */
class debug_mesh_pressure_update_flag_color_policy : public debug_mesh_color_policy {
  public:
    debug_mesh_pressure_update_flag_color_policy() {}

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        const boost::uint8_t value =
            *reinterpret_cast<const boost::uint8_t*>( srcChannelAccessor.data( iterator.get_center_data_index() ) );

        // int reorderedFaces[6] = {1,0,3,2,5,4};
        frantic::graphics::color3f col( 0.f );

        switch( value ) {
        case 0: // frantic::fluids::pressure_solve_none:
            return false;
        case 1: // frantic::fluids::pressure_solve_free_surface_ghost:
        case 2: // frantic::fluids::pressure_solve_free_surface_discrete:
            col = frantic::graphics::color3f( 0, 1.f, 0 );
            break;
        case 3: // frantic::fluids::pressure_solve_air:
            col = frantic::graphics::color3f( 1.f, 0, 0 );
            break;
        case 4: // frantic::fluids::pressure_solve_fluid:
            col = frantic::graphics::color3f( 0, 0, 1.f );
            break;
        case 5: // frantic::fluids::pressure_solve_dirichlet:
            col = frantic::graphics::color3f( 1.f, 0, 1.f );
            break;
        case 6: // frantic::fluids::pressure_solve_neumann:
            col = frantic::graphics::color3f( 0.6f, 0.6f, 0.6f );
            break;
        }
        for( std::size_t i = 0; i < 24; ++i ) {
            outColors[i] = col;
        }
        return true;
    }
};

/**
 *  A color policy intended to visualize linear unsigned data
 * in a channel of arity 1.
 *
 *  This policy isn't strictly for linear unsigned data,
 * but zero and positive and negative numbers don't get
 * any special treatment in the color policy.
 *
 *  The color scheme is similar to ( or maybe the same as ? )
 * the default colour scheme used by Matlab.
 *
 * ( small values ) blue - red ( large values )
 */
class debug_mesh_linear_unsigned_color_policy : public debug_mesh_color_policy {
  private:
    double m_minValue, m_maxValue, m_deltaValue, m_invDeltaValue;
    debug_mesh_linear_unsigned_color_policy() {}

    void init( const double minValue, const double maxValue ) {
        m_minValue = minValue;
        m_maxValue = maxValue;

        if( m_minValue > m_maxValue ) {
            std::swap( m_minValue, m_maxValue );
        }

        if( m_minValue > 0 )
            m_minValue = 0;

        m_deltaValue = m_maxValue - m_minValue;

        if( m_deltaValue > 0 )
            m_invDeltaValue = 1.0 / m_deltaValue;
    }

    double
    get_channel_value_as_double( frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                 const std::size_t i ) {
        if( srcChannelAccessor.arity() != 1 ) {
            throw std::runtime_error( "debug_mesh_linear_unsigned_color_policy::get_channel_value_as_double Error: the "
                                      "channel arity must be 1, but instead the channel is " +
                                      frantic::strings::to_string( frantic::channels::channel_data_type_str(
                                          srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }
        switch( srcChannelAccessor.data_type() ) {
        case frantic::channels::data_type_int8:
            return static_cast<double>( *reinterpret_cast<const boost::int8_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_int16:
            return static_cast<double>( *reinterpret_cast<const boost::int16_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_int32:
            return static_cast<double>( *reinterpret_cast<const boost::int32_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_int64:
            return static_cast<double>( *reinterpret_cast<const boost::int64_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_uint8:
            return static_cast<double>( *reinterpret_cast<const boost::uint8_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_uint16:
            return static_cast<double>( *reinterpret_cast<const boost::uint16_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_uint32:
            return static_cast<double>( *reinterpret_cast<const boost::uint32_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_uint64:
            return static_cast<double>( *reinterpret_cast<const boost::uint64_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_float16:
            return static_cast<double>( *reinterpret_cast<const half*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_float32:
            return static_cast<double>( *reinterpret_cast<const float*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_float64:
            return static_cast<double>( *reinterpret_cast<const double*>( srcChannelAccessor.data( i ) ) );
        default:
            throw std::runtime_error(
                "debug_mesh_linear_unsigned_color_policy::get_channel_value_as_double Error: the "
                "channel type must have a numeric data type and arity 1, but instead the channel is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }
    }

  public:
    /**
     *  The color policy is scaled to accomodate numbers minValue
     * through maxValue.
     */
    debug_mesh_linear_unsigned_color_policy( const double minValue, const double maxValue ) {
        init( minValue, maxValue );
    }

    /**
     *  The color policy is scaled to accomodate the minimum and
     * maximum values in the channel.
     */
    debug_mesh_linear_unsigned_color_policy(
        frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor ) {
        if( srcChannelAccessor.arity() != 1 ) {
            throw std::runtime_error(
                "debug_mesh_linear_unsigned_color_policy Error: the channel type must be float32, but instead it is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }

        double minValue = 0;
        double maxValue = 0;

        for( std::size_t i = 0; i < srcChannelAccessor.size(); ++i ) {
            const double value = get_channel_value_as_double( srcChannelAccessor, i );

            if( i == 0 ) {
                minValue = value;
                maxValue = value;
            }

            if( value < minValue )
                minValue = value;
            if( value > maxValue )
                maxValue = value;
        }

        init( minValue, maxValue );
    }

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        if( srcChannelAccessor.arity() != 1 ) {
            throw std::runtime_error(
                "debug_mesh_linear_unsigned_color_policy Error: the channel type must be float32, but instead it is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }

        const double value = get_channel_value_as_double( srcChannelAccessor, iterator.get_center_data_index() );

        frantic::graphics::color3f col( 0.f );

        if( m_deltaValue == 0 ) {
            col = frantic::graphics::color3f( 0, 0, 1.f );
        } else {
            const float valueFraction =
                frantic::math::clamp( static_cast<float>( ( value - m_minValue ) * m_invDeltaValue ), 0.f, 1.f );

            float r = 1.f;
            float g = 1.f;
            float b = 1.f;

            switch( static_cast<int>( floor( 8.f * valueFraction ) ) ) {
            case 0: // [0,1/8)
                b = 0.5f + 4.f * valueFraction;
                g = 0;
                r = 0;
                break;
            case 1: // [1/8,3/8)
            case 2:
                b = 1.f;
                g = -0.5f + 4.f * valueFraction;
                r = 0;
                break;
            case 3: // [3/8,5/8)
            case 4:
                b = 2.5f - 4.f * valueFraction;
                g = 1.f;
                r = -1.5f + 4.f * valueFraction;
                break;
            case 5: // [5/8,6/8)
            case 6:
                b = 0;
                g = 3.5f - 4.f * valueFraction;
                r = 1.f;
                break;
            case 7: // [7/8,1]
            case 8:
                b = 0;
                g = 0;
                r = 4.5f - 4.f * valueFraction;
                break;
            }

            col = frantic::graphics::color3f( r, g, b );
        }

        for( std::size_t i = 0; i < 24; ++i ) {
            outColors[i] = col;
        }

        return true;
    }
};

/**
 *  A color policy intended to visualize linear signed data
 * in a channel of arity 1.
 *
 * negative values: blue
 * zero: black
 * positive values: red
 */
class debug_mesh_linear_signed_color_policy : public debug_mesh_color_policy {
  private:
    double m_maxMagnitude, m_invMaxMagnitude;
    boost::function<bool( double )> m_displayCondition;

    debug_mesh_linear_signed_color_policy() {}
    void init( const double maxMagnitude ) {
        m_maxMagnitude = maxMagnitude;
        if( m_maxMagnitude != 0 )
            m_invMaxMagnitude = 1.0 / maxMagnitude;
    }

    double
    get_channel_value_as_double( frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                 const std::size_t i ) {
        if( srcChannelAccessor.arity() != 1 ) {
            throw std::runtime_error(
                "debug_mesh_linear_signed_color_policy::get_channel_value_as_double Error: the channel "
                "arity must be 1, but instead the channel is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }
        switch( srcChannelAccessor.data_type() ) {
        case frantic::channels::data_type_int8:
            return static_cast<double>( *reinterpret_cast<const boost::int8_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_int16:
            return static_cast<double>( *reinterpret_cast<const boost::int16_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_int32:
            return static_cast<double>( *reinterpret_cast<const boost::int32_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_int64:
            return static_cast<double>( *reinterpret_cast<const boost::int64_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_uint8:
            return static_cast<double>( *reinterpret_cast<const boost::uint8_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_uint16:
            return static_cast<double>( *reinterpret_cast<const boost::uint16_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_uint32:
            return static_cast<double>( *reinterpret_cast<const boost::uint32_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_uint64:
            return static_cast<double>( *reinterpret_cast<const boost::uint64_t*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_float16:
            return static_cast<double>( *reinterpret_cast<const half*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_float32:
            return static_cast<double>( *reinterpret_cast<const float*>( srcChannelAccessor.data( i ) ) );
        case frantic::channels::data_type_float64:
            return static_cast<double>( *reinterpret_cast<const double*>( srcChannelAccessor.data( i ) ) );
        default:
            throw std::runtime_error(
                "debug_mesh_linear_signed_color_policy::get_channel_value_as_double Error: the channel "
                "type must have a numeric data type and arity 1, but instead the channel is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }
    }

  public:
    /**
     *  The color policy is scaled to accomodate numbers minValue
     * through maxValue.
     */
    debug_mesh_linear_signed_color_policy( const double maxMagnitude,
                                           boost::function<bool( double )> displayCondition = 0 )
        : m_displayCondition( displayCondition ) {
        init( maxMagnitude );
    }

    /**
     *  The color policy is scaled to accomodate the minimum and
     * maximum values in the channel.
     */
    debug_mesh_linear_signed_color_policy(
        frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
        boost::function<bool( double )> displayCondition = 0 )
        : m_displayCondition( displayCondition ) {
        if( srcChannelAccessor.arity() != 1 ) {
            throw std::runtime_error( "debug_mesh_linear_signed_color_policy Error: the channel arity must be 1, but "
                                      "instead the channel is " +
                                      frantic::strings::to_string( frantic::channels::channel_data_type_str(
                                          srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }

        double maxMagnitude = 0;

        for( std::size_t i = 0; i < srcChannelAccessor.size(); ++i ) {
            // const double value = static_cast<double>( *reinterpret_cast<const float*>( srcChannelAccessor.data( i ) )
            // );
            const double absValue = fabs( get_channel_value_as_double( srcChannelAccessor, i ) );

            if( i == 0 ) {
                maxMagnitude = absValue;
            }

            if( absValue > maxMagnitude )
                maxMagnitude = absValue;
        }

        init( maxMagnitude );
    }

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        if( srcChannelAccessor.arity() != 1 ) {
            throw std::runtime_error( "debug_mesh_linear_signed_color_policy Error: the channel arity must be 1, but "
                                      "instead the channel is " +
                                      frantic::strings::to_string( frantic::channels::channel_data_type_str(
                                          srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }

        const double value = get_channel_value_as_double( srcChannelAccessor, iterator.get_center_data_index() );

        if( !m_displayCondition.empty() && !m_displayCondition( value ) )
            return false;

        frantic::graphics::color3f col( 0.f );

        if( m_maxMagnitude == 0 ) {
            col = frantic::graphics::color3f( 0, 0, 0 );
        } else {
            // [-maxMagnitude,+maxMagnitude] -> [0,1]
            const float valueFraction =
                frantic::math::clamp( static_cast<float>( 0.5 * value * m_invMaxMagnitude + 0.5 ), 0.f, 1.f );

            /*
            const float r = frantic::math::clamp( -1.f + 5.f * fabs( valueFraction - 0.3f ), 0.f, 1.f );
            const float g = frantic::math::clamp( -1.f + 5.f * fabs( valueFraction - 0.5f ), 0.f, 1.f );
            const float b = frantic::math::clamp( -1.f + 5.f * fabs( valueFraction - 0.7f ), 0.f, 1.f );
            */

            const float r = ( valueFraction > 0.25f && valueFraction <= 0.5f )
                                ? 0
                                : frantic::math::clamp( -0.5f + 4.f * fabsf( valueFraction - 0.25f ), 0.f, 1.f );
            const float g = frantic::math::clamp( -0.5f + 4.f * fabsf( valueFraction - 0.5f ), 0.f, 1.f );
            const float b = ( valueFraction >= 0.5f && valueFraction < 0.75f )
                                ? 0
                                : frantic::math::clamp( -0.5f + 4.f * fabsf( valueFraction - 0.75f ), 0.f, 1.f );

            col = frantic::graphics::color3f( r, g, b );
        }

        for( std::size_t i = 0; i < 24; ++i ) {
            outColors[i] = col;
        }

        return true;
    }
};

/**
 * This color policy class specifies a coloring suitable for a vector3f channel
 */
class debug_mesh_staggered_linear_signed_color_policy : public debug_mesh_color_policy {
    double m_maxMagnitude, m_invMaxMagnitude;

    debug_mesh_staggered_linear_signed_color_policy() {}

    void init( const double maxMagnitude ) {
        m_maxMagnitude = maxMagnitude;
        if( m_maxMagnitude != 0 )
            m_invMaxMagnitude = 1.0 / maxMagnitude;
    }

    frantic::graphics::vector3f get_channel_value_as_vector3f(
        frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor, const std::size_t i ) {
        if( srcChannelAccessor.arity() != 3 ) {
            throw std::runtime_error(
                "debug_mesh_staggered_linear_signed_color_policy::get_channel_value_as_vector3f Error: "
                "the channel arity must be 3, but instead the channel is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }

        frantic::channels::channel_type_convertor_function_t convert = get_channel_type_convertor_function(
            srcChannelAccessor.data_type(), frantic::channels::channel_data_type_traits<float>::data_type(),
            _T("debug_mesh_staggered_linear_signed_color_policy") );

        frantic::graphics::vector3f value;
        convert( reinterpret_cast<char*>( &value[0] ), srcChannelAccessor.data( i ), 3 );

        return value;
    }

  public:
    /**
     *  The color policy is scaled to accomodate numbers minValue
     * through maxValue.
     */
    debug_mesh_staggered_linear_signed_color_policy( const double maxMagnitude ) { init( maxMagnitude ); }

    /**
     *  The color policy is scaled to accomodate the minimum and
     * maximum values in the channel.
     */
    debug_mesh_staggered_linear_signed_color_policy(
        frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor ) {
        if( srcChannelAccessor.arity() != 3 ) {
            throw std::runtime_error(
                "debug_mesh_staggered_linear_signed_color_policy Error: the channel must be of arity "
                "3, but instead the channel is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }

        double maxMagnitude = 0;

        for( std::size_t i = 0; i < srcChannelAccessor.size(); ++i ) {
            // const double value = static_cast<double>( *reinterpret_cast<const float*>( srcChannelAccessor.data( i ) )
            // );
            const frantic::graphics::vector3f value = get_channel_value_as_vector3f( srcChannelAccessor, i );
            const double absValue = value.max_abs_component();

            if( i == 0 ) {
                maxMagnitude = absValue;
            }

            if( absValue > maxMagnitude )
                maxMagnitude = absValue;
        }

        init( maxMagnitude );
    }

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        int reorderedFaces[6] = { 1, 0, 3, 2, 5, 4 };
        //			int invalidCount = 0;
        frantic::graphics::color3f col( 0.f );

        const bool include = true;
        // color the 3 locally stored faces
        for( int meshFace = 0; meshFace < 6; ++meshFace ) {
            const int voxelFace = reorderedFaces[meshFace];
            const int v = meshFace * 4;

            const boost::int32_t dataIndex =
                frantic::volumetrics::levelset::is_neighbor_index_direction_positive( voxelFace )
                    ? iterator.get_adjacent_data_index( voxelFace )
                    : iterator.get_center_data_index();

            if( dataIndex >= 0 ) {
                const int axis = frantic::volumetrics::levelset::neighbor_index_axis( voxelFace );
                const float faceValue = ( get_channel_value_as_vector3f( srcChannelAccessor, dataIndex ) )[axis];

                if( m_maxMagnitude == 0 ) {
                    col = frantic::graphics::color3f( 0, 0, 0 );
                } else {
                    // [-maxMagnitude,+maxMagnitude] -> [0,1]
                    const float valueFraction = frantic::math::clamp(
                        static_cast<float>( 0.5 * faceValue * m_invMaxMagnitude + 0.5 ), 0.f, 1.f );

                    /*
                    const float r = frantic::math::clamp( -1.f + 5.f * fabs( valueFraction - 0.3f ), 0.f, 1.f );
                    const float g = frantic::math::clamp( -1.f + 5.f * fabs( valueFraction - 0.5f ), 0.f, 1.f );
                    const float b = frantic::math::clamp( -1.f + 5.f * fabs( valueFraction - 0.7f ), 0.f, 1.f );
                    */

                    const float r =
                        ( valueFraction > 0.25f && valueFraction <= 0.5f )
                            ? 0
                            : frantic::math::clamp( -0.5f + 4.f * fabsf( valueFraction - 0.25f ), 0.f, 1.f );
                    const float g = frantic::math::clamp( -0.5f + 4.f * fabsf( valueFraction - 0.5f ), 0.f, 1.f );
                    const float b =
                        ( valueFraction >= 0.5f && valueFraction < 0.75f )
                            ? 0
                            : frantic::math::clamp( -0.5f + 4.f * fabsf( valueFraction - 0.75f ), 0.f, 1.f );

                    col = frantic::graphics::color3f( r, g, b );
                }
            } else {
                col = frantic::graphics::color3f( 1.f, 1.f, 1.f ); // white - none
            }

            for( size_t i = v; i != (size_t)v + 4; ++i ) {
                outColors[i] = col;
            }
        }

        return include;
    }
};

/**
 * This color policy class specifies a coloring suitable for a vector3f channel
 * with unsigned values.
 */
class debug_mesh_staggered_linear_unsigned_color_policy : public debug_mesh_color_policy {
    debug_mesh_staggered_linear_unsigned_color_policy() {}

    double m_minValue, m_maxValue, m_deltaValue, m_invDeltaValue;

    void init( const double minValue, const double maxValue ) {
        m_minValue = minValue;
        m_maxValue = maxValue;

        if( m_minValue > m_maxValue ) {
            std::swap( m_minValue, m_maxValue );
        }

        if( m_minValue > 0 )
            m_minValue = 0;

        m_deltaValue = m_maxValue - m_minValue;

        if( m_deltaValue > 0 )
            m_invDeltaValue = 1.0 / m_deltaValue;
    }

    frantic::graphics::vector3f get_channel_value_as_vector3f(
        frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor, const std::size_t i ) {
        if( srcChannelAccessor.arity() != 3 ) {
            throw std::runtime_error(
                "debug_mesh_staggered_linear_unsigned_color_policy::get_channel_value_as_vector3f "
                "Error: the channel arity must be 3, but instead the channel is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }

        frantic::channels::channel_type_convertor_function_t convert = get_channel_type_convertor_function(
            srcChannelAccessor.data_type(), frantic::channels::channel_data_type_traits<float>::data_type(),
            _T("debug_mesh_staggered_linear_unsigned_color_policy") );

        frantic::graphics::vector3f value;
        convert( reinterpret_cast<char*>( &value[0] ), srcChannelAccessor.data( i ), 3 );

        return value;
    }

  public:
    /**
     *  The color policy is scaled to accomodate numbers minValue
     * through maxValue.
     */
    debug_mesh_staggered_linear_unsigned_color_policy( const double minValue, const double maxValue ) {
        init( minValue, maxValue );
    }

    /**
     *  The color policy is scaled to accomodate the minimum and
     * maximum values in the channel.
     */
    debug_mesh_staggered_linear_unsigned_color_policy(
        frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor ) {
        if( srcChannelAccessor.arity() != 3 ) {
            throw std::runtime_error(
                "debug_mesh_staggered_linear_unsigned_color_policy Error: the channel must be of arity "
                "3, but instead the channel is " +
                frantic::strings::to_string( frantic::channels::channel_data_type_str(
                    srcChannelAccessor.arity(), srcChannelAccessor.data_type() ) ) );
        }

        double globalMinValue = 0;
        double globalMaxValue = 0;

        for( std::size_t i = 0; i < srcChannelAccessor.size(); ++i ) {
            // const double value = static_cast<double>( *reinterpret_cast<const float*>( srcChannelAccessor.data( i ) )
            // );
            const frantic::graphics::vector3f value = get_channel_value_as_vector3f( srcChannelAccessor, i );
            const double minValue = std::min( value.x, std::min( value.y, value.z ) );
            const double maxValue = std::max( value.x, std::max( value.y, value.z ) );

            if( i == 0 ) {
                globalMaxValue = maxValue;
                globalMaxValue = minValue;
            }

            if( maxValue > globalMaxValue )
                globalMaxValue = maxValue;
            if( minValue < globalMinValue )
                globalMinValue = minValue;
        }

        init( globalMinValue, globalMaxValue );
    }

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        int reorderedFaces[6] = { 1, 0, 3, 2, 5, 4 };
        //			int invalidCount = 0;
        frantic::graphics::color3f col( 0.f );

        const bool include = true;
        // color the 3 locally stored faces
        for( int meshFace = 0; meshFace < 6; ++meshFace ) {
            const int voxelFace = reorderedFaces[meshFace];
            const int v = meshFace * 4;

            const boost::int32_t dataIndex =
                frantic::volumetrics::levelset::is_neighbor_index_direction_positive( voxelFace )
                    ? iterator.get_adjacent_data_index( voxelFace )
                    : iterator.get_center_data_index();

            if( dataIndex >= 0 ) {
                const int axis = frantic::volumetrics::levelset::neighbor_index_axis( voxelFace );
                const float faceValue = ( get_channel_value_as_vector3f( srcChannelAccessor, dataIndex ) )[axis];

                if( m_deltaValue == 0 ) {
                    col = frantic::graphics::color3f( 0, 0, 0 );
                } else {
                    // [-maxMagnitude,+maxMagnitude] -> [0,1]
                    const float valueFraction = frantic::math::clamp(
                        static_cast<float>( ( faceValue - m_minValue ) * m_invDeltaValue ), 0.f, 1.f );

                    float r = 1.f;
                    float g = 1.f;
                    float b = 1.f;

                    switch( static_cast<int>( floor( 8.f * valueFraction ) ) ) {
                    case 0: // [0,1/8)
                        b = 0.5f + 4.f * valueFraction;
                        g = 0;
                        r = 0;
                        break;
                    case 1: // [1/8,3/8)
                    case 2:
                        b = 1.f;
                        g = -0.5f + 4.f * valueFraction;
                        r = 0;
                        break;
                    case 3: // [3/8,5/8)
                    case 4:
                        b = 2.5f - 4.f * valueFraction;
                        g = 1.f;
                        r = -1.5f + 4.f * valueFraction;
                        break;
                    case 5: // [5/8,6/8)
                    case 6:
                        b = 0;
                        g = 3.5f - 4.f * valueFraction;
                        r = 1.f;
                        break;
                    case 7: // [7/8,1]
                    case 8:
                        b = 0;
                        g = 0;
                        r = 4.5f - 4.f * valueFraction;
                        break;
                    }

                    col = frantic::graphics::color3f( r, g, b );
                }
            } else {
                col = frantic::graphics::color3f( 1.f, 1.f, 1.f ); // white - none
            }

            for( size_t i = v; i != (size_t)v + 4; ++i ) {
                outColors[i] = col;
            }
        }

        return include;
    }
};

class debug_mesh_is_non_zero_color_policy : public debug_mesh_color_policy {
  public:
    debug_mesh_is_non_zero_color_policy() {}

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        bool isZero = true;

        const char* data = srcChannelAccessor.data( iterator.get_center_data_index() );
        for( size_t i = 0; i < srcChannelAccessor.primitive_size(); ++i ) {
            if( data[i] != 0 ) {
                isZero = false;
                break;
            }
        }

        if( isZero ) {
            return false;
        }

        const frantic::graphics::color3f col( 0, 0, 1.f );

        for( std::size_t i = 0; i < 24; ++i ) {
            outColors[i] = col;
        }

        return true;
    }
};

class debug_mesh_is_zero_color_policy : public debug_mesh_color_policy {
  public:
    debug_mesh_is_zero_color_policy() {}

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        bool isZero = true;

        const char* data = srcChannelAccessor.data( iterator.get_center_data_index() );
        for( size_t i = 0; i < srcChannelAccessor.size(); ++i ) {
            if( data[i] != 0 ) {
                isZero = false;
                break;
            }
        }

        if( !isZero ) {
            return false;
        }

        const frantic::graphics::color3f col( 0, 0, 1.f );

        for( std::size_t i = 0; i < 24; ++i ) {
            outColors[i] = col;
        }

        return true;
    }
};

/**
 * This color policy class is intended to display the sign of float32[1]
 * data
 */
class debug_mesh_float_sign_color_policy : public debug_mesh_color_policy {
  public:
    debug_mesh_float_sign_color_policy() {}

    bool compute_vertex_colors( const frantic::volumetrics::levelset::rle_defined_and_adj_iterator& iterator,
                                frantic::volumetrics::levelset::const_rle_channel_general_accessor& srcChannelAccessor,
                                frantic::graphics::color3f* outColors ) {
        const float value =
            *reinterpret_cast<const float*>( srcChannelAccessor.data( iterator.get_center_data_index() ) );

        frantic::graphics::color3f col( 0.f );

        if( value < 0 ) {
            col = frantic::graphics::color3f( 0, 0, 1.f );
        } else if( value == 0 ) {
            col = frantic::graphics::color3f( 0, 0, 0 );
        } else { // ( value > 0 )
            col = frantic::graphics::color3f( 1.f, 0, 0 );
        }

        for( std::size_t i = 0; i < 24; ++i ) {
            outColors[i] = col;
        }

        return true;
    }
};

/**
 * This function is similar to the RLE debug meshing function. It will use the rle_index_spec to create a
 * mesh with a cube for every voxel and color the cubes based on the supplied color policy. This function
 * can be used to create a debug mesh for a given channel by simply implementing the appropriate color policy.
 *
 *@param vcs the voxel coord system of the containing field (or level set)
 *@param ris the rle_index_spec of the field
 *@param cubeSize the scaling size for the cubes
 *@param channelAcc the channel accessor to the data to use to color the cubes
 *@param colorPolicy the color policy that is called to color the cubes
 *@param outMesh the resulting mesh
 */
void convert_voxel_channel_to_debug_mesh(
    const frantic::volumetrics::voxel_coord_system& vcs, const frantic::volumetrics::levelset::rle_index_spec& ris,
    float cubeSize, frantic::volumetrics::levelset::const_rle_channel_general_accessor& channelAcc,
    debug_mesh_color_policy& colorPolicy, frantic::geometry::trimesh3& outMesh );

} // namespace visualization
} // namespace volumetrics
} // namespace frantic
