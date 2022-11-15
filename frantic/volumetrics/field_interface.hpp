// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_map.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frantic {
namespace volumetrics {

class field_interface {
  public:
    virtual ~field_interface() {}

    /**
     * Evaluates the field at the specified location, storing the result in the supplied buffer. The buffer must adhere
     * to the layout described by the channel_map retrieved via get_channel_map().
     *
     * @note This member function is const in order to suggest a thread-safe approach to its implementation.
     *
     * @param dest The results of evaluating the field at 'pos' will be stored in this buffer.
     * @param pos The location to evaluate the field at. This location is local to the field and has no global meaning.
     * @return Returns false if the field is undefined at the given location.
     */
    virtual bool evaluate_field( void* dest, const frantic::graphics::vector3f& pos ) const = 0;

    /**
     * Returns the description of the data produced by each evaluation of the field. The pointer passed to
     * evaluate_field() must point to a memory region adhering to the layout described by this map.
     */
    virtual const frantic::channels::channel_map& get_channel_map() const = 0;

    // TODO: We can support these consider extending the interface to support derivative calculations. I imagine we will
    // want to
    //       support derivative calculation on a per-channel basis instead of the entire list of channels.

    /**
     * Declares whether the field can calculate a better derivative than the usual finite differences approach. This
     * means the implementations of evaluate_jacobian(), evaluate_gradient(), evaluate_divergence() & evaluate_curl()
     * are more meaningful than if calculated via finite differences at some small, reasonable offset.
     * @return True if explicit derivatives will be calculated.
     */
    /*virtual bool supports_derivatives() const = 0;

    virtual bool evaluate_gradient( void* dest, const frantic::graphics::vector3f& pos, std::size_t channelNum = 0 )
    const = 0;

    virtual bool evaluate_divergence( void* dest, const frantic::graphics::vector3f& pos, std::size_t channelNum = 0 )
    const = 0;

    virtual bool evaluate_curl( void* dest, const frantic::graphics::vector3f& pos, std::size_t channelNum = 0 ) const =
    0;

    virtual bool evaluate_jacobian( void* dest, const frantic::graphics::vector3f& pos, std::size_t channelNum = 0 )
    const = 0;*/
};

} // namespace volumetrics
} // namespace frantic
