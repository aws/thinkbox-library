// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <boost/shared_ptr.hpp>
#include <frantic/particles/streams/particle_container_particle_istream.hpp>

namespace frantic {
namespace particles {
namespace streams {

template <class Container>
class shared_particle_container_particle_istream
    : public frantic::particles::streams::particle_container_particle_istream<typename Container::iterator> {
  private:
    boost::shared_ptr<Container> m_particles;

  public:
    shared_particle_container_particle_istream( boost::shared_ptr<Container> particles )
        : particle_container_particle_istream<typename Container::iterator>( particles->begin(), particles->end(),
                                                                             particles->get_channel_map() )
        , m_particles( particles ) {}

    virtual ~shared_particle_container_particle_istream() {}
};

} // namespace streams
} // namespace particles
} // namespace frantic
