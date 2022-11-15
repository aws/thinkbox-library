// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/diagnostics/profiling_section.hpp>

namespace frantic {
namespace diagnostics {

class profiling_manager {
  private:
    std::vector<profiling_section> m_profs;
    template <class CharType>
    friend std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>&, const profiling_manager& );

  public:
    profiling_manager() {}
    int new_profiling_section( const frantic::tstring& name ) {
        int result = (int)m_profs.size();
        m_profs.push_back( profiling_section( name ) );
        return result;
    }

    void enter_section( int i ) { m_profs[i].enter(); }

    void exit_section( int i ) { m_profs[i].exit(); }
};

template <class CharType>
inline std::basic_ostream<CharType>& operator<<( std::basic_ostream<CharType>& o, const profiling_manager& pfm ) {
    for( int i = 0; i < (int)pfm.m_profs.size(); ++i )
        o << pfm.m_profs[i] << "\n";
    return o;
}

} // namespace diagnostics
} // namespace frantic
