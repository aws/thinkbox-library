// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#ifndef RADIAN_H
#define RADIAN_H

#include <frantic/math/utils.hpp>
#include <iostream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_2PI
#define M_2PI 6.28318530717958647692
#endif

#ifndef DEG_PER_RAD
#define DEG_PER_RAD 57.29577951308232
#endif

#ifndef RAD_PER_DEG
#define RAD_PER_DEG 0.017453292519943
#endif

namespace frantic {
namespace math {
// a class representing a radian angle between
// -PI and PI.
class Radian {
  public:
    Radian( double angle = 0.0 )
        : m_angle( Radian::angleMod( angle ) ) {
        // Nothing
    }

    // from a normal vector

    Radian( double nx, double ny )
        : m_angle( atan2( ny, nx ) ) {
        // Nothing
    }

    ~Radian() {
        // Nothing
    }

    // assignment
    Radian& operator=( Radian const& rhs );
    Radian& operator=( double const rhs );

    bool operator<( Radian const& rhs ) const;
    bool operator<( double const rhs ) const;
    bool operator<=( Radian const& rhs ) const;
    bool operator<=( double const rhs ) const;
    bool operator>( Radian const& rhs ) const;
    bool operator>( double const rhs ) const;
    bool operator>=( Radian const& rhs ) const;
    bool operator>=( double const rhs ) const;
    bool operator==( Radian const& rhs ) const;
    bool operator==( double const rhs ) const;
    bool operator!=( Radian const& rhs ) const;
    bool operator!=( double const rhs ) const;

    // add or subtract an angle.
    Radian& operator+=( Radian const& rhs );
    Radian& operator+( Radian const& rhs ) const;
    Radian& operator-=( Radian const& rhs );
    Radian& operator-( Radian const& rhs ) const;
    // times a scalar.
    Radian& operator*=( double const rhs );
    Radian& operator*( double const rhs ) const;
    Radian& operator/=( double const rhs );
    Radian& operator/( double const rhs ) const;

    operator double() const;
    double magnitude() const;
    double dbl() const;

    // conversion to degrees.
    double degrees() const;
    // to a normal vector
    void normalVector( double& nx, double& ny ) const;
    // to a normal vector.  PointType requires the interface:
    // constructor(double x, double y)
    template <typename VectorType>
    VectorType normalVector() const;

    // static utility functions
    static double angleMod( double angle );
    static double angleDiff( double angle1, double angle2 );
    static double radiansToDegrees( Radian rad );
    static Radian degreesToRadians( double degrees );
    static Radian vectorToAngle( double vx, double vy );

  private:
    double m_angle;
};

bool operator<( double const lhs, Radian const& rhs );
bool operator<=( double const lhs, Radian const& rhs );
bool operator>( double const lhs, Radian const& rhs );
bool operator>=( double const lhs, Radian const& rhs );
bool operator==( double const lhs, Radian const& rhs );
bool operator!=( double const lhs, Radian const& rhs );

template <typename VectorType>
inline VectorType Radian::normalVector() const {
    double nx, ny;
    normalVector( nx, ny );
    return VectorType( nx, ny );
}

std::ostream& operator<<( std::ostream& stream, Radian const& rhs );

inline Radian::operator double() const { return m_angle; }

// add or subtract an angle.
inline Radian& Radian::operator+=( Radian const& rhs ) {
    m_angle = angleMod( m_angle + rhs.m_angle );
    return *this;
}
inline Radian& Radian::operator+( Radian const& rhs ) const {
    Radian angle( m_angle );
    return angle += rhs;
}
inline Radian& Radian::operator-=( Radian const& rhs ) {
    m_angle = angleMod( m_angle - rhs.m_angle );
    return *this;
}
inline Radian& Radian::operator-( Radian const& rhs ) const {
    Radian angle( m_angle );
    return angle -= rhs;
}
// times a scalar.
inline Radian& Radian::operator*=( double const rhs ) {
    m_angle = angleMod( m_angle * rhs );
    return *this;
}
inline Radian& Radian::operator*( double const rhs ) const {
    Radian angle( m_angle );
    return angle *= rhs;
}
inline Radian& Radian::operator/=( double const rhs ) {
    m_angle = angleMod( m_angle / rhs );
    return *this;
}
inline Radian& Radian::operator/( double const rhs ) const {
    Radian angle( m_angle );
    return angle /= rhs;
}

// conversion to degrees.
inline double Radian::degrees() const { return radiansToDegrees( m_angle ); }

// static utility functions
inline double Radian::angleMod( double angle ) { return angle - ( M_2PI * frantic::math::round( angle / ( M_2PI ) ) ); }
inline double Radian::angleDiff( double angle1, double angle2 ) { return angleMod( angle1 - angle2 ); }
inline double Radian::radiansToDegrees( Radian rad ) { return rad.dbl() * DEG_PER_RAD; }
inline Radian Radian::degreesToRadians( double degrees ) { return degrees / DEG_PER_RAD; }

inline double Radian::magnitude() const { return m_angle > 0 ? m_angle : -m_angle; }

inline double Radian::dbl() const { return m_angle; }

inline void Radian::normalVector( double& nx, double& ny ) const {
    nx = cos( m_angle );
    ny = sin( m_angle );
}

inline Radian Radian::vectorToAngle( double vx, double vy ) { return atan2( vy, vx ); }

inline std::ostream& operator<<( std::ostream& stream, Radian const& rhs ) {
    stream << rhs.dbl();
    return stream;
}

inline Radian& Radian::operator=( Radian const& rhs ) {
    m_angle = rhs.m_angle;
    return *this;
}
inline Radian& Radian::operator=( double const rhs ) {
    m_angle = rhs;
    return *this;
}

inline bool Radian::operator<( Radian const& rhs ) const { return m_angle < rhs.m_angle; }
inline bool Radian::operator<( double const rhs ) const { return m_angle < rhs; }
inline bool Radian::operator<=( Radian const& rhs ) const { return m_angle <= rhs.m_angle; }
inline bool Radian::operator<=( double const rhs ) const { return m_angle <= rhs; }
inline bool Radian::operator>( Radian const& rhs ) const { return m_angle > rhs.m_angle; }
inline bool Radian::operator>( double const rhs ) const { return m_angle > rhs; }
inline bool Radian::operator>=( Radian const& rhs ) const { return m_angle >= rhs.m_angle; }
inline bool Radian::operator>=( double const rhs ) const { return m_angle >= rhs; }
inline bool Radian::operator==( Radian const& rhs ) const { return m_angle == rhs.m_angle; }
inline bool Radian::operator==( double const rhs ) const { return m_angle == rhs; }
inline bool Radian::operator!=( Radian const& rhs ) const { return m_angle != rhs.m_angle; }
inline bool Radian::operator!=( double const rhs ) const { return m_angle != rhs; }

inline bool operator<( double const lhs, Radian const& rhs ) { return lhs < rhs.dbl(); }
inline bool operator<=( double const lhs, Radian const& rhs ) { return lhs <= rhs.dbl(); }
inline bool operator>( double const lhs, Radian const& rhs ) { return lhs > rhs.dbl(); }
inline bool operator>=( double const lhs, Radian const& rhs ) { return lhs >= rhs.dbl(); }
inline bool operator==( double const lhs, Radian const& rhs ) { return lhs == rhs.dbl(); }
inline bool operator!=( double const lhs, Radian const& rhs ) { return lhs != rhs.dbl(); }

} // namespace math
} // namespace frantic
#endif
