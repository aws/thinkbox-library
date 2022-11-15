// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/math/perlin_noise_generator.hpp>

namespace frantic {
namespace math {
// clang-format off
	static int grad2[][2] = {{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,1},{1,-1},{-1,-1}};

	//although we only need 12 gradients for 3d we pad to 16 to avoid division
	//the 4 duplicates form a regular tetrahedron so they produce no visual bias
	//described in Perlin's Improved Noise Paper, link in header
	static int grad3[][3] = {{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},
							{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1},
							{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
							{1,1,0},{-1,1,0},{0,-1,1},{0,-1,-1}};

	static int grad4[][4]= {{0,1,1,1}, {0,1,1,-1}, {0,1,-1,1}, {0,1,-1,-1},
							{0,-1,1,1}, {0,-1,1,-1}, {0,-1,-1,1}, {0,-1,-1,-1},
							{1,0,1,1}, {1,0,1,-1}, {1,0,-1,1}, {1,0,-1,-1},
							{-1,0,1,1}, {-1,0,1,-1}, {-1,0,-1,1}, {-1,0,-1,-1},
							{1,1,0,1}, {1,1,0,-1}, {1,-1,0,1}, {1,-1,0,-1},
							{-1,1,0,1}, {-1,1,0,-1}, {-1,-1,0,1}, {-1,-1,0,-1},
							{1,1,1,0}, {1,1,-1,0}, {1,-1,1,0}, {1,-1,-1,0},
							{-1,1,1,0}, {-1,1,-1,0}, {-1,-1,1,0}, {-1,-1,-1,0}};
// clang-format on

float perlin_noise_generator::get_noise( float x ) { return pnoise( x ); }

float perlin_noise_generator::get_dnoise( float x, float* dx ) { return pnoise( x, dx ); }

float perlin_noise_generator::get_noise( float x, float y ) { return pnoise( x, y ); }

float perlin_noise_generator::get_dnoise( float x, float y, float* dx, float* dy ) { return pnoise( x, y, dx, dy ); }

float perlin_noise_generator::get_noise( float x, float y, float z ) { return pnoise( x, y, z ); }

float perlin_noise_generator::get_dnoise( float x, float y, float z, float* dx, float* dy, float* dz ) {
    return pnoise( x, y, z, dx, dy, dz );
}

float perlin_noise_generator::noise4d( float x, float y, float z, float t ) {
    return dnoise4d( x, y, z, t, 0, 0, 0, 0 );
}

float perlin_noise_generator::dnoise4d( float x, float y, float z, float t, float* dx, float* dy, float* dz,
                                        float* dt ) {
    float total = 0.0f;
    float dx_total = 0.0f;
    float dy_total = 0.0f;
    float dz_total = 0.0f;
    float dt_total = 0.0f;
    float frequency = powf( m_lacunarity, static_cast<float>( m_startingOctave ) );

    for( int i = 0; i < m_octaves; i++ ) {
        float amplitude = powf( m_persistence, static_cast<float>( i ) );

        total += pnoise( x * frequency, y * frequency, z * frequency, t * frequency, dx, dy, dz, dt ) * amplitude;

        if( dx != 0 && dy != 0 && dz != 0 && dt != 0 ) {
            dx_total += *dx;
            dy_total += *dy;
            dz_total += *dz;
            dt_total += *dt;
        }

        frequency *= m_lacunarity;
    }

    if( dx != 0 && dy != 0 && dz != 0 && dt != 0 ) {
        *dx = dx_total;
        *dy = dy_total;
        *dz = dz_total;
        *dt = dt_total;
    }

    return total;
}

void perlin_noise_generator::generate_permutations( unsigned int count, int seed ) {
    srand( seed );

    // round up to the nearest power of 2
    count = round_up_to_power_of_two( count );
    m_permCount = count;

    m_permutations.clear();
    m_permutations.resize( count * 2 );

    // generate numbers
    for( unsigned int i = 0; i < count; i++ )
        m_permutations[i] = i;

    // randomly arrange numbers
    int j, temp;
    for( unsigned int i = count; i > 0; i-- ) {
        temp = m_permutations[i];
        j = rand() & ( count - 1 );
        m_permutations[i] = m_permutations[j];
        m_permutations[j] = temp;
    }

    // extend vector
    for( unsigned int i = 0; i < count; i++ )
        m_permutations[i + count] = m_permutations[i];
}

float perlin_noise_generator::lerp( float m, float a, float b ) { return a + m * ( b - a ); }

float perlin_noise_generator::grad( int hash ) {
    // gradients  -8.0,-7.0 .... 7.0,8.0
    hash &= 15;
    float g = 1.0f + ( hash & 7 );

    return ( hash & 8 ) ? g : -g;
}

float perlin_noise_generator::grad( int hash, float x, float y, float z, float t ) {
    int* grad = grad4[hash & 31];

    return grad[0] * x + grad[1] * y + grad[2] * z + grad[3] * t;
}

float perlin_noise_generator::dot( int* grad, float x, float y ) { return grad[0] * x + grad[1] * y; }

float perlin_noise_generator::dot( int* grad, float x, float y, float z ) {
    return grad[0] * x + grad[1] * y + grad[2] * z;
}

float perlin_noise_generator::dot( int* grad, float x, float y, float z, float t ) {
    return grad[0] * x + grad[1] * y + grad[2] * z + grad[3] * t;
}

float perlin_noise_generator::fade( float t ) { return t * t * t * ( t * ( t * 6 - 15 ) + 10 ); }

// the derivative of the above fade function
float perlin_noise_generator::dfade( float t ) { return t * t * ( t * ( t * 30 - 60 ) + 30 ); }

float perlin_noise_generator::pnoise( float x, float* dx ) {
    int X = FASTFLOOR( x ) & ( m_permCount - 1 );
    x -= floor( x );

    float u = fade( x );

    float g0 = grad( m_permutations[X] );
    float g1 = grad( m_permutations[X + 1] );

    float a = g0 * x;
    float b = g1 * ( x - 1 );

    float noise = lerp( u, a, b );

    // lerp expands to
    // a + u(b-a)
    if( dx != 0 ) {
        float du = dfade( x );

        *dx = g0;
        *dx += u * ( g1 - g0 );
        *dx += du * ( b - a );
        *dx *= 0.25f;
    }

    // can generate a max value of 4
    return noise * 0.25f; // limit to [-1,1]
}

float perlin_noise_generator::pnoise( float x, float y, float* dx, float* dy ) {
    int X = FASTFLOOR( x ) & ( m_permCount - 1 );
    int Y = FASTFLOOR( y ) & ( m_permCount - 1 );

    x -= FASTFLOOR( x );
    y -= FASTFLOOR( y );

    float u = fade( x );
    float v = fade( y );

    int* g0 = grad2[m_permutations[m_permutations[X] + Y] & 7];
    int* g1 = grad2[m_permutations[m_permutations[X + 1] + Y] & 7];
    int* g2 = grad2[m_permutations[m_permutations[X] + Y + 1] & 7];
    int* g3 = grad2[m_permutations[m_permutations[X + 1] + Y + 1] & 7];

    float a = dot( g0, x, y );
    float b = dot( g1, x - 1, y );
    float c = dot( g2, x, y - 1 );
    float d = dot( g3, x - 1, y - 1 );

    float noise = lerp( v, lerp( u, a, b ), lerp( u, c, d ) );

    // expanding the lerp gives us
    // a + u(b-a) + v(c-a) + uv(a-b-c+d)
    if( dx != 0 && dy != 0 ) {
        float du = dfade( x );
        float dv = dfade( y );

        *dx = (float)g0[0];
        *dx += u * ( g1[0] - g0[0] );
        *dx += du * ( b - a );
        *dx += v * ( g2[0] - g0[0] );
        *dx += u * v * ( g0[0] - g1[0] - g2[0] + g3[0] );
        *dx += du * v * ( a - b - c + d );

        *dy = (float)g0[1];
        *dy += u * ( g1[1] - g0[1] );
        *dy += v * ( g2[1] - g0[1] );
        *dy += dv * ( c - a );
        *dy += u * v * ( g0[1] - g1[1] - g2[1] + g3[1] );
        *dy += dv * u * ( a - b - c + d );
    }

    return noise;
}

float perlin_noise_generator::pnoise( float x, float y, float z, float* dx, float* dy, float* dz ) {
    int X = FASTFLOOR( x ) & ( m_permCount - 1 ), Y = FASTFLOOR( y ) & ( m_permCount - 1 ),
        Z = FASTFLOOR( z ) & ( m_permCount - 1 );

    x -= FASTFLOOR( x );
    y -= FASTFLOOR( y );
    z -= FASTFLOOR( z );

    float u = fade( x ), v = fade( y ), w = fade( z );

    int A = m_permutations[X] + Y, AA = m_permutations[A] + Z, AB = m_permutations[A + 1] + Z,
        B = m_permutations[X + 1] + Y, BA = m_permutations[B] + Z, BB = m_permutations[B + 1] + Z;

    int* g0 = grad3[m_permutations[AA] & 15];
    int* g1 = grad3[m_permutations[BA] & 15];
    int* g2 = grad3[m_permutations[AB] & 15];
    int* g3 = grad3[m_permutations[BB] & 15];
    int* g4 = grad3[m_permutations[AA + 1] & 15];
    int* g5 = grad3[m_permutations[BA + 1] & 15];
    int* g6 = grad3[m_permutations[AB + 1] & 15];
    int* g7 = grad3[m_permutations[BB + 1] & 15];

    float a = dot( g0, x, y, z );
    float b = dot( g1, x - 1, y, z );
    float c = dot( g2, x, y - 1, z );
    float d = dot( g3, x - 1, y - 1, z );
    float e = dot( g4, x, y, z - 1 );
    float f = dot( g5, x - 1, y, z - 1 );
    float g = dot( g6, x, y - 1, z - 1 );
    float h = dot( g7, x - 1, y - 1, z - 1 );

    float noise = lerp( w, lerp( v, lerp( u, a, b ), lerp( u, c, d ) ), lerp( v, lerp( u, e, f ), lerp( u, g, h ) ) );

    // the above lerping can be expanded to
    // a + (b-a)u + (c-a)v + (e-a)w + (a-b-c+d)uv + (a-c-e+g)vw + (a-b-e+f)uw + (-a+b+c-d+e-f-g+h)uvw
    if( dx != 0 && dy != 0 && dz != 0 ) {
        float du = dfade( x );
        float dv = dfade( y );
        float dw = dfade( z );

        *dx = 0.0f;
        *dy = 0.0f;

        *dx = (float)g0[0];
        *dx += u * ( g1[0] - g0[0] );
        *dx += du * ( b - a );
        *dx += v * ( g2[0] - g0[0] );
        *dx += w * ( g4[0] - g0[0] );
        *dx += u * v * ( g0[0] - g1[0] - g2[0] + g3[0] );
        *dx += du * v * ( a - b - c + d );
        *dx += v * w * ( g0[0] - g2[0] - g4[0] + g6[0] );
        *dx += u * w * ( g0[0] - g1[0] - g4[0] + g5[0] );
        *dx += du * w * ( a - b - e + f );
        *dx += u * v * w * ( -g0[0] + g1[0] + g2[0] - g3[0] + g4[0] - g5[0] - g6[0] + g7[0] );
        *dx += du * v * w * ( -a + b + c - d + e - f - g + h );

        *dy = (float)g0[1];
        *dy += u * ( g1[1] - g0[1] );
        *dy += v * ( g2[1] - g0[1] );
        *dy += dv * ( c - a );
        *dy += w * ( g4[1] - g0[1] );
        *dy += u * v * ( g0[1] - g1[1] - g2[1] + g3[1] );
        *dy += u * dv * ( a - b - c + d );
        *dy += v * w * ( g0[1] - g2[1] - g4[1] + g6[1] );
        *dy += dv * w * ( a - c - e + g );
        *dy += u * w * ( g0[1] - g1[1] - g4[1] + g5[1] );
        *dy += u * v * w * ( -g0[1] + g1[1] + g2[1] - g3[1] + g4[1] - g5[1] - g6[1] + g7[1] );
        *dy += u * dv * w * ( -a + b + c - d + e - f - g + h );

        *dz = (float)g0[2];
        *dz += u * ( g1[2] - g0[2] );
        *dz += v * ( g2[2] - g0[2] );
        *dz += w * ( g4[2] - g0[2] );
        *dz += dw * ( e - a );
        *dz += u * v * ( g0[2] - g1[2] - g2[2] + g3[2] );
        *dz += v * w * ( g0[2] - g2[2] - g4[2] + g6[2] );
        *dz += dw * v * ( a - c - e + g );
        *dz += u * w * ( g0[2] - g1[2] - g4[2] + g5[2] );
        *dz += dw * u * ( a - b - e + f );
        *dz += u * v * w * ( -g0[2] + g1[2] + g2[2] - g3[2] + g4[2] - g5[2] - g6[2] + g7[2] );
        *dz += u * v * dw * ( -a + b + c - d + e - f - g + h );
    }

    return noise;
}

float perlin_noise_generator::pnoise( float x, float y, float z, float t, float* dx, float* dy, float* dz, float* dt ) {
    int X = FASTFLOOR( x ) & ( m_permCount - 1 );
    int Y = FASTFLOOR( y ) & ( m_permCount - 1 );
    int Z = FASTFLOOR( z ) & ( m_permCount - 1 );
    int T = FASTFLOOR( t ) & ( m_permCount - 1 );

    x -= FASTFLOOR( x );
    y -= FASTFLOOR( y );
    z -= FASTFLOOR( z );
    t -= FASTFLOOR( t );

    float fadeX = fade( x );
    float fadeY = fade( y );
    float fadeZ = fade( z );
    float fadeT = fade( t );

    int* g0 = grad4[( m_permutations[m_permutations[m_permutations[X] + Y] + Z] + T ) & 31];
    int* g1 = grad4[( m_permutations[m_permutations[m_permutations[X + 1] + Y] + Z] + T ) & 31];
    int* g2 = grad4[( m_permutations[m_permutations[m_permutations[X] + Y + 1] + Z] + T ) & 31];
    int* g3 = grad4[( m_permutations[m_permutations[m_permutations[X + 1] + Y + 1] + Z] + T ) & 31];
    int* g4 = grad4[( m_permutations[m_permutations[m_permutations[X] + Y] + Z + 1] + T ) & 31];
    int* g5 = grad4[( m_permutations[m_permutations[m_permutations[X + 1] + Y] + Z + 1] + T ) & 31];
    int* g6 = grad4[( m_permutations[m_permutations[m_permutations[X] + Y + 1] + Z + 1] + T ) & 31];
    int* g7 = grad4[( m_permutations[m_permutations[m_permutations[X + 1] + Y + 1] + Z + 1] + T ) & 31];
    int* g8 = grad4[( m_permutations[m_permutations[m_permutations[X] + Y] + Z] + T + 1 ) & 31];
    int* g9 = grad4[( m_permutations[m_permutations[m_permutations[X + 1] + Y] + Z] + T + 1 ) & 31];
    int* g10 = grad4[( m_permutations[m_permutations[m_permutations[X] + Y + 1] + Z] + T + 1 ) & 31];
    int* g11 = grad4[( m_permutations[m_permutations[m_permutations[X + 1] + Y + 1] + Z] + T + 1 ) & 31];
    int* g12 = grad4[( m_permutations[m_permutations[m_permutations[X] + Y] + Z + 1] + T + 1 ) & 31];
    int* g13 = grad4[( m_permutations[m_permutations[m_permutations[X + 1] + Y] + Z + 1] + T + 1 ) & 31];
    int* g14 = grad4[( m_permutations[m_permutations[m_permutations[X] + Y + 1] + Z + 1] + T + 1 ) & 31];
    int* g15 = grad4[( m_permutations[m_permutations[m_permutations[X + 1] + Y + 1] + Z + 1] + T + 1 ) & 31];

    float a = dot( g0, x, y, z, t );
    float b = dot( g1, x - 1, y, z, t );
    float c = dot( g2, x, y - 1, z, t );
    float d = dot( g3, x - 1, y - 1, z, t );
    float e = dot( g4, x, y, z - 1, t );
    float f = dot( g5, x - 1, y, z - 1, t );
    float g = dot( g6, x, y - 1, z - 1, t );
    float h = dot( g7, x - 1, y - 1, z - 1, t );
    float i = dot( g8, x, y, z, t - 1 );
    float j = dot( g9, x - 1, y, z, t - 1 );
    float k = dot( g10, x, y - 1, z, t - 1 );
    float l = dot( g11, x - 1, y - 1, z, t - 1 );
    float m = dot( g12, x, y, z - 1, t - 1 );
    float n = dot( g13, x - 1, y, z - 1, t - 1 );
    float o = dot( g14, x, y - 1, z - 1, t - 1 );
    float p = dot( g15, x - 1, y - 1, z - 1, t - 1 );

    // lerp t
    float noise = lerp( fadeT,
                        lerp( fadeZ, lerp( fadeY, lerp( fadeX, a, b ), lerp( fadeX, c, d ) ),
                              lerp( fadeY, lerp( fadeX, e, f ), lerp( fadeX, g, h ) ) ),
                        lerp( fadeZ, lerp( fadeY, lerp( fadeX, i, j ), lerp( fadeX, k, l ) ),
                              lerp( fadeY, lerp( fadeX, m, n ), lerp( fadeX, o, p ) ) ) );

    // this can be expanded where u=fadeX, v=fadeY, w=fadeZ, t=fadeT to
    // a + u(b-a) + v(c-a) + w(e-a) + s(i-a) +
    // uv(a-b-c+d) + uw(a-b-e-+f) + ut(a-b-i+j) +
    // vw(a-c-e+g) + vt(a-c-i+k) + wt(a-e-i+m) +
    // uvw(-a+b+c-d+e-f-g+h) + uvt(-a+b+c-d+i-j-k+l) +
    // vwt(-a+c+e-g+i-k-m+o) + uwt(-a+b+e-f+i-j-m+n) +
    // uvwt(a-b-c+d-e+f+g-h-i+j+k-l+m-n-o+p)

    if( dx != 0 && dy != 0 && dz != 0 && dt != 0 ) {
        float dfadeX = dfade( x );
        float dfadeY = dfade( y );
        float dfadeZ = dfade( z );
        float dfadeT = dfade( t );

        *dx = (float)g0[0];
        *dx += fadeX * ( g1[0] - g0[0] ) + dfadeX * ( b - a );
        *dx += fadeY * ( g2[0] - g0[0] );
        *dx += fadeZ * ( g4[0] - g0[0] );
        *dx += fadeT * ( g8[0] - g0[0] );
        *dx += fadeX * fadeY * ( g0[0] - g1[0] - g2[0] + g3[0] ) + dfadeX * fadeY * ( a - b - c + d );
        *dx += fadeX * fadeZ * ( g0[0] - g1[0] - g4[0] + g5[0] ) + dfadeX * fadeZ * ( a - b - e + f );
        *dx += fadeX * fadeT * ( g0[0] - g1[0] - g8[0] + g9[0] ) + dfadeX * fadeT * ( a - b - i + j );
        *dx += fadeY * fadeZ * ( g0[0] - g2[0] - g4[0] + g6[0] );
        *dx += fadeY * fadeT * ( g0[0] - g2[0] - g8[0] + g9[0] );
        *dx += fadeZ * fadeT * ( g0[0] - g4[0] - g8[0] + g12[0] );
        *dx += fadeX * fadeY * fadeZ * ( -g0[0] + g1[0] + g2[0] - g3[0] + g4[0] - g5[0] - g6[0] + g7[0] ) +
               dfadeX * fadeY * fadeZ * ( -a + b + c - d + e - f - g + h );
        *dx += fadeX * fadeY * fadeT * ( -g0[0] + g1[0] + g2[0] - g3[0] + g8[0] - g9[0] - g10[0] + g11[0] ) +
               dfadeX * fadeY * fadeT * ( -a + b + c - d + i - j - k + l );
        *dx += fadeY * fadeZ * fadeT * ( -g0[0] + g2[0] + g4[0] - g6[0] + g8[0] - g10[0] - g12[0] + g14[0] );
        *dx += fadeX * fadeZ * fadeT * ( -g0[0] + g1[0] + g4[0] - g5[0] + g8[0] - g9[0] - g12[0] + g13[0] ) +
               dfadeX * fadeZ * fadeT * ( -a + b + e - f + i - j - m + n );
        *dx += fadeX * fadeY * fadeZ * fadeT *
               ( g0[0] - g1[0] - g2[0] + g3[0] - g4[0] + g5[0] + g6[0] - g7[0] - g8[0] + g9[0] + g10[0] - g11[0] +
                 g12[0] - g13[0] - g14[0] + g15[0] );
        *dx += dfadeX * fadeY * fadeZ * fadeT * ( a - b - c + d - e + f + g - h - i + j + k - l + m - n - o + p );

        *dy = (float)g0[1];
        *dy += fadeX * ( g1[1] - g0[1] );
        *dy += fadeY * ( g2[1] - g0[1] ) + dfadeY * ( c - a );
        *dy += fadeZ * ( g4[1] - g0[1] );
        *dy += fadeT * ( g8[1] - g0[1] );
        *dy += fadeX * fadeY * ( g0[1] - g1[1] - g2[1] + g3[1] ) + fadeX * dfadeY * ( a - b - c + d );
        *dy += fadeX * fadeZ * ( g0[1] - g1[1] - g4[1] + g5[1] );
        *dy += fadeX * fadeT * ( g0[1] - g1[1] - g8[1] + g9[1] );
        *dy += fadeY * fadeZ * ( g0[1] - g2[1] - g4[1] + g6[1] ) + dfadeY * fadeZ * ( a - c - e + g );
        *dy += fadeY * fadeT * ( g0[1] - g2[1] - g8[1] + g9[1] ) + dfadeY * fadeT * ( a - b - i + j );
        *dy += fadeZ * fadeT * ( g0[1] - g4[1] - g8[1] + g12[1] );
        *dy += fadeX * fadeY * fadeZ * ( -g0[1] + g1[1] + g2[1] - g3[1] + g4[1] - g5[1] - g6[1] + g7[1] ) +
               fadeX * dfadeY * fadeZ * ( -a + b + c - d + e - f - g + h );
        *dy += fadeX * fadeY * fadeT * ( -g0[1] + g1[1] + g2[1] - g3[1] + g8[1] - g9[1] - g10[1] + g11[1] ) +
               fadeX * dfadeY * fadeT * ( -a + b + c - d + i - j - k + l );
        *dy += fadeY * fadeZ * fadeT * ( -g0[1] + g2[1] + g4[1] - g6[1] + g8[1] - g10[1] - g12[1] + g14[1] ) +
               dfadeY * fadeZ * fadeT * ( -a + c + e - g + i - k - m + o );
        *dy += fadeX * fadeZ * fadeT * ( -g0[1] + g1[1] + g4[1] - g5[1] + g8[1] - g9[1] - g12[1] + g13[1] );
        *dy += fadeX * fadeY * fadeZ * fadeT *
               ( g0[1] - g1[1] - g2[1] + g3[1] - g4[1] + g5[1] + g6[1] - g7[1] - g8[1] + g9[1] + g10[1] - g11[1] +
                 g12[1] - g13[1] - g14[1] + g15[1] );
        *dy += fadeX * dfadeY * fadeZ * fadeT * ( a - b - c + d - e + f + g - h - i + j + k - l + m - n - o + p );

        *dz = (float)g0[2];
        *dz += fadeX * ( g1[2] - g0[2] );
        *dz += fadeY * ( g2[2] - g0[2] );
        *dz += fadeZ * ( g4[2] - g0[2] ) + dfadeZ * ( e - a );
        *dz += fadeT * ( g8[2] - g0[2] );
        *dz += fadeX * fadeY * ( g0[2] - g1[2] - g2[2] + g3[2] );
        *dz += fadeX * fadeZ * ( g0[2] - g1[2] - g4[2] + g5[2] ) + fadeX * dfadeZ * ( a - b - e + f );
        *dz += fadeX * fadeT * ( g0[2] - g1[2] - g8[2] + g9[2] );
        *dz += fadeY * fadeZ * ( g0[2] - g2[2] - g4[2] + g6[2] ) + fadeY * dfadeZ * ( a - c - e + g );
        *dz += fadeY * fadeT * ( g0[2] - g2[2] - g8[2] + g9[2] );
        *dz += fadeZ * fadeT * ( g0[2] - g4[2] - g8[2] + g12[2] ) + dfadeZ * fadeT * ( a - e - i + m );
        *dz += fadeX * fadeY * fadeZ * ( -g0[2] + g1[2] + g2[2] - g3[2] + g4[2] - g5[2] - g6[2] + g7[2] ) +
               fadeX * fadeY * dfadeZ * ( -a + b + c - d + e - f - g + h );
        *dz += fadeX * fadeY * fadeT * ( -g0[2] + g1[2] + g2[2] - g3[2] + g8[2] - g9[2] - g10[2] + g11[2] );
        *dz += fadeY * fadeZ * fadeT * ( -g0[2] + g2[2] + g4[2] - g6[2] + g8[2] - g10[2] - g12[2] + g14[2] ) +
               fadeY * dfadeZ * fadeT * ( -a + c + e - g + i - k - m + o );
        *dz += fadeX * fadeZ * fadeT * ( -g0[2] + g1[2] + g4[2] - g5[2] + g8[2] - g9[2] - g12[2] + g13[2] ) +
               fadeX * dfadeZ * fadeT * ( -a + b + e - f + i - j - m + n );
        *dz += fadeX * fadeY * fadeZ * fadeT *
               ( g0[2] - g1[2] - g2[2] + g3[2] - g4[2] + g5[2] + g6[2] - g7[2] - g8[2] + g9[2] + g10[2] - g11[2] +
                 g12[2] - g13[2] - g14[2] + g15[2] );
        *dz += fadeX * fadeY * dfadeZ * fadeT * ( a - b - c + d - e + f + g - h - i + j + k - l + m - n - o + p );

        *dt = (float)g0[3];
        *dt += fadeX * ( g1[3] - g0[3] );
        *dt += fadeY * ( g2[3] - g0[3] );
        *dt += fadeZ * ( g4[3] - g0[3] );
        *dt += fadeT * ( g8[3] - g0[3] ) + dfadeT * ( i - a );
        *dt += fadeX * fadeY * ( g0[3] - g1[3] - g2[3] + g3[3] );
        *dt += fadeX * fadeZ * ( g0[3] - g1[3] - g4[3] + g5[3] );
        *dt += fadeX * fadeT * ( g0[3] - g1[3] - g8[3] + g9[3] ) + fadeX * dfadeT * ( a - b - i + k );
        *dt += fadeY * fadeZ * ( g0[3] - g2[3] - g4[3] + g6[3] );
        *dt += fadeY * fadeT * ( g0[3] - g2[3] - g8[3] + g9[3] ) + fadeY * dfadeT * ( a - b - i + j );
        *dt += fadeZ * fadeT * ( g0[3] - g4[3] - g8[3] + g12[3] ) + fadeZ * dfadeT * ( a - e - i + m );
        *dt += fadeX * fadeY * fadeZ * ( -g0[3] + g1[3] + g2[3] - g3[3] + g4[3] - g5[3] - g6[3] + g7[3] );
        *dt += fadeX * fadeY * fadeT * ( -g0[3] + g1[3] + g2[3] - g3[3] + g8[3] - g9[3] - g10[3] + g11[3] ) +
               fadeX * fadeY * dfadeT * ( -a + b + c - d + i - j - k + l );
        *dt += fadeY * fadeZ * fadeT * ( -g0[3] + g2[3] + g4[3] - g6[3] + g8[3] - g10[3] - g12[3] + g14[3] ) +
               fadeY * fadeZ * dfadeT * ( -a + c + e - g + i - k - m + o );
        *dt += fadeX * fadeZ * fadeT * ( -g0[3] + g1[3] + g4[3] - g5[3] + g8[3] - g9[3] - g12[3] + g13[3] ) +
               fadeX * fadeZ * dfadeT * ( -a + b + e - f + i - j - m + n );
        *dt += fadeX * fadeY * fadeZ * fadeT *
               ( g0[3] - g1[3] - g2[3] + g3[3] - g4[3] + g5[3] + g6[3] - g7[3] - g8[3] + g9[3] + g10[3] - g11[3] +
                 g12[3] - g13[3] - g14[3] + g15[3] );
        *dt += fadeX * fadeY * fadeZ * dfadeT * ( a - b - c + d - e + f + g - h - i + j + k - l + m - n - o + p );
    }

    return noise;
}

} // namespace math
} // namespace frantic
