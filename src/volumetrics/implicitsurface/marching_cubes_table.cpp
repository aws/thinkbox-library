// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @author Mike Yurick
 *
 * Adapted from the paper/implementation:
 * "Efficient Implementation of Marching Cubes' cases with topological guarantees".
 *	Lewinar, et al.
 *
 * He uses a different sign conventions, cube corner and edge setup, so the implementation has
 * been reworked  to be made usable with our cube conventions.  There were also a couple of errors
 * in his tables/tests which I have corrected.  There is a unit test in the VolumetricsTestSuite
 * called testMarchingCubesTopology which verifies the correctness of the topology of a mesh
 * created from a randomly generated level set.
 */

// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/volumetrics/implicitsurface/marching_cubes_table.hpp>

#pragma warning( push )
#pragma warning( disable : 4512 )
#include <boost/assign/std/vector.hpp>
#pragma warning( pop )

namespace frantic {
namespace volumetrics {
using frantic::graphics::vector3;

marching_cubes_table::marching_cubes_table() {

    using namespace boost::assign;

    /**
     * Classic Marching Cubes table.
     *
     * For each of the possible cube cases that a cube can take (determined by the cube
     * corner values as indicated above), a correct triangulation of the edge vertices
     * is listed in the appropriate entry of this table.  The table entries are in the form
     * of triplets of edge vertices which describe faces, where the edge vertices are
     * numbered as indicated above.
     */

    m_classicCubeCases = std::vector<std::vector<vector3>>( 256 );

    /*   0:                          */
    /*   1: 0,                       */ m_classicCubeCases[1] += vector3( 2, 0, 1 );
    /*   2:    1,                    */ m_classicCubeCases[2] += vector3( 3, 0, 4 );
    /*   3: 0, 1,                    */ m_classicCubeCases[3] += vector3( 2, 3, 1 ), vector3( 2, 4, 3 );
    /*   4:       2,                 */ m_classicCubeCases[4] += vector3( 6, 1, 5 );
    /*   5: 0,    2,                 */ m_classicCubeCases[5] += vector3( 6, 0, 5 ), vector3( 6, 2, 0 );
    /*   6:    1, 2,                 */ m_classicCubeCases[6] += vector3( 4, 3, 0 ), vector3( 1, 5, 6 );
    /*   7: 0, 1, 2,                 */ m_classicCubeCases[7] += vector3( 6, 3, 5 ), vector3( 4, 3, 6 ),
        vector3( 2, 4, 6 );
    /*   8:          3,              */ m_classicCubeCases[8] += vector3( 5, 3, 7 );
    /*   9: 0,       3,              */ m_classicCubeCases[9] += vector3( 2, 0, 1 ), vector3( 5, 3, 7 );
    /*  10:    1,    3,              */ m_classicCubeCases[10] += vector3( 5, 4, 7 ), vector3( 5, 0, 4 );
    /*  11: 0, 1,    3,              */ m_classicCubeCases[11] += vector3( 2, 5, 1 ), vector3( 7, 5, 2 ),
        vector3( 4, 7, 2 );
    /*  12:       2, 3,              */ m_classicCubeCases[12] += vector3( 7, 1, 3 ), vector3( 7, 6, 1 );
    /*  13: 0,    2, 3,              */ m_classicCubeCases[13] += vector3( 7, 0, 3 ), vector3( 2, 0, 7 ),
        vector3( 6, 2, 7 );
    /*  14:    1, 2, 3,              */ m_classicCubeCases[14] += vector3( 4, 1, 0 ), vector3( 6, 1, 4 ),
        vector3( 7, 6, 4 );
    /*  15: 0, 1, 2, 3,              */ m_classicCubeCases[15] += vector3( 2, 4, 7 ), vector3( 2, 7, 6 );
    /*  16:             4,           */ m_classicCubeCases[16] += vector3( 9, 8, 2 );
    /*  17: 0,          4,           */ m_classicCubeCases[17] += vector3( 1, 8, 0 ), vector3( 1, 9, 8 );
    /*  18:    1,       4,           */ m_classicCubeCases[18] += vector3( 3, 0, 4 ), vector3( 8, 2, 9 );
    /*  19: 0, 1,       4,           */ m_classicCubeCases[19] += vector3( 3, 8, 4 ), vector3( 9, 8, 3 ),
        vector3( 1, 9, 3 );
    /*  20:       2,    4,           */ m_classicCubeCases[20] += vector3( 8, 2, 9 ), vector3( 6, 1, 5 );
    /*  21: 0,    2,    4,           */ m_classicCubeCases[21] += vector3( 8, 6, 9 ), vector3( 5, 6, 8 ),
        vector3( 0, 5, 8 );
    /*  22:    1, 2,    4,           */ m_classicCubeCases[22] += vector3( 0, 4, 3 ), vector3( 8, 2, 9 ),
        vector3( 1, 5, 6 );
    /*  23: 0, 1, 2,    4,           */ m_classicCubeCases[23] += vector3( 9, 8, 6 ), vector3( 8, 4, 6 ),
        vector3( 6, 4, 5 ), vector3( 5, 4, 3 );
    /*  24:          3, 4,           */ m_classicCubeCases[24] += vector3( 5, 3, 7 ), vector3( 8, 2, 9 );
    /*  25: 0,       3, 4,           */ m_classicCubeCases[25] += vector3( 8, 1, 9 ), vector3( 0, 1, 8 ),
        vector3( 5, 3, 7 );
    /*  26:    1,    3, 4,           */ m_classicCubeCases[26] += vector3( 5, 4, 7 ), vector3( 0, 4, 5 ),
        vector3( 8, 2, 9 );
    /*  27: 0, 1,    3, 4,           */ m_classicCubeCases[27] += vector3( 7, 5, 4 ), vector3( 4, 5, 9 ),
        vector3( 9, 5, 1 ), vector3( 4, 9, 8 );
    /*  28:       2, 3, 4,           */ m_classicCubeCases[28] += vector3( 7, 1, 3 ), vector3( 6, 1, 7 ),
        vector3( 2, 9, 8 );
    /*  29: 0,    2, 3, 4,           */ m_classicCubeCases[29] += vector3( 6, 3, 7 ), vector3( 8, 3, 6 ),
        vector3( 0, 3, 8 ), vector3( 6, 9, 8 );
    /*  30:    1, 2, 3, 4,           */ m_classicCubeCases[30] += vector3( 9, 8, 2 ), vector3( 0, 4, 6 ),
        vector3( 6, 4, 7 ), vector3( 0, 6, 1 );
    /*  31: 0, 1, 2, 3, 4,           */ m_classicCubeCases[31] += vector3( 9, 8, 6 ), vector3( 6, 8, 4 ),
        vector3( 6, 4, 7 );
    /*  32:                5,        */ m_classicCubeCases[32] += vector3( 10, 4, 8 );
    /*  33: 0,             5,        */ m_classicCubeCases[33] += vector3( 10, 4, 8 ), vector3( 2, 0, 1 );
    /*  34:    1,          5,        */ m_classicCubeCases[34] += vector3( 10, 0, 8 ), vector3( 10, 3, 0 );
    /*  35: 0, 1,          5,        */ m_classicCubeCases[35] += vector3( 10, 2, 8 ), vector3( 1, 2, 10 ),
        vector3( 3, 1, 10 );
    /*  36:       2,       5,        */ m_classicCubeCases[36] += vector3( 10, 4, 8 ), vector3( 1, 5, 6 );
    /*  37: 0,    2,       5,        */ m_classicCubeCases[37] += vector3( 6, 0, 5 ), vector3( 2, 0, 6 ),
        vector3( 4, 8, 10 );
    /*  38:    1, 2,       5,        */ m_classicCubeCases[38] += vector3( 10, 0, 8 ), vector3( 3, 0, 10 ),
        vector3( 1, 5, 6 );
    /*  39: 0, 1, 2,       5,        */ m_classicCubeCases[39] += vector3( 3, 5, 10 ), vector3( 10, 5, 2 ),
        vector3( 2, 5, 6 ), vector3( 2, 8, 10 );
    /*  40:          3,    5,        */ m_classicCubeCases[40] += vector3( 5, 3, 7 ), vector3( 10, 4, 8 );
    /*  41: 0,       3,    5,        */ m_classicCubeCases[41] += vector3( 0, 1, 2 ), vector3( 5, 3, 7 ),
        vector3( 4, 8, 10 );
    /*  42:    1,    3,    5,        */ m_classicCubeCases[42] += vector3( 5, 10, 7 ), vector3( 8, 10, 5 ),
        vector3( 0, 8, 5 );
    /*  43: 0, 1,    3,    5,        */ m_classicCubeCases[43] += vector3( 7, 5, 10 ), vector3( 5, 1, 10 ),
        vector3( 10, 1, 8 ), vector3( 8, 1, 2 );
    /*  44:       2, 3,    5,        */ m_classicCubeCases[44] += vector3( 1, 7, 6 ), vector3( 3, 7, 1 ),
        vector3( 10, 4, 8 );
    /*  45: 0,    2, 3,    5,        */ m_classicCubeCases[45] += vector3( 4, 8, 10 ), vector3( 2, 0, 3 ),
        vector3( 7, 2, 3 ), vector3( 6, 2, 7 );
    /*  46:    1, 2, 3,    5,        */ m_classicCubeCases[46] += vector3( 8, 10, 0 ), vector3( 0, 10, 6 ),
        vector3( 6, 10, 7 ), vector3( 0, 6, 1 );
    /*  47: 0, 1, 2, 3,    5,        */ m_classicCubeCases[47] += vector3( 8, 10, 2 ), vector3( 2, 10, 7 ),
        vector3( 2, 7, 6 );
    /*  48:             4, 5,        */ m_classicCubeCases[48] += vector3( 9, 4, 2 ), vector3( 9, 10, 4 );
    /*  49: 0,          4, 5,        */ m_classicCubeCases[49] += vector3( 1, 4, 0 ), vector3( 10, 4, 1 ),
        vector3( 9, 10, 1 );
    /*  50:    1,       4, 5,        */ m_classicCubeCases[50] += vector3( 9, 0, 2 ), vector3( 3, 0, 9 ),
        vector3( 10, 3, 9 );
    /*  51: 0, 1,       4, 5,        */ m_classicCubeCases[51] += vector3( 10, 3, 1 ), vector3( 10, 1, 9 );
    /*  52:       2,    4, 5,        */ m_classicCubeCases[52] += vector3( 4, 9, 10 ), vector3( 2, 9, 4 ),
        vector3( 6, 1, 5 );
    /*  53: 0,    2,    4, 5,        */ m_classicCubeCases[53] += vector3( 10, 4, 9 ), vector3( 9, 4, 5 ),
        vector3( 5, 4, 0 ), vector3( 9, 5, 6 );
    /*  54:    1, 2,    4, 5,        */ m_classicCubeCases[54] += vector3( 1, 5, 6 ), vector3( 3, 0, 2 ),
        vector3( 9, 3, 2 ), vector3( 10, 3, 9 );
    /*  55: 0, 1, 2,    4, 5,        */ m_classicCubeCases[55] += vector3( 5, 6, 3 ), vector3( 3, 6, 9 ),
        vector3( 3, 9, 10 );
    /*  56:          3, 4, 5,        */ m_classicCubeCases[56] += vector3( 9, 4, 2 ), vector3( 10, 4, 9 ),
        vector3( 3, 7, 5 );
    /*  57: 0,       3, 4, 5,        */ m_classicCubeCases[57] += vector3( 3, 7, 5 ), vector3( 10, 4, 0 ),
        vector3( 1, 10, 0 ), vector3( 9, 10, 1 );
    /*  58:    1,    3, 4, 5,        */ m_classicCubeCases[58] += vector3( 0, 2, 5 ), vector3( 5, 2, 10 ),
        vector3( 10, 2, 9 ), vector3( 10, 7, 5 );
    /*  59: 0, 1,    3, 4, 5,        */ m_classicCubeCases[59] += vector3( 7, 5, 10 ), vector3( 10, 5, 1 ),
        vector3( 10, 1, 9 );
    /*  60:       2, 3, 4, 5,        */ m_classicCubeCases[60] += vector3( 10, 4, 2 ), vector3( 10, 2, 9 ),
        vector3( 3, 7, 1 ), vector3( 1, 7, 6 );
    /*  61: 0,    2, 3, 4, 5,        */ m_classicCubeCases[61] += vector3( 9, 10, 0 ), vector3( 0, 10, 4 ),
        vector3( 6, 9, 0 ), vector3( 0, 3, 7 ), vector3( 7, 6, 0 );
    /*  62:    1, 2, 3, 4, 5,        */ m_classicCubeCases[62] += vector3( 7, 6, 0 ), vector3( 0, 6, 1 ),
        vector3( 10, 7, 0 ), vector3( 0, 2, 9 ), vector3( 9, 10, 0 );
    /*  63: 0, 1, 2, 3, 4, 5,        */ m_classicCubeCases[63] += vector3( 7, 6, 10 ), vector3( 6, 9, 10 );
    /*  64:                   6,     */ m_classicCubeCases[64] += vector3( 11, 9, 6 );
    /*  65: 0,                6,     */ m_classicCubeCases[65] += vector3( 0, 1, 2 ), vector3( 9, 6, 11 );
    /*  66:    1,             6,     */ m_classicCubeCases[66] += vector3( 3, 0, 4 ), vector3( 9, 6, 11 );
    /*  67: 0, 1,             6,     */ m_classicCubeCases[67] += vector3( 3, 2, 4 ), vector3( 1, 2, 3 ),
        vector3( 9, 6, 11 );
    /*  68:       2,          6,     */ m_classicCubeCases[68] += vector3( 5, 9, 1 ), vector3( 5, 11, 9 );
    /*  69: 0,    2,          6,     */ m_classicCubeCases[69] += vector3( 0, 9, 2 ), vector3( 11, 9, 0 ),
        vector3( 5, 11, 0 );
    /*  70:    1, 2,          6,     */ m_classicCubeCases[70] += vector3( 9, 5, 11 ), vector3( 1, 5, 9 ),
        vector3( 3, 0, 4 );
    /*  71: 0, 1, 2,          6,     */ m_classicCubeCases[71] += vector3( 11, 3, 5 ), vector3( 2, 3, 11 ),
        vector3( 4, 3, 2 ), vector3( 9, 2, 11 );
    /*  72:          3,       6,     */ m_classicCubeCases[72] += vector3( 3, 7, 5 ), vector3( 6, 11, 9 );
    /*  73: 0,       3,       6,     */ m_classicCubeCases[73] += vector3( 5, 3, 7 ), vector3( 0, 1, 2 ),
        vector3( 6, 11, 9 );
    /*  74:    1,    3,       6,     */ m_classicCubeCases[74] += vector3( 4, 5, 0 ), vector3( 7, 5, 4 ),
        vector3( 6, 11, 9 );
    /*  75: 0, 1,    3,       6,     */ m_classicCubeCases[75] += vector3( 6, 11, 9 ), vector3( 7, 5, 1 ),
        vector3( 2, 7, 1 ), vector3( 4, 7, 2 );
    /*  76:       2, 3,       6,     */ m_classicCubeCases[76] += vector3( 9, 7, 11 ), vector3( 3, 7, 9 ),
        vector3( 1, 3, 9 );
    /*  77: 0,    2, 3,       6,     */ m_classicCubeCases[77] += vector3( 9, 7, 11 ), vector3( 9, 3, 7 ),
        vector3( 2, 3, 9 ), vector3( 0, 3, 2 );
    /*  78:    1, 2, 3,       6,     */ m_classicCubeCases[78] += vector3( 1, 0, 9 ), vector3( 9, 0, 7 ),
        vector3( 7, 0, 4 ), vector3( 7, 11, 9 );
    /*  79: 0, 1, 2, 3,       6,     */ m_classicCubeCases[79] += vector3( 11, 9, 7 ), vector3( 7, 9, 2 ),
        vector3( 7, 2, 4 );
    /*  80:             4,    6,     */ m_classicCubeCases[80] += vector3( 2, 11, 8 ), vector3( 2, 6, 11 );
    /*  81: 0,          4,    6,     */ m_classicCubeCases[81] += vector3( 11, 1, 6 ), vector3( 0, 1, 11 ),
        vector3( 8, 0, 11 );
    /*  82:    1,       4,    6,     */ m_classicCubeCases[82] += vector3( 11, 2, 6 ), vector3( 8, 2, 11 ),
        vector3( 0, 4, 3 );
    /*  83: 0, 1,       4,    6,     */ m_classicCubeCases[83] += vector3( 8, 4, 11 ), vector3( 11, 4, 1 ),
        vector3( 1, 4, 3 ), vector3( 1, 6, 11 );
    /*  84:       2,    4,    6,     */ m_classicCubeCases[84] += vector3( 5, 2, 1 ), vector3( 8, 2, 5 ),
        vector3( 11, 8, 5 );
    /*  85: 0,    2,    4,    6,     */ m_classicCubeCases[85] += vector3( 8, 0, 5 ), vector3( 11, 8, 5 );
    /*  86:    1, 2,    4,    6,     */ m_classicCubeCases[86] += vector3( 4, 3, 0 ), vector3( 1, 5, 8 ),
        vector3( 8, 5, 11 ), vector3( 1, 8, 2 );
    /*  87: 0, 1, 2,    4,    6,     */ m_classicCubeCases[87] += vector3( 4, 3, 8 ), vector3( 8, 3, 5 ),
        vector3( 8, 5, 11 );
    /*  88:          3, 4,    6,     */ m_classicCubeCases[88] += vector3( 2, 11, 8 ), vector3( 6, 11, 2 ),
        vector3( 7, 5, 3 );
    /*  89: 0,       3, 4,    6,     */ m_classicCubeCases[89] += vector3( 5, 3, 7 ), vector3( 0, 1, 6 ),
        vector3( 11, 0, 6 ), vector3( 8, 0, 11 );
    /*  90:    1,    3, 4,    6,     */ m_classicCubeCases[90] += vector3( 6, 8, 2 ), vector3( 11, 8, 6 ),
        vector3( 5, 0, 4 ), vector3( 7, 5, 4 );
    /*  91: 0, 1,    3, 4,    6,     */ m_classicCubeCases[91] += vector3( 4, 7, 1 ), vector3( 1, 7, 5 ),
        vector3( 8, 4, 1 ), vector3( 1, 6, 11 ), vector3( 11, 8, 1 );
    /*  92:       2, 3, 4,    6,     */ m_classicCubeCases[92] += vector3( 3, 2, 1 ), vector3( 11, 2, 3 ),
        vector3( 8, 2, 11 ), vector3( 7, 11, 3 );
    /*  93: 0,    2, 3, 4,    6,     */ m_classicCubeCases[93] += vector3( 3, 7, 0 ), vector3( 0, 7, 11 ),
        vector3( 0, 11, 8 );
    /*  94:    1, 2, 3, 4,    6,     */ m_classicCubeCases[94] += vector3( 11, 8, 1 ), vector3( 1, 8, 2 ),
        vector3( 7, 11, 1 ), vector3( 1, 0, 4 ), vector3( 4, 7, 1 );
    /*  95: 0, 1, 2, 3, 4,    6,     */ m_classicCubeCases[95] += vector3( 4, 7, 8 ), vector3( 7, 11, 8 );
    /*  96:                5, 6,     */ m_classicCubeCases[96] += vector3( 4, 8, 10 ), vector3( 11, 9, 6 );
    /*  97: 0,             5, 6,     */ m_classicCubeCases[97] += vector3( 2, 0, 1 ), vector3( 4, 8, 10 ),
        vector3( 9, 6, 11 );
    /*  98:    1,          5, 6,     */ m_classicCubeCases[98] += vector3( 0, 10, 3 ), vector3( 8, 10, 0 ),
        vector3( 11, 9, 6 );
    /*  99: 0, 1,          5, 6,     */ m_classicCubeCases[99] += vector3( 9, 6, 11 ), vector3( 1, 2, 8 ),
        vector3( 10, 1, 8 ), vector3( 3, 1, 10 );
    /* 100:       2,       5, 6,     */ m_classicCubeCases[100] += vector3( 5, 9, 1 ), vector3( 11, 9, 5 ),
        vector3( 8, 10, 4 );
    /* 101: 0,    2,       5, 6,     */ m_classicCubeCases[101] += vector3( 10, 4, 8 ), vector3( 2, 0, 11 ),
        vector3( 11, 0, 5 ), vector3( 2, 11, 9 );
    /* 102:    1, 2,       5, 6,     */ m_classicCubeCases[102] += vector3( 11, 1, 5 ), vector3( 9, 1, 11 ),
        vector3( 10, 3, 0 ), vector3( 8, 10, 0 );
    /* 103: 0, 1, 2,       5, 6,     */ m_classicCubeCases[103] += vector3( 5, 11, 2 ), vector3( 2, 11, 9 ),
        vector3( 3, 5, 2 ), vector3( 2, 8, 10 ), vector3( 10, 3, 2 );
    /* 104:          3,    5, 6,     */ m_classicCubeCases[104] += vector3( 10, 4, 8 ), vector3( 3, 7, 5 ),
        vector3( 11, 9, 6 );
    /* 105: 0,       3,    5, 6,     */ m_classicCubeCases[105] += vector3( 6, 11, 9 ), vector3( 5, 3, 7 ),
        vector3( 2, 0, 1 ), vector3( 4, 8, 10 );
    /* 106:    1,    3,    5, 6,     */ m_classicCubeCases[106] += vector3( 11, 9, 6 ), vector3( 8, 10, 7 ),
        vector3( 5, 8, 7 ), vector3( 0, 8, 5 );
    /* 107: 0, 1,    3,    5, 6,     */ m_classicCubeCases[107] += vector3( 8, 1, 2 ), vector3( 10, 1, 8 ),
        vector3( 5, 1, 10 ), vector3( 10, 7, 5 ), vector3( 9, 6, 11 );
    /* 108:       2, 3,    5, 6,     */ m_classicCubeCases[108] += vector3( 10, 4, 8 ), vector3( 3, 7, 11 ),
        vector3( 9, 3, 11 ), vector3( 1, 3, 9 );
    /* 109: 0,    2, 3,    5, 6,     */ m_classicCubeCases[109] += vector3( 11, 3, 7 ), vector3( 9, 3, 11 ),
        vector3( 0, 3, 9 ), vector3( 9, 2, 0 ), vector3( 10, 4, 8 );
    /* 110:    1, 2, 3,    5, 6,     */ m_classicCubeCases[110] += vector3( 0, 8, 7 ), vector3( 7, 8, 10 ),
        vector3( 1, 0, 7 ), vector3( 7, 11, 9 ), vector3( 9, 1, 7 );
    /* 111: 0, 1, 2, 3,    5, 6,     */ m_classicCubeCases[111] += vector3( 11, 9, 7 ), vector3( 7, 9, 2 ),
        vector3( 8, 10, 7 ), vector3( 2, 8, 7 );
    /* 112:             4, 5, 6,     */ m_classicCubeCases[112] += vector3( 4, 11, 10 ), vector3( 6, 11, 4 ),
        vector3( 2, 6, 4 );
    /* 113: 0,          4, 5, 6,     */ m_classicCubeCases[113] += vector3( 11, 1, 6 ), vector3( 11, 0, 1 ),
        vector3( 10, 0, 11 ), vector3( 4, 0, 10 );
    /* 114:    1,       4, 5, 6,     */ m_classicCubeCases[114] += vector3( 6, 0, 2 ), vector3( 10, 0, 6 ),
        vector3( 3, 0, 10 ), vector3( 11, 10, 6 );
    /* 115: 0, 1,       4, 5, 6,     */ m_classicCubeCases[115] += vector3( 6, 11, 1 ), vector3( 1, 11, 10 ),
        vector3( 1, 10, 3 );
    /* 116:       2,    4, 5, 6,     */ m_classicCubeCases[116] += vector3( 2, 10, 4 ), vector3( 5, 10, 2 ),
        vector3( 11, 10, 5 ), vector3( 2, 1, 5 );
    /* 117: 0,    2,    4, 5, 6,     */ m_classicCubeCases[117] += vector3( 10, 4, 11 ), vector3( 11, 4, 0 ),
        vector3( 11, 0, 5 );
    /* 118:    1, 2,    4, 5, 6,     */ m_classicCubeCases[118] += vector3( 10, 3, 2 ), vector3( 2, 3, 0 ),
        vector3( 11, 10, 2 ), vector3( 2, 1, 5 ), vector3( 5, 11, 2 );
    /* 119: 0, 1, 2,    4, 5, 6,     */ m_classicCubeCases[119] += vector3( 10, 3, 11 ), vector3( 3, 5, 11 );
    /* 120:          3, 4, 5, 6,     */ m_classicCubeCases[120] += vector3( 5, 3, 7 ), vector3( 10, 4, 6 ),
        vector3( 6, 4, 2 ), vector3( 10, 6, 11 );
    /* 121: 0,       3, 4, 5, 6,     */ m_classicCubeCases[121] += vector3( 6, 0, 1 ), vector3( 11, 0, 6 ),
        vector3( 4, 0, 11 ), vector3( 11, 10, 4 ), vector3( 5, 3, 7 );
    /* 122:    1,    3, 4, 5, 6,     */ m_classicCubeCases[122] += vector3( 2, 6, 10 ), vector3( 10, 6, 11 ),
        vector3( 0, 2, 10 ), vector3( 10, 7, 5 ), vector3( 5, 0, 10 );
    /* 123: 0, 1,    3, 4, 5, 6,     */ m_classicCubeCases[123] += vector3( 6, 11, 1 ), vector3( 1, 11, 10 ),
        vector3( 7, 5, 1 ), vector3( 10, 7, 1 );
    /* 124:       2, 3, 4, 5, 6,     */ m_classicCubeCases[124] += vector3( 1, 3, 11 ), vector3( 11, 3, 7 ),
        vector3( 2, 1, 11 ), vector3( 11, 10, 4 ), vector3( 4, 2, 11 );
    /* 125: 0,    2, 3, 4, 5, 6,     */ m_classicCubeCases[125] += vector3( 3, 7, 0 ), vector3( 0, 7, 11 ),
        vector3( 10, 4, 0 ), vector3( 11, 10, 0 );
    /* 126:    1, 2, 3, 4, 5, 6,     */ m_classicCubeCases[126] += vector3( 1, 0, 2 ), vector3( 11, 10, 7 );
    /* 127: 0, 1, 2, 3, 4, 5, 6,     */ m_classicCubeCases[127] += vector3( 10, 7, 11 );
    /* 128:                      7,  */ m_classicCubeCases[128] += vector3( 11, 7, 10 );
    /* 129: 0,                   7,  */ m_classicCubeCases[129] += vector3( 2, 0, 1 ), vector3( 7, 10, 11 );
    /* 130:    1,                7,  */ m_classicCubeCases[130] += vector3( 0, 4, 3 ), vector3( 7, 10, 11 );
    /* 131: 0, 1,                7,  */ m_classicCubeCases[131] += vector3( 2, 3, 1 ), vector3( 4, 3, 2 ),
        vector3( 7, 10, 11 );
    /* 132:       2,             7,  */ m_classicCubeCases[132] += vector3( 1, 5, 6 ), vector3( 11, 7, 10 );
    /* 133: 0,    2,             7,  */ m_classicCubeCases[133] += vector3( 0, 6, 2 ), vector3( 5, 6, 0 ),
        vector3( 11, 7, 10 );
    /* 134:    1, 2,             7,  */ m_classicCubeCases[134] += vector3( 3, 0, 4 ), vector3( 1, 5, 6 ),
        vector3( 7, 10, 11 );
    /* 135: 0, 1, 2,             7,  */ m_classicCubeCases[135] += vector3( 7, 10, 11 ), vector3( 4, 3, 5 ),
        vector3( 6, 4, 5 ), vector3( 2, 4, 6 );
    /* 136:          3,          7,  */ m_classicCubeCases[136] += vector3( 11, 3, 10 ), vector3( 11, 5, 3 );
    /* 137: 0,       3,          7,  */ m_classicCubeCases[137] += vector3( 11, 3, 10 ), vector3( 5, 3, 11 ),
        vector3( 0, 1, 2 );
    /* 138:    1,    3,          7,  */ m_classicCubeCases[138] += vector3( 11, 4, 10 ), vector3( 0, 4, 11 ),
        vector3( 5, 0, 11 );
    /* 139: 0, 1,    3,          7,  */ m_classicCubeCases[139] += vector3( 4, 10, 2 ), vector3( 2, 10, 5 ),
        vector3( 5, 10, 11 ), vector3( 5, 1, 2 );
    /* 140:       2, 3,          7,  */ m_classicCubeCases[140] += vector3( 1, 11, 6 ), vector3( 10, 11, 1 ),
        vector3( 3, 10, 1 );
    /* 141: 0,    2, 3,          7,  */ m_classicCubeCases[141] += vector3( 2, 0, 6 ), vector3( 6, 0, 10 ),
        vector3( 10, 0, 3 ), vector3( 6, 10, 11 );
    /* 142:    1, 2, 3,          7,  */ m_classicCubeCases[142] += vector3( 6, 1, 11 ), vector3( 1, 0, 11 ),
        vector3( 11, 0, 10 ), vector3( 10, 0, 4 );
    /* 143: 0, 1, 2, 3,          7,  */ m_classicCubeCases[143] += vector3( 10, 11, 4 ), vector3( 4, 11, 6 ),
        vector3( 4, 6, 2 );
    /* 144:             4,       7,  */ m_classicCubeCases[144] += vector3( 7, 10, 11 ), vector3( 9, 8, 2 );
    /* 145: 0,          4,       7,  */ m_classicCubeCases[145] += vector3( 1, 8, 0 ), vector3( 9, 8, 1 ),
        vector3( 10, 11, 7 );
    /* 146:    1,       4,       7,  */ m_classicCubeCases[146] += vector3( 4, 3, 0 ), vector3( 7, 10, 11 ),
        vector3( 8, 2, 9 );
    /* 147: 0, 1,       4,       7,  */ m_classicCubeCases[147] += vector3( 11, 7, 10 ), vector3( 4, 3, 9 ),
        vector3( 9, 3, 1 ), vector3( 4, 9, 8 );
    /* 148:       2,    4,       7,  */ m_classicCubeCases[148] += vector3( 6, 1, 5 ), vector3( 2, 9, 8 ),
        vector3( 11, 7, 10 );
    /* 149: 0,    2,    4,       7,  */ m_classicCubeCases[149] += vector3( 7, 10, 11 ), vector3( 9, 8, 5 ),
        vector3( 5, 8, 0 ), vector3( 9, 5, 6 );
    /* 150:    1, 2,    4,       7,  */ m_classicCubeCases[150] += vector3( 3, 0, 4 ), vector3( 9, 8, 2 ),
        vector3( 1, 5, 6 ), vector3( 7, 10, 11 );
    /* 151: 0, 1, 2,    4,       7,  */ m_classicCubeCases[151] += vector3( 5, 4, 3 ), vector3( 6, 4, 5 ),
        vector3( 8, 4, 6 ), vector3( 6, 9, 8 ), vector3( 7, 10, 11 );
    /* 152:          3, 4,       7,  */ m_classicCubeCases[152] += vector3( 3, 11, 5 ), vector3( 10, 11, 3 ),
        vector3( 9, 8, 2 );
    /* 153: 0,       3, 4,       7,  */ m_classicCubeCases[153] += vector3( 5, 3, 10 ), vector3( 5, 10, 11 ),
        vector3( 0, 1, 8 ), vector3( 8, 1, 9 );
    /* 154:    1,    3, 4,       7,  */ m_classicCubeCases[154] += vector3( 8, 2, 9 ), vector3( 0, 4, 10 ),
        vector3( 11, 0, 10 ), vector3( 5, 0, 11 );
    /* 155: 0, 1,    3, 4,       7,  */ m_classicCubeCases[155] += vector3( 1, 9, 4 ), vector3( 4, 9, 8 ),
        vector3( 5, 1, 4 ), vector3( 4, 10, 11 ), vector3( 11, 5, 4 );
    /* 156:       2, 3, 4,       7,  */ m_classicCubeCases[156] += vector3( 8, 2, 9 ), vector3( 6, 1, 10 ),
        vector3( 10, 1, 3 ), vector3( 6, 10, 11 );
    /* 157: 0,    2, 3, 4,       7,  */ m_classicCubeCases[157] += vector3( 3, 10, 6 ), vector3( 6, 10, 11 ),
        vector3( 0, 3, 6 ), vector3( 6, 9, 8 ), vector3( 8, 0, 6 );
    /* 158:    1, 2, 3, 4,       7,  */ m_classicCubeCases[158] += vector3( 10, 0, 4 ), vector3( 11, 0, 10 ),
        vector3( 1, 0, 11 ), vector3( 11, 6, 1 ), vector3( 8, 2, 9 );
    /* 159: 0, 1, 2, 3, 4,       7,  */ m_classicCubeCases[159] += vector3( 10, 11, 4 ), vector3( 4, 11, 6 ),
        vector3( 9, 8, 4 ), vector3( 6, 9, 4 );
    /* 160:                5,    7,  */ m_classicCubeCases[160] += vector3( 8, 7, 4 ), vector3( 8, 11, 7 );
    /* 161: 0,             5,    7,  */ m_classicCubeCases[161] += vector3( 7, 8, 11 ), vector3( 4, 8, 7 ),
        vector3( 2, 0, 1 );
    /* 162:    1,          5,    7,  */ m_classicCubeCases[162] += vector3( 0, 7, 3 ), vector3( 11, 7, 0 ),
        vector3( 8, 11, 0 );
    /* 163: 0, 1,          5,    7,  */ m_classicCubeCases[163] += vector3( 1, 2, 3 ), vector3( 3, 2, 11 ),
        vector3( 11, 2, 8 ), vector3( 3, 11, 7 );
    /* 164:       2,       5,    7,  */ m_classicCubeCases[164] += vector3( 8, 7, 4 ), vector3( 11, 7, 8 ),
        vector3( 5, 6, 1 );
    /* 165: 0,    2,       5,    7,  */ m_classicCubeCases[165] += vector3( 2, 0, 5 ), vector3( 2, 5, 6 ),
        vector3( 4, 8, 7 ), vector3( 7, 8, 11 );
    /* 166:    1, 2,       5,    7,  */ m_classicCubeCases[166] += vector3( 6, 1, 5 ), vector3( 3, 0, 11 ),
        vector3( 11, 0, 8 ), vector3( 3, 11, 7 );
    /* 167: 0, 1, 2,       5,    7,  */ m_classicCubeCases[167] += vector3( 8, 11, 3 ), vector3( 3, 11, 7 ),
        vector3( 2, 8, 3 ), vector3( 3, 5, 6 ), vector3( 6, 2, 3 );
    /* 168:          3,    5,    7,  */ m_classicCubeCases[168] += vector3( 8, 3, 4 ), vector3( 5, 3, 8 ),
        vector3( 11, 5, 8 );
    /* 169: 0,       3,    5,    7,  */ m_classicCubeCases[169] += vector3( 0, 1, 2 ), vector3( 5, 3, 4 ),
        vector3( 8, 5, 4 ), vector3( 11, 5, 8 );
    /* 170:    1,    3,    5,    7,  */ m_classicCubeCases[170] += vector3( 5, 0, 8 ), vector3( 5, 8, 11 );
    /* 171: 0, 1,    3,    5,    7,  */ m_classicCubeCases[171] += vector3( 1, 2, 5 ), vector3( 5, 2, 8 ),
        vector3( 5, 8, 11 );
    /* 172:       2, 3,    5,    7,  */ m_classicCubeCases[172] += vector3( 11, 4, 8 ), vector3( 1, 4, 11 ),
        vector3( 3, 4, 1 ), vector3( 11, 6, 1 );
    /* 173: 0,    2, 3,    5,    7,  */ m_classicCubeCases[173] += vector3( 6, 2, 3 ), vector3( 3, 2, 0 ),
        vector3( 11, 6, 3 ), vector3( 3, 4, 8 ), vector3( 8, 11, 3 );
    /* 174:    1, 2, 3,    5,    7,  */ m_classicCubeCases[174] += vector3( 6, 1, 11 ), vector3( 11, 1, 0 ),
        vector3( 11, 0, 8 );
    /* 175: 0, 1, 2, 3,    5,    7,  */ m_classicCubeCases[175] += vector3( 8, 11, 2 ), vector3( 11, 6, 2 );
    /* 176:             4, 5,    7,  */ m_classicCubeCases[176] += vector3( 7, 9, 11 ), vector3( 2, 9, 7 ),
        vector3( 4, 2, 7 );
    /* 177: 0,          4, 5,    7,  */ m_classicCubeCases[177] += vector3( 9, 0, 1 ), vector3( 7, 0, 9 ),
        vector3( 4, 0, 7 ), vector3( 9, 11, 7 );
    /* 178:    1,       4, 5,    7,  */ m_classicCubeCases[178] += vector3( 11, 7, 9 ), vector3( 7, 3, 9 ),
        vector3( 9, 3, 2 ), vector3( 2, 3, 0 );
    /* 179: 0, 1,       4, 5,    7,  */ m_classicCubeCases[179] += vector3( 11, 7, 9 ), vector3( 9, 7, 3 ),
        vector3( 9, 3, 1 );
    /* 180:       2,    4, 5,    7,  */ m_classicCubeCases[180] += vector3( 1, 5, 6 ), vector3( 11, 7, 2 ),
        vector3( 2, 7, 4 ), vector3( 11, 2, 9 );
    /* 181: 0,    2,    4, 5,    7,  */ m_classicCubeCases[181] += vector3( 0, 5, 9 ), vector3( 9, 5, 6 ),
        vector3( 4, 0, 9 ), vector3( 9, 11, 7 ), vector3( 7, 4, 9 );
    /* 182:    1, 2,    4, 5,    7,  */ m_classicCubeCases[182] += vector3( 2, 3, 0 ), vector3( 9, 3, 2 ),
        vector3( 7, 3, 9 ), vector3( 9, 11, 7 ), vector3( 1, 5, 6 );
    /* 183: 0, 1, 2,    4, 5,    7,  */ m_classicCubeCases[183] += vector3( 5, 6, 3 ), vector3( 3, 6, 9 ),
        vector3( 11, 7, 3 ), vector3( 9, 11, 3 );
    /* 184:          3, 4, 5,    7,  */ m_classicCubeCases[184] += vector3( 5, 3, 11 ), vector3( 11, 3, 2 ),
        vector3( 2, 3, 4 ), vector3( 11, 2, 9 );
    /* 185: 0,       3, 4, 5,    7,  */ m_classicCubeCases[185] += vector3( 11, 5, 4 ), vector3( 4, 5, 3 ),
        vector3( 9, 11, 4 ), vector3( 4, 0, 1 ), vector3( 1, 9, 4 );
    /* 186:    1,    3, 4, 5,    7,  */ m_classicCubeCases[186] += vector3( 2, 9, 0 ), vector3( 0, 9, 11 ),
        vector3( 0, 11, 5 );
    /* 187: 0, 1,    3, 4, 5,    7,  */ m_classicCubeCases[187] += vector3( 1, 9, 5 ), vector3( 9, 11, 5 );
    /* 188:       2, 3, 4, 5,    7,  */ m_classicCubeCases[188] += vector3( 4, 2, 11 ), vector3( 11, 2, 9 ),
        vector3( 3, 4, 11 ), vector3( 11, 6, 1 ), vector3( 1, 3, 11 );
    /* 189: 0,    2, 3, 4, 5,    7,  */ m_classicCubeCases[189] += vector3( 4, 0, 3 ), vector3( 11, 6, 9 );
    /* 190:    1, 2, 3, 4, 5,    7,  */ m_classicCubeCases[190] += vector3( 2, 9, 0 ), vector3( 0, 9, 11 ),
        vector3( 6, 1, 0 ), vector3( 11, 6, 0 );
    /* 191: 0, 1, 2, 3, 4, 5,    7,  */ m_classicCubeCases[191] += vector3( 6, 9, 11 );
    /* 192:                   6, 7,  */ m_classicCubeCases[192] += vector3( 10, 6, 7 ), vector3( 10, 9, 6 );
    /* 193: 0,                6, 7,  */ m_classicCubeCases[193] += vector3( 10, 6, 7 ), vector3( 9, 6, 10 ),
        vector3( 1, 2, 0 );
    /* 194:    1,             6, 7,  */ m_classicCubeCases[194] += vector3( 6, 10, 9 ), vector3( 7, 10, 6 ),
        vector3( 4, 3, 0 );
    /* 195: 0, 1,             6, 7,  */ m_classicCubeCases[195] += vector3( 9, 7, 10 ), vector3( 6, 7, 9 ),
        vector3( 2, 4, 3 ), vector3( 1, 2, 3 );
    /* 196:       2,          6, 7,  */ m_classicCubeCases[196] += vector3( 10, 5, 7 ), vector3( 1, 5, 10 ),
        vector3( 9, 1, 10 );
    /* 197: 0,    2,          6, 7,  */ m_classicCubeCases[197] += vector3( 5, 2, 0 ), vector3( 10, 2, 5 ),
        vector3( 9, 2, 10 ), vector3( 5, 7, 10 );
    /* 198:    1, 2,          6, 7,  */ m_classicCubeCases[198] += vector3( 0, 4, 3 ), vector3( 7, 10, 1 ),
        vector3( 1, 10, 9 ), vector3( 7, 1, 5 );
    /* 199: 0, 1, 2,          6, 7,  */ m_classicCubeCases[199] += vector3( 2, 4, 5 ), vector3( 5, 4, 3 ),
        vector3( 9, 2, 5 ), vector3( 5, 7, 10 ), vector3( 10, 9, 5 );
    /* 200:          3,       6, 7,  */ m_classicCubeCases[200] += vector3( 3, 6, 5 ), vector3( 9, 6, 3 ),
        vector3( 10, 9, 3 );
    /* 201: 0,       3,       6, 7,  */ m_classicCubeCases[201] += vector3( 2, 0, 1 ), vector3( 5, 3, 9 ),
        vector3( 9, 3, 10 ), vector3( 5, 9, 6 );
    /* 202:    1,    3,       6, 7,  */ m_classicCubeCases[202] += vector3( 9, 4, 10 ), vector3( 5, 4, 9 ),
        vector3( 0, 4, 5 ), vector3( 6, 5, 9 );
    /* 203: 0, 1,    3,       6, 7,  */ m_classicCubeCases[203] += vector3( 10, 9, 5 ), vector3( 5, 9, 6 ),
        vector3( 4, 10, 5 ), vector3( 5, 1, 2 ), vector3( 2, 4, 5 );
    /* 204:       2, 3,       6, 7,  */ m_classicCubeCases[204] += vector3( 1, 3, 10 ), vector3( 9, 1, 10 );
    /* 205: 0,    2, 3,       6, 7,  */ m_classicCubeCases[205] += vector3( 2, 0, 9 ), vector3( 9, 0, 3 ),
        vector3( 9, 3, 10 );
    /* 206:    1, 2, 3,       6, 7,  */ m_classicCubeCases[206] += vector3( 0, 4, 1 ), vector3( 1, 4, 10 ),
        vector3( 1, 10, 9 );
    /* 207: 0, 1, 2, 3,       6, 7,  */ m_classicCubeCases[207] += vector3( 2, 4, 9 ), vector3( 4, 10, 9 );
    /* 208:             4,    6, 7,  */ m_classicCubeCases[208] += vector3( 2, 10, 8 ), vector3( 7, 10, 2 ),
        vector3( 6, 7, 2 );
    /* 209: 0,          4,    6, 7,  */ m_classicCubeCases[209] += vector3( 0, 10, 8 ), vector3( 6, 10, 0 ),
        vector3( 7, 10, 6 ), vector3( 1, 6, 0 );
    /* 210:    1,       4,    6, 7,  */ m_classicCubeCases[210] += vector3( 3, 0, 4 ), vector3( 8, 2, 7 ),
        vector3( 7, 2, 6 ), vector3( 8, 7, 10 );
    /* 211: 0, 1,       4,    6, 7,  */ m_classicCubeCases[211] += vector3( 6, 7, 8 ), vector3( 8, 7, 10 ),
        vector3( 1, 6, 8 ), vector3( 8, 4, 3 ), vector3( 3, 1, 8 );
    /* 212:       2,    4,    6, 7,  */ m_classicCubeCases[212] += vector3( 10, 5, 7 ), vector3( 10, 1, 5 ),
        vector3( 8, 1, 10 ), vector3( 2, 1, 8 );
    /* 213: 0,    2,    4,    6, 7,  */ m_classicCubeCases[213] += vector3( 7, 10, 5 ), vector3( 5, 10, 8 ),
        vector3( 5, 8, 0 );
    /* 214:    1, 2,    4,    6, 7,  */ m_classicCubeCases[214] += vector3( 7, 1, 5 ), vector3( 10, 1, 7 ),
        vector3( 2, 1, 10 ), vector3( 10, 8, 2 ), vector3( 3, 0, 4 );
    /* 215: 0, 1, 2,    4,    6, 7,  */ m_classicCubeCases[215] += vector3( 7, 10, 5 ), vector3( 5, 10, 8 ),
        vector3( 4, 3, 5 ), vector3( 8, 4, 5 );
    /* 216:          3, 4,    6, 7,  */ m_classicCubeCases[216] += vector3( 10, 5, 3 ), vector3( 2, 5, 10 ),
        vector3( 6, 5, 2 ), vector3( 10, 8, 2 );
    /* 217: 0,       3, 4,    6, 7,  */ m_classicCubeCases[217] += vector3( 8, 0, 6 ), vector3( 6, 0, 1 ),
        vector3( 10, 8, 6 ), vector3( 6, 5, 3 ), vector3( 3, 10, 6 );
    /* 218:    1,    3, 4,    6, 7,  */ m_classicCubeCases[218] += vector3( 5, 0, 10 ), vector3( 10, 0, 4 ),
        vector3( 6, 5, 10 ), vector3( 10, 8, 2 ), vector3( 2, 6, 10 );
    /* 219: 0, 1,    3, 4,    6, 7,  */ m_classicCubeCases[219] += vector3( 8, 4, 10 ), vector3( 6, 5, 1 );
    /* 220:       2, 3, 4,    6, 7,  */ m_classicCubeCases[220] += vector3( 8, 2, 10 ), vector3( 10, 2, 1 ),
        vector3( 10, 1, 3 );
    /* 221: 0,    2, 3, 4,    6, 7,  */ m_classicCubeCases[221] += vector3( 8, 0, 10 ), vector3( 0, 3, 10 );
    /* 222:    1, 2, 3, 4,    6, 7,  */ m_classicCubeCases[222] += vector3( 8, 2, 10 ), vector3( 10, 2, 1 ),
        vector3( 0, 4, 10 ), vector3( 1, 0, 10 );
    /* 223: 0, 1, 2, 3, 4,    6, 7,  */ m_classicCubeCases[223] += vector3( 8, 4, 10 );
    /* 224:                5, 6, 7,  */ m_classicCubeCases[224] += vector3( 6, 8, 9 ), vector3( 4, 8, 6 ),
        vector3( 7, 4, 6 );
    /* 225: 0,             5, 6, 7,  */ m_classicCubeCases[225] += vector3( 2, 0, 1 ), vector3( 4, 8, 9 ),
        vector3( 6, 4, 9 ), vector3( 7, 4, 6 );
    /* 226:    1,          5, 6, 7,  */ m_classicCubeCases[226] += vector3( 7, 3, 6 ), vector3( 6, 3, 8 ),
        vector3( 8, 3, 0 ), vector3( 8, 9, 6 );
    /* 227: 0, 1,          5, 6, 7,  */ m_classicCubeCases[227] += vector3( 3, 1, 8 ), vector3( 8, 1, 2 ),
        vector3( 7, 3, 8 ), vector3( 8, 9, 6 ), vector3( 6, 7, 8 );
    /* 228:       2,       5, 6, 7,  */ m_classicCubeCases[228] += vector3( 4, 5, 7 ), vector3( 9, 5, 4 ),
        vector3( 1, 5, 9 ), vector3( 8, 9, 4 );
    /* 229: 0,    2,       5, 6, 7,  */ m_classicCubeCases[229] += vector3( 7, 4, 9 ), vector3( 9, 4, 8 ),
        vector3( 5, 7, 9 ), vector3( 9, 2, 0 ), vector3( 0, 5, 9 );
    /* 230:    1, 2,       5, 6, 7,  */ m_classicCubeCases[230] += vector3( 9, 1, 7 ), vector3( 7, 1, 5 ),
        vector3( 8, 9, 7 ), vector3( 7, 3, 0 ), vector3( 0, 8, 7 );
    /* 231: 0, 1, 2,       5, 6, 7,  */ m_classicCubeCases[231] += vector3( 7, 3, 5 ), vector3( 9, 2, 8 );
    /* 232:          3,    5, 6, 7,  */ m_classicCubeCases[232] += vector3( 6, 8, 9 ), vector3( 6, 4, 8 ),
        vector3( 5, 4, 6 ), vector3( 3, 4, 5 );
    /* 233: 0,       3,    5, 6, 7,  */ m_classicCubeCases[233] += vector3( 9, 4, 8 ), vector3( 6, 4, 9 ),
        vector3( 3, 4, 6 ), vector3( 6, 5, 3 ), vector3( 2, 0, 1 );
    /* 234:    1,    3,    5, 6, 7,  */ m_classicCubeCases[234] += vector3( 9, 6, 8 ), vector3( 8, 6, 5 ),
        vector3( 8, 5, 0 );
    /* 235: 0, 1,    3,    5, 6, 7,  */ m_classicCubeCases[235] += vector3( 9, 6, 8 ), vector3( 8, 6, 5 ),
        vector3( 1, 2, 8 ), vector3( 5, 1, 8 );
    /* 236:       2, 3,    5, 6, 7,  */ m_classicCubeCases[236] += vector3( 4, 8, 3 ), vector3( 3, 8, 9 ),
        vector3( 3, 9, 1 );
    /* 237: 0,    2, 3,    5, 6, 7,  */ m_classicCubeCases[237] += vector3( 4, 8, 3 ), vector3( 3, 8, 9 ),
        vector3( 2, 0, 3 ), vector3( 9, 2, 3 );
    /* 238:    1, 2, 3,    5, 6, 7,  */ m_classicCubeCases[238] += vector3( 0, 8, 1 ), vector3( 8, 9, 1 );
    /* 239: 0, 1, 2, 3,    5, 6, 7,  */ m_classicCubeCases[239] += vector3( 2, 8, 9 );
    /* 240:             4, 5, 6, 7,  */ m_classicCubeCases[240] += vector3( 7, 4, 2 ), vector3( 6, 7, 2 );
    /* 241: 0,          4, 5, 6, 7,  */ m_classicCubeCases[241] += vector3( 0, 1, 4 ), vector3( 4, 1, 6 ),
        vector3( 4, 6, 7 );
    /* 242:    1,       4, 5, 6, 7,  */ m_classicCubeCases[242] += vector3( 3, 0, 7 ), vector3( 7, 0, 2 ),
        vector3( 7, 2, 6 );
    /* 243: 0, 1,       4, 5, 6, 7,  */ m_classicCubeCases[243] += vector3( 3, 1, 7 ), vector3( 1, 6, 7 );
    /* 244:       2,    4, 5, 6, 7,  */ m_classicCubeCases[244] += vector3( 1, 5, 2 ), vector3( 2, 5, 7 ),
        vector3( 2, 7, 4 );
    /* 245: 0,    2,    4, 5, 6, 7,  */ m_classicCubeCases[245] += vector3( 7, 4, 5 ), vector3( 4, 0, 5 );
    /* 246:    1, 2,    4, 5, 6, 7,  */ m_classicCubeCases[246] += vector3( 1, 5, 2 ), vector3( 2, 5, 7 ),
        vector3( 3, 0, 2 ), vector3( 7, 3, 2 );
    /* 247: 0, 1, 2,    4, 5, 6, 7,  */ m_classicCubeCases[247] += vector3( 7, 3, 5 );
    /* 248:          3, 4, 5, 6, 7,  */ m_classicCubeCases[248] += vector3( 5, 3, 6 ), vector3( 6, 3, 4 ),
        vector3( 6, 4, 2 );
    /* 249: 0,       3, 4, 5, 6, 7,  */ m_classicCubeCases[249] += vector3( 0, 1, 4 ), vector3( 4, 1, 6 ),
        vector3( 5, 3, 4 ), vector3( 6, 5, 4 );
    /* 250:    1,    3, 4, 5, 6, 7,  */ m_classicCubeCases[250] += vector3( 5, 0, 6 ), vector3( 0, 2, 6 );
    /* 251: 0, 1,    3, 4, 5, 6, 7,  */ m_classicCubeCases[251] += vector3( 5, 1, 6 );
    /* 252:       2, 3, 4, 5, 6, 7,  */ m_classicCubeCases[252] += vector3( 1, 3, 2 ), vector3( 3, 4, 2 );
    /* 253: 0,    2, 3, 4, 5, 6, 7,  */ m_classicCubeCases[253] += vector3( 4, 0, 3 );
    /* 254:    1, 2, 3, 4, 5, 6, 7,  */ m_classicCubeCases[254] += vector3( 1, 0, 2 );
    /* 255: 0, 1, 2, 3, 4, 5, 6, 7,  */
}

bool marching_cubes_table::test_face( char face, const float cubeCorners[8] ) const {
    double A, B, C, D;

    switch( face ) {
    case -1:
    case 1:
        A = cubeCorners[0];
        B = cubeCorners[1];
        C = cubeCorners[5];
        D = cubeCorners[4];
        break;
    case -2:
    case 2:
        A = cubeCorners[1];
        B = cubeCorners[3];
        C = cubeCorners[7];
        D = cubeCorners[5];
        break;
    case -3:
    case 3:
        A = cubeCorners[3];
        B = cubeCorners[2];
        C = cubeCorners[6];
        D = cubeCorners[7];
        break;
    case -4:
    case 4:
        A = cubeCorners[2];
        B = cubeCorners[0];
        C = cubeCorners[4];
        D = cubeCorners[6];
        break;
    case -5:
    case 5:
        A = cubeCorners[0];
        B = cubeCorners[1];
        C = cubeCorners[3];
        D = cubeCorners[2];
        break;
    case -6:
    case 6:
        A = cubeCorners[4];
        B = cubeCorners[5];
        C = cubeCorners[7];
        D = cubeCorners[6];
        break;
    default:
        throw std::runtime_error( "marching_cubes_table::test_face() - Invalid face test code (" +
                                  boost::lexical_cast<std::string>( (int)face ) + ")." );
        break;
    };

    const double result = ( A * C - B * D );
    if( result == 0 ) {
        return face >= 0;
    }

    // 'A' should flip the sign appropriately for the face test being performed
    const double sign = ( A >= 0 ? 1.0 : -1.0 );
    return face * sign * ( result ) >= 0;
}

bool marching_cubes_table::test_interior( char flag, char cubeCase, char config, char subconfig,
                                          const float cubeCorners[8] ) const {

    // temp values for testing
    float t, a, b, At = 0, Bt = 0, Ct = 0, Dt = 0;
    float A0, A1, B0, B1, C0, C1, D0, D1;
    float A1A0, B1B0, C1C0, D1D0;
    char test = 0;

    A0 = cubeCorners[0];
    A1 = cubeCorners[4];
    B0 = cubeCorners[2];
    B1 = cubeCorners[6];
    C0 = cubeCorners[3];
    C1 = cubeCorners[7];
    D0 = cubeCorners[1];
    D1 = cubeCorners[5];

    // for some cases, an edge is required which is stored in the test table.  the edge determines the orientation
    // of the plane used to determine if corner vertices are "connected"
    char edge = -1;

    switch( cubeCase ) {

    // case 4 and 10 have enough symmetry to always use the same corners to test
    case 4:
    case 10:

        // set up temp values for testing
        A1A0 = A1 - A0;
        B1B0 = B1 - B0;
        C1C0 = C1 - C0;
        D1D0 = D1 - D0;

        // find the max/min of Chernyaev's 2nd order equation
        // to determine the height of the plane which connects
        // the corners
        a = A1A0 * C1C0 - B1B0 * D1D0;
        b = C0 * A1A0 + A0 * C1C0 - D0 * B1B0 - B0 * D1D0;
        t = -b / ( 2 * a );

        // if t lies outside the cube, return
        if( t < 0.f || t > 1.f )
            return flag > 0;

        // otherwise, set up the test values for
        At = A0 + A1A0 * t;
        Bt = B0 + B1B0 * t;
        Ct = C0 + C1C0 * t;
        Dt = D0 + D1D0 * t;
        break;

    case 6:
    case 7:
    case 12:
    case 13:
        switch( cubeCase ) {
        case 6:
            edge = test6[config][2];
            break;
        case 7:
            edge = test7[config][4];
            break;
        case 12:
            edge = test12[config][3];
            break;
        case 13:
            edge = triangles13_5_1[config][subconfig][0];
            break;
        }
        switch( edge ) {
        case 0:
            t = A0 / ( A0 - D0 );
            At = 0;
            Bt = B0 + ( C0 - B0 ) * t;
            Ct = B1 + ( C1 - B1 ) * t;
            Dt = A1 + ( D1 - A1 ) * t;
            break;
        case 1:
            t = B0 / ( B0 - A0 );
            At = 0;
            Bt = C0 + ( D0 - C0 ) * t;
            Ct = C1 + ( D1 - C1 ) * t;
            Dt = B1 + ( A1 - B1 ) * t;
            break;
        case 2:
            t = A0 / ( A0 - A1 );
            At = 0;
            Bt = B0 + ( B1 - B0 ) * t;
            Ct = C0 + ( C1 - C0 ) * t;
            Dt = D0 + ( D1 - D0 ) * t;
            break;
        case 3:
            t = D0 / ( D0 - C0 );
            At = 0;
            Bt = A0 + ( B0 - A0 ) * t;
            Ct = A1 + ( B1 - A1 ) * t;
            Dt = D1 + ( C1 - D1 ) * t;
            break;
        case 4:
            t = D0 / ( D0 - D1 );
            At = 0;
            Bt = A0 + ( A1 - A0 ) * t;
            Ct = B0 + ( B1 - B0 ) * t;
            Dt = C0 + ( C1 - C0 ) * t;
            break;
        case 5:
            t = C0 / ( C0 - B0 );
            At = 0;
            Bt = D0 + ( A0 - D0 ) * t;
            Ct = D1 + ( A1 - D1 ) * t;
            Dt = C1 + ( B1 - C1 ) * t;
            break;
        case 6:
            t = B0 / ( B0 - B1 );
            At = 0;
            Bt = C0 + ( C1 - C0 ) * t;
            Ct = D0 + ( D1 - D0 ) * t;
            Dt = A0 + ( A1 - A0 ) * t;
            break;
        case 7:
            t = C0 / ( C0 - C1 );
            At = 0;
            Bt = D0 + ( D1 - D0 ) * t;
            Ct = A0 + ( A1 - A0 ) * t;
            Dt = B0 + ( B1 - B0 ) * t;
            break;
        case 8:
            t = A1 / ( A1 - D1 );
            At = 0;
            Bt = B1 + ( C1 - B1 ) * t;
            Ct = B0 + ( C0 - B0 ) * t;
            Dt = A0 + ( D0 - A0 ) * t;
            break;
        case 9:
            t = B1 / ( B1 - A1 );
            At = 0;
            Bt = C1 + ( D1 - C1 ) * t;
            Ct = C0 + ( D0 - C0 ) * t;
            Dt = B0 + ( A0 - B0 ) * t;
            break;
        case 10:
            t = D1 / ( D1 - C1 );
            At = 0;
            Bt = A1 + ( B1 - A1 ) * t;
            Ct = A0 + ( B0 - A0 ) * t;
            Dt = D0 + ( C0 - D0 ) * t;
            break;
        case 11:
            t = C1 / ( C1 - B1 );
            At = 0;
            Bt = D1 + ( A1 - D1 ) * t;
            Ct = D0 + ( A0 - D0 ) * t;
            Dt = C0 + ( B0 - C0 ) * t;
            break;
        default:
            throw std::runtime_error( "marching_cubes_table::test_interior() - Invalid edge code (" +
                                      boost::lexical_cast<std::string>( (int)edge ) + ")." );
        }
        break;

    default:
        throw std::runtime_error( "marching_cubes_table::test_interior() - Invalid ambiguous case (" +
                                  boost::lexical_cast<std::string>( (int)cubeCase ) + ")." );
        break;
    }

    if( At >= 0 )
        ++test;
    if( Bt >= 0 )
        test += 2;
    if( Ct >= 0 )
        test += 4;
    if( Dt >= 0 )
        test += 8;

    // normalize the result using the largest input value
    float result =
        ( At * Ct - Bt * Dt ) / std::max( std::max( fabs( At ), fabs( Bt ) ), std::max( fabs( Ct ), fabs( Dt ) ) );

    // the only cases which remain ambiguous after the intermediate plane corner vertex values are
    // examined are those where opposing vertices of the cube share the same sign.  ie cases 5, 10.
    // the rest can be explained away using the plane corner vertex configuration.
    switch( test ) {
    case 0:
        return flag > 0;
    case 1:
        return flag > 0;
    case 2:
        return flag > 0;
    case 3:
        return flag > 0;
    case 4:
        return flag > 0;
    case 5:
        if( result < 10e-5f )
            return flag > 0;
        break;
    case 6:
        return flag > 0;
    case 7:
        return flag < 0;
    case 8:
        return flag > 0;
    case 9:
        return flag > 0;
    case 10:
        if( result >= 10e-5f )
            return flag > 0;
        break;
    case 11:
        return flag < 0;
    case 12:
        return flag > 0;
    case 13:
        return flag < 0;
    case 14:
        return flag < 0;
    case 15:
        return flag < 0;
    }
    return flag < 0;
}

void marching_cubes_table::get_faces_from_table( const char* table, int n, std::vector<vector3>& outFaces ) const {
    outFaces.resize( n );
    for( int i = 0; i < n; ++i ) {
        outFaces[i].set( table[3 * i], table[3 * i + 1], table[3 * i + 2] );
    }
}

void marching_cubes_table::dump_cube( const float cubeCorners[8], std::ostream& o ) const {
    o << cubeCorners[0] << " " << cubeCorners[1] << " " << cubeCorners[2] << " " << cubeCorners[3] << " "
      << cubeCorners[4] << " " << cubeCorners[5] << " " << cubeCorners[6] << " " << cubeCorners[7] << std::endl;
}

const std::vector<vector3>& marching_cubes_table::get_cubecase_faces( unsigned char cubeCase ) const {
    return m_classicCubeCases[cubeCase];
}

// bool marching_cubes_table::get_cubecase_faces( unsigned char cubeCase, std::vector<float>& cubeCorners,
// std::vector<frantic::graphics::vector3>& faces, frantic::graphics::vector3& //cubeCaseDebug ) const {
bool marching_cubes_table::get_cubecase_faces( unsigned char cubeCase, const float cubeCorners[8],
                                               std::vector<frantic::graphics::vector3>& faces ) const {

    // cubeCaseDebug.set(-1,-1,-1);

    bool extraVert = false;
    const unsigned char uniqueCubeCase = disambiguatedCubeCases[cubeCase][0];
    const unsigned char config = disambiguatedCubeCases[cubeCase][1];
    unsigned char subconfig = 0;

    switch( uniqueCubeCase ) {
    case 1:
        // cubeCaseCounter["1"]++;
        // cubeCaseDebug.set(1,0,0);
        get_faces_from_table( triangles1[config], 1, faces );
        break;

    case 2:
        // cubeCaseCounter["2"]++;
        // cubeCaseDebug.set(2,0,0);
        get_faces_from_table( triangles2[config], 2, faces );
        break;

    case 3:
        if( test_face( test3[config], cubeCorners ) ) {
            // cubeCaseCounter["3.2"]++;
            // cubeCaseDebug.set(3,2,0);
            get_faces_from_table( triangles3_2[config], 4, faces ); // 3.2
        } else {
            // cubeCaseCounter["3.1"]++;
            // cubeCaseDebug.set(3,2,0);
            get_faces_from_table( triangles3_1[config], 2, faces ); // 3.1
        }
        break;

    case 4:
        if( test_interior( test4[config], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
            // cubeCaseCounter["4.1"]++;
            // cubeCaseDebug.set(4,1,0);
            get_faces_from_table( triangles4_1[config], 2, faces ); // 4.1
        } else {
            // cubeCaseCounter["4.2"]++;
            // cubeCaseDebug.set(4,2,0);
            get_faces_from_table( triangles4_2[config], 6, faces ); // 4.2
        }
        break;

    case 5:
        // cubeCaseCounter["5"]++;
        // cubeCaseDebug.set(5,0,0);
        get_faces_from_table( triangles5[config], 3, faces );
        break;

    case 6:
        if( test_face( test6[config][0], cubeCorners ) ) {
            // cubeCaseCounter["6.2"]++;
            // cubeCaseDebug.set(6,2,0);
            // std::cout << "6.2" << std::endl;
            get_faces_from_table( triangles6_2[config], 5, faces ); // 6.2
        } else {
            if( test_interior( test6[config][1], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
                // cubeCaseCounter["6.1.1"]++;
                // cubeCaseDebug.set(6,1,1);
                // std::cout << "6.1.1" << std::endl;
                get_faces_from_table( triangles6_1_1[config], 3, faces ); // 6.1.1
            } else {
                // cubeCaseCounter["6.1.2"]++;
                // cubeCaseDebug.set(6,1,2);
                // std::cout << "6.1.2" << std::endl;
                extraVert = true;
                get_faces_from_table( triangles6_1_2[config], 9, faces ); // 6.1.2
            }
        }
        break;

    case 7:
        if( test_face( test7[config][0], cubeCorners ) )
            subconfig += 1;
        if( test_face( test7[config][1], cubeCorners ) )
            subconfig += 2;
        if( test_face( test7[config][2], cubeCorners ) )
            subconfig += 4;

        switch( subconfig ) {
        case 0:
            // cubeCaseCounter["7.1"]++;
            // cubeCaseDebug.set(7,1,0);
            get_faces_from_table( triangles7_1[config], 3, faces );
            break;
        case 1:
            // cubeCaseCounter["7.2"]++;
            // cubeCaseDebug.set(7,2,0);
            get_faces_from_table( triangles7_2[config][0], 5, faces );
            break;
        case 2:
            // cubeCaseCounter["7.2"]++;
            // cubeCaseDebug.set(7,2,0);
            get_faces_from_table( triangles7_2[config][1], 5, faces );
            break;
        case 3:
            // cubeCaseCounter["7.3"]++;
            // cubeCaseDebug.set(7,3,0);
            extraVert = true;
            get_faces_from_table( triangles7_3[config][0], 9, faces );
            break;
        case 4:
            // cubeCaseCounter["7.2"]++;
            // cubeCaseDebug.set(7,2,0);
            get_faces_from_table( triangles7_2[config][2], 5, faces );
            break;
        case 5:
            // cubeCaseCounter["7.3"]++;
            // cubeCaseDebug.set(7,3,0);
            extraVert = true;
            get_faces_from_table( triangles7_3[config][1], 9, faces );
            break;
        case 6:
            // cubeCaseCounter["7.3"]++;
            // cubeCaseDebug.set(7,3,0);
            extraVert = true;
            get_faces_from_table( triangles7_3[config][2], 9, faces );
            break;
        case 7:
            if( test_interior( test7[config][3], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
                // cubeCaseCounter["7.4.2"]++;
                // cubeCaseDebug.set(7,4,2);
                get_faces_from_table( triangles7_4_2[config], 9, faces );
            } else {
                // cubeCaseCounter["7.4.1"]++;
                // cubeCaseDebug.set(7,4,1);
                get_faces_from_table( triangles7_4_1[config], 5, faces );
            }
            break;
        };
        break;

    case 8:
        // cubeCaseCounter["8"]++;
        // cubeCaseDebug.set(8,0,0);
        get_faces_from_table( triangles8[config], 2, faces );
        break;

    case 9:
        // cubeCaseCounter["9"]++;
        // cubeCaseDebug.set(9,0,0);
        get_faces_from_table( triangles9[config], 4, faces );
        break;

    case 10:
        if( test_face( test10[config][0], cubeCorners ) ) {
            if( test_face( test10[config][1], cubeCorners ) ) {
                // cubeCaseCounter["10.1.1"]++;
                // cubeCaseDebug.set(10,1,1);
                get_faces_from_table( triangles10_1_1_[config], 4, faces ); // 10.1.1
            } else {
                // cubeCaseCounter["10.2"]++;
                // cubeCaseDebug.set(10,2,0);
                extraVert = true;
                get_faces_from_table( triangles10_2[config], 8, faces ); // 10.2
            }
        } else {
            if( test_face( test10[config][1], cubeCorners ) ) {
                // cubeCaseCounter["10.2"]++;
                // cubeCaseDebug.set(10,2,0);
                extraVert = true;
                get_faces_from_table( triangles10_2_[config], 8, faces ); // 10.2
            } else {
                if( test_interior( test10[config][2], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
                    // cubeCaseCounter["10.1.1"]++;
                    // cubeCaseDebug.set(10,1,1);
                    get_faces_from_table( triangles10_1_1[config], 4, faces ); // 10.1.1
                } else {
                    // cubeCaseCounter["10.1.2"]++;
                    // cubeCaseDebug.set(10,1,2);
                    get_faces_from_table( triangles10_1_2[config], 8, faces ); // 10.1.2
                }
            }
        }
        break;

    case 11:
        // cubeCaseCounter["11"]++;
        // cubeCaseDebug.set(11,0,0);
        get_faces_from_table( triangles11[config], 4, faces );
        break;

    case 12:
        if( test_face( test12[config][0], cubeCorners ) ) {
            if( test_face( test12[config][1], cubeCorners ) ) {
                // cubeCaseCounter["12.1.1"]++;
                // cubeCaseDebug.set(12,1,1);
                // std::cout << "12.1.1" << std::endl;
                get_faces_from_table( triangles12_1_1_[config], 4, faces ); // 12.1.1
            } else {
                // cubeCaseCounter["12.2"]++;
                // cubeCaseDebug.set(12,2,0);
                // std::cout << "12.2" << std::endl;
                extraVert = true;
                get_faces_from_table( triangles12_2[config], 8, faces ); // 12.2
            }
        } else {
            if( test_face( test12[config][1], cubeCorners ) ) {
                // cubeCaseCounter["12.2"]++;
                // cubeCaseDebug.set(12,2,0);
                // std::cout << "12.2" << std::endl;
                extraVert = true;
                get_faces_from_table( triangles12_2_[config], 8, faces ); // 12.2
            } else {
                if( test_interior( test12[config][2], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
                    // cubeCaseCounter["12.1.1"]++;
                    // cubeCaseDebug.set(12,1,1);
                    // std::cout << "12.1.1" << std::endl;
                    get_faces_from_table( triangles12_1_1[config], 4, faces ); // 12.1.1
                } else {
                    // cubeCaseCounter["12.1.2"]++;
                    // cubeCaseDebug.set(12,1,2);
                    // std::cout << "12.1.2" << std::endl;
                    get_faces_from_table( triangles12_1_2[config], 8, faces ); // 12.1.2
                }
            }
        }
        break;

    case 13:
        if( test_face( test13[config][0], cubeCorners ) )
            subconfig += 1;
        if( test_face( test13[config][1], cubeCorners ) )
            subconfig += 2;
        if( test_face( test13[config][2], cubeCorners ) )
            subconfig += 4;
        if( test_face( test13[config][3], cubeCorners ) )
            subconfig += 8;
        if( test_face( test13[config][4], cubeCorners ) )
            subconfig += 16;
        if( test_face( test13[config][5], cubeCorners ) )
            subconfig += 32;

        switch( subconfig ) {
        // subcases for 13.1
        case 0:
            // cubeCaseCounter["13.1"]++;
            // cubeCaseDebug.set(13,1,0);
            get_faces_from_table( triangles13_1[config], 4, faces );
            break;
        case 63:
            // cubeCaseCounter["13.1"]++;
            // cubeCaseDebug.set(13,1,0);
            get_faces_from_table( triangles13_1_[config], 4, faces );
            break;

        // subcases for 13.2, with 1 positive face test
        case 1:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2[config][0], 6, faces );
            break;
        case 2:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2[config][1], 6, faces );
            break;
        case 4:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2[config][2], 6, faces );
            break;
        case 8:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2[config][3], 6, faces );
            break;
        case 16:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2[config][4], 6, faces );
            break;
        case 32:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2[config][5], 6, faces );
            break;

        // subcases for 13.2, with all but 1 positive face test
        case 62:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2_[config][0], 6, faces );
            break;
        case 61:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2_[config][1], 6, faces );
            break;
        case 59:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2_[config][2], 6, faces );
            break;
        case 55:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2_[config][3], 6, faces );
            break;
        case 47:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2_[config][4], 6, faces );
            break;
        case 31:
            // cubeCaseCounter["13.2"]++;
            // cubeCaseDebug.set(13,2,0);
            get_faces_from_table( triangles13_2_[config][5], 6, faces );
            break;

        // subcases for 13.3, with 2 positive face tests
        case 3:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][0], 10, faces );
            break;
        case 6:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][4], 10, faces );
            break;
        case 9:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][1], 10, faces );
            break;
        case 12:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][7], 10, faces );
            break;
        case 17:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][2], 10, faces );
            break;
        case 18:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][5], 10, faces );
            break;
        case 20:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][8], 10, faces );
            break;
        case 24:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][10], 10, faces );
            break;
        case 33:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][3], 10, faces );
            break;
        case 34:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][6], 10, faces );
            break;
        case 36:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][9], 10, faces );
            break;
        case 40:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3[config][11], 10, faces );
            break;

        // subcases for 13.3, with all but 2 positive face tests
        case 60:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][0], 10, faces );
            break;
        case 57:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][4], 10, faces );
            break;
        case 54:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][1], 10, faces );
            break;
        case 51:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][7], 10, faces );
            break;
        case 46:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][2], 10, faces );
            break;
        case 45:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][5], 10, faces );
            break;
        case 43:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][8], 10, faces );
            break;
        case 39:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][10], 10, faces );
            break;
        case 30:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][3], 10, faces );
            break;
        case 29:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][6], 10, faces );
            break;
        case 27:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][9], 10, faces );
            break;
        case 23:
            // cubeCaseCounter["13.3"]++;
            // cubeCaseDebug.set(13,3,0);
            extraVert = true;
            get_faces_from_table( triangles13_3_[config][11], 10, faces );
            break;

        // subcases for 13.4, 3 positive face tests
        case 22:
            // cubeCaseCounter["13.4"]++;
            // cubeCaseDebug.set(13,4,0);
            extraVert = true;
            get_faces_from_table( triangles13_4[config][2], 12, faces );
            break;
        case 25:
            // cubeCaseCounter["13.4"]++;
            // cubeCaseDebug.set(13,4,0);
            extraVert = true;
            get_faces_from_table( triangles13_4[config][1], 12, faces );
            break;
        case 44:
            // cubeCaseCounter["13.4"]++;
            // cubeCaseDebug.set(13,4,0);
            extraVert = true;
            get_faces_from_table( triangles13_4[config][3], 12, faces );
            break;
        case 35:
            // cubeCaseCounter["13.4"]++;
            // cubeCaseDebug.set(13,4,0);
            extraVert = true;
            get_faces_from_table( triangles13_4[config][0], 12, faces );
            break;

        // subcases for case 13.5, 3 positive face tests configurations requiring an internal test
        case 19:
            subconfig = 0;
            if( test_interior( test13[config][6], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
                // cubeCaseCounter["13.5.1"]++;
                // cubeCaseDebug.set(13,5,1);
                get_faces_from_table( triangles13_5_1[config][0], 6, faces );
            } else {
                // cubeCaseCounter["13.5.2"]++;
                // cubeCaseDebug.set(13,5,2);
                get_faces_from_table( triangles13_5_2[config][0], 10, faces );
            }
            break;
        case 28:
            subconfig = 3;
            if( test_interior( test13[config][6], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
                // cubeCaseCounter["13.5.1"]++;
                // cubeCaseDebug.set(13,5,1);
                get_faces_from_table( triangles13_5_1[config][3], 6, faces );
            } else {
                // cubeCaseCounter["13.5.2"]++;
                // cubeCaseDebug.set(13,5,2);
                get_faces_from_table( triangles13_5_2[config][3], 10, faces );
            }
            break;
        case 38:
            subconfig = 2;
            if( test_interior( test13[config][6], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
                // cubeCaseCounter["13.5.1"]++;
                // cubeCaseDebug.set(13,5,1);
                get_faces_from_table( triangles13_5_1[config][2], 6, faces );
            } else {
                // cubeCaseCounter["13.5.2"]++;
                // cubeCaseDebug.set(13,5,2);
                get_faces_from_table( triangles13_5_2[config][2], 10, faces );
            }
            break;
        case 41:
            subconfig = 1;
            if( test_interior( test13[config][6], uniqueCubeCase, config, subconfig, cubeCorners ) ) {
                // cubeCaseCounter["13.5.1"]++;
                // cubeCaseDebug.set(13,5,1);
                get_faces_from_table( triangles13_5_1[config][1], 6, faces );
            } else {
                // cubeCaseCounter["13.5.2"]++;
                // cubeCaseDebug.set(13,5,2);
                get_faces_from_table( triangles13_5_2[config][1], 10, faces );
            }
            break;

        default: {
            std::stringstream ss;
            ss << "marching_cubes_table::get_cubecase_faces(): Encountered an unknown cubeCorners case.  Corner values "
                  "as "
                  "follows:"
               << std::endl;
            dump_cube( cubeCorners, ss );
            throw std::runtime_error( ss.str() );
        }
        }
        break;

    case 14:
        // cubeCaseCounter["14"]++;
        // cubeCaseDebug.set(14,0,0);
        get_faces_from_table( triangles14[config], 4, faces );
        break;

    default: {
        std::stringstream ss;
        ss << "marching_cubes_table::get_cubecase_faces(): Encountered an unknown cubeCorners case.  Corner values as "
              "follows:"
           << std::endl;
        dump_cube( cubeCorners, ss );
        throw std::runtime_error( ss.str() );
    }
    };

    return extraVert;
}

} // namespace volumetrics
} // namespace frantic
