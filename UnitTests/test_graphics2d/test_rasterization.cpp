// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <Eigen/Dense>
#include <frantic/graphics/vector3.hpp>
#include <frantic/graphics2d/rasterization.hpp>
#include <frantic/graphics2d/vector2f.hpp>

using namespace frantic::graphics2d;
using namespace frantic::graphics;

static void rasterize_triangle_onto_buffer_simple( const vector2f vertices[3], const bool boundaryEdges[3],
                                                   Eigen::MatrixXd& outRasterizedImage ) {

    const int bufferWidth = (int)outRasterizedImage.cols();
    const int bufferHeight = (int)outRasterizedImage.rows();
    int startingRow;
    std::vector<frantic::graphics2d::simple_run> runs;
    frantic::graphics2d::rasterize_triangle( vertices[0], vertices[1], vertices[2], bufferWidth, bufferHeight,
                                             boundaryEdges, startingRow, runs );
    for( int rowNum = 0; rowNum < (int)runs.size(); ++rowNum ) {
        const simple_run& run = runs[rowNum];
        const int ri = startingRow + rowNum;
        for( int ci = run.columnRange.first; ci < run.columnRange.second; ++ci ) {
            outRasterizedImage( ri, ci ) = 1;
        }
    }
}

static void rasterize_triangle_onto_buffer_bary( const vector2f vertices[3], const bool boundaryEdges[3],
                                                 Eigen::MatrixXd& outRasterizedImage, Eigen::MatrixXd& outBary0,
                                                 Eigen::MatrixXd& outBary1, Eigen::MatrixXd& outBary2 ) {

    const int bufferWidth = (int)outRasterizedImage.cols();
    const int bufferHeight = (int)outRasterizedImage.rows();
    int startingRow;
    std::vector<frantic::graphics2d::barycentric_run> runs;
    frantic::graphics2d::rasterize_triangle( vertices[0], vertices[1], vertices[2], bufferWidth, bufferHeight,
                                             boundaryEdges, startingRow, runs );
    for( int rowNum = 0; rowNum < (int)runs.size(); ++rowNum ) {
        const barycentric_run& run = runs[rowNum];
        vector3f barycentricCoordinates = run.firstBarycentricCoord;
        const vector3f barycentricCoordinatesDiff = run.deltaBarycentricCoord;
        const int ri = startingRow + rowNum;
        for( int ci = run.columnRange.first; ci < run.columnRange.second; ++ci ) {
            outRasterizedImage( ri, ci ) = 1;
            outBary0( ri, ci ) = barycentricCoordinates[0];
            outBary1( ri, ci ) = barycentricCoordinates[1];
            outBary2( ri, ci ) = barycentricCoordinates[2];
            barycentricCoordinates += barycentricCoordinatesDiff;
        }
    }
}

static Eigen::MatrixXd mat_from_string( const std::string& s ) {
    std::vector<std::string> rows;
    size_t start = 0;
    size_t end = s.find( ' ' );
    size_t numCols = end - start;
    while( end <= s.length() ) {
        rows.push_back( s.substr( start, numCols ) );
        start = end + 1;
        end = start + numCols;
    }
    Eigen::MatrixXd matrix( rows.size(), rows[0].length() );
    for( int r = 0; r < (int)rows.size(); ++r ) {
        for( int c = 0; c < (int)rows[r].length(); ++c ) {
            const char val = rows[r].at( c );
            matrix( r, c ) = val - '0';
        }
    }
    return matrix;
}

struct test_case {
    vector2f vertices[3];
    Eigen::MatrixXd expectedImage;
    bool boundaryEdges[3];
};

TEST( Rasterization, Run ) {
    // For the following test cases, the y coordinates increase from top to bottom
    test_case testCases[] = {
        /*
         * The following set of tests are for triangles that look like the following
         *      X
         *    X X
         */
        // the left edge of the triangle goes through the middle of all diagonal pixels
        { { vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 3.0f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "001 011 111" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 3.0f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "001 011 111" ),
          { false, false, true } },
        // the left edge of the triangle goes slightly below the middle of all diagonal pixels
        { { vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 3.0f ), vector2f( 3.0f, 0.0001f ) },
          mat_from_string( "000 001 011" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 3.0f ), vector2f( 3.0f, 0.0001f ) },
          mat_from_string( "001 011 111" ),
          { false, false, true } },
        // the left edge of the triangle goes through the middle of pixel at row 1, col 2
        { { vector2f( 0.0f, 3.0f ), vector2f( 2.5f, 3.0f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "000 001 011" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 3.0f ), vector2f( 2.5f, 3.0f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "000 011 111" ),
          { false, false, true } },
        // the left edge of the triangle goes slightly below the middle of pixel at row 1, col 2
        { { vector2f( 0.0f, 3.0f ), vector2f( 2.5f, 3.0f ), vector2f( 2.5f, 1.5001f ) },
          mat_from_string( "000 000 011" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 3.0f ), vector2f( 2.5f, 3.0f ), vector2f( 2.5f, 1.5001f ) },
          mat_from_string( "000 011 111" ),
          { false, false, true } },
        // the left edge of the triangle goes through the middle of the pixel at row 2, col 2
        { { vector2f( 0.0f, 3.0f ), vector2f( 2.5f, 3.0f ), vector2f( 2.5f, 2.5000f ) },
          mat_from_string( "000 000 001" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 3.0f ), vector2f( 2.5f, 3.0f ), vector2f( 2.5f, 2.5000f ) },
          mat_from_string( "000 000 111" ),
          { false, false, true } },
        // the left edge of the triangle goes slightly below the middle of the pixel at row 2, col 2
        { { vector2f( 0.0f, 3.0f ), vector2f( 2.5f, 3.0f ), vector2f( 2.5f, 2.5001f ) },
          mat_from_string( "000 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 3.0f ), vector2f( 2.5f, 3.0f ), vector2f( 2.5f, 2.5001f ) },
          mat_from_string( "000 000 111" ),
          { false, false, true } },
        /*
         * The following set of tests are for triangles that look like the following
         *     X
         *     X X
         */
        // the right edge of the triangle goes through the middle of all diagonal pixels
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "100 110 111" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "100 110 111" ),
          { false, false, true } },
        { { vector2f( 3.0f, 3.0f ), vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ) },
          mat_from_string( "100 110 111" ),
          { true, true, true } },
        // the right edge of the triangle goes slightly below the middle of all diagonal pixels
        { { vector2f( 0.0f, 0.0001f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 100 110" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0001f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "100 110 111" ),
          { false, false, true } },
        // the right edge of the triangle goes through the middle of pixel at row 1, col 0
        { { vector2f( 0.5f, 1.5f ), vector2f( 0.5f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 100 110" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 1.5f ), vector2f( 0.5f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 110 111" ),
          { false, false, true } },
        // the right edge of the triangle goes slightly below the middle of pixel at row 1, col 0
        { { vector2f( 0.5f, 1.5001f ), vector2f( 0.5f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 000 110" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 1.5001f ), vector2f( 0.5f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 110 111" ),
          { false, false, true } },
        // the right edge of the triangle goes through the middle of the pixel at row 2, col 0
        { { vector2f( 0.5f, 2.5f ), vector2f( 0.5f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 000 100" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 2.5f ), vector2f( 0.5f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 000 111" ),
          { false, false, true } },
        // the right edge of the triangle goes slightly below the middle of the pixel at row 2, col 0
        { { vector2f( 0.5f, 2.5001f ), vector2f( 0.5f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 2.5001f ), vector2f( 0.5f, 3.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "000 000 111" ),
          { false, false, true } },
        /*
         * The following set of tests are for triangles that look like the following
         *    X X
         *      X
         */
        // the left edge of the triangle goes through the middle of all diagonal pixels
        { { vector2f( 0.0f, 0.0f ), vector2f( 3.0f, 0.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "111 011 001" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 3.0f, 0.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "111 011 001" ),
          { false, false, true } },
        // same as above but switch the order of vertex 0 and 1
        { { vector2f( 3.0f, 0.0f ), vector2f( 0.0f, 0.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "111 011 001" ),
          { false, false, false } },
        //
        { { vector2f( 3.0f, 0.0f ), vector2f( 0.0f, 0.0f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "111 011 001" ),
          { false, false, true } },
        // the left edge of the triangle goes slightly below the middle of all diagonal pixels
        { { vector2f( 0.0f, 0.0f ), vector2f( 3.0f, 0.0f ), vector2f( 3.0f, 2.9999f ) },
          mat_from_string( "011 001 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 3.0f, 0.0f ), vector2f( 3.0f, 2.9999f ) },
          mat_from_string( "111 011 001" ),
          { false, false, true } },
        // same as above but switch the order of vertex 0 and 1
        { { vector2f( 3.0f, 0.0f ), vector2f( 0.0f, 0.0f ), vector2f( 3.0f, 2.9999f ) },
          mat_from_string( "011 001 000" ),
          { false, false, false } },
        // the left edge of the triangle goes through the middle of pixel at row 1, col 2
        { { vector2f( 0.0f, 0.0f ), vector2f( 2.5f, 0.0f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "011 001 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 2.5f, 0.0f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "111 011 000" ),
          { false, false, true } },
        // the left edge of the triangle goes slightly below the middle of pixel at row 1, col 2
        { { vector2f( 0.0f, 0.0f ), vector2f( 2.5f, 0.0f ), vector2f( 2.5f, 1.4999f ) },
          mat_from_string( "011 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 2.5f, 0.0f ), vector2f( 2.5f, 1.4999f ) },
          mat_from_string( "111 011 000" ),
          { false, false, true } },
        // the left edge of the triangle goes through the middle of the pixel at row 0, col 2
        { { vector2f( 0.0f, 0.0f ), vector2f( 2.5f, 0.0f ), vector2f( 2.5f, 0.5f ) },
          mat_from_string( "001 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 2.5f, 0.0f ), vector2f( 2.5f, 0.5f ) },
          mat_from_string( "111 000 000" ),
          { false, false, true } },
        // the left edge of the triangle goes slightly below the middle of the pixel at row 0, col 2
        { { vector2f( 0.0f, 0.0f ), vector2f( 2.5f, 0.0f ), vector2f( 2.5f, 0.4999f ) },
          mat_from_string( "000 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 2.5f, 0.0f ), vector2f( 2.5f, 0.4999f ) },
          mat_from_string( "111 000 000" ),
          { false, false, true } },
        /*
         * The following set of tests are for triangles that look like the following
         *     X X
         *     X
         */
        // the right edge of the triangle goes through the middle of all diagonal pixels
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "111 110 100" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "111 110 100" ),
          { false, true, false } },
        // the right edge of the triangle goes slightly below the middle of all diagonal pixels
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 2.9999f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "110 100 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 2.9999f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "111 110 100" ),
          { false, true, false } },
        // the right edge of the triangle goes through the middle of pixel at row 1, col 0
        { { vector2f( 0.5f, 0.0f ), vector2f( 0.5f, 1.5f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "110 100 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.0f ), vector2f( 0.5f, 1.5f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "111 110 000" ),
          { false, true, false } },
        // the right edge of the triangle goes slightly below the middle of pixel at row 1, col 0
        { { vector2f( 0.5f, 0.0f ), vector2f( 0.5f, 1.4999f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "110 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.0f ), vector2f( 0.5f, 1.4999f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "111 110 000" ),
          { false, true, false } },
        // the right edge of the triangle goes through the middle of the pixel at row 0, col 0
        { { vector2f( 0.5f, 0.0f ), vector2f( 0.5f, 0.5f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "100 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.0f ), vector2f( 0.5f, 0.5f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "111 000 000" ),
          { false, true, false } },
        // same as above but switch the order of vertex 0 and 1
        { { vector2f( 0.5f, 0.5f ), vector2f( 0.5f, 0.0f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "100 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 0.5f, 0.0f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "111 000 000" ),
          { false, false, true } },
        // the right edge of the triangle goes slightly below the middle of the pixel at row 0, col 0
        { { vector2f( 0.5f, 0.0f ), vector2f( 0.5f, 0.4999f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "000 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.0f ), vector2f( 0.5f, 0.4999f ), vector2f( 3.0f, 0.0f ) },
          mat_from_string( "111 000 000" ),
          { false, true, false } },
        /*
         * The following set of tests are for triangles that look like the following
         *     X X X
         *       X
         */
        // all the vertices are in the middle of a pixel
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "111 010 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "111 010 000" ),
          { false, false, true } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "111 010 000" ),
          { false, true, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "111 010 000" ),
          { true, false, false } },
        // lowest vertex is slightly right from the middle of pixel of the pixel, others are in the middle
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5001f, 1.5f ) },
          mat_from_string( "111 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5001f, 1.5f ) },
          mat_from_string( "111 001 000" ),
          { false, true, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5001f, 1.5f ) },
          mat_from_string( "111 010 000" ),
          { false, false, true } },
        // lowest vertex is slightly left from the middle of pixel of the pixel, others are in the middle
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.4999f, 1.5f ) },
          mat_from_string( "111 000 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.4999f, 1.5f ) },
          mat_from_string( "111 010 000" ),
          { false, true, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.4999f, 1.5f ) },
          mat_from_string( "111 100 000" ),
          { false, false, true } },
        // left vertex is slightly right from the middle of pixel of the pixel, others are in the middle
        { { vector2f( 0.5001f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "011 010 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5001f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "011 010 000" ),
          { false, true, false } },
        //
        { { vector2f( 0.5001f, 0.5f ), vector2f( 2.5f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "111 010 000" ),
          { false, false, true } },
        // right vertex is slightly left from the middle of pixel of the pixel, others are in the middle
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.4999f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "110 010 000" ),
          { false, false, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.4999f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "111 010 000" ),
          { false, true, false } },
        //
        { { vector2f( 0.5f, 0.5f ), vector2f( 2.4999f, 0.5f ), vector2f( 1.5f, 1.5f ) },
          mat_from_string( "110 010 000" ),
          { false, false, true } },
        /*
         * The following set of tests are for triangles that look like the following
         *        X
         *      X X X
         */
        // all the vertices are in the middle of a pixel
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 111 000" ),
          { false, false, false } },
        //
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 111 000" ),
          { false, false, true } },
        //
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 111 000" ),
          { false, true, false } },
        //
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 111 000" ),
          { true, false, false } },
        // highest vertex is slightly right from the middle of pixel of the pixel, others are in the middle
        { { vector2f( 1.5001f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "000 111 000" ),
          { false, false, false } },
        //
        { { vector2f( 1.5001f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 111 000" ),
          { true, false, false } },
        //
        { { vector2f( 1.5001f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "001 111 000" ),
          { false, false, true } },
        //
        { { vector2f( 1.5001f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "011 111 000" ),
          { true, false, true } },
        // highest vertex is slightly left from the middle of pixel of the pixel, others are in the middle
        { { vector2f( 1.4999f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "000 111 000" ),
          { false, false, false } },
        //
        { { vector2f( 1.4999f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "100 111 000" ),
          { true, false, false } },
        //
        { { vector2f( 1.4999f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 111 000" ),
          { false, false, true } },
        //
        { { vector2f( 1.4999f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "110 111 000" ),
          { true, false, true } },
        // left vertex is slightly right from the middle of pixel of the pixel, others are in the middle
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5001f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 011 000" ),
          { false, false, false } },
        //
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5001f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 111 000" ),
          { true, false, false } },
        //
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5001f, 1.5f ), vector2f( 2.5f, 1.5f ) },
          mat_from_string( "010 011 000" ),
          { false, false, true } },
        // right vertex is slightly left from the middle of pixel of the pixel, others are in the middle
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.4999f, 1.5f ) },
          mat_from_string( "010 110 000" ),
          { false, false, false } },
        //
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.4999f, 1.5f ) },
          mat_from_string( "010 110 000" ),
          { true, false, false } },
        //
        { { vector2f( 1.5f, 0.5f ), vector2f( 0.5f, 1.5f ), vector2f( 2.4999f, 1.5f ) },
          mat_from_string( "010 111 000" ),
          { false, false, true } },
        /*
         * The following set of tests are for triangles that look like the following
         *        X
         *      X X
         *        X
         */
        { { vector2f( 3.0f, 0.0f ), vector2f( 0.0f, 1.5f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "001 111 001" ),
          { false, false, false } },
        //
        { { vector2f( 3.0f, 0.0f ), vector2f( 0.0f, 1.5f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "011 111 001" ),
          { true, false, false } },
        //
        { { vector2f( 3.0f, 0.0f ), vector2f( 0.0f, 1.5f ), vector2f( 3.0f, 3.0f ) },
          mat_from_string( "001 111 011" ),
          { false, true, false } },
        // with both the top and bottom lines almost horizontal
        { { vector2f( 4.5f, 0.0f ), vector2f( 4.5f, 0.3f ), vector2f( 1.5f, 0.3f ) },
          mat_from_string( "011110 000000" ),
          { true, true, true } },
        /*
         * The following set of tests are for triangles that look like the following
         *      X
         *      X X
         *      X
         */
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 1.5f ) },
          mat_from_string( "100 111 100" ),
          { false, false, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 1.5f ) },
          mat_from_string( "100 111 110" ),
          { false, true, false } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 1.5f ) },
          mat_from_string( "110 111 100" ),
          { false, false, true } },
        //
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 1.5f ) },
          mat_from_string( "110 111 110" ),
          { true, true, true } },
        // with the top line almost horizontal
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 3.0f ), vector2f( 3.0f, 0.001f ) },
          mat_from_string( "111 111 110" ),
          { true, true, true } },
        // with both the top and bottom lines almost horizontal
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 0.001f ), vector2f( 3.0f, 0.001f ) },
          mat_from_string( "111 000 000" ),
          { true, true, true } },
        // with both the top and bottom lines almost horizontal
        { { vector2f( 3.0f, 1.0f ), vector2f( 0.5f, 0.5f ), vector2f( 0.0f, 0.0f ) },
          mat_from_string( "111 000 000" ),
          { true, true, true } },
        // with both the top and bottom lines almost horizontal
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 0.3f ), vector2f( 4.5f, 0.3f ) },
          mat_from_string( "111110 000000" ),
          { true, true, true } },
        //
        { { vector2f( 1.5f, 0.0f ), vector2f( 1.5f, 0.2f ), vector2f( 4.5f, 0.2f ) },
          mat_from_string( "011110 000000" ),
          { true, true, true } },
        // with both the top and bottom lines almost horizontal
        { { vector2f( 0.0f, 0.0f ), vector2f( 0.0f, 0.001f ), vector2f( 3.0f, 0.0001f ) },
          mat_from_string( "000000 000000" ),
          { false, false, false } } };

    for( size_t i = 0; i < sizeof( testCases ) / sizeof( testCases[0] ); ++i ) {
        const test_case& testCase = testCases[i];
        const int bufferWidth = (int)testCase.expectedImage.cols();
        const int bufferHeight = (int)testCase.expectedImage.rows();

        // Simple run
        Eigen::MatrixXd rasterizedImageSimple = Eigen::MatrixXd::Zero( bufferHeight, bufferWidth );
        rasterize_triangle_onto_buffer_simple( testCase.vertices, testCase.boundaryEdges, rasterizedImageSimple );

        std::stringstream failBitmapSimpleMsg;
        failBitmapSimpleMsg << "testCase# " << i << " failed" << std::endl;
        failBitmapSimpleMsg << "rasterized image is " << std::endl << rasterizedImageSimple << std::endl;
        failBitmapSimpleMsg << "expected image is " << std::endl << testCase.expectedImage << std::endl;
        ASSERT_EQ( ( rasterizedImageSimple - testCase.expectedImage ).norm(), 0.0 ) << failBitmapSimpleMsg.str();

        // Barycentric run
        Eigen::MatrixXd rasterizedImageBary = Eigen::MatrixXd::Zero( bufferHeight, bufferWidth );
        Eigen::MatrixXd baryCoord0 = Eigen::MatrixXd::Zero( bufferHeight, bufferWidth );
        Eigen::MatrixXd baryCoord1 = Eigen::MatrixXd::Zero( bufferHeight, bufferWidth );
        Eigen::MatrixXd baryCoord2 = Eigen::MatrixXd::Zero( bufferHeight, bufferWidth );
        rasterize_triangle_onto_buffer_bary( testCase.vertices, testCase.boundaryEdges, rasterizedImageBary, baryCoord0,
                                             baryCoord1, baryCoord2 );

        std::stringstream failBitmapBaryMsg;
        failBitmapBaryMsg << "test# " << i << " failed" << std::endl;
        failBitmapBaryMsg << "rasterized image is " << std::endl << rasterizedImageBary << std::endl;
        failBitmapBaryMsg << "expected image is " << std::endl << testCase.expectedImage << std::endl;
        ASSERT_EQ( ( rasterizedImageBary - testCase.expectedImage ).norm(), 0.0 ) << failBitmapBaryMsg.str();

        // Check that the barycentric coordinates make sense
        for( int r = 0; r < bufferHeight; ++r ) {
            for( int c = 0; c < bufferWidth; ++c ) {
                if( testCase.expectedImage( r, c ) != 0 ) {
                    // We can get the cartesian coordinates for a pixel from the barycentric coordinates by a linear
                    // combination of the barycentric coordinates with the corresponding cartesian coordinates of the
                    // vertices of the triangle
                    vector2f coord = (float)baryCoord0( r, c ) * testCase.vertices[0] +
                                     (float)baryCoord1( r, c ) * testCase.vertices[1] +
                                     (float)baryCoord2( r, c ) * testCase.vertices[2];

                    // If the following 3 assertions are true, baryCoord satisfies the definition
                    // of barycentric coordinates
                    std::stringstream failBarycoordMsg;
                    failBarycoordMsg << "Rasterized barycentric coordinates for test# " << i << " at pixel ( " << c
                                     << ", " << r << " )" << std::endl;
                    failBarycoordMsg << "which are ( " << baryCoord0( r, c ) << ", " << baryCoord1( r, c ) << ", "
                                     << baryCoord2( r, c ) << " ) ";
                    failBarycoordMsg << "do not satisfy the definition of barycentric coordinates " << std::endl;

                    // The cartesian coordinates for the pixel is expected to be located at the center of the pixel.
                    vector2f expectedCoord( c + 0.5f, r + 0.5f );
                    EXPECT_LT( ( coord - expectedCoord ).get_magnitude(), 0.001 )
                        << failBarycoordMsg.str() << "Expected cartesian coord " << expectedCoord << ". Got " << coord
                        << std::endl;

                    // The barycentric coordinates must sum up to 1
                    double baryCoordSum = baryCoord0( r, c ) + baryCoord1( r, c ) + baryCoord2( r, c );
                    EXPECT_NEAR( baryCoordSum, 1.0, 0.001 )
                        << failBarycoordMsg.str() << "Expected sum of barycentric coordinates 1.0. Got " << baryCoordSum
                        << std::endl;
                }
            }
        }
    }
}
