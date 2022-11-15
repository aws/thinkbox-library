// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( E57_AVAILABLE )

#include "utilities/e57_generator.hpp"
#include <E57Format/E57Format.h>

#include <frantic/channels/channel_accessor.hpp>
#include <frantic/channels/channel_map.hpp>
#include <frantic/particles/particle_file_metadata.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/prt_metadata.hpp>
#include <frantic/particles/streams/prt2_particle_ostream.hpp>

#include <cmath>
#include <frantic/graphics/vector3.hpp>
#include <string>

using namespace std;
using namespace e57;
using namespace frantic::channels;
using namespace frantic::strings;
using namespace frantic::particles::streams;
using namespace frantic::particles;
using namespace frantic::graphics;

void generate_files( const string& fileName, bool hasIntensity, bool hasColor, bool hasNormal, bool posIsDouble,
                     int numScans, bool posIsScaledInt, bool hasTransform, bool colorIsInt ) {

    try {
        string outfile = fileName + ".e57";
        /// Open new file for writing, get the initialized root node (a Structure).
        ImageFile imf( outfile, "w" );
        StructureNode root = imf.root();

        // TODO: Modify these variables when producing a new file ===================
        /*
            // Boolean variables to determine which fields should be included in the file.
            bool hasIntensity = true;
            bool hasColor = true;
            bool hasNormal = false;
            int numScans = 1;
            //bool hasExtents = true;

            // Determine if position is double or single precision
            bool posIsDouble = false;
        */
        //============================================================================

        /// Register extension with URI=www.example.com/DemoExtension and prefix=demo
        imf.extensionsAdd( "nor", "http://www.libe57.org/E57_NOR_surface_normals.txt" );

        /// Set per-file properties.
        /// Path names: "/formatName", "/majorVersion", "/minorVersion", "/coordinateMetadata"
        root.set( "formatName", StringNode( imf, "ASTM E57 3D Imaging Data File" ) );
        root.set( "guid", StringNode( imf, "3F2504E0-4F89-11D3-9A0C-0305E82C3300" ) );

        /// Get ASTM version number supported by library, so can write it into file
        int astmMajor, astmMinor;
        ustring libraryId;
        e57::Utilities::getVersions( astmMajor, astmMinor, libraryId );
        root.set( "versionMajor", IntegerNode( imf, astmMajor ) );
        root.set( "versionMinor", IntegerNode( imf, astmMinor ) );

        /// Save a dummy string for coordinate system.
        /// Really should be a valid WKT string identifying the coordinate reference system (CRS).
        root.set( "coordinateMetadata", StringNode( imf, "...A WKT string here..." ) );

        /// Create creationDateTime structure
        /// Path name: "/creationDateTime
        StructureNode creationDateTime = StructureNode( imf );
        root.set( "creationDateTime", creationDateTime );
        creationDateTime.set( "dateTimeValue", FloatNode( imf, 123.456 ) ); //!!! convert time() to GPStime

        /// Create 3D data area.
        /// Path name: "/data3D"
        VectorNode data3D = VectorNode( imf, true );
        root.set( "data3D", data3D );
        vector<StructureNode> scans;
        vector<CompressedVectorNode> groupNodes; // one for each scan
        vector<CompressedVectorNode> pointNodes; // one for each scan

        /// Add first scan
        /// Path name: "/data3D/0"
        for( int i = 0; i < numScans; ++i ) {
            StructureNode scan = StructureNode( imf );
            data3D.append( scan );
            /// Add guid to scan
            /// Path name: "/data3D/0/guid".
            static const char* guids[] = {
                "3F2504E0-4F89-11D3-9A0C-0305E82C3301", "04FDBC32-C52A-491D-8248-D8DE54469E12",
                "E47BBBF3-C4A5-43E4-BC59-963545A6AD8A", "F27EE8D1-CA66-4BD5-8D19-667AF139DA89",
                "499D194A-B013-4F6D-8950-34BE0A69D7E7", "8C151AB6-D4C1-43F3-BE6B-70AA8674938D",
                "36DE968F-C21D-43AB-8093-83ABFA7F84A8", "C6731D21-B59D-4FD2-B2AF-696FFEE5E04D",
                "AB3EF4D4-9082-4FA1-8C44-30F6C0EC4727", "E7538FA9-B939-4D35-83D1-B5E7B9621731" };
            scan.set( "guid", StringNode( imf, guids[i % ( sizeof( guids ) / sizeof( guids[0] ) )] ) );
            scans.push_back( scan );
        }

        /// Make a prototype of datatypes that will be stored in points record.
        /// This prototype will be used in creating the points CompressedVector.
        /// Using this proto in a CompressedVector will define path names like:
        ///      "/data3D/0/points/0/cartesianX"
        StructureNode proto = StructureNode( imf );

        if( posIsScaledInt ) {
            proto.set( "cartesianX", ScaledIntegerNode( imf, 0, 0, 32767, 0.001 ) );
            proto.set( "cartesianY", ScaledIntegerNode( imf, 0, 0, 32767, 0.001 ) );
            proto.set( "cartesianZ", ScaledIntegerNode( imf, 0, 0, 32767, 0.001 ) );
        } else {
            proto.set( "cartesianX", FloatNode( imf ) );
            proto.set( "cartesianY", FloatNode( imf ) );
            proto.set( "cartesianZ", FloatNode( imf ) );
        }
        proto.set( "cartesianInvalidState", IntegerNode( imf, 0, 0, 2 ) );
        proto.set( "rowIndex", IntegerNode( imf, 0, 0, 1 ) );
        proto.set( "columnIndex", IntegerNode( imf, 0, 0, 4 ) );
        proto.set( "returnIndex", IntegerNode( imf, 0, 0, 0 ) );
        proto.set( "returnCount", IntegerNode( imf, 1, 1, 1 ) );
        proto.set( "timeStamp", FloatNode( imf, 0.0, E57_DOUBLE ) );
        proto.set( "intensity", IntegerNode( imf, 0, 0, 255 ) );
        if( !colorIsInt ) {
            proto.set( "colorRed", FloatNode( imf, 0.0, E57_SINGLE, 0.0, 1.0 ) );
            proto.set( "colorGreen", FloatNode( imf, 0.0, E57_SINGLE, 0.0, 1.0 ) );
            proto.set( "colorBlue", FloatNode( imf, 0.0, E57_SINGLE, 0.0, 1.0 ) );
        } else {
            proto.set( "colorRed", IntegerNode( imf, 0, 0, 255 ) );
            proto.set( "colorGreen", IntegerNode( imf, 0, 0, 255 ) );
            proto.set( "colorBlue", IntegerNode( imf, 0, 0, 255 ) );
        }
        /// Add normal fields
        if( hasNormal ) {
            proto.set( "nor:normalX", FloatNode( imf, 0.0, E57_SINGLE ) );
            proto.set( "nor:normalY", FloatNode( imf, 0.0, E57_SINGLE ) );
            proto.set( "nor:normalZ", FloatNode( imf, 0.0, E57_SINGLE ) );
        }

        /// Make empty codecs vector for use in creating points CompressedVector.
        /// If this vector is empty, it is assumed that all fields will use the BitPack codec.
        VectorNode codecs = VectorNode( imf, true );

        /// Create CompressedVector for storing points.  Path Name: "/data3D/0/points".
        /// We use the prototype and empty codecs tree from above.
        /// The CompressedVector will be filled by code below.
        for( int i = 0; i < scans.size(); ++i ) {
            StructureNode& scan = scans[i];
            CompressedVectorNode points = CompressedVectorNode( imf, proto, codecs );
            scan.set( "points", points );
            pointNodes.push_back( points );
        }

        const quat4fd rotationTransform = quat4fd();
        const vector3fd translationTransform = vector3fd();
        /// translation increased by translationTransformStep for each subsequent scan
        const vector3fd translationTransformStep = ( hasTransform ? vector3fd( 1.0, 0.2, 3.0 ) : vector3fd() );

        for( int i = 0; i < numScans; ++i ) {
            StructureNode& scan = scans[i];
            /// Create pose structure for scan.
            /// Path names: "/data3D/0/pose/rotation/w", etc...
            ///             "/data3D/0/pose/translation/x", etc...
            StructureNode pose = StructureNode( imf );
            scan.set( "pose", pose );
            StructureNode rotation = StructureNode( imf );
            pose.set( "rotation", rotation );
            rotation.set( "w", FloatNode( imf, rotationTransform.w ) );
            rotation.set( "x", FloatNode( imf, rotationTransform.x ) );
            rotation.set( "y", FloatNode( imf, rotationTransform.y ) );
            rotation.set( "z", FloatNode( imf, rotationTransform.z ) );
            StructureNode translation = StructureNode( imf );
            pose.set( "translation", translation );
            translation.set( "x", FloatNode( imf, ( translationTransform + i * translationTransformStep ).x ) );
            translation.set( "y", FloatNode( imf, ( translationTransform + i * translationTransformStep ).y ) );
            translation.set( "z", FloatNode( imf, ( translationTransform + i * translationTransformStep ).z ) );

            /// Add grouping scheme area
            /// Path name: "/data3D/0/pointGroupingSchemes"
            StructureNode pointGroupingSchemes = StructureNode( imf );
            scan.set( "pointGroupingSchemes", pointGroupingSchemes );

            /// Add a line grouping scheme
            /// Path name: "/data3D/0/pointGroupingSchemes/groupingByLine"
            StructureNode groupingByLine = StructureNode( imf );
            pointGroupingSchemes.set( "groupingByLine", groupingByLine );

            /// Add idElementName to groupingByLine, specify a line is column oriented
            /// Path name: "/data3D/0/pointGroupingSchemes/groupingByLine/idElementName"
            groupingByLine.set( "idElementName", StringNode( imf, "columnIndex" ) );

            /// Make a prototype of datatypes that will be stored in LineGroupRecord.
            /// This prototype will be used in creating the groups CompressedVector.
            /// Will define path names like:
            ///     "/data3D/0/pointGroupingSchemes/groupingByLine/groups/0/idElementValue"
            StructureNode lineGroupProto = StructureNode( imf );
            lineGroupProto.set( "idElementValue", IntegerNode( imf, 0, 0, 4 ) );
            lineGroupProto.set( "startPointIndex", IntegerNode( imf, 0, 0, 9 ) );
            lineGroupProto.set( "pointCount", IntegerNode( imf, 1, 1, 2 ) );

            /// Add cartesian bounds to line group prototype
            /// Will define path names like:
            ///     "/data3D/0/pointGroupingSchemes/groupingByLine/groups/0/cartesianBounds/xMinimum"
            StructureNode lineGroupBbox = StructureNode( imf );
            lineGroupProto.set( "cartesianBounds", lineGroupBbox );
            lineGroupBbox.set( "xMinimum", FloatNode( imf, 0.0 ) );
            lineGroupBbox.set( "xMaximum", FloatNode( imf, 0.0 ) );
            lineGroupBbox.set( "yMinimum", FloatNode( imf, 0.0 ) );
            lineGroupBbox.set( "yMaximum", FloatNode( imf, 0.0 ) );
            lineGroupBbox.set( "zMinimum", FloatNode( imf, 0.0 ) );
            lineGroupBbox.set( "zMaximum", FloatNode( imf, 0.0 ) );

            /// Make empty codecs vector for use in creating groups CompressedVector.
            /// If this vector is empty, it is assumed that all fields will use the BitPack codec.
            VectorNode lineGroupCodecs = VectorNode( imf, true );

            /// Create CompressedVector for storing groups.
            /// Path Name: "/data3D/0/pointGroupingSchemes/groupingByLine/groups".
            /// We use the prototype and empty codecs tree from above.
            /// The CompressedVector will be filled by code below.
            CompressedVectorNode groups = CompressedVectorNode( imf, lineGroupProto, lineGroupCodecs );
            groupingByLine.set( "groups", groups );
            groupNodes.push_back( groups );

            /// Add name and description to scan
            /// Path names: "/data3D/0/name", "/data3D/0/description".
            scan.set( "name", StringNode( imf, "Point Cloud" ) );
            scan.set( "description",
                      StringNode( imf, "Point cloud generated by FranticLibrary for testing purposes." ) );
        }

        // Create bound box in local coordinate system for each scan. BoundBox should be transformed in the istream
        // along with point data.
        boundbox3fd boundBox;
        if( posIsDouble )
            boundBox = boundbox3fd(
                vector3fd( 100000000000000.01110, 1.1 * 100000000000000.01110, 1.2 * 100000000000000.01110 ),
                vector3fd( 10.0 * 100000000000000.01110, 10.1 * 100000000000000.01110, 10.2 * 100000000000000.01110 ) );
        else
            boundBox = boundbox3fd( vector3fd( 1.0, 1.1, 1.2 ), vector3fd( 10.0, 10.1, 10.2 ) );

        boundbox3fd boundBox1 = boundBox;

        /// Add Cartesian bounding box to scan.
        /// Path names: "/data3D/0/cartesianBounds/xMinimum", etc...
        for( int i = 0; i < scans.size(); ++i ) {
            StructureNode& scan = scans[i];
            StructureNode bbox = StructureNode( imf );
            bbox.set( "xMinimum", FloatNode( imf, boundBox.xminimum() ) );
            bbox.set( "xMaximum", FloatNode( imf, boundBox.xmaximum() ) );
            bbox.set( "yMinimum", FloatNode( imf, boundBox.yminimum() ) );
            bbox.set( "yMaximum", FloatNode( imf, boundBox.ymaximum() ) );
            bbox.set( "zMinimum", FloatNode( imf, boundBox.zminimum() ) );
            bbox.set( "zMaximum", FloatNode( imf, boundBox.zmaximum() ) );
            scan.set( "cartesianBounds", bbox );

            /// Add start/stop acquisition times to scan.
            /// Path names: "/data3D/0/acquisitionStart/dateTimeValue",
            ///             "/data3D/0/acquisitionEnd/dateTimeValue"
            StructureNode acquisitionStart = StructureNode( imf );
            scan.set( "acquisitionStart", acquisitionStart );
            acquisitionStart.set( "dateTimeValue", FloatNode( imf, 1235. ) );
            StructureNode acquisitionEnd = StructureNode( imf );
            scan.set( "acquisitionEnd", acquisitionEnd );
            acquisitionEnd.set( "dateTimeValue", FloatNode( imf, 1235. ) );
        }

        /// ---------------------------------------------------------------------------------------
        /// Prepare vector of source buffers for writing in the CompressedVector of points
        const int N = 10;
        double cartesianX[N] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
        double cartesianY[N] = { 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1 };
        double cartesianZ[N] = { 1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2, 9.1, 10.2 };
        if( posIsDouble ) {
            for( int i = 0; i < N; i++ ) {
                cartesianX[i] = 100000000000000.01110 * cartesianX[i];
                cartesianY[i] = 100000000000000.01110 * cartesianY[i];
                cartesianZ[i] = 100000000000000.01110 * cartesianZ[i];
            }
        }
        int32_t cartesianInvalidState[N] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        int32_t rowIndex[N] = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };
        int32_t columnIndex[N] = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4 };
        int32_t returnIndex[N] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        int32_t returnCount[N] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        double timeStamp[N] = { .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0 };
        int32_t intensity[N] = { 1, 2, 3, 2, 1, 1, 2, 3, 2, 1 };
        float colorRed[N] = { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f };
        float colorGreen[N] = { 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f };
        float colorBlue[N] = { 0.9f, 0.8f, 0.7f, 0.6f, 0.5f, 0.4f, 0.3f, 0.2f, 0.1f, 0.0f };
        int intColorRed[N] = { 0, 1, 2, 4, 8, 16, 32, 64, 255, 255 };
        int intColorGreen[N] = { 128, 128, 128, 128, 128, 128, 128, 128, 255, 128 };
        int intColorBlue[N] = { 255, 127, 63, 31, 15, 7, 3, 1, 255, 155 };
        float normalX[N] = { 1.0f, 0.0f, 0.0f, 0.6f, 0.0f, 0.8f, (float)std::sqrt( 1 / 3.0 ), 0.8f, 0.0f, 0.6f };
        float normalY[N] = { 0.0f, 1.0f, 0.0f, 0.8f, 0.6f, 0.0f, (float)std::sqrt( 1 / 3.0 ), 0.6f, 0.8f, 0.0f };
        float normalZ[N] = { 0.0f, 0.0f, 1.0f, 0.0f, 0.8f, 0.6f, (float)std::sqrt( 1 / 3.0 ), 0.0f, 0.6f, 0.8f };
        vector<SourceDestBuffer> sourceBuffers;
        sourceBuffers.push_back( SourceDestBuffer( imf, "cartesianX", cartesianX, N, true, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "cartesianY", cartesianY, N, true, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "cartesianZ", cartesianZ, N, true, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "cartesianInvalidState", cartesianInvalidState, N, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "rowIndex", rowIndex, N, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "columnIndex", columnIndex, N, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "returnIndex", returnIndex, N, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "returnCount", returnCount, N, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "timeStamp", timeStamp, N, true ) );
        sourceBuffers.push_back( SourceDestBuffer( imf, "intensity", intensity, N, true ) );
        if( !colorIsInt ) {
            sourceBuffers.push_back( SourceDestBuffer( imf, "colorRed", colorRed, N, true ) );
            sourceBuffers.push_back( SourceDestBuffer( imf, "colorGreen", colorGreen, N, true ) );
            sourceBuffers.push_back( SourceDestBuffer( imf, "colorBlue", colorBlue, N, true ) );
        } else {
            sourceBuffers.push_back( SourceDestBuffer( imf, "colorRed", intColorRed, N, true ) );
            sourceBuffers.push_back( SourceDestBuffer( imf, "colorGreen", intColorGreen, N, true ) );
            sourceBuffers.push_back( SourceDestBuffer( imf, "colorBlue", intColorBlue, N, true ) );
        }
        if( hasNormal ) {
            sourceBuffers.push_back( SourceDestBuffer( imf, "nor:normalX", normalX, N, true ) );
            sourceBuffers.push_back( SourceDestBuffer( imf, "nor:normalY", normalY, N, true ) );
            sourceBuffers.push_back( SourceDestBuffer( imf, "nor:normalZ", normalZ, N, true ) );
        }

        /// Write source buffers into CompressedVector
        for( int i = 0; i < pointNodes.size(); ++i ) {
            CompressedVectorWriter writer = pointNodes[i].writer( sourceBuffers );
            writer.write( N );
            writer.close();
        }

        /// Prepare vector of source buffers for writing in the CompressedVector of groups
        const int NG = 5;
        int32_t idElementValue[NG] = { 0, 1, 2, 3, 4 };
        int32_t startPointIndex[NG] = { 0, 2, 4, 6, 8 };
        int32_t pointCount[NG] = { 2, 2, 2, 2, 2 };
        double xMinimum[NG] = { 1.0, 3.0, 5.0, 7.0, 9.0 };
        double xMaximum[NG] = { 2.0, 4.0, 6.0, 8.0, 10.0 };
        double yMinimum[NG] = { 1.1, 3.1, 5.1, 7.1, 9.1 };
        double yMaximum[NG] = { 2.1, 4.1, 6.1, 8.1, 10.1 };
        double zMinimum[NG] = { 1.2, 3.2, 5.2, 7.2, 9.2 };
        double zMaximum[NG] = { 2.2, 4.2, 6.2, 8.2, 10.2 };
        vector<SourceDestBuffer> groupSDBuffers;
        groupSDBuffers.push_back( SourceDestBuffer( imf, "idElementValue", idElementValue, NG, true ) );
        groupSDBuffers.push_back( SourceDestBuffer( imf, "startPointIndex", startPointIndex, NG, true ) );
        groupSDBuffers.push_back( SourceDestBuffer( imf, "pointCount", pointCount, NG, true ) );
        groupSDBuffers.push_back( SourceDestBuffer( imf, "cartesianBounds/xMinimum", xMinimum, NG, true ) );
        groupSDBuffers.push_back( SourceDestBuffer( imf, "cartesianBounds/xMaximum", xMaximum, NG, true ) );
        groupSDBuffers.push_back( SourceDestBuffer( imf, "cartesianBounds/yMinimum", yMinimum, NG, true ) );
        groupSDBuffers.push_back( SourceDestBuffer( imf, "cartesianBounds/yMaximum", yMaximum, NG, true ) );
        groupSDBuffers.push_back( SourceDestBuffer( imf, "cartesianBounds/zMinimum", zMinimum, NG, true ) );
        groupSDBuffers.push_back( SourceDestBuffer( imf, "cartesianBounds/zMaximum", zMaximum, NG, true ) );

        /// Write source buffers into CompressedVector
        for( int i = 0; i < groupNodes.size(); ++i ) {
            CompressedVectorWriter writer = groupNodes[i].writer( groupSDBuffers );
            writer.write( NG );
            writer.close();
        }

        imf.close();

        //=====================================================================
        // Write corresponding prt file

        // Get filename
        string prtFileName = fileName + ".prt";

        // Set metadata
        particle_file_metadata metadata;
        property_map& generalMetadata = metadata.get_general_metadata();
        prt::length_unit_in_micrometers::add_channel( generalMetadata );
        prt::length_unit_in_micrometers::set_value( generalMetadata, 1.0e6 );

        // Set scanner transforms
        vector<transform4fd> scannerTransforms;
        if( numScans > 1 ) {
            for( int i = 0; i < numScans; ++i ) {
                transform4fd transform;
                rotationTransform.as_transform4f( transform );
                transform.set_translation( translationTransform + i * translationTransformStep );
                scannerTransforms.push_back( transform );
            }
            prt::set_scanner_transforms( generalMetadata, scannerTransforms );
        }

        property_map positionMeta;
        prt::add_channel_extents( positionMeta, data_type_float64 );
        prt::set_extents( positionMeta, boundBox );
        metadata.append_channel_metadata( _T( "Position" ), positionMeta );

        // Generate channel_map
        channel_map cm = channel_map();
        if( posIsDouble ) {
            cm.define_channel<vector3fd>( _T("Position") );
        } else {
            cm.define_channel<vector3f>( _T("Position") );
        }
        if( hasIntensity )
            cm.define_channel<float>( _T("Intensity") );
        if( hasColor )
            cm.define_channel<vector3f>( _T("Color") );
        if( hasNormal )
            cm.define_channel<vector3f>( _T("Normal") );
        if( numScans > 1 )
            cm.define_channel<uint32_t>( _T("ScannerIndex") );
        cm.end_channel_definition();

        // Create prt ostream
        frantic::particles::particle_file_stream_factory_object factory;
        factory.enable_prt2_saving();
        boost::shared_ptr<particle_ostream> ostream =
            factory.create_ostream( to_tstring( prtFileName ), cm, cm, &metadata );

        // Create channel accessors for writing particles
        channel_cvt_accessor<vector3fd> posAccessor;
        channel_accessor<float> intensityAccessor;
        channel_accessor<vector3f> colorAccessor;
        channel_accessor<vector3f> normalAccessor;
        channel_accessor<uint32_t> scanIndexAccessor;
        posAccessor = cm.get_cvt_accessor<vector3fd>( _T("Position") );
        if( hasIntensity )
            intensityAccessor = cm.get_accessor<float>( _T("Intensity") );
        if( hasColor )
            colorAccessor = cm.get_accessor<vector3f>( _T("Color") );
        if( hasNormal )
            normalAccessor = cm.get_accessor<vector3f>( _T("Normal") );
        if( numScans > 1 )
            scanIndexAccessor = cm.get_accessor<uint32_t>( _T("ScannerIndex") );

        // Copy particle data into buffer
        vector<char> rawParticleData( ostream->particle_size() );
        for( int j = 1; j <= numScans; j++ ) {

            for( int i = 0; i < N; i++ ) {

                if( numScans == 1 ) {
                    posAccessor.set( rawParticleData, vector3fd( cartesianX[i], cartesianY[i], cartesianZ[i] ) );
                } else {
                    posAccessor.set( rawParticleData, scannerTransforms[j - 1].projection_transform(
                                                          vector3fd( cartesianX[i], cartesianY[i], cartesianZ[i] ) ) );
                }
                if( hasIntensity ) {
                    intensityAccessor.get( rawParticleData ) = static_cast<float>( intensity[i] );
                }

                if( hasColor ) {
                    vector3f& colChannel = colorAccessor.get( rawParticleData );
                    if( !colorIsInt ) {
                        colChannel = vector3f( colorRed[i], colorGreen[i], colorBlue[i] );
                    } else {
                        colChannel =
                            vector3f( intColorRed[i] / 255.f, intColorGreen[i] / 255.f, intColorBlue[i] / 255.f );
                    }
                }

                if( hasNormal ) {
                    normalAccessor.get( rawParticleData ) = vector3f( normalX[i], normalY[i], normalZ[i] );
                }
                if( numScans > 1 ) {
                    uint32_t& scanIndexChannel = scanIndexAccessor.get( rawParticleData );
                    scanIndexChannel = j;
                }

                // Stream particle into the prt file
                ostream->put_particle( rawParticleData );
            }
        }

        ostream->close();
    } catch( E57Exception& ex ) {
        ex.report( __FILE__, __LINE__, __FUNCTION__ );
    } catch( std::exception& ex ) {
        cerr << "Got an std::exception, what=" << ex.what() << endl;
    } catch( ... ) {
        cerr << "Got an unknown exception" << endl;
    }
}

#endif
