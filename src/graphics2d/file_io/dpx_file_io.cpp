// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>

#include "boost/cstdint.hpp"
#include "frantic/channels/channel_buffer.hpp"
#include "frantic/graphics/color3f.hpp"
#include "frantic/graphics/color_converter.hpp"
#include "frantic/graphics/color_rgb_f.hpp"
#include "frantic/graphics/color_rgba_f.hpp"
#include "frantic/graphics2d/file_io/dpx_file_io.hpp"
#include "frantic/graphics2d/framebuffer.hpp"
#include "frantic/graphics2d/framebufferiterator.hpp"
#include "frantic/math/utils.hpp"

using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::graphics2d::file_io;
using namespace frantic::math;

void dpx_file_io::file_information::reset() { memset( this, 0, sizeof( file_information ) ); }

void dpx_file_io::file_information::print() {
    std::cout << "file_information" << std::endl;
    std::cout << "\tmagic_num:           " << std::hex << magic_num << std::dec << std::endl;
    std::cout << "\toffset:              " << offset << std::endl;
    std::cout << "\tvers:                " << (char*)vers << std::endl;
    std::cout << "\tfile_size:           " << file_size << std::endl;
    std::cout << "\tditto_key:           " << ditto_key << std::endl;
    std::cout << "\tgen_hdr_size:        " << gen_hdr_size << std::endl;
    std::cout << "\tind_hdr_size:        " << ind_hdr_size << std::endl;
    std::cout << "\tuser_data_size:      " << user_data_size << std::endl;
    std::cout << "\tfile_name:           " << (char*)file_name << std::endl;
    std::cout << "\tcreate_time:         " << (char*)create_time << std::endl;
    std::cout << "\tcreator:             " << (char*)creator << std::endl;
    std::cout << "\tproject:             " << (char*)project << std::endl;
    std::cout << "\tcopyright:           " << (char*)copyright << std::endl;
    std::cout << "\tkey:                 " << std::hex << key << std::dec << std::endl;
}

void dpx_file_io::file_information::convert_to_local_endian( int magic_num ) {
    if( magic_num == 0x58504453 ) {
        math::byte_swap( offset );
        math::byte_swap( file_size );
        math::byte_swap( ditto_key );
        math::byte_swap( gen_hdr_size );
        math::byte_swap( ind_hdr_size );
        math::byte_swap( user_data_size );
        math::byte_swap( key );
    }
}

void dpx_file_io::file_information::convert_to_file_endian( int magic_num ) { convert_to_local_endian( magic_num ); }

void dpx_file_io::image_element::print() {
    std::cout << "\timage_element" << std::endl;
    std::cout << "\t\tdata_sign:           " << data_sign << std::endl;
    std::cout << "\t\tref_low_data:        " << ref_low_data << std::endl;
    std::cout << "\t\tref_low_quantity:    " << ref_low_quantity << std::endl;
    std::cout << "\t\tref_high_data:       " << ref_high_data << std::endl;
    std::cout << "\t\tref_high_quantity:   " << ref_high_quantity << std::endl;
    std::cout << "\t\tdescriptor:          " << (int)descriptor << std::endl;
    std::cout << "\t\ttransfer:            " << (int)transfer << std::endl;
    std::cout << "\t\tcolorimetric:        " << (int)colorimetric << std::endl;
    std::cout << "\t\tbit_size:            " << (int)bit_size << std::endl;
    std::cout << "\t\tpacking:             " << packing << std::endl;
    std::cout << "\t\tencoding:            " << encoding << std::endl;
    std::cout << "\t\tdata_offset:         " << data_offset << std::endl;
    std::cout << "\t\teol_padding:         " << eol_padding << std::endl;
    std::cout << "\t\teo_image_padding:    " << eo_image_padding << std::endl;
    std::cout << "\t\tdescription:         " << (char*)description << std::endl;
}

void dpx_file_io::image_element::convert_to_local_endian( int magic_num ) {
    if( magic_num == 0x58504453 ) {
        math::byte_swap( data_sign );
        math::byte_swap( ref_low_data );
        math::byte_swap( ref_low_quantity );
        math::byte_swap( ref_high_data );
        math::byte_swap( ref_high_quantity );
        math::byte_swap( packing );
        math::byte_swap( encoding );
        math::byte_swap( data_offset );
        math::byte_swap( eol_padding );
        math::byte_swap( eo_image_padding );
    }
}

void dpx_file_io::image_element::convert_to_file_endian( int magic_num ) { convert_to_local_endian( magic_num ); }

void dpx_file_io::image_information::print() {
    std::cout << "image_information" << std::endl;
    std::cout << "\torientation:         " << orientation << std::endl;
    std::cout << "\telement_number:      " << element_number << std::endl;
    std::cout << "\tpixels_per_line:     " << pixels_per_line << std::endl;
    std::cout << "\tlines_per_image_ele: " << lines_per_image_ele << std::endl;
    for( int i = 0; i < element_number; i++ )
        element[i].print();
}

void dpx_file_io::image_information::convert_to_local_endian( int magic_num ) {
    if( magic_num == 0x58504453 ) {
        math::byte_swap( orientation );
        math::byte_swap( element_number );
        math::byte_swap( pixels_per_line );
        math::byte_swap( lines_per_image_ele );
        for( int i = 0; i < element_number; i++ )
            element[i].convert_to_local_endian( magic_num );
    }
}

void dpx_file_io::image_information::convert_to_file_endian( int magic_num ) {
    if( magic_num == 0x58504453 ) {
        math::byte_swap( orientation );
        math::byte_swap( pixels_per_line );
        math::byte_swap( lines_per_image_ele );
        for( int i = 0; i < element_number; i++ )
            element[i].convert_to_file_endian( magic_num );
        math::byte_swap( element_number );
    }
}

void dpx_file_io::image_orientation::print() {
    std::cout << "image_orientation" << std::endl;
    std::cout << "\tx_offset:            " << x_offset << std::endl;
    std::cout << "\ty_offset:            " << y_offset << std::endl;
    std::cout << "\tx_center:            " << x_center << std::endl;
    std::cout << "\ty_center:            " << y_center << std::endl;
    std::cout << "\tx_orig_size:         " << x_orig_size << std::endl;
    std::cout << "\ty_orig_size:         " << y_orig_size << std::endl;
    std::cout << "\tfile_name:           " << (char*)file_name << std::endl;
    std::cout << "\tcreation_time:       " << (char*)creation_time << std::endl;
    std::cout << "\tinput_dev:           " << (char*)input_dev << std::endl;
    std::cout << "\tinput_serial:        " << (char*)input_serial << std::endl;
    for( int i = 0; i < 4; i++ )
        std::cout << "\tborder " << i << ":            " << border[i] << std::endl;
    for( int i = 0; i < 2; i++ )
        std::cout << "\tpixel_aspect " << i << ":      " << pixel_aspect[i] << std::endl;
}

void dpx_file_io::image_orientation::convert_to_local_endian( int magic_num ) {
    if( magic_num == 0x58504453 ) {
        math::byte_swap( x_offset );
        math::byte_swap( y_offset );
        math::byte_swap( x_center );
        math::byte_swap( y_center );
        math::byte_swap( x_orig_size );
        math::byte_swap( y_orig_size );
        for( int i = 0; i < 4; i++ )
            math::byte_swap( border[i] );
        for( int i = 0; i < 2; i++ )
            math::byte_swap( pixel_aspect[i] );
    }
}

void dpx_file_io::image_orientation::convert_to_file_endian( int magic_num ) { convert_to_local_endian( magic_num ); }

void dpx_file_io::motion_picture_film_header::print() {
    std::cout << "motion_picture_film_header" << std::endl;
    std::cout << "\tfilm_mfg_id:         ( ";
    for( int i = 0; i < 2; i++ )
        std::cout << (int)film_mfg_id[i] << " ";
    std::cout << ")" << std::endl;
    std::cout << "\tfilm_type:           ( ";
    for( int i = 0; i < 2; i++ )
        std::cout << (int)film_type[i] << " ";
    std::cout << ")" << std::endl;
    std::cout << "\toffset:              ( ";
    for( int i = 0; i < 2; i++ )
        std::cout << (int)offset[i] << " ";
    std::cout << ")" << std::endl;
    std::cout << "\tprefix:              ( ";
    for( int i = 0; i < 2; i++ )
        std::cout << (int)prefix[i] << " ";
    std::cout << ")" << std::endl;
    std::cout << "\tcount:               ( ";
    for( int i = 0; i < 2; i++ )
        std::cout << (int)count[i] << " ";
    std::cout << ")" << std::endl;
    std::cout << "\tformat:              " << (char*)format << std::endl;
    std::cout << "\tframe_position:      " << frame_position << std::endl;
    std::cout << "\tsequence_len:        " << sequence_len << std::endl;
    std::cout << "\theld_count:          " << held_count << std::endl;
    std::cout << "\tframe_rate:          " << frame_rate << std::endl;
    std::cout << "\tshutter_angle:       " << shutter_angle << std::endl;
    std::cout << "\tframe_id:            " << (char*)frame_id << std::endl;
    std::cout << "\tslate_info:          " << (char*)slate_info << std::endl;
}

void dpx_file_io::motion_picture_film_header::convert_to_local_endian( int magic_num ) {
    if( magic_num == 0x58504453 ) {
        math::byte_swap( frame_position );
        math::byte_swap( sequence_len );
        math::byte_swap( held_count );
        math::byte_swap( frame_rate );
        math::byte_swap( shutter_angle );
    }
}

void dpx_file_io::motion_picture_film_header::convert_to_file_endian( int magic_num ) {
    convert_to_local_endian( magic_num );
}

void dpx_file_io::television_header::print() {
    std::cout << "television_header" << std::endl;
    std::cout << "\ttim_code:            " << tim_code << std::endl;
    std::cout << "\tuserBits:            " << userBits << std::endl;
    std::cout << "\tinterlace:           " << (int)interlace << std::endl;
    std::cout << "\tfield_num:           " << (int)field_num << std::endl;
    std::cout << "\tvideo_signal:        " << (int)video_signal << std::endl;
    std::cout << "\thor_sample_rate:     " << hor_sample_rate << std::endl;
    std::cout << "\tver_sample_rate:     " << ver_sample_rate << std::endl;
    std::cout << "\tframe_rate:          " << frame_rate << std::endl;
    std::cout << "\ttime_offset:         " << time_offset << std::endl;
    std::cout << "\tgamma:               " << gamma << std::endl;
    std::cout << "\tblack_level:         " << black_level << std::endl;
    std::cout << "\tblack_gain:          " << black_gain << std::endl;
    std::cout << "\tbreak_point:         " << break_point << std::endl;
    std::cout << "\twhite_level:         " << white_level << std::endl;
    std::cout << "\tintegration_times:   " << integration_times << std::endl;
}

void dpx_file_io::television_header::convert_to_local_endian( int magic_num ) {
    if( magic_num == 0x58504453 ) {
        math::byte_swap( tim_code );
        math::byte_swap( userBits );
        math::byte_swap( hor_sample_rate );
        math::byte_swap( ver_sample_rate );
        math::byte_swap( frame_rate );
        math::byte_swap( time_offset );
        math::byte_swap( gamma );
        math::byte_swap( black_level );
        math::byte_swap( black_gain );
        math::byte_swap( break_point );
        math::byte_swap( white_level );
        math::byte_swap( integration_times );
    }
}

void dpx_file_io::television_header::convert_to_file_endian( int magic_num ) { convert_to_local_endian( magic_num ); }

void dpx_file_io::user_defined_data::print( int userSize ) {
    std::cout << "user_defined_data" << std::endl;
    std::cout << "\tuserId:              " << (char*)userId << std::endl;
    std::cout << "\tdata:" << std::endl;
    for( int i = 0; i < userSize; i++ )
        std::cout << "\t\t" << std::hex << i << "\t" << data[i] << "\t" << std::hex << (int)data[i] << "\t" << std::dec
                  << (int)data[i] << std::endl;
}

void dpx_file_io::image_data_element::convert_to_local_endian( int magic_num, int length ) {
    if( magic_num == 0x58504453 ) {
        for( int i = 0; i < length; i++ ) {
            math::byte_swap( *reinterpret_cast<boost::uint32_t*>( &data[4 * i] ) );
        }
    }
}

void dpx_file_io::image_data_element::convert_to_file_endian( int magic_num, int length ) {
    convert_to_local_endian( magic_num, length );
}

void dpx_file_io::clean() {
    clean_header();
    clean_data();
}

void dpx_file_io::clean_header() {
    memset( &m_fileInformation, 0, sizeof( file_information ) );
    memset( &m_imageInformation, 0, sizeof( image_information ) );
    memset( &m_imageOrientation, 0, sizeof( image_orientation ) );
    memset( &m_motionPictureFilmHeader, 0, sizeof( motion_picture_film_header ) );
    memset( &m_televisionHeader, 0, sizeof( television_header ) );
    if( m_userDefinedData.data != NULL )
        delete[] m_userDefinedData.data;
    memset( &m_userDefinedData, 0, sizeof( user_defined_data ) );
}

void dpx_file_io::clean_data() {
    if( m_imageDataElement.data != NULL ) {
        delete[] m_imageDataElement.data;
        m_imageDataElement.data = NULL;
    }
}

void dpx_file_io::fill_default() {
    clean();

    // Verify the image format header
    if( ( m_fileInformation.magic_num != 0x53445058 ) || ( m_fileInformation.magic_num != 0x58504453 ) )
        m_fileInformation.magic_num = 0x53445058;

    memcpy( m_fileInformation.vers, "V1.0", 5 );

    m_fileInformation.offset = 8192;
    m_fileInformation.gen_hdr_size =
        sizeof( m_fileInformation ) + sizeof( m_imageInformation ) + sizeof( m_imageOrientation );
    m_fileInformation.ind_hdr_size = sizeof( m_motionPictureFilmHeader ) + sizeof( m_televisionHeader );

    // Set the file size
    m_fileInformation.file_size =
        m_fileInformation.offset + m_imageInformation.pixels_per_line * m_imageInformation.lines_per_image_ele * 4;

    m_imageInformation.element_number = 1;
    m_imageInformation.element[0].data_sign = 0;
    m_imageInformation.element[0].ref_low_data = 0;
    m_imageInformation.element[0].ref_high_data = 1023;
    m_imageInformation.element[0].descriptor = 50;
    m_imageInformation.element[0].packing = 1;
    m_imageInformation.element[0].bit_size = 10;
    m_imageInformation.element[0].encoding = 0;
    m_imageInformation.element[0].data_offset = m_fileInformation.offset;
    m_imageInformation.element[0].eol_padding = 0;
}

void dpx_file_io::resize_data( int width, int height ) {
    if( m_imageDataElement.data != NULL ) {
        delete[] m_imageDataElement.data;
        m_imageDataElement.data = NULL;
    }

    m_imageDataElement.data = new boost::uint8_t[width * height * 4];

    m_imageInformation.pixels_per_line = width;
    m_imageInformation.lines_per_image_ele = height;
}

/** An overloaded function call that reads the file into memory, copys
 * the image into a channel_buffer, then free the memory and cleans the header
 *
 * @arg	fileName		The file name
 * @arg	channelBuffer	The channel_buffer where the image is placed.  If need it will
 *							be resized
 * @return					Returns true if the file was successfully read.
 */
bool dpx_file_io::read_file( const std::string& fileName, channel_buffer* channelBuffer ) {
    bool ret = this->read_file( fileName );
    if( ret ) {
        this->copy_to_channel_buffer( channelBuffer );
    }
    this->clean();
    return ret;
}

/** An overloaded function call that copys the channel_buffer into memory, writes
 * the file, then free the memory and cleans the header
 *
 * @arg	fileName		The file name
 * @arg	channelBuffer	The channel_buffer where the image is placed.  If need it will
 *							be resized
 * @arg	resetHeader		If true a default header will be created, false will preserve
 *							existing header data and only modify img size values
 * @return					Returns true if the file was successfully written.
 */
bool dpx_file_io::write_file( const std::string& fileName, channel_buffer& channelBuffer, bool resetHeader ) {
    bool ret;
    this->copy_from_channel_buffer( channelBuffer, resetHeader );
    ret = this->write_file( fileName );
    this->clean();
    return ret;
}

bool dpx_file_io::read_file( const std::string& fileName ) {
    std::ifstream inFile;
    clean_data();

    inFile.open( fileName.c_str(), std::ifstream::in | std::ifstream::binary );

    if( inFile.good() ) {
        inFile.read( (char*)&m_fileInformation, sizeof( m_fileInformation ) );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to open dpx file " << fileName;
        throw std::runtime_error( strstm.str() );
    }
    if( inFile.good() ) {
        // Verify the image format and size
        if( ( m_fileInformation.magic_num == 0x53445058 ) || ( m_fileInformation.magic_num == 0x58504453 ) ) {
            m_fileInformation.convert_to_local_endian( m_fileInformation.magic_num );
        } else {
            std::stringstream strstm;
            strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Dpx magic number failed for "
                   << fileName;
            throw std::runtime_error( strstm.str() );
        }
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to read file information from "
               << fileName;
        throw std::runtime_error( strstm.str() );
    }
    if( inFile.good() ) {
        inFile.read( (char*)&m_imageInformation, sizeof( m_imageInformation ) );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to read image information from "
               << fileName;
        throw std::runtime_error( strstm.str() );
    }
    if( inFile.good() ) {
        m_imageInformation.convert_to_local_endian( m_fileInformation.magic_num );
        inFile.read( (char*)&m_imageOrientation, sizeof( m_imageOrientation ) );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to read image orientation from "
               << fileName;
        throw std::runtime_error( strstm.str() );
    }
    if( inFile.good() ) {
        m_imageOrientation.convert_to_local_endian( m_fileInformation.magic_num );
        inFile.read( (char*)&m_motionPictureFilmHeader, sizeof( m_motionPictureFilmHeader ) );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to read motion picture info from "
               << fileName;
        throw std::runtime_error( strstm.str() );
    }
    if( inFile.good() ) {
        m_motionPictureFilmHeader.convert_to_local_endian( m_fileInformation.magic_num );
        inFile.read( (char*)&m_televisionHeader, sizeof( m_televisionHeader ) );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to read television info from "
               << fileName;
        throw std::runtime_error( strstm.str() );
    }
    if( inFile.good() ) {
        m_televisionHeader.convert_to_local_endian( m_fileInformation.magic_num );
        inFile.read( (char*)&m_userDefinedData.userId, 32 );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__
               << " - Failed to read user defined data header from " << fileName;
        throw std::runtime_error( strstm.str() );
    }
    if( inFile.good() ) {
        if( ( (int)m_fileInformation.user_data_size ) - 32 > 0 ) {
            m_userDefinedData.data = new boost::uint8_t[m_fileInformation.user_data_size - 32];
            inFile.read( (char*)m_userDefinedData.data, m_fileInformation.user_data_size - 32 );
        }
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to read user defined data from "
               << fileName;
        throw std::runtime_error( strstm.str() );
    }
    inFile.seekg( m_fileInformation.offset );
    if( inFile.good() ) {
        m_imageDataElement.data =
            new boost::uint8_t[m_imageInformation.pixels_per_line * m_imageInformation.lines_per_image_ele * 4];
        inFile.read( (char*)m_imageDataElement.data,
                     m_imageInformation.pixels_per_line * m_imageInformation.lines_per_image_ele * 4 );
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to seek to image from " << fileName;
        throw std::runtime_error( strstm.str() );
    }
    if( inFile.good() ) {
        m_imageDataElement.convert_to_local_endian(
            m_fileInformation.magic_num, m_imageInformation.pixels_per_line * m_imageInformation.lines_per_image_ele );
        return true;
    } else {
        clean_data();

        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to read image from " << fileName;
        throw std::runtime_error( strstm.str() );

        return false;
    }
}

bool dpx_file_io::write_file( const std::string& fileName ) {
    std::ofstream outFile;

    outFile.open( fileName.c_str(), std::ifstream::out | std::ofstream::binary );

    if( outFile.good() ) {
        m_fileInformation.convert_to_file_endian( m_fileInformation.magic_num );
        outFile.write( (char*)&m_fileInformation, sizeof( m_fileInformation ) );
        m_fileInformation.convert_to_local_endian( m_fileInformation.magic_num );
    }
    if( outFile.good() ) {
        m_imageInformation.convert_to_file_endian( m_fileInformation.magic_num );
        outFile.write( (char*)&m_imageInformation, sizeof( m_imageInformation ) );
        m_imageInformation.convert_to_local_endian( m_fileInformation.magic_num );
    }
    if( outFile.good() ) {
        m_imageOrientation.convert_to_file_endian( m_fileInformation.magic_num );
        outFile.write( (char*)&m_imageOrientation, sizeof( m_imageOrientation ) );
        m_imageOrientation.convert_to_local_endian( m_fileInformation.magic_num );
    }
    if( outFile.good() ) {
        m_motionPictureFilmHeader.convert_to_file_endian( m_fileInformation.magic_num );
        outFile.write( (char*)&m_motionPictureFilmHeader, sizeof( m_motionPictureFilmHeader ) );
        m_motionPictureFilmHeader.convert_to_local_endian( m_fileInformation.magic_num );
    }
    if( outFile.good() ) {
        m_televisionHeader.convert_to_file_endian( m_fileInformation.magic_num );
        outFile.write( (char*)&m_televisionHeader, sizeof( m_televisionHeader ) );
        m_televisionHeader.convert_to_local_endian( m_fileInformation.magic_num );
    }
    if( m_fileInformation.user_data_size > 0 ) {
        if( outFile.good() ) {
            outFile.write( (char*)&m_userDefinedData.userId, 32 );
        }
        if( outFile.good() ) {
            if( m_userDefinedData.data != NULL )
                outFile.write( (char*)m_userDefinedData.data, m_fileInformation.user_data_size - 32 );
        }
    }
    while( outFile.tellp() < (int)m_fileInformation.offset )
        outFile.put( 0 );
    if( outFile.good() ) {
        if( ( m_imageInformation.pixels_per_line > 0 ) && ( m_imageInformation.lines_per_image_ele > 0 ) &&
            ( m_imageDataElement.data != NULL ) ) {
            m_imageDataElement.convert_to_file_endian( m_fileInformation.magic_num,
                                                       m_imageInformation.pixels_per_line *
                                                           m_imageInformation.lines_per_image_ele );
            outFile.write( (char*)m_imageDataElement.data,
                           m_imageInformation.pixels_per_line * m_imageInformation.lines_per_image_ele * 4 );
            m_imageDataElement.convert_to_local_endian( m_fileInformation.magic_num,
                                                        m_imageInformation.pixels_per_line *
                                                            m_imageInformation.lines_per_image_ele );
        }
    }
    if( outFile.good() ) {
        return true;
    } else {
        clean_data();
        return false;
    }
}

void dpx_file_io::print() {
    if( m_imageDataElement.data != NULL ) {
        m_fileInformation.print();
        m_imageInformation.print();
        m_imageOrientation.print();
        m_motionPictureFilmHeader.print();
        m_televisionHeader.print();
        m_userDefinedData.print( m_fileInformation.user_data_size - 32 );
    } else
        std::cout << __FILE__ << ":" << __LINE__ << " - Dpx image contains no information." << std::endl;
}

void dpx_file_io::copy_to_channel_buffer( channel_buffer* frame ) {
    if( m_imageDataElement.data != NULL ) {
        try {
            graphics::color_rgba_f pixel;

            frame->resize( size2( m_imageInformation.pixels_per_line, m_imageInformation.lines_per_image_ele ) );
            generic_channel_buffer_iterator<color_rgba_f> iter =
                frame->get_generic_iterator<color_rgba_f, color_converter>();

            iter.go_y( m_imageInformation.lines_per_image_ele - 1 );
            for( unsigned int row = 0, i = 0; row < m_imageInformation.lines_per_image_ele; row++ ) {
                for( unsigned int col = 0; col < m_imageInformation.pixels_per_line; col++, i++ ) {
                    boost::uint32_t read = ( (boost::uint32_t*)m_imageDataElement.data )[i];

                    pixel.set_r( ( ( read >> 22 ) & 0x3ff ) / 1023.0f );
                    pixel.set_g( ( ( read >> 12 ) & 0x3ff ) / 1023.0f );
                    pixel.set_b( ( ( read >> 2 ) & 0x3ff ) / 1023.0f );

                    iter.set( pixel );

                    iter.next_x();
                }
                iter.go_x( 0 );
                iter.prev_y();
            }
        } catch( const std::exception& e ) {
            std::cout << e.what() << std::endl;
            std::stringstream strstm;
            strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to copy data";
            throw std::runtime_error( strstm.str() );
        }
    } else {
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Dpx image contains no information";
        throw std::runtime_error( strstm.str() );
    }
}

void dpx_file_io::copy_from_channel_buffer( channel_buffer& frame, bool resetHeader ) {
    graphics::color_rgb_f pixel;
    boost::uint32_t read = 0;

    try {
        if( resetHeader )
            fill_default();

        resize_data( frame.get_size2().xsize, frame.get_size2().ysize );

        generic_channel_buffer_iterator<color_rgb_f> iter = frame.get_generic_iterator<color_rgb_f, color_converter>();

        iter.go_y( m_imageInformation.lines_per_image_ele - 1 );
        for( unsigned int row = 0, i = 0; row < m_imageInformation.lines_per_image_ele; row++ ) {
            for( unsigned int col = 0; col < m_imageInformation.pixels_per_line; col++, i++ ) {
                pixel = iter.get();

                read = ( ( boost::uint32_t )( pixel.get_r() * 1023.0f ) ) << 22;
                read |= ( ( boost::uint32_t )( pixel.get_g() * 1023.0f ) ) << 12;
                read |= ( ( boost::uint32_t )( pixel.get_b() * 1023.0f ) ) << 2;

                ( (boost::uint32_t*)m_imageDataElement.data )[i] = read;

                iter.next_x();
            }
            iter.go_x( 0 );
            iter.prev_y();
        }
    } catch( const std::exception& e ) {
        std::cout << e.what() << std::endl;
        std::stringstream strstm;
        strstm << __FILE__ << "@" << __LINE__ << "." << __FUNCTION__ << " - Failed to copy data";
        throw std::runtime_error( strstm.str() );
    }
}
