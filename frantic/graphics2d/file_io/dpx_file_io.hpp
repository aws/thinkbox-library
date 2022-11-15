// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "boost/cstdint.hpp"
#include <fstream>
#include <iostream>
#include "frantic/channels/channel_buffer.hpp"
#include "frantic/math/utils.hpp"

// Fwd decl to break the include cycle dpx_file_io.hpp -> framebuffer.hpp -> repeat.
namespace frantic {
namespace graphics2d {
template <class T>
class framebuffer;
}
} // namespace frantic

namespace frantic {
namespace graphics2d {
namespace file_io {

/**
 * The DPX file format encodes (by default) RGB data using 10-bit's per channel.
*/
class dpx_file_io {
  public:
    dpx_file_io() { clean(); }
    ~dpx_file_io() { clean(); }

    class file_information {
      public:
        boost::uint32_t magic_num;      /* magic number 0x53445058 (SDPX) or 0x58504453 (XPDS) */
        boost::uint32_t offset;         /* offset to image data in bytes */
        boost::uint8_t vers[8];         /* which header format version is being used (v1.0)*/
        boost::uint32_t file_size;      /* file size in bytes */
        boost::uint32_t ditto_key;      /* read time short cut - 0 = same, 1 = new */
        boost::uint32_t gen_hdr_size;   /* generic header length in bytes */
        boost::uint32_t ind_hdr_size;   /* industry header length in bytes */
        boost::uint32_t user_data_size; /* user-defined data length in bytes */
        boost::uint8_t file_name[100];  /* iamge file name */
        boost::uint8_t create_time[24]; /* file creation date "yyyy:mm:dd:hh:mm:ss:LTZ" */
        boost::uint8_t creator[100];    /* file creator's name */
        boost::uint8_t project[200];    /* project name */
        boost::uint8_t copyright[200];  /* right to use or copyright info */
        boost::uint32_t key;            /* encryption ( FFFFFFFF = unencrypted ) */
        boost::uint8_t reserved[104];   /* reserved field TBD (need to pad) */

        void reset();
        void print();
        void convert_to_local_endian( int magic_num );
        void convert_to_file_endian( int magic_num );
    };

    class image_element {
      public:
        boost::uint32_t data_sign;        /* data sign (0 = unsigned, 1 = signed ) */
                                          /* "Core set images are unsigned" */
        boost::uint32_t ref_low_data;     /* reference low data code value */
        float ref_low_quantity;           /* reference low quantity represented */
        boost::uint32_t ref_high_data;    /* reference high data code value */
        float ref_high_quantity;          /* reference high quantity represented */
        boost::uint8_t descriptor;        /* descriptor for image element */
        boost::uint8_t transfer;          /* transfer characteristics for element */
        boost::uint8_t colorimetric;      /* colormetric specification for element */
        boost::uint8_t bit_size;          /* bit size for element */
        boost::uint16_t packing;          /* packing for element */
        boost::uint16_t encoding;         /* encoding for element */
        boost::uint32_t data_offset;      /* offset to data of element */
        boost::uint32_t eol_padding;      /* end of line padding used in element */
        boost::uint32_t eo_image_padding; /* end of image padding used in element */
        boost::uint8_t description[32];   /* description of element */

        void print();
        void convert_to_local_endian( int magic_num );
        void convert_to_file_endian( int magic_num );
    };

    class image_information {
      public:
        boost::uint16_t orientation;         /* image orientation */
        boost::uint16_t element_number;      /* number of image elements */
        boost::uint32_t pixels_per_line;     /* or x value */
        boost::uint32_t lines_per_image_ele; /* or y value, per element */
        image_element element[8];            /* NOTE THERE ARE EIGHT OF THESE */
        boost::uint8_t reserved[52];         /* reserved for future use (padding) */

        void print();
        void convert_to_local_endian( int magic_num );
        void convert_to_file_endian( int magic_num );
    };

    class image_orientation {
      public:
        boost::uint32_t x_offset;         /* X offset */
        boost::uint32_t y_offset;         /* Y offset */
        float x_center;                   /* X center */
        float y_center;                   /* Y center */
        boost::uint32_t x_orig_size;      /* X original size */
        boost::uint32_t y_orig_size;      /* Y original size */
        boost::uint8_t file_name[100];    /* source image file name */
        boost::uint8_t creation_time[24]; /* source image creation date and time */
        boost::uint8_t input_dev[32];     /* input device name */
        boost::uint8_t input_serial[32];  /* input device serial number */
        boost::uint16_t border[4];        /* border validity (XL, XR, YT, YB) */
        boost::uint32_t pixel_aspect[2];  /* pixel aspect ratio (H:V) */
        boost::uint8_t reserved[28];      /* reserved for future use (padding) */

        void print();
        void convert_to_local_endian( int magic_num );
        void convert_to_file_endian( int magic_num );
    };

    class motion_picture_film_header {
      public:
        boost::uint8_t film_mfg_id[2];  /* film manufacturer ID code (2 digits from film edge code) */
        boost::uint8_t film_type[2];    /* file type (2 digits from film edge code) */
        boost::uint8_t offset[2];       /* offset in perfs (2 digits from film edge code)*/
        boost::uint8_t prefix[6];       /* prefix (6 digits from film edge code) */
        boost::uint8_t count[4];        /* count (4 digits from film edge code)*/
        boost::uint8_t format[32];      /* format (i.e. academy) */
        boost::uint32_t frame_position; /* frame position in sequence */
        boost::uint32_t sequence_len;   /* sequence length in frames */
        boost::uint32_t held_count;     /* held count (1 = default) */
        float frame_rate;               /* frame rate of original in frames/sec */
        float shutter_angle;            /* shutter angle of camera in degrees */
        boost::uint8_t frame_id[32];    /* frame identification (i.e. keyframe) */
        boost::uint8_t slate_info[100]; /* slate information */
        boost::uint8_t reserved[56];    /* reserved for future use (padding) */

        void print();
        void convert_to_local_endian( int magic_num );
        void convert_to_file_endian( int magic_num );
    };

    class television_header {
      public:
        boost::uint32_t tim_code;    /* SMPTE time code */
        boost::uint32_t userBits;    /* SMPTE user bits */
        boost::uint8_t interlace;    /* interlace ( 0 = noninterlaced, 1 = 2:1 interlace*/
        boost::uint8_t field_num;    /* field number */
        boost::uint8_t video_signal; /* video signal standard (table 4)*/
        boost::uint8_t unused;       /* used for byte alignment only */
        float hor_sample_rate;       /* horizontal sampling rate in Hz */
        float ver_sample_rate;       /* vertical sampling rate in Hz */
        float frame_rate;            /* temporal sampling rate or frame rate in Hz */
        float time_offset;           /* time offset from sync to first pixel */
        float gamma;                 /* gamma value */
        float black_level;           /* black level code value */
        float black_gain;            /* black gain */
        float break_point;           /* breakpoint */
        float white_level;           /* reference white level code value */
        float integration_times;     /* integration time(s) */
        boost::uint8_t reserved[76]; /* reserved for future use (padding) */

        void print();
        void convert_to_local_endian( int magic_num );
        void convert_to_file_endian( int magic_num );
    };

    class user_defined_data {
      public:
        boost::uint8_t userId[32]; /* User-defined identification string */
        boost::uint8_t* data;      /* User-defined data */

        user_defined_data()
            : data( NULL ) {}

        void print( int userSize = 0 );
    };

    class image_data_element {
      public:
        image_data_element()
            : data( NULL ) {}
        boost::uint8_t* data;

        void convert_to_local_endian( int magic_num, int length );
        void convert_to_file_endian( int magic_num, int length );
    };

    void clean();
    void clean_header();
    void clean_data();
    void fill_default();
    void resize_data( int width, int height );
    bool read_file( const std::string& fileName );
    bool read_file( const std::string& fileName, frantic::channels::channel_buffer* channelBuffer );
    bool write_file( const std::string& fileName );
    bool write_file( const std::string& fileName, frantic::channels::channel_buffer& channelBuffer,
                     bool resetHeader = true );
    void print();

    template <typename T>
    void copy_to_framebuffer( framebuffer<T>* frame );
    template <typename T>
    void copy_from_framebuffer( framebuffer<T>& frame, bool resetHeader = true );

    void copy_to_channel_buffer( frantic::channels::channel_buffer* frame );
    void copy_from_channel_buffer( frantic::channels::channel_buffer& frame, bool resetHeader = true );

    file_information& get_file_information() { return m_fileInformation; }
    image_information& get_image_information() { return m_imageInformation; }
    image_orientation& get_image_orientation() { return m_imageOrientation; }
    motion_picture_film_header& get_motion_picture_film_header() { return m_motionPictureFilmHeader; }
    television_header& get_television_header() { return m_televisionHeader; }
    user_defined_data& get_user_defined_data() { return m_userDefinedData; }
    image_data_element& get_image_data_element() { return m_imageDataElement; }

  private:
    file_information m_fileInformation;
    image_information m_imageInformation;
    image_orientation m_imageOrientation;
    motion_picture_film_header m_motionPictureFilmHeader;
    television_header m_televisionHeader;
    user_defined_data m_userDefinedData;
    image_data_element m_imageDataElement;
};

} // namespace file_io
} // namespace graphics2d
} // namespace frantic

#include <frantic/graphics2d/framebuffer.hpp>
#include <frantic/graphics2d/framebufferiterator.hpp>

namespace frantic {
namespace graphics2d {
namespace file_io {

template <typename T>
void dpx_file_io::copy_to_framebuffer( framebuffer<T>* frame ) {
    if( m_imageDataElement.data != NULL ) {
        graphics::color3f pixel;
        unsigned int i = 0;

        frame->set_size( size2( m_imageInformation.pixels_per_line, m_imageInformation.lines_per_image_ele ) );
        framebufferiterator<T> iter( frame );
        iter.goRow( frame->height() - 1 );

        while( iter.isRowValid() &&
               ( i < m_imageInformation.pixels_per_line * m_imageInformation.lines_per_image_ele ) ) {
            while( iter.isColValid() ) {
                boost::uint32_t read = ( (boost::uint32_t*)m_imageDataElement.data )[i];

                pixel.r = ( ( read >> 22 ) & 0x3ff ) / 1023.0f;
                pixel.g = ( ( read >> 12 ) & 0x3ff ) / 1023.0f;
                pixel.b = ( ( read >> 2 ) & 0x3ff ) / 1023.0f;

                iter.setData( pixel );

                i++;
                iter.nextCol();
            }
            iter.prevRowStart();
        }
    } else
        std::cout << __FILE__ << ":" << __LINE__ << " - Dpx image contains no information." << std::endl;
}

template <typename T>
void dpx_file_io::copy_from_framebuffer( framebuffer<T>& frame, bool resetHeader ) {
    graphics::color3f pixel;
    unsigned int i = 0;
    boost::uint32_t read = 0;

    if( resetHeader )
        fill_default();

    resize_data( frame.width(), frame.height() );

    framebufferiterator<T> iter( &frame );
    iter.goRow( frame.height() - 1 );

    while( iter.isRowValid() && ( i < m_imageInformation.pixels_per_line * m_imageInformation.lines_per_image_ele ) ) {
        while( iter.isColValid() ) {
            iter.getData( &pixel );

            read = ( ( boost::uint32_t )( pixel.r * 1023.0f ) ) << 22;
            read |= ( ( boost::uint32_t )( pixel.g * 1023.0f ) ) << 12;
            read |= ( ( boost::uint32_t )( pixel.b * 1023.0f ) ) << 2;

            ( (boost::uint32_t*)m_imageDataElement.data )[i] = read;

            i++;
            iter.nextCol();
        }
        iter.prevRowStart();
    }
}

} // namespace file_io
} // namespace graphics2d
} // namespace frantic
