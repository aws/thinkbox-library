// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/boundbox3f.hpp>
#include <frantic/strings/tstring.hpp>

#include <tinyxml2.h>

#include <boost/filesystem/path.hpp>
#include <boost/rational.hpp>

#include <map>

namespace frantic {
namespace geometry {

class xmesh_metadata {
  public:
    enum transform_type_t {
        transform_type_none,
        transform_type_point,
        transform_type_vector,
        transform_type_normal,
        transform_type_invalid
    };

    enum length_unit_t {
        length_unit_unitless,
        length_unit_inches,
        length_unit_feet,
        length_unit_miles,
        length_unit_millimeters,
        length_unit_centimeters,
        length_unit_meters,
        length_unit_kilometers,
        length_unit_invalid
    };

  private:
    bool hasLengthUnit;
    double lengthUnitScale;
    length_unit_t lengthUnit;

    bool hasFramesPerSecond;
    boost::rational<boost::int64_t> framesPerSecond;

    bool hasBoundbox;
    frantic::graphics::boundbox3f boundbox;

    typedef std::map<frantic::tstring, frantic::tstring> user_data_collection_t;
    user_data_collection_t userData;

    // typedef std::map<std::string,mesh_channel_metadata> channel_metadata_collection_t;
    // channel_metadata_collection_t meshChannelMetadata;

    typedef std::map<frantic::tstring, transform_type_t> channel_transform_type_collection_t;
    channel_transform_type_collection_t channelTransformType;

    // TODO: this shouldn't be a friend
    friend void insert_channel_transform_type_tags( const xmesh_metadata&, tinyxml2::XMLHandle, const std::string&, bool );

  public:
    xmesh_metadata();

    /**
     *  Remove all defined metadata.
     */
    void clear( void );

    // Frames per second
    /**
     *  Remove the frames per second metadata.
     */
    void clear_frames_per_second();
    /**
     *  Set the frames per second.
     *
     *  The frames per second is defined as a rational numerator/denominator
     * because that allows us to express frame rates such as NTSC exactly
     * (30000/1001).
     *
     * @param numerator the numerator part of the frames per second.
     * @param denominator the denominator part of the frames per second.
     */
    void set_frames_per_second( boost::int64_t numerator, boost::int64_t denominator = 1 );
    /**
     *  Set the frames per second.
     *
     * @param framesPerSecond the frames per second.
     */
    void set_frames_per_second( const boost::rational<boost::int64_t>& framesPerSecond );
    /**
     *  Return true if the frames per second is defined, and false otherwise.
     *
     * @return true if the frames per second is defined, and false otherwise.
     */
    bool has_frames_per_second() const;
    /**
     *  Return the frames per second.  If has_frames_per_second() is false,
     * then this value is undefined.
     *
     * @return the frames per second.
     */
    boost::rational<boost::int64_t> get_frames_per_second() const;

    // Boundbox
    /**
     *  Remove the boundbox metadata.
     */
    void clear_boundbox();
    /**
     * Set the boundbox.
     *
     * @param boundbox the boundbox value.
     */
    void set_boundbox( const frantic::graphics::boundbox3f& boundbox );
    /**
     * Return true if the boundbox is defined, and false otherwise.
     *
     * @return true if the boundbox is defined, and false otherwise.
     */
    bool has_boundbox() const;
    /**
     *  Return the boundbox.  If has_boundbox() is false, then this value is
     * undefined.
     *
     * @return the boundbox.
     */
    frantic::graphics::boundbox3f get_boundbox() const;

    // Length unit
    /**
     *  Remove the length unit metadata.
     */
    void clear_length_unit();
    /**
     *  Set the length unit.  The length unit comprises a unit and a multiplier.
     * For example, a unit or 10 cm corresponds to lengthUnitScale 10 and
     * lengthUnit length_unit_centimeters.
     *
     * @param lengthUnitScale the multiple of lengthUnit.
     * @param lengthUnit the length unit.
     */
    void set_length_unit( double lengthUnitScale, length_unit_t lengthUnit );
    /**
     * Return true if the length unit is defined, and false otherwise.
     *
     * @return true if the length unit is defined, and false otherwise.
     */
    bool has_length_unit() const;
    /**
     *  Return the length unit multiplier.
     *
     * @return the length unit multiplier.
     */
    double get_length_unit_scale() const;
    /**
     *  Return the length unit.
     *
     * @return the length unit.
     */
    length_unit_t get_length_unit() const;

    // Channel transform types
    /**
     *  Remove the transform type metadata for all channels.
     */
    void clear_channel_transform_type();
    /**
     *  Set the transform type for the specified channel.
     *
     * @param channelName the channel name to specify a transform type for.
     * @param transformType the transform type for the specified channel.
     */
    void set_channel_transform_type( const frantic::tstring& channelName, transform_type_t transformType );
    /**
     *  Return the transform type for the specified channel.  Defaults to
     * transform_type_none.
     *
     * @param channelName the channel name to find.
     * @return the transform_type_t of the specified channel.
     */
    transform_type_t get_channel_transform_type( const frantic::tstring& channelName ) const;

    // User data
    /**
     *  Remove all user data.
     */
    void clear_user_data();
    /**
     *  Set the user data for the specified key to the specified value.  If the key already
     * exists, its old value is overwritten.
     *
     * @param key the key for the user data.
     * @param value the value for the user data.
     */
    void set_user_data( const frantic::tstring& key, const frantic::tstring& value );
    /**
     * Remove the user data with the specified key.
     *
     * @param key the key for the user data to remove.
     */
    void erase_user_data( const frantic::tstring& key );
    /**
     *  Return true if user data with the specified key exists, and false otherwise.
     *
     * @param key the user data key to find.
     * @return true if user data with the specified key exists, and false otherwise.
     */
    bool has_user_data( const frantic::tstring& key ) const;
    /**
     *  Set the out vector to a list of the defined user data keys.
     *
     * @param[out] out a vector to populate with a list of the defined user data keys.
     */
    void get_user_data_keys( std::vector<frantic::tstring>& out ) const;
    /**
     *  Return the user data value for the specified key.  Throws an exception if
     * there is no user data for the specified key.
     *
     * @param key the key to find.
     * @return the value associated with the key.
     */
    frantic::tstring get_user_data( const frantic::tstring& key ) const;
};

/**
 *  Populate the metadata using values found in the XML document.
 *
 * @param doc the XML document to read.
 * @param[out] outMetadata the metadata to populate.
 */
void read_xmesh_metadata( tinyxml2::XMLDocument& doc, xmesh_metadata& outMetadata );

/**
 *  Populate the metadata using values found in the XML document.
 *
 * @param path the XML file to read.
 * @param[out] outMetadata the metadata to populate.
 */
void read_xmesh_metadata( const boost::filesystem::path& path, xmesh_metadata& outMetadata );

/**
 *  Append metadata to the specified XML document.
 *  This is intended for internal use by classes such as xmesh_sequence_saver
 * and xmesh_writer.
 *  You should normally call this after the document is already populated with
 * xmesh channels and files, etc.
 *
 * @param[out] outDocument append metadata onto this document.
 * @param metadata the metadata to append.
 * @param standalone if true, create a standalone XML file without any
 *        real xmesh data.  This should normally be false.
 */
void write_xmesh_metadata( tinyxml2::XMLDocument& outDocument, const xmesh_metadata& metadata, bool standalone = false );
/**
 *  Create a standalone XML metadata file without any real xmesh data.
 *  This is intended for internal use.
 *
 * @param path the XML filename to create.
 * @param metadata the metadata to write.
 */
void write_xmesh_metadata( const boost::filesystem::path& path, const xmesh_metadata& metadata );

/**
 * Converts a meters scale to an appropriate length unit.
 * The length unit comprises a unit and a multiplier.
 * For example, a unit of 10 cm corresponds to lengthUnitScale 10 and
 * lengthUnit length_unit_centimeters.
 *
 * @param metersScale the unit size in meters.
 * @return a pair containing < the multiple of lengthUnit, the length unit >
 */
std::pair<double, xmesh_metadata::length_unit_t> get_xmesh_length_unit_from_meters( double metersScale );
/**
 * Converts an xmesh length unit to meters
 */
double get_meters_from_xmesh_length_unit( double, frantic::geometry::xmesh_metadata::length_unit_t );
} // namespace geometry
} // namespace frantic
