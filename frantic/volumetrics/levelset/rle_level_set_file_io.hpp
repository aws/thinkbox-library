// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/files/filename_sequence.hpp>
#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {
/**
 * This function loads from the file a rle voxel field.
 *
 * @param  file           The filename of the file to read.
 * @param  outVoxelField  The level set to populate with the data in the file.
 */
void read_rle_voxel_field_file( const frantic::tstring& file, frantic::fluids::rle_voxel_field& outVoxelField );

/**
 * This function loads from the file a rle voxel field with a channel propagation policy.
 *
 * @param  file           The filename of the file to read.
 * @param  policy         The channel propagation policy
 * @param  outVoxelField  The level set to populate with the data in the file.
 */
void read_rle_voxel_field_file( const frantic::tstring& file,
                                const frantic::channels::channel_propagation_policy& policy,
                                frantic::fluids::rle_voxel_field& outVoxelField );

/**
 * This function writes to an rle voxel field file in the .rls file format.
 *
 * @param  file        The filename of the file to write.
 * @param  voxelField  The voxel field to write to the file.
 */
void write_rle_voxel_field_file( const frantic::tstring& file, const frantic::fluids::rle_voxel_field& voxelField );

/**
 * This function writes to an rle voxel field file in the .rls file format with a channel
 * propagation policy.
 *
 * @param  file        The filename of the file to write.
 * @param  policy      The channel propagation policy
 * @param  voxelField  The lvoxel field to write to the file.
 */
void write_rle_voxel_field_file( const frantic::tstring& file,
                                 const frantic::channels::channel_propagation_policy& policy,
                                 const frantic::fluids::rle_voxel_field& voxelField );

/**
 * This function uses the file extension of the input file to call the appropriate function for
 * loading the level set file.
 *
 * @param  file         The filename of the file to read.
 * @param  outLevelSet  The level set to populate with the data in the file.
 */
void read_rle_level_set_file( const frantic::tstring& file, levelset::rle_level_set& outLevelSet );

/**
 * This function reads in an rle level set file in the .rls file format.
 *
 * @param  file         The filename of the file to read.
 * @param  outLevelSet  The level set to populate with the data in the file.
 */
void read_rls_rle_level_set_file( const frantic::tstring& file, levelset::rle_level_set& outLevelSet );

/**
 * This function uses the file extension of the input file to call the appropriate function for
 * loading the level set "header" info.  This includes voxel coord system, interfacewidths and
 * channel info.
 *
 * @param  file         The filename of the file to read.
 * @param  outLevelSet  The level set to populate with the header data in the file.
 */
void read_rle_level_set_file_header( const frantic::tstring& file, levelset::rle_level_set& outLevelSet );

/**
 * This function reads rle level set file header info into a level set in the .rls file format.
 *
 * @param  file         The filename of the file to read.
 * @param  outLevelSet  The level set to populate with the header data in the file.
 */
void read_rls_rle_level_set_file_header( const frantic::tstring& file, rle_level_set& outLevelSet );

/**
 * This function writes to an rle level set file in the .rls file format.
 *
 * @param  file      The filename of the file to write.
 * @param  levelSet  The level set to write to the file.
 */
void write_rls_rle_level_set_file( const frantic::tstring& file, const levelset::rle_level_set& levelSet );

/**
 * This function interpolates between two level set files at the given fractional alpha value.  If the
 * file names are the same no interpolation takes place and the first file is simple loaded and returned.
 * If the alpha value is not between 0 and 1, it throws an exception.
 *
 * @param  firstFile	The filename of the first level set file
 * @param  secondFile	The filename of the second level set file
 * @param  alpha		The float alpha value between 0 and 1
 * @param  outLevelSet	The interpolated level set.
 */
void interpolate_rls_rle_level_set_files( const frantic::tstring& firstFile, const frantic::tstring& secondFile,
                                          const float alpha, levelset::rle_level_set& outLevelSet );

/**
 * This class contains functionality to transparently manage synching of
 * rls files from a network location to the local temp folder.
 *
 * @todo DEPRECATED
 * This was for the old version of the max scene source to help with caching network files during journey.
 * it is still being used in the old flood:spray and should be removed when the python version is completed
 */
class rls_network_cache {

  private:
    frantic::files::filename_sequence m_fsq;
    frantic::tstring m_tempDir;
    std::vector<frantic::tstring> filesCopiedLocally;
    bool m_initialized;

  public:
    /// DEPRECATED
    rls_network_cache()
        : m_initialized( false ) {}

    /**
     * Constructor that initializes the member variables.
     *
     * @todo DEPRECATED
     *
     * @param  tempDir		the location of the temp folder to copy the rls files to
     * @param  cacheDir		the location of one of the files in the source sequence
     */
    rls_network_cache( const frantic::tstring& tempDir, const frantic::tstring& cacheDir );

    /**
     * Initializes member variables. These variables need to be initalized here, or in the overloaded constuctor before
     * using the object.
     *
     * @todo DEPRECATED
     *
     * @param  tempDir		the location of the temp folder to copy the rls files to
     * @param  cacheDir		the location of one of the files in the source sequence
     */
    void initialize( const frantic::tstring& tempDir, const frantic::tstring& cacheDir );

    /**
     * Function to see if level set cache is initialized
     *
     * @todo DEPRECATED
     *
     * @return bool		true if initialized
     */
    bool is_initialized() const { return m_initialized; }

    ~rls_network_cache();

    /**
     * Retrieves a level from the cache at the provided path at the time requested, interpolating if required
     * between the two closest frames.  The level set must match the provided coordinate system or exceptions
     * will be thrown.
     *
     * @todo DEPRECATED
     *
     * @param  frame			the frame at which to get the level set
     * @param  interfaceVoxelWidthInside	the inside width to match
     * @param  interfaceVoxelWidthOutside	the outside width to match
     * @param[out] outLevelSet	the output level set, already with a vcs to match
     * @param  useCacheSettings  flag for the new level set use the cache voxel coord system and interface voxel widths,
     * even though it is different than the provided settings
     * @return bool			true if everything matched and the cache was loaded from, false otherwise
     */
    bool get_level_set( double frame, float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside,
                        frantic::volumetrics::levelset::rle_level_set& outLevelSet, bool useCacheSettings = false );

    frantic::files::filename_pattern& get_filename_pattern() { return m_fsq.get_filename_pattern(); }

  private:
    /**
     * Retrieves the subframe numbers from the given filename sequence which are the closest bracket
     * to the provided timevalue.
     *
     * @todo DEPRECATED
     *
     * @param  fsq			the filename sequence to be checked
     * @param  frame			the frame at which to get the closest neighbouring frames
     * @param[out] interval		a standard pair of doubles indicating the frame numbers in the fsq of the
     * closest files
     * @param[out] alpha		the fractional alpha value distance between the frame numbers the requested time was
     * as (e.g. for interpolation)
     * @return bool			true if it could find appropriate frames, false otherwise
     */
    bool get_nearest_subframe_interval( frantic::files::filename_sequence fsq, double frame,
                                        std::pair<double, double>& interval, float& alpha );

    /** Copies a cache file at a given location to the local temp directory.  Uses the extension to
     *  determine whether the file is a legacy level set file and then copies all the other channel
     *  files if it is.
     *
     * @todo DEPRECATED
     *
     * @param  sourceFile	the location of the source file
     * @param  destFile 	the location of the destination file
     * @param  failFlag		whether to fail if the file exists at the source location, false will overwrite
     *
     */
    void copy_level_set_file( const frantic::tstring& sourceFile, const frantic::tstring& destFile,
                              const bool failFlag );
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
