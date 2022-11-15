// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/named_channel_data.hpp>
#include <frantic/particles/particle_file_metadata.hpp>
#include <frantic/particles/streams/particle_istream.hpp>
#include <frantic/strings/tstring.hpp>

#include <boost/filesystem/path.hpp>

#include <utility>

namespace frantic {
namespace files {
namespace aln {

/**
 * Get the transforms in an ALN file. Returns a vector of (filename, transform) pairs.
 * This parser is stricter than MeshLab's. We define a valid ALN as one satisfying the following regular expression
 * (broken across lines for clarity):
 *
 *     \s*\d+\s*\n
 *     (
 *     [^\n]+\n
 *     (\s*#\s*\n)*
 *     (\s*(-?(\d+(\.\d*)?|\.\d+))(\s+(-?(\d+(\.\d*)?|\.\d+))){3}\s*\n){4}
 *     )*
 *     (\s*0\s*(\n.*)?)?
 *
 * where \s is one of: ' ', '\t', '\v', '\f', '\r'
 * In addition, the number in the first nonempty line must correspond to the number of scans.
 * However, we will allow an empty line or a line containing only whitespace at any point in the file.
 * For a concrete example, see "UnitTests/TestInputs/align.aln" and "UnitTests/TestInputs/variety.aln".
 * Note that the parser permits the following oddities that would not be produced by a MeshLab export:
 *
 *  + There can be any amount of "#" separators between a filename and its transform (even 0).
 *  + Anything can appear after the line indicating the end of file character "0".
 *  + An empty line or a line containing only whitespace at any point in the file.
 *  + Multiple whitespace characters will be truncated to one space, except at the start or end of a line, in which it
 *    will be removed.
 *  + Numbers can be specified in any form which can be boost:lexical_cast to MatrixT::float_type.
 *
 * @param <MatrixT> One of frantic::graphics::transform4f or frantic::graphics::transform4fd.
 * @param file The ALN file to get transforms from.
 * @param[out] outTransforms A std::map mapping scan file names to their transformation matrices.
 * @throws std::runtime_error If the file cannot be opened, or if the ALN is invalid.
 */
template <typename MatrixT>
void get_transforms( const frantic::tstring& file, std::vector<std::pair<frantic::tstring, MatrixT>>& outTransforms );

/**
 * Create a combined particle istream for the scans referenced in the ALN file.
 * @param file The ALN file containing registration data.
 * @param[out] outMetadata Metadata containing scanner transforms.
 * @param positionTypeHint The type hint to use when loading the scans.
 * @throws std::runtime_error
 */
boost::shared_ptr<frantic::particles::streams::particle_istream>
create_aln_particle_istream( const boost::filesystem::path& file,
                             frantic::particles::particle_file_metadata& outMetadata,
                             const frantic::channels::data_type_t positionTypeHint );
} // namespace aln
} // namespace files
} // namespace frantic
