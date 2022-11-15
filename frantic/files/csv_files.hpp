// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/strings/tstring.hpp>
#include <iostream>

namespace frantic {
namespace files {

// This function splits a single line of a .csv file into an array of entries.
// The streamName parameter should generally be the name of the input file, for error message purposes.
void read_csv_file_line( FILE* in, const frantic::tstring& streamName, std::vector<std::string>& outColumnData );

class line_reader_interface;

/**
 *  Split a CSV file into its records and entries.
 *
 *  NOTE: The input file may be buffered, so eof( file ) may return
 * true before the last record is returned.
 */
class csv_reader {
    char m_delimiter;
    std::string m_line;  // current line
    std::string m_token; // current token inside quotes
    bool m_emitEmptyColumns;

    boost::int64_t m_beginRowLineNumber;
    boost::int64_t m_lineNumber;

    boost::int64_t m_fileSize;

    bool m_useFileBuffer;
    boost::shared_ptr<line_reader_interface> m_lineReader;
    FILE* m_file;

    frantic::tstring m_streamName;

    bool is_delimiter( char c );
    void init( FILE* file, const frantic::tstring& streamName, bool useFileBuffer );
    bool getline( std::string& line );
    void do_emit( const std::string& s, std::vector<std::string>& outColumnData, std::size_t& outColumn );
    void maybe_emit( const std::string& s, std::vector<std::string>& outColumnData, std::size_t& outColumn );
    void maybe_emit( const std::string& s, std::string::size_type startIndex, std::string::size_type endIndex,
                     std::vector<std::string>& outColumnData, std::size_t& outColumn );

  public:
    /**
     *  Create an object that will read CSV files, with
     * optional buffering.
     *
     * @param file the file to read.
     * @param streamName the name of the stream, used for error
     *		messages.
     * @param useBuffer if true, the csv_reader will load the
     *		file into an internal buffer.
     */
    csv_reader( FILE* file = 0, const frantic::tstring& streamName = _T("<NONE>"), bool useBuffer = false );
    virtual ~csv_reader();
    /**
     *  Set the delimiter that splits columns in a record.
     */
    void set_delimiter( char c );
    /**
     *  Choose whether read_line should output empty columns that
     * are found outside of quotes.
     *
     * @param emitEmptyColumns if true, read_line may output
     *		empty columns.  This may happen if there are two
     *		delimiters in a row, or if there is a delimiter at
     *		the end of a line.
     */
    void set_emit_empty_columns( bool emitEmptyColumns );
    /**
     *  Read a line of the CSV file and split it into
     * columns.
     *
     * @param outColumnData if the function returns true,
     *		this is the columns in a record of the CSV file.
     * @return true if a line was read.  false if no line could
     *		be read, due to an error or end of file.
     */
    bool read_line( std::vector<std::string>& outColumnData );

    /**
     * @return the zero-based line number corresponding to the beginning of
     *		the most recent line returned by read_line(), or -1 if no line
     *		has been read from the file yet.
     */
    boost::int64_t get_line_number() const;

    /**
     * @return returns the size of the file
     */
    boost::int64_t file_progress_count() const;

    /**
     * @return returns the current position in the file
     */
    boost::int64_t file_progress_index() const;
};

/**
 *  Attempt to determine the CSV delimiter, that is, the character
 * that separates the columns.
 *
 * @param in a file opened for reading.  This must not be NULL.
 * @param streamName the name of the input file, used for error messages.
 * @param outDelimiter if the function returns true, then this is set to
 *		the delimiter.
 * @return true if a delimiter could be found, and outDelimiter is set.
 *		false otherwise.
 */
bool guess_csv_delimiter( FILE* in, const frantic::tstring& streamName, char& outDelimiter );

// This function builds a .csv file line out of an array of entries.
// NOTE: You might be tempted to just print them out with commas in between, but if any of the columns have certain
//       special characters, such as a comma, that will not work.  This function deals with such cases.
// void write_csv_file_line( const std::vector<std::string>& columnData, std::ostream& out );

class line_reader_interface {
  public:
    virtual ~line_reader_interface();

    /*
     * returns the number of bytes that have been read from the file
     */
    virtual boost::int64_t get_file_progress() = 0;

    /**
     *  Get a line of text from the input file.
     *
     * @param line if the function returns true, then this is the next
     *		line of text in the input file.
     * @return true is a line could be read, or false if an error occurred
     *		or the end of file was reached.
     */
    virtual bool getline( std::string& line ) = 0;
};

boost::shared_ptr<frantic::files::line_reader_interface> create_line_reader( FILE* file, bool useFileBuffer );

} // namespace files
} // namespace frantic
