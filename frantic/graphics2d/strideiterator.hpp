// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/* strideiterator.hpp
 *
 * This class provides a 2D iterator for a vector.
 *
 *  Written by: Brian McKinnon	Nov 7, 2006
 */

#ifndef STRIDE_ITERATOR_H
#define STRIDE_ITERATOR_H

#include <iostream>
#include <vector>

namespace frantic {
namespace graphics2d {

template <typename T>
class strideiterator {
  public:
    strideiterator( std::vector<T>* buf, int width, int height, int borderWidth = 0, int borderHeight = 0 );

    bool prevCol() {
        _colOffset--;
        return isColValid();
    }
    bool prevCol( int step ) {
        _colOffset -= step;
        return isColValid();
    }

    bool nextCol() {
        _colOffset++;
        return isColValid();
    }
    bool nextCol( int step ) {
        _colOffset += step;
        return isColValid();
    }

    bool prevRow() {
        _rowOffset -= _width;
        return isRowValid();
    }
    bool prevRow( int step ) {
        _rowOffset -= _width * step;
        return isRowValid();
    }

    bool nextRow() {
        _rowOffset += _width;
        return isRowValid();
    }
    bool nextRow( int step ) {
        _rowOffset += _width * step;
        return isRowValid();
    }

    bool nextRowStart() {
        _colOffset = _colOffsetFromHead;
        return nextRow();
    }
    bool prevRowStart() {
        _colOffset = _colOffsetFromHead;
        return prevRow();
    }

    bool nextColStart() {
        _rowOffset = _rowOffsetFromHead;
        return nextCol();
    }
    bool prevColStart() {
        _rowOffset = _rowOffsetFromHead;
        return prevCol();
    }

    bool goRow( int row ) {
        _rowOffset = row * _width;
        return isRowValid();
    }
    bool goCol( int col ) {
        _colOffset = col;
        return isColValid();
    }

    int getRowOffset() { return _rowOffset / _width; }
    int getColOffset() { return _colOffset; }

    bool isColValid() const;
    bool isRowValid() const;

    void gotoStart() {
        _rowOffset = _rowOffsetFromHead;
        _colOffset = _colOffsetFromHead;
    }

    void getData( T* data ) const { *data = ( *_data )[_rowOffset + _colOffset]; }
    void setData( T data ) { ( *_data )[_rowOffset + _colOffset] = data; }

    void setBuffer( T data ) {
        for( unsigned int i = 0; i < _data->size(); i++ )
            ( *_data )[i] = data;
    }
    void flipHorizontal();

    virtual void print( std::string header = "strideiterator info" ) const;

  private:
    std::vector<T>* _data;
    int _width;
    int _height;
    int _rowOffsetFromHead;
    int _rowOffsetFromTail;
    int _colOffsetFromHead;
    int _colOffsetFromTail;
    int _rowOffset;
    int _colOffset;
};

template <typename T>
strideiterator<T>::strideiterator( std::vector<T>* data, int width, int height, int borderWidth, int borderHeight ) {
    assert( data != NULL );
    assert( width > 0 );
    assert( height > 0 );
    _data = data;
    _width = width;
    _height = height;
    _rowOffsetFromHead = borderHeight * width;
    _rowOffsetFromTail = height * width - _rowOffsetFromHead;
    _colOffsetFromHead = borderWidth;
    _colOffsetFromTail = width - _colOffsetFromHead;
    _rowOffset = _rowOffsetFromHead;
    _colOffset = _colOffsetFromHead;
}

template <typename T>
bool strideiterator<T>::isColValid() const {
    if( _colOffset < _colOffsetFromHead )
        return false;
    else if( _colOffset < _colOffsetFromTail )
        return true;
    else
        return false;
}

template <typename T>
bool strideiterator<T>::isRowValid() const {
    if( _rowOffset < _rowOffsetFromHead )
        return false;
    else if( _rowOffset < _rowOffsetFromTail )
        return true;
    else
        return false;
}

template <typename T>
void strideiterator<T>::flipHorizontal() {
    T temp;

    for( int row = _rowOffsetFromHead; row < _rowOffsetFromTail; row += _width ) {
        for( int head = _colOffsetFromHead, tail = _colOffsetFromTail - 1; head < tail; head++, tail-- ) {
            temp = ( *_data )[row + head];
            ( *_data )[row + head] = ( *_data )[row + tail];
            ( *_data )[row + tail] = temp;
        }
    }
}

template <typename T>
void strideiterator<T>::print( std::string header ) const {
    std::cout << header << ":" << std::endl;
    std::cout << "\tData addr:        " << std::hex << _data << std::dec << std::endl;
    std::cout << "\t\tData size:          " << (unsigned int)_data->size() << std::endl;
    std::cout << "\tRow pointer:      " << _rowOffsetFromHead << " <- " << _rowOffset << " -> " << _rowOffsetFromTail
              << std::endl;
    std::cout << "\tColumn pointer:   " << _colOffsetFromHead << " <- " << _colOffset << " -> " << _colOffsetFromTail
              << std::endl;
    if( isRowValid() )
        std::cout << "\tRow is valid" << std::endl;
    else
        std::cout << "\tRow is not valid" << std::endl;

    if( isColValid() )
        std::cout << "\tColumn is valid" << std::endl;
    else
        std::cout << "\tColumn is not valid" << std::endl;
}
} // namespace graphics2d
} // namespace frantic

#endif
