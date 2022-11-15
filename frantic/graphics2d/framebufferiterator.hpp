// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// This class provides a 2D iterator for an image.

#pragma once

#include "frantic/graphics2d/strideiterator.hpp"
#include <iostream>

namespace frantic {
namespace graphics2d {

template <class T>
class framebuffer;

template <typename T>
class framebufferiterator : public strideiterator<T> {
  public:
    framebufferiterator( framebuffer<T>* img, int borderWidth = 0, int borderHeight = 0 );

    void print( std::string header = "framebufferiterator info" ) const;

  private:
    framebuffer<T>* _img;
};

} // namespace graphics2d
} // namespace frantic

#include "frantic/graphics2d/framebuffer.hpp"

namespace frantic {
namespace graphics2d {

template <typename T>
framebufferiterator<T>::framebufferiterator( framebuffer<T>* img, int borderWidth, int borderHeight )
    : strideiterator<T>( &( img->data() ), img->xsize(), img->ysize(), borderWidth, borderHeight ) {
    assert( img != NULL );
    _img = img;
}

template <typename T>
void framebufferiterator<T>::print( std::string header ) const {
    std::cout << header << ":" << std::endl;
    std::cout << "\tframebuffer addr: " << std::hex << _img << std::dec << std::endl;
    std::cout << "\t\tDimensions:         ( " << _img->width() << ", " << _img->height() << " )" << std::endl;
    strideiterator<T>::print();
}

} // namespace graphics2d
} // namespace frantic
