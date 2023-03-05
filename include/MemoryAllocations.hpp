/*
* ==============================================================================
*
*  Copyright (C) 2021 Elliot Cartee, April Nellis, Jacob Van Hook
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
* ------------------------------------------------------------------------------
*
* File: MemoryAllocations.hpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This file contains helper functions for allocating and resizing
* multi-dimensional Boost arrays.
*
* ==============================================================================
*/

#ifndef BOOST_HELPER_HPP
#define BOOST_HELPER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include "boost/multi_array.hpp"

namespace memory {

/** Define template for 2D arrays */
template<class T>
using array2D_t = boost::multi_array<T, 2>;

/** Define template for 3D arrays */
template<class T>
using array3D_t = boost::multi_array<T, 3>;

/** Method for allocating 2D arrays */
template <class T>
array2D_t<T> allocateArray2D(const unsigned long n1, const unsigned long n2) {
  array2D_t<T> array(boost::extents[n1][n2]);
  return array;
}

/** Method for allocating 3D arrays */
template <class T>
array3D_t<T> allocateArray3D(const unsigned long n1, const unsigned long n2,
                             const unsigned long n3) {
  array3D_t<T> array(boost::extents[n1][n2][n3]);
  return array;
}

/** Method for initializing 2D arrays */
template <class T, std::size_t size1, std::size_t size2>
array2D_t<T> initializeArray2D(const T init_array[size1][size2]) {
  array2D_t<T> array(boost::extents[size1][size2]);
  for (unsigned long i = 0; i < size1; i++) {
    for (unsigned long j = 0; j < size2; j++) {
      array[i][j] = init_array[i][j];
    }
  }
  return array;
}

/** Method for initializing 3D arrays */
template <class T, std::size_t size1, std::size_t size2, std::size_t size3>
array3D_t<T> initializeArray3D(const T init_array[size1][size2][size3]) {
  array3D_t<T> array(boost::extents[size1][size2][size3]);
  for (unsigned long i = 0; i < size1; i++) {
    for (unsigned long j = 0; j < size2; j++) {
      for (unsigned long k = 0; k < size3; k++) {
          array[i][j][k] = init_array[i][j][k];
      }
    }
  }
  return array;
}

} // namespace memory

#endif