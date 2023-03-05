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
*
* ------------------------------------------------------------------------------
*
* File: WriteToFile.hpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This file contains helper functions for writing multi-dimensional
* Boost arrays and vectors to file, and for printing nicely-formatted vectors to
* the command line.
*
* ==============================================================================
*/

#ifndef WRITE_TO_FILE_HPP
#define WRITE_TO_FILE_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <string>
#include <vector>
#include <fstream>

/** ----- Project-specific header files --------------------------------------*/
#include "MemoryAllocations.hpp"

namespace io {

/**
* This function writes the 2D Boost array aArray to a file with name aFilename
*/
template <class T>
void writeToFile2D(std::string aFilename, memory::array2D_t<T> aArray) {
  const unsigned long dim0 = aArray.shape()[0];
  const unsigned long dim1 = aArray.shape()[1];
  std::ofstream out("output/" + aFilename, std::ios::binary);

  for (unsigned long i = 0; i < dim0; i++) {
    for (unsigned long j = 0; j < dim1; j++) {
      out.write((char*) &aArray[i][j], sizeof(T));
    }
  }
}

/**
* This function writes the 3D Boost array aArray to a file with name aFilename
*/
template <class T>
void writeToFile3D(std::string aFilename, memory::array3D_t<T> aArray) {
  const unsigned long dim0 = aArray.shape()[0];
  const unsigned long dim1 = aArray.shape()[1];
  const unsigned long dim2 = aArray.shape()[2];
  std::ofstream out("output/" + aFilename, std::ios::binary);

  for (unsigned long i = 0; i < dim0; i++) {
    for (unsigned long j = 0; j < dim1; j++) {
      for (unsigned long k = 0; k < dim2; k++) {
        out.write((char*) &aArray[i][j][k], sizeof(T));
      }
    }
  }
}

/**
* This function writes the 1D std::vector aVec to a file with name aFilename
*/
template <class T>
void writeVectorToFile(std::string aFilename, std::vector<T> aVec) {
  std::ofstream out("output/" + aFilename);
  for (unsigned long i = 0; i < aVec.size(); i++) {
    out.write((char*) &aVec[i], sizeof(T));
  }
}

} // namespace io

#endif