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
* File: InterpolationMethods.cpp
*
* Authors: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This file contains generic helper methods for interpolation.
*    Given a physical location and the grid sizes, this function interpolates
*    boost multi-arrays using multi-linear interpolation.
*
* ==============================================================================
*/

#include "InterpolationMethods.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>

/** Interpolate 2D boost array using bilinear interpolation
  *    @param: aX, aY        physical location (aX,aY)
  *    @param: aDx, aDy      grid stepsizes in x and y directions
  *    @param: aMinX, aMinY  minimum physical x and values
  *    @param: aArray        pointer to 2D boost array to be interpolated */
double interpolate2D(const double aX, const double aY,
                     const double aDx, const double aDy,
                     const double aMinX, const double aMinY,
                     std::shared_ptr<memory::array2D_t<double>> aArray) {
  /** Find coordinates of "lower-left corner" */
  const int i = floor((aX - aMinX) / aDx);
  const int j = floor((aY - aMinY) / aDy);

  /** Find interpolation coefficients */
  // @TODO: Delete max step ??
  /** The max step avoids some issues with rounding making these coefficients
   *  negative when using very aggressive compiler optimizations */
  const double r_x = std::max(((aX - aMinX) / aDx) - i, 0.0);
  const double r_y = std::max(((aY - aMinY) / aDy) - j, 0.0);

  /** Find gridpoints for neighbors, checking bounds */
  const int nx = (*aArray).shape()[0];
  const int ny = (*aArray).shape()[1];
  const int i1 = std::min(i + 1, nx - 1);
  const int j1 = std::min(j + 1, ny - 1);

  /** Interpolate along x-axis */
  double top, bottom;
  /** If statements necessary to avoid rounding issues when optimizing controls */
  if ((*aArray)[i][j1] == (*aArray)[i1][j1]) {
    top = (*aArray)[i][j1];
  } else {
    top = (1-r_x) * (*aArray)[i][j1] + r_x * (*aArray)[i1][j1];
  }
  if ((*aArray)[i][j] == (*aArray)[i1][j]) {
    bottom = (*aArray)[i][j];
  } else {
    bottom = (1-r_x) * (*aArray)[i][j]  + r_x * (*aArray)[i1][j];
  }

  /** Interpolate along y-axis */
  return (1-r_y)*bottom + r_y*top;
}

/** Interpolate 3D boost array using trilinear interpolation
  *    @param: aX, aY, aZ           physical location (aX,aY,aZ)
  *    @param: aDx, aDy, aDz        grid stepsizes in x, y, and z directions
  *    @param: aMinX, aMinY, aMinZ  minimum physical x and values
  *    @param: aArray               pointer to 3D boost array to be interpolated */
double interpolate3D(const double aX, const double aY, const double aZ,
                     const double aDx, const double aDy, const double aDz,
                     const double aMinX, const double aMinY, const double aMinZ,
                     std::shared_ptr<memory::array3D_t<double>> aArray) {
  /* Find grid coordinates of "lower-left corner" */
  const int i = floor((aX - aMinX) / aDx);
  const int j = floor((aY - aMinY) / aDy);
  const int k = floor((aZ - aMinZ) / aDz);

  /* Find interpolation coefficients */
  // @TODO: Delete max step ??
  /** The max step avoids some issues with rounding making these coefficients
   *  negative when using very aggressive compiler optimizations */
  const double r_x = std::max(((aX - aMinX) / aDx) - i, 0.0);
  const double r_y = std::max(((aY - aMinY) / aDy) - j, 0.0);
  const double r_z = std::max(((aZ - aMinZ) / aDz) - k, 0.0);

  /* Find gridpoints for neighbors, checking bounds */
  const int nx = (*aArray).shape()[0];
  const int ny = (*aArray).shape()[1];
  const int nz = (*aArray).shape()[2];
  const int i1 = std::min(i + 1, nx - 1);
  const int j1 = std::min(j + 1, ny - 1);
  const int k1 = std::min(k + 1, nz - 1);

  /** Trilinear interpolation */
  /* Interpolate along s-axis */
  const double u1 = (1-r_z) * (*aArray)[i][j][k]   + (r_z) * (*aArray)[i][j][k1];
  const double u2 = (1-r_z) * (*aArray)[i1][j][k]  + (r_z) * (*aArray)[i1][j][k1];
  const double u3 = (1-r_z) * (*aArray)[i][j1][k]  + (r_z) * (*aArray)[i][j1][k1];
  const double u4 = (1-r_z) * (*aArray)[i1][j1][k] + (r_z) * (*aArray)[i1][j1][k1];

  /* Interpolate along y-axis */
  const double u5 = (1-r_y)*u1 + (r_y)*u3;
  const double u6 = (1-r_y)*u2 + (r_y)*u4;

  /* Interpolate along x-axis */
  return (1-r_x)*u5 + (r_x)*u6;
}