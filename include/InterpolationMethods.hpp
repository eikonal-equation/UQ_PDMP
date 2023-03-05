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
* File: InterpolationMethods.hpp
*
* Authors: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This file contains generic helper methods for interpolation.
*    Given a physical location and the grid sizes, this function interpolates
*    boost multi-arrays using multi-linear interpolation.
*
* ==============================================================================
*/

#ifndef INTERP_HPP
#define INTERP_HPP

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"

/** Interpolate 2D boost array using bilinear interpolation
  *    @param: aX, aY        physical location (aX,aY)
  *    @param: aDx, aDy      grid stepsizes in x and y directions
  *    @param: aMinX, aMinY  minimum physical x and values
  *    @param: aArray        pointer to 2D boost array to be interpolated */
double interpolate2D(const double aX, const double aY,
                     const double aDx, const double aDy,
                     const double aMinX, const double aMinY,
                     std::shared_ptr<memory::array2D_t<double>> aArray);

/** Interpolate 3D boost array using trilinear interpolation
  *    @param: aX, aY, aZ           physical location (aX,aY,aZ)
  *    @param: aDx, aDy, aDz        grid stepsizes in x, y, and z directions
  *    @param: aMinX, aMinY, aMinZ  minimum physical x and values
  *    @param: aArray               pointer to 3D boost array to be interpolated */
double interpolate3D(const double aX, const double aY, const double aZ,
                     const double aDx, const double aDy, const double aDz,
                     const double aMinX, const double aMinY, const double aMinZ,
                     std::shared_ptr<memory::array3D_t<double>> aArray);

#endif