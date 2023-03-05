/*
* ==============================================================================
*
*  Copyright (C) 2021  Elliot Cartee, April Nellis, Jacob Van Hook
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
* File: CMode2D.cpp
*
* Authors: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This class is used to hold all the data associated with an
* individual mode with two space and one time dimension.
* This class is an interface used to access the underlying data structures.
* The CDF, controls, and expected values are stored as boost::multi_arrays,
* while the cost and speed functions are stored as function pointers.
*
* (See also CMode2D.hpp)
*
* ==============================================================================
*/
#include "CMode2D.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <string>
#include <iostream>

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"
#include "WriteToFile.hpp"
#include "InterpolationMethods.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace memory;

/*==============================================================================
  Constructor
==============================================================================*/
/** Constructor for the CMode2D class.
 * @param aNx number of gridpoints in x-dimension (space)
 * @param aNy number of gridpoints in y-dimension (space)
 * @param aNs number of gridpoints in s-dimension (cost)
 * @param aMinX, aMaxX minimum and maximum physical value in x-dimension
 * @param aMinY, aMaxY minimum and maximum physical value in y-dimension
 * @param aMaxS maximum physical value in s-dimension
 * @param aCostFunction  running cost function C(x,y,a)
 * @param aVelocity      velocity function f(x,y,a) */
CMode2D::CMode2D(const int aNx, const int aNy, const int aNs,
                 const double aMinX, const double aMaxX, 
                 const double aMinY, const double aMaxY, const double aMaxS,
                 std::function<double(double,double,double)> aCostFunction,
                 std::function<double(double,double,double)> aVelocityX,
                 std::function<double(double,double,double)> aVelocityY) {
  /** Set grid parameters */
  fNx = aNx;
  fNy = aNy;
  fNs = aNs;
  fMinX = aMinX;
  fMaxX = aMaxX;
  fMinY = aMinY;
  fMaxY = aMaxY;
  fMaxS = aMaxS;
  fDx = (fMaxX - fMinX) / (fNx - 1);
  fDy = (fMaxY - fMinY) / (fNy - 1);
  fDs = fMaxS / (fNs - 1);

  /** Allocate arrays */
  fCDF           = std::make_shared<array3D_t<double>>(allocateArray3D<double>(fNx,fNy,fNs));
  fCDF_Control   = std::make_shared<array3D_t<double>>(allocateArray3D<double>(fNx,fNy,fNs));
  fEV_Control    = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));
  fExpectedValue = std::make_shared<array3D_t<double>>(allocateArray3D<double>(fNx,fNy,fNs));

  /** Speed and cost functions */
  fCost  = aCostFunction;
  fVelocityX = aVelocityX;
  fVelocityY = aVelocityY;
}

/*==============================================================================
  Write grid to file
==============================================================================*/
/** An I/O function which writes the grid variables to file.
 * @param aVariables vector of variables to write to file.
 *      Choices are "CDF", "EV", "Velocity", and "Cost"
 * @param aFilename a string which contains the prefix to the names
 *      of the files to which the grids will be printed. */
void CMode2D::writeModeToFile(const std::vector<std::string> aVariables,
                              std::string aFilename) const {
  for (unsigned long n = 0; n < aVariables.size(); n++) {
    const std::string variable = aVariables[n];
    if (variable == "CDF" && fCDF != NULL && fCDF_Control != NULL) {
      io::writeToFile3D<double>(aFilename + "_CDF", *fCDF);
      io::writeToFile3D<double>(aFilename + "_CDF_Control", *fCDF_Control);
    } else if (variable == "EV" && fExpectedValue != NULL && fEV_Control != NULL) {
      io::writeToFile3D<double>(aFilename + "_ExpectedValue", *fExpectedValue);
      io::writeToFile2D<double>(aFilename + "_EV_Control", *fEV_Control);
    } else if (variable == "Velocity") {
      array3D_t<double> velocity_x_array = allocateArray3D<double>(fNx,fNy,fNs);
      array3D_t<double> velocity_y_array = allocateArray3D<double>(fNx,fNy,fNs);
      for (int i = 0; i < fNx; i++) {
        for (int j = 0; j < fNy; j++) {
          for (int k = 0; k < fNs; k++) {
            const double x = xGridToPhysical(i);
            const double y = yGridToPhysical(j);
            const double a = getCDF_Control(i,j,k);
            velocity_x_array[i][j][k] = getVelocityX(x,y,a);
            velocity_y_array[i][j][k] = getVelocityY(x,y,a);
          }
        }
      }
      io::writeToFile3D<double>(aFilename + "_VelocityX", velocity_x_array);
      io::writeToFile3D<double>(aFilename + "_VelocityY", velocity_y_array);
    } else if (variable == "Cost") {
      array3D_t<double> cost_array = allocateArray3D<double>(fNx,fNy,fNs);
      for (int i = 0; i < fNx; i++) {
        for (int j = 0; j < fNy; j++) {
          for (int k = 0; k < fNs; k++) {
            const double x = xGridToPhysical(i);
            const double y = yGridToPhysical(j);
            const double a = getCDF_Control(i,j,k);
            cost_array[i][j][k] = getCost(x,y,a);
          }
        }
      }
      io::writeToFile3D<double>(aFilename + "_Cost", cost_array);
    }
  }

  return;
}

/*==============================================================================
  Interpolation
==============================================================================*/
/** Interpolate CDF at physical coordinates (x,y,s) using trilinear interpolation */
double CMode2D::interpolateCDF(const double aX, const double aY, const double aS) const {
  return interpolate3D(aX, aY, aS, fDx, fDy, fDs, fMinX, fMinY, 0.0, fCDF);
}

/** Interpolate Expected Value at physical coordinates (x,y,s) using trilinear interpolation */
double CMode2D::interpolateEV(const double aX, const double aY, const double aS) const {
  return interpolate3D(aX, aY, aS, fDx, fDy, fDs, fMinX, fMinY, 0.0, fExpectedValue);
}

/** Interpolate CDF Upper Bound at physical coordinates (x,y,s) using trilinear interpolation */
double CMode2D::interpolateCDFUpperBound(const double aX, const double aY, const double aS) const {
  return interpolate3D(aX, aY, aS, fDx, fDy, fDs, fMinX, fMinY, 0.0, fCDFUpperBound);
}

/** Interpolate CDF Lower Bound at physical coordinates (x,y,s) using trilinear interpolation */
double CMode2D::interpolateCDFLowerBound(const double aX, const double aY, const double aS) const {
  return interpolate3D(aX, aY, aS, fDx, fDy, fDs, fMinX, fMinY, 0.0, fCDFLowerBound);
}