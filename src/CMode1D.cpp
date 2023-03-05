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
* File: CMode1D.cpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This class is used to hold all the data associated with an
* individual mode with one space and one time dimension.
* This class is an interface used to access the underlying data structures.
* The CDF, controls, and expected values are stored as boost::multi_arrays,
* while the cost and speed functions are stored as function pointers.
*
* (See also CMode1D.hpp)
*
* ==============================================================================
*/
#include "CMode1D.hpp"

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
/** Constructor for the CMode1D class.
 * @param aNx number of gridpoints in x-dimension (space)
 * @param aNs number of gridpoints in s-dimension (cost)
 * @param aMinX, aMaxX minimum and maximum physical value in x-dimension
 * @param aMaxS maximum physical value in s-dimension
 * @param aCostFunction  running cost function C(x,a)
 * @param aVelocity      velocity function f(x,a) */
CMode1D::CMode1D(const int aNx, const int aNs,
                 const double aMinX, const double aMaxX, const double aMaxS,
                 std::function<double(double,double)> aCostFunction,
                 std::function<double(double,double)> aVelocity) {
  /** Set grid parameters */
  fNx = aNx;
  fNs = aNs;
  fMinX = aMinX;
  fMaxX = aMaxX;
  fMaxS = aMaxS;
  fDx = (fMaxX - fMinX) / (fNx - 1);
  fDs = fMaxS / (fNs - 1);

  /** Allocate arrays */
  fEV_Control    = std::vector<double>(fNx);
  fExpectedValue = std::make_unique<array2D_t<double>>(allocateArray2D<double>(fNx,fNs));

  /** Speed and cost functions */
  fCost     = aCostFunction;
  fVelocity = aVelocity;
}

/*==============================================================================
  Write grid to file
==============================================================================*/
/** An I/O function which writes the grid variables to file.
 * @param aVariables vector of variables to write to file.
 *      Choices are "CDF", "EV", "CDFBounds", "Velocity", and "Cost"
 * @param aFilename a string which contains the prefix to the names
 *      of the files to which the grids will be printed. */
void CMode1D::writeModeToFile(const std::vector<std::string> aVariables, 
                              const std::string aFilename) const {
  for (unsigned long n = 0; n < aVariables.size(); n++) {
    const std::string variable = aVariables[n];
    if (variable == "CDF" && fCDF != NULL && fCDF_Control != NULL) {
      io::writeToFile2D<double>(aFilename + "_CDF", *fCDF);
      io::writeToFile2D<double>(aFilename + "_CDF_Control", *fCDF_Control);
    } else if (variable == "EV" && fExpectedValue != NULL) {
      io::writeToFile2D<double>(aFilename + "_ExpectedValue", *fExpectedValue);
      io::writeVectorToFile<double>(aFilename + "_EV_Control", fEV_Control);
    } else if (variable == "CDFBounds" && fCDFUpperBound != NULL && fCDFLowerBound != NULL) {
      io::writeToFile2D<double>(aFilename + "_CDF_UpperBound", *fCDFUpperBound);
      io::writeToFile2D<double>(aFilename + "_CDF_LowerBound", *fCDFLowerBound);
    } else if (variable == "Velocity") {
      array2D_t<double> velocity_array = allocateArray2D<double>(fNx,fNs);
      for (int k = 0; k < fNs; k++) {
        for (int i = 0; i < fNx; i++) {
          const double x = xGridToPhysical(i);
          const double a = getCDF_Control(i,k);
          velocity_array[i][k] = getVelocity(x,a);
        }
      }
      io::writeToFile2D<double>(aFilename + "_Velocity", velocity_array);
    } else if (variable == "Cost") {
      array2D_t<double> cost_array  = allocateArray2D<double>(fNx,fNs);
      for (int k = 0; k < fNs; k++) {
        for (int i = 0; i < fNx; i++) {
          const double x = xGridToPhysical(i);
          const double a = getCDF_Control(i,k);
          cost_array[i][k] = getCost(x,a);
        }
      }
      io::writeToFile2D<double>(aFilename + "_Cost",  cost_array);
    } else {
      std::cout << "Warning: invalid argument used when writing to file" << std::endl;
    }
  }

  return;
}

/*==============================================================================
    Allocating arrays
==============================================================================*/
/** Helper function used for allocating arrays for variables
 * @param aVariables vector of variables to allocate arrays for.
 *      Choices are "CDF" and "CDFBounds" */
void CMode1D::allocateVariables(const std::vector<std::string> aVariables) {
  for (unsigned long n = 0; n < aVariables.size(); n++) {
    const std::string variable = aVariables[n];
    if (variable == "CDF") {
      fCDF         = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNs));
      fCDF_Control = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNs));
    } else if (variable == "CDFBounds") {
      fCDFUpperBound = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNs));
      fCDFLowerBound = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNs));
      if (fCDF_Control == NULL) {
        fCDF_Control   = std::make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNs));
      }
    }
  }
  return;
}

/*==============================================================================
  Interpolation
==============================================================================*/
/** Interpolate CDF at physical coordinates (aX,aS) */
double CMode1D::interpolateCDF(const double aX, const double aS) const {
  return interpolate2D(aX, aS, fDx, fDs, fMinX, 0.0, fCDF);
}

/** Interpolate Expected Value at physical coordinates (aX,aS) */
double CMode1D::interpolateEV(const double aX, const double aS) const {
  return interpolate2D(aX, aS, fDx, fDs, fMinX, 0.0, fExpectedValue);
}

/** Interpolate CDF Upper Bound at physical coordinates (aX,aS) */
double CMode1D::interpolateCDFUpperBound(const double aX, const double aS) const {
  return interpolate2D(aX, aS, fDx, fDs, fMinX, 0.0, fCDFUpperBound);
}

/** Interpolate CDF LowerBound at physical coordinates (aX,aS) */
double CMode1D::interpolateCDFLowerBound(const double aX, const double aS) const {
  return interpolate2D(aX, aS, fDx, fDs, fMinX, 0.0, fCDFLowerBound);
}