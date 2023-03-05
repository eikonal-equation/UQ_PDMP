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
* File: CMode1D.hpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This class is used to hold all the data associated with an
* individual mode with one space and one time dimension.
* This class is an interface used to access the underlying data structures.
* The CDF, controls, and expected values are stored as boost::multi_arrays,
* while the cost and speed functions are stored as function pointers.
*
* ==============================================================================
*/

#ifndef MODE1D_HPP
#define MODE1D_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"

class CMode1D
{
  protected:
    /** Grid parameters */
    int fNx;
    int fNs;
    double fDx;
    double fDs;
    double fMinX;
    double fMaxX;
    double fMaxS;

    /* Grid approximating CDF and CDF bounds.
     * Stored as 2D boost::multi_array */
    std::shared_ptr<memory::array2D_t<double>> fCDF;
    std::shared_ptr<memory::array2D_t<double>> fCDFUpperBound;
    std::shared_ptr<memory::array2D_t<double>> fCDFLowerBound;

    /** Arrays for storing controls */
    std::shared_ptr<memory::array2D_t<double>> fCDF_Control;
    std::vector<double>                        fEV_Control;

    /** Vector for holding expected value */
    std::shared_ptr<memory::array2D_t<double>> fExpectedValue;

    /** Running cost function, stored as a 2D function pointer
     *  Signature is C(x,a) */
    std::function<double(double,double)> fCost;

    /** Velocity function, stored as a 2D function pointer
     *  Signature is f(x,a) */
    std::function<double(double,double)> fVelocity;

  public:
    /** ========================================================================
    *    Constructor
    * ========================================================================*/
    /** Constructor for the CMode1D class.
     * @param aNx number of gridpoints in x-dimension (space)
     * @param aNs number of gridpoints in s-dimension (cost)
     * @param aMinX, aMaxX minimum and maximum physical value in x-dimension
     * @param aMaxS maximum physical value in s-dimension
     * @param aCostFunction  running cost function C(x,a)
     * @param aVelocity      velocity function f(x,a) */
    CMode1D(const int aNx, const int aNs,
            const double aMinX, const double aMaxX, const double aMaxS,
            std::function<double(double,double)> aCostFunction,
            std::function<double(double,double)> aVelocity);

    /** ========================================================================
    *    Setters
    *=========================================================================*/
    void setCDF(const int aK, const int aN, const double aValue);
    void setCDFUpperBound(const int aK, const int aN, const double aValue);
    void setCDFLowerBound(const int aK, const int aN, const double aValue);
    void setCDF_Control(const int aK, const int aN, const double aControl);
    void setEV_Control(const int aK, const double aControl);
    void setExpectedValue(const int aK, const int aN, const double aValue);

    void setCost(std::function<double(double,double)> aCostFunction);
    void setVelocity(std::function<double(double,double)> aVelocity);

    /** ========================================================================
    *    Getters
    * ========================================================================*/
    /** Logical-coordinate functions */
    double getCDF(const int aK, const int aN) const;
    double getCDFLowerBound(const int aK, const int aN) const;
    double getCDFUpperBound(const int aK, const int aN) const;
    double getCDF_Control(const int aK, const int aN) const;
    double getEV_Control(const int aK) const;
    double getExpectedValue(const int aK, const int aN) const;

    /** Physical-coordinate functions */
    double getCost(const double aX, const double aA) const;
    double getVelocity(const double aX, const double aA) const;

    /** Grid parameters */
    int getGridSizeX() const;
    int getGridSizeS() const;

    double getMinX() const;
    double getMaxX() const;
    double getMaxS() const;

    double getDx() const;
    double getDs() const;

    /** ========================================================================
    *    Interpolation Methods
    * ========================================================================*/
    /** Interpolate CDF at physical coordinates (aX,aS) */
    double interpolateCDF(const double aX, const double aS) const;

    /** Interpolate expected value at physical coordinates (aX,aS) */
    double interpolateEV(const double aX, const double aS) const;

    /** Interpolate CDF Upper Bound at physical coordinates (aX,aS) */
    double interpolateCDFUpperBound(const double aX, const double aS) const;

    /** Interpolate CDF Lower Bound at physical coordinates (aX,aS) */
    double interpolateCDFLowerBound(const double aX, const double aS) const;

    /** ========================================================================
    *    Writing to file
    * ========================================================================*/
    /** An I/O function which writes the grid variables to file.
     * @param aVariables vector of variables to write to file.
     *      Choices are "CDF", "EV", "CDFBounds", "Velocity", and "Cost"
     * @param aFilename a string which contains the prefix for the filenames
     *      where the grids will be printed. */
    void writeModeToFile(const std::vector<std::string> aVariables, 
                         const std::string aFilename) const;

    /** ========================================================================
    *    Other
    * ========================================================================*/
    /** Mapping from logical coordinates to physical coordinates */
    double xGridToPhysical(const int aK) const;
    double sGridToPhysical(const int aN) const;

    /** Allocating arrays
     * @param aVariables vector of variables to allocate arrays for.
     *      Choices are "CDF" and "CDFBounds" */
    void allocateVariables(const std::vector<std::string> aVariables);
};

/** ============================================================================
*    Inline function definitions
* ============================================================================*/
/** Inline definitions must be in the header file.
 * These functions are used frequently. */

/** ------ Inline definition of setters ------------------------------*/
inline void CMode1D::setCDF(const int aK, const int aN, const double aValue) {
  (*fCDF)[aK][aN] = aValue;
}

inline void CMode1D::setCDFUpperBound(const int aK, const int aN, const double aValue) {
  (*fCDFUpperBound)[aK][aN] = aValue;
}

inline void CMode1D::setCDFLowerBound(const int aK, const int aN, const double aValue) {
  (*fCDFLowerBound)[aK][aN] = aValue;
}

inline void CMode1D::setCDF_Control(const int aK, const int aN, const double aControl) {
  (*fCDF_Control)[aK][aN] = aControl;
}

inline void CMode1D::setEV_Control(const int aK, const double aControl) {
  fEV_Control[aK] = aControl;
}

inline void CMode1D::setExpectedValue(const int aK, const int aN, const double aValue) {
  (*fExpectedValue)[aK][aN] = aValue;
}

inline void CMode1D::setCost(std::function<double(double,double)> aCostFunction) {
  fCost = aCostFunction;
}

inline void CMode1D::setVelocity(std::function<double(double,double)> aVelocity) {
  fVelocity = aVelocity;
}


/** ------ Inline definition of getters --------------------------------------*/
inline double CMode1D::getCDF(const int aK, const int aN) const {
  return (*fCDF)[aK][aN];
}

inline double CMode1D::getCDFUpperBound(const int aK, const int aN) const {
  return (*fCDFUpperBound)[aK][aN];
}

inline double CMode1D::getCDFLowerBound(const int aK, const int aN) const {
  return (*fCDFLowerBound)[aK][aN];
}

inline double CMode1D::getCDF_Control(const int aK, const int aN) const {
  return (*fCDF_Control)[aK][aN];
}

inline double CMode1D::getEV_Control(const int aK) const {
  return fEV_Control[aK];
}

inline double CMode1D::getExpectedValue(const int aK, const int aN) const {
  return (*fExpectedValue)[aK][aN];
}

inline double CMode1D::getCost(const double aX, const double aA) const {
  return fCost(aX,aA);
}

inline double CMode1D::getVelocity(const double aX, const double aA) const {
  return fVelocity(aX,aA);
}

inline int CMode1D::getGridSizeX() const {
  return fNx;
}

inline int CMode1D::getGridSizeS() const {
  return fNs;
}

inline double CMode1D::getMinX() const {
  return fMinX;
}

inline double CMode1D::getMaxX() const {
  return fMaxX;
}

inline double CMode1D::getMaxS() const {
  return fMaxS;
}

inline double CMode1D::getDx() const {
  return fDx;
}

inline double CMode1D::getDs() const {
  return fDs;
}

/** Logical-to-physical coordinate transformation functions */
inline double CMode1D::xGridToPhysical(const int aK) const {
  return fMinX + (double)aK * fDx;
}

inline double CMode1D::sGridToPhysical(const int aN) const {
  return (double)aN * fDs;
}

#endif