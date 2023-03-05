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
* File: CMode2D.hpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This class is used to hold all the data associated with an
* individual mode with two space and one time dimension.
* This class is an interface used to access the underlying data structures.
* The CDF, controls, and expected values are stored as boost::multi_arrays,
* while the cost and speed functions are stored as function pointers.
*
* ==============================================================================
*/

#ifndef MODE2D_HPP
#define MODE2D_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"

class CMode2D
{
  private:
    /** Grid parameters */
    int fNx;
    int fNy;
    int fNs;
    double fDx;
    double fDy;
    double fDs;
    double fMinX;
    double fMinY;
    double fMaxX;
    double fMaxY;
    double fMaxS;

    /* Grid of CDF values, bounds, and expected values. 
     * Stored as 3D boost::multi_arrays */
    std::shared_ptr<memory::array3D_t<double>> fCDF;
    std::shared_ptr<memory::array3D_t<double>> fCDFUpperBound;
    std::shared_ptr<memory::array3D_t<double>> fCDFLowerBound;
    std::shared_ptr<memory::array3D_t<double>> fExpectedValue;

    /** Arrays for storing optimal controls */
    std::shared_ptr<memory::array3D_t<double>> fCDF_Control;
    std::shared_ptr<memory::array2D_t<double>> fEV_Control;

    /** Cost function, stored as a 3D function C(x,y,a) */
    std::function<double(double,double,double)> fCost;

    /** Horizontal and Vertical velocities, stored as 3D functions f(x,y,a) */
    std::function<double(double,double,double)> fVelocityX;
    std::function<double(double,double,double)> fVelocityY;

  public:
    /** ========================================================================
    *    Constructor
    * ========================================================================*/
    /** Constructor for the CMode2D class.
     * @param aNx number of gridpoints in x-dimension (space)
     * @param aNy number of gridpoints in y-dimension (space)
     * @param aNs number of gridpoints in s-dimension (cost)
     * @param aMinX, aMaxX minimum and maximum physical value in x-dimension
     * @param aMinY, aMaxY minimum and maximum physical value in y-dimension
     * @param aMaxS maximum physical value in s-dimension
     * @param aCostFunction  running cost function C(x,y,a)
     * @param aVelocity      velocity function f(x,y,a) */
    CMode2D(const int aNx, const int aNy, const int aNs,
            const double aMinX, const double aMaxX, 
            const double aMinY, const double aMaxY, const double aMaxS,
            std::function<double(double,double,double)> aCostFunction,
            std::function<double(double,double,double)> aVelocityX,
            std::function<double(double,double,double)> aVelocityY);

    /** ========================================================================
    *    Setters
    *=========================================================================*/
    void setCDF(const int aI, const int aJ, const int aK, const double aValue);
    void setCDFUpperBound(const int aI, const int aJ, const int aK, const double aValue);
    void setCDFLowerBound(const int aI, const int aJ, const int aK, const double aValue);
    void setCDF_Control(const int aI, const int aJ, const int aK, const double aValue);
    void setEV_Control(const int aI, const int aJ, const double aValue);
    void setExpectedValue(const int aI, const int aJ, const int aK, const double aValue);

    void setCost(std::function<double(double,double,double)> aCostFunction);
    void setVelocityX(std::function<double(double,double,double)> aVelocityX);
    void setVelocityY(std::function<double(double,double,double)> aVelocityY);

    /** ========================================================================
    *    Getters
    * ========================================================================*/
    /** Logical-coordinate functions */
    double getCDF(const int aI, const int aJ, const int aK) const;
    double getCDFUpperBound(const int aI, const int aJ, const int aK) const;
    double getCDFLowerBound(const int aI, const int aJ, const int aK) const;
    double getCDF_Control(const int aI, const int aJ, const int aK) const;
    double getEV_Control(const int aI, const int aJ) const;
    double getExpectedValue(const int aI, const int aJ, const int aK) const;

    /** Physical-coordinate functions */
    double getCost(const double aX, const double aY, const double aControl) const;
    double getVelocityX(const double aX, const double aY, const double aControl) const;
    double getVelocityY(const double aX, const double aY, const double aControl) const;

    /** Grid parameters */
    int getGridSizeX() const;
    int getGridSizeY() const;
    int getGridSizeS() const;

    double getMinX() const;
    double getMinY() const;
    double getMaxX() const;
    double getMaxY() const;
    double getMaxS() const;

    double getDx() const;
    double getDy() const;
    double getDs() const;

    /** ========================================================================
    *    Interpolation Methods
    * ========================================================================*/
    /** Interpolate CDF at physical coordinates (aX,aY,aS) */
    double interpolateCDF(const double aX, const double aY, const double aS) const;

    /** Interpolate Expected Value at physical coordinates (aX,aY,aS) */
    double interpolateEV(const double aX, const double aY, const double aS) const;

    /** Interpolate CDF Upper Bound at physical coordinates (aX,aY,aS) */
    double interpolateCDFUpperBound(const double aX, const double aY, const double aS) const;

    /** Interpolate CDF Lower Bound at physical coordinates (aX,aY,aS) */
    double interpolateCDFLowerBound(const double aX, const double aY, const double aS) const;

    /** ========================================================================
    *    Writing to file
    * ========================================================================*/
    /** An I/O function which writes the grid variables to file.
     * @param aVariables vector of variables to write to file.
     *      Choices are "CDF", "EV", "Velocity", and "Cost"
     * @param aFilename a string which contains the prefix for the filenames
     *      where the grids will be printed. */
    void writeModeToFile(const std::vector<std::string> aVariables, 
                         const std::string aFilename) const;

    /** ========================================================================
    *    Other
    * ========================================================================*/
    /** Mapping from logical coordinates to physical coordinates */
    double xGridToPhysical(const int aI) const;
    double yGridToPhysical(const int aJ) const;
    double sGridToPhysical(const int aK) const;
};

/** ============================================================================
*    Inline function definitions
* ============================================================================*/
/** Inline definitions must be in the header file.
 * These functions are used frequently. */

/** ------ Inline definition of setters ------------------------------*/
inline void CMode2D::setCDF(const int aI, const int aJ, const int aK, const double aValue) {
  (*fCDF)[aI][aJ][aK] = aValue;
}

inline void CMode2D::setCDFUpperBound(const int aI, const int aJ, const int aK, const double aValue) {
  (*fCDFUpperBound)[aI][aJ][aK] = aValue;
}

inline void CMode2D::setCDFLowerBound(const int aI, const int aJ, const int aK, const double aValue) {
  (*fCDFLowerBound)[aI][aJ][aK] = aValue;
}

inline void CMode2D::setCDF_Control(const int aI, const int aJ, const int aK, const double aControl) {
  (*fCDF_Control)[aI][aJ][aK] = aControl;
}

inline void CMode2D::setEV_Control(const int aI, const int aJ, const double aControl) {
  (*fEV_Control)[aI][aJ] = aControl;
}

inline void CMode2D::setExpectedValue(const int aI, const int aJ, const int aK, const double aValue) {
  (*fExpectedValue)[aI][aJ][aK] = aValue;
}

inline void CMode2D::setCost(std::function<double(double,double,double)> aCostFunction) {
  fCost = aCostFunction;
}

inline void CMode2D::setVelocityX(std::function<double(double,double,double)> aVelocityX) {
  fVelocityX = aVelocityX;
}

inline void CMode2D::setVelocityY(std::function<double(double,double,double)> aVelocityY) {
  fVelocityY = aVelocityY;
}

/** ------ Inline definition of getters --------------------------------------*/
inline double CMode2D::getCDF(const int aI, const int aJ, const int aK) const {
  return (*fCDF)[aI][aJ][aK];
}

inline double CMode2D::getCDFUpperBound(const int aI, const int aJ, const int aK) const {
  return (*fCDFUpperBound)[aI][aJ][aK];
}

inline double CMode2D::getCDFLowerBound(const int aI, const int aJ, const int aK) const {
  return (*fCDFLowerBound)[aI][aJ][aK];
}

inline double CMode2D::getCDF_Control(const int aI, const int aJ, const int aK) const {
  return (*fCDF_Control)[aI][aJ][aK];
}

inline double CMode2D::getEV_Control(const int aI, const int aJ) const {
  return (*fEV_Control)[aI][aJ];
}

inline double CMode2D::getExpectedValue(const int aI, const int aJ, const int aK) const {
  return (*fExpectedValue)[aI][aJ][aK];
}

inline double CMode2D::getCost(const double aX, const double aY, const double aControl) const {
  return fCost(aX,aY,aControl);
}

inline double CMode2D::getVelocityX(const double aX, const double aY, const double aControl) const {
  return fVelocityX(aX,aY,aControl);
}

inline double CMode2D::getVelocityY(const double aX, const double aY, const double aControl) const {
  return fVelocityY(aX,aY,aControl);
}

inline int CMode2D::getGridSizeX() const {
  return fNx;
}

inline int CMode2D::getGridSizeY() const {
  return fNy;
}

inline int CMode2D::getGridSizeS() const {
  return fNs;
}

inline double CMode2D::getMinX() const {
  return fMinX;
}

inline double CMode2D::getMinY() const {
  return fMinY;
}

inline double CMode2D::getMaxX() const {
  return fMaxX;
}

inline double CMode2D::getMaxY() const {
  return fMaxY;
}

inline double CMode2D::getMaxS() const {
  return fMaxS;
}

inline double CMode2D::getDx() const {
  return fDx;
}

inline double CMode2D::getDy() const {
  return fDy;
}

inline double CMode2D::getDs() const {
  return fDs;
}

/** Logical-to-physical coordinate transformation functions */
inline double CMode2D::xGridToPhysical(const int aI) const {
  return fMinX + (double)aI * fDx;
}

inline double CMode2D::yGridToPhysical(const int aJ) const {
  return fMinY + (double)aJ * fDy;
}

inline double CMode2D::sGridToPhysical(const int aK) const {
  return (double)aK * fDs;
}

#endif