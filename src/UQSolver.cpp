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
* File: UQSolver.cpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description:  This header file is for the abstract class UQSolver. This is an
*   abstract class that provides an interface for both 1D and 2D solvers.
*
* ==============================================================================
*/

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>

/** ------ Project-specific header files -------------------------------------*/
#include "UQSolver.hpp"

/** Destructor */
UQSolver::~UQSolver() = default;

/** Calculates the probability of transitioning from mode aI to mode aJ
  * after time aTau using fixed matrix of transition rates */
double UQSolver::transitionProbability(const int aI, const int aJ,
                                       const double aTau) const {
  if (aI != aJ) {
    return (*fLambdaMatrix)[aI][aJ] * aTau;
  } else {
    double prob = 1.0;
    for (int j1 = 0; j1 < fM; j1++) {
      if (j1 != aI) {
        prob = prob - transitionProbability(aI,j1,aTau);
      }
    }
    return prob;
  }
}

/** Calculates the probability of transitioning from mode aI to mode aJ
  * after time aTau using estimates of transition rates
  * @param aUpperBound: Determines which bound on lambdas is used
  *                       to calculate transition probability
              aUpperBound = True:  Use upper bound on lambda value
              aUpperBound = False: Use lower bound on lambda value */
double UQSolver::transitionProbabilityEstimate(const int aI, const int aJ,
                                               const double aTau,
                                               const bool aUpperBound) const {
  /** Choose whether to use upper or lower bound, 
   *  converting to int to use as index in array */
  int b;
  if (aUpperBound) {
    b = 1;
  } else {
    b = 0;
  }

  /** Calculate transition probability */
  if (aI != aJ) {
    return (*fLambdaEstimates)[b][aI][aJ] * aTau;
  } else {
    double prob = 1.0;
    for (int j1 = 0; j1 < fM; j1++) {
      if (j1 != aI) {
        prob = prob - transitionProbabilityEstimate(aI,j1,aTau,b);
      }
    }
    return prob;
  }
}