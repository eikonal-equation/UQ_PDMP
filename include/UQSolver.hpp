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
* File: UQSolver.hpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description:  This header file is for the abstract class UQSolver. This is an
*   abstract class that provides an interface for both 1D and 2D solvers.
*
* ==============================================================================
*/

#ifndef ABSTRACT_UQ_HPP
#define ABSTRACT_UQ_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <vector>
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"

class UQSolver
{
  public:
    /*==========================================================================
        Transition Probabilities
    ==========================================================================*/
    /** Calculates the probability of transitioning from mode aI to mode aJ
      * after time aTau using fixed matrix of transition rates */
    double transitionProbability(const int aI, const int aJ,
                                 const double aTau) const;

    /** Calculates the probability of transitioning from mode aI to mode aJ
      * after time aTau using estimates of transition rates
      * @param aUpperBound: Determines which bound on lambdas is used
      *                       to calculate transition probability
                  aUpperBound = True:  Use upper bound on lambda value
                  aUpperBound = False: Use lower bound on lambda value */
    double transitionProbabilityEstimate(const int aI, const int aJ,
                                         const double aTau, 
                                         const bool aUpperBound) const;

    /*==========================================================================
        Monte Carlo Methods
    ==========================================================================*/
    /** Compute aNTrials Monte Carlo trials starting from aInitialLocation.
      *  @param aInitialPosition: vector containing starting location
      *  @param aInitialMode:     starting mode
      *  @param aNtrials:         number of Monte Carlo Trials
      *  @param aControlStrategy: choice of control policy to use in MC trials.
                      choices are "CDF", "EV", and "none" */
    virtual void monteCarlo(const std::vector<double> aInitialPosition,
                            const int aInitialMode, const int aNTrials,
                            const std::string aControlStrategy) {};

    /*==========================================================================
        PDE Numerics
    ==========================================================================*/
    /** Compute expectation-optimal policy and corresponding expected values
     *      for the piecewise deterministic Markov process.
     *  @param aUseFixedControl: flag for deciding whether to use fixed controls
     *                True:  use uniform zero controls
     *                False: chooses controls to optimize EV 
     *                         and places these controls in fModes */
    virtual void computeExpectedValue(const bool aUseFixedControl) {};

    /** Compute the CDF of the piecewise deterministic Markov process.
     *  @param aUseMinCost: flag for deciding whether to initialize with
     *                      minimum cost computations
     *  @param aControlStrategy: choice of control policy to use:
                      choices are "CDF", "EV", and "none"
                        CDF: choose controls to optimize CDF
                        EV: use EV-optimal controls
                        none: use uniformly zero controls */
    virtual void computeCDF(const bool aUseMinCost, 
                            const std::string aControlStrategy) {};

    /** Compute bounds on the CDFs of the piecewise deterministic Markov process.
     *  @param aControlStrategy: choice of control policy to use:
                      choices are "CDF", "EV", and "none"
                        CDF: choose controls to optimize CDF
                        EV: use EV-optimal controls
                        none: use uniformly zero controls */
    virtual void computeCDFBounds(const bool aUseMinCost,
                                  const std::string aControlStrategy,
                                  const bool aNumBound) {};

    /** Destructor */
    virtual ~UQSolver();

  protected:
    /** ------ Fields --------------------------------------------------------*/
    /** Number of modes */
    int fM;

    /** Lambda matrix of transition rates */
    std::shared_ptr<memory::array2D_t<double>> fLambdaMatrix;

    /** Estimated ranges for transition rates lambda */
    std::shared_ptr<memory::array3D_t<double>> fLambdaEstimates;

    /** A string which contains the prefix for the filenames 
     *  where the grids will be printed. */
    std::string fFilename;
};

#endif