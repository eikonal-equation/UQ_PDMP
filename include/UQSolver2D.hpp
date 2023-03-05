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
* File: UQSolver2D.hpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description:  This header file is for the 2-dimensional UQ Solver.
*  This class is used for computing the CDFs of PDMPs.
*
* ==============================================================================
*/

#ifndef UQ_2D_HPP
#define UQ_2D_HPP

/** ------ Superclass --------------------------------------------------------*/
#include "UQSolver.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <vector>
#include <boost/heap/binomial_heap.hpp>
#include <random>

/** ------ Project-specific header files -------------------------------------*/
#include "WriteToFile.hpp"
#include "MemoryAllocations.hpp"
#include "CMode2D.hpp"

class UQSolver2D : public UQSolver
{
  public:
    /*==========================================================================
        Constructor
    ==========================================================================*/
    /** Default constructor */
    UQSolver2D() = default;

    /** Main constructor 
     *  @param aModes: vector containing object for each mode 
     *  @param aLambdaMatrix: 2D matrix of transition rates (lambda) 
     *  @param aFilename: prefix to be used when writing data to file */
    UQSolver2D(std::vector<CMode2D> aModes,
               std::shared_ptr<memory::array2D_t<double>> aLambdaMatrix,
               const std::string aFilename);

    /*==========================================================================
        Monte Carlo Method
    ==========================================================================*/
    /** Compute aNTrials Monte Carlo trials starting from aInitialLocation.
      *  @param aInitialPosition: vector containing starting location
      *  @param aInitialMode:     starting mode
      *  @param aNtrials:         number of Monte Carlo Trials
      *  @param aControlStrategy: choice of control policy to use in MC trials.
                      choices are "CDF", "EV", and "none" */
    void monteCarlo(const std::vector<double> aInitialPosition,
                    const int aInitialMode, const int aNTrials,
                    const std::string aControlStrategy) override;

    /*==========================================================================
        PDE Numerics
    ==========================================================================*/
    /** Compute the CDF of the piecewise deterministic Markov process in 1D.
     *  @param aUseMinCost: flag for deciding whether to initialize with
     *                      minimum cost computations
     *  @param aControlStrategy: choice of control policy to use .
                      choices are "CDF", "EV", and "none"
                        CDF: choose controls to optimize CDF
                        EV: use EV-optimal controls
                        none: use uniformly zero controls */
    void computeCDF(const bool aUseMinCost, 
                    const std::string aControlStrategy) override;

    /** Compute expectation-optimal policy and corresponding expected value.
     *  @param aUseFixedControl: flag for deciding whether to use fixed controls
     *                             True:  use uniform zero controls
     *                             False: chooses controls to optimize EV 
     *                                    and places these controls in fModes */
    void computeExpectedValue(const bool aUseFixedControl) override;

  private:
    /**========================================================================/
        Helper methods for CDF
    /=========================================================================*/
    /** Compute minimum cost and associated probabilities using Djikstra-like method
      *  @param aControlStrategy: choice of control policy to use.
                      choices are "CDF", "EV", and "none" */
    void compute_min_cost(const std::string aControlStrategy);

    /** Recursive implementation of golden section search (for CDF maximization)
     *    This is designed to include tiebreakers,
     *  @param f:    function to be maximized 
     *              should take in a double and output vector containing two doubles
     *              this method will then maximize first output
     *              while using secound output as a tiebreaker
     *  @param a,b:  endpoints of interval on which to use golden section search
     *  returns: maximizer x */
    double cdf_golden_section_search(const std::function<std::vector<double>(double)> f,
                                     double a, double b) const;

    /** Choose control to optimize CDF value
      *  @param aI: logical x-coordinate
      *  @param aJ: logical y-coordinate
      *  @param aN: logical s-coordinate
      *  @param aMode: current mode
      *  output: optimal (cdf,ev,control) */
    std::vector<double> optimize_cdf(const int aX, const int aY, const int aS,
                                     const int aMode) const;

    /** Update a single CDF value
      *  @param aX: physical x-coordinate
      *  @param aY: physical y-coordinate
      *  @param aS: physical s-coordinate
      *  @param aMode: current mode
      *  @param aControl: control value
      *  output: (cdf_value, ev_value) */
    std::vector<double> update_cdf(const int aI, const int aJ, const int aN,
                                   const int aMode, const double aControl) const;

    /**========================================================================/
        Helper methods for Expected Value
    /=========================================================================*/
    /** Chooses control to optimize EV at gridpoint (aI,aJ) and mode aMode */
    std::vector<double> optimize_ev(const int aI, const int aJ, 
                                    const int aMode) const;

    /** Compute EV update at gridpoint (aI,aJ) and mode aMode using aControl */
    double update_ev(const int aI, const int aJ, const int aMode, 
                     const double aControl) const;

    /** Recursive implementation of golden section search (for EV minimization)
     *  @param f:    function to be maximized
     *  @param a,b:  endpoints of interval
     *  returns: maximizer x */
    double ev_golden_section_search(const std::function<double(double)> f,
                                    double a, double b) const;

    /**========================================================================/
        Helper methods for Monte Carlo
    /=========================================================================*/
    /** Individual Monte Carlo trial.
     *  Given starting position and mode, returns accumulated cost. 
     *   @param aInitialPosition: vector containing starting location
     *   @param aInitialMode:     starting mode
     *   @param aDeadline:        control strategy will maximize probability of arriving
     *                            before aDeadline, then switch to minimizing EV */
    double monteCarloTrial(const std::vector<double> aInitialPosition,
                           const int aInitialMode,
                           const double aDeadline) const;

    /** Choose control to optimize CDF value
      *  @param aX: physical x-coordinate
      *  @param aY: physical y-coordinate
      *  @param aS: physical s-coordinate
      *  @param aMode:    current mode
      *  @param aTau:     timestep
      *  output: optimal (cdf,ev,control) */
    std::vector<double> mc_optimize(const double aX, const double aY,
                                    const double aS, const int aMode,
                                    const double aTau) const;

    /** Update a single CDF value at specified physical coordinates
      *  @param aX: physical x-coordinate
      *  @param aY: physical y-coordinate
      *  @param aS: physical s-coordinate
      *  @param aMode:    current mode
      *  @param aControl: control value
      *  @param aTau:     timestep
      *  output: (cdf_value, ev_value) */
    std::vector<double> mc_update(const double aX, const double aY,
                                  const double aS, const int aMode,
                                  const double aControl, const double aTau) const;

    /**========================================================================/
        Typedefs and helper methods for Djikstra's implementation in 2D
    /=========================================================================*/
    /** Helper function for updating gridpoint.
     *  Update EV at gridpoint (aI,aJ) using strategy aControlStrategy
     *    Returns True if value is updated,
     *    Returns False if value is not updated. */
    bool update_min_cost(const int aI, const int aJ, 
                         const std::string aControlStrategy);

    /* A GridPoint Class to be used by the heap. */
    class HeapGP {
      public:
        int i;
        int j;
        double value;
        HeapGP(int aI, int aJ, double aValue): i(aI), j(aJ), value(aValue) {};
    };

    /* Struct comparison which is required by boost::heap */
    struct compare_HeapGP
    {
      bool operator()(const HeapGP& point1, const HeapGP& point2) const
      {
        return point1.value > point2.value;
      }
    };

    /* Typedef Heap types to make it easier to read. */
    enum status_t {FAR, CONSIDERED, ACCEPTED};
    typedef boost::heap::binomial_heap<HeapGP, boost::heap::compare<compare_HeapGP> > heap_t;
    typedef typename boost::heap::binomial_heap<HeapGP, boost::heap::compare<compare_HeapGP> >::handle_type handle_t;

    /*==========================================================================
        Write to file
    ==========================================================================*/
    /** This method writes the grid sizes, step sizes, and lambda values to file */
    void write_config_to_file() const;

    /** An I/O member which writes the grid functions for each mode to file.
     * @param aVariables vector of variables to write to file.
     *      Choices are "CDF", "EV", "Velocity", and "Cost" */
    void write_modes_to_file(const std::vector<std::string> aVariables) const;

    /**========================================================================/
        Fields
    /=========================================================================*/
    /** Modes */
    std::vector<CMode2D> fModes;

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

    /** Arrays for s^0(x,i) and w^0(x,i) */
    std::shared_ptr<memory::array2D_t<status_t>> fStatus;
    std::shared_ptr<memory::array2D_t<handle_t>> fHeapPointers;
    std::shared_ptr<memory::array2D_t<double>> fMinCost;
    std::shared_ptr<memory::array3D_t<double>> fMinCostProb;
    std::shared_ptr<memory::array2D_t<double>> fMinCost_Controls;

    /** Array for u(x,i) */
    std::unique_ptr<memory::array3D_t<double>> fExpectedValue;

    /** Monte Carlo max steps */
    int fMaxSteps = 100000;

    /** Golden section search tolerance */
    double fTolGSS = 1.0e-3;
};

#endif