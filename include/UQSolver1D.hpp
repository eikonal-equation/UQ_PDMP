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
* File: UQSolver1D.hpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description:  This header file is for the 1-dimensional UQ Solver.
*  This class implements all the numerical methods for solving PDEs, as well as
*  all Monte Carlo methods in one space dimension.
*
* ==============================================================================
*/

#ifndef UQ_1D_HPP
#define UQ_1D_HPP

/** ------ Superclass --------------------------------------------------------*/
#include "UQSolver.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <vector>
#include <boost/heap/fibonacci_heap.hpp>

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"
#include "CMode1D.hpp"

class UQSolver1D : public UQSolver
{
  public:
    /*==========================================================================
        Constructor
    ==========================================================================*/
    /** Default constructor */
    UQSolver1D() = default;

    /** Main constructor 
     *  @param aModes: vector containing object for each mode 
     *  @param aLambdaMatrix: 2D matrix of transition rates (lambda) 
     *  @param aLambdaEstimates: 3D matrix of bounds on transition rates (lambda) 
     *  @param aFilename: prefix to be used when writing data to file 
     *  @param aFishing:  boolean used to set whether solving fishing example or not */
    UQSolver1D(std::vector<CMode1D> aModes,
               std::shared_ptr<memory::array2D_t<double>> aLambdaMatrix,
               std::shared_ptr<memory::array3D_t<double>> aLambdaEstimates,
               const std::string aFilename, const bool aFishing);

    /*==========================================================================
        Monte Carlo Methods
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

    /** Compute bounds on the CDFs of the piecewise deterministic Markov process.
     *  @param aUseMinCost: flag for deciding whether to initialize with
     *                      minimum cost computations
     *  @param aControlStrategy: choice of control policy to use .
                      choices are "CDF", "EV", and "none"
                        CDF: choose controls to optimize CDF bounds
                            Note: Controls are chosen separately to optimize both
                                  upper and lower bounds, but only controls for
                                  optimizing lower bound are saved
                        EV: use EV-optimal controls
                        none: use uniformly zero controls
     *  @param aNumBound: flag for computing numerical or conservative bounds
     *                    True: use numerical bounds
     *                    False: use conservative analytical bounds */
    void computeCDFBounds(const bool aUseMinCost, 
                          const std::string aControlStrategy,
                          const bool aNumBound) override;

  private:
    /**========================================================================/
        Helper methods for CDF
    /=========================================================================*/
    /** Update a single CDF value at specified logical coordinates
      *  @param aX: physical x-coordinate
      *  @param aS: physical s-coordinate
      *  @param aMode: current mode
      *  @param aControl: control value
      *  output: (cdf_value, ev_value) */
    std::vector<double> update_cdf(const int aK, const int aN,
                                   const int aMode, 
                                   const double aControl) const;

    /** Update bounds for a single CDF value
      *  @param aK: logical x-coordinate
      *  @param aN: logical s-coordinate
      *  @param aMode: current mode
      *  @param aControl: control value
      *  @param aB: determines whether updating lower (0) or upper (1) bound
      *  @param aNumBound: flag for deciding whether to compute numerical or conservative bounds
      *                    True: use numerical bounds
      *                    False: use conservative analytical bounds */
    double update_cdf_bounds(const int aK, const int aN, const int aMode,
                             const double aControl, const int aB,
                             const bool aNumBound) const;

    /**========================================================================/
        Typedefs and helper methods for Minimum cost and
          Djikstra's implementation in 1D
    /=========================================================================*/
    /** Compute minimum cost and corresponding probabilities using Djikstra-like method
      *  @param aBounds - notes whether or not we are calculating bounds on CDF 
      *  @param aControlStrategy: choice of control policy to use.
                      choices are "CDF", "EV", and "none" */
    void compute_min_cost(const bool aBounds, const std::string aControlStrategy);

    /* A GridPoint Class to be used by the heap. */
    class HeapGP {
      public:
        int k;
        double value;
        HeapGP(int aK, double aValue): k(aK), value(aValue) {};
    };

    struct Neighbor {
      public:
        int k;
        int i;
        Neighbor(int aK, int aI): k(aK), i(aI) {};
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
    typedef boost::heap::fibonacci_heap<HeapGP, boost::heap::compare<compare_HeapGP> > heap_t;
    typedef typename boost::heap::fibonacci_heap<HeapGP, boost::heap::compare<compare_HeapGP> >::handle_type handle_t;

    /**========================================================================/
        Helper methods for Expected Value and Minimum cost
    /=========================================================================*/
    /** Update expected value at gridpoint (aK) and mode aI
     *    Returns True if value is updated,
     *    Returns False if value is not updated. */
    bool update_ev(const int aK, const int aI, const bool aUseFixedControl);

    /** Calculate update for minimum cost at gridpoint (aK)
      *  @param aK - logical x-coordinate of gridpoint to update
      *  @param bounds - notes whether or not we are calculating bounds on CDF
      *  @param aControlStrategy: choice of control policy to use.
                      choices are "CDF", "EV", and "none" */
    bool update_min_cost(const int aK, const bool aBounds, 
                         const std::string aControlStrategy);

    /**========================================================================/
        Helper methods for Monte Carlo
    /=========================================================================*/
    /** Individual Monte Carlo trial.
     *  Given starting position and mode, returns accumulated cost. 
     *   @param aInitialPosition: vector containing starting location
     *   @param aInitialMode:     starting mode
     *   @param aDeadline:        control strategy will maximize probability of arriving
     *                            before aDeadline, then switch to minimizing EV */
    double monteCarloTrial(const double aInitialPosition,
                           const int aInitialMode,
                           const double aDeadline) const;

    /** Update a single CDF value at specified physical coordinates
      *  @param aX: physical x-coordinate
      *  @param aS: physical s-coordinate
      *  @param aMode: current mode
      *  @param aControl: control value
      *  @param aTau: timestep
      *  output: (cdf_value, ev_value) */
    std::vector<double> mc_update(const double aX, const double aS, 
                                  const int aMode, const double aControl,
                                  const double aTau) const;

    /*==========================================================================
        Write to file
    ==========================================================================*/
    /** This method writes grid sizes, step sizes, and lambda values to file */
    void write_config_to_file() const;

    /** An I/O member which writes the grid functions for each mode to file.
     * @param aVariables vector of variables to write to file.
     *      Choices are "CDF", "EV", "CDFBounds", "Velocity", and "Cost" */
    void write_modes_to_file(const std::vector<std::string> aVariables) const;

    /*==========================================================================
        Allocating arrays
    ==========================================================================*/
    /** An I/O member which writes the grid functions for each mode to file.
     * @param aVariables vector of variables to write to file.
     *      Choices are "CDF" and "CDFBounds"*/
    void allocate_variables(const std::vector<std::string> aVariables);

    /**========================================================================/
        Fields
    /=========================================================================*/
    /** Modes */
    std::vector<CMode1D> fModes;

    /** Grid parameters */
    int fNx;
    int fNs;
    double fDx;
    double fDs;
    double fMinX;
    double fMaxX;
    double fMaxS;

    /** Arrays for s^0(x,i) and w^0(x,i) */
    std::vector<double> fMinCost;
    std::vector<double> fMinCost_Controls;
    std::shared_ptr<memory::array2D_t<double>> fMinCostProb;
    std::shared_ptr<memory::array2D_t<double>> fMinCostProbLower;
    std::shared_ptr<memory::array2D_t<double>> fMinCostProbUpper;

    /** Monte Carlo simulation results */
    int fMaxSteps = 100000;
    
    /** whether or not the fish harvesting example is being run */
    bool fFishing;
};

#endif